"""
Author: Eugene Tan
Date Created: 4/11/2024
Last Updated: 12/11/2024
Script to plot posterior draws for each country and their subnational region to check fit
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import Public Packages
using JLD2
using CSV
using Plots
using DataFrames
using ProgressBars
using LinearAlgebra
using StatsBase

# # %% Import custom packages
using DateConversions
# using NetLoss
# using NetAccessModel
# using Subnat_NetCropModel

##############################################
# %% GLOBAL SETTINGS (NOT COUNTRY DEPENDENT)
##############################################
# %% Define paths
dataset_dir = "datasets/subnational/"
subnat_reg_dir = "outputs/regressions/subnational/"
nat_netcrop_post_dir = "outputs/draws/national/crop_access/"
nat_netage_post_dir = "outputs/draws/national/demography/"
nat_access_ext_dir = "outputs/extractions/access/"
nat_access_reg_dir = "outputs/regressions/access/"
save_dir = "output_plots/subnational/"

# %% Year bounds
REG_YEAR_START_NAT = 2000 # Start year for national model
REG_YEAR_START = 2010#2011 # Start year for subnational model
REG_YEAR_END = 2023

##############################################
# %% BATCH RUN CODE BLOCK!
##############################################
# %% Get ISO List
ISO_list = String.(CSV.read(raw"C:\Users\ETan\Documents\Prototype Analyses\itn-updated\datasets\ISO_list.csv", DataFrame)[:,1])
exclusion_ISOs = ["CPV","ZAF"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)


# %%
for ISO_i in ProgressBar(1:length(filt_ISOs))
    ##############################
    # %% Load Regressed Data
    ##############################
    ISO = filt_ISOs[ISO_i]
    merged_outputs = load("outputs/draws/subnational/"*"$(ISO)_SUBNAT_draws.jld2")["merged_outputs"]
    n_admin1 = length(merged_outputs)
    n_samples, n_months = size(merged_outputs[1]["NPC_MONTHLY_TOTAL_samples"])[1:2]

    ##############################################
    # %% COUNTRY SPECIFIC SETTINGS
    ##############################################
    # %% Data filenames
    # Reg filename
    subnat_reg_filename = "$(ISO)_SUBNAT_NETCROP_$(REG_YEAR_START_NAT)_$(REG_YEAR_START)_$(REG_YEAR_END)_regression.jld2"
    # National draws filename
    nat_netcrop_post_filename = "$(ISO)_2000_2023_post_crop_access.jld2"
    # Net Distribution data
    distributions_filename = "net_distributions_admin1_dummy_combined_amp.csv"
    # Net age demography posterior
    net_age_filename = "$(ISO)_net_age_demography_samples.csv"

    # %% National MCMC chain (for getting monthly disaggregation ratios)
    nat_cropchain_dir = "outputs/regressions/crop/Compact Regressions/$(REG_YEAR_START_NAT)_$(REG_YEAR_END)/"
    nat_cropchain_filename = "$(ISO)_$(REG_YEAR_START_NAT)_$(REG_YEAR_END)_cropchains.jld2"

    # %% Access Model MCMC Chain
    net_access_input_dict_filename = "$(REG_YEAR_START_NAT)_$(REG_YEAR_END)/$(ISO)_$(REG_YEAR_START_NAT)_$(REG_YEAR_END)_accessextract.jld2"
    net_access_chain_filename = "netaccesschains.jld2"

    ##############################################
    # %% Open relevant datasets
    ##############################################
    # %% Load Data
    # National net crop draw data
    nat_netcrop_post_data = JLD2.load(nat_netcrop_post_dir*nat_netcrop_post_filename)
    subnat_reg_data = JLD2.load(subnat_reg_dir*subnat_reg_filename)
    # Load distribution data
    master_distributions = CSV.read(dataset_dir*distributions_filename, DataFrame)
    # Load posterior net demography data
    master_net_age = CSV.read(nat_netage_post_dir*net_age_filename, DataFrame)
    # Access Model data
    net_access_input_dict = JLD2.load(nat_access_ext_dir*net_access_input_dict_filename)
    net_access_chain = JLD2.load(nat_access_reg_dir*net_access_chain_filename)
    # Household Survey Data
    nat_npc_monthly = CSV.read("datasets/npc_monthly_data.csv", DataFrame)
    subnat_npc_monthly = CSV.read("datasets/subnational/subnat_npc_monthly_data.csv", DataFrame)
    # National MCMC Chain (for getting monthly disaggregation ratios)
    nat_cropchain = load(nat_cropchain_dir*nat_cropchain_filename)

    ##############################################
    # %% Calculate relevant loop bounds (TIME INDEX)
    ##############################################

    # %% Get metadata
    # List of admin1 names
    admin1_names = subnat_reg_data["admin1_names"]

    # %% Calculate indexing bounds
    YEAR_START_NAT = 2000 #subnat_reg_data["YEAR_START_NAT"]
    YEAR_START = 2010 #subnat_reg_data["YEAR_START"]
    YEAR_END = 2023 #subnat_reg_data["YEAR_END"]

    YEARS_ANNUAL = Vector(YEAR_START:1:(YEAR_END))
    MONTHS_MONTHLY = Vector(1:(YEAR_END-YEAR_START+1)*12)
    FULL_MONTHS_MONTHLY = Vector(1:(YEAR_END-YEAR_START_NAT+1)*12)

    n_admin1 = length(admin1_names)
    n_months = length(FULL_MONTHS_MONTHLY)
    n_net_types = length(subnat_reg_data["admin1_outputs"][1]["NET_NAMES"])

    ##############################
    # %% COMPILE POPULATION DATA
    ##############################
    SUBNAT_POPULATION_MONTHLY = zeros(n_admin1, n_months)
    for admin1_name_i in 1:n_admin1
        SUBNAT_POPULATION_MONTHLY[admin1_name_i,:] = subnat_reg_data["admin1_outputs"][admin1_name_i]["FULL_POPULATION_MONTHLY"]
    end
    NAT_POPULATION_MONTHLY = sum(SUBNAT_POPULATION_MONTHLY, dims = 1)[:]

    ##############################
    # %% NATIONAL DATA PREP
    ##############################
    # Define Arrays to extract required Net Crop values
    NATIONAL_Γ_MONTHLY_TOTAL_samples = sum(nat_netcrop_post_data["Γ_MONTHLY_samples_BYNET"], dims = 3)[:,:,1]

    # Extract relevant data from household surveys
    filt_nat_npc_monthly = nat_npc_monthly[(nat_npc_monthly.ISO .== ISO) .&
                                            (nat_npc_monthly.year .>= YEAR_START_NAT) .&
                                            (nat_npc_monthly.year .<= YEAR_END), :]

    NAT_NET_CROP_MONTHLY = missings(Float64, n_months)
    NAT_NET_CROP_STD_MONTHLY = missings(Float64, n_months)

    for row_i in 1:size(filt_nat_npc_monthly)[1]
        month_val = filt_nat_npc_monthly[row_i,"month"]
        year_val = filt_nat_npc_monthly[row_i,"year"]
        monthidx = monthyear_to_monthidx(month_val, year_val, YEAR_START = YEAR_START_NAT)

        μ_est = filt_nat_npc_monthly[row_i,"NPC_mean"]*NAT_POPULATION_MONTHLY[monthidx]
        σ_est = filt_nat_npc_monthly[row_i,"NPC_adj_se"]*NAT_POPULATION_MONTHLY[monthidx]
        if (σ_est < μ_est) && σ_est>0
            NAT_NET_CROP_MONTHLY[monthidx] = μ_est
            NAT_NET_CROP_STD_MONTHLY[monthidx] = σ_est
        end
    end

    ##############################
    # %% SUBNATIONAL DATA PREP
    ##############################
    # Define Arrays to extract required Net Crop values, and household surveys
    ALLREGIONS_Γ_MONTHLY_TOTAL_samples = zeros(n_admin1, n_samples, n_months)

    SUBNAT_NET_CROP_MONTHLY = missings(Float64, n_admin1, n_months)
    SUBNAT_NET_CROP_STD_MONTHLY = missings(Float64, n_admin1, n_months)

    admin1_outputs = subnat_reg_data["admin1_outputs"]
    for admin1_name_i in 1:n_admin1

        # Extract relevant data from household surveys
        area_id = admin1_outputs[admin1_name_i]["area_id"]
        filt_subnat_npc_monthly = subnat_npc_monthly[(subnat_npc_monthly.ISO .== ISO) .&
                                        (subnat_npc_monthly.area_id .== area_id) .&
                                        (subnat_npc_monthly.year .>= YEAR_START_NAT) .&
                                        (subnat_npc_monthly.year .<= YEAR_END), :]
        
        for row_i in 1:size(filt_subnat_npc_monthly)[1]
            month_val = filt_subnat_npc_monthly[row_i,"month"]
            year_val = filt_subnat_npc_monthly[row_i,"year"]
            monthidx = monthyear_to_monthidx(month_val, year_val, YEAR_START = YEAR_START_NAT)
            μ_est = filt_subnat_npc_monthly[row_i,"NPC_mean"]*SUBNAT_POPULATION_MONTHLY[admin1_name_i, monthidx]
            σ_est = filt_subnat_npc_monthly[row_i,"NPC_adj_se"]*SUBNAT_POPULATION_MONTHLY[admin1_name_i, monthidx]
            if (σ_est < μ_est) && σ_est>0
                SUBNAT_NET_CROP_MONTHLY[admin1_name_i,monthidx] = μ_est
                SUBNAT_NET_CROP_STD_MONTHLY[admin1_name_i,monthidx] = σ_est
            end
        end
    end

    ##############################
    # %% Calculate NPC values
    ##############################

    SUBNAT_NPC_MONTHLY = missings(Float64, n_admin1, n_months)
    SUBNAT_NPC_STD_MONTHLY = missings(Float64, n_admin1, n_months)

    # Subnational: Get NPC for individual regions
    nonmissing_idxs = findall(.!ismissing.(SUBNAT_NET_CROP_MONTHLY))
    SUBNAT_NPC_MONTHLY[nonmissing_idxs] = SUBNAT_NET_CROP_MONTHLY[nonmissing_idxs]./SUBNAT_POPULATION_MONTHLY[nonmissing_idxs]
    SUBNAT_NPC_STD_MONTHLY[nonmissing_idxs] = SUBNAT_NET_CROP_STD_MONTHLY[nonmissing_idxs]./SUBNAT_POPULATION_MONTHLY[nonmissing_idxs]

    # Define Arrays to extract required Net Crop values
    ADJ_ALLREGIONS_Γ_MONTHLY_TOTAL_samples = zeros(n_admin1, n_samples, n_months)

    for admin1_name_i in 1:n_admin1
        # Look up posterior draws that were generated
        ADJ_ALLREGIONS_Γ_MONTHLY_TOTAL_samples[admin1_name_i,:,:,:] .= sum(merged_outputs[admin1_name_i]["Γ_MONTHLY_BYNET_samples"], dims = 3)[:,:,1]
    end

    ADJ_ALLREGIONS_Γ_MONTHLY_TOTAL_mean = mean(ADJ_ALLREGIONS_Γ_MONTHLY_TOTAL_samples, dims = 2)[:,1,:]
    ADJ_SUBNAT_NPC_MONTHLY = ADJ_ALLREGIONS_Γ_MONTHLY_TOTAL_mean./SUBNAT_POPULATION_MONTHLY
    # National
    NAT_REG_NAT_Γ_MONTHLY = mean(NATIONAL_Γ_MONTHLY_TOTAL_samples, dims = 1)[:]
    SUBNAT_REG_NAT_Γ_MONTHLY = sum(ADJ_ALLREGIONS_Γ_MONTHLY_TOTAL_mean, dims = 1)[:]

    NAT_REG_NAT_NPC_MONTHLY = NAT_REG_NAT_Γ_MONTHLY./NAT_POPULATION_MONTHLY
    SUBNAT_REG_NAT_NPC_MONTHLY = SUBNAT_REG_NAT_Γ_MONTHLY./NAT_POPULATION_MONTHLY

    ##############################
    # %% Plot and Check Adjustments （NET CROP)
    ##############################
    # General plot settings
    pythonplot()
    Plots.theme(:vibrant)
    YEARS_ANNUAL = REG_YEAR_START_NAT:REG_YEAR_END
    MONTHS_MONTHLY = 1:length(YEARS_ANNUAL)*12
    YEARS_ANNUAL[1]:YEARS_ANNUAL[end]

    # %% Check National Sums

    fig = plot(xlabel = "Years", ylabel = "Net Crop (mil)", title = "$(ISO) (National)", legend = :topleft,
    xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90)
    plot!(fig, NAT_REG_NAT_Γ_MONTHLY./1e6, label = "National Regression", linewidth = 1.5)
    plot!(fig, SUBNAT_REG_NAT_Γ_MONTHLY./1e6, label = "Subnational Regression", linewidth = 1.5)
    scatter!(fig, NAT_NET_CROP_MONTHLY./1e6, label = "National Survey", markercolor = 5, markersize = 4)

    save_subdir = save_dir*"netcrop_fits/"
    mkpath(save_subdir)
    savefig(fig, save_subdir*"$(ISO)_NAT_SUBNAT_TOTALS.pdf")

    ##############################
    # %% Plot and Check Adjustments （NPC)
    ##############################

    # # %% Check National Sums
    # fig = plot(xlabel = "Months", ylabel = "NPC", title = "$(ISO) (National)", legend = :topleft,
    #         xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90)
    # plot!(fig, NAT_REG_NAT_NPC_MONTHLY, label = "National Regression")
    # plot!(fig, SUBNAT_REG_NAT_NPC_MONTHLY, label = "Subnational Regression")

    # %% Check Subnational Sums
    admin1_figs_collection = []

    for admin1_name_i in 1:n_admin1
        fig = plot(xlabel = "Months", ylabel = "NPC", title = "$(ISO) ($(admin1_names[admin1_name_i]))", legend = false,
        xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90)
        plot!(fig, ADJ_SUBNAT_NPC_MONTHLY[admin1_name_i, :], label = "Subnational Regression", linewidth = 1.5)
        scatter!(fig, SUBNAT_NPC_MONTHLY[admin1_name_i,:], label = "Subnational Survey",
                    markerstrokewidth = 0, markercolor = 5)
        scatter!(fig, SUBNAT_NPC_MONTHLY[admin1_name_i,:].+2 .*SUBNAT_NPC_STD_MONTHLY[admin1_name_i,:], label = "Subnational Survey", 
                    markersize = 3.5, markerstrokewidth = 0, markercolor = 3)
        scatter!(fig, max.(0,SUBNAT_NPC_MONTHLY[admin1_name_i,:].-2 .*SUBNAT_NPC_STD_MONTHLY[admin1_name_i,:]), label = "Subnational Survey", 
                    markersize = 3.5, markerstrokewidth = 0, markercolor = 3)
        push!(admin1_figs_collection, fig)
    end

    # %%

    # Helper function to calculate ideal plot layout given number of subplots n_plots
    function calc_layout(n_plots)
        # n1 = floor(Int,sqrt(n_plots))
        # n2 = ceil(Int,sqrt(n_plots))
    
        # if n1*(n1-1) >= n_admin1
        #     return (n1-1, n1)
        # else
        #     return (n2-1, n2)
        # end
        n = ceil(Int, (1+sqrt(1+4*n_plots))/2)
        return return (n-1, n)
    end


    # %% Save plot
    fig = plot(admin1_figs_collection..., layout = calc_layout(n_admin1), size = (1920,1080))
    save_subdir = save_dir*"subnat_npc/"
    mkpath(save_subdir)
    savefig(fig, save_subdir*"$(ISO)_SUBNAT_NPC.pdf")

end
