"""
Author: Eugene Tan
Date Created: 4/11/2024
Last Updated: 11/11/2024
Script to plot comparison between BV and MITN model outputs
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %%
using GeoIO
using GeoStats
using GeoInterface
using ProgressBars
using CSV
using DataFrames
using Missings
using JLD2
using DateConversions
using ProgressBars
using Distributions

# %% Plotting packages
using Plots
pythonplot()
theme(:vibrant)

# %% Custom Packages
using Theil

# %% Import datasets
# Map geometries
nat_geometries_dataset = GeoIO.load(raw"Z:\master_geometries\Admin_Units\Global\MAP\2023\MG_5K\admin2023_0_MG_5K.shp")
subnat_geometries_dataset = GeoIO.load(raw"Z:\master_geometries\Admin_Units\Global\MAP\2023\MG_5K\admin2023_1_MG_5K.shp")

# BV NPC Estimates
BV_dataset = CSV.read("datasets/subnational/BV_inla_npc_estimates.csv", DataFrame)
# %% Get ISO List

# Perform draws and save outputs. Filter out unwanted countries
ISO_list = String.(CSV.read(raw"C:\Users\ETan\Documents\Prototype Analyses\itn-updated\datasets\ISO_list.csv", DataFrame)[:,1])
exclusion_ISOs = ["CPV","ZAF"]#["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Make Storage variables
NAT_NPC_collection = Vector{Any}(undef, length(filt_ISOs))
NAT_ACCESS_collection = Vector{Any}(undef, length(filt_ISOs))
SUBNAT_NPC_collection = Vector{Any}(undef, length(filt_ISOs))
SUBNAT_ACCESS_collection = Vector{Any}(undef, length(filt_ISOs))
NPC_GAP_collection = Vector{Any}(undef, length(filt_ISOs))
NPC_GAP_normalised_collection = Vector{Any}(undef, length(filt_ISOs))
ACCESS_GAP_collection = Vector{Any}(undef, length(filt_ISOs))
ACCESS_GAP_normalised_collection = Vector{Any}(undef, length(filt_ISOs))
# NPC_VAR_collection = Vector{Any}(undef, length(filt_ISOs))
ACCESS_STD_collection = Vector{Any}(undef, length(filt_ISOs))
NAT_geometries_collection = Vector{Any}(undef, length(filt_ISOs))
SUBNAT_geometries_collection = Vector{Any}(undef, length(filt_ISOs))
THEIL_collection =  Vector{Any}(undef, length(filt_ISOs))

BV_SUBNAT_NPC_collection = Vector{Any}(undef, length(filt_ISOs))
BV_MITN_SUBNAT_NPC_GAP_collection = Vector{Any}(undef, length(filt_ISOs))
BV_MITN_SUBNAT_NPC_GAP_normalised_collection = Vector{Any}(undef, length(filt_ISOs))

SURVEY_NPC_collection = Vector{Any}(undef, length(filt_ISOs))
SURVEY_NPC_STD_collection = Vector{Any}(undef, length(filt_ISOs))

# %%
for ISO_i in ProgressBar(1:length(filt_ISOs))
    # %%
    ISO = filt_ISOs[ISO_i]
    YEAR_START = 2000
    YEAR_END = 2023

    # %% Load Draw data
    # National Estimates
    nat_dataset = JLD2.load("outputs/draws/national/crop_access/$(ISO)_$(YEAR_START)_$(YEAR_END)_post_crop_access.jld2")
    # Subnational Estimates
    subnat_dataset = JLD2.load("outputs/draws/subnational/$(ISO)_SUBNAT_draws.jld2")

    # %% Get number of regions
    n_admin1 = length(subnat_dataset["admin1_names"])

    # %% Get population data
    subnat_reg_dir = "outputs/regressions/subnational/"
    REG_YEAR_START = 2010 # YEAR CHOSEN FOR START OF SUBNATIONAL REGRESSION WHEN DOING MCMC DRAWS
    subnat_reg_filename = "$(ISO)_SUBNAT_NETCROP_$(YEAR_START)_$(REG_YEAR_START)_$(YEAR_END)_regression.jld2"
    subnat_reg_data = JLD2.load(subnat_reg_dir*subnat_reg_filename)

    # Extract population values
    # Get month index w.r.t to trained indices
    TRAINED_YEAR_START = subnat_dataset["YEAR_START_NAT"]
    year_idx_ref_1 = YEAR_START-TRAINED_YEAR_START+1
    year_idx_ref_2 = YEAR_END-TRAINED_YEAR_START+1

    monthidx_ref_1 = (year_idx_ref_1-1)*12 + 1
    monthidx_ref_2 = (year_idx_ref_2-1)*12 + 12

    subnat_monthidxs = monthidx_ref_1:monthidx_ref_2
    n_months = length(subnat_monthidxs)

    # Extract admin1 populations
    FULL_POPULATION_MONTHLY = zeros(n_admin1, n_months)

    for i in 1:n_admin1
        FULL_POPULATION_MONTHLY[i,:] = subnat_reg_data["admin1_outputs"][i]["FULL_POPULATION_MONTHLY"][subnat_monthidxs]
    end

    SUBNAT_POPULATIONS = zeros(n_admin1, n_months)
    for i in 1:n_admin1
        SUBNAT_POPULATIONS[i,:] = subnat_reg_data["admin1_outputs"][i]["FULL_POPULATION_MONTHLY"]
    end
    
    # Extracted adjusted subnational means and calculate gap

    SUBNAT_NPC_mean = zeros(n_admin1, n_months)
    SUBNAT_λ_access_mean = zeros(n_admin1, n_months)

    for i in 1:n_admin1
        # SUBNAT_NPC_mean[i,:] = mean(subnat_dataset["admin1_outputs"][i]["NPC_MONTHLY_TOTAL_samples"], dims = 1)[:]
        # SUBNAT_λ_access_mean[i,:] = mean(subnat_dataset["admin1_outputs"][i]["λ_ACCESS_samples"], dims = 1)[:]

        SUBNAT_NPC_mean[i,:] = mean(subnat_dataset["merged_outputs"][i]["ADJ_NPC_MONTHLY_TOTAL_samples"], dims = 1)[:]
        SUBNAT_λ_access_mean[i,:] = mean(subnat_dataset["merged_outputs"][i]["ADJ_λ_ACCESS_samples"], dims = 1)[:]
    end

    # %% Calculate National NPC using Subnational estimates
    NAT_NPC_mean = (sum(SUBNAT_NPC_mean .* FULL_POPULATION_MONTHLY, dims = 1)[:])./(sum(FULL_POPULATION_MONTHLY, dims = 1)[:])
    NAT_λ_access_mean = (sum(SUBNAT_λ_access_mean .* FULL_POPULATION_MONTHLY, dims = 1)[:])./(sum(FULL_POPULATION_MONTHLY, dims = 1)[:])
    
    # %% Extract gap values
    # NAT_A_NPC_mean_BYNET = nat_dataset["A_NPC_mean_BYNET"]
    # NAT_NPC_mean = sum(NAT_A_NPC_mean_BYNET, dims = (2,3))[:]
    # NAT_λ_access_mean = nat_dataset["λ_access_mean"]

    SUBNAT_NPC_GAP = SUBNAT_NPC_mean .- repeat(NAT_NPC_mean, 1,n_admin1)'
    SUBNAT_ACCESS_GAP = SUBNAT_λ_access_mean .- repeat(NAT_λ_access_mean, 1,n_admin1)'

    # Calculate normalised gap variation coefficient    
    NPC_GAP_mean = mean(SUBNAT_NPC_GAP, dims = 1)[:]
    ACCESS_GAP_mean = mean(SUBNAT_NPC_GAP, dims = 1)[:]

    # Calculate normalised Gap values
    SUBNAT_NPC_GAP_normalised = SUBNAT_NPC_GAP./repeat(NAT_NPC_mean', n_admin1, 1)
    SUBNAT_ACCESS_GAP_normalised = SUBNAT_ACCESS_GAP./repeat(NAT_λ_access_mean', n_admin1, 1)
    
   

    # Theil values
    THEIL_values = zeros(n_months)
    
    for month_i in 1:n_months
        THEIL_values[month_i] = NPC_to_theil(SUBNAT_NPC_mean[:,month_i], FULL_POPULATION_MONTHLY[:,month_i])
    end


    # NPC_VAR = (sqrt.(sum(((SUBNAT_NPC_GAP.- repeat(NPC_GAP_mean', n_admin1)).^2).*SUBNAT_POPULATIONS, dims = 1)./sum(SUBNAT_POPULATIONS, dims = 1))[:])
    ACCESS_STD = (sqrt.(sum(((SUBNAT_ACCESS_GAP.- repeat(ACCESS_GAP_mean', n_admin1)).^2).*SUBNAT_POPULATIONS, dims = 1)./sum(SUBNAT_POPULATIONS, dims = 1))[:])

    # %% Extract geometries
    # National
    NAT_geometries_collection[ISO_i] = nat_geometries_dataset[nat_geometries_dataset.ISO .== ISO,:][1,:].geometry

    # Subnational
    SUBNAT_geometries = Vector{Any}(undef, n_admin1)

    for i in 1:n_admin1
        area_id = subnat_dataset["merged_outputs"][i]["area_id"]
        admin1_geometry = subnat_geometries_dataset[subnat_geometries_dataset.area_id .== area_id,:][1,:].geometry
        SUBNAT_geometries[i] = admin1_geometry
    end

    # Get BV estimates for NPC
    # Storage variable
    BV_SUBNAT_NPC_mean = zeros(size(SUBNAT_NPC_mean))
    for i in 1:n_admin1
        area_id = subnat_dataset["merged_outputs"][i]["area_id"]
        # Import data from BV extraction
        BV_SUBNAT_NPC_ANNUAL = BV_dataset[BV_dataset.area_id .== area_id,:].NPC
        for year_i in 1:(size(BV_SUBNAT_NPC_ANNUAL)[1]-1)
            BV_SUBNAT_NPC_mean[i,12*(year_i-1)+1:12*year_i] = LinRange(BV_SUBNAT_NPC_ANNUAL[year_i],BV_SUBNAT_NPC_ANNUAL[year_i+1],13)[1:12]
        end
    end

    BV_MITN_SUBNAT_NPC_GAP = (SUBNAT_NPC_mean .- BV_SUBNAT_NPC_mean)
    BV_MITN_SUBNAT_NPC_GAP_normalised = BV_MITN_SUBNAT_NPC_GAP./BV_SUBNAT_NPC_mean

    # Extract survey estimates for npc
    survey_npcs = subnat_dataset["merged_outputs"][1]["SUBNAT_NPC_MONTHLY"]
    survey_npcs_std = subnat_dataset["merged_outputs"][1]["SUBNAT_NPC_STD_MONTHLY"]

    # %% Store in variable
    NAT_NPC_collection[ISO_i] = NAT_NPC_mean
    NAT_ACCESS_collection[ISO_i] = NAT_λ_access_mean
    SUBNAT_NPC_collection[ISO_i] = SUBNAT_NPC_mean
    SUBNAT_ACCESS_collection[ISO_i] = SUBNAT_λ_access_mean
    NPC_GAP_collection[ISO_i] = SUBNAT_NPC_GAP
    NPC_GAP_normalised_collection[ISO_i] = SUBNAT_NPC_GAP_normalised
    ACCESS_GAP_collection[ISO_i] = SUBNAT_ACCESS_GAP
    ACCESS_GAP_normalised_collection[ISO_i] = SUBNAT_ACCESS_GAP_normalised
    THEIL_collection[ISO_i] = THEIL_values
    ACCESS_STD_collection[ISO_i] = ACCESS_STD
    SUBNAT_geometries_collection[ISO_i] = SUBNAT_geometries

    BV_SUBNAT_NPC_collection[ISO_i] = BV_SUBNAT_NPC_mean
    BV_MITN_SUBNAT_NPC_GAP_collection[ISO_i] = BV_MITN_SUBNAT_NPC_GAP
    BV_MITN_SUBNAT_NPC_GAP_normalised_collection[ISO_i] = BV_MITN_SUBNAT_NPC_GAP_normalised

    SURVEY_NPC_collection[ISO_i] = survey_npcs
    SURVEY_NPC_STD_collection[ISO_i] = survey_npcs_std
end

# %%
BV_MITN_NPC_ests = zeros(0,2)

for i in 1:length(BV_SUBNAT_NPC_collection)
    NPC_ests = hcat(BV_SUBNAT_NPC_collection[i][:,1:12:end-12][:,11:end][:],SUBNAT_NPC_collection[i][:,1:12:end-12][:,11:end][:])
    BV_MITN_NPC_ests = vcat(BV_MITN_NPC_ests, NPC_ests)
end


# %%
xlims = (-0.05, 1.05)
ylims = (-0.05, 1.05)
fig1 = scatter(BV_MITN_NPC_ests[:,1], BV_MITN_NPC_ests[:,2], markersize = 4,
            xlims = xlims, ylims = ylims)

# %% NPC RMSE plots by Country
rmse_fig_collection = []
rmses = Matrix{Float64}(undef, length(filt_ISOs),2)
log_pdfs = Matrix{Float64}(undef, length(filt_ISOs),2)

for ISO_i in ProgressBar(1:length(filt_ISOs))
    ISO = filt_ISOs[ISO_i]
    println("$(ISO),$(ISO_i)")
    # extract required data for country
    survey_npc = SURVEY_NPC_collection[ISO_i]
    survey_npc_std = SURVEY_NPC_STD_collection[ISO_i]
    bv_npc = BV_SUBNAT_NPC_collection[ISO_i]
    subnat_npc = SUBNAT_NPC_collection[ISO_i]

    # Find all nonmissing entries from survey data
    nonmissing_idxs = findall(.!ismissing.(survey_npc))
    
    # Make plots
    fig = plot(0:0.01:1, 0:0.01:1, color = :grey, legend = :topleft, label = nothing,
                title = "Subnational\nNPC Fit RMSE ($(ISO))", xlabel = "Survey Estimate",
                ylabel = "Model Estimate",
                titlefontsize = 10, labelfontsize = 7, legendfontsize = 7, tickfontsize = 5)
    if !isempty(nonmissing_idxs)
        # Calculate RMSE
        mitn_rmse = sqrt(mean((subnat_npc[nonmissing_idxs] .- survey_npc[nonmissing_idxs]).^2))
        bv_rmse = sqrt(mean((bv_npc[nonmissing_idxs] .- survey_npc[nonmissing_idxs]).^2))

        # Calculate Log Likelihood
        mitn_logpdf = logpdf(MvNormal(Float64.(survey_npc[nonmissing_idxs]), Diagonal(survey_npc_std[nonmissing_idxs].^2)), subnat_npc[nonmissing_idxs])
        bv_logpdf = logpdf(MvNormal(Float64.(survey_npc[nonmissing_idxs]), Diagonal(survey_npc_std[nonmissing_idxs].^2)), bv_npc[nonmissing_idxs])

        # Save to variables
        rmses[ISO_i,:] = [mitn_rmse, bv_rmse]
        log_pdfs[ISO_i,:] = [mitn_logpdf, bv_logpdf]

        # Make Plots
        
        scatter!(fig, survey_npc[nonmissing_idxs], subnat_npc[nonmissing_idxs],
                    markerstrokewidth = 0, markercolor = 3, markeralpha = 0.8, markersize = 2.5,
                    xlims = xlims, ylims = ylims, label = "MITN $(round(mitn_rmse, digits  =3))")
        scatter!(fig, survey_npc[nonmissing_idxs], bv_npc[nonmissing_idxs],
                    markerstrokewidth = 0, markercolor = 4, markeralpha = 0.8, markersize = 2.5,
                    xlims = xlims, ylims = ylims, label = "BV $(round(bv_rmse, digits  =3))")
    end
    push!(rmse_fig_collection, fig)
end

# %%
country_rmse_plot = plot(rmse_fig_collection..., layout = (5,9), size = (1920,1080))
savefig(country_rmse_plot, "output_plots/subnational/all_countries_rmse.pdf")

# %%
sep = 0.15
ms = 4.5
sortidx = sortperm(rmses[:,1])
fig1 = plot(xlabel = "Country", ylabel = "NPC RMSE", title = "Subnational Model RMSE",
            xticks = (1:length(filt_ISOs), filt_ISOs[sortidx]),
            xtickfontrotation = 90, size = (900,300))
scatter!(fig1, (1:length(sortidx)) .- sep, rmses[sortidx,1], label = "MITN", markercolor = 3, markersize = ms)
scatter!(fig1, (1:length(sortidx)) .+ sep,rmses[sortidx,2], label = "BV", markercolor = 4, markersize = ms)

# %%
sep = 0.15
ms = 4.5
sortidx = sortperm(rmses[:,1])

fig2 = plot(xlabel = "Country", ylabel = "NPC RMSE", title = "Subnational Log PDFs",
            xticks = (1:length(filt_ISOs), filt_ISOs[sortidx]),
            xtickfontrotation = 90, size = (900,300), ylims = (-10, 1.05*maximum(log_pdfs[findall(.!isnan.(log_pdfs))])))
scatter!(fig2, (1:length(sortidx)) .- sep, log_pdfs[sortidx,1], label = "MITN", markercolor = 3, markersize = ms)
scatter!(fig2, (1:length(sortidx)) .+ sep,log_pdfs[sortidx,2], label = "BV", markercolor = 4, markersize = ms)

# %%
error_fig = plot(fig1, fig2, layout = (2,1), size = (900,500))
savefig(error_fig, "output_plots/subnational/model_errors.pdf")