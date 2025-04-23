"""
Author: Eugene Tan
Date Created: 4/11/2024
Last Updated: 12/11/2024
Generate subnational crop and access draws after training subnational crop model
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Import Public Packages
using JLD2
using CSV
using DataFrames
using ProgressBars
using LinearAlgebra
using StatsBase
using Turing
using AdvancedMH
using StatsPlots
using Distributions

# %% Import custom packages
using DateConversions
using NetLoss
using NetAccessModel
using NetAccessPrediction
using Subnat_NetCropModel

##############################################
# %% GLOBAL SETTINGS (NOT COUNTRY DEPENDENT)
##############################################
# %% Define paths
dataset_dir = RAW_SUBNAT_DATASET_DIR
subnat_reg_dir = OUTPUT_REGRESSIONS_DIR*"subnational/"
nat_netcrop_post_dir = OUTPUT_DRAWS_DIR*"national/crop_access/"
nat_netage_post_dir = OUTPUT_DRAWS_DIR*"national/demography/"
nat_access_ext_dir = OUTPUT_EXTRACTIONS_DIR*"access/pred_data/"
nat_access_reg_dir = OUTPUT_REGRESSIONS_DIR*"access/"
save_dir = OUTPUT_DRAWS_DIR*"subnational/"

# %% Sampling settings
n_samples = SUBNAT_CROP_ACCESS_N_DRAWS

# %% MCMC NAT-SUBNAT Adjustment ratio Sampler settings
iterations = NAT_SUBNAT_ADJ_MCMC_ITERATIONS
burn_in = NAT_SUBNAT_ADJ_MCMC_BURNIN
sample_var = NAT_SUBNAT_ADJ_MCMC_SAMPLING_VAR

# %% Year bounds
REG_YEAR_START_NAT = YEAR_NAT_START # Start year for national model
REG_YEAR_START = YEAR_SUBNAT_TRANS #2011 # Start year for subnational model
REG_YEAR_END = YEAR_NAT_END

##############################################
# %% BATCH RUN CODE BLOCK!
##############################################

# %% Perform draws and save outputs. Filter out unwanted countries
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)


# %%
for ISO in filt_ISOs
    ##############################################
    # %% COUNTRY SPECIFIC SETTINGS
    ##############################################
    # %% Data filenames
    # Reg filename
    subnat_reg_filename = "$(ISO)_SUBNAT_NETCROP_$(REG_YEAR_START_NAT)_$(REG_YEAR_START)_$(REG_YEAR_END)_regression.jld2"
    # National draws filename
    nat_netcrop_post_filename = "$(ISO)_2000_2023_post_crop_access.jld2"
    # Net Distribution data
    distributions_filename = SUBNAT_DISTRIBUTION_DATA_FILENAME
    # Net age demography posterior
    net_age_filename = "$(ISO)_net_age_demography_samples.csv"

    # %% National MCMC chain (for getting monthly disaggregation ratios)
    nat_cropchain_dir = OUTPUT_REGRESSIONS_DIR*"crop/$(REG_YEAR_START_NAT)_$(REG_YEAR_END)/"
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
    nat_npc_monthly = CSV.read(OUTPUT_DATAPREP_DIR*HOUSEHOLD_NAT_SUMMARY_DATA_FILENAME, DataFrame)
    subnat_npc_monthly = CSV.read(OUTPUT_SUBNAT_DATAPREP_DIR*HOUSEHOLD_SUBNAT_SUMMARY_DATA_FILENAME, DataFrame)
    # National MCMC Chain (for getting monthly disaggregation ratios)
    nat_cropchain = load(nat_cropchain_dir*nat_cropchain_filename)

    ##############################################
    # %% Calculate relevant loop bounds (TIME INDEX)
    ##############################################

    # %% Get metadata
    # List of admin1 names
    admin1_names = subnat_reg_data["admin1_names"]
    n_admin1 = length(admin1_names)

    # %% Calculate indexing bounds
    YEAR_START_NAT = YEAR_NAT_START #subnat_reg_data["YEAR_START_NAT"]
    YEAR_START = YEAR_SUBNAT_TRANS #subnat_reg_data["YEAR_START"]
    YEAR_END = YEAR_NAT_END #subnat_reg_data["YEAR_END"]

    YEARS_ANNUAL = Vector(YEAR_START:1:(YEAR_END))
    MONTHS_MONTHLY = Vector(1:(YEAR_END-YEAR_START+1)*12)
    FULL_MONTHS_MONTHLY = Vector(1:(YEAR_END-YEAR_START_NAT+1)*12)

    # %% Get required dimension sizes for sampling posterior and adjusted estimates
    n_admin1 = length(admin1_names)
    n_months = length(FULL_MONTHS_MONTHLY)
    n_net_types = length(subnat_reg_data["admin1_outputs"][1]["NET_NAMES"])

    
    ##############################################
    # %% Extract relevant variables
    ##############################################
    # %% Get monthly disaggregation ratios and re-use for subnational regression
    full_monthly_p = nat_cropchain["monthly_p"]

    # %% Import National Access model parameters
    ρ_chain_df = net_access_chain["ρ_chain_df"]
    μ_chain_df = net_access_chain["μ_chain_df"]
    p_h = mean(net_access_input_dict["p_h_aggregated"], dims = 1)[:]

    # %% Extract Net Distribution Numbers
    ALLREGIONS_DISTRIBUTION_ANNUAL_BYNET = zeros(n_admin1, length(YEARS_ANNUAL), n_net_types)
    for i in 1:n_admin1
        # ALLREGIONS_DISTRIBUTION_ANNUAL_BYNET[i,:,:] = subnat_reg_data["admin1_outputs"][i]["DISTRIBUTION_ANNUAL_BYNET"]
        ALLREGIONS_DISTRIBUTION_ANNUAL_BYNET[i,:,:] = subnat_reg_data["admin1_outputs"][i]["ADJ_DISTRIBUTION_ANNUAL_BYNET"]
    end
    

    ##############################################
    # %% Subnational Crop and Access Draws
    ##############################################

    # %% Get number of Admin 1 Subnational regions
    n_admin1 = length(admin1_names)

    # %% Storage variable for results
    admin1_outputs = Vector{Any}(undef, n_admin1)

    ##############################################
    # UNADJUSTED POSTERIOR DRAWS FOR EACH SUBNATIONAL REGION
    # %% Do draws for forward evolution
    ##############################################
    for admin1_name_i in 1:length(admin1_names)
        println("Sampling for $(ISO), region $(admin1_name_i) of $(length(admin1_names))...")

        # Get admin1 parameters
        admin1_name = admin1_names[admin1_name_i]
        area_id = subnat_reg_data["admin1_outputs"][admin1_name_i]["area_id"]
        NET_NAMES = subnat_reg_data["admin1_outputs"][admin1_name_i]["NET_NAMES"]

        # Get trimmed month index range where subnational is relevant
        
        idx_year_ref_1 = YEAR_START - REG_YEAR_START_NAT + 1
        idx_year_ref_2 = YEAR_END - REG_YEAR_START_NAT + 1

        idx_month_ref_1 = (idx_year_ref_1-1)*12+1
        idx_month_ref_2 = (idx_year_ref_2*12)
        subnat_monthidxs = idx_month_ref_1:idx_month_ref_2 

        monthly_p = full_monthly_p[subnat_monthidxs]

        # %% Forward simulation

        FULL_POPULATION_MONTHLY = subnat_reg_data["admin1_outputs"][admin1_name_i]["FULL_POPULATION_MONTHLY"]
        FULL_POPULATION_ANNUAL = FULL_POPULATION_MONTHLY[1:12:end]
        YEAR_START_idx = YEAR_START-YEAR_START_NAT+1
        POPULATION_ANNUAL = FULL_POPULATION_ANNUAL[YEAR_START_idx:end]

        # Get required net distribution data ( from subnational regression)


        # filt_distributions = master_distributions[(master_distributions.iso .== ISO) .& 
        #                                     (master_distributions.area_id .== area_id) .& 
        #                                     (master_distributions.year .>= YEAR_START) .&
        #                                     (master_distributions.year .<= YEAR_END),:]
        # DISTRIBUTION_ANNUAL_BYNET = Matrix(filt_distributions[sortperm(filt_distributions.year),NET_NAMES])

        DISTRIBUTION_ANNUAL_BYNET =  ALLREGIONS_DISTRIBUTION_ANNUAL_BYNET[admin1_name_i,:,:]
        
        # Get required net demography
        initial_demography_condition = master_net_age[(master_net_age.ISO .== ISO) .&
                                                        (master_net_age.YEAR .== YEAR_START) .&
                                                        (master_net_age.MONTH .== 1),:]

        # Load national demography data for combining final net crop draws
        # Extract original net age demography and scale by required population proportion size
        FULL_A_NPC_BYNET = nat_netcrop_post_data["A_NPC_mean_BYNET"]
        FULL_A_BYNET = zeros(size(FULL_A_NPC_BYNET))

        for birth_j in 1:size(FULL_A_BYNET)[2]
            for net_type_i in 1:size(FULL_A_BYNET)[3]
                FULL_A_BYNET[:,birth_j,net_type_i] = FULL_A_NPC_BYNET[:,birth_j,net_type_i].*FULL_POPULATION_MONTHLY
            end
        end

        # Also load raw net-crop for use in de-meaning net subnational estimates prior to YEAR_START
        Γ_MONTHLY_samples_BYNET = nat_netcrop_post_data["Γ_MONTHLY_samples_BYNET"]
        Γ_MONTHLY_mean_BYNET = mean(Γ_MONTHLY_samples_BYNET, dims = 1)[1,:,:]
        
        # Get number of net types and missing net distribution entries
        n_net_types = length(NET_NAMES)
        n_missing = sum(ismissing.(DISTRIBUTION_ANNUAL_BYNET[:,1]))*n_net_types

        # Import MCMC sampled chains
        τ_chain = subnat_reg_data["admin1_outputs"][admin1_name_i]["τ_chain"]
        κ_chain = subnat_reg_data["admin1_outputs"][admin1_name_i]["κ_chain"]
        missing_nets_chain = subnat_reg_data["admin1_outputs"][admin1_name_i]["missing_nets_chain"]
        subnat_reg_data["admin1_outputs"][admin1_name_i]

        # Sample from MCMC joint posterior
        sampled_idxs = rand(1:size(τ_chain)[1], n_samples)

        # Define storage variables for sample draws
        COMBINED_A_BYNET_samples = zeros(n_samples, size(FULL_A_BYNET)...)
        Γ_MONTHLY_BYNET_samples = zeros(n_samples, size(FULL_A_BYNET)[1], n_net_types)
        NPC_MONTHLY_BYNET_samples = zeros(n_samples, size(FULL_A_BYNET)[1], n_net_types)

        COMBINED_A_TOTAL_samples = zeros(n_samples, size(FULL_A_BYNET)[1:2]...)
        Γ_MONTHLY_TOTAL_samples = zeros(n_samples, size(FULL_A_BYNET)[1])
        NPC_MONTHLY_TOTAL_samples = zeros(n_samples, size(FULL_A_BYNET)[1])

        # Draw from subnational crop posterior
        for i in ProgressBar(1:n_samples, leave = false)
            τ_est = Vector(τ_chain[sampled_idxs[i],:])
            κ_est = Vector(κ_chain[sampled_idxs[i],:])

            missing_nets_est = []

            if n_missing > 0
                missing_nets_est = Vector(missing_nets_chain[sampled_idxs[i],:])
            end 

            # Forward simulation using drawn posterior parameters
            rand_init_demography_idx = rand(unique(initial_demography_condition[:,4]))

            rand_init_demography_idx
            init_demography_sample = initial_demography_condition[initial_demography_condition.sample_idx .== rand_init_demography_idx,:]
            initial_demography_condition    
            
            # Use posterior mean final net demography condition at YEAR_START from national model to initialise
            A_BYNET = SUBNAT_monthly_itn_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                                            POPULATION_ANNUAL,
                                                            DISTRIBUTION_ANNUAL_BYNET,
                                                            NET_NAMES,
                                                            init_demography_sample,
                                                            τ_est, κ_est,
                                                            missing_nets_est,
                                                            monthly_p = monthly_p)[2]

            # Join national and subnational draws at YEAR_START
            # Need to first "de-mean" FULL_A_BYNET by multiplying the ratio of Γ_MONTHLY_samples_BYNET/Γ_MONTHLY_mean_BYNET
            Γ_MONTHLY_sample = Γ_MONTHLY_samples_BYNET[rand(1:size(Γ_MONTHLY_samples_BYNET)[1]),:,:]
            ADJ_ratio = Γ_MONTHLY_sample./Γ_MONTHLY_mean_BYNET
            ADJ_ratio[findall(isnan.(ADJ_ratio))] .= 1 # Clean cases where there are NaNs
            
            ADJ_FULL_A_BYNET = copy(FULL_A_BYNET)
            for month_i in 1:size(ADJ_FULL_A_BYNET)[1]
                for net_type_i in 1:n_net_types
                    ADJ_FULL_A_BYNET[month_i,:,net_type_i] = FULL_A_BYNET[month_i,:,net_type_i]*ADJ_ratio[month_i,net_type_i]
                end
            end
            
            # Replace national values with subnational draws
            COMBINED_A_BYNET = copy(ADJ_FULL_A_BYNET)
            
            COMBINED_A_BYNET[end-size(A_BYNET)[1]+2:end,end-size(A_BYNET)[2]+2:end,:] .= copy(A_BYNET[1:end-1,1:end-1,:])

            # Calculate summaries
            Γ_MONTHLY_BYNET = sum(COMBINED_A_BYNET, dims = 2)[:,1,:]
            NPC_MONTHLY_BYNET = zeros(size(Γ_MONTHLY_BYNET))
            for net_type_i in 1:n_net_types
                NPC_MONTHLY_BYNET[:,net_type_i] = Γ_MONTHLY_BYNET[:,net_type_i]./FULL_POPULATION_MONTHLY
            end
            
            COMBINED_A_TOTAL = sum(COMBINED_A_BYNET, dims = 3)[:,:,1]
            Γ_MONTHLY_TOTAL = sum(COMBINED_A_TOTAL, dims = 2)[:,1]
            NPC_MONTHLY_TOTAL = Γ_MONTHLY_TOTAL./FULL_POPULATION_MONTHLY
            Plots.plot(NPC_MONTHLY_TOTAL)
            # Store sample in storage variable
            COMBINED_A_BYNET_samples[i,:,:,:] = COMBINED_A_BYNET
            Γ_MONTHLY_BYNET_samples[i,:,:] = Γ_MONTHLY_BYNET
            NPC_MONTHLY_BYNET_samples[i,:,:] = NPC_MONTHLY_BYNET
            COMBINED_A_TOTAL_samples[i,:,:] = COMBINED_A_TOTAL
            Γ_MONTHLY_TOTAL_samples[i,:] = Γ_MONTHLY_TOTAL
            NPC_MONTHLY_TOTAL_samples[i,:] = NPC_MONTHLY_TOTAL
        end

        # %% Sample posterior access, give subnational NPC draws
        λ_ACCESS_samples = sample_net_access(ρ_chain_df, μ_chain_df, p_h,
                                    FULL_POPULATION_MONTHLY, Γ_MONTHLY_TOTAL_samples;
                                    n_max = 20)

        # %% Add to collection of outputs
        # Only keep mean for net sampled demography due to file size constraints
        COMBINED_A_BYNET_mean = mean(COMBINED_A_BYNET_samples, dims = 1)[1,:,:,:]
        COMBINED_A_TOTAL_mean = mean(COMBINED_A_TOTAL_samples,dims = 1)[1,:,:]
        admin1_output = Dict(   "area_id" => area_id,
                                "COMBINED_A_BYNET_mean" => COMBINED_A_BYNET_mean,
                                "Γ_MONTHLY_BYNET_samples" => Γ_MONTHLY_BYNET_samples,
                                "NPC_MONTHLY_BYNET_samples" => NPC_MONTHLY_BYNET_samples,
                                "COMBINED_A_TOTAL_mean" => COMBINED_A_TOTAL_mean,
                                "Γ_MONTHLY_TOTAL_samples" => Γ_MONTHLY_TOTAL_samples,
                                "NPC_MONTHLY_TOTAL_samples" => NPC_MONTHLY_TOTAL_samples,
                                "λ_ACCESS_samples" => λ_ACCESS_samples)

        # Store in variabl
        admin1_outputs[admin1_name_i] = copy(admin1_output)
        # push!(admin1_outputs, admin1_output)

    end

    ##############################################
    # COMPILE UNADJUSTED POSTERIOR NET CROP DRAWS TO A NICE ARRAY FOR LATER REGRESSION
    # %% Do draws for forward evolution
    ##############################################

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
    NAT_NPC_MONTHLY = missings(Float64, n_months)
    NAT_NPC_STD_MONTHLY = missings(Float64, n_months)

    for row_i in 1:size(filt_nat_npc_monthly)[1]
        month_val = filt_nat_npc_monthly[row_i,"month"]
        year_val = filt_nat_npc_monthly[row_i,"year"]
        monthidx = monthyear_to_monthidx(month_val, year_val, YEAR_START = YEAR_START_NAT)

        μ_est = filt_nat_npc_monthly[row_i,"NPC_mean"]*NAT_POPULATION_MONTHLY[monthidx]
        σ_est = filt_nat_npc_monthly[row_i,"NPC_adj_se"]*NAT_POPULATION_MONTHLY[monthidx]
        if (σ_est < μ_est) && σ_est>0
            NAT_NET_CROP_MONTHLY[monthidx] = μ_est
            NAT_NET_CROP_STD_MONTHLY[monthidx] = σ_est
            NAT_NPC_MONTHLY[monthidx] = filt_nat_npc_monthly[row_i,"NPC_mean"]
            NAT_NPC_STD_MONTHLY[monthidx] = filt_nat_npc_monthly[row_i,"NPC_adj_se"]
        end
    end


    ##############################
    # %% SUBNATIONAL DATA PREP
    ##############################
    # Define Arrays to extract required Net Crop values, and household surveys
    ALLREGIONS_Γ_MONTHLY_TOTAL_samples = zeros(n_admin1, n_samples, n_months)

    SUBNAT_NET_CROP_MONTHLY = missings(Float64, n_admin1, n_months)
    SUBNAT_NET_CROP_STD_MONTHLY = missings(Float64, n_admin1, n_months)
    SUBNAT_NPC_MONTHLY = missings(Float64, n_admin1, n_months)
    SUBNAT_NPC_STD_MONTHLY = missings(Float64, n_admin1, n_months)

    for admin1_name_i in 1:n_admin1
        # Look up posterior draws that were generated
        ALLREGIONS_Γ_MONTHLY_TOTAL_samples[admin1_name_i,:,:,:] .= sum(admin1_outputs[admin1_name_i]["Γ_MONTHLY_BYNET_samples"], dims = 3)[:,:,1]

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
                SUBNAT_NPC_MONTHLY[admin1_name_i,monthidx] = filt_subnat_npc_monthly[row_i,"NPC_mean"]
                SUBNAT_NPC_STD_MONTHLY[admin1_name_i,monthidx] = filt_subnat_npc_monthly[row_i,"NPC_adj_se"]
            end
        end
    end

    ALLREGIONS_Γ_MONTHLY_TOTAL_mean = mean(ALLREGIONS_Γ_MONTHLY_TOTAL_samples, dims = 2)[:,1,:]

    ##############################
    # %% MCMC to reconcile National and Subnational Estimates
    ##############################
    # %% Define RWMH Sampler
    rwmh_sampler = externalsampler(RWMH(MvNormal(zeros(n_admin1), I.*sample_var)))

    DISTRIBUTION_ANNUAL_TOTAL = sum(ALLREGIONS_DISTRIBUTION_ANNUAL_BYNET, dims = 3)[:,:,1]
    # %% MCMC Sampling for adjustment weights``
    model = nat_subnat_calibration(NAT_NET_CROP_MONTHLY,
                                        NAT_NET_CROP_STD_MONTHLY,
                                        SUBNAT_NET_CROP_MONTHLY,
                                        SUBNAT_NET_CROP_STD_MONTHLY,
                                        ALLREGIONS_Γ_MONTHLY_TOTAL_mean)
                                        
    output = sample(model, rwmh_sampler, iterations,
                        progress = true, discard_initial = burn_in)

    # %% Extract Admin1 Weights
    ω = mean(Matrix(DataFrame(output)[:,3:end-1]), dims = 1)[:]
    α = (ω./sum(ω)).*n_admin1 # Adjustment Weights

    ##############################
    # %% Update Posterior Net Crop Estimates by Adjustment Weights
    ##############################
    merged_outputs = Vector{Any}(undef, n_admin1)
    for admin1_name_i in 1:n_admin1
        println("Adjusted Access Sampling for $(ISO), region $(admin1_name_i) of $(n_admin1)...")
        admin1_output = admin1_outputs[admin1_name_i]

        # Get Area Id
        area_id = admin1_output["area_id"]

        # Calculate Adjusted Estimates
        ADJ_COMBINED_A_BYNET_mean = admin1_output["COMBINED_A_BYNET_mean"].*α[admin1_name_i]
        ADJ_Γ_MONTHLY_BYNET_samples = admin1_output["Γ_MONTHLY_BYNET_samples"].*α[admin1_name_i]
        ADJ_COMBINED_A_TOTAL_mean = sum(ADJ_COMBINED_A_BYNET_mean, dims = 3)[:,:,1]
        ADJ_Γ_MONTHLY_TOTAL_samples = sum(ADJ_Γ_MONTHLY_BYNET_samples, dims = 3)[:,:,1]
        ADJ_NPC_MONTHLY_BYNET_samples = copy(ADJ_Γ_MONTHLY_BYNET_samples)
        for i in 1:n_samples
            for net_type_i in 1:n_net_types
                ADJ_NPC_MONTHLY_BYNET_samples[i,:,net_type_i] = ADJ_NPC_MONTHLY_BYNET_samples[i,:,net_type_i]./SUBNAT_POPULATION_MONTHLY[admin1_name_i,:]
            end
        end
        ADJ_NPC_MONTHLY_TOTAL_samples = ADJ_Γ_MONTHLY_TOTAL_samples./repeat(SUBNAT_POPULATION_MONTHLY[admin1_name_i,:]', n_samples,1)


        # Re-sample Access

        FULL_POPULATION_MONTHLY = subnat_reg_data["admin1_outputs"][admin1_name_i]["FULL_POPULATION_MONTHLY"]
        ADJ_λ_ACCESS_samples = sample_net_access(ρ_chain_df, μ_chain_df, p_h,
                                                    FULL_POPULATION_MONTHLY, ADJ_Γ_MONTHLY_TOTAL_samples;
                                                n_max = 20)

        adj_admin1_output = Dict(   "area_id" => area_id,
                                    "ADJ_COMBINED_A_BYNET_mean" => ADJ_COMBINED_A_BYNET_mean,
                                    "ADJ_Γ_MONTHLY_BYNET_samples" => ADJ_Γ_MONTHLY_BYNET_samples,
                                    "ADJ_NPC_MONTHLY_BYNET_samples" => ADJ_NPC_MONTHLY_BYNET_samples,
                                    "ADJ_COMBINED_A_TOTAL_mean" => ADJ_COMBINED_A_TOTAL_mean,
                                    "ADJ_Γ_MONTHLY_TOTAL_samples" => ADJ_Γ_MONTHLY_TOTAL_samples,
                                    "ADJ_NPC_MONTHLY_TOTAL_samples" => ADJ_NPC_MONTHLY_TOTAL_samples,
                                    "ADJ_λ_ACCESS_samples" => ADJ_λ_ACCESS_samples,
                                    "NAT_NPC_MONTHLY" => NAT_NPC_MONTHLY,
                                    "NAT_NPC_STD_MONTHLY" => NAT_NPC_STD_MONTHLY,
                                    "SUBNAT_NPC_MONTHLY" => SUBNAT_NPC_MONTHLY,
                                    "SUBNAT_NPC_STD_MONTHLY" => SUBNAT_NPC_STD_MONTHLY
                                    )

        merged_output = merge(admin1_output, adj_admin1_output)
        
        merged_outputs[admin1_name_i] = copy(merged_output)
    end
    
    ##############################
    # %% Save outputs
    ##############################
    mkpath(save_dir)
    save_dir*"$(ISO)_SUBNAT_draws.jld2"
    jldsave(save_dir*"$(ISO)_SUBNAT_draws.jld2";
                    YEAR_START_NAT, YEAR_START, YEAR_END,
                    merged_outputs,
                    admin1_names,
                    α_weights = α)
end
