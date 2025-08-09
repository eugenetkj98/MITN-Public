"""
Author: Eugene Tan
Date Created: 19/5/2024
Last Updated: 19/5/2024
AWS version of the subnational MITN regression code to be run using batch.
"""
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from TOML file
include(pwd()*"/scripts/read_toml.jl")

# %% Import Public Packages
using JLD2
using CSV
using DataFrames
using Missings
using ProgressBars
using LinearAlgebra
using StatsBase
using Distributions
using Turing
using AdvancedMH
using StatsPlots

# %% Import custom packages
using DateConversions
using NetLoss
using Subnat_NetCropModel

##############################################
# %% GLOBAL SETTINGS (NOT COUNTRY DEPENDENT)
##############################################
# %% Define paths
dataset_dir = RAW_SUBNAT_DATASET_DIR
dataprep_dir = OUTPUT_SUBNAT_DATAPREP_DIR
nat_netcrop_post_dir = OUTPUT_DRAWS_DIR*"national/crop_access/"
net_age_dir = OUTPUT_DRAWS_DIR*"national/demography/"
save_dir = OUTPUT_REGRESSIONS_DIR*"subnational/"
posterior_dir = OUTPUT_DIR

# %% MCMC and Regression Settings
n_chains = min(SUBNAT_CROP_ADJ_N_CHAINS, Threads.nthreads())

# Settings for calculating adjustment weights
adj_iterations = SUBNAT_CROP_ADJ_MCMC_ITERATIONS
adj_burn_in = SUBNAT_CROP_ADJ_MCMC_BURNIN

# Settings for fitting subnational attrition parameters
scaling_constant = SUBNAT_CROP_ATR_SCALING_CONSTANT
proposal_resampling_variances = SUBNAT_CROP_ATR_PROPOSAL_SAMPLING_VAR # Compact support sampling
subnat_iterations = SUBNAT_CROP_ATR_MCMC_ITERATIONS
subnat_burn_in = SUBNAT_CROP_ATR_MCMC_BURNIN

# %% Year bounds
YEAR_START_NAT = YEAR_NAT_START # Start year for national model
YEAR_START = YEAR_SUBNAT_TRANS # Start year for subnational model
YEAR_END = YEAR_NAT_END

# %% Perform draws and save outputs. Filter out unwanted countries
exclusion_ISOs = EXCLUSION_ISOS
ISO = ARGS[1]

# %%

if ISO ∈ exclusion_ISOs
    println("$(ISO) is on exclusion list. Skipping regression.")
    flush(stdout)
else
    ##############################################
    # %% COUNTRY SPECIFIC SETTINGS
    ##############################################
    println("Currently fitting subnational model for $(ISO)...")

    # %% Data Filenames
    # Population Data
    subnat_population_filename = SUBNAT_POPULATION_FILENAME
    # Net Distribution Data
    distributions_filename = SUBNAT_DISTRIBUTION_DATA_FILENAME
    # Admin 1 ID Legend
    id_legend_filename = ADMIN1_AREAID_LEGEND_FILENAME
    # Net age demography posterior
    net_age_filename = "$(ISO)_net_age_demography_mean.csv"
    # Net attrition posteriors (from national model)
    net_attrition_filename = "net_attrition_posteriors.csv"
    # NPC monthly dataset for regression
    subnat_npc_monthly_data_filename = HOUSEHOLD_SUBNAT_SUMMARY_DATA_FILENAME
    # National input dict used for the national net crop regression (used to get monthly population numbers)
    nat_netcrop_input_dict_filename = "$(ISO)_$(YEAR_START_NAT)_$(YEAR_END)_cropextract.jld2"
    # National net crop draws
    nat_netcrop_post_filename = "$(ISO)_$(YEAR_START_NAT)_$(YEAR_END)_post_crop_access.jld2"

    # %% National MCMC chain (for getting monthly disaggregation ratios)
    nat_cropchain_dir = OUTPUT_REGRESSIONS_DIR*"crop/$(YEAR_START_NAT)_$(YEAR_END)/"
    nat_cropchain_filename = "$(ISO)_$(YEAR_START_NAT)_$(YEAR_END)_cropchains.jld2"

    # %% File save settings
    output_filename = "$(ISO)_SUBNAT_NETCROP_$(YEAR_START_NAT)_$(YEAR_START)_$(YEAR_END)_regression.jld2"

    ##############################################
    # %% Open relevant datasets
    ##############################################

    # %% Import master datasets
    master_populations = CSV.read(dataset_dir*subnat_population_filename, DataFrame)
    master_distributions = CSV.read(dataset_dir*distributions_filename, DataFrame)
    master_id_legend = CSV.read(dataset_dir*id_legend_filename, DataFrame)
    master_net_age = CSV.read(net_age_dir*net_age_filename, DataFrame)
    master_net_attrition = CSV.read(posterior_dir*net_attrition_filename, DataFrame)

    # %% Load required country specific national level data
    # National input dict used for the national net crop regression (used to get monthly population numbers)
    nat_netcrop_input_dict = load(OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_START_NAT)_$(YEAR_END)/"*nat_netcrop_input_dict_filename)
    # Subnational filtered household survey data
    subnat_npc_monthly_data = CSV.read(dataprep_dir*subnat_npc_monthly_data_filename, DataFrame)
    # Load national net crop draws
    nat_netcrop_post_data = load(nat_netcrop_post_dir*nat_netcrop_post_filename)
    # Load National regression MCMC chain (for monthly disaggregation ratios)
    nat_cropchain = load(nat_cropchain_dir*nat_cropchain_filename)

    ##############################################
    # %% Extract relevant variables
    ##############################################
    # %% Get monthly disaggregation ratios and re-use for subnational regression
    full_monthly_p = nat_cropchain["monthly_p"] # Also extract the full set of monthly disaggregation ratios

    # %% National monthly populations
    FULL_NAT_POPULATION_MONTHLY = nat_netcrop_input_dict["POPULATION_MONTHLY"]

    ##############################################
    # %% Calculate relevant loop bounds (TIME INDEX)
    ##############################################
    # %% Calculate indexing bounds
    YEARS_ANNUAL = Vector(YEAR_START:1:(YEAR_END))
    FULL_YEARS_ANNUAL = Vector(YEAR_START_NAT:1:(YEAR_END))
    MONTHS_MONTHLY = Vector(1:(YEAR_END-YEAR_START+1)*12)
    FULL_MONTHS_MONTHLY = Vector(1:(YEAR_END-YEAR_START_NAT+1)*12)

    # Indexing calculations to filter full time period from national, to subnational time period
    idx_year_ref_1 = findfirst(FULL_YEARS_ANNUAL .== YEAR_START)
    idx_year_ref_2 = findfirst(FULL_YEARS_ANNUAL .== YEAR_END)

    idx_month_ref_1 = (idx_year_ref_1-1)*12+1
    idx_month_ref_2 = (idx_year_ref_2*12)
    subnat_monthidxs = idx_month_ref_1:idx_month_ref_2 

    # %% Get admin1 metadata
    country_legend = master_id_legend[master_id_legend.ISO.==ISO,:]
    admin1_names = country_legend.Name_1

    # %% Get number of Admin 1 Subnational regions
    n_admin1 = length(admin1_names)

    ##########################################
    # Get national model dynamics parameters
    ##########################################
    # Net Attrition Parameters
    net_attrition_params = master_net_attrition[master_net_attrition.ISO .== ISO,:]

    # Get Initial Condition Net demography
    initial_demography_condition = master_net_age[(master_net_age.ISO .== ISO) .&
                                                (master_net_age.YEAR .== YEAR_START) .&
                                                (master_net_age.MONTH .== 1),:]
    INIT_NET_NAMES = unique(initial_demography_condition.NET_TYPE)
    max_age_months = size(initial_demography_condition[:,6:end])[2]
    YEAR_HISTORY_START = YEAR_START-Int(max_age_months/12)

    ##############################################
    # %% Data Prep
    ##############################################
    ALLREGIONS_DISTRIBUTIONS_EXTRACT = []
    ALLREGIONS_POPULATION_ANNUAL_EXTRACT = []
    ALLREGIONS_FULL_POPULATION_MONTHLY_EXTRACT = []
    ALLREGIONS_NET_NAMES = []
    ALLREGIONS_NET_CROP_MONTHLY_EXTRACT = []
    ALLREGIONS_NET_CROP_STD_MONTHLY_EXTRACT = []

    
    # # Extract data for all regions
    for admin1_name_i in 1:n_admin1
        ##############################################
        # %% Select Region
        ##############################################
        # Get area_id
        admin1_name = admin1_names[admin1_name_i]
        area_id = country_legend[country_legend.Name_1 .== admin1_name,"area_id"][1]

        ##############################################
        # %% Population Data
        ##############################################
        admin1_master_populations = master_populations[(master_populations.area_id .== area_id) .&
                                                        (master_populations.year .>= YEAR_START) .&
                                                        (master_populations.year .<= YEAR_END+1),:]
        FULL_admin1_master_populations = master_populations[(master_populations.area_id .== area_id) .&
                                                        (master_populations.year .>= YEAR_START_NAT) .&
                                                        (master_populations.year .<= YEAR_END+1),:]

        POPULATION_ANNUAL = admin1_master_populations[sortperm(admin1_master_populations.year),"SUM"]
        FULL_POPULATION_ANNUAL = FULL_admin1_master_populations[sortperm(FULL_admin1_master_populations.year),"SUM"]                                    

        # POPULATION_MONTHLY = zeros(length(MONTHS_MONTHLY))
        FULL_POPULATION_MONTHLY = zeros(length(FULL_MONTHS_MONTHLY))

        # for i in 1:(length(YEARS_ANNUAL))
        #     POPULATION_MONTHLY[(12*(i-1)+1):(12*i)] = LinRange(POPULATION_ANNUAL[i], POPULATION_ANNUAL[i+1],13)[1:12]
        # end

        for i in 1:(length(FULL_YEARS_ANNUAL))
            FULL_POPULATION_MONTHLY[(12*(i-1)+1):(12*i)] = LinRange(FULL_POPULATION_ANNUAL[i], FULL_POPULATION_ANNUAL[i+1],13)[1:12]
        end
        push!(ALLREGIONS_POPULATION_ANNUAL_EXTRACT, POPULATION_ANNUAL)
        push!(ALLREGIONS_FULL_POPULATION_MONTHLY_EXTRACT, FULL_POPULATION_MONTHLY)

        ##########################################
        # Distribution Data
        ##########################################
        # %% Raw Distribution Data + Disaggregated from reported national reported distributions + AMP Tracker
        # Get list of net n_net_types
        filt_distributions = master_distributions[(master_distributions.iso .== ISO) .& 
                                        (master_distributions.area_id .== area_id) .& 
                                        (master_distributions.year .>= YEAR_START) .&
                                        (master_distributions.year .<= YEAR_END),:]

        FULL_filt_distributions = master_distributions[(master_distributions.iso .== ISO) .& 
                                        (master_distributions.area_id .== area_id) .& 
                                        (master_distributions.year .>= YEAR_START_NAT) .&
                                        (master_distributions.year .<= YEAR_END),:]

        present_net_ids = findall(sum(Array(FULL_filt_distributions[.!ismissing.(FULL_filt_distributions[:,"Total Nets"]), 10:end]), dims = 1)[:] .> 0)
        DIST_NET_NAMES = names(filt_distributions)[10:end][present_net_ids]

        # %% Get total list of net types
        NET_NAMES = union(DIST_NET_NAMES, INIT_NET_NAMES)
        n_net_types = length(NET_NAMES)

        # %% Annual distributions learned from national model
        NAT_model_est_DISTRIBUTION_ANNUAL_BYNET_mean = mean(nat_netcrop_post_data["DISTRIBUTION_ANNUAL_samples_BYNET"][:, idx_year_ref_1:idx_year_ref_2,:], dims = 1)[1,:,:]

        # %% Construct Matrix of net distributions
        SUBNAT_RAW_DISTRIBUTION_ANNUAL_BYNET = Matrix(filt_distributions[sortperm(filt_distributions.year),NET_NAMES])

        # %% Find all missing distribution data, and replace with population based disaggregation from national model estimates
        missing_dist_idx = findall(ismissing.(sum(SUBNAT_RAW_DISTRIBUTION_ANNUAL_BYNET, dims = 2)[:]))

        # Calculate population disaggregation ratios
        population_disagg_ratios = POPULATION_ANNUAL[1:end-1]./FULL_NAT_POPULATION_MONTHLY[subnat_monthidxs][1:12:end]

        # Estimates of net distribution by type for required years based on posterior national model estimates
        SUBNAT_model_est_DISTRIBUTION_ANNUAL_BYNET_mean = NAT_model_est_DISTRIBUTION_ANNUAL_BYNET_mean.*repeat(population_disagg_ratios,1, n_net_types)
        
        # Fill in missing data with posterior national estimates
        DISTRIBUTION_ANNUAL_BYNET = copy(SUBNAT_RAW_DISTRIBUTION_ANNUAL_BYNET)
        DISTRIBUTION_ANNUAL_BYNET[missing_dist_idx,:] = SUBNAT_model_est_DISTRIBUTION_ANNUAL_BYNET_mean[missing_dist_idx,:]

        push!(ALLREGIONS_DISTRIBUTIONS_EXTRACT, DISTRIBUTION_ANNUAL_BYNET)
        push!(ALLREGIONS_NET_NAMES, NET_NAMES)


        ##########################################
        # Household Survey Data Prep: Import subnat crop data for regression
        ##########################################

        # Extract observed netcrop estimates for survey aggregate (subnat_npc_monthly_data.csv)
        FULL_SUBNAT_NET_CROP_MONTHLY = missings(Float64,length(FULL_MONTHS_MONTHLY))
        FULL_SUBNAT_NET_CROP_STD_MONTHLY = missings(Float64,length(FULL_MONTHS_MONTHLY))
        
        subnat_npc_entries = subnat_npc_monthly_data[subnat_npc_monthly_data.area_id .== area_id,:]
        for row_i in 1:size(subnat_npc_entries)[1]
            month_val = subnat_npc_entries[row_i,"month"]
            year_val = subnat_npc_entries[row_i,"year"]
            monthidx = monthyear_to_monthidx(month_val, year_val, YEAR_START = YEAR_START_NAT)
            FULL_SUBNAT_NET_CROP_MONTHLY[monthidx] = subnat_npc_entries[row_i,"NPC_mean"]*FULL_POPULATION_MONTHLY[monthidx]
            FULL_SUBNAT_NET_CROP_STD_MONTHLY[monthidx] = subnat_npc_entries[row_i,"NPC_adj_se"]*FULL_POPULATION_MONTHLY[monthidx]
        end

        # IMPORTANT!! => For years where there are no surveys, use values from National regression
        # This is a simple way to include some delivery information, make models match a little better
        # without making the system too complicated

        # Find years where there are no surveys and need to be filled in from national estimates
        missing_years = setdiff(FULL_YEARS_ANNUAL,unique(subnat_npc_entries.year))

        # Convert missing years to month idx (specificall choose estimate at the start of the year)
        missing_monthidxs = monthyear_to_monthidx.(1,missing_years, YEAR_START = YEAR_START_NAT)

        # Get National posterior net_crop estimates
        nat_npc_samples = sum(nat_netcrop_post_data["Γ_MONTHLY_samples_BYNET"], dims = 3)[:,:,1]./repeat(FULL_NAT_POPULATION_MONTHLY',size(nat_netcrop_post_data["Γ_MONTHLY_samples_BYNET"])[1])
        nat_npc_mean = mean(nat_npc_samples, dims = 1)[1,:]
        nat_npc_std = std(nat_npc_samples, dims = 1)[1,:]
        
        # Fill data vector used for regression with national model estimates scaled by population for years with no survey data
        FULL_SUBNAT_NET_CROP_MONTHLY[missing_monthidxs] .= (FULL_POPULATION_MONTHLY.*nat_npc_mean)[missing_monthidxs]
        FULL_SUBNAT_NET_CROP_STD_MONTHLY[missing_monthidxs] .= (FULL_POPULATION_MONTHLY.*nat_npc_std)[missing_monthidxs]

        # Remove all cases where the STD is 0 as these will cause problems in the regression
        nonmissing_idxs = findall(.!ismissing.(FULL_SUBNAT_NET_CROP_STD_MONTHLY))
        nonmissing_zero_idxs = findall(FULL_SUBNAT_NET_CROP_STD_MONTHLY[nonmissing_idxs].==0)
        zero_std_idxs = nonmissing_idxs[nonmissing_zero_idxs]
        FULL_SUBNAT_NET_CROP_MONTHLY[zero_std_idxs] .= missing
        FULL_SUBNAT_NET_CROP_STD_MONTHLY[zero_std_idxs] .= missing

        push!(ALLREGIONS_NET_CROP_MONTHLY_EXTRACT, FULL_SUBNAT_NET_CROP_MONTHLY)
        push!(ALLREGIONS_NET_CROP_STD_MONTHLY_EXTRACT, FULL_SUBNAT_NET_CROP_STD_MONTHLY)
    end 

    # Compile extracted data into array
    n_subnat_years, n_net_types = size(ALLREGIONS_DISTRIBUTIONS_EXTRACT[1])[1:2]
    ALLREGIONS_FULL_POPULATION_MONTHLY_EXTRACT[1]
    ALLREGIONS_DISTRIBUTIONS_ARRAY = zeros(n_admin1, n_subnat_years, n_net_types)
    ALLREGIONS_FULL_POPULATION_ARRAY = zeros(n_admin1, size(ALLREGIONS_FULL_POPULATION_MONTHLY_EXTRACT[1])[1])
    ALLREGIONS_FULL_NET_CROP_ARRAY = missings(Float64, n_admin1, size(ALLREGIONS_FULL_POPULATION_MONTHLY_EXTRACT[1])[1])
    ALLREGIONS_FULL_NET_CROP_STD_ARRAY = missings(Float64, n_admin1, size(ALLREGIONS_FULL_POPULATION_MONTHLY_EXTRACT[1])[1])
    for i in 1:n_admin1
        ALLREGIONS_DISTRIBUTIONS_ARRAY[i,:,:] .= ALLREGIONS_DISTRIBUTIONS_EXTRACT[i]
        ALLREGIONS_FULL_POPULATION_ARRAY[i,:] .= ALLREGIONS_FULL_POPULATION_MONTHLY_EXTRACT[i]
        ALLREGIONS_FULL_NET_CROP_ARRAY[i,:] .= ALLREGIONS_NET_CROP_MONTHLY_EXTRACT[i]
        ALLREGIONS_FULL_NET_CROP_STD_ARRAY[i,:] .= ALLREGIONS_NET_CROP_STD_MONTHLY_EXTRACT[i]
    end
    
    ##############################################
    # %% Subnational Regression
    ##############################################

    ##############################################
    # %% Step 1: Get Distribution Adjustment Weights by subnational based on household survey. 
    # Assume same attrition parameters as national model
    ##############################################

    # Trim monthly_p to desired subnational regression period
    monthly_p = full_monthly_p[subnat_monthidxs]

    # Get required household net crop data
    ALLREGIONS_NET_CROP_MONTHLY = ALLREGIONS_FULL_NET_CROP_ARRAY[:,subnat_monthidxs]
    ALLREGIONS_NET_CROP_STD_MONTHLY = ALLREGIONS_FULL_NET_CROP_STD_ARRAY[:,subnat_monthidxs]

    # Calculate total national net crop to optimise MCMC regression algorithm
    NAT_DISTRIBUTIONS_TOTAL = sum(ALLREGIONS_DISTRIBUTIONS_ARRAY, dims = (1,3))[:]
    
    # %% Do MCMC sampling to get posterior adjustment weights
    # Define proposal sampler distribution
    adj_weights_rwmh = externalsampler(RWMH(MvNormal(zeros(n_admin1), I.*0.0002)))
    
    # Sample MCMC
    output = []
    if n_chains > 1
        output = sample(subnat_weighted_monthly_itn(YEARS_ANNUAL, MONTHS_MONTHLY,
                                            ALLREGIONS_POPULATION_ANNUAL_EXTRACT,
                                            ALLREGIONS_DISTRIBUTIONS_ARRAY,
                                            NAT_DISTRIBUTIONS_TOTAL,
                                            ALLREGIONS_NET_NAMES,
                                            net_attrition_params,
                                            initial_demography_condition,
                                            ALLREGIONS_NET_CROP_MONTHLY,
                                            ALLREGIONS_NET_CROP_STD_MONTHLY; 
                                            monthly_p = monthly_p,
                                            ϵ = 0.003), 
                                                adj_weights_rwmh, MCMCThreads(), adj_iterations, n_chains,
                                                progress = true, discard_initial = adj_burn_in)
    else
        output = sample(subnat_weighted_monthly_itn(YEARS_ANNUAL, MONTHS_MONTHLY,
                                            ALLREGIONS_POPULATION_ANNUAL_EXTRACT,
                                            ALLREGIONS_DISTRIBUTIONS_ARRAY,
                                            NAT_DISTRIBUTIONS_TOTAL,
                                            ALLREGIONS_NET_NAMES,
                                            net_attrition_params,
                                            initial_demography_condition,
                                            ALLREGIONS_NET_CROP_MONTHLY,
                                            ALLREGIONS_NET_CROP_STD_MONTHLY; 
                                            monthly_p = monthly_p,
                                            ϵ = 0.003), 
                                                adj_weights_rwmh, adj_iterations, 
                                                progress = true, discard_initial = adj_burn_in)
    end
    
    ##############################################
    # %% Step 2: Rescale distribution values by regressed weights for each subnational admin region
    ##############################################

    # Calculate Weight Adjusted Net Distributions for Subnational Regression
    ω = mean(Matrix(DataFrame(output)[:,3:3+n_admin1-1]), dims = 1)[:]
    
    ADJ_ALLREGIONS_DISTRIBUTIONS_ARRAY = zeros(size(ALLREGIONS_DISTRIBUTIONS_ARRAY))
    for i in 1:n_admin1
        ADJ_ALLREGIONS_DISTRIBUTIONS_ARRAY[i,:,:] = ALLREGIONS_DISTRIBUTIONS_ARRAY[i,:,:].*ω[i]
    end
 

    ##############################################
    # %% Step 3: Perform subnational level regression of attrition parameters
    ##############################################

    # %% Storage variable for results
    admin1_outputs = Vector{Any}(undef, n_admin1)

    # %% Analyse all admin regions and save in required list variables
    Threads.@threads for admin1_name_i in 1:n_admin1

        println("Analysing subnational region for $(ISO). $(admin1_name_i) of $(length(admin1_names))...")
        
        ##########################################
        # Prepare required data
        ##########################################

        # Get area_id
        admin1_name = admin1_names[admin1_name_i]
        area_id = country_legend[country_legend.Name_1 .== admin1_name,"area_id"][1]

        # Get required annual population data for subnational region
        POPULATION_ANNUAL = ALLREGIONS_POPULATION_ANNUAL_EXTRACT[admin1_name_i]

        # Get required list of net names
        NET_NAMES = ALLREGIONS_NET_NAMES[admin1_name_i]

        # Get household net crop values to regres against
        SUBNAT_NET_CROP_MONTHLY = ALLREGIONS_NET_CROP_MONTHLY[admin1_name_i,:]
        SUBNAT_NET_CROP_STD_MONTHLY = ALLREGIONS_NET_CROP_STD_MONTHLY[admin1_name_i,:]

        # Get annual distribution numbers using adjusted values from pervious section
        DISTRIBUTION_ANNUAL_BYNET = ADJ_ALLREGIONS_DISTRIBUTIONS_ARRAY[admin1_name_i,:,:]

        # Get population for entrie year period
        FULL_POPULATION_MONTHLY = ALLREGIONS_FULL_POPULATION_ARRAY[admin1_name_i,:]

        ##########################################
        # Try running MCMC Turing regression
        ##########################################
        
        # Calculate number of r.v. added to impute missing data
        n_missing = sum(ismissing.(DISTRIBUTION_ANNUAL_BYNET[:,1]))*n_net_types

        # Construct final resampling variance vector and covariance matrix
        resampling_variances = vcat(repeat(proposal_resampling_variances[1:2], n_net_types),
                                        ones(n_missing).*(proposal_resampling_variances[3]/scaling_constant))
        subnat_rwmh = externalsampler(RWMH(MvNormal(zeros(length(resampling_variances)), Diagonal(resampling_variances))))

        # Define Random Walk Metropolis Hastings external sampler. Taken from AdvancedMH package.
        output = sample(SUBNAT_monthly_itn(YEARS_ANNUAL, MONTHS_MONTHLY,
                                    POPULATION_ANNUAL,
                                    DISTRIBUTION_ANNUAL_BYNET,
                                    SUBNAT_NET_CROP_MONTHLY,
                                    SUBNAT_NET_CROP_STD_MONTHLY,
                                    NET_NAMES,
                                    initial_demography_condition,
                                    net_attrition_params,
                                    monthly_p = monthly_p), 
                                    subnat_rwmh, subnat_iterations,
                                    progress = true, discard_initial = subnat_burn_in)

        # %% Get mean posterior estimates
        τ_chain = DataFrame(output)[:,3:2:(2+n_net_types*2)]
        κ_chain = DataFrame(output)[:,4:2:(2+n_net_types*2)]
        missing_nets_chain = DataFrame(output)[:,(2+n_net_types*2+1):(2+n_net_types*2+n_missing)]

        τ_est = mean(Matrix(τ_chain), dims = 1)[:]
        κ_est = mean(Matrix(κ_chain), dims = 1)[:]

        missing_nets_est = []

        if n_missing > 0
            missing_nets_est = mean(Matrix(DataFrame(output)[:,(2+n_net_types*2+1):(2+n_net_types*2+n_missing)]), dims = 1)[:]
        end 

        # Define dict of results to save
        output_dict = Dict("area_id" => area_id, "admin1_name" => admin1_name,
                            "NET_NAMES" => NET_NAMES,
                            "FULL_POPULATION_MONTHLY" => FULL_POPULATION_MONTHLY,
                            "RAW_DISTRIBUTION_ANNUAL_BYNET" => ALLREGIONS_DISTRIBUTIONS_ARRAY,
                            "ADJ_DISTRIBUTION_ANNUAL_BYNET" => DISTRIBUTION_ANNUAL_BYNET,
                            "τ_chain" => τ_chain, "κ_chain" => κ_chain, "missing_nets_chain" => missing_nets_chain,
                            "τ_est" => τ_est, "κ_est" => κ_est, "missing_nets_est" => missing_nets_est,
                            "monthly_p" => monthly_p)

        # Store in variable
        admin1_outputs[admin1_name_i] = copy(output_dict)
    end

    ##############################################
    # %% Step 4: Save Outputs
    ##############################################
    mkpath(save_dir)
    jldsave(save_dir*output_filename; ISO,
                        YEAR_START_NAT,
                        YEAR_START,
                        YEAR_END,
                        admin1_names,
                        admin1_outputs,
                        adj_weights = ω)
end

