"""
Author: Eugene Tan
Date Created: 15/4/2025
Last Updated: 15/4/2025
Script to use existing aggregated LLIN SNF model to produce estimates for multitype case for all countries and save outputs
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Prep environment and subdirectories
include(pwd()*"/scripts/dir_configs.jl")

# %% Load packages
using JLD2
using CSV
using DataFrames

# %% Other useful packages
using Dates
using ProgressBars
using StatsBase

# %% Custom modules
using DateConversions
using NetCropModel
using NetLoss
using NetCropPrediction

# %% Get ISO List
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Define countries where cITNs are forced to have the same decay curve as LLINs due to lack of surveys pre 2010
excl_citn_ISOs = EXCLUSION_CITN_ISOS

# %% LOOP TO GENERATE PREDICTIONS
# for ISO in filt_ISOs
    ########################################################
    # %% Define required directories and prediction settings
    ########################################################
    input_timeseries_file = "datasets/forward_prediction_inputs/$(ISO)_prediction_input.csv"
    mitn_posterior_chain = "outputs/regressions/crop/2000_2023/$(ISO)_2000_2023_cropchains.jld2"
    n_samples = 100 # Number of posterior draws of time series

    output_dir = "outputs/predictions/"
    mkpath(output_dir)

    output_filename = "$(ISO)_netcrop_prediction.jld2" # NEED TO THINK OF A WAY TO AUTOMATE ADDING FILES

    ########################################################
    # %% Load Input Time Series and Scrape required meta-data
    ########################################################
    # Meta Data
    input_timeseries = CSV.read(input_timeseries_file, DataFrame)
    NET_NAMES = names(input_timeseries)[4:end]
    n_net_types = length(NET_NAMES)
    YEAR_START, YEAR_END = input_timeseries.year[[1,end]]
    YEARS_ANNUAL = YEAR_START:YEAR_END

    # Construct Delivery and Distribution Time Series
    DELIVERIES_ANNUAL = Array(input_timeseries.LLIN_delivery)
    DISTRIBUTION_ANNUAL = Array(input_timeseries[:,NET_NAMES])

    ########################################################
    # %% Load regression chain and construct lookup datafile for posterior draws
    ########################################################
    # Regressed YEAR RANGES
    reg_YEAR_START = YEAR_NAT_START
    reg_YEAR_END = YEAR_NAT_END

    # Load required base data
    mitn_model = JLD2.load(mitn_posterior_chain)
    reg_monthly_p = mitn_model["monthly_p"]
    reg_NET_NAMES = mitn_model["NET_NAMES"]
    reg_ϕ_chain = mitn_model["chain"].ϕ
    reg_α_init_chain = mitn_model["chain"].α_init
    reg_α_LLIN_chain = mitn_model["chain"].α_LLIN
    reg_τ_chain = rename!(mitn_model["chain"][:,5:2:5+2*(length(reg_NET_NAMES).-1)], reg_NET_NAMES)
    reg_κ_chain = rename!(mitn_model["chain"][:,6:2:6+2*(length(reg_NET_NAMES).-1)], reg_NET_NAMES)

    # Make cITN-LLIN attrition rate adjustment if country is one where cITNs are to be treated like LLINs due to lack of surveys pre-2010
    if ISO ∈ excl_citn_ISOs
        reg_τ_chain[:, "cITN"] .= reg_τ_chain[:, "LLIN"]
        reg_κ_chain[:, "cITN"] .= reg_κ_chain[:, "LLIN"]
    end

    # Extrapolate monthly_p to required annual ranges
    # Make storage variable
    n_months = (YEAR_END-YEAR_START+1)*12
    monthly_p = zeros(n_months)

    # Fill in known range from data
    reg_idx_start = monthyear_to_monthidx(1, reg_YEAR_START, YEAR_START = YEAR_START)
    reg_idx_end = monthyear_to_monthidx(12, reg_YEAR_END, YEAR_START = YEAR_START)
    monthly_p[reg_idx_start:reg_idx_end] = mitn_model["monthly_p"]

    # For unknown range, copy from the last 12 months of the regression which is based on the average annual trend
    annual_distribution_pattern = reg_monthly_p[end-11:end] # Get monthly pattern from last 12 months of regression
    fill_in_years = setdiff(YEAR_START:YEAR_END, reg_YEAR_START:reg_YEAR_END)
    for year in fill_in_years
        monthidx_start, monthidx_end = monthyear_to_monthidx.([1,12], year, YEAR_START = YEAR_START)
        monthly_p[monthidx_start:monthidx_end] = annual_distribution_pattern
    end

    # Construct posterior chains for parameter samples to align with desired NET_NAMES
    # If prediction net type is not found in existing list, then just pretend it behaves similarly to an LLIN

    τ_chain = DataFrame()
    κ_chain = DataFrame()

    for net_type in NET_NAMES
        if net_type ∈ reg_NET_NAMES
            τ_chain = hcat(τ_chain, reg_τ_chain[:, net_type], makeunique = true)
            κ_chain = hcat(κ_chain, reg_κ_chain[:, net_type], makeunique = true)
        else
            net_type_idx_LLIN = findfirst(reg_NET_NAMES .== "LLIN")
            τ_chain = hcat(τ_chain, reg_τ_chain[:, net_type_idx_LLIN], makeunique = true)
            κ_chain = hcat(κ_chain, reg_κ_chain[:, net_type_idx_LLIN], makeunique = true)
        end
    end

    rename!(τ_chain, NET_NAMES)
    rename!(κ_chain, NET_NAMES)

    ########################################################
    # %% Simulate SNF
    ########################################################
    # Storage variable for sampled posterior prediction
    Γ_BYNET_pred = zeros(n_samples, n_months, n_net_types)
    A_BYNET_pred = zeros(n_samples, n_months, n_months, n_net_types)

    Γ_BYNET_dist_pred = zeros(n_samples, n_months, n_net_types)
    A_BYNET_dist_pred = zeros(n_samples, n_months, n_months, n_net_types)

    # Select idx values to sample from
    sampled_idxs = rand(1:length(reg_ϕ_chain),n_samples)

    # Perform predictive sampling for each selected idx in MCMC chain 
    for i in ProgressBar(1:n_samples, leave = false)
        # Select sampled index
        idx = sampled_idxs[i]

        # Extract relevant chain values
        ϕ_est = reg_ϕ_chain[idx]
        α_LLIN_est = reg_α_LLIN_chain[idx]
        τ_net_est = τ_chain[idx,:]
        κ_net_est = κ_chain[idx,:]

        # Perform prediction
        Γ_MONTHLY_BYNET_sample, A_BYNET_sample = mitn_national_predict(YEARS_ANNUAL,
                                                            DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                                            ϕ_est, τ_net_est, κ_net_est, 
                                                            α_LLIN_est;
                                                            monthly_p = monthly_p)
                        
        Γ_MONTHLY_BYNET_dist_sample, A_BYNET_dist_sample = mitn_national_noredist_predict(YEARS_ANNUAL,
                                                                            DISTRIBUTION_ANNUAL,
                                                                            τ_net_est, κ_net_est;
                                                                            monthly_p = monthly_p)

        # Save prediction to storage variable
        Γ_BYNET_pred[i,:,:] = Γ_MONTHLY_BYNET_sample
        A_BYNET_pred[i,:,:,:] = A_BYNET_sample

        Γ_BYNET_dist_pred[i,:,:] = Γ_MONTHLY_BYNET_dist_sample
        A_BYNET_dist_pred[i,:,:,:] = A_BYNET_dist_sample
    end

    # %% Calculate 95% CI
    Γ_TOTAL_pred_sample = sum(Γ_BYNET_pred, dims = 3)[:,:,1]
    Γ_TOTAL_dist_pred_sample = sum(Γ_BYNET_dist_pred, dims = 3)[:,:,1]

    # Storage variable
    Γ_BYNET_pred_CI = zeros(size(Γ_TOTAL_pred_sample)[2],3, n_net_types)
    Γ_TOTAL_pred_CI = zeros(size(Γ_TOTAL_pred_sample)[2],3)

    Γ_BYNET_dist_pred_CI = zeros(size(Γ_TOTAL_dist_pred_sample)[2],3, n_net_types)
    Γ_TOTAL_dist_pred_CI = zeros(size(Γ_TOTAL_dist_pred_sample)[2],3)

    for i in 1:(size(Γ_TOTAL_pred_CI)[1])
        Γ_TOTAL_pred_CI[i,:] = quantile(Γ_TOTAL_pred_sample[:,i],[0.025, 0.5, 0.975])
        Γ_TOTAL_dist_pred_CI[i,:] = quantile(Γ_TOTAL_dist_pred_sample[:,i],[0.025, 0.5, 0.975])
        for j in 1:n_net_types
            Γ_BYNET_pred_CI[i,:,j] = quantile(Γ_BYNET_pred[:,i,j],[0.025, 0.5, 0.975])
            Γ_BYNET_dist_pred_CI[i,:,j] = quantile(Γ_BYNET_dist_pred[:,i,j],[0.025, 0.5, 0.975])
        end
    end

    # %%

    ########################################################
    # %% Save Data
    ########################################################

    output_dict = Dict(   "Γ_BYNET_pred_samples" => Γ_BYNET_pred,
            "A_BYNET_pred_samples" => A_BYNET_pred,
            "Γ_BYNET_pred" => Γ_BYNET_pred_CI,
            "Γ_TOTAL_pred" => Γ_TOTAL_pred_CI,
            "YEAR_VALS" => YEARS_ANNUAL,
            "NET_NAMES" => NET_NAMES)

    save(output_dir*output_filename, output_dict)
end

