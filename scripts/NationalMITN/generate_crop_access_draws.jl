"""
Author: Eugene Tan
Date Created: 23/7/2024
Last Updated: 26/8/2024
Code to generate the mean net crop posterior estimate. Saves as a CSV file
"""


# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import Public Packages
using JLD2
using CSV
using DataFrames

# %% Import Custom Packages
using DataExtractions
using DateConversions
using NetCropModel
using NetCropRegression
using NetAccessModel
using NetAccessRegression

# Maths packages
using LinearAlgebra
using StatsBase
using ProgressBars

# %% Output directory
output_dir = "outputs/draws/national/crop_access/"

# %% Sampling Settings
n_samples = 50#1000

# %% Get ISO List
ISO_list = String.(CSV.read("datasets/ISO_list.csv", DataFrame)[:,1])
exclusion_ISOs = []#["CPV","BWA","GNQ","DJI","ETH","SOM","ZAF","SSD"]
# ["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD"]
# GAB excluded
# %% Run Analysis
YEAR_START = 2000
YEAR_END = 2023 

# %% Perform draws and save outputs
# ISO_list = String.(CSV.read(raw"C:\Users\ETan\Documents\Prototype Analyses\itn-updated\datasets\ISO_list.csv", DataFrame)[:,1])
# exclusion_ISOs = ["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

for ISO in filt_ISOs
    println("Generating draws of national crop and access for $(ISO)...")
    # Import Data
    input_dict = load("outputs/extractions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
    # regression_dict = load("outputs/regressions/crop/Compact Regressions/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")
    regression_dict = load("outputs/regressions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")
    net_access_input_dict = load("outputs/extractions/access/pred_data/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")
    net_access_chain = load("outputs/regressions/access/netaccesschains.jld2")

    # Crop Regression input data
    ISO = input_dict["ISO"]
    YEARS_ANNUAL = input_dict["YEARS_ANNUAL"]
    MONTHS_MONTHLY = input_dict["MONTHS_MONTHLY"]
    DELIVERIES_ANNUAL = input_dict["DELIVERIES_ANNUAL"]
    DISTRIBUTION_ANNUAL = input_dict["DISTRIBUTION_ANNUAL"]
    NET_CROP_MONTHLY = input_dict["NET_CROP_MONTHLY"]
    POPULATION_MONTHLY = input_dict["POPULATION_MONTHLY"]
    NET_NAMES = input_dict["NET_NAMES"]
    n_net_types = length(NET_NAMES)

    # Get number of missing net values
    n_missing_nets_vals = sum(ismissing.(DISTRIBUTION_ANNUAL[:,1]))

    # Crop Regression output data
    chain = regression_dict["chain"]
    monthly_p = regression_dict["monthly_p"]

    # Access Regression output data
    ρ_chain_df = net_access_chain["ρ_chain_df"]
    μ_chain_df = net_access_chain["μ_chain_df"]
    p_h = mean(net_access_input_dict["p_h_aggregated"], dims = 1)[:]
    
    ### Draw posterior samples of model parameters
    # Randomly generate sample indexes to sample from chain
    chain_length = size(chain)[1]
    sample_idxs = sample(1:chain_length, n_samples, replace = false)

    # Extract MCMC draws for parameters
    ϕ_posterior_draws = chain[sample_idxs, :ϕ]
    α_init_posterior_draws = chain[sample_idxs,:α_init]
    α_LLIN_posterior_draws = chain[sample_idxs,:α_LLIN]

    b_net_posterior_draws = Matrix(DataFrame(chain)[:,5:2:5+2*(n_net_types-1)])[sample_idxs, :]
    k_net_posterior_draws = Matrix(DataFrame(chain)[:,6:2:6+2*(n_net_types-1)])[sample_idxs, :]
    if n_missing_nets_vals > 0
        n_missing_nets_posterior_draws = Matrix(DataFrame(chain)[:,(4+2*(n_net_types)+1):end])[sample_idxs, :]
    else
        # No missing data, so just feed in a zero matrix with no columns
        n_missing_nets_posterior_draws = zeros(size(chain)[1], 0)
    end


    ##### Generate Net Crop Trajectories
    Γ_MONTHLY_samples_BYNET = zeros(n_samples, length(MONTHS_MONTHLY), n_net_types)
    A_samples_BYNET = zeros(n_samples, length(MONTHS_MONTHLY), length(MONTHS_MONTHLY), n_net_types)
    DISTRIBUTION_ANNUAL_samples_BYNET = zeros(n_samples, length(YEARS_ANNUAL), n_net_types)

    Threads.@threads for i in ProgressBar(1:n_samples, leave = false)
        # Select MCMC posterior draw parameters
        ϕ = ϕ_posterior_draws[i]
        α_init = α_init_posterior_draws[i]
        α_LLIN = α_LLIN_posterior_draws[i]
        b_nets = b_net_posterior_draws[i,:]
        k_nets = k_net_posterior_draws[i,:]
        missing_nets = n_missing_nets_posterior_draws[i,:]

        
        
        Γ_MONTHLY_BYNET, A_BYNET, 
        COUNTRY_LLIN_STOCK_ANNUAL_EOY, 
        UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET, 
        ADJUSTED_DISTRIBUTION_ANNUAL_BYNET = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                                                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                                                                    ϕ, b_nets, k_nets,
                                                                                    α_init, α_LLIN,
                                                                                    missing_nets; 
                                                                                    monthly_p = monthly_p,
                                                                                    return_age = true, return_stock = true)

        Γ_MONTHLY_samples_BYNET[i,:,:,:] = Γ_MONTHLY_BYNET
        A_samples_BYNET[i,:,:,:] = A_BYNET
        DISTRIBUTION_ANNUAL_samples_BYNET[i,:,:] = ADJUSTED_DISTRIBUTION_ANNUAL_BYNET
    end

    # Get totals
    Γ_MONTHLY_samples_TOTAL = sum(Γ_MONTHLY_samples_BYNET, dims = 3)[:,:,1]

    # Get calculate per capita net demography matrix
    A_NPC_samples_BYNET = zeros(size(A_samples_BYNET))
    for i in ProgressBar(1:n_samples, leave = false)
        for month_i in 1:length(MONTHS_MONTHLY)
            A_NPC_samples_BYNET[i,month_i,:,:] = A_samples_BYNET[i,month_i,:,:]./POPULATION_MONTHLY[month_i]
        end
    end

    # Calculate means
    Γ_MONTHLY_mean_TOTAL = mean(Γ_MONTHLY_samples_TOTAL, dims = 1)[1,:]
    NPC_MONTHLY_mean_TOTAL = Γ_MONTHLY_mean_TOTAL./POPULATION_MONTHLY
    
    A_mean_BYNET = mean(A_samples_BYNET, dims = 1)[1,:,:,:]
    A_NPC_mean_BYNET = zeros(size(A_mean_BYNET))
    for month_i in 1:length(MONTHS_MONTHLY)
        A_NPC_mean_BYNET[month_i,:,:] = A_mean_BYNET[month_i,:,:]./POPULATION_MONTHLY[month_i]
    end

    ##### Generate Net Access Trajectories
    λ_access_samples = sample_net_access(ρ_chain_df, μ_chain_df, p_h,
                                        POPULATION_MONTHLY, Γ_MONTHLY_samples_TOTAL) # Need to temporarily unscale input by mil for calculating access

    λ_access_mean = mean(λ_access_samples, dims = 1)[:]

    #### Save posterior estimate as a JLD2 file
    filename = "$(ISO)_$(YEAR_START)_$(YEAR_END)_post_crop_access.jld2"
    mkpath(output_dir)
    jldsave(output_dir*filename; 
            YEAR_START = YEAR_START,
            YEAR_END = YEAR_END, 
            DISTRIBUTION_ANNUAL_samples_BYNET= DISTRIBUTION_ANNUAL_samples_BYNET,
            Γ_MONTHLY_samples_BYNET = Γ_MONTHLY_samples_BYNET,
            Γ_MONTHLY_mean_TOTAL = Γ_MONTHLY_mean_TOTAL, 
            A_mean_BYNET = A_mean_BYNET,
            A_NPC_mean_BYNET = A_NPC_mean_BYNET, 
            λ_access_samples = λ_access_samples,
            λ_access_mean = λ_access_mean,
            NET_NAMES = NET_NAMES,
            POPULATION_MONTHLY = POPULATION_MONTHLY)
end
