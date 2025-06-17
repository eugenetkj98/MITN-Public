"""
Author: Eugene Tan
Date Created: 6/5/2025
Last Updated: 6/5/2025
Make plots for convergence behaviour of EM algorithm (Nigeria)
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/read_toml.jl")

# %% Import packages
using JLD2
using CairoMakie
using LaTeXStrings
using DataFrames
using StatsBase
using ProgressBars

# %% Import Custom Packages
# using DataExtractions
# using DateConversions
using NetCropModel
# using NetCropRegression
# using NetAccessModel
# using NetAccessPrediction
# using NetAccessRegression

# %% Get Year Bounds
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

# %% Number of samples to draw
n_samples = 100

# %% Country to make plot for
ISO = "SEN"

####################################################
# %% Make prior/posterior predictions
####################################################

# %% Load Data
input_dict = load(OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
regression_dict = load(OUTPUT_REGRESSIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")

# %% Crop Regression input data
ISO = input_dict["ISO"]
YEARS_ANNUAL = input_dict["YEARS_ANNUAL"]
MONTHS_MONTHLY = input_dict["MONTHS_MONTHLY"]
DELIVERIES_ANNUAL = input_dict["DELIVERIES_ANNUAL"]
DISTRIBUTION_ANNUAL = input_dict["DISTRIBUTION_ANNUAL"]
NET_CROP_MONTHLY = input_dict["NET_CROP_MONTHLY"]
POPULATION_MONTHLY = input_dict["POPULATION_MONTHLY"]
NET_NAMES = input_dict["NET_NAMES"]
n_net_types = length(NET_NAMES)

# %% Get number of missing net values
n_missing_nets_vals = sum(ismissing.(DISTRIBUTION_ANNUAL[:,1]))

# # %% Crop Regression output data
# chain = regression_dict["chain"]
# monthly_p = regression_dict["monthly_p"]
# regression_dict["chain_epochs"]
# # %% Draw posterior samples of model parameters
# # Randomly generate sample indexes to sample from chain
# chain_length = size(chain)[1]
# sample_idxs = sample(1:chain_length, n_samples, replace = false)

# # Extract MCMC draws for parameters
# ϕ_posterior_draws = chain[sample_idxs, :ϕ]
# α_init_posterior_draws = chain[sample_idxs,:α_init]
# α_LLIN_posterior_draws = chain[sample_idxs,:α_LLIN]

# b_net_posterior_draws = Matrix(DataFrame(chain)[:,5:2:5+2*(n_net_types-1)])[sample_idxs, :]
# k_net_posterior_draws = Matrix(DataFrame(chain)[:,6:2:6+2*(n_net_types-1)])[sample_idxs, :]
# if n_missing_nets_vals > 0
#     n_missing_nets_posterior_draws = Matrix(DataFrame(chain)[:,(4+2*(n_net_types)+1):end])[sample_idxs, :]
# else
#     # No missing data, so just feed in a zero matrix with no columns
#     n_missing_nets_posterior_draws = zeros(size(chain)[1], 0)
# end

######################
# %% SCENARIO 1: Attrition parameter posteriors only
######################

# Extract chain and attrition parameters
chain = regression_dict["chain_epochs"][end]
monthly_p = regression_dict["monthly_weights"][end]

# Randomly sample idxs
chain_length = size(chain)[1]
sample_idxs = sample(1:chain_length, n_samples, replace = false)

# No missing init imputation, conversion or redist
# ϕ_posterior_draws = zeros(n_samples) #chain[sample_idxs, :ϕ]
# α_init_posterior_draws = zeros(n_samples) #chain[sample_idxs,:α_init]
# α_LLIN_posterior_draws = (ones(n_samples) .= 0.0001) #chain[sample_idxs,:α_LLIN]

ϕ_posterior_draws = chain[sample_idxs, :ϕ]
α_init_posterior_draws = chain[sample_idxs,:α_init]
α_LLIN_posterior_draws = chain[sample_idxs,:α_LLIN]



# No missing data imputation
if n_missing_nets_vals > 0
    n_missing_nets_posterior_draws = (zeros(n_samples,n_missing_nets_vals) .= 0)
else
    # No missing data, so just feed in a zero matrix with no columns
    n_missing_nets_posterior_draws = zeros(size(chain)[1], 0)
end

# Attrition Parameters
b_net_posterior_draws = Matrix(DataFrame(chain)[:,5:2:5+2*(n_net_types-1)])[sample_idxs, :]
k_net_posterior_draws = Matrix(DataFrame(chain)[:,6:2:6+2*(n_net_types-1)])[sample_idxs, :]

# Generate draws
Γ_MONTHLY_samples_BYNET_SCENARIO_1 = zeros(n_samples, length(MONTHS_MONTHLY), n_net_types)

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

                                                    

    Γ_MONTHLY_samples_BYNET_SCENARIO_1[i,:,:,:] = Γ_MONTHLY_BYNET
end

# Get totalsa and CI
Γ_MONTHLY_samples_TOTAL_SCENARIO_1 = sum(Γ_MONTHLY_samples_BYNET_SCENARIO_1, dims = 3)[:,:,1]

Γ_MONTHLY_CI_SCENARIO_1 = Matrix{Float64}(undef, size(Γ_MONTHLY_samples_TOTAL_SCENARIO_1)[2],3)
for i in 1:size(Γ_MONTHLY_CI_SCENARIO_1)[1]
    Γ_MONTHLY_CI_SCENARIO_1[i,[1,3]] .= quantile(Γ_MONTHLY_samples_TOTAL_SCENARIO_1[:,i], [0.025, 0.975])
    Γ_MONTHLY_CI_SCENARIO_1[i,2] = mean(Γ_MONTHLY_samples_TOTAL_SCENARIO_1[:,i])
end

# %% Test plot
fig = Figure(size = (800,300))
ax = Axis(fig[1,1])

band!(ax, 1:size(Γ_MONTHLY_CI_SCENARIO_1)[1], 
            Γ_MONTHLY_CI_SCENARIO_1[:,1],
            Γ_MONTHLY_CI_SCENARIO_1[:,3],
            alpha = 0.5)
lines!(ax, 1:size(Γ_MONTHLY_CI_SCENARIO_1)[1], 
                Γ_MONTHLY_CI_SCENARIO_1[:,2])

nonmissingidxs = findall(.!ismissing.(NET_CROP_MONTHLY))
scatter!(ax, (1:size(Γ_MONTHLY_CI_SCENARIO_1)[1])[nonmissingidxs], 
                Float64.(NET_CROP_MONTHLY[nonmissingidxs]))
fig