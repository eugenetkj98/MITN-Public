"""
Author: Eugene Tan
Date Created: 4/11/2024
Last Updated: 4/11/2024
Extract posterior net crop and attrition estimates from regressed national crop model to interface with subnational model.
"""
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Import Public Packages
using JLD2
using CSV
using DataFrames
using DateConversions
using ProgressBars
using LinearAlgebra
using StatsBase
using KernelDensity
using Plots


# %% Extraction Settings
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END


# %% Define Posterior Models data location
netcrop_reg_dir = OUTPUT_REGRESSIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/"

"""
    kernelmode_gaussian(samples::Vector{Float64})
Calculates the kernel density estimate of the sample draws. Returns the mode and standard deviation based on the KDE.
"""

function kernelmode_gaussian(samples)
    # Calculate KDE of samples
    kernel = kde(samples)

    # Get Gaussian parameters with mean = mode of the KDE, std = std(sample-mode)
    μ = kernel.x[argmax(kernel.density)]
    σ = sqrt(sum((samples.-μ).^2)/(length(samples)-1))
    
    return μ, σ
end

# %%
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = ["CPV","ZAF"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Extract net attrition parameters for all countries
net_attrition_summary = DataFrame()

for ISO in filt_ISOs
    # Get Net Attrition Posteriors MCMC chain
    post_crop_filename = ISO*"_$(YEAR_START)_$(YEAR_END)_cropchains.jld2"
    crop_chain = load(netcrop_reg_dir*post_crop_filename)

    # Get net profiles
    NET_NAMES = crop_chain["NET_NAMES"]
    n_net_types = length(NET_NAMES)
    τ_net_chain = crop_chain["chain"][:,5:2:5+(n_net_types-1)*2]
    κ_net_chain = crop_chain["chain"][:,6:2:6+(n_net_types-1)*2]

    # Calculate Gaussian parameters based on the mode and std (centered around the mode)
    net_summary_values = zeros(n_net_types, 6)

    for net_type_i in 1:n_net_types
        net_summary_values[net_type_i,1] = mean(τ_net_chain[:,net_type_i])
        net_summary_values[net_type_i,2] = mean(κ_net_chain[:,net_type_i])
        net_summary_values[net_type_i,3:4] .= kernelmode_gaussian(τ_net_chain[:,net_type_i])
        net_summary_values[net_type_i,5:6] .= kernelmode_gaussian(κ_net_chain[:,net_type_i])
    end

    # Construct DataFrame of processed posterior data
    country_net_attrition_summary = DataFrame(ISO = ISO, NET_TYPE = NET_NAMES, 
                                                mean_tau = net_summary_values[:,1],
                                                mean_kappa = net_summary_values[:,2],
                                                kde_gaussian_tau_mu = net_summary_values[:,3],
                                                kde_gaussian_tau_sigma = net_summary_values[:,4],
                                                kde_gaussian_kappa_mu = net_summary_values[:,5],
                                                kde_gaussian_kappa_sigma = net_summary_values[:,6])
    
    # Update combined data frame
    net_attrition_summary = vcat(net_attrition_summary, country_net_attrition_summary)
end

# Save summary as CSV
CSV.write(OUTPUT_DIR*"net_attrition_posteriors.csv",net_attrition_summary)
