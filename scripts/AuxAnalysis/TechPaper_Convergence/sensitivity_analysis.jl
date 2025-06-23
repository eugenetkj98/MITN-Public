"""
Author: Eugene Tan
Date Created: 17/6/2025
Last Updated: 17/6/2025
Sensitivity analysis on amount of data required for convergence (i.e. sensitivity analysis)
Tests three different types of data omission
- Random Removal
- Chronological
- Reverse Chronological
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/read_toml.jl")

# %% Import Public Packages
using JLD2
using CSV
using DataFrames

# %% Import Custom Packages
using DataExtractions
using DateConversions
using NetCropModel
using NetCropRegression

# Maths packages
using LinearAlgebra
using StatsBase

# %% Get ISO to analyse from argument input
ISO = "NGA"
exclusion_ISOs = EXCLUSION_ISOS

# %% Run Analysis
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

# %% 
output_dir =  OUTPUT_REGRESSIONS_DIR*"crop/convergence_test/sensitivity_analysis/"

# %% Regression
# Load extracted data
input_dict = JLD2.load(OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")

nonmissingidx_llin = findall(.!ismissing.(input_dict["LLIN_CROP_MONTHLY"]))
nonmissingidx_nets = findall(.!ismissing.(input_dict["NET_CROP_MONTHLY"]))

# Define number of data to preserve for regression

# n_data = min(parse(Int, ARGS[1]), length(nonmissingidx_llin))
n_data = min(10, length(nonmissingidx_llin))
println("Number of data points in regression: $(n_data)/$(length(nonmissingidx_llin))...")


# Get alternative input_dicts
## Chronological
CHRONOLOGICAL_LLIN_CROP_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["LLIN_CROP_MONTHLY"])) .= missing)
CHRONOLOGICAL_NET_CROP_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["NET_CROP_MONTHLY"])) .= missing)
CHRONOLOGICAL_NET_CROP_STD_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["NET_CROP_STD_MONTHLY"])) .= missing)
CHRONOLOGICAL_cITN_CROP_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["cITN_CROP_MONTHLY"])) .= missing)
CHRONOLOGICAL_HOUSEHOLD_NPC_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["HOUSEHOLD_NPC_MONTHLY"])) .= missing)
CHRONOLOGICAL_HOUSEHOLD_NPC_STD_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["HOUSEHOLD_NPC_STD_MONTHLY"])) .= missing)

CHRONOLOGICAL_LLIN_CROP_MONTHLY[nonmissingidx_llin[1:n_data]] .= input_dict["LLIN_CROP_MONTHLY"][nonmissingidx_llin[1:n_data]]
CHRONOLOGICAL_NET_CROP_MONTHLY[nonmissingidx_nets[1:n_data]] .= input_dict["NET_CROP_MONTHLY"][nonmissingidx_nets[1:n_data]]
CHRONOLOGICAL_NET_CROP_STD_MONTHLY[nonmissingidx_nets[1:n_data]] .= input_dict["NET_CROP_STD_MONTHLY"][nonmissingidx_nets[1:n_data]]
CHRONOLOGICAL_cITN_CROP_MONTHLY[nonmissingidx_nets[1:n_data]] .= input_dict["cITN_CROP_MONTHLY"][nonmissingidx_nets[1:n_data]]
CHRONOLOGICAL_HOUSEHOLD_NPC_MONTHLY[nonmissingidx_nets[1:n_data]] .= input_dict["HOUSEHOLD_NPC_MONTHLY"][nonmissingidx_nets[1:n_data]]
CHRONOLOGICAL_HOUSEHOLD_NPC_STD_MONTHLY[nonmissingidx_nets[1:n_data]] .= input_dict["HOUSEHOLD_NPC_STD_MONTHLY"][nonmissingidx_nets[1:n_data]]


## Reverse Chronological
REV_CHRONOLOGICAL_LLIN_CROP_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["LLIN_CROP_MONTHLY"])) .= missing)
REV_CHRONOLOGICAL_NET_CROP_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["NET_CROP_MONTHLY"])) .= missing)
REV_CHRONOLOGICAL_NET_CROP_STD_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["NET_CROP_STD_MONTHLY"])) .= missing)
REV_CHRONOLOGICAL_cITN_CROP_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["cITN_CROP_MONTHLY"])) .= missing)
REV_CHRONOLOGICAL_HOUSEHOLD_NPC_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["HOUSEHOLD_NPC_MONTHLY"])) .= missing)
REV_CHRONOLOGICAL_HOUSEHOLD_NPC_STD_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["HOUSEHOLD_NPC_STD_MONTHLY"])) .= missing)

REV_CHRONOLOGICAL_LLIN_CROP_MONTHLY[nonmissingidx_llin[(end-n_data+1):end]] .= input_dict["LLIN_CROP_MONTHLY"][nonmissingidx_llin[(end-n_data+1):end]]
REV_CHRONOLOGICAL_NET_CROP_MONTHLY[nonmissingidx_llin[(end-n_data+1):end]] .= input_dict["NET_CROP_MONTHLY"][nonmissingidx_llin[(end-n_data+1):end]]
REV_CHRONOLOGICAL_NET_CROP_STD_MONTHLY[nonmissingidx_llin[(end-n_data+1):end]] .= input_dict["NET_CROP_STD_MONTHLY"][nonmissingidx_llin[(end-n_data+1):end]]
REV_CHRONOLOGICAL_cITN_CROP_MONTHLY[nonmissingidx_llin[(end-n_data+1):end]] .= input_dict["cITN_CROP_MONTHLY"][nonmissingidx_llin[(end-n_data+1):end]]
REV_CHRONOLOGICAL_HOUSEHOLD_NPC_MONTHLY[nonmissingidx_llin[(end-n_data+1):end]] .= input_dict["HOUSEHOLD_NPC_MONTHLY"][nonmissingidx_llin[(end-n_data+1):end]]
REV_CHRONOLOGICAL_HOUSEHOLD_NPC_STD_MONTHLY[nonmissingidx_llin[(end-n_data+1):end]] .= input_dict["HOUSEHOLD_NPC_STD_MONTHLY"][nonmissingidx_llin[(end-n_data+1):end]]

## Random Removal
n_samples = 10

RANDOM_LLIN_CROP_MONTHLY_SAMPLES = Vector{Any}(undef, n_samples)
RANDOM_NET_CROP_MONTHLY_SAMPLES = Vector{Any}(undef, n_samples)
RANDOM_NET_CROP_STD_MONTHLY_SAMPLES = Vector{Any}(undef, n_samples)
RANDOM_cITN_CROP_MONTHLY_SAMPLES = Vector{Any}(undef, n_samples)
RANDOM_HOUSEHOLD_NPC_MONTHLY_SAMPLES = Vector{Any}(undef, n_samples)
RANDOM_HOUSEHOLD_NPC_STD_MONTHLY_SAMPLES = Vector{Any}(undef, n_samples)

for i in 1:n_samples
    RANDOM_LLIN_CROP_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["LLIN_CROP_MONTHLY"])) .= missing)
    RANDOM_NET_CROP_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["NET_CROP_MONTHLY"])) .= missing)
    RANDOM_NET_CROP_STD_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["NET_CROP_STD_MONTHLY"])) .= missing)
    RANDOM_cITN_CROP_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["cITN_CROP_MONTHLY"])) .= missing)
    RANDOM_HOUSEHOLD_NPC_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["HOUSEHOLD_NPC_MONTHLY"])) .= missing)
    RANDOM_HOUSEHOLD_NPC_STD_MONTHLY = (Vector{Union{Missing, Float64}}(undef, length(input_dict["HOUSEHOLD_NPC_STD_MONTHLY"])) .= missing)

    sampled_idx = sample(nonmissingidx_llin, n_data, replace = false)

    RANDOM_LLIN_CROP_MONTHLY[sampled_idx] .= input_dict["LLIN_CROP_MONTHLY"][sampled_idx]
    RANDOM_NET_CROP_MONTHLY[sampled_idx] .= input_dict["NET_CROP_MONTHLY"][sampled_idx]
    RANDOM_NET_CROP_STD_MONTHLY[sampled_idx] .= input_dict["NET_CROP_STD_MONTHLY"][sampled_idx]
    RANDOM_cITN_CROP_MONTHLY[sampled_idx] .= input_dict["cITN_CROP_MONTHLY"][sampled_idx]
    RANDOM_HOUSEHOLD_NPC_MONTHLY[sampled_idx] .= input_dict["HOUSEHOLD_NPC_MONTHLY"][sampled_idx]
    RANDOM_HOUSEHOLD_NPC_STD_MONTHLY[sampled_idx] .= input_dict["HOUSEHOLD_NPC_STD_MONTHLY"][sampled_idx]
    

    RANDOM_LLIN_CROP_MONTHLY_SAMPLES[i] = deepcopy(RANDOM_LLIN_CROP_MONTHLY)
    RANDOM_NET_CROP_MONTHLY_SAMPLES[i] = deepcopy(RANDOM_NET_CROP_MONTHLY)
    RANDOM_NET_CROP_STD_MONTHLY_SAMPLES[i] = deepcopy(RANDOM_NET_CROP_STD_MONTHLY)
    RANDOM_cITN_CROP_MONTHLY_SAMPLES[i] = deepcopy(RANDOM_cITN_CROP_MONTHLY)
    RANDOM_HOUSEHOLD_NPC_MONTHLY_SAMPLES[i] = deepcopy(RANDOM_HOUSEHOLD_NPC_MONTHLY)
    RANDOM_HOUSEHOLD_NPC_STD_MONTHLY_SAMPLES[i] = deepcopy(RANDOM_HOUSEHOLD_NPC_STD_MONTHLY)
end

########################################################
# %% Perform Sensitivity Analysis Regression for Chronological Case
########################################################
# Construct altered input_dict
CHRONOLOGICAL_input_dict = deepcopy(input_dict)
CHRONOLOGICAL_input_dict["LLIN_CROP_MONTHLY"] = CHRONOLOGICAL_LLIN_CROP_MONTHLY
CHRONOLOGICAL_input_dict["NET_CROP_MONTHLY"] = CHRONOLOGICAL_NET_CROP_MONTHLY
CHRONOLOGICAL_input_dict["NET_CROP_STD_MONTHLY"] = CHRONOLOGICAL_NET_CROP_STD_MONTHLY
CHRONOLOGICAL_input_dict["cITN_CROP_MONTHLY"] = CHRONOLOGICAL_cITN_CROP_MONTHLY
CHRONOLOGICAL_input_dict["HOUSEHOLD_NPC_MONTHLY"] = CHRONOLOGICAL_HOUSEHOLD_NPC_MONTHLY
CHRONOLOGICAL_input_dict["HOUSEHOLD_NPC_STD_MONTHLY"] = CHRONOLOGICAL_HOUSEHOLD_NPC_STD_MONTHLY

# Do regression
findall(.!ismissing.(CHRONOLOGICAL_input_dict["HOUSEHOLD_NPC_STD_MONTHLY"]))

scatter(input_dict["HOUSEHOLD_NPC_STD_MONTHLY"])

#

println("Commencing SA for Chronological Case...")
bayes_GD(CHRONOLOGICAL_input_dict;
            chain_output_dir = output_dir*"chronological/",
            filename = "$(ISO)_chronological_$(n_data).jld2")
println("Completed SA for Chronological Case.")

########################################################
# %% Perform Sensitivity Analysis Regression for Reverse Chronological
########################################################
# Construct altered input_dict
REV_CHRONOLOGICAL_input_dict = deepcopy(input_dict)
REV_CHRONOLOGICAL_input_dict["LLIN_CROP_MONTHLY"] = REV_CHRONOLOGICAL_LLIN_CROP_MONTHLY
REV_CHRONOLOGICAL_input_dict["NET_CROP_MONTHLY"] = REV_CHRONOLOGICAL_NET_CROP_MONTHLY
REV_CHRONOLOGICAL_input_dict["NET_CROP_STD_MONTHLY"] = REV_CHRONOLOGICAL_NET_CROP_STD_MONTHLY
REV_CHRONOLOGICAL_input_dict["cITN_CROP_MONTHLY"] = REV_CHRONOLOGICAL_cITN_CROP_MONTHLY
REV_CHRONOLOGICAL_input_dict["HOUSEHOLD_NPC_MONTHLY"] = REV_CHRONOLOGICAL_HOUSEHOLD_NPC_MONTHLY
REV_CHRONOLOGICAL_input_dict["HOUSEHOLD_NPC_STD_MONTHLY"] = REV_CHRONOLOGICAL_HOUSEHOLD_NPC_STD_MONTHLY

println("Commencing SA for Reverse Chronological Case...")
bayes_GD(REV_CHRONOLOGICAL_input_dict;
            chain_output_dir = output_dir*"rev_chronological/",
            filename = "$(ISO)_rev_chronological_$(n_data).jld2")
println("Completed SA for Reverse Chronological Case.")

########################################################
# %% Perform Sensitivity Analysis Regression for Reverse Chronological
########################################################
for i in 1:n_samples
    # Construct altered input_dict
    RANDOM_input_dict = deepcopy(input_dict)
    RANDOM_input_dict["LLIN_CROP_MONTHLY"] = RANDOM_LLIN_CROP_MONTHLY_SAMPLES[i]
    RANDOM_input_dict["NET_CROP_MONTHLY"] = RANDOM_NET_CROP_MONTHLY_SAMPLES[i]
    RANDOM_input_dict["NET_CROP_STD_MONTHLY"] = RANDOM_NET_CROP_STD_MONTHLY_SAMPLES[i]
    RANDOM_input_dict["cITN_CROP_MONTHLY"] = RANDOM_cITN_CROP_MONTHLY_SAMPLES[i]
    RANDOM_input_dict["HOUSEHOLD_NPC_MONTHLY"] = RANDOM_HOUSEHOLD_NPC_MONTHLY_SAMPLES[i]
    RANDOM_input_dict["HOUSEHOLD_NPC_STD_MONTHLY"] = RANDOM_HOUSEHOLD_NPC_STD_MONTHLY_SAMPLES[i]

    println("Commencing SA for Random Case $(i)/$(n_samples)...")
    bayes_GD(RANDOM_input_dict;
                chain_output_dir = output_dir*"random/",
                filename = "$(ISO)_random_$(n_data)_sample_$(i).jld2")
    println("Completed SA for Random Case $(i)/$(n_samples).")
end