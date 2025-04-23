"""
Author: Eugene Tan
Date Created: 2/4/2025
Last Updated: 2/4/2025
Script to run regression for net access
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Import Public Packages
using JLD2
using CSV
using DataFrames

# %% Import Custom Packages
using DataExtractions
using DateConversions
using NetAccessModel
using NetAccessRegression

# Maths packages
using LinearAlgebra
using StatsBase

# %% Get ISO List
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Run Analysis
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

# %%
for i in 1:length(filt_ISOs)
    # Select ISO
    ISO = filt_ISOs[i]

    println("Extracting access data for Country $(i) of $(length(filt_ISOs)) → $(ISO)...")
    
    if ISO ∈ exclusion_ISOs
        println("$(ISO) is on exclusion list. Moving to next country.")
        continue
    else
        # Import net crop data
        input_dict = load(OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
        regression_dict = load(OUTPUT_REGRESSIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")
        
        # Extract net access data
        extract_data_netaccess(ISO, YEAR_START, YEAR_END,
                                                    input_dict, regression_dict, reg_mode = true)
        extract_data_netaccess(ISO, YEAR_START, YEAR_END,
                                                    input_dict, regression_dict, reg_mode = false)
    end
end

println("Extracted all required net access data.")

# %% Aggregate survey data for regressing net access parameters and save .jld2
global init_variables = false

p_h_globaldata = zeros(0,10)
ρ_h_globaldata = zeros(0,10)
μ_h_globaldata = zeros(0,10)
γ_globaldata = zeros(0)

for i in 1:length(filt_ISOs)
    # Select ISO
    ISO = filt_ISOs[i]

    if ISO ∈ exclusion_ISOs
        continue
    else
        # Import extracted data
        net_access_input_dict = load(OUTPUT_EXTRACTIONS_DIR*"access/reg_data/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")
        
        γ_surveys = net_access_input_dict["γ_aggregated"]
        p_h = net_access_input_dict["p_h_aggregated"]
        ρ_h = net_access_input_dict["ρ_h_aggregated"]
        μ_h = net_access_input_dict["μ_h_aggregated"]
        
        h_max = size(p_h)[2]
        # Check if storage variable for aggregated global data has been made.
        if init_variables == false
            global p_h_globaldata = zeros(0,h_max)
            global ρ_h_globaldata = zeros(0,h_max)
            global μ_h_globaldata = zeros(0,h_max)
            global γ_globaldata = zeros(0)

            global init_variables = true
        end

        # Add to aggregate storage variable
        p_h_globaldata = vcat(p_h_globaldata, p_h)
        ρ_h_globaldata = vcat(ρ_h_globaldata, ρ_h)
        μ_h_globaldata = vcat(μ_h_globaldata, μ_h)
        γ_globaldata = vcat(γ_globaldata, γ_surveys)
    end
end

# Save aggregated data
mkpath(OUTPUT_EXTRACTIONS_DIR*"access/reg_data/aggregated_inputs/")
jldsave(OUTPUT_EXTRACTIONS_DIR*"access/reg_data/aggregated_inputs/netaccess_allsurvey_inputs.jld2";
            p_h_globaldata, ρ_h_globaldata, μ_h_globaldata, γ_globaldata)



# %% MCMC Regression for Net Access model
access_survey_globaldata = load(OUTPUT_EXTRACTIONS_DIR*"access/reg_data/aggregated_inputs/netaccess_allsurvey_inputs.jld2")
bayes_access(access_survey_globaldata)

println("Net access regression complete.")

