"""
Author: Eugene Tan
Date Created: 23/7/2024
Last Updated: 26/8/2024
Need to write documentation
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

# %% Get ISO List
ISO_list = String.(CSV.read("datasets/ISO_list.csv", DataFrame)[:,1])
exclusion_ISOs = []
# ["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD"]
# GAB excluded
# %% Run Analysis
YEAR_START = 2000
YEAR_END = 2023 # Inclusive until end of 2021 December, excludes 2022

# %%
for i in 1:length(ISO_list)
    # Select ISO
    ISO = ISO_list[i]

    println("Fitting model for Country $(i) of $(length(ISO_list)) → $(ISO)...")
    
    if ISO ∈ exclusion_ISOs
        println("$(ISO) is on exclusion list. Moving to next country.")
        continue
    else
        # Stock and Flow regression
        extract_data_netcrop(ISO, YEAR_START, YEAR_END)
        input_dict = load("outputs/extractions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")

        bayes_GD(input_dict, save_output = true, N_EPOCHS = 5)
        regression_dict = load("outputs/regressions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")
        
        extract_data_netaccess(ISO, YEAR_START, YEAR_END,
                                                    input_dict, regression_dict, reg_mode = true)
        extract_data_netaccess(ISO, YEAR_START, YEAR_END,
                                                    input_dict, regression_dict, reg_mode = false)
        # net_access_input_dict = load("outputs/extractions/access/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")

        # %%
        println("Net Crop Regression complete for $(ISO). Data saved")
    end
end

# %% Aggregate survey data for regressing net access parameters and save .jld2
init_variables = false

p_h_globaldata = zeros(0,10)
ρ_h_globaldata = zeros(0,10)
μ_h_globaldata = zeros(0,10)
γ_globaldata = zeros(0)

for i in 1:length(ISO_list)
    # Select ISO
    ISO = ISO_list[i]

    if ISO ∈ exclusion_ISOs
        continue
    else
        # Import extracted data
        net_access_input_dict = load("outputs/extractions/access/reg_data/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")
        
        γ_surveys = net_access_input_dict["γ_aggregated"]
        p_h = net_access_input_dict["p_h_aggregated"]
        ρ_h = net_access_input_dict["ρ_h_aggregated"]
        μ_h = net_access_input_dict["μ_h_aggregated"]
        
        h_max = size(p_h)[2]
        # Check if storage variable for aggregated global data has been made.
        if init_variables == false
            p_h_globaldata = zeros(0,h_max)
            ρ_h_globaldata = zeros(0,h_max)
            μ_h_globaldata = zeros(0,h_max)
            γ_globaldata = zeros(0)

            init_variables = true
        end

        # Add to aggregate storage variable
        p_h_globaldata = vcat(p_h_globaldata, p_h)
        ρ_h_globaldata = vcat(ρ_h_globaldata, ρ_h)
        μ_h_globaldata = vcat(μ_h_globaldata, μ_h)
        γ_globaldata = vcat(γ_globaldata, γ_surveys)
    end
end

# Save aggregated data
mkpath("outputs/extractions/access/reg_data/aggregated_inputs/")
jldsave("outputs/extractions/access/reg_data/aggregated_inputs/netaccess_allsurvey_inputs.jld2";
            p_h_globaldata, ρ_h_globaldata, μ_h_globaldata, γ_globaldata)

# %% MCMC Regression for Net Access model
access_survey_globaldata = load("outputs/extractions/access/reg_data/aggregated_inputs/netaccess_allsurvey_inputs.jld2")
bayes_access(access_survey_globaldata)

# println("Analysis Complete. Check directory for outputs.")

