"""
Author: Eugene Tan
Date Created: 12/9/2024
Last Updated: 12/9/2024
Script to run regression for sensitivity analysis across all countries.
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
ISO_list = String.(CSV.read(raw"C:\Users\ETan\Documents\Prototype Analyses\itn-updated\datasets\ISO_list.csv", DataFrame)[:,1])
exclusion_ISOs = ["CPV","BWA","GNQ","DJI","ETH","SOM","ZAF","SSD"]
# ["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD"]
# GAB excluded
# %% Run Analysis
YEAR_START = 2000
YEAR_END = 2023 # Inclusive until end of 2021 December, excludes 2022

# %% Extract all instances for sensitivity analysis
for i in 1:length(ISO_list)
    # Select ISO
    ISO = ISO_list[i]

    println("SA MODE! Extracting data for Country $(i) of $(length(ISO_list)) → $(ISO)...")
    
    if ISO ∈ exclusion_ISOs
        println("$(ISO) is on exclusion list. Moving to next country.")
        continue
    else
        # Data Extraction
        SA_extract_data_netcrop(ISO, YEAR_START, YEAR_END)

        println("Net Crop data extraction complete for $(ISO). Data saved")
    end
end

# %%
mode_names = Dict("chronological" => "CSA",
                    "random" => "RSA",
                    "cross" => "CVA")

# Base directory for crop data
crop_data_dir = "outputs/extractions/crop/sensitivity_analysis/$(YEAR_START)_$(YEAR_END)/"

analysis_modes = ["chronological", "random", "cross"]

# %% Run sensitivity analysis net crop regression for each setting of number of surveys etc.
for i in 1:length(ISO_list)
    # Select ISO
    ISO = ISO_list[i]

    println("SA_MODE! Fitting model for Country $(i) of $(length(ISO_list)) → $(ISO)...")
    
    if ISO ∈ exclusion_ISOs
        println("$(ISO) is on exclusion list. Moving to next country.")
        continue
    else
        # Stock and Flow regression for each mode

        for mode in analysis_modes
            mode_suffix = mode_names[mode]
        
            # Get list and number of files
            n_analyses = length(readdir(crop_data_dir*mode*"/$(ISO)"))
            n_vals = parse.(Int64, readdir(crop_data_dir*mode*"/$(ISO)"))
        
            # Compute MCMC regression for each input file
            for j in 1:length(n_vals)
                
                # Get n_survey value
                n = n_vals[j]
        
                println("Country: $(ISO) [$(j)/$(length(ISO_list))].  Current Mode: $(mode_suffix). n_surveys $(n)/$(maximum(n_vals))")
        
                # Locate extracted data and Import
                dict_dir = crop_data_dir*mode*"/$(ISO)/$(n)/"
                dict_filename = readdir(dict_dir)[1]
                dict_path = dict_dir*dict_filename
                input_dict = load(dict_path)
        
                # Define location to save data to and corresponding file name for MCMC chain
                save_filename = "$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains_$(mode_suffix).jld2"
                chain_output_dir = "outputs/regressions/crop/sensitivity_analysis/$(YEAR_START)_$(YEAR_END)/"*mode*"/$(ISO)/$(n)/"
        
                # MCMC regression
                bayes_GD(input_dict; chain_output_dir = chain_output_dir, filename = save_filename)
            end
        end

        # %%
        println("Net Crop Regression complete for $(ISO). Data saved")
    end
end
