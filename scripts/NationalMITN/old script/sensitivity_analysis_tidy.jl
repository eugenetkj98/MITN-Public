"""
Author: Eugene Tan
Date Created: 23/7/2024
Last Updated: 26/8/2024
Script to compile and tidy outputs from sensitivity analysis into a single file 
for each country
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import Public Packages
using JLD2
using CSV
using DataFrames
using ProgressBars

# %% Import Custom Packages
using PlottingFunctions

# Maths packages
using LinearAlgebra
using StatsBase

# %% Get ISO List
ISO_list = String.(CSV.read(raw"C:\Users\ETan\Documents\Prototype Analyses\itn-updated\datasets\ISO_list.csv", DataFrame)[:,1])
exclusion_ISOs = ["CPV","BWA","GNQ","DJI","ETH","SOM","ZAF","SSD"]
# ["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD"]
# GAB excluded
# %% Run Analysis configs
YEAR_START = 2000
YEAR_END = 2023 # Inclusive until end of 2021 December, excludes 2022

SA_reg_dir = "outputs/regressions/crop/sensitivity_analysis/" # Directory for sensitivity analysis regression files
mode_names = Dict("chronological" => "CSA",
                    "random" => "RSA",
                    "cross" => "CVA")
analysis_modes = ["chronological", "random", "cross"]

# %%
for i in 1:length(ISO_list)
    # Select ISO
    ISO = ISO_list[i]

    println("Compiling SA results for Country $(i) of $(length(ISO_list)) → $(ISO)...")
    
    if ISO ∈ exclusion_ISOs
        println("$(ISO) is on exclusion list. Moving to next country.")
        continue
    else
        # Data Cleaning and Tidying
        output_dict = Dict() # Storage variable for output later
        for sa_mode in analysis_modes # SA mode to extract from options: random, chronological, cross (cross validation)
            country_sa_dir = SA_reg_dir*"$(YEAR_START)_$(YEAR_END)/$sa_mode/$ISO" # Directory of analysis files containing SA results for selected country ISO
            dir_entries = readdir(country_sa_dir)
            n_entries = length(dir_entries)
            id_vals = parse.(Int64, dir_entries)

            println("ISO: $(ISO). Analysis Mode: $(sa_mode). #Entries $(n_entries)")

            halflife_quantiles_entries = []
            for j in ProgressBar(1:n_entries, leave = false)

                sa_filename = readdir(country_sa_dir*"/$(dir_entries[j])")[1]
                sa_reg_output = load(country_sa_dir*"/$(dir_entries[j])/$sa_filename")

                # Extract halflife estimates
                chain = sa_reg_output["chain"]
                lifecurve_quantiles_BYNET, halflife_quantiles_BYNET = netlife_posterior_draws(sa_reg_output; t_vals = 0:0.01:5)
                push!(halflife_quantiles_entries,halflife_quantiles_BYNET)
            end

            # Compress halflife entries into a tensor for neatness
            halflife_results_tidy = zeros(n_entries, size(halflife_quantiles_entries[1])...)
            for j in ProgressBar(1:n_entries, leave = false)
                halflife_results_tidy[j,:,:] = halflife_quantiles_entries[j]
            end

            # Compile results into a single dict for plotting
            mode_name = mode_names[sa_mode]
            output_dict[mode_name] = halflife_results_tidy
        end

        # Save extracted results into JLD2 file
        output_filename = "$(ISO)_$(YEAR_START)_$(YEAR_END)_sensitivity_analysis_results.jld2"
        output_dir = "outputs/regressions/crop/sensitivity_analysis/$(YEAR_START)_$(YEAR_END)/compiled_results/"
        save(output_dir*output_filename, output_dict)

        println("SA summary complete for $(ISO). Data saved!")
    end
end

# %%


