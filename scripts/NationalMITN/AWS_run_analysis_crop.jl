"""
Author: Eugene Tan
Date Created: 2/4/2025
Last Updated: 2/4/2025
Script to run regression for net crop
"""
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from TOML file
include(pwd()*"/scripts/read_toml.jl")

###########################################
# %% Load Packages
###########################################
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

###########################################
# %% Define Directories, paths and parameters
###########################################
# %% Get ISO to analyse from argument input
ISO = ARGS[1]
exclusion_ISOs = EXCLUSION_ISOS

# %% Run Analysis
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

###########################################
# %% Perform data extraction and preparation + Regression
###########################################

# %% Net Crop Data Extraction
println("Extracting Data for Country $(ISO).")
flush(stdout)
extract_data_netcrop(ISO, YEAR_START, YEAR_END)
println("Net Crop Data Extraction complete for $(ISO). Data saved")
flush(stdout)

# %% Regression
if ISO ∈ exclusion_ISOs
    println("$(ISO) is on exclusion list. Ending analysis.")
    flush(stdout)
else
    println("Fitting model for Country $(ISO) with $(Threads.nthreads())")
    flush(stdout)

    # Load extracted data
    input_dict = load(OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")

    # Net crop regression
    bayes_GD(input_dict, save_output = true, N_EPOCHS = 5, n_chains = Threads.nthreads())

    println("Net Crop Regression complete for $(ISO). Data saved")
    flush(stdout)
end


