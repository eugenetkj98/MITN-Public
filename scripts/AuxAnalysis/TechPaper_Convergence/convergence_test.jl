"""
Author: Eugene Tan
Date Created: 16/6/2024
Last Updated: 16/6/2024
Code to do multiple randomly initiated National SNF runs to show convergence behaviour
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Global Constants
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
ISO = "SEN"
exclusion_ISOs = EXCLUSION_ISOS

# %% Run Analysis
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

# %% 
output_dir =  OUTPUT_REGRESSIONS_DIR*"crop/convergence_test/"

# %% Regression
# Load extracted data
input_dict = JLD2.load(OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")

# %%
# Sample index
sample_idx = ARGS[1]
bayes_GD(input_dict;
            chain_output_dir = output_dir,
            filename = "$(ISO)_convergence_sample_$(sample_idx).jld2",
            verbose = true, N_EPOCHS = 10)
