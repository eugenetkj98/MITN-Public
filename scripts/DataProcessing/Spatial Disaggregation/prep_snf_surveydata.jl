"""
Author: Eugene Tan
Date Created: 2/12/2024
Last Updated: 5/12/2024
Script to normalise and post process extracted spatial household data points and reduce dimension for collinearity
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Import relevant packages
using ProgressBars
using DataFrames
using CSV
using StatsBase
using Statistics
using MultivariateStats
using LinearAlgebra
using SparseArrays

# %% Directories
input_dir = OUTPUT_DATAPREP_DIR*"INLA/"
input_filename = "inla_dataset_snf_adj.csv"
cov_legend_dir = RAW_DATASET_DIR*"INLA/"
cov_legend_filename = COV_LEGEND_FILENAME

output_dir = OUTPUT_DATAPREP_DIR*"INLA/"
output_filename = "inla_dataset_snf_reduced.csv"

# %% Load extracted household data
data = CSV.read(input_dir*input_filename, DataFrame)

# %% Do some covariate transformations
data[:,"cov_ACCESS"] .= log.(data[:,"cov_ACCESS"] .+ 1e-5)
data[:,"cov_ARID"] .= sqrt.(data[:,"cov_ARID"])

cov_legend = CSV.read(cov_legend_dir*cov_legend_filename, DataFrame)

# %% Import normalisation constants and PCA transformation
norm_constants = CSV.read(output_dir*"covariate_normalisation_constants.csv", DataFrame)
proj_matrix = CSV.read(output_dir*"proj_matrix.csv", DataFrame)

# %% Create normalised version of dataset
norm_data = copy(data)
cov_names = col_names[findall(contains.(col_names, "cov_"))]
for cov_name in cov_names
    norm_data[:, cov_name] = data[:, cov_name]./norm_constants[norm_constants.cov .== cov_name,"std"][1]
end

# %% Find all columns with "cov_" prefix (i.e. the covariates)
col_names = names(data)
cov_names = col_names[findall(contains.(col_names, "cov_"))]
n_covs = length(cov_names)
proj_matrix

# %% Calculate projection
proj_M = Matrix(proj_matrix[:,2:end])
proj_vals = (proj_M*(Matrix(norm_data[:, (end-n_covs+1):end])'))'
proj_names = proj_matrix[:,1]

# %% Construct reduced INLA dataset for snf adj
inla_gen_data = hcat(norm_data[:,1:15], DataFrame(proj_vals, proj_names))

# %% Filter out problematic data

filt_inla_gen_data = inla_gen_data[findall(.!isnan.(inla_gen_data.npc) .&& .!isnan.(inla_gen_data.access) .&& 
                                    .!isnan.(inla_gen_data.npc_gap) .&& .!isnan.(inla_gen_data.access_gap)),:]

# %% Save dataset
CSV.write(output_dir*"inla_dataset_snf_adj_reduced.csv", filt_inla_gen_data)
