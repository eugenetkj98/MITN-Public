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
input_dir = OUTPUT_DATAPREP_DIR
input_filename = INLA_DATAPREP_FILENAME
cov_legend_filename = COV_LEGEND_FILENAME

output_dir = OUTPUT_DATAPREP_DIR
output_filename = INLA_REDUCED_DATAPREP_FILENAME

# %% Load extracted household data
data = CSV.read(input_dir*input_filename, DataFrame)

# %% Do some covariate transformations
data[:,"cov_ACCESS"] .= log.(data[:,"cov_ACCESS"] .+ 1e-5)
data[:,"cov_ARID"] .= sqrt.(data[:,"cov_ARID"])

cov_legend = CSV.read(input_dir*cov_legend_filename, DataFrame)

# %% Find all columns with "cov_" prefix (i.e. the covariates)
col_names = names(data)
cov_names = col_names[findall(contains.(col_names, "cov_"))]
n_covs = length(cov_names)

# %% Categorise covariates into static, annual and monthly
cov_types = Vector{String}(undef, n_covs)
for i in 1:n_covs
    cov_types[i] = cov_legend[findfirst(cov_legend.cov_name .== cov_names[i]), "type"]
end

# Static Covariates
idx_static = findall(cov_types .== "static")
cov_static = cov_names[idx_static]

# Annual Covariates
idx_annual = findall(cov_types .== "annual")
cov_annual = cov_names[idx_annual]

# Monthly Covariates
idx_monthly = findall(cov_types .== "monthly")
cov_monthly = cov_names[idx_monthly]

# %% Calculate correlation matrices
function zeromean_unitvar_norm(mat)
    output = zeros(size(mat))
    μ_vals = zeros(size(mat)[2])
    σ_vals = zeros(size(mat)[2])
    
    for i in 1:size(mat)[2]
        μ = mean(mat[:,i])
        σ = std(mat[:,i])
        output[:,i] .= (mat[:,i])./σ#(mat[:,i] .- μ)./σ
        μ_vals[i] = μ
        σ_vals[i] = σ
    end

    return output, μ_vals, σ_vals
end

# Normalise Data
data_static, μ_static, σ_static = zeromean_unitvar_norm(data[:, cov_static])
data_annual, μ_annual, σ_annual = zeromean_unitvar_norm(data[:, cov_annual])
data_monthly, μ_monthly, σ_monthly = zeromean_unitvar_norm(data[:, cov_monthly])

# fix NaNs in zero std data
data_static[findall(isnan.(data_static))] .= 0
data_annual[findall(isnan.(data_annual))] .= 0
data_monthly[findall(isnan.(data_monthly))] .= 0

# Calculate correlation matrices
corr_static = cor(Matrix(data_static))
corr_annual = cor(Matrix(data_annual))
corr_monthly = cor(Matrix(data_monthly))

# Fix NaNs in correlation matrices
corr_static[findall(isnan.(corr_static))] .= 0
corr_annual[findall(isnan.(corr_annual))] .= 0
corr_monthly[findall(isnan.(corr_monthly))] .= 0


# %% Find transformation to get normalised, reduced covariate set
function get_transform_matrix(data_mat, corr_mat; ρ = 0.4, pratio = 0.7)
    # Find all variables whose correlation is above threshold
    ci = findall((abs.(corr_mat)-I).>ρ)

    # Note down index of variables
    corr_idx = unique(vcat(getindex.(ci, 1), getindex.(ci, 2)))
    uncorr_idx = setdiff(1:size(corr_mat)[1], corr_idx)

    if length(corr_idx) > 0
        # Calculate principal components
        pca = fit(PCA, Matrix(data_mat[:, corr_idx])', maxoutdim = length(corr_idx), pratio = pratio)

        # Construct transformation matrix from raw covariates to components
        n_components = length(uncorr_idx) + size(pca.proj)[2]
        M = zeros(n_components, size(corr_mat)[1])

        for i in 1:length(uncorr_idx)
            M[i, uncorr_idx[i]] = 1
        end
        for i in 1:size(pca.proj)[2]
            M[length(uncorr_idx)+i,corr_idx] .= pca.proj[:,i]
        end

        return M
    else # There were no significant correlated components
        return Float64.(Matrix(I(length(uncorr_idx))))
    end
end

M_static = get_transform_matrix(data_static, corr_static)
M_annual = get_transform_matrix(data_annual, corr_annual)
M_monthly = get_transform_matrix(data_monthly, corr_monthly)

M = blockdiag(sparse(M_static), sparse(M_annual), sparse(M_monthly))

# %% Construct output datasets


# Calculate reduced components
red_comps = (M * hcat(data_static, data_annual, data_monthly)')'

# Construct list of covariate names
old_cov_names = vcat(cov_static, cov_annual, cov_monthly)
new_cov_names = []
for i in 1:size(M_static)[1]
    push!(new_cov_names, "static_$(i)")
end

for i in 1:size(M_annual)[1]
    push!(new_cov_names, "annual_$(i)")
end

for i in 1:size(M_monthly)[1]
    push!(new_cov_names, "monthly_$(i)")
end

# Reduced inla dataset
output_dataset = hcat(data[:,1:15],DataFrame(red_comps, new_cov_names))

# Reference matrix for component projections
M_matrix = hcat(DataFrame(raw_cov = new_cov_names), DataFrame(M, old_cov_names))
M_matrix
# %% Save results
# Write normalised covariate values
CSV.write(output_dir*output_filename, output_dataset)
# Save normalised projection matrices
CSV.write(output_dir*"proj_matrix.csv", M_matrix)

# %% Normalisation constants
norm_constants = DataFrame(cov = old_cov_names, mean = vcat(μ_static, μ_annual, μ_monthly), std = vcat(σ_static, σ_annual, σ_monthly))
CSV.write(output_dir*"covariate_normalisation_constants.csv", norm_constants)
