"""
Author: Eugene Tan
Date Created: 2/12/2024
Last Updated: 5/12/2024
Script to include adjusted access scores to normalise observed use against. 
NEED TO WRITE A BETTER DESCRIPTION
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
input_filename = INLA_REDUCED_DATAPREP_FILENAME

output_dir = OUTPUT_DATAPREP_DIR
output_filename = INLA_USE_REDUCED_DATAPREP_FILENAME

adj_access_raster_dir = OUTPUT_RASTERS_DIR*"final_access/pmodel_access/"

# %% Import data
survey_data = CSV.read(input_dir*input_filename, DataFrame)

# %% Get month and year values
# Define storage variable
adj_access_vals = zeros(size(survey_data)[1])
YEAR_VALS = unique(survey_data.interview_year)

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]
    # Get list of required months
    months = unique(survey_data[survey_data.interview_year .== year,"interview_month"])

    # Select a month
    for month_idx in ProgressBar(1:length(months), leave = false)
        month = months[month_idx]
      
        # Load raster
        adj_access_raster_filename = "NA"
        if month < 10
            adj_access_raster_filename = "adj_access_$(year)_0$(month)_mean.tif"
        else
            adj_access_raster_filename =  "adj_access_$(year)_$(month)_mean.tif"
        end
        adj_access_raster = Raster(adj_access_raster_dir*adj_access_raster_filename)

        entry_idxs = findall((survey_data.interview_year .== year) .& 
                                (survey_data.interview_month .== month))

        temp_adj_access_vals = Vector{Float64}(undef, length(entry_idxs))

        for i in 1:length(entry_idxs)
            idx = entry_idxs[i]

            # Get lat, lon valuyes
            lat, lon = survey_data[idx,["latitude", "longitude"]]

            # Find index of raster corresponding to latlon
            lats = lookup(adj_access_raster, Y)
            lons = lookup(adj_access_raster, X)

            lat_idx = argmin(abs.(lats .- lat))
            lon_idx = argmin(abs.(lons .- lon))

            # Extract required value
            if ismissing(adj_access_raster.data[lon_idx,lat_idx]) # if missing, then just take current access value from survey
                temp_adj_access_vals[i] = survey_data[idx,"access"]
            else
                temp_adj_access_vals[i] = Float64(adj_access_raster.data[lon_idx,lat_idx])
            end
        end

        # Save values to cov storage variable
        adj_access_vals[entry_idxs] .= copy(temp_adj_access_vals)
    end
end

# %% Construct extended data

df_insert = DataFrame(Matrix(Matrix(adj_access_vals')'), ["adj_access"])
output_df = hcat(survey_data[:,1:10],df_insert,survey_data[:,11:end])

CSV.write(output_dir*output_filename, output_df)