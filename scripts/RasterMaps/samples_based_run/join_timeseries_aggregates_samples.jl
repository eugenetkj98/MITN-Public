"""
Author: Eugene Tan
Date Created: 20/5/2025
Last Updated: 20/5/2025
Combines generate fragments from raster_timeseries_aggregations.jl, sorts and saves dataframe
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from TOML file
include(pwd()*"/scripts/read_toml.jl")

# %% Import relevant packages
using ProgressBars
using DataFrames
using CSV

# %% Navigate to directory and find list of all filenames
aggregate_parts_dir_continent = OUTPUT_DIR*"coverage_timeseries/master_extractions_parts/continent/"
aggregate_parts_dir_admin0 = OUTPUT_DIR*"coverage_timeseries/master_extractions_parts/admin0/"
output_dir = OUTPUT_DIR*"coverage_timeseries/"

#############################################
# %% Join datasets
#############################################
# Get list of filenames
filenames_continent = readdir(aggregate_parts_dir_continent)
filenames_admin0 = readdir(aggregate_parts_dir_admin0)

# Read each file and add to DataFrame list
df_collection = []

for filename in filenames_continent
    push!(df_collection, CSV.read(aggregate_parts_dir_continent*filename, DataFrame))
end

for filename in filenames_admin0
    push!(df_collection, CSV.read(aggregate_parts_dir_admin0*filename, DataFrame))
end

# Concatenate all dataframes into a master frame and save
println("Compiling data extraction parts into a single CSV.")
master_df = vcat(df_collection...)
sort!(master_df, [order(:category), order(:year), order(:month), order(:ISO), order(:admin_name)])
CSV.write(output_dir*"master_extraction.csv", master_df)
println("Saved continent time series extraction at: $(OUTPUT_DIR*"coverage_timeseries/master_extraction_continent.csv")")
