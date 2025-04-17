"""
Author: Eugene Tan
Date Created: 8/4/2025
Last Updated: 8/4/2025
Combines generate fragments from raster_timeseries_aggregations.jl, sorts and saves dataframe
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Import relevant packages
using ProgressBars
using DataFrames
using CSV

# %% Navigate to directory and find list of all filenames
aggregate_parts_dir = OUTPUT_DIR*"coverage_timeseries/master_extractions_parts/"
output_dir = OUTPUT_DIR*"coverage_timeseries/"

# Get list of filenames
filenames = readdir(aggregate_parts_dir)

# Read each file and add to DataFrame list
df_collection = []

for filename in filenames
    push!(df_collection, CSV.read(aggregate_parts_dir*filename, DataFrame))
end

# Concatenate all dataframes into a master frame 
master_df = vcat(df_collection...)
sort!(master_df, [order(:year), order(:month), order(:ISO), order(:category),order(:admin_name)])
CSV.write(OUTPUT_DIR*"coverage_timeseries/master_extraction.csv", master_df)