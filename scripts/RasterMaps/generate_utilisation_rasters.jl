"""
Author: Eugene Tan
Date Created: 27/3/2025
Last Updated: 27/3/2025
Combine Use and Access rasters to get overall utilisation rasters on a monthly basis
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Import relevant packages
using ProgressBars
using DataFrames
using CSV
using Rasters
using Shapefile
using GeoInterface
using GeoIO
using JLD2
using StatsBase

# Input and output directory for rasters
raster_dir = OUTPUT_RASTERS_DIR
output_dir = OUTPUT_RASTERS_DIR

# Make paths tosave location if doesn't already exist
mkpath(output_dir*"final_utilisation/monthly/")

# %% Time bounds
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

# %%
for year in ProgressBar(YEAR_START:YEAR_END)
    for month in 1:12
        month_str = ""
        if month < 10
            month_str = "0$(month)"
        else
            month_str = "$(month)"
        end

        # Get filename
        access_filename = "access_$(year)_$(month_str)_mean.tif"
        use_filename = "use_$(year)_$(month_str)_mean.tif"

        # Import Rasters
        access_raster = replace_missing(Raster(raster_dir*"final_access/pmodel_access/"*access_filename), missingval = NaN)
        use_raster = use_rasters = replace_missing(Raster(raster_dir*"final_use/logis_use/"*use_filename), missingval = NaN)

        # Calculate utilisation raster
        mean_util_raster = use_raster./access_raster

        # Save utilisation raster
        write(output_dir*"final_utilisation/monthly/utilisation_$(year)_$(month_str)_mean.tif", mean_util_raster, force = true)
    end
end