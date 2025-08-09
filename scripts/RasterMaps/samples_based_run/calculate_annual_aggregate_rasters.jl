"""
Author: Eugene Tan
Date Created: 19/5/2025
Last Updated: 19/5/2025
Combines monthly rasters into mean annual rasters.
"""
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from TOML file
include(pwd()*"/scripts/read_toml.jl")

# %% Import relevant packages
using ProgressBars
using DataFrames
using CSV
using Rasters
using GeoIO
using StatsBase

# %% File paths
raster_dir = OUTPUT_RASTERS_DIR
output_dir = OUTPUT_RASTERS_DIR

# Output directories
mkpath(raster_dir*"final_npc/mean/annual")
mkpath(raster_dir*"final_access/mean/annual")
mkpath(raster_dir*"final_use/mean/annual")
mkpath(raster_dir*"final_utilisation/mean/annual")

# %% Num of samples to import from INLA outputs for rasters
n_samples = INLA_UNCERTAINTY_N_SAMPLES

# %% Year Bounds
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

# %%
year = parse(Int, ARGS[1])

# Define all required filenames for each metric
netage_raster_monthly_dir = raster_dir*"final_netage/snf_netage/"
npc_raster_monthly_dir = raster_dir*"final_npc/mean/monthly/"
access_raster_monthly_dir = raster_dir*"final_access/mean/monthly/"
use_raster_monthly_dir = raster_dir*"final_use/mean/monthly/"
util_raster_monthly_dir = raster_dir*"final_utilisation/mean/monthly/"

netage_all_filenames = readdir(netage_raster_monthly_dir)
npc_all_filenames = readdir(npc_raster_monthly_dir)
access_all_filenames = readdir(access_raster_monthly_dir)
use_all_filenames = readdir(use_raster_monthly_dir)
util_all_filenames = readdir(util_raster_monthly_dir)

netage_monthly_filenames = netage_all_filenames[occursin.("$(year)", npc_all_filenames)]
npc_monthly_filenames = npc_all_filenames[occursin.("$(year)", npc_all_filenames)]
access_monthly_filenames = access_all_filenames[occursin.("$(year)", access_all_filenames)]
use_monthly_filenames = use_all_filenames[occursin.("$(year)", use_all_filenames)]
util_monthly_filenames = util_all_filenames[occursin.("$(year)", util_all_filenames)]

# Check to make sure nonempty filename lists
if (length(npc_monthly_filenames) < 1) || (length(access_monthly_filenames) < 1) || (length(use_monthly_filenames) < 1) || (length(util_monthly_filenames) < 1)
    println("No complete sample raster sets found for $(year)...skipping. :3")
else
    # Import each set of sample rasters
    netage_monthly_rasters = replace_missing.(Raster.(netage_raster_monthly_dir.*netage_monthly_filenames), missingval = NaN)
    npc_monthly_rasters = replace_missing.(Raster.(npc_raster_monthly_dir.*npc_monthly_filenames), missingval = NaN)
    access_monthly_rasters = replace_missing.(Raster.(access_raster_monthly_dir.*access_monthly_filenames), missingval = NaN)
    use_monthly_rasters = replace_missing.(Raster.(use_raster_monthly_dir.*use_monthly_filenames), missingval = NaN)
    util_monthly_rasters = replace_missing.(Raster.(util_raster_monthly_dir.*util_monthly_filenames), missingval = NaN)

    # Calculate mean raster for each metric
    netage_mean_raster = replace_missing(mean(netage_monthly_rasters), missingval = NaN)
    npc_mean_raster = replace_missing(mean(npc_monthly_rasters), missingval = NaN)
    access_mean_raster = replace_missing(mean(access_monthly_rasters), missingval = NaN)
    use_mean_raster = replace_missing(mean(use_monthly_rasters), missingval = NaN)
    util_mean_raster = replace_missing(mean(util_monthly_rasters), missingval = NaN)

    # Write rasters
    write(output_dir*"final_netage/annual/netage_$(year)_mean.tif", netage_mean_raster, force = true)
    write(output_dir*"final_npc/mean/annual/npc_$(year)_mean.tif", npc_mean_raster, force = true)
    write(output_dir*"final_access/mean/annual/access_$(year)_mean.tif", access_mean_raster, force = true)
    write(output_dir*"final_use/mean/annual/use_$(year)_mean.tif", use_mean_raster, force = true)
    write(output_dir*"final_utilisation/mean/annual/utilisation_$(year)_mean.tif", util_mean_raster, force = true)

    # Progress message
    println("Constructed mean rasters for $(year). Woohoo!")
end