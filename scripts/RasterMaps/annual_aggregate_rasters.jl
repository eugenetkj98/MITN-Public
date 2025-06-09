"""
Author: Eugene Tan
Date Created: 16/1/2025
Last Updated: 8/5/2025
Script to combine monthly resolution rasters for ITN metrics into annual rasters
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/read_toml.jl")

# %% Import relevant packages
using ProgressBars
using Rasters
using GeoIO
using JLD2
using StatsBase

# %% Directories
input_dir = OUTPUT_RASTERS_DIR
output_dir = OUTPUT_RASTERS_DIR
mkpath(output_dir*"final_netage/annual/")
mkpath(output_dir*"final_npc/annual/")
mkpath(output_dir*"final_access/annual/")
mkpath(output_dir*"final_use/annual/")
mkpath(output_dir*"final_utilisation/annual/")

# %% Year ranges
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

# %%
for year in ProgressBar(YEAR_START:YEAR_END, leave = false)
    # Get list of month strings
    month_vals = 1:12
    netage_filename_strings = Vector{String}(undef, length(month_vals))
    npc_filename_strings = Vector{String}(undef, length(month_vals))
    access_filename_strings = Vector{String}(undef, length(month_vals))
    use_filename_strings = Vector{String}(undef, length(month_vals))
    utilisation_filename_strings = Vector{String}(undef, length(month_vals))
    for i in 1:length(month_vals)
        month = month_vals[i]
        month_str = ""
        if month < 10
            month_str = "0$(month)"
        else
            month_str = "$(month)"
        end
        netage_filename_strings[i] = "netage_$(year)_$(month_str)_mean.tif"
        npc_filename_strings[i] = "npc_$(year)_$(month_str)_mean.tif"
        access_filename_strings[i] = "access_$(year)_$(month_str)_mean.tif"
        use_filename_strings[i] = "use_$(year)_$(month_str)_mean.tif"
        utilisation_filename_strings[i] = "utilisation_$(year)_$(month_str)_mean.tif"
    end

    # Get list of raster filenames for current year
    netage_rasters = replace_missing.(Raster.(input_dir.*"final_netage/snf_netage/".*netage_filename_strings), missingval = NaN)
    npc_rasters = replace_missing.(Raster.(input_dir.*"final_npc/mean/monthly/".*npc_filename_strings), missingval = NaN)
    access_rasters = replace_missing.(Raster.(input_dir.*"final_access/pmodel_access/".*access_filename_strings), missingval = NaN)
    use_rasters = replace_missing.(Raster.(input_dir.*"final_use/logis_use/".*use_filename_strings), missingval = NaN)
    util_rasters = replace_missing.(Raster.(input_dir.*"final_utilisation/monthly/".*utilisation_filename_strings), missingval = NaN)

    mean_netage_raster = mean(netage_rasters)
    mean_npc_raster = mean(npc_rasters)
    mean_access_raster = mean(access_rasters)
    mean_use_raster = mean(use_rasters)
    mean_util_raster = mean(util_rasters)

    # Save annual raster
    write(output_dir*"final_netage/annual/netage_$(year)_mean.tif", mean_netage_raster, force = true)
    write(output_dir*"final_npc/annual/npc_$(year)_mean.tif", mean_npc_raster, force = true)
    write(output_dir*"final_access/annual/access_$(year)_mean.tif", mean_access_raster, force = true)
    write(output_dir*"final_use/annual/use_$(year)_mean.tif", mean_use_raster, force = true)
    write(output_dir*"final_utilisation/annual/utilisation_$(year)_mean.tif", mean_util_raster, force = true)
end