"""
Author: Eugene Tan
Date Created: 19/5/2025
Last Updated: 19/5/2025
Combines sampled rasters from posterior samples rasters of ITN metrics
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
using GeoIO
using StatsBase

# %% File paths
raster_dir = OUTPUT_RASTERS_DIR
output_dir = OUTPUT_RASTERS_DIR

# Output directories
mkpath(raster_dir*"final_npc/mean/monthly")
mkpath(raster_dir*"final_access/mean/monthly")
mkpath(raster_dir*"final_use/mean/monthly")
mkpath(raster_dir*"final_utilisation/mean/monthly")


# %% Year Bounds
YEAR_START = 2004 #YEAR_NAT_START
YEAR_END = YEAR_NAT_END

for year in YEAR_START:YEAR_END
    for month in 1:12
        # Define month substring (Convention)
        month_str = "$(month)"
        if month < 10
            month_str = "0$(month)"
        end

        # Define all required filenames for each metric
        npc_raster_post_dir = raster_dir*"final_npc/posterior_samples/"
        access_raster_post_dir = raster_dir*"final_access/posterior_samples/"
        use_raster_post_dir = raster_dir*"final_use/posterior_samples/"
        util_raster_post_dir = raster_dir*"final_utilisation/posterior_samples/"

        npc_all_filenames = readdir(npc_raster_post_dir)
        access_all_filenames = readdir(access_raster_post_dir)
        use_all_filenames = readdir(use_raster_post_dir)
        util_all_filenames = readdir(util_raster_post_dir)

        npc_sample_filenames = npc_all_filenames[occursin.("$(year)_$(month_str)", npc_all_filenames)]
        access_sample_filenames = access_all_filenames[occursin.("$(year)_$(month_str)", access_all_filenames)]
        use_sample_filenames = use_all_filenames[occursin.("$(year)_$(month_str)", use_all_filenames)]
        util_sample_filenames = util_all_filenames[occursin.("$(year)_$(month_str)", util_all_filenames)]

        # Check to make sure nonempty filename lists
        if (length(npc_sample_filenames) < 1) || (length(access_sample_filenames) < 1) || (length(use_sample_filenames) < 1) || (length(util_sample_filenames) < 1)
            println("No complete sample raster sets found for $(year)-$(month_str)...skipping. :3")
            continue
        else
            # Import each set of sample rasters
            npc_sample_rasters = replace_missing.(Raster.(npc_raster_post_dir.*npc_sample_filenames), missingval = NaN)
            access_sample_rasters = replace_missing.(Raster.(access_raster_post_dir.*access_sample_filenames), missingval = NaN)
            use_sample_rasters = replace_missing.(Raster.(use_raster_post_dir.*use_sample_filenames), missingval = NaN)
            util_sample_rasters = replace_missing.(Raster.(util_raster_post_dir.*util_sample_filenames), missingval = NaN)

            # Calculate mean raster for each metric
            npc_mean_raster = replace_missing(mean(npc_sample_rasters), missingval = NaN)
            access_mean_raster = replace_missing(mean(access_sample_rasters), missingval = NaN)
            use_mean_raster = replace_missing(mean(use_sample_rasters), missingval = NaN)
            util_mean_raster = replace_missing(mean(util_sample_rasters), missingval = NaN)

            # Write rasters
            write(output_dir*"final_npc/mean/monthly/npc_$(year)_$(month_str)_mean.tif", npc_mean_raster, force = true)
            write(output_dir*"final_access/mean/monthly/access_$(year)_$(month_str)_mean.tif", access_mean_raster, force = true)
            write(output_dir*"final_use/mean/monthly/use_$(year)_$(month_str)_mean.tif", use_mean_raster, force = true)
            write(output_dir*"final_utilisation/mean/monthly/utilisation_$(year)_$(month_str)_mean.tif", util_mean_raster, force = true)

            # Progress message
            println("Constructed mean rasters for $(year)-$(month_str). Woohoo!")
        end
    end
end
