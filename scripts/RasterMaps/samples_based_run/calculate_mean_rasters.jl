"""
Author: Eugene Tan
Date Created: 19/5/2025
Last Updated: 19/5/2025
Combines sampled rasters from posterior samples rasters of ITN metrics
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
mkpath(raster_dir*"final_npc/mean/monthly")
mkpath(raster_dir*"final_access/mean/monthly")
mkpath(raster_dir*"final_use/mean/monthly")
mkpath(raster_dir*"final_utilisation/mean/monthly")

# %% Year Bounds
# YEAR_START = 2004 #YEAR_NAT_START
# YEAR_END = YEAR_NAT_END

# for year in YEAR_START:YEAR_END
#     for month in 1:12

year = parse(Int, ARGS[1])
month = parse(Int, ARGS[2])


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
        else
            # Create tally rasters
            global npc_sum_raster = replace_missing(Raster(npc_raster_post_dir*npc_sample_filenames[1]), missingval = NaN)
            global access_sum_raster = replace_missing(Raster(access_raster_post_dir.*access_sample_filenames[1]), missingval = NaN)
            global use_sum_raster = replace_missing(Raster(use_raster_post_dir.*use_sample_filenames[1]), missingval = NaN)
            global util_sum_raster = replace_missing(Raster(util_raster_post_dir.*util_sample_filenames[1]), missingval = NaN)

            for sample_i in 2:length(npc_sample_filenames)
                println("Loading NPC raster $(sample_i)/$(length(npc_sample_filenames))")
                global npc_sum_raster = npc_sum_raster .+ replace_missing(Raster(npc_raster_post_dir*npc_sample_filenames[sample_i]), missingval = NaN)
            end

            for sample_i in 2:length(access_sample_filenames)
                println("Loading Access raster $(sample_i)/$(length(access_sample_filenames))")
                global access_sum_raster = access_sum_raster .+ replace_missing(Raster(access_raster_post_dir*access_sample_filenames[sample_i]), missingval = NaN)
            end

            for sample_i in 2:length(use_sample_filenames)
                println("Loading Use raster $(sample_i)/$(length(use_sample_filenames))")
                global use_sum_raster = use_sum_raster .+ replace_missing(Raster(use_raster_post_dir*use_sample_filenames[sample_i]), missingval = NaN)
            end

            for sample_i in 2:length(util_sample_filenames)
                println("Loading Utilisation raster $(sample_i)/$(length(util_sample_filenames))")
                global util_sum_raster = util_sum_raster .+ replace_missing(Raster(util_raster_post_dir*util_sample_filenames[sample_i]), missingval = NaN)
            end
            

            # Calculate mean raster for each metric
            npc_mean_raster = replace_missing(npc_sum_raster./length(npc_sample_filenames), missingval = NaN)
            access_mean_raster = replace_missing(access_sum_raster./length(access_sample_filenames), missingval = NaN)
            use_mean_raster = replace_missing(use_sum_raster./length(use_sample_filenames), missingval = NaN)
            util_mean_raster = replace_missing(util_sum_raster./length(util_sample_filenames), missingval = NaN)

            # Write rasters
            write(output_dir*"final_npc/mean/monthly/npc_$(year)_$(month_str)_mean.tif", npc_mean_raster, force = true)
            write(output_dir*"final_access/mean/monthly/access_$(year)_$(month_str)_mean.tif", access_mean_raster, force = true)
            write(output_dir*"final_use/mean/monthly/use_$(year)_$(month_str)_mean.tif", use_mean_raster, force = true)
            write(output_dir*"final_utilisation/mean/monthly/utilisation_$(year)_$(month_str)_mean.tif", util_mean_raster, force = true)

            # Progress message
            println("Constructed mean rasters for $(year)-$(month_str). Woohoo!")
        end


#     end
# end
