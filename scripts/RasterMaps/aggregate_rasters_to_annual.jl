"""
Author: Eugene Tan
Date Created: 28/1/2025
Last Updated: 28/1/2025
Script to aggregate/summarise monthly raster estimates into annual resolution.
"""
# %% Import packages
using Rasters

# %% File paths
raster_dir = "outputs/rasters/"
output_dir = "outputs/rasters/"

# Make required paths to save in (if they don't already exist)
mkpath(raster_dir*"final_npc/annual/")
mkpath(raster_dir*"final_access/annual/")
mkpath(raster_dir*"final_use/annual/")

# %% Time bounds
YEAR_START = 2020
YEAR_END = 2023#2023

# %%
for year in ProgressBar(YEAR_START:YEAR_END)
    println("Compiling annual rasters for year $(year)...")

    # Define temporary storage variables to store monthly rasters in a given year
    npc_monthly_mean_rasters = Vector{Raster}(undef, 12)
    npc_monthly_upper_rasters = Vector{Raster}(undef, 12)
    npc_monthly_lower_rasters = Vector{Raster}(undef, 12)

    access_monthly_mean_rasters = Vector{Raster}(undef, 12)
    access_monthly_upper_rasters = Vector{Raster}(undef, 12)
    access_monthly_lower_rasters = Vector{Raster}(undef, 12)

    use_monthly_mean_rasters = Vector{Raster}(undef, 12)
    use_monthly_upper_rasters = Vector{Raster}(undef, 12)
    use_monthly_lower_rasters = Vector{Raster}(undef, 12)

    # Import rasters for each year
    for month in 1:12
        # Get required month string based on filenames convention
        month_str = "$(month)"
        if month < 10
            month_str = "0$(month)"
        end

        # Import required rasters and store in appropriate storage variable
        npc_monthly_mean_rasters[month] = Raster(raster_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_mean.tif")
        npc_monthly_upper_rasters[month] = Raster(raster_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_upper.tif")
        npc_monthly_lower_rasters[month] = Raster(raster_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_lower.tif")

        access_monthly_mean_rasters[month] = Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_mean.tif")
        access_monthly_upper_rasters[month] = Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_upper.tif")
        access_monthly_lower_rasters[month] = Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_lower.tif")

        use_monthly_mean_rasters[month] = Raster(raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_mean.tif")
        use_monthly_upper_rasters[month] = Raster(raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_upper.tif")
        use_monthly_lower_rasters[month] = Raster(raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_lower.tif")
    end

    # Calculate annual raster as average across all monthly raster in a given year
    npc_annual_mean_raster = mean(npc_monthly_mean_rasters)
    npc_annual_upper_raster = mean(npc_monthly_upper_rasters)
    npc_annual_lower_raster = mean(npc_monthly_lower_rasters)

    access_annual_mean_raster = mean(access_monthly_mean_rasters)
    access_annual_upper_raster = mean(access_monthly_upper_rasters)
    access_annual_lower_raster = mean(access_monthly_lower_rasters)

    use_annual_mean_raster = mean(use_monthly_mean_rasters)
    use_annual_upper_raster = mean(use_monthly_upper_rasters)
    use_annual_lower_raster = mean(use_monthly_lower_rasters)

    # Save/Write rasters
    write(raster_dir*"final_npc/annual/adj_npc_$(year)_mean.tif", npc_annual_mean_raster)
    write(raster_dir*"final_npc/annual/adj_npc_$(year)_upper.tif", npc_annual_upper_raster)
    write(raster_dir*"final_npc/annual/adj_npc_$(year)_lower.tif", npc_annual_lower_raster)

    write(raster_dir*"final_access/annual/adj_access_$(year)_mean.tif", access_annual_mean_raster)
    write(raster_dir*"final_access/annual/adj_access_$(year)_upper.tif", access_annual_upper_raster)
    write(raster_dir*"final_access/annual/adj_access_$(year)_lower.tif", access_annual_lower_raster)

    write(raster_dir*"final_use/annual/use_$(year)_mean.tif", use_annual_mean_raster)
    write(raster_dir*"final_use/annual/use_$(year)_upper.tif", use_annual_upper_raster)
    write(raster_dir*"final_use/annual/use_$(year)_lower.tif", use_annual_lower_raster)
end


