"""
Author: Eugene Tan
Date Created: 10/4/2025
Last Updated: 10/4/2025
Get BV and MITN pixel level predictions of ITN metrics and compile into dataset for plotting RMSE
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %%
using DataFrames
using CSV
using ProgressMeter
using GeoIO
using Rasters
using ProgressBars
using RasterLookup

# %% Data Directories
# Rasters
mitn_raster_dir = OUTPUT_RASTERS_DIR
bv_raster_dir = BV_OUTPUT_DIR

# Training data
inla_data_dir = OUTPUT_DATAPREP_DIR*"INLA/"
inla_data_filename = "inla_dataset_reduced.csv"

# %% Load survey training data and sort by time
inla_data = CSV.read(inla_data_dir*inla_data_filename, DataFrame)
sort!(inla_data, [order(:interview_year), order(:interview_month)])

# %% Define storage variable to compile results in
df_collection = []

# %% Get number of unique timestamps
YEAR_VALS = unique(inla_data.interview_year)

# %% Go through each year and scrape required data and lookup model predictions
for year_i in ProgressBar(1:length(YEAR_VALS), leave = false)
    year = YEAR_VALS[year_i]

    # Filter data
    filt_data_annual = inla_data[inla_data.interview_year .== year,:]

    # unique months to lookup
    months = unique(filt_data_annual.interview_month)

    # Import annual resolution rasters
    bv_npc_raster = replace_missing(Raster(bv_raster_dir*"nets_per_capita/ITN_$(year)_percapita_nets_mean.tif"), missingval = NaN)
    bv_access_raster = replace_missing(Raster(bv_raster_dir*"access/ITN_$(year)_access_mean.tif"), missingval = NaN)
    bv_use_raster = replace_missing(Raster(bv_raster_dir*"ITN_$(year)_use_mean.tif"), missingval = NaN)

    for month_i in ProgressBar(1:length(months), leave = false)
        month = months[month_i]

        month_str = "$(month)"
        if month < 10
            month_str = "0$(month)"
        end

        # Import monthly resolution rasters
        mitn_npc_raster = replace_missing(Raster(mitn_raster_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        mitn_access_raster = replace_missing(Raster(mitn_raster_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        mitn_use_raster = replace_missing(Raster(mitn_raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_mean.tif"), missingval = NaN)

        # Further filter data to desired month
        filt_data_monthly = filt_data_annual[filt_data_annual.interview_month .== month,:]

        bv_npc_vals = zeros(size(filt_data_monthly)[1])
        bv_access_vals = zeros(size(filt_data_monthly)[1])
        bv_use_vals = zeros(size(filt_data_monthly)[1])

        mitn_npc_vals = zeros(size(filt_data_monthly)[1])
        mitn_access_vals = zeros(size(filt_data_monthly)[1])
        mitn_use_vals = zeros(size(filt_data_monthly)[1])

        for i in 1:size(filt_data_monthly)[1]
            lon, lat = filt_data_monthly[i,["longitude", "latitude"]]

            # BV model values
            
            bv_npc_vals[i] = interp_lookup(bv_npc_raster, lat, lon)
            bv_access_vals[i] = interp_lookup(bv_access_raster, lat, lon)
            bv_use_vals[i] = interp_lookup(bv_use_raster, lat, lon)

            # MITN model values
            # try mitn_npc_raster[X = Near(lon),Y = Near(lat)] # see if available in dimension
            mitn_npc_vals[i] = interp_lookup(mitn_npc_raster, lat, lon)
            mitn_access_vals[i] = interp_lookup(mitn_access_raster, lat, lon)
            mitn_use_vals[i] = interp_lookup(mitn_use_raster, lat, lon)
        end

        # Compile data frame for current month
        df = hcat(filt_data_monthly[:,["ISO","admin1_name", "area_id", "latitude", "longitude", "interview_month", "interview_year", "cluster_sample_wt"]],
                DataFrame(      npc = filt_data_monthly.npc,
                                access = filt_data_monthly.access,
                                use = filt_data_monthly.use,
                                bv_npc = bv_npc_vals,
                                mitn_npc = mitn_npc_vals,
                                bv_access = bv_access_vals,
                                mitn_access = mitn_access_vals,
                                bv_use = bv_use_vals,
                                mitn_use = mitn_use_vals))

        # Add to collection of dataframes
        push!(df_collection, df)
    end
end

# %% Compile data frame into a 
compiled_dataframe = vcat(df_collection...)

# %% Save into outputs
CSV.write(OUTPUT_DIR*"coverage_timeseries/bv_mitn_validation_values.csv", compiled_dataframe)