"""
Author: Eugene Tan
Date Created: 6/4/2025
Last Updated: 6/4/2025
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

# Custom packages
using DateConversions

# %% File paths
# Save file name
output_dir = OUTPUT_DATAPREP_DIR*"INLA/"
output_filename = "inla_dataset_snf_adj.csv"

# Population Rasters directory
pop_dir = POPULATION_RASTER_DIR

# Rasters directory
raster_dir = OUTPUT_RASTERS_DIR

# Region boundaries
admin0_shapes_geoIO = GeoIO.load(ADMIN0_SHAPEFILE)
admin1_shapes_geoIO = GeoIO.load(ADMIN1_SHAPEFILE)

# %% Perform draws and save outputs. Filter out unwanted countries
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = ["CPV","ZAF"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Construct spatial sampling distribution based on population in each pixel relative to national.

# %% Loop to first construct SNF block map up to subnational resolution
# Import log model npc rasters (any calibrated raster will do 5x5km resolution)
raster_base = replace_missing(Raster(OUTPUT_RASTERS_DIR*"inla_logmodel_npc/NPC_logmodel_$(2000)_mean.tif"), missingval = NaN)

# %% Import population raster
population_raster = replace_missing(Raster(pop_dir*"WorldPop_UNAdj_v3_DRC_fix.$(2020).Annual.Data.5km.sum.tif"), missingval = NaN)

# %% Storage variable for weight rasters
weight_rasters = []

for ISO_i in ProgressBar(1:length(filt_ISOs))
    # %% For each country calculate population distribution probability
    ISO = filt_ISOs[ISO_i]

    # Get Country geometry to mask raster
    admin0_geometry = admin0_shapes_geoIO[admin0_shapes_geoIO.ISO .== ISO,:].geometry

    # Get masked population raster
    pop_masked = Rasters.trim(mask(population_raster, with = admin0_geometry); pad=0)

    # Calculate total population of country
    total_pop = sum(pop_masked[.!isnan.(pop_masked)])

    # Construct probability weight matrix and save to storage variable
    country_weight_raster = pop_masked./total_pop
    push!(weight_rasters, country_weight_raster)
end

# %%
combined_weight_raster = mosaic(first, weight_rasters..., atol = 0.01)

# %% Construct sampling settings and weights
n_samples = 3000

# Extract all nonnan locations
nonnan_idxs = findall(.!isnan.(combined_weight_raster))
sampling_weights = Weights(combined_weight_raster[nonnan_idxs])

# %% Sample and extract data, save to data frame
df_collection = []
for year in ProgressBar(YEAR_NAT_START:YEAR_NAT_END)
    for month in ProgressBar(1:12)
        month_str = "$(month)"
        if month < 10
            month_str = "0$(month)"
        end

        # Import npc and access rasters
        npc_snf_raster = replace_missing(Raster(raster_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        access_snf_raster = replace_missing(Raster(raster_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        adj_npc_mean_raster = replace_missing(Raster(raster_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        adj_access_mean_raster = replace_missing(Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_mean.tif"), missingval = NaN)

        # Sample from weighted distribution
        sampled_idxs = sample(nonnan_idxs, sampling_weights, n_samples)

        # Extract all lat and lon values
        lat_vals = zeros(length(sampled_idxs))
        lon_vals = zeros(length(sampled_idxs))
        combined_weight_raster
        for i in 1:length(sampled_idxs)
            lat_vals[i] = lookup(combined_weight_raster, Y)[sampled_idxs[i][2]]
            lon_vals[i] = lookup(combined_weight_raster, X)[sampled_idxs[i][1]]
        end

        # Storage variables for npc and access
        npc_snf_vals = zeros(n_samples)
        access_snf_vals = zeros(n_samples)
        npc_vals = zeros(n_samples)
        access_vals = zeros(n_samples)

        # Extract required data values from rasters
        # Get list of lat lons from npc and access rasters
        npc_model_lats = lookup(adj_npc_mean_raster, Y)
        npc_model_lons = lookup(adj_npc_mean_raster, X)
        access_model_lats = lookup(adj_access_mean_raster, Y)
        access_model_lons = lookup(adj_access_mean_raster, X)

        # Retrieve npc and access from each sampled points
        for i in 1:n_samples
            lat = lat_vals[i]
            lon = lon_vals[i]

            # Lookup raster idx location
            npc_model_lat_idx = argmin(abs.(npc_model_lats .- lat))
            npc_model_lon_idx = argmin(abs.(npc_model_lons .- lon))
            access_model_lat_idx = argmin(abs.(access_model_lats .- lat))
            access_model_lon_idx = argmin(abs.(access_model_lons .- lon))

            # Use idx to sample value
            npc_snf_vals[i] = Float64(npc_snf_raster[npc_model_lon_idx, npc_model_lat_idx])
            access_snf_vals[i] = Float64(access_snf_raster[access_model_lon_idx, access_model_lat_idx])
            npc_vals[i] = Float64(adj_npc_mean_raster[npc_model_lon_idx, npc_model_lat_idx])
            access_vals[i] = Float64(adj_access_mean_raster[access_model_lon_idx, access_model_lat_idx])
        end


        df = DataFrame(ISO = "XXX", admin1_name = "generated", area_id = 123456789,
                                latitude = lat_vals, longitude = lon_vals,
                                interview_month = month, interview_year = year,
                                cluster_sample_wt = 1,
                                npc = npc_vals, access = access_vals, use = NaN, 
                                npc_gap = npc_vals .- npc_snf_vals, 
                                access_gap = access_vals .- access_snf_vals,
                                use_gap = NaN)
        push!(df_collection, df)
    end
end

# Concatenate all dataframe fragments
hh_data_summary = vcat(df_collection...)

####################################
# %% Get Year bounds to optimise data extraction time
####################################
YEAR_VALS = sort(unique(hh_data_summary.interview_year))

####################################
# %% Get Covariates from Raster at each generated survey value
####################################

#######
# %% Accessibility to Cities
#######
global cov_raster_filepath = ACCESSIBILITY_RASTER_DIR
global cov_raster_filename = ACCESSIBILITY_COV_FILENAME

# Load Rasters
global cov_raster = RasterStack(cov_raster_filepath*cov_raster_filename)

# Change default missing values to NaN so that summation is not affected
global cov_raster = replace_missing(cov_raster, missingval=NaN)

# get original array resolution
# raster_res = size(cov_raster)
downsampling_ratio = 5
# new_raster_res = round.(Int, raster_res./downsampling_ratio)

# resample and overwrite into new raster
global cov_raster_m = Raster(Rasters.aggregate(Rasters.Center(), cov_raster, downsampling_ratio; skipmissingval=true, progress=false))[:,:,1]
global cov_raster_m = replace_missing(cov_raster_m, missingval = -9999)

# Define storage variable
ACCESS_cov = zeros(size(hh_data_summary)[1])

# %% for each row in hh_data_summary
for row_i in 1:size(hh_data_summary)[1]
    # Get lat, lon valuyes
    lat, lon = hh_data_summary[row_i,["latitude", "longitude"]]

    # Find index of raster corresponding to latlon
    lats = lookup(cov_raster_m, Y)
    lons = lookup(cov_raster_m, X)

    lat_idx = argmin(abs.(lats .- lat))
    lon_idx = argmin(abs.(lons .- lon))

    # Extract required value
    ACCESS_cov[row_i] = Float64(cov_raster_m[lon_idx, lat_idx])
end


#######
# %% PET
#######
global cov_raster_filepath = COV_PET_DIR
global cov_raster_filename = COV_PET_FILENAME

# Load Rasters
global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

# Define storage variable
PET_cov = zeros(size(hh_data_summary)[1])

# %% for each row in hh_data_summary
for row_i in 1:size(hh_data_summary)[1]
    # Get lat, lon valuyes
    lat, lon = hh_data_summary[row_i,["latitude", "longitude"]]

    # Find index of raster corresponding to latlon
    lats = lookup(cov_raster, Y)
    lons = lookup(cov_raster, X)

    lat_idx = argmin(abs.(lats .- lat))
    lon_idx = argmin(abs.(lons .- lon))

    # Extract required value
    PET_cov[row_i] = Float64(cov_raster[lon_idx, lat_idx])
end

PET_cov = max.(0, PET_cov)

#######
# %% Aridity Covariate
#######
global cov_raster_filepath = COV_ARID_DIR
global cov_raster_filename = COV_ARID_FILENAME

# Load Rasters
global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

# Define storage variable
ARID_cov = zeros(size(hh_data_summary)[1])

# %% for each row in hh_data_summary
for row_i in 1:size(hh_data_summary)[1]
    # Get lat, lon valuyes
    lat, lon = hh_data_summary[row_i,["latitude", "longitude"]]

    # Find index of raster corresponding to latlon
    lats = lookup(cov_raster, Y)
    lons = lookup(cov_raster, X)

    lat_idx = argmin(abs.(lats .- lat))
    lon_idx = argmin(abs.(lons .- lon))

    # Extract required value
    ARID_cov[row_i] = Float64(cov_raster[lon_idx, lat_idx])
end

ARID_cov = max.(0, ARID_cov)

#######
# %% Night Time Lights Covariate
#######
global cov_raster_filepath = COV_NTL_DIR
global cov_raster_filename = COV_NTL_FILENAME

# Load Rasters
global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

# Define storage variable
NTL_cov = zeros(size(hh_data_summary)[1])

# %% for each row in hh_data_summary
for row_i in 1:size(hh_data_summary)[1]
    # Get lat, lon valuyes
    lat, lon = hh_data_summary[row_i,["latitude", "longitude"]]

    # Find index of raster corresponding to latlon
    lats = lookup(cov_raster, Y)
    lons = lookup(cov_raster, X)

    lat_idx = argmin(abs.(lats .- lat))
    lon_idx = argmin(abs.(lons .- lon))

    # Extract required value
    NTL_cov[row_i] = Float64(cov_raster[lon_idx, lat_idx])
end

NTL_cov = max.(0, NTL_cov)

#######
# %% Elevation
#######
global cov_raster_filepath = COV_ELEV_DIR
global cov_raster_filename = COV_ELEV_FILENAME

# Load Rasters
global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

# Define storage variable
ELEV_cov = zeros(size(hh_data_summary)[1])

# %% for each row in hh_data_summary
for row_i in 1:size(hh_data_summary)[1]
    # Get lat, lon valuyes
    lat, lon = hh_data_summary[row_i,["latitude", "longitude"]]

    # Find index of raster corresponding to latlon
    lats = lookup(cov_raster, Y)
    lons = lookup(cov_raster, X)

    lat_idx = argmin(abs.(lats .- lat))
    lon_idx = argmin(abs.(lons .- lon))

    # Extract required value
    ELEV_cov[row_i] = Float64(cov_raster[lon_idx, lat_idx])
end

ELEV_cov = max.(0, ELEV_cov)

#######
# %% Slope
#######
global cov_raster_filepath = COV_SLP_DIR
global cov_raster_filename = COV_SLP_FILENAME

# Load Rasters
global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

# Define storage variable
SLP_cov = zeros(size(hh_data_summary)[1])

# %% for each row in hh_data_summary
for row_i in 1:size(hh_data_summary)[1]
    # Get lat, lon valuyes
    lat, lon = hh_data_summary[row_i,["latitude", "longitude"]]

    # Find index of raster corresponding to latlon
    lats = lookup(cov_raster, Y)
    lons = lookup(cov_raster, X)

    lat_idx = argmin(abs.(lats .- lat))
    lon_idx = argmin(abs.(lons .- lon))

    # Extract required value
    SLP_cov[row_i] = Float64(cov_raster[lon_idx, lat_idx])
end

SLP_cov = max.(0, SLP_cov)

################################################
# Monthly varying covariates
################################################

#######
# %% EVI
#######

# Define storage variable
EVI_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_EVI_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # # Temporary truncation of data (covariate data not avaialable after 2016)
    # if year > 2022
    #     year = 2022
    # end

    # Get list of required months
    months = unique(hh_data_summary[hh_data_summary.interview_year .== year,"interview_month"])

    # Select a month
    for month_idx in ProgressBar(1:length(months), leave = false)
        month = months[month_idx]

        # If required time period is not available, just take closest one
        if (year == 2000) && (month < 2)
            month = 2
        end

        # Define raster filename
        global cov_raster_filename = "NA"
        
        if month < 10
            global cov_raster_filename = "EVI_v061.$(year).0$(month).max.5km.max.tif"
        else
            global cov_raster_filename =  "EVI_v061.$(year).$(month).max.5km.max.tif"
        end

        # Load raster
        global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

        entry_idxs = findall((hh_data_summary.interview_year .== year) .& 
                                (hh_data_summary.interview_month .== month))

        temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

        for i in 1:length(entry_idxs)
            idx = entry_idxs[i]

            # Get lat, lon valuyes
            lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

            # Find index of raster corresponding to latlon
            lats = lookup(cov_raster, Y)
            lons = lookup(cov_raster, X)

            lat_idx = argmin(abs.(lats .- lat))
            lon_idx = argmin(abs.(lons .- lon))

            # Extract required value
            temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
        end

        # Save values to cov storage variable
        EVI_cov[entry_idxs] .= copy(temp_cov_vals)
    end
end
EVI_cov = max.(0, EVI_cov)

#######
# %% LST Day
#######

# Define storage variable
LSTD_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LSTD_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # # Temporary truncation of data (covariate data not avaialable after 2016)
    # if year > 2021
    #     year = 2021
    # end

    # Get list of required months
    months = unique(hh_data_summary[hh_data_summary.interview_year .== year,"interview_month"])

    # Select a month
    for month_idx in ProgressBar(1:length(months), leave = false)
        month = months[month_idx]

        # If required time period is not available, just take closest one
        if (year == 2000) && (month < 2)
            month = 2
        end

        # Define raster filename
        global cov_raster_filename = "NA"
        
        if month < 10
            global cov_raster_filename = "LST_Day_v061.$(year).0$(month).max.5km.max.tif"
        else
            global cov_raster_filename =  "LST_Day_v061.$(year).$(month).max.5km.max.tif"
        end

        # Load raster
        global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

        entry_idxs = findall((hh_data_summary.interview_year .== year) .& 
                                (hh_data_summary.interview_month .== month))

        temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

        for i in 1:length(entry_idxs)
            idx = entry_idxs[i]

            # Get lat, lon valuyes
            lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

            # Find index of raster corresponding to latlon
            lats = lookup(cov_raster, Y)
            lons = lookup(cov_raster, X)

            lat_idx = argmin(abs.(lats .- lat))
            lon_idx = argmin(abs.(lons .- lon))

            # Extract required value
            temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
        end

        # Save values to cov storage variable
        LSTD_cov[entry_idxs] .= copy(temp_cov_vals)
    end
end
LSTD_cov = max.(0, LSTD_cov)

#######
# %% LST Night
#######

# Define storage variable
LSTN_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LSTN_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # # Temporary truncation of data (covariate data not avaialable after 2016)
    # if year > 2021
    #     year = 2021
    # end

    # Get list of required months
    months = unique(hh_data_summary[hh_data_summary.interview_year .== year,"interview_month"])

    # Select a month
    for month_idx in ProgressBar(1:length(months), leave = false)
        month = months[month_idx]

        # If required time period is not available, just take closest one
        if (year == 2000) && (month < 2)
            month = 2
        end

        # Define raster filename
        global cov_raster_filename = "NA"
        
        if month < 10
            global cov_raster_filename = "LST_Night_v061.$(year).0$(month).max.5km.max.tif"
        else
            global cov_raster_filename =  "LST_Night_v061.$(year).$(month).max.5km.max.tif"
        end

        # Load raster
        global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

        entry_idxs = findall((hh_data_summary.interview_year .== year) .& 
                                (hh_data_summary.interview_month .== month))

        temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

        for i in 1:length(entry_idxs)
            idx = entry_idxs[i]

            # Get lat, lon valuyes
            lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

            # Find index of raster corresponding to latlon
            lats = lookup(cov_raster, Y)
            lons = lookup(cov_raster, X)

            lat_idx = argmin(abs.(lats .- lat))
            lon_idx = argmin(abs.(lons .- lon))

            # Extract required value
            temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
        end

        # Save values to cov storage variable
        LSTN_cov[entry_idxs] .= copy(temp_cov_vals)
    end
end
LSTN_cov = max.(0, LSTN_cov)

#######
# %% LST DELTA
#######

# Define storage variable
LSTDELTA_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LSTDELTA_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Get list of required months
    months = unique(hh_data_summary[hh_data_summary.interview_year .== year,"interview_month"])

    # Select a month
    for month_idx in ProgressBar(1:length(months), leave = false)
        month = months[month_idx]

        # If required time period is not available, just take closest one
        if (year == 2000) && (month < 2)
            month = 2
        end

        # Define raster filename
        global cov_raster_filename = "NA"
        
        if month < 10
            global cov_raster_filename = "LST_DiurnalDifference.$(year).0$(month).max.5km.max.tif"
        else
            global cov_raster_filename =  "LST_DiurnalDifference.$(year).$(month).max.5km.max.tif"
        end

        # Load raster
        global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

        entry_idxs = findall((hh_data_summary.interview_year .== year) .& 
                                (hh_data_summary.interview_month .== month))

        temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

        for i in 1:length(entry_idxs)
            idx = entry_idxs[i]

            # Get lat, lon valuyes
            lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

            # Find index of raster corresponding to latlon
            lats = lookup(cov_raster, Y)
            lons = lookup(cov_raster, X)

            lat_idx = argmin(abs.(lats .- lat))
            lon_idx = argmin(abs.(lons .- lon))

            # Extract required value
            temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
        end

        # Save values to cov storage variable
        LSTDELTA_cov[entry_idxs] .= copy(temp_cov_vals)
    end
end
LSTDELTA_cov = max.(0, LSTDELTA_cov)

#######
# %% TCW
#######

# Define storage variable
TCW_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_TCW_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Get list of required months
    months = unique(hh_data_summary[hh_data_summary.interview_year .== year,"interview_month"])

    # Select a month
    for month_idx in ProgressBar(1:length(months), leave = false)
        month = months[month_idx]

        # If required time period is not available, just take closest one
        if (year == 2000) && (month < 2)
            month = 2
        end

        # Define raster filename
        global cov_raster_filename = "NA"
        
        if month < 10
            global cov_raster_filename = "TCW_v061.$(year).0$(month).max.5km.max.tif"
        else
            global cov_raster_filename =  "TCW_v061.$(year).$(month).max.5km.max.tif"
        end

        # Load raster
        global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

        entry_idxs = findall((hh_data_summary.interview_year .== year) .& 
                                (hh_data_summary.interview_month .== month))

        temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

        for i in 1:length(entry_idxs)
            idx = entry_idxs[i]

            # Get lat, lon valuyes
            lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

            # Find index of raster corresponding to latlon
            lats = lookup(cov_raster, Y)
            lons = lookup(cov_raster, X)

            lat_idx = argmin(abs.(lats .- lat))
            lon_idx = argmin(abs.(lons .- lon))

            # Extract required value
            temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
        end

        # Save values to cov storage variable
        TCW_cov[entry_idxs] .= copy(temp_cov_vals)
    end
end
TCW_cov = max.(0, TCW_cov)

#######
# %% TSI
#######

# Define storage variable
TSI_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_TSI_DIR

# Select a year
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # # Temporary truncation of data (covariate data not avaialable after 2016)
    # if year > 2022
    #     year = 2022
    # end

    # Get list of required months
    months = unique(hh_data_summary[hh_data_summary.interview_year .== year,"interview_month"])

    # Select a month
    for month_idx in ProgressBar(1:length(months), leave = false)
        month = months[month_idx]

        # If required time period is not available, just take closest one
        if (year == 2000) && (month < 4)
            month = 4
        end

        # Define raster filename
        global cov_raster_filename = "NA"
        
        if month < 10
            global cov_raster_filename = "TSI-Martens2-Pf.$(year).0$(month).Data.5km.Data.tif"
        else
            global cov_raster_filename =  "TSI-Martens2-Pf.$(year).$(month).Data.5km.Data.tif"
        end

        # Load raster
        global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

        entry_idxs = findall((hh_data_summary.interview_year .== year) .& 
                                (hh_data_summary.interview_month .== month))

        temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

        for i in 1:length(entry_idxs)
            idx = entry_idxs[i]

            # Get lat, lon valuyes
            lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

            # Find index of raster corresponding to latlon
            lats = lookup(cov_raster, Y)
            lons = lookup(cov_raster, X)

            lat_idx = argmin(abs.(lats .- lat))
            lon_idx = argmin(abs.(lons .- lon))

            # Extract required value
            temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
        end

        # Save values to cov storage variable
        TSI_cov[entry_idxs] .= copy(temp_cov_vals)
    end
end
TSI_cov = max.(0, TSI_cov)

#######
# %% TCB
#######

# Define storage variable
TCB_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_TCB_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # # Temporary truncation of data (covariate data not avaialable after 2016)
    # if year > 2022
    #     year = 2022
    # end

    # Get list of required months
    months = unique(hh_data_summary[hh_data_summary.interview_year .== year,"interview_month"])

    # Select a month
    for month_idx in ProgressBar(1:length(months), leave = false)
        month = months[month_idx]

        # If required time period is not available, just take closest one
        if (year == 2000) && (month < 2)
            month = 2
        end

        # Define raster filename
        global cov_raster_filename = "NA"
        
        if month < 10
            global cov_raster_filename = "TCB_v061.$(year).0$(month).max.5km.max.tif"
        else
            global cov_raster_filename =  "TCB_v061.$(year).$(month).max.5km.max.tif"
        end

        # Load raster
        global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

        entry_idxs = findall((hh_data_summary.interview_year .== year) .& 
                                (hh_data_summary.interview_month .== month))

        temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

        for i in 1:length(entry_idxs)
            idx = entry_idxs[i]

            # Get lat, lon valuyes
            lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

            # Find index of raster corresponding to latlon
            lats = lookup(cov_raster, Y)
            lons = lookup(cov_raster, X)

            lat_idx = argmin(abs.(lats .- lat))
            lon_idx = argmin(abs.(lons .- lon))

            # Extract required value
            temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
        end

        # Save values to cov storage variable
        TCB_cov[entry_idxs] .= copy(temp_cov_vals)
    end
end
TCB_cov = max.(0, TCB_cov)

################################################
# Annual varying covariates
################################################

#######
# %% Landcover 00
#######

# Define storage variable
LAND00_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND00_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-00_Unclassified.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND00_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND00_cov = max.(0, LAND00_cov)

#######
# %% Landcover 01
#######

# Define storage variable
LAND01_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND01_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-01_Evergreen_Needleleaf_Forest.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND01_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND01_cov = max.(0, LAND01_cov)

#######
# %% Landcover 02
#######

# Define storage variable
LAND02_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND02_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-02_Evergreen_Broadleaf_Forest.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND02_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND02_cov = max.(0, LAND02_cov)

#######
# %% Landcover 03
#######

# Define storage variable
LAND03_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND03_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-03_Deciduous_Needleleaf_Forest.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND03_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND03_cov = max.(0, LAND03_cov)

#######
# %% Landcover 04
#######

# Define storage variable
LAND04_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND04_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-04_Deciduous_Broadleaf_Forest.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND04_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND04_cov = max.(0, LAND04_cov)


#######
# %% Landcover 05
#######

# Define storage variable
LAND05_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND05_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-05_Mixed_Forest.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND05_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND05_cov = max.(0, LAND05_cov)

#######
# %% Landcover 06
#######

# Define storage variable
LAND06_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND06_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-06_Closed_Shrublands.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND06_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND06_cov = max.(0, LAND06_cov)

#######
# %% Landcover 07
#######

# Define storage variable
LAND07_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND07_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-07_Open_Shrublands.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND07_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND07_cov = max.(0, LAND07_cov)

#######
# %% Landcover 08
#######

# Define storage variable
LAND08_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND08_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-08_Woody_Savannas.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND08_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND08_cov = max.(0, LAND08_cov)

#######
# %% Landcover 09
#######

# Define storage variable
LAND09_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND09_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-09_Savannas.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND09_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND09_cov = max.(0, LAND09_cov)

#######
# %% Landcover 10
#######

# Define storage variable
LAND10_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND10_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-10_Grasslands.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND10_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND10_cov = max.(0, LAND10_cov)

#######
# %% Landcover 11
#######

# Define storage variable
LAND11_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND11_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-11_Permanent_Wetlands.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND11_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND11_cov = max.(0, LAND11_cov)

#######
# %% Landcover 12
#######

# Define storage variable
LAND12_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND12_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-12_Croplands.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND12_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND12_cov = max.(0, LAND12_cov)

#######
# %% Landcover 13
#######

# Define storage variable
LAND13_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND13_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-13_Urban_And_Built_Up.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND13_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND13_cov = max.(0, LAND13_cov)

#######
# %% Landcover 14
#######

# Define storage variable
LAND14_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND14_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-14_Cropland_Natural_Vegetation_Mosaic.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND14_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND14_cov = max.(0, LAND14_cov)

#######
# %% Landcover 15
#######

# Define storage variable
LAND15_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND15_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-15_Snow_And_Ice.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND15_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND15_cov = max.(0, LAND15_cov)

#######
# %% Landcover 16
#######

# Define storage variable
LAND16_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND16_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-16_Barren_Or_Sparsely_Populated.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND16_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND16_cov = max.(0, LAND16_cov)

#######
# %% Landcover 17
#######

# Define storage variable
LAND17_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND17_DIR

# Select ayear
for year_idx in ProgressBar(1:length(YEAR_VALS))
    year = YEAR_VALS[year_idx]

    # Temporary truncation of data (covariate data not avaialable after 2016)
    if year < 2001
        year = 2001
    elseif year > 2023
        year = 2023
    end

    # Define raster filename
    global cov_raster_filename = "IGBP_Landcover_Class-17_Water.$(year).Annual.Data.5km.fraction.tif"
        
    # Load raster
    global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

    entry_idxs = findall((hh_data_summary.interview_year .== year))

    temp_cov_vals = Vector{Float64}(undef, length(entry_idxs))

    for i in 1:length(entry_idxs)
        idx = entry_idxs[i]

        # Get lat, lon valuyes
        lat, lon = hh_data_summary[idx,["latitude", "longitude"]]

        # Find index of raster corresponding to latlon
        lats = lookup(cov_raster, Y)
        lons = lookup(cov_raster, X)

        lat_idx = argmin(abs.(lats .- lat))
        lon_idx = argmin(abs.(lons .- lon))

        # Extract required value
        temp_cov_vals[i] = Float64(cov_raster.data[lon_idx,lat_idx])
    end

    # Save values to cov storage variable
    LAND17_cov[entry_idxs] .= copy(temp_cov_vals)

end
LAND17_cov = max.(0, LAND17_cov)

#######
# %% Temperature suitability covariate
#######
monthidxs = monthyear_to_monthidx.(hh_data_summary.interview_month,hh_data_summary.interview_year, YEAR_START = 2000)


####################################
# %% Concatenate Processed Data with Covarariates
####################################

INLA_dataset = hcat(hh_data_summary, DataFrame(monthidx = monthidxs, 
                                                cov_ACCESS = ACCESS_cov,
                                                cov_PET = PET_cov,
                                                cov_ARID = ARID_cov,
                                                cov_NTL = NTL_cov,
                                                cov_ELEV = ELEV_cov,
                                                cov_SLP = SLP_cov,
                                                cov_EVI = EVI_cov,
                                                cov_LSTD = LSTD_cov,
                                                cov_LSTN = LSTN_cov,
                                                cov_LSTDELTA = LSTDELTA_cov,
                                                # cov_TCW = TCW_cov,
                                                cov_TSI = TSI_cov,
                                                cov_TCB = TCB_cov,
                                                # cov_LAND00 = LAND00_cov,
                                                cov_LAND01 = LAND01_cov,
                                                cov_LAND02 = LAND02_cov,
                                                # cov_LAND03 = LAND03_cov,
                                                cov_LAND04 = LAND04_cov,
                                                cov_LAND05 = LAND05_cov,
                                                cov_LAND06 = LAND06_cov,
                                                cov_LAND07 = LAND07_cov,
                                                cov_LAND08 = LAND08_cov,
                                                cov_LAND09 = LAND09_cov,
                                                cov_LAND10 = LAND10_cov,
                                                cov_LAND11 = LAND11_cov,
                                                cov_LAND12 = LAND12_cov,
                                                cov_LAND13 = LAND13_cov,
                                                cov_LAND14 = LAND14_cov,
                                                # cov_LAND15 = LAND15_cov,
                                                cov_LAND16 = LAND16_cov,
                                                cov_LAND17 = LAND17_cov))

#######
# %% Find list of all entries with missing covariate data and remove/filter out
#######
cov_values = Matrix(INLA_dataset[:,15:end])
missing_cov_idxs = findall(sum(cov_values .== -9999, dims = 2)[:].>0)
valid_entries_idx = setdiff(1:size(INLA_dataset)[1], missing_cov_idxs)
filt_INLA_dataset = INLA_dataset[valid_entries_idx,:]

#######
# %% Save dataset
#######
CSV.write(output_dir*output_filename, filt_INLA_dataset)