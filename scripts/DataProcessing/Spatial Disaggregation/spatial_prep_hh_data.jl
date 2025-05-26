"""
Author: Eugene Tan
Date Created: 2/12/2024
Last Updated: 20/5/2024
Script to go through list of all household survey entries and extract relevant covariate values at each survey
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

# using Shapefile
using LinearAlgebra
using GeoInterface
using GeoIO
using JLD2
using StatsBase

# Custom packages
using DateConversions

# %% File paths
# Output save dir
output_dir = OUTPUT_DATAPREP_DIR*"INLA/"
output_filename = INLA_DATAPREP_FILENAME

# Household Survey Data
hh_dir = OUTPUT_DATAPREP_DIR#"datasets/subnational/"
hh_filename = HOUSEHOLD_SURVEY_DATA_FILENAME

# MITN Posterior Estimates
snf_post_dir = OUTPUT_DRAWS_DIR*"subnational/"

# Region Admin 1 area id legend
admin1_legend_dir = RAW_SUBNAT_DATASET_DIR
admin1_legend_filename = ADMIN1_AREAID_LEGEND_FILENAME

# Region boundaries
admin1_shapes_geoIO = GeoIO.load(ADMIN1_SHAPEFILE)

# Covariates
####################################
# %% Pre-process and Summarise Household Data
####################################
# Load admin1 area id legend
admin1_legend =  CSV.read(admin1_legend_dir*admin1_legend_filename, DataFrame)

# Load raw survey data
full_hh_data = CSV.read(hh_dir*hh_filename, DataFrame)

# filter for only survey entries where there is lat-lon info
hh_data = full_hh_data[.!ismissing.(full_hh_data.latitude),:]

# Get list of unique ISOs
ISO_list = ISO_LIST
exclusion_ISOs = ["CPV","ZAF"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# Storage variable for collection of dataframe rows
df_collection = []

# Go through each ISO and calculate/extract required data
for ISO_i in 1:length(filt_ISOs)
    
    ISO = filt_ISOs[ISO_i]

    println("Extracting houshold survey data for country $(ISO_i) of $(length(filt_ISOs)) → $(ISO)")

    #####################################
    # Load MITN Posteriors data
    #####################################
    snf_post_filename = "$(ISO)_SUBNAT_draws.jld2"
    snf_post_draws = JLD2.load(snf_post_dir*snf_post_filename)

    # Extract relevant household survey data based on current ISO    
    country_hh_data = hh_data[(hh_data.ISO .== ISO) .&
                            (.!ismissing.(hh_data.area_id)) .&
                            (hh_data.interview_year .>= snf_post_draws["YEAR_START_NAT"]),:]
    # Compress raw survey data into weighted average estimates at unique latlon entries
    sid_lat_lons = unique([(country_hh_data.SurveyId[i], 
                        country_hh_data.latitude[i], 
                        country_hh_data.longitude[i],
                        country_hh_data.interview_month[i],
                        country_hh_data.interview_year[i]) for i in 1:size(country_hh_data)[1]])

    # Define storage variables
    n_unique_obs = length(sid_lat_lons)
    surveyids = Vector{String}(undef, n_unique_obs)
    ISOs = Vector{String}(undef, n_unique_obs)
    admin1_names = Vector{String}(undef, n_unique_obs)
    area_ids = Vector{Int64}(undef, n_unique_obs)
    latitude = Vector{Float64}(undef, n_unique_obs)
    longitude = Vector{Float64}(undef, n_unique_obs)
    interview_month = Vector{Int64}(undef, n_unique_obs)
    interview_year = Vector{Int64}(undef, n_unique_obs)
    cluster_sample_wt = Vector{Float64}(undef, n_unique_obs)
    npc_vals = Vector{Float64}(undef, n_unique_obs)
    access_vals = Vector{Float64}(undef, n_unique_obs)
    use_vals = Vector{Float64}(undef, n_unique_obs)
    npc_gap_vals = Vector{Float64}(undef, n_unique_obs)
    access_gap_vals = Vector{Float64}(undef, n_unique_obs)
    use_gap_vals = Vector{Float64}(undef, n_unique_obs)

    Threads.@threads for i in 1:n_unique_obs
        surveyid, lat, lon, mth, year = sid_lat_lons[i]

        filt_ISOs = country_hh_data[(country_hh_data.SurveyId .== surveyid) .& (country_hh_data.latitude .== lat) .& (country_hh_data.longitude .== lon),:]

        #####################################
        # Extract required metadata
        #####################################
        area_id = filt_ISOs[1,"area_id"]
        admin1_name = admin1_legend[findfirst(admin1_legend.area_id .== area_id),"Name_1"]
        snf_dict = snf_post_draws["merged_outputs"][findfirst(snf_post_draws["admin1_names"] .== admin1_name)]

        # Get required month index
        monthidx = monthyear_to_monthidx(mth, year, YEAR_START = snf_post_draws["YEAR_START_NAT"])

        # Get MITN snf estimates of admin1 region NPC and access
        snf_npc = mean(snf_dict["ADJ_NPC_MONTHLY_TOTAL_samples"][:,monthidx])
        snf_access = mean(snf_dict["ADJ_λ_ACCESS_samples"][:,monthidx])

        #####################################
        # Calculate survey estimates of NPC, access, use and cluster weight
        #####################################
        npc_entries = filt_ISOs.n_itn./filt_ISOs.hh_size
        access_entries = min.(0.5, npc_entries)./0.5
        # use_entries = filt_ISOs.n_itn_used./filt_ISOs.hh_size
        weights = filt_ISOs.hh_sample_wt
        npc = sum(filt_ISOs.n_itn.*weights)/sum(filt_ISOs.hh_size.*weights) #sum(weights.*npc_entries)/sum(weights)
        access = sum(min.((2 .*filt_ISOs.n_itn)./(filt_ISOs.hh_size),1).*(filt_ISOs.hh_size).*weights)/sum((filt_ISOs.hh_size).*weights) #sum(weights.*access_entries)/sum(weights)
        use = sum(min.(2 .*filt_ISOs.n_itn_used ./ filt_ISOs.hh_size, 1).*(filt_ISOs.hh_size).*weights)/sum(weights .* filt_ISOs.hh_size) #sum(weights.*use_entries)/sum(weights)
        
        # use = dot(weights, filt_ISOs.n_slept_under_itn)/dot(weights, filt_ISOs.hh_size)
        sample_weight = sum(weights)

        #####################################
        # Calculate NPC and access gap
        #####################################
        npc_gap = npc-snf_npc
        access_gap = access-snf_access
        use_gap = access-use

        #####################################
        # Save results to storage variable
        #####################################
        surveyids[i] = surveyid
        ISOs[i] = ISO
        admin1_names[i] = admin1_name
        area_ids[i] = area_id
        latitude[i] = lat
        longitude[i] = lon
        interview_month[i] = mth
        interview_year[i] = year
        cluster_sample_wt[i] = sample_weight
        npc_vals[i] = npc
        access_vals[i] = access
        use_vals[i] = use
        npc_gap_vals[i] = npc_gap
        access_gap_vals[i] = access_gap
        use_gap_vals[i] = use_gap
    end

    # Construct Data Frame
    df = DataFrame(ISO = ISOs, admin1_name = admin1_names, area_id = area_ids,
                        latitude = latitude, longitude = longitude,
                        interview_month = interview_month, interview_year = interview_year,
                        cluster_sample_wt = cluster_sample_wt,
                        npc = npc_vals, access = access_vals, use = use_vals, 
                        npc_gap = npc_gap_vals, access_gap = access_gap_vals,
                        use_gap = use_gap_vals)

    push!(df_collection, df)
end

GC.gc()

# %% Concatenate all dataframe fragments
hh_data_summary = vcat(df_collection...)

####################################
# %% Get Year bounds to optimise data extraction time
####################################
YEAR_VALS = sort(unique(hh_data_summary.interview_year))

####################################
# %% Get Covariates from Raster at each survey value
####################################

#######
# %% Accessibility to Cities
#######

println("Extracting static covariate: Accessibility")
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
Threads.@threads for row_i in 1:size(hh_data_summary)[1]
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
GC.gc()

#######
# %% PET
#######
println("Extracting static covariate: PET")

global cov_raster_filepath = COV_PET_DIR
global cov_raster_filename = COV_PET_FILENAME

# Load Rasters
global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

# Define storage variable
PET_cov = zeros(size(hh_data_summary)[1])

# %% for each row in hh_data_summary
Threads.@threads for row_i in 1:size(hh_data_summary)[1]
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
GC.gc()

#######
# %% Aridity Covariate
#######
println("Extracting static covariate: ARID")

global cov_raster_filepath = COV_ARID_DIR
global cov_raster_filename = COV_ARID_FILENAME

# Load Rasters
global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

# Define storage variable
ARID_cov = zeros(size(hh_data_summary)[1])

# %% for each row in hh_data_summary
Threads.@threads for row_i in 1:size(hh_data_summary)[1]
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
GC.gc()

#######
# %% Night Time Lights Covariate
#######
println("Extracting static covariate: NTL")

global cov_raster_filepath = COV_NTL_DIR
global cov_raster_filename = COV_NTL_FILENAME

# Load Rasters
global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

# Define storage variable
NTL_cov = zeros(size(hh_data_summary)[1])

# %% for each row in hh_data_summary
Threads.@threads for row_i in 1:size(hh_data_summary)[1]
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
GC.gc()

#######
# %% Elevation
#######
println("Extracting static covariate: Elevation")

global cov_raster_filepath = COV_ELEV_DIR
global cov_raster_filename = COV_ELEV_FILENAME

# Load Rasters
global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

# Define storage variable
ELEV_cov = zeros(size(hh_data_summary)[1])

# %% for each row in hh_data_summary
Threads.@threads for row_i in 1:size(hh_data_summary)[1]
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
GC.gc()

#######
# %% Slope
#######
println("Extracting static covariate: Slope")

global cov_raster_filepath = COV_SLP_DIR
global cov_raster_filename = COV_SLP_FILENAME

# Load Rasters
global cov_raster = replace_missing(Raster(cov_raster_filepath*cov_raster_filename), missingval = -9999)

# Define storage variable
SLP_cov = zeros(size(hh_data_summary)[1])

# %% for each row in hh_data_summary
Threads.@threads for row_i in 1:size(hh_data_summary)[1]
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
GC.gc()

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
for year_idx in 1:length(YEAR_VALS)
    
    year = YEAR_VALS[year_idx]

    println("Extracting monthly covariate: EVI, Year $(year)")

    # # Temporary truncation of data (covariate data not avaialable after 2016)
    # if year > 2022
    #     year = 2022
    # end

    # Get list of required months
    months = unique(hh_data_summary[hh_data_summary.interview_year .== year,"interview_month"])

    # Select a month
    for month_idx in 1:length(months)
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

        Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% LST Day
#######

# Define storage variable
LSTD_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LSTD_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting monthly covariate: LSTD, Year $(year)")

    # # Temporary truncation of data (covariate data not avaialable after 2016)
    # if year > 2021
    #     year = 2021
    # end

    # Get list of required months
    months = unique(hh_data_summary[hh_data_summary.interview_year .== year,"interview_month"])

    # Select a month
    for month_idx in 1:length(months)
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

        Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% LST Night
#######

# Define storage variable
LSTN_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LSTN_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting monthly covariate: LSTN, Year $(year)")

    # # Temporary truncation of data (covariate data not avaialable after 2016)
    # if year > 2021
    #     year = 2021
    # end

    # Get list of required months
    months = unique(hh_data_summary[hh_data_summary.interview_year .== year,"interview_month"])

    # Select a month
    for month_idx in 1:length(months)
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

        Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% LST DELTA
#######

# Define storage variable
LSTDELTA_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LSTDELTA_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting monthly covariate: LSTDELTA, Year $(year)")

    # Get list of required months
    months = unique(hh_data_summary[hh_data_summary.interview_year .== year,"interview_month"])

    # Select a month
    for month_idx in 1:length(months)
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

        Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% TCW
#######

# Define storage variable
TCW_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_TCW_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting monthly covariate: TCW, Year $(year)")

    # Get list of required months
    months = unique(hh_data_summary[hh_data_summary.interview_year .== year,"interview_month"])

    # Select a month
    for month_idx in 1:length(months)
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

        Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% TSI
#######

# Define storage variable
TSI_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_TSI_DIR

# Select a year
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting monthly covariate: TSI, Year $(year)")

    # # Temporary truncation of data (covariate data not avaialable after 2016)
    # if year > 2022
    #     year = 2022
    # end

    # Get list of required months
    months = unique(hh_data_summary[hh_data_summary.interview_year .== year,"interview_month"])

    # Select a month
    for month_idx in 1:length(months)
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

        Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% TCB
#######

# Define storage variable
TCB_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_TCB_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting monthly covariate: TCB, Year $(year)")

    # # Temporary truncation of data (covariate data not avaialable after 2016)
    # if year > 2022
    #     year = 2022
    # end

    # Get list of required months
    months = unique(hh_data_summary[hh_data_summary.interview_year .== year,"interview_month"])

    # Select a month
    for month_idx in 1:length(months)
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

        Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

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
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND00, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 01
#######

# Define storage variable
LAND01_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND01_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND01, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 02
#######

# Define storage variable
LAND02_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND02_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND02, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 03
#######

# Define storage variable
LAND03_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND03_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND03, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 04
#######

# Define storage variable
LAND04_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND04_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND04, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()


#######
# %% Landcover 05
#######

# Define storage variable
LAND05_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND05_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND05, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 06
#######

# Define storage variable
LAND06_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND06_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND06, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 07
#######

# Define storage variable
LAND07_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND07_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND07, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 08
#######

# Define storage variable
LAND08_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND08_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND08, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 09
#######

# Define storage variable
LAND09_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND09_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND09, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 10
#######

# Define storage variable
LAND10_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND10_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND10, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 11
#######

# Define storage variable
LAND11_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND11_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND11, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 12
#######

# Define storage variable
LAND12_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND12_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND12, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 13
#######

# Define storage variable
LAND13_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND13_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND13, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 14
#######

# Define storage variable
LAND14_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND14_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND14, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 15
#######

# Define storage variable
LAND15_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND15_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND15, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 16
#######

# Define storage variable
LAND16_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND16_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND16, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

#######
# %% Landcover 17
#######

# Define storage variable
LAND17_cov = zeros(size(hh_data_summary)[1])

# Define raster filepath
global cov_raster_filepath = COV_LAND17_DIR

# Select ayear
for year_idx in 1:length(YEAR_VALS)
    year = YEAR_VALS[year_idx]

    println("Extracting annual covariate: LAND17, Year $(year)")

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

    Threads.@threads for i in 1:length(entry_idxs)
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
GC.gc()

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
# %% Save dataset
#######
mkpath(output_dir)
CSV.write(output_dir*"unfiltered_inla_dataset.csv", INLA_dataset)

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

