"""
Author: Eugene Tan
Date Created: 6/4/2025
Last Updated: 6/4/2025
Calculate prediction error of the entire summarised survey dataset w.r.t to the final regressed rasters
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

# %% Dataset Directories
raster_dir = OUTPUT_RASTERS_DIR
dataprep_dir = OUTPUT_DATAPREP_DIR*"INLA/"
bv_output_dir = BV_OUTPUT_DIR
pop_rasters_dir = POPULATION_RASTER_DIR

# %% Filenames
inla_data_filename = "inla_dataset_reduced.csv"
country_summaries_filename = HOUSEHOLD_NAT_SUMMARY_DATA_FILENAME

# %% Load Datasets
inla_data = CSV.read(dataprep_dir*inla_data, DataFrame)
country_summaries = CSV.read(OUTPUT_DATAPREP_DIR*country_summaries_filename, DataFrame)
hh_survey_data = CSV.read(OUTPUT_DATAPREP_DIR*HOUSEHOLD_SURVEY_DATA_FILENAME, DataFrame)

# Region boundaries
admin0_shapes_geoIO = GeoIO.load(ADMIN0_SHAPEFILE)
################################
# %% LOCAL LEVEL
################################
# %% Compile new Dataframe with comparisons
df_local_collection = []

# Get list of unique years
YEAR_LIST = unique(inla_data.interview_year)

# Extract model estimates at local cluster level and compile dataframe
for year in YEAR_LIST
    # Select Year
    year = YEAR_LIST[1]

    # Import annual resolution final rasters (i.e. BV model) and get require lat lon collection
    bv_npc_raster = replace_missing(Raster(BV_OUTPUT_DIR*"nets_per_capita/ITN_$(year)_percapita_nets_mean.tif"), missing_val = NaN)
    bv_use_raster = replace_missing(Raster(BV_OUTPUT_DIR*"ITN_$(year)_use_mean.tif"), missingval = NaN)
    bv_model_lats = lookup(bv_npc_raster, Y)
    bv_model_lons = lookup(bv_npc_raster, X)

    # Filter survey data and get list of unique month values
    data_yearslice = inla_data[inla_data.interview_year .== year,:]
    monthidxs = unique(data_yearslice.interview_month)

    for month in monthidxs

        data_monthslice = data_yearslice[data_yearslice.interview_month .== month,:]

        lats = data_monthslice.latitude
        lons = data_monthslice.longitude

        # Scrape BV local prediction values
        bv_npc_vals = zeros(size(data_monthslice)[1])
        bv_use_vals = zeros(size(data_monthslice)[1])

        for j in 1:size(data_monthslice)[1]
            lat = lats[j]
            lon = lons[j]

            # Lookup raster idx location
            bv_model_lat_idx = argmin(abs.(bv_model_lats .- lat))
            bv_model_lon_idx = argmin(abs.(bv_model_lons .- lon))

            # Use raster idx to lookup value
            bv_npc_vals[j] = bv_npc_raster[bv_model_lon_idx, bv_model_lat_idx]
            bv_use_vals[j] = bv_use_raster[bv_model_lon_idx, bv_model_lat_idx]
        end

        # Scrape MITN local prediction values
        month_str = "$(month)"
        if month < 10
            month_str = "0$(month)"
        end

        # Import MITN rasters
        mitn_npc_raster = replace_missing(Raster(raster_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        mitn_adj_npc_raster = replace_missing(Raster(raster_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        # mitn_final_npc_raster = replace_missing(Raster(raster_dir*"final_npc/logmodel_npc/final_npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)

        mitn_access_raster = replace_missing(Raster(raster_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        mitn_adj_access_raster = replace_missing(Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        # mitn_final_access_raster = replace_missing(Raster(raster_dir*"final_access/pmodel_access/final_access_$(year)_$(month_str)_mean.tif"), missingval = NaN)

        mitn_use_raster = replace_missing(Raster(raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        mitn_final_use_raster = replace_missing(Raster(raster_dir*"final_use/logis_use/final_use_$(year)_$(month_str)_mean.tif"), missingval = NaN)

        mitn_model_lats = lookup(mitn_npc_raster, Y)
        mitn_model_lons = lookup(mitn_npc_raster, X)

        mitn_npc_vals = zeros(size(data_monthslice)[1])
        mitn_access_vals = zeros(size(data_monthslice)[1])
        mitn_use_vals = zeros(size(data_monthslice)[1])

        mitn_adj_npc_vals = zeros(size(data_monthslice)[1])
        mitn_adj_access_vals = zeros(size(data_monthslice)[1])

        mitn_final_npc_vals = zeros(size(data_monthslice)[1])
        mitn_final_access_vals = zeros(size(data_monthslice)[1])
        mitn_final_use_vals = zeros(size(data_monthslice)[1])

        for j in 1:size(data_monthslice)[1]
            lat = lats[j]
            lon = lons[j]

            # Lookup raster idx location
            mitn_model_lat_idx = argmin(abs.(mitn_model_lats .- lat))
            mitn_model_lon_idx = argmin(abs.(mitn_model_lons .- lon))

            # Use raster idx to lookup value
            mitn_npc_vals[j] = mitn_npc_raster[mitn_model_lon_idx, mitn_model_lat_idx]
            mitn_access_vals[j] = mitn_access_raster[mitn_model_lon_idx, mitn_model_lat_idx]
            # mitn_use_vals[j] = mitn_use_raster[mitn_model_lon_idx, mitn_model_lat_idx]

            mitn_adj_npc_vals[j] = mitn_adj_npc_raster[mitn_model_lon_idx, mitn_model_lat_idx]
            mitn_adj_access_vals[j] = mitn_adj_access_raster[mitn_model_lon_idx, mitn_model_lat_idx]

            # mitn_final_npc_vals[j] = mitn_final_npc_raster[mitn_model_lon_idx, mitn_model_lat_idx]
            # mitn_final_access_vals[j] = mitn_final_access_raster[mitn_model_lon_idx, mitn_model_lat_idx]
            # mitn_final_use_vals[j] = mitn_final_use_raster[mitn_model_lon_idx, mitn_model_lat_idx]
        end

        # Construct Dataframe for year-month survey entries
        df = DataFrame(Type = "local",
                    ISO = data_monthslice.ISO,
                    admin1_name = data_monthslice.admin1_name,
                    area_id = data_monthslice.area_id,
                    latitude = data_monthslice.latitude,
                    longitude = data_monthslice.longitude,
                    interview_month = data_monthslice.interview_month,
                    interview_year = data_monthslice.interview_year,
                    cluster_sample_wt = data_monthslice.cluster_sample_wt,
                    npc = data_monthslice.npc,
                    bv_npc = bv_npc_vals,
                    mitn_npc = mitn_npc_vals,
                    mitn_adj_npc = mitn_adj_npc_vals,
                    mitn_final_npc = mitn_final_npc_vals,
                    access = data_monthslice.access,
                    bv_access = NaN,
                    mitn_access = mitn_access_vals,
                    mitn_adj_access = mitn_adj_access_vals,
                    mitn_final_access = mitn_final_access_vals,
                    use = data_monthslice.use,
                    bv_use = bv_use_vals,
                    mitn_use = mitn_use_vals,
                    mitn_final_use = mitn_final_use_vals)
        push!(df_local_collection, df)
    end
end

################################
# %% COUNTRY LEVEL
################################


# %%
admin0_df_entries = []
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = ["CPV", "ZAF"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

using NetAccessModel
using DateConversions

for ISO_i in ProgressBar(1:length(filt_ISOs))
    ISO = filt_ISOs[ISO_i]
    # Import access data and calculate household access
    net_access_input_dict = load(OUTPUT_EXTRACTIONS_DIR*"access/reg_data/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")

    # Calculate reference values for National household Access from Surveys
    national_H_aggregated = net_access_input_dict["H_aggregated"]
    national_access_aggregated = zeros(size(national_H_aggregated)[1])
    for survey_i in 1:length(national_access_aggregated)
        national_access_aggregated[survey_i] = sum(H_to_access(national_H_aggregated[survey_i,:,:]))
    end
    

    # Get list of monthidx to extract values from (those present in surveys)
    nat_monthidxs = monthyear_to_monthidx.(country_summaries[(country_summaries.ISO .== ISO) .& 
                                                            (country_summaries.year .>= YEAR_START) .&
                                                            (country_summaries.year .<= YEAR_END),:].month, 
                                                            country_summaries[(country_summaries.ISO .== ISO) .& 
                                                    (country_summaries.year .>= YEAR_START) .&
                                                    (country_summaries.year .<= YEAR_END),:].year, YEAR_START = YEAR_START)
    access_monthidxs = net_access_input_dict["survey_monthidx_aggregated"]
    use_monthidxs = monthyear_to_monthidx.(country_summaries[(country_summaries.ISO .== ISO) .& 
                                                            (country_summaries.year .>= YEAR_START) .&
                                                            (country_summaries.year .<= YEAR_END),:].month, 
                                                            country_summaries[(country_summaries.ISO .== ISO) .& 
                                                            (country_summaries.year .>= YEAR_START) .&
                                                            (country_summaries.year .<= YEAR_END),:].year, YEAR_START = YEAR_START)

    monthidxs = union(nat_monthidxs, use_monthidxs, access_monthidxs)

    for i in ProgressBar(1:length(monthidxs))
        monthidx = monthidxs[i]

        # Get corresponding month and year
        month = monthidx_to_monthyear(monthidx)[1]
        year = (monthidx_to_monthyear(monthidx)[2]-1)+YEAR_START

        month_str = "$(month)"
        if month < 10
            month_str = "0"*month_str
        end

        # Import population raster
        population_raster = replace_missing(Raster(pop_rasters_dir*"WorldPop_UNAdj_v3_DRC_fix.$(min(year,2020)).Annual.Data.5km.sum.tif"), missingval = NaN)

        # Import BV model rasters
        bv_npc_raster = replace_missing(Raster(BV_OUTPUT_DIR*"nets_per_capita/ITN_$(year)_percapita_nets_mean.tif"), missingval = NaN)
        bv_use_raster = replace_missing(Raster(BV_OUTPUT_DIR*"ITN_$(year)_use_mean.tif"), missingval = NaN)

        # Import MITN model rasters
        mitn_npc_raster = replace_missing(Raster(raster_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        mitn_adj_npc_raster = replace_missing(Raster(raster_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        # mitn_final_npc_raster = replace_missing(Raster(raster_dir*"final_npc/logmodel_npc/final_npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)

        mitn_access_raster = replace_missing(Raster(raster_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        mitn_adj_access_raster = replace_missing(Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        # mitn_final_access_raster = replace_missing(Raster(raster_dir*"final_access/pmodel_access/final_access_$(year)_$(month_str)_mean.tif"), missingval = NaN)

        mitn_use_raster = replace_missing(Raster(raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        # mitn_final_use_raster = replace_missing(Raster(raster_dir*"final_use/logis_use/final_use_$(year)_$(month_str)_mean.tif"), missingval = NaN)

        # Extract summary for admin0 region based on population weights
        admin0_npc = NaN
        admin0_access = NaN
        admin0_use = NaN

        

        if monthidx ∈ nat_monthidxs
            admin0_npc = country_summaries[(country_summaries.ISO .== ISO) .& 
                                            (country_summaries.month .== month) .&
                                            (country_summaries.year .== year),"NPC_mean"][1]
        end

        hh_survey_slice = hh_survey_data[(hh_survey_data.ISO .== ISO) .& 
                                            (hh_survey_data.interview_month .== month) .&
                                            (hh_survey_data.interview_year .== year),:]

        if monthidx ∈ access_monthidxs
            admin0_access = national_access_aggregated[findfirst(access_monthidxs .== monthidx)]
        else
            hh_survey_slice = hh_survey_data[(hh_survey_data.ISO .== ISO) .& 
                                            (hh_survey_data.interview_month .== month) .&
                                            (hh_survey_data.interview_year .== year),:]
            admin0_access = sum(min.(2 .*hh_survey_slice.n_itn./hh_survey_slice.hh_size,1).*(hh_survey_slice.hh_size).*(hh_survey_slice.hh_sample_wt))/sum((hh_survey_slice.hh_size).*(hh_survey_slice.hh_sample_wt))
        end


        if monthidx ∈ use_monthidxs
            admin0_use = country_summaries[(country_summaries.ISO .== ISO) .& 
                                            (country_summaries.month .== month) .&
                                            (country_summaries.year .== year),"use_mean"][1]
        end

        admin0_id = admin0_shapes_geoIO[admin0_shapes_geoIO.ISO .== ISO,"area_id"][1]
        admin0_name = admin0_shapes_geoIO[admin0_shapes_geoIO.ISO .== ISO,"Name_0"][1]
        admin0_geometry = admin0_shapes_geoIO[admin0_shapes_geoIO.ISO .== ISO,:].geometry

        pop_masked = Rasters.trim(mask(population_raster, with = admin0_geometry), pad = 0)

        admin0_population = sum(pop_masked[findall(.!isnan.(pop_masked))])

        idx = findfirst(model_nat_est_timeseries["filt_ISOs"] .== ISO)
        
        admin0_bv_npc = model_nat_est_timeseries["bv_nat_npc"][idx, monthidx,2]
        admin0_bv_access = model_nat_est_timeseries["bv_nat_access"][idx, monthidx,2]
        admin0_bv_use = model_nat_est_timeseries["bv_nat_use"][idx, monthidx,2]

        admin0_mitn_npc =  model_nat_est_timeseries["mitn_nat_npc"][idx,monthidx,2]
        admin0_mitn_access =  model_nat_est_timeseries["mitn_nat_access"][idx,monthidx,2]
        admin0_mitn_use =  model_nat_est_timeseries["mitn_nat_use"][idx,monthidx,2]

        # Construct country level entry
        output_df_entry = DataFrame(ISO = ISO,
                                        admin_level = 0,
                                        admin_name = admin0_name,
                                        admin_id = admin0_id,
                                        interview_month = month,
                                        interview_year = year,
                                        monthidx = monthidx,
                                        population = admin0_population,
                                        npc = admin0_npc,
                                        access = admin0_access,
                                        use = admin0_use,
                                        bv_npc = admin0_bv_npc,
                                        bv_access = admin0_bv_access,
                                        bv_use = admin0_bv_use,
                                        mitn_npc = admin0_mitn_npc,
                                        mitn_access = admin0_mitn_access,
                                        mitn_use = admin0_mitn_use)

        # Construct Dataframe for year-month survey entries
        df = DataFrame(Type = "country",
                    ISO = data_monthslice.ISO,
                    admin1_name = NaN,
                    area_id = admin0_id,
                    latitude = NaN,
                    longitude = NaN,
                    interview_month = month,
                    interview_year = year,
                    cluster_sample_wt = 1,
                    npc = admin0_npc,
                    bv_npc = admin0_bv_npc,
                    mitn_npc = admin0_mitn_npc,
                    mitn_adj_npc = admin0_mitn_adj_npc,
                    mitn_final_npc = admin0_mitn_final_npc,
                    access = admin0_access,
                    bv_access = admin0_bv_access,
                    mitn_access = admin0_mitn_access,
                    mitn_adj_access = admin0_mitn_adj_access,
                    mitn_final_access = admin0_mitn_final_access,
                    use = admin0_use,
                    bv_use = admin0_bv_use,
                    mitn_use = admin0_mitn_use,
                    mitn_final_use = admin0_mitn_final_use,)

        # Combine all entries and add to collection)
        push!(admin0_df_entries, df)
    end
end
