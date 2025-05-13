"""
Author: Eugene Tan
Date Created: 8/4/2025
Last Updated: 8/4/2025
Extracts time series from rasters and conmbines with Stock and flow estimates to yield a large dataframe summary of the time series.
Saves as a CSV
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
using UsefulTransformations
using RasterLookup

# %% File paths
# MITN Posterior Estimates
snf_post_dir = OUTPUT_DRAWS_DIR*"subnational/"

# Region Admin 1 area id legend
admin1_legend_dir = RAW_DATASET_DIR*"subnational/"
admin1_legend_filename = ADMIN1_AREAID_LEGEND_FILENAME

# Population Rasters directory
pop_dir = POPULATION_RASTER_DIR

# Input and output directory for ITN rasters
raster_dir = OUTPUT_RASTERS_DIR

# Region boundaries
admin0_shapes_geoIO = GeoIO.load(ADMIN0_SHAPEFILE)
admin1_shapes_geoIO = GeoIO.load(ADMIN1_SHAPEFILE)

# INLA Spatial Disaggregation survey data
survey_data_dir = OUTPUT_DATAPREP_DIR
survey_data_filename = INLA_REDUCED_DATAPREP_FILENAME
# Num of samples to import from INLA outputs for use
n_samples = INLA_N_SAMPLES # Number of samples that were saved in the INLA raster process

# %% Get list of countries to analyse
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Time bounds
YEAR_START = 2023#YEAR_NAT_START
YEAR_END = 2023#YEAR_NAT_END

# %% Location to save data to
output_dir = OUTPUT_DIR*"coverage_timeseries/master_extractions_parts/"
output_filename = "extraction_$(YEAR_START)_$(YEAR_END).csv"

# Loop extraction code for each year and month
for year in YEAR_START:YEAR_END
    # Import population raster
    pop_year = min(max(year, 2000), 2020)
    population_raster = replace_missing(Raster(pop_dir*"WorldPop_UNAdj_v3_DRC_fix.$(pop_year).Annual.Data.5km.sum.tif"), missingval = NaN)

    for month in 4:12
        println("Extracting time series for $(year)-$(month)...")
        df0_full_collection = []
        df1_full_collection = []

        # Construct month string for reading rasters
        month_str = "$(month)"
        if month < 10
            month_str = "0$(month)"
        end

        # Calculate month idx
        monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEAR_NAT_START)

        # Import ITN mean rasters
        npc_mean_raster = replace_missing(Raster(raster_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)

        access_mean_raster = replace_missing(Raster(raster_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_mean.tif"), missingval = NaN)

        use_mean_raster = replace_missing(Raster(raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_mean.tif"), missingval = NaN)

        # Import ITN sample rasters
        npc_sample_raster_filenames = readdir(raster_dir*"final_npc/joint_posterior_samples/$(year)_$(month_str)/")
        npc_sample_rasters = replace_missing.(Raster.(raster_dir*"final_npc/joint_posterior_samples/$(year)_$(month_str)/".*npc_sample_raster_filenames), missingval = NaN)

        access_sample_raster_filenames = readdir(raster_dir*"final_access/joint_posterior_samples/$(year)_$(month_str)/")
        access_sample_rasters = replace_missing.(Raster.(raster_dir*"final_access/joint_posterior_samples/$(year)_$(month_str)/".*access_sample_raster_filenames), missingval = NaN)

        raster_dir*"final_use/joint_posterior_samples/$(year)_$(month_str)/"
        use_sample_raster_filenames = readdir(raster_dir*"final_use/joint_posterior_samples/$(year)_$(month_str)/")
        use_sample_rasters = replace_missing.(Raster.(raster_dir*"final_use/joint_posterior_samples/$(year)_$(month_str)/".*use_sample_raster_filenames), missingval = NaN)

        for ISO_i in 1:length(filt_ISOs)
            # Select Country
            ISO = filt_ISOs[ISO_i]
            println("Raster Extraction Country $(ISO_i)/$(length(filt_ISOs)): $(ISO)...")

            # %% Load Subnational Stock and Flow Data
            subnat_snf_filename = "$(ISO)_SUBNAT_draws.jld2"
            subnat_snf_post_draws = JLD2.load(snf_post_dir*subnat_snf_filename)

            # %% Get list of admin 1 regions
            admin1_names = subnat_snf_post_draws["admin1_names"]
            n_admin1 = length(admin1_names)

            # %% Dataframe storage variable for data
            df_admin1_collection = []

            # %% Extract data for each admin1 region
            for admin1_i in 1:n_admin1
                println("Subnational Extraction ($(ISO)): Admin region $(admin1_i) of $(n_admin1).")

                # Metadata
                admin1_name = admin1_names[admin1_i]
                admin1_id = subnat_snf_post_draws["merged_outputs"][admin1_i]["area_id"]
                admin1_geometry = admin1_shapes_geoIO[findfirst(admin1_shapes_geoIO.area_id .== admin1_id),"geometry"]

                ########################
                # SNF ESTIMATES
                ########################
                # Stock and flow NPC values and CI
                snf_npc_lower, snf_npc_upper = quantile(subnat_snf_post_draws["merged_outputs"][admin1_i]["ADJ_NPC_MONTHLY_TOTAL_samples"][:,monthidx], [0.025, 0.975])
                snf_npc_mean = mean(subnat_snf_post_draws["merged_outputs"][admin1_i]["ADJ_NPC_MONTHLY_TOTAL_samples"][:,monthidx])

                # Stock and flow Access values and CI
                snf_access_lower, snf_access_upper = quantile(subnat_snf_post_draws["merged_outputs"][admin1_i]["ADJ_λ_ACCESS_samples"][:,monthidx], [0.025, 0.975])
                snf_access_mean = mean(subnat_snf_post_draws["merged_outputs"][admin1_i]["ADJ_λ_ACCESS_samples"][:,monthidx])

                # Get admin1 population raster
                pop_masked = Rasters.trim(mask(population_raster, with = admin1_geometry); pad=0)
                pop_total = sum(pop_masked[findall(.!isnan.(pop_masked))])

                # ########################
                # # GEOSPATIAL INLA ESTIMATES
                # ########################
                # # Construct masked mean rasters
                # npc_mean_masked = resample(Rasters.trim(mask(npc_mean_raster, with = admin1_geometry); pad = 0), to = pop_masked)
                # access_mean_masked = resample(Rasters.trim(mask(access_mean_raster, with = admin1_geometry); pad = 0), to = pop_masked)
                # use_mean_masked = resample(Rasters.trim(mask(use_mean_raster, with = admin1_geometry); pad = 0), to = pop_masked)

                # # Calculate population weighted means
                # raster_npc_mean = aggregate_raster_weighted_mean(npc_mean_masked, pop_masked)
                # raster_access_mean = aggregate_raster_weighted_mean(access_mean_masked, pop_masked)
                # raster_use_mean = aggregate_raster_weighted_mean(use_mean_masked, pop_masked)

                # # Calculate ITN sample summaries for each sample raster
                # n_npc_samples = length(npc_sample_rasters)
                # npc_raster_vals = zeros(n_npc_samples)
                # println("Calculating NPC sample rasters...")
                # for sample_i in ProgressBar(1:n_npc_samples, leave = false)
                #     npc_sample_raster_masked = resample(Rasters.trim(mask(npc_sample_rasters[sample_i], with = admin1_geometry); pad = 0), to = pop_masked)
                #     npc_raster_vals[sample_i] = aggregate_raster_weighted_mean(npc_sample_raster_masked, pop_masked)
                # end

                # n_access_samples = length(access_sample_rasters)
                # access_raster_vals = zeros(n_access_samples)
                # println("Calculating Access sample rasters...")
                # for sample_i in ProgressBar(1:n_access_samples, leave = false)
                #     access_sample_raster_masked = resample(Rasters.trim(mask(access_sample_rasters[sample_i], with = admin1_geometry); pad = 0), to = pop_masked)
                #     access_raster_vals[sample_i] = aggregate_raster_weighted_mean(access_sample_raster_masked, pop_masked)
                # end

                # n_use_samples = length(use_sample_rasters)
                # use_raster_vals = zeros(n_use_samples)
                # println("Calculating Use sample rasters...")
                # for sample_i in ProgressBar(1:n_use_samples, leave = false)
                #     use_sample_raster_masked = resample(Rasters.trim(mask(use_sample_rasters[sample_i], with = admin1_geometry); pad = 0), to = pop_masked)
                #     use_raster_vals[sample_i] = aggregate_raster_weighted_mean(use_sample_raster_masked, pop_masked)
                # end

                # # Calculate ITN quantiles
                # npc_quantiles = quantile(npc_raster_vals, [0.025, 0.975])
                # access_quantiles = quantile(access_raster_vals, [0.025, 0.975])
                # use_quantiles = quantile(use_raster_vals, [0.025, 0.975])

                # Compile into Dataframe entry
                df1_i = DataFrame(ISO = ISO, category = "Admin1", admin_name = admin1_name, area_id = admin1_id,
                            month = month, year = year,
                            population = pop_total,
                            snf_npc_95lower = snf_npc_lower,
                            snf_npc_mean = snf_npc_mean,
                            snf_npc_95upper = snf_npc_upper,
                            # raster_npc_95lower = npc_quantiles[1],
                            # raster_npc_mean = raster_npc_mean,
                            # raster_npc_95upper = npc_quantiles[2],
                            snf_access_95lower = snf_access_lower,
                            snf_access_mean = snf_access_mean,
                            snf_access_95upper = snf_access_upper,
                            # raster_access_95lower = access_quantiles[1],
                            # raster_access_mean = raster_access_mean,
                            # raster_access_95upper = access_quantiles[2],
                            # raster_use_95lower = use_quantiles[1],
                            # raster_use_mean = raster_use_mean,
                            # raster_use_95upper = use_quantiles[2]
                            )

                push!(df_admin1_collection, df1_i)
            end

            # Concatenate admin1 entries for current country
            df1_i_entries = vcat(df_admin1_collection...)
            push!(df1_full_collection, df1_i_entries)


            #######################################
            # %% ADMIN0: Extract data for admin0 region
            #######################################
            println("National Extraction ($(ISO))")
            admin0_name = admin0_shapes_geoIO[findfirst(admin0_shapes_geoIO.ISO .== ISO),"Name_0"]
            admin0_id = admin0_shapes_geoIO[findfirst(admin0_shapes_geoIO.ISO .== ISO),"area_id"]
            admin0_geometry = admin0_shapes_geoIO[findfirst(admin0_shapes_geoIO.ISO .== ISO),"geometry"]

            # Tally estimates
            pop_total = sum(df1_i_entries.population)

            ########################
            # ADMIN0 SNF ESTIMATES
            ########################
            snf_npc_mean = sum(df1_i_entries.snf_npc_mean .* df1_i_entries.population)/sum(df1_i_entries.population)
            snf_npc_upper = sum(df1_i_entries.snf_npc_95upper .* df1_i_entries.population)/sum(df1_i_entries.population)
            snf_npc_lower = sum(df1_i_entries.snf_npc_95lower .* df1_i_entries.population)/sum(df1_i_entries.population)

            snf_access_mean = sum(df1_i_entries.snf_access_mean .* df1_i_entries.population)/sum(df1_i_entries.population)
            snf_access_upper = sum(df1_i_entries.snf_access_95upper .* df1_i_entries.population)/sum(df1_i_entries.population)
            snf_access_lower = sum(df1_i_entries.snf_access_95lower .* df1_i_entries.population)/sum(df1_i_entries.population)
            
            ########################
            # GEOSPATIAL INLA ESTIMATES
            ########################
            # Get admin1 population raster
            pop_masked = Rasters.trim(mask(population_raster, with = admin0_geometry); pad=0)
            pop_total = sum(pop_masked[findall(.!isnan.(pop_masked))])

            # Construct masked mean rasters
            npc_mean_masked = resample(Rasters.trim(mask(npc_mean_raster, with = admin0_geometry); pad = 0), to = pop_masked)
            access_mean_masked = resample(Rasters.trim(mask(access_mean_raster, with = admin0_geometry); pad = 0), to = pop_masked)
            use_mean_masked = resample(Rasters.trim(mask(use_mean_raster, with = admin0_geometry); pad = 0), to = pop_masked)

            # Calculate population weighted means
            raster_npc_mean = aggregate_raster_weighted_mean(npc_mean_masked, pop_masked)
            raster_access_mean = aggregate_raster_weighted_mean(access_mean_masked, pop_masked)
            raster_use_mean = aggregate_raster_weighted_mean(use_mean_masked, pop_masked)

            # Calculate ITN sample summaries for each sample raster
            n_npc_samples = length(npc_sample_rasters)
            npc_raster_vals = zeros(n_npc_samples)
            println("Calculating NPC sample rasters...")
            for sample_i in ProgressBar(1:n_npc_samples, leave = false)
                npc_sample_raster_masked = resample(Rasters.trim(mask(npc_sample_rasters[sample_i], with = admin0_geometry); pad = 0), to = pop_masked)
                npc_raster_vals[sample_i] = aggregate_raster_weighted_mean(npc_sample_raster_masked, pop_masked)
            end

            n_access_samples = length(access_sample_rasters)
            access_raster_vals = zeros(n_access_samples)
            println("Calculating Access sample rasters...")
            for sample_i in ProgressBar(1:n_access_samples, leave = false)
                access_sample_raster_masked = resample(Rasters.trim(mask(access_sample_rasters[sample_i], with = admin0_geometry); pad = 0), to = pop_masked)
                access_raster_vals[sample_i] = aggregate_raster_weighted_mean(access_sample_raster_masked, pop_masked)
            end

            n_use_samples = length(use_sample_rasters)
            use_raster_vals = zeros(n_use_samples)
            println("Calculating Use sample rasters...")
            for sample_i in ProgressBar(1:n_use_samples, leave = false)
                use_sample_raster_masked = resample(Rasters.trim(mask(use_sample_rasters[sample_i], with = admin0_geometry); pad = 0), to = pop_masked)
                use_raster_vals[sample_i] = aggregate_raster_weighted_mean(use_sample_raster_masked, pop_masked)
            end
            
            raster_npc_lower, raster_npc_upper = quantile(npc_raster_vals, [0.025, 0.975])
            raster_access_lower, raster_access_upper = quantile(access_raster_vals, [0.025, 0.975])
            raster_use_lower, raster_use_upper = quantile(use_raster_vals, [0.025, 0.975])

            # Construct Dataframe entry
            df0_i = DataFrame(ISO = ISO, category = "Admin0", admin_name = admin0_name, area_id = admin0_id,
                            month = month, year = year,
                            population = pop_total,
                            snf_npc_95lower = snf_npc_lower,
                            snf_npc_mean = snf_npc_mean,
                            snf_npc_95upper = snf_npc_upper,
                            raster_npc_95lower = raster_npc_lower,
                            raster_npc_mean = raster_npc_mean,
                            raster_npc_95upper = raster_npc_upper,
                            snf_access_95lower = snf_access_lower,
                            snf_access_mean = snf_access_mean,
                            snf_access_95upper = snf_access_upper,
                            raster_access_95lower = raster_access_lower,
                            raster_access_mean = raster_access_mean,
                            raster_access_95upper = raster_access_upper,
                            raster_use_95lower = raster_use_lower,
                            raster_use_mean = raster_use_mean,
                            raster_use_95upper = raster_use_upper)

            # Add data frame to admin0 df collection
            push!(df0_full_collection, df0_i)
        end

        # %% Combine collections together and save
        # compiled_dataframe = vcat(df0_full_collection...,df1_full_collection...)
        compiled_dataframe = vcat(df0_full_collection...)
        mkpath(output_dir)
        CSV.write(output_dir*"extraction_$(year)_$(month).csv", compiled_dataframe)
    end
end

