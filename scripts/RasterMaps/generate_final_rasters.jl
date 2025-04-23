"""
Author: Eugene Tan
Date Created: 16/1/2025
Last Updated: 16/3/2025
Combines all raster components to get final maps of ITN coverage and saves rasters. Makes the following components
1. Creates raster of subnational stock and flow outputs
2. Combines with INLA outputs of deviation metrics to produce raw disaggregated values of npc and access
3. Rakes raw npc and access disaggregated rasters against stock and flow NATIONAL estimates (attempts to preserve 95% CI)
4. Calculates a global use raster using raked (adjusted) rasters
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

# %% File paths
# MITN Posterior Estimates
snf_post_dir = OUTPUT_DRAWS_DIR*"subnational/"

# Region Admin 1 area id legend
admin1_legend_dir = RAW_DATASET_DIR*"subnational/"
admin1_legend_filename = ADMIN1_AREAID_LEGEND_FILENAME

# Population Rasters directory
pop_dir = POPULATION_RASTER_DIR

# Region boundaries
admin0_shapes_geoIO = GeoIO.load(ADMIN0_SHAPEFILE)
admin1_shapes_geoIO = GeoIO.load(ADMIN1_SHAPEFILE)

# INLA Spatial Disaggregation survey data
survey_data_dir = OUTPUT_DATAPREP_DIR
survey_data_filename = INLA_REDUCED_DATAPREP_FILENAME

# Num of samples to import from INLA outputs for rasters
n_samples = INLA_UNCERTAINTY_N_SAMPLES

# Get base raster with required resolution to build from
raster_base = replace_missing(Raster("outputs/rasters/inla_logmodel_npc/NPC_logmodel_$(2000)_mean.tif"), missingval = -NaN)

# Input and output directory for rasters
inla_dir = OUTPUT_RASTERS_DIR
raster_dir = OUTPUT_RASTERS_DIR
input_dir = OUTPUT_RASTERS_DIR
output_dir = OUTPUT_RASTERS_DIR

# Make paths tosave location if doesn't already exist
mkpath(output_dir*"final_npc/snf_npc/")
mkpath(output_dir*"final_npc/logmodel_npc/")
mkpath(output_dir*"final_access/snf_access/")
mkpath(output_dir*"final_access/pmodel_access/")
mkpath(output_dir*"final_use/logis_use/")

# %% Perform draws and save outputs. Filter out unwanted countries
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Time bounds
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

# %%

# %% Loop to first construct SNF block map upto subnational resolution
for year in ProgressBar(YEAR_START:YEAR_END)
    for month in 1:12
        monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEAR_START)

        println("Constructing stock and flow rasters year [$(year)/$(YEAR_END)], month [$(month)/12]")

        # %% Declare storage variables
        npc_nat_snf_mean_rasters = Vector{Any}(undef, length(filt_ISOs))
        npc_nat_snf_upper_rasters = Vector{Any}(undef, length(filt_ISOs))
        npc_nat_snf_lower_rasters = Vector{Any}(undef, length(filt_ISOs))

        access_nat_snf_mean_rasters = Vector{Any}(undef, length(filt_ISOs))
        access_nat_snf_upper_rasters = Vector{Any}(undef, length(filt_ISOs))
        access_nat_snf_lower_rasters = Vector{Any}(undef, length(filt_ISOs))

        println("Constructing rasters using subnational SNF values...")
        for ISO_i in ProgressBar(1:length(filt_ISOs), leave = false)
            ISO = filt_ISOs[ISO_i]
            snf_post_filename = "$(ISO)_SUBNAT_draws.jld2"
            snf_post_draws = load(snf_post_dir*snf_post_filename)
            n_admin1 = length(snf_post_draws["merged_outputs"])

            # Declare storage variables
            npc_subnat_snf_mean_rasters = []
            npc_subnat_snf_upper_rasters = []
            npc_subnat_snf_lower_rasters = []

            access_subnat_snf_mean_rasters = []
            access_subnat_snf_upper_rasters = []
            access_subnat_snf_lower_rasters = []

            for subnat_i in ProgressBar(1:n_admin1, leave = false)
                admin1_id = snf_post_draws["merged_outputs"][subnat_i]["area_id"]
                npc_subnat_mean_estimate = mean(snf_post_draws["merged_outputs"][subnat_i]["ADJ_NPC_MONTHLY_TOTAL_samples"][:,monthidx])
                npc_subnat_upper_estimate = quantile(snf_post_draws["merged_outputs"][subnat_i]["ADJ_NPC_MONTHLY_TOTAL_samples"][:,monthidx], 0.975)
                npc_subnat_lower_estimate = quantile(snf_post_draws["merged_outputs"][subnat_i]["ADJ_NPC_MONTHLY_TOTAL_samples"][:,monthidx], 0.025)

                # Remove NaNs for this draw possibly due to transition period in problematic countries (MIGHT WANT TO CHECK)
                access_subnat_draws = snf_post_draws["merged_outputs"][subnat_i]["ADJ_λ_ACCESS_samples"][:,monthidx]
                access_subnat_draws = access_subnat_draws[findall(.!isnan.(access_subnat_draws))]
                access_subnat_mean_estimate = mean(access_subnat_draws)
                access_subnat_upper_estimate = quantile(access_subnat_draws, 0.975)
                access_subnat_lower_estimate = quantile(access_subnat_draws, 0.025)
                
                if isnan(access_subnat_mean_estimate) # i.e. no data available, so no valid access. Assume 0 access
                    access_subnat_mean_estimate = 0
                    access_subnat_std_estimate = 0
                end
                
                # Get required geometry
                admin1_geometry = admin1_shapes_geoIO[admin1_shapes_geoIO.area_id .== admin1_id,:].geometry

                # Mask raster base and trim to desired subregion and set values to SNF estimates of NPC and access
                subnat_masked = Rasters.trim(mask(raster_base, with = admin1_geometry); pad=0)
                
                 # check if able to find region in default tiling. If region is too small, then skip in analysis
                if size(raster_base) == size(subnat_masked)
                    continue
                end

                npc_subnat_mean_raster = copy(subnat_masked)
                npc_subnat_upper_raster = copy(subnat_masked)
                npc_subnat_lower_raster = copy(subnat_masked)
                access_subnat_mean_raster = copy(subnat_masked)
                access_subnat_upper_raster = copy(subnat_masked)
                access_subnat_lower_raster = copy(subnat_masked)

                # Check if regions were valid to begin with
                if isempty(findall(.!isnan.(npc_subnat_mean_raster))) # For some reason if not accounted for, causes problems with mosaic()
                    continue
                end

                npc_subnat_mean_raster[findall(.!isnan.(npc_subnat_mean_raster))] .= npc_subnat_mean_estimate
                npc_subnat_upper_raster[findall(.!isnan.(npc_subnat_upper_raster))] .= npc_subnat_upper_estimate
                npc_subnat_lower_raster[findall(.!isnan.(npc_subnat_lower_raster))] .= npc_subnat_lower_estimate
                access_subnat_mean_raster[findall(.!isnan.(access_subnat_mean_raster))] .= access_subnat_mean_estimate
                access_subnat_upper_raster[findall(.!isnan.(access_subnat_upper_raster))] .= access_subnat_upper_estimate
                access_subnat_lower_raster[findall(.!isnan.(access_subnat_lower_raster))] .= access_subnat_lower_estimate

                push!(npc_subnat_snf_mean_rasters, npc_subnat_mean_raster)
                push!(npc_subnat_snf_upper_rasters, npc_subnat_upper_raster)
                push!(npc_subnat_snf_lower_rasters, npc_subnat_lower_raster)

                push!(access_subnat_snf_mean_rasters, access_subnat_mean_raster)
                push!(access_subnat_snf_upper_rasters, access_subnat_upper_raster)
                push!(access_subnat_snf_lower_rasters, access_subnat_lower_raster)
            end

            npc_nat_snf_mean_rasters[ISO_i] = copy(npc_subnat_snf_mean_rasters)
            npc_nat_snf_upper_rasters[ISO_i] = copy(npc_subnat_snf_upper_rasters)
            npc_nat_snf_lower_rasters[ISO_i] = copy(npc_subnat_snf_lower_rasters)

            access_nat_snf_mean_rasters[ISO_i] = copy(access_subnat_snf_mean_rasters)
            access_nat_snf_upper_rasters[ISO_i] = copy(access_subnat_snf_upper_rasters)
            access_nat_snf_lower_rasters[ISO_i] = copy(access_subnat_snf_lower_rasters)
        end

        # Combine Subnational rasters into national level rasters
        Threads.@threads for ISO_i in ProgressBar(1:length(filt_ISOs), leave = false)
            npc_nat_snf_mean_rasters[ISO_i] = mosaic(first, npc_nat_snf_mean_rasters[ISO_i]..., atol = 0.01)
            npc_nat_snf_upper_rasters[ISO_i] =  mosaic(first, npc_nat_snf_upper_rasters[ISO_i]..., atol = 0.01)
            npc_nat_snf_lower_rasters[ISO_i] =  mosaic(first, npc_nat_snf_lower_rasters[ISO_i]..., atol = 0.01)

            access_nat_snf_mean_rasters[ISO_i] = mosaic(first, access_nat_snf_mean_rasters[ISO_i]..., atol = 0.01)
            access_nat_snf_upper_rasters[ISO_i] =  mosaic(first, access_nat_snf_upper_rasters[ISO_i]..., atol = 0.01)
            access_nat_snf_lower_rasters[ISO_i] =  mosaic(first, access_nat_snf_lower_rasters[ISO_i]..., atol = 0.01)
        end
        
        # Combine national level rasters into Africa raster
        # Month string
        month_str = "$(month)"
        if month < 10
            month_str = "0$(month)"
        end

        # Mosaic and save to disk.
        println("Combining all national rasters into Africa level raster and write to disk...")
        for thread_i in 1:6
            if thread_i == 1
                println("Constructing NPC SNF mean raster...")
                combined_npc_snf_mean_raster = resample(mosaic(first, npc_nat_snf_mean_rasters..., atol = 0.01), to = raster_base)
                write(output_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_mean.tif", combined_npc_snf_mean_raster, force = true)
                println("Constructed NPC SNF mean raster...")
            elseif thread_i == 2
                println("Constructing NPC SNF upper raster...")
                combined_npc_snf_upper_raster = resample(mosaic(first, npc_nat_snf_upper_rasters..., atol = 0.01), to = raster_base)
                write(output_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_upper.tif", combined_npc_snf_upper_raster, force = true)
                println("Constructed NPC SNF upper raster...")
            elseif thread_i == 3
                println("Constructing NPC SNF lower raster...")
                combined_npc_snf_lower_raster = resample(mosaic(first, npc_nat_snf_lower_rasters..., atol = 0.01), to = raster_base)
                write(output_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_lower.tif", combined_npc_snf_lower_raster, force = true)
                println("Constructed NPC SNF lower raster...")
            elseif thread_i == 4
                println("Constructing Access SNF mean raster...")
                combined_access_snf_mean_raster = resample(mosaic(first, access_nat_snf_mean_rasters..., atol = 0.01), to = raster_base)
                write(output_dir*"final_access/snf_access/access_$(year)_$(month_str)_mean.tif", combined_access_snf_mean_raster, force = true)
                println("Constructed Access SNF mean raster...")
            elseif thread_i == 5
                println("Constructing Access SNF upper raster...")
                combined_access_snf_upper_raster = resample(mosaic(first, access_nat_snf_upper_rasters..., atol = 0.01), to = raster_base)
                write(output_dir*"final_access/snf_access/access_$(year)_$(month_str)_upper.tif", combined_access_snf_upper_raster, force = true)
                println("Constructed Access SNF upper raster...")
            else
                println("Constructing Access SNF lower raster...")
                combined_access_snf_lower_raster = resample(mosaic(first, access_nat_snf_lower_rasters..., atol = 0.01), to = raster_base)
                write(output_dir*"final_access/snf_access/access_$(year)_$(month_str)_lower.tif", combined_access_snf_lower_raster, force = true)
                println("Constructed Access SNF lower raster...")
            end
        end

        println("Raster construction complete.")
    end
end

# %% Loop to construct unadjusted NPC, Access and Use rasters using INLA deviation regression outputs
for year in YEAR_START:YEAR_END
    for month in 1:12
        println("Constructing spatial disaggregated rasters year [$(year)/$(YEAR_END)], month [$(month)/12]")

        # Get month string for importing files
        month_str = "$(month)"
        if month < 10
            month_str = "0$(month)"
        end
        
        println("Importing rasters...")
        # Import stock and flow rasters
        npc_snf_mean_raster = replace_missing(Raster(input_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        npc_snf_upper_raster = replace_missing(Raster(input_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_upper.tif"), missingval = NaN)
        npc_snf_lower_raster = replace_missing(Raster(input_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_lower.tif"), missingval = NaN)
        access_snf_mean_raster = replace_missing(Raster(input_dir*"final_access/snf_access/access_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        access_snf_upper_raster = replace_missing(Raster(input_dir*"final_access/snf_access/access_$(year)_$(month_str)_upper.tif"), missingval = NaN)
        access_snf_lower_raster = replace_missing(Raster(input_dir*"final_access/snf_access/access_$(year)_$(month_str)_lower.tif"), missingval = NaN)

        #################################
        # Construct mean rasters
        #################################
        # Import log model npc rasters
        logmodel_npc_mean_raster = replace_missing(Raster(inla_dir*"inla_logmodel_npc/NPC_logmodel_$(year)_mean.tif"), missingval = NaN)

        # Import pmodel access rasters
        pmodel_access_mean_raster =  replace_missing(Raster(inla_dir*"inla_pmodel_access/ACCESS_pmodel_$(year)_mean.tif"), missingval = NaN)

        # Import logismodel use rasters
        logis_use_mean_raster = replace_missing(Raster(inla_dir*"inla_use_logis/USE_logismodel_$(year)_$(month)_mean.tif"), missingval = NaN)

        # TEMP - Used just to make sure rasters are all aligned
        raster_base = copy(logmodel_npc_mean_raster)

        # Consruct rasters
        println("Calculating mean NPC rasters...")
        logmodel_npc_ratio_mean_raster = exp.(logmodel_npc_mean_raster)
        npc_map_mean_raster = logmodel_npc_ratio_mean_raster.* npc_snf_mean_raster
        println("Calculating mean Access rasters...")
        access_map_mean_raster = inv_p_transform.(pmodel_access_mean_raster, access_snf_mean_raster, n=2)
        println("Calculating mean Use rasters...")
        use_mean_raster = inv_p_transform.(logis_use_mean_raster, access_map_mean_raster, n=2)

        # Save mean rasters
        println("Saving mean rasters...")
        write(output_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_mean.tif", npc_map_mean_raster, force = true)
        write(output_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_mean.tif", access_map_mean_raster, force = true)
        write(output_dir*"final_use/logis_use/use_$(year)_$(month_str)_mean.tif", use_mean_raster, force = true)

        #################################
        # Construct posterior sample rasters
        #################################
        mkpath(output_dir*"final_npc/joint_posterior_samples/$(year)_$(month_str)/")
        mkpath(output_dir*"final_access/joint_posterior_samples/$(year)_$(month_str)/")
        mkpath(output_dir*"final_use/joint_posterior_samples/$(year)_$(month_str)/")
        

        # Calculate spatial disaggregated raster for each INLA posterior draw
        println("Constructing NPC and Access joint posterior sample rasters...")
        for sample_i in ProgressBar(1:n_samples, leave = false)
            # NPC
            logmodel_npc_sample_raster = replace_missing(Raster(inla_dir*"inla_logmodel_npc/NPC_logmodel_$(year)_sample_$(sample_i).tif"), missingval = NaN)
            logmodel_npc_ratio_sample_raster = exp.(logmodel_npc_sample_raster)
            npc_map_mean_sample_raster = logmodel_npc_ratio_sample_raster .* npc_snf_mean_raster
            npc_map_upper_sample_raster = logmodel_npc_ratio_sample_raster .* npc_snf_upper_raster
            npc_map_lower_sample_raster = logmodel_npc_ratio_sample_raster .* npc_snf_lower_raster

            write(output_dir*"final_npc/joint_posterior_samples/$(year)_$(month_str)/npc_mean_snf_sample_$(sample_i).tif", npc_map_mean_sample_raster, force = true)
            write(output_dir*"final_npc/joint_posterior_samples/$(year)_$(month_str)/npc_upper_snf_sample_$(sample_i).tif", npc_map_upper_sample_raster, force = true)
            write(output_dir*"final_npc/joint_posterior_samples/$(year)_$(month_str)/npc_lower_snf_sample_$(sample_i).tif", npc_map_lower_sample_raster, force = true)

            # Access
            pmodel_access_sample_raster = replace_missing(Raster(inla_dir*"inla_pmodel_access/ACCESS_pmodel_$(year)_sample_$(sample_i).tif"), missingval = NaN)
            access_map_mean_sample_raster = inv_p_transform.(pmodel_access_sample_raster, access_snf_mean_raster, n=2)
            access_map_upper_sample_raster = inv_p_transform.(pmodel_access_sample_raster, access_snf_upper_raster, n=2)
            access_map_lower_sample_raster = inv_p_transform.(pmodel_access_sample_raster, access_snf_lower_raster, n=2)
            
            write(output_dir*"final_access/joint_posterior_samples/$(year)_$(month_str)/access_mean_snf_sample_$(sample_i).tif", access_map_mean_sample_raster, force = true)
            write(output_dir*"final_access/joint_posterior_samples/$(year)_$(month_str)/access_upper_snf_sample_$(sample_i).tif", access_map_upper_sample_raster, force = true)
            write(output_dir*"final_access/joint_posterior_samples/$(year)_$(month_str)/access_lower_snf_sample_$(sample_i).tif", access_map_lower_sample_raster, force = true)
        end

        println("Constructing Use joint posterior sample rasters...")
        for sample_i in ProgressBar(1:n_samples, leave = false)
            access_sample_i = rand(1:n_samples)
            use_sample_i = rand(1:n_samples)

            # Import Use model rasters
            logis_use_sample_raster = replace_missing(Raster(inla_dir*"inla_use_logis/USE_logismodel_$(year)_$(month)_sample_$(use_sample_i).tif"), missingval = NaN)
            
            # Import constructed access sample rasters
            access_mean_sample_raster = replace_missing(Raster(inla_dir*"final_access/joint_posterior_samples/$(year)_$(month_str)/access_mean_snf_sample_$(access_sample_i).tif"), missingval = NaN)
            access_upper_sample_raster = replace_missing(Raster(inla_dir*"final_access/joint_posterior_samples/$(year)_$(month_str)/access_upper_snf_sample_$(access_sample_i).tif"), missingval = NaN)
            access_lower_sample_raster = replace_missing(Raster(inla_dir*"final_access/joint_posterior_samples/$(year)_$(month_str)/access_lower_snf_sample_$(access_sample_i).tif"), missingval = NaN)
            
            # Construct sample use rasters
            use_map_mean_sample_raster = inv_p_transform.(logis_use_sample_raster, access_mean_sample_raster, n=2)
            use_map_upper_sample_raster = inv_p_transform.(logis_use_sample_raster, access_upper_sample_raster, n=2)
            use_map_lower_sample_raster = inv_p_transform.(logis_use_sample_raster, access_lower_sample_raster, n=2)

            # Save rasters
            write(output_dir*"final_use/joint_posterior_samples/$(year)_$(month_str)/use_mean_snf_sample_$(sample_i).tif", use_map_mean_sample_raster, force = true)
            write(output_dir*"final_use/joint_posterior_samples/$(year)_$(month_str)/use_upper_snf_sample_$(sample_i).tif", use_map_upper_sample_raster, force = true)
            write(output_dir*"final_use/joint_posterior_samples/$(year)_$(month_str)/use_lower_snf_sample_$(sample_i).tif", use_map_lower_sample_raster, force = true)
        end

        println("Raster construction complete.")
    end
end

# # %% # Calculate adjusted NPC and Access Rasters (Raking on a country level)
# # Subnat estimates of NPC should be close to national estimates when summed by population (SNF Calibration step)
# # Subnat estimates of access may be considered unreliable for raking as it assumes household demographic distribution to be consistent across whole country

# # Load INLA dataset
# survey_data = CSV.read(OUTPUT_DATAPREP_DIR*HOUSEHOLD_SURVEY_DATA_FILENAME, DataFrame)
# survey_data = survey_data[.!ismissing.(survey_data.latitude),:]

# for year in ProgressBar(YEAR_START:YEAR_END, leave = false)
#     # Import population raster
#     pop_year = min(max(year, 2000), 2020)
    
#     population_raster = replace_missing(Raster(pop_dir*"WorldPop_UNAdj_v3_DRC_fix.$(pop_year).Annual.Data.5km.sum.tif"), missingval = NaN)

#     for month in ProgressBar(1:12, leave = false)
#         println("Processing raster year [$(year)/$(YEAR_END)], month [$(month)/12]")

#         # Calculate reference monthidx to access data
#         monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEAR_NAT_START)

#         # Get correct year/month string for importing file
#         year_str = "$(year)"
#         month_str = "$(month)"
#         if month < 10
#             month_str = "0"*month_str
#         end

#         println("Importing rasters...")
#         # Import calculated NPC rasters
#         npc_mean_raster = replace_missing(Raster(input_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)
#         npc_upper_raster = replace_missing(Raster(input_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_upper.tif"), missingval = NaN)
#         npc_lower_raster = replace_missing(Raster(input_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_lower.tif"), missingval = NaN)
        
#         # Import calculated Access rasters 
#         access_mean_raster = replace_missing(Raster(input_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_mean.tif"), missingval = NaN)
#         access_upper_raster = replace_missing(Raster(input_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_upper.tif"), missingval = NaN)
#         access_lower_raster = replace_missing(Raster(input_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_lower.tif"), missingval = NaN)

#         # Storage variables for country rasters of ITN coverage with CI
#         adj_npc_nat_mean_rasters = Vector{Raster}(undef, length(filt_ISOs))
#         adj_npc_nat_upper_rasters = Vector{Raster}(undef, length(filt_ISOs))
#         adj_npc_nat_lower_rasters = Vector{Raster}(undef, length(filt_ISOs))

#         adj_access_nat_mean_rasters = Vector{Raster}(undef, length(filt_ISOs))
#         adj_access_nat_upper_rasters = Vector{Raster}(undef, length(filt_ISOs))
#         adj_access_nat_lower_rasters = Vector{Raster}(undef, length(filt_ISOs))

#         println("Raking rasters using SNF National estimates...")
#         # Do extractions and calculations for each country
#         for ISO_i in ProgressBar(1:length(filt_ISOs), leave = false)
#             ISO = filt_ISOs[ISO_i]

#             # Import subnational SNF draws and construct raked rasters
#             subnat_snf_post_draws = load(snf_post_dir*"$(ISO)_SUBNAT_draws.jld2")
#             n_admin1 = length(subnat_snf_post_draws["admin1_names"])

#             adj_npc_mean_rasters = Vector{Any}(undef, n_admin1)
#             adj_npc_upper_rasters = Vector{Any}(undef, n_admin1)
#             adj_npc_lower_rasters = Vector{Any}(undef, n_admin1)

#             adj_access_mean_rasters = Vector{Any}(undef, n_admin1)
#             adj_access_upper_rasters = Vector{Any}(undef, n_admin1)
#             adj_access_lower_rasters = Vector{Any}(undef, n_admin1)

#             for admin1_i in ProgressBar(1:n_admin1, leave = false)

#                 subnat_snf_npc_mean = mean(subnat_snf_post_draws["merged_outputs"][admin1_i]["ADJ_NPC_MONTHLY_TOTAL_samples"][:,monthidx])
#                 subnat_snf_npc_std = std(subnat_snf_post_draws["merged_outputs"][admin1_i]["ADJ_NPC_MONTHLY_TOTAL_samples"][:,monthidx])
                
#                 subnat_snf_access_mean = mean(subnat_snf_post_draws["merged_outputs"][admin1_i]["ADJ_λ_ACCESS_samples"][:,monthidx])
#                 subnat_snf_access_std = std(subnat_snf_post_draws["merged_outputs"][admin1_i]["ADJ_λ_ACCESS_samples"][:,monthidx])

#                 # Get region geometry to mask raster
#                 area_id = subnat_snf_post_draws["merged_outputs"][admin1_i]["area_id"]
#                 admin1_geometry = admin1_shapes_geoIO[findfirst(admin1_shapes_geoIO.area_id .== area_id),"geometry"]

#                 # Get masked + trimmed versions of each required component raster
#                 ## Population
#                 pop_masked = Rasters.trim(mask(population_raster, with = admin1_geometry); pad=0)
                
#                 ## NPC
#                 npc_mean_masked = resample(Rasters.trim(mask(npc_mean_raster, with = admin1_geometry); pad=0), to = pop_masked)
#                 npc_upper_masked = resample(Rasters.trim(mask(npc_upper_raster, with = admin1_geometry); pad=0), to = pop_masked)
#                 npc_lower_masked = resample(Rasters.trim(mask(npc_lower_raster, with = admin1_geometry); pad=0), to = pop_masked)

#                 ## Access
#                 access_mean_masked = resample(Rasters.trim(mask(access_mean_raster, with = admin1_geometry); pad=0), to = pop_masked)
#                 access_upper_masked = resample(Rasters.trim(mask(access_upper_raster, with = admin1_geometry); pad=0), to = pop_masked)
#                 access_lower_masked = resample(Rasters.trim(mask(access_lower_raster, with = admin1_geometry); pad=0), to = pop_masked)
                
#                 #####
#                 # Check if there is survey data to rake against. If yes get survey data slice to rake data against
#                 #####

#                 # Base estimate of scaling k (i.e. take fully from SNF estimate)
#                 # for NPC
#                 nonmissing_idx_npc_mean = intersect(findall(.!isnan.(pop_masked)), findall(.!isnan.(npc_mean_masked)))
#                 npc_spatial_mean_estimate = sum(pop_masked[nonmissing_idx_npc_mean].*npc_mean_masked[nonmissing_idx_npc_mean])/sum(pop_masked[nonmissing_idx_npc_mean])

#                 npc_scaling_k = subnat_snf_npc_mean/npc_spatial_mean_estimate
#                 if (npc_spatial_mean_estimate == 0) && isnan(npc_scaling_k)
#                     npc_scaling_k = 1
#                 end

#                 # for Access
#                 nonmissing_idx_access_mean = intersect(findall(.!isnan.(pop_masked)), findall(.!isnan.(access_mean_masked)))
#                 access_spatial_mean_estimate = sum(pop_masked[nonmissing_idx_npc_mean].*access_mean_masked[nonmissing_idx_npc_mean])/sum(pop_masked[nonmissing_idx_npc_mean])

#                 access_scaling_k = subnat_snf_access_mean/access_spatial_mean_estimate

#                 if (access_spatial_mean_estimate == 0) && isnan(access_scaling_k)
#                     access_scaling_k = 1
#                 end

#                 # Check survey data if there is any raw surveys present (need at least 20. Arbitrarily chosen)
#                 survey_data_slice = survey_data[survey_data.ISO .== ISO .&&
#                                                 survey_data.interview_year .== year .&&
#                                                 survey_data.interview_month .== month .&&
#                                                 .!ismissing.(survey_data.area_id),:]
#                 survey_data_slice = survey_data_slice[survey_data_slice.area_id .== area_id,:]

#                 n_survey_datapoints = size(survey_data_slice)[1]

#                 if n_survey_datapoints > 20
#                     # Get number of unique lat longs
#                     sid_lat_lons = unique([(survey_data_slice.latitude[i], 
#                                             survey_data_slice.longitude[i]) for i in 1:size(survey_data_slice)[1]])
#                     n_unique_lat_lons = length(sid_lat_lons)

#                     # Create storage variable with sample size, standard error, lat long estimate from INLA model, survey point estimate
#                     npc_rake_data = zeros(n_unique_lat_lons,4)
#                     access_rake_data = zeros(n_unique_lat_lons,4)

#                     for i in 1:n_unique_lat_lons
#                         # Target lat, lon
#                         lat, lon = sid_lat_lons[i]
                        
#                         # Get list of all points that relate to target lat-lon
#                         survey_points = survey_data_slice[survey_data_slice.latitude .== lat .&&
#                                                         survey_data_slice.longitude .== lon,:]
#                         names(survey_points)
#                         # Calculate size of sample
#                         n_sample = size(survey_points)[1]

#                         # Calculate local observation estimate
#                         n_itn = survey_points.n_itn
#                         hh_size = survey_points.hh_size
#                         hh_wt = survey_points.hh_sample_wt
#                         obs_npc = sum(n_itn.*hh_wt)/sum(hh_size.*hh_wt)
#                         obs_access = sum(min.(2 .* n_itn./hh_size, 1) .* hh_size .* hh_wt)/sum(hh_size .* hh_wt)
#                         std_npc = 1
#                         std_access = 1 # Baseline uninformative values if sample size onle = 1
#                         if n_sample > 10 # If there was more than 1 sample, std must be defined
#                             std_npc = sum((((n_itn./hh_size) .- obs_npc).^2) .* hh_wt)/sum(hh_wt)
#                             std_access = sum((min.(2 .* n_itn./hh_size, 1) .- obs_access).^2 .* hh_size .* hh_wt)/sum(hh_size .* hh_wt)
#                         end
#                         std_npc
#                         std_err_npc = std_npc/sqrt(n_sample)
#                         std_err_access = std_access/sqrt(n_sample)

#                         # Find index of raster corresponding to latlon
#                         npc_model_lats = lookup(npc_mean_masked, Y)
#                         npc_model_lons = lookup(npc_mean_masked, X)
#                         npc_model_lat_idx = argmin(abs.(npc_model_lats .- lat))
#                         npc_model_lon_idx = argmin(abs.(npc_model_lons .- lon))

#                         access_model_lats = lookup(access_mean_masked, Y)
#                         access_model_lons = lookup(access_mean_masked, X)
#                         access_model_lat_idx = argmin(abs.(access_model_lats .- lat))
#                         access_model_lon_idx = argmin(abs.(access_model_lons .- lon))

#                         # Extract required value
#                         npc_rake_data[i,:] .= n_sample, std_err_npc, obs_npc, Float64(npc_mean_masked[npc_model_lon_idx, npc_model_lat_idx])
#                         access_rake_data[i,:] .= n_sample, std_err_access, obs_access, Float64(access_mean_masked[access_model_lon_idx, access_model_lat_idx])
#                     end

#                     # Get rid of all NaN entries
#                     npc_nanidx = union(findall(isnan.(npc_rake_data[:,3])),findall(isnan.(npc_rake_data[:,4])))
#                     access_nanidx = union(findall(isnan.(access_rake_data[:,3])),findall(isnan.(access_rake_data[:,4])))

#                     npc_rake_data = npc_rake_data[setdiff(1:size(npc_rake_data)[1],npc_nanidx),:]
#                     access_rake_data = access_rake_data[setdiff(1:size(access_rake_data)[1],access_nanidx),:]

#                     # Calculate modified constant that minimises the error between both raw survey data and SNF estimates. Assuming errors are gaussian.
#                     npc_RMSE = 0.1 # Taken from from model fit values
#                     npc_A = npc_spatial_mean_estimate*subnat_snf_npc_mean/((subnat_snf_npc_std^2))
#                     npc_B = sum(npc_rake_data[:,3].*npc_rake_data[:,4]./(npc_rake_data[:,2].^2 .+ npc_RMSE^2))
#                     npc_C = (npc_spatial_mean_estimate^2)/(subnat_snf_npc_std^2)
#                     npc_D = sum((npc_rake_data[:,4].^2)./(npc_rake_data[:,2].^2 .+ npc_RMSE^2))

#                     npc_scaling_k = (npc_A+npc_B)/(npc_C+npc_D)

#                     access_RMSE = 0.1 # Taken from from model fit values
#                     access_A = access_spatial_mean_estimate*subnat_snf_access_mean/((subnat_snf_access_std^2))
#                     access_B = sum(access_rake_data[:,3].*access_rake_data[:,4]./(access_rake_data[:,2].^2 .+ access_RMSE^2))
#                     access_C = (access_spatial_mean_estimate^2)/(subnat_snf_access_std^2)
#                     access_D = sum((access_rake_data[:,4].^2)./(access_rake_data[:,2].^2 .+ access_RMSE^2))

#                     access_scaling_k = (access_A+access_B)/(access_C+access_D)
#                 end

#                 # Calculate adjusted maps/rasters using the scaling constant
#                 adj_npc_mean_masked = npc_scaling_k.*npc_mean_masked
#                 adj_npc_upper_masked = npc_scaling_k.*npc_upper_masked
#                 adj_npc_lower_masked = npc_scaling_k.*npc_lower_masked

#                 adj_access_mean_masked = access_scaling_k.*access_mean_masked
#                 adj_access_upper_masked = access_scaling_k.*access_upper_masked
#                 adj_access_lower_masked = access_scaling_k.*access_lower_masked
                
#                 # Store raked rasters into collection
#                 adj_npc_mean_rasters[admin1_i] = adj_npc_mean_masked
#                 adj_npc_upper_rasters[admin1_i] = adj_npc_upper_masked
#                 adj_npc_lower_rasters[admin1_i] = adj_npc_lower_masked
#                 plot(npc_mean_masked)
#                 adj_access_mean_rasters[admin1_i] = adj_access_mean_masked
#                 adj_access_upper_rasters[admin1_i] = adj_access_upper_masked
#                 adj_access_lower_rasters[admin1_i] = adj_access_lower_masked
#             end

#             # Combine subnational adjusted rasters into national and then save in collection variable
            
#             adj_npc_nat_mean_rasters[ISO_i] = mosaic(first, adj_npc_mean_rasters..., atol = 0.01)
#             adj_npc_nat_upper_rasters[ISO_i] = mosaic(first, adj_npc_upper_rasters..., atol = 0.01)
#             adj_npc_nat_lower_rasters[ISO_i] = mosaic(first, adj_npc_lower_rasters..., atol = 0.01)

#             adj_access_nat_mean_rasters[ISO_i] = mosaic(first, adj_access_mean_rasters..., atol = 0.01)
#             adj_access_nat_upper_rasters[ISO_i] = mosaic(first, adj_access_upper_rasters..., atol = 0.01)
#             adj_access_nat_lower_rasters[ISO_i] = mosaic(first, adj_access_lower_rasters..., atol = 0.01)

#         end

#         # Combine national rasters into single raster and write (Multi threaded version - Hardcoded)
#         println("Mosaic national rasters to Africa level, save and write rasters")
        
#         for thread_i in 1:6
#             if thread_i == 1
#                 println("Constructing adjusted national NPC mean raster...")
#                 combined_npc_adj_mean_raster = mosaic(first, adj_npc_nat_mean_rasters, atol = 0.01)
#                 write(output_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_mean.tif", combined_npc_adj_mean_raster, force = true)
#                 println("Constructed adjusted national NPC mean raster...")
#             elseif thread_i == 2
#                 println("Constructing adjusted national NPC upper CI raster...")
#                 combined_npc_adj_upper_raster = mosaic(first, adj_npc_nat_upper_rasters, atol = 0.01)
#                 write(output_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_upper.tif", combined_npc_adj_upper_raster, force = true)
#                 println("Constructed adjusted national NPC upper CI raster...")
#             elseif thread_i == 3
#                 println("Constructing adjusted national NPC lower CI raster...")
#                 combined_npc_adj_lower_raster = mosaic(first, adj_npc_nat_lower_rasters, atol = 0.01)
#                 write(output_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_lower.tif", combined_npc_adj_lower_raster, force = true)
#                 println("Constructed adjusted national NPC lower CI raster...")
#             elseif thread_i == 4
#                 println("Constructing adjusted national access mean raster...")
#                 combined_access_adj_mean_raster = mosaic(first, adj_access_nat_mean_rasters, atol = 0.01)
#                 write(output_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_mean.tif", combined_access_adj_mean_raster, force = true)
#                 println("Constructed adjusted national access mean raster...")
#             elseif thread_i == 5
#                 println("Constructing adjusted national access upper CI raster...")
#                 combined_access_adj_upper_raster = mosaic(first, adj_access_nat_upper_rasters, atol = 0.01)
#                 write(output_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_upper.tif", combined_access_adj_upper_raster, force = true)
#                 println("Constructed adjusted national access upper CI raster...")
#             elseif thread_i == 6
#                 println("Constructing adjusted national access lower CI raster...")
#                 combined_access_adj_lower_raster = mosaic(first, adj_access_nat_lower_rasters, atol = 0.01)
#                 write(output_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_lower.tif", combined_access_adj_lower_raster, force = true)
#                 println("Constructed adjusted national access upper CI raster...")
#             end
#         end 

#         println("Raster construction complete.")
#     end
# end

# # %% # Calculate Use Rasters and country level quantiles
# for year in ProgressBar(YEAR_START:YEAR_END, leave = false)
#     # Import population raster
#     pop_year = min(max(year, 2000), 2020)
#     population_raster = replace_missing(Raster(pop_dir*"WorldPop_UNAdj_v3_DRC_fix.$(pop_year).Annual.Data.5km.sum.tif"), missingval = NaN)

#     for month in 1:12
#         println("Processing raster year [$(year)/$(YEAR_END)], month [$(month)/12]")

#         # Calculate reference monthidx to access data
#         monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEAR_START)

#         # Get correct year/month string for importing file
#         year_str = "$(year)"
#         month_str = "$(month)"
#         if month < 10
#             month_str = "0"*month_str
#         end

#         println("Importing rasters...")
#         # Import calculated Access rasters 
#         adj_access_mean_raster = replace_missing(Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_mean.tif"), missingval = NaN)
#         adj_access_upper_raster = replace_missing(Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_upper.tif"), missingval = NaN)
#         adj_access_lower_raster = replace_missing(Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_lower.tif"), missingval = NaN)

#         # adj_access_mean_raster = replace_missing(Raster(raster_dir*"final_access/snf_access/access_$(year)_$(month_str)_mean.tif"), missingval = NaN)
#         # adj_access_upper_raster = replace_missing(Raster(raster_dir*"final_access/snf_access/access_$(year)_$(month_str)_upper.tif"), missingval = NaN)
#         # adj_access_lower_raster = replace_missing(Raster(raster_dir*"final_access/snf_access/access_$(year)_$(month_str)_lower.tif"), missingval = NaN)
        
#         # Import INLA regression of use rasters
#         # logis_use_mean_raster = resample(replace_missing(Raster(inla_dir*"inla_use_logis/USE_logismodel_$(year)_$(month)_mean.tif"), missingval = NaN), to = adj_access_mean_raster) # OLD
#         # logis_use_mean_raster = resample(replace_missing(Raster(inla_dir*"use_logis/USE_logismodel_$(year)_$(month)_mean.tif"), missingval = NaN), to = adj_access_mean_raster) # OLD
#         # logis_use_mean_raster = resample(replace_missing(Raster("outputs/INLA/rasters/inla_use_logis/USE_logismodel_$(year)_$(month)_mean.tif"), missingval = NaN), to = adj_access_mean_raster)
#         logis_use_mean_raster = resample(replace_missing(Raster("outputs/INLA/rasters/inla_adj_use_logis/adj_USE_logismodel_$(year)_$(month)_mean.tif"), missingval = NaN), to = adj_access_mean_raster)
        
#         logis_use_sample_rasters = Vector{Raster}(undef, n_samples)
#         for sample_i in 1:n_samples
#             # logis_use_sample_rasters[sample_i] = resample(replace_missing(Raster(inla_dir*"inla_use_logis/USE_logismodel_$(year)_$(month)_sample_$(sample_i).tif"), missingval = NaN), to = adj_access_mean_raster) # OLD
#             # logis_use_sample_rasters[sample_i] = resample(replace_missing(Raster(inla_dir*"use_logis/USE_logismodel_$(year)_$(month)_sample_$(sample_i).tif"), missingval = NaN), to = adj_access_mean_raster) # OLD
#             # logis_use_sample_rasters[sample_i] = resample(replace_missing(Raster("outputs/INLA/rasters/inla_use_logis/USE_logismodel_$(year)_$(month)_sample_$(sample_i).tif"), missingval = NaN), to = adj_access_mean_raster)
#             logis_use_sample_rasters[sample_i] = resample(replace_missing(Raster("outputs/INLA/rasters/inla_adj_use_logis/adj_USE_logismodel_$(year)_$(month)_sample_$(sample_i).tif"), missingval = NaN), to = adj_access_mean_raster)
#         end

#         println("Calculating use rasters for each country (with CI)")
#         # Pre-calculate continent level use rasters for mean and samples
#         use_mean_raster = inv_p_transform.(logis_use_mean_raster, adj_access_mean_raster, n=2)

#         use_sample_upper_rasters = Vector{Raster}(undef, n_samples)
#         use_sample_lower_rasters = Vector{Raster}(undef, n_samples)
#         for sample_i in 1:n_samples
#             use_sample_upper_rasters[sample_i] = inv_p_transform.(logis_use_sample_rasters[sample_i], adj_access_upper_raster, n=2)
#             use_sample_lower_rasters[sample_i] = inv_p_transform.(logis_use_sample_rasters[sample_i], adj_access_lower_raster, n=2)
#         end

#         # Storage variables for country rasters of ITN coverage with CI
#         use_nat_upper_rasters = Vector{Raster}(undef, length(filt_ISOs))
#         use_nat_lower_rasters = Vector{Raster}(undef, length(filt_ISOs))

#         # Do extractions and calculations for each country
#         for ISO_i in ProgressBar(1:length(filt_ISOs), leave = false)
#             # Get ISO
#             ISO = filt_ISOs[ISO_i]

#             # Get required national geometry
#             admin0_geometry = admin0_shapes_geoIO[admin0_shapes_geoIO.ISO .== ISO,:].geometry

#             # Get trimmed subset of population raster for target country
#             pop_nat_masked = Rasters.trim(mask(population_raster, with = admin0_geometry); pad=0)

#             # Get trimmed subset of use rasters for target country
#             # This line is not needed, but really just used to make code more readable in population weighting line
#             use_mean_raster_masked = resample(Rasters.trim(mask(use_mean_raster, with = admin0_geometry); pad=0), to = pop_nat_masked)

#             use_sample_upper_rasters_masked = Vector{Raster}(undef, n_samples)
#             use_sample_lower_rasters_masked = Vector{Raster}(undef, n_samples)
#             for sample_i in 1:n_samples
#                 use_sample_upper_rasters_masked[sample_i] = resample(Rasters.trim(mask(use_sample_upper_rasters[sample_i], with = admin0_geometry); pad=0), to = pop_nat_masked)
#                 use_sample_lower_rasters_masked[sample_i] = resample(Rasters.trim(mask(use_sample_lower_rasters[sample_i], with = admin0_geometry); pad=0), to = pop_nat_masked)
#             end

#             # Calculate the population weighted average the upper and lower use rasters to get joint CI

#             # Find non-missing entries in raster
#             nonmissing_idxs = intersect(findall(.!isnan.(pop_nat_masked)),findall(.!isnan.(use_mean_raster_masked)))

#             # Calculated weighted means
#             use_rates_samples_upper = Vector{Float64}(undef, n_samples)
#             use_rates_samples_lower = Vector{Float64}(undef, n_samples)
#             for sample_i in 1:n_samples
#                 use_rates_samples_upper[sample_i] = sum(pop_nat_masked[nonmissing_idxs].*use_sample_upper_rasters_masked[sample_i][nonmissing_idxs])/sum(pop_nat_masked[nonmissing_idxs])
#                 use_rates_samples_lower[sample_i] = sum(pop_nat_masked[nonmissing_idxs].*use_sample_lower_rasters_masked[sample_i][nonmissing_idxs])/sum(pop_nat_masked[nonmissing_idxs])
#             end

#             # Get the quantile rasters based on population weighted means on national level
#             upper_rank_idx = round(Int,0.95*n_samples)
#             lower_rank_idx = max(round(Int, 0.05*n_samples),1)

#             use_nat_upper_rasters[ISO_i] = use_sample_upper_rasters_masked[sortperm(use_rates_samples_upper)][upper_rank_idx]
#             use_nat_lower_rasters[ISO_i] = use_sample_lower_rasters_masked[sortperm(use_rates_samples_lower)][lower_rank_idx]
#         end

#         # Mosaic together the upper and lower national rasters
#         use_upper_raster = resample(mosaic(first, use_nat_upper_rasters, atol = 0.01), to = use_mean_raster)
#         use_lower_raster = resample(mosaic(first, use_nat_lower_rasters, atol = 0.01), to = use_mean_raster)
        
#         # # Write rasters
#         println("Saving rasters...")
#         write(output_dir*"final_use/logis_use/use_$(year)_$(month_str)_mean.tif", use_mean_raster, force = true)
#         write(output_dir*"final_use/logis_use/use_$(year)_$(month_str)_upper.tif", use_upper_raster, force = true)
#         write(output_dir*"final_use/logis_use/use_$(year)_$(month_str)_lower.tif", use_lower_raster, force = true)

#         println("Raster construction complete.")
#     end
# end