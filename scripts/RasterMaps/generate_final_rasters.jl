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
snf_post_dir = "outputs/draws/subnational/"

# Region Admin 1 area id legend
admin1_legend_dir = "datasets/subnational/"
admin1_legend_filename = "admin2023_1_MG_5K_config.csv"

# Population Rasters directory
pop_dir = "Z:/mastergrids/Other_Global_Covariates/Population/WorldPop/v3/PopulationCounts_DRC_fixed/5km/"

# Region boundaries
admin0_shapes_geoIO = GeoIO.load("Z:/master_geometries/Admin_Units/Global/MAP/2023/MG_5K/admin2023_0_MG_5K.shp")
admin1_shapes_geoIO = GeoIO.load("Z:/master_geometries/Admin_Units/Global/MAP/2023/MG_5K/admin2023_1_MG_5K.shp")

# Num of samples to import from INLA outputs for use
n_samples = 10 

# Get base raster with required resolution to build from
# raster_base = replace_missing(Raster("outputs/rasters/inla_logmodel_npc/NPC_logmodel_$(2000)_mean.tif"), missingval = -NaN)

# Input and output directory for rasters
inla_dir = "Z:/eugene/INLA Outputs/"
raster_dir = "Z:/eugene/Final Rasters/rasters/"#"outputs/rasters/"
input_dir = "Z:/eugene/Final Rasters/rasters/"#"outputs/rasters/"
output_dir = "Z:/eugene/Final Rasters/rasters/"#"outputs/rasters/"

# Make paths tosave location if doesn't already exist
mkpath(output_dir*"final_npc/snf_npc/")
mkpath(output_dir*"final_npc/logmodel_npc/")
mkpath(output_dir*"final_access/snf_access/")
mkpath(output_dir*"final_access/pmodel_access/")
mkpath(output_dir*"final_use/logis_use/")

# %% Perform draws and save outputs. Filter out unwanted countries
ISO_list = String.(CSV.read("/datasets/ISO_list.csv", DataFrame)[:,1])
exclusion_ISOs = ["CPV","ZAF"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Time bounds
YEAR_START = 2000
YEAR_END = 2023

# %% Loop to first construct SNF block map upto subnational resolution
# Import log model npc rasters (any calibrated raster will do 5x5km resolution)
raster_base = replace_missing(Raster("Z:/eugene/INLA Outputs/inla_logmodel_npc/NPC_logmodel_$(2000)_mean.tif"), missingval = NaN)

# %% Loop to construct raster
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


# %% Loop to construct unadjusted NPC and Access rasters using INLA deviation regression outputs
for year in ProgressBar(YEAR_START:YEAR_END)
    for month in 1:12
        println("Constructing spatial disaggregated rasters year [$(year)/$(YEAR_END)], month [$(month)/12]")

        # Get month string for importing files
        month_str = "$(month)"
        if month < 10
            month_str = "0$(month)"
        end

        # Get monthidx for SNF lookup
        monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEAR_START)
        
        println("Importing rasters...")
        # Import stock and flow rasters
        npc_snf_mean_raster = replace_missing(Raster(input_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        npc_snf_upper_raster = replace_missing(Raster(input_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_upper.tif"), missingval = NaN)
        npc_snf_lower_raster = replace_missing(Raster(input_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_lower.tif"), missingval = NaN)
        access_snf_mean_raster = replace_missing(Raster(input_dir*"final_access/snf_access/access_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        access_snf_upper_raster = replace_missing(Raster(input_dir*"final_access/snf_access/access_$(year)_$(month_str)_upper.tif"), missingval = NaN)
        access_snf_lower_raster = replace_missing(Raster(input_dir*"final_access/snf_access/access_$(year)_$(month_str)_lower.tif"), missingval = NaN)

        # Import log model npc rasters
        logmodel_npc_mean_raster = replace_missing(Raster(inla_dir*"inla_logmodel_npc/NPC_logmodel_$(year)_mean.tif"), missingval = NaN)

        # Import pmodel access rasters
        pmodel_access_mean_raster =  replace_missing(Raster(inla_dir*"inla_pmodel_access/ACCESS_pmodel_$(year)_mean.tif"), missingval = NaN)

        # TEMP - Used just to make sure rasters are all aligned
        raster_base = copy(logmodel_npc_mean_raster)

        
        # Calculate spatial disaggregated rasters
        # Apply logmodel ratios to get estimate of npc map
        println("Calculating NPC rasters...")
        logmodel_npc_ratio_mean_raster = exp.(logmodel_npc_mean_raster)

        npc_map_mean_raster = logmodel_npc_ratio_mean_raster.* npc_snf_mean_raster
        npc_map_upper_raster = logmodel_npc_ratio_mean_raster .* npc_snf_upper_raster
        npc_map_lower_raster = logmodel_npc_ratio_mean_raster .* npc_snf_lower_raster

        # Inverse transform inla estimates to get access map estimates
        println("Calculating Access rasters...")
        access_map_mean_raster = inv_p_transform.(pmodel_access_mean_raster, access_snf_mean_raster, n=2)
        access_map_upper_raster = inv_p_transform.(pmodel_access_mean_raster, access_snf_upper_raster, n=2)
        access_map_lower_raster = inv_p_transform.(pmodel_access_mean_raster, access_snf_lower_raster, n=2)

        # Save and write rasters
        println("Save and write rasters")
        write(output_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_mean.tif", npc_map_mean_raster, force = true)
        write(output_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_upper.tif", npc_map_upper_raster, force = true)
        write(output_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_lower.tif", npc_map_lower_raster, force = true)

        write(output_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_mean.tif", access_map_mean_raster, force = true)
        write(output_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_upper.tif", access_map_upper_raster, force = true)
        write(output_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_lower.tif", access_map_lower_raster, force = true)

        println("Raster construction complete.")
    end
end

# %% # Calculate adjusted NPC and Access Rasters (Raking on a country level)
# Subnat estimates of NPC should be close to national estimates when summed by population (SNF Calibration step)
# Subnat estimates of access may be considered unreliable for raking as it assumes household demographic distribution to be consistent across whole country
for year in ProgressBar(YEAR_START:YEAR_END, leave = false)
    # Import population raster
    pop_year = min(max(year, 2000), 2020)
    population_raster = replace_missing(Raster(pop_dir*"WorldPop_UNAdj_v3_DRC_fix.$(pop_year).Annual.Data.5km.sum.tif"), missingval = NaN)

    for month in 1:12
        println("Processing raster year [$(year)/$(YEAR_END)], month [$(month)/12]")

        # Calculate reference monthidx to access data
        monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEAR_START)

        # Get correct year/month string for importing file
        year_str = "$(year)"
        month_str = "$(month)"
        if month < 10
            month_str = "0"*month_str
        end

        println("Importing rasters...")
        # Import calculated NPC rasters
        npc_mean_raster = Raster(input_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_mean.tif")
        npc_upper_raster = Raster(input_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_upper.tif")
        npc_lower_raster = Raster(input_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_lower.tif")
        
        # Import calculated Access rasters 
        access_mean_raster = Raster(input_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_mean.tif")
        access_upper_raster = Raster(input_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_upper.tif")
        access_lower_raster = Raster(input_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_lower.tif")

    
        # Storage variables for country rasters of ITN coverage with CI
        adj_npc_nat_mean_rasters = Vector{Raster}(undef, length(filt_ISOs))
        adj_npc_nat_upper_rasters = Vector{Raster}(undef, length(filt_ISOs))
        adj_npc_nat_lower_rasters = Vector{Raster}(undef, length(filt_ISOs))

        adj_access_nat_mean_rasters = Vector{Raster}(undef, length(filt_ISOs))
        adj_access_nat_upper_rasters = Vector{Raster}(undef, length(filt_ISOs))
        adj_access_nat_lower_rasters = Vector{Raster}(undef, length(filt_ISOs))

        println("Raking rasters using SNF National estimates...")
        # Do extractions and calculations for each country
        for ISO_i in ProgressBar(1:length(filt_ISOs), leave = false)
            ISO = filt_ISOs[ISO_i]

            # Extract national SNF draws to rake rasters
            nat_snf_post_draws = load("outputs/draws/national/crop_access/$(ISO)_2000_2023_post_crop_access.jld2")
            nat_population = nat_snf_post_draws["POPULATION_MONTHLY"][monthidx]
            nat_npc_draws = sum(nat_snf_post_draws["Γ_MONTHLY_samples_BYNET"], dims = 3)[:,monthidx,1]./nat_population
            nat_access_draws = nat_snf_post_draws["λ_access_samples"][:, monthidx]

            # Calculate npc and access estimates to rake country raster by
            npc_snf_mean_estimate = mean(nat_npc_draws)
            npc_snf_upper_estimate = quantile(nat_npc_draws, 0.95)
            npc_snf_lower_estimate = quantile(nat_npc_draws, 0.05)

            access_snf_mean_estimate = mean(nat_access_draws)
            access_snf_upper_estimate = quantile(nat_access_draws, 0.95)
            access_snf_lower_estimate = quantile(nat_access_draws, 0.05)

            # Get Country geometry to mask raster
            admin0_geometry = admin0_shapes_geoIO[admin0_shapes_geoIO.ISO .== ISO,:].geometry

            # Get masked + trimmed versions of each required component raster
            ## Population
            pop_masked = Rasters.trim(mask(population_raster, with = admin0_geometry); pad=0)

            ## NPC
            npc_mean_masked = resample(Rasters.trim(mask(npc_mean_raster, with = admin0_geometry); pad=0), to = pop_masked)
            npc_upper_masked = resample(Rasters.trim(mask(npc_upper_raster, with = admin0_geometry); pad=0), to = pop_masked)
            npc_lower_masked = resample(Rasters.trim(mask(npc_lower_raster, with = admin0_geometry); pad=0), to = pop_masked)

            ## Access
            access_mean_masked = resample(Rasters.trim(mask(access_mean_raster, with = admin0_geometry); pad=0), to = pop_masked)
            access_upper_masked = resample(Rasters.trim(mask(access_upper_raster, with = admin0_geometry); pad=0), to = pop_masked)
            access_lower_masked = resample(Rasters.trim(mask(access_lower_raster, with = admin0_geometry); pad=0), to = pop_masked)

            #####
            # NPC Raster Calculations
            #####

            # TEMP FIX: Calculate separate nonmissing_idxs for npc, access and use (INLA transformation function didn't deal with access = 0 and 1 very well)
            # Calculate population weighted mean npc for subnational region
            nonmissing_idx_npc_mean = intersect(findall(.!isnan.(pop_masked)), findall(.!isnan.(npc_mean_masked)))
            nonmissing_idx_npc_upper = intersect(findall(.!isnan.(pop_masked)), findall(.!isnan.(npc_upper_masked)))
            nonmissing_idx_npc_lower = intersect(findall(.!isnan.(pop_masked)), findall(.!isnan.(npc_lower_masked)))
            npc_spatial_mean_estimate = sum(pop_masked[nonmissing_idx_npc_mean].*npc_mean_masked[nonmissing_idx_npc_mean])/sum(pop_masked[nonmissing_idx_npc_mean])
            npc_spatial_upper_estimate = sum(pop_masked[nonmissing_idx_npc_upper].*npc_upper_masked[nonmissing_idx_npc_upper])/sum(pop_masked[nonmissing_idx_npc_upper])
            npc_spatial_lower_estimate = sum(pop_masked[nonmissing_idx_npc_lower].*npc_lower_masked[nonmissing_idx_npc_lower])/sum(pop_masked[nonmissing_idx_npc_lower])

            # Calculate scaling constant required to adjust map estimates of NPC to match SNF estimates of NPC
            npc_mean_scaling_k = npc_snf_mean_estimate/npc_spatial_mean_estimate
            npc_upper_scaling_k = npc_snf_upper_estimate/npc_spatial_upper_estimate
            npc_lower_scaling_k = npc_snf_lower_estimate/npc_spatial_lower_estimate

            # Calculate adjusted maps/rasters using the scaling constant
            adj_npc_mean_masked = npc_mean_scaling_k.*npc_mean_masked
            adj_npc_upper_masked = npc_upper_scaling_k.*npc_upper_masked
            adj_npc_lower_masked = npc_lower_scaling_k.*npc_lower_masked

            #####
            # Access Raster Calculations
            #####
            # TEMP FIX
            nonmissing_idx_access_mean = intersect(findall(.!isnan.(pop_masked)), findall(.!isnan.(access_mean_masked)))
            nonmissing_idx_access_upper = intersect(findall(.!isnan.(pop_masked)), findall(.!isnan.(access_upper_masked)))
            nonmissing_idx_access_lower = intersect(findall(.!isnan.(pop_masked)), findall(.!isnan.(access_lower_masked)))

            # Calculate population weighted mean access for subnational region
            access_spatial_mean_estimate = sum(pop_masked[nonmissing_idx_access_mean].*access_mean_masked[nonmissing_idx_access_mean])/sum(pop_masked[nonmissing_idx_access_mean])
            access_spatial_upper_estimate = sum(pop_masked[nonmissing_idx_access_upper].*access_upper_masked[nonmissing_idx_access_upper])/sum(pop_masked[nonmissing_idx_access_upper])
            access_spatial_lower_estimate = sum(pop_masked[nonmissing_idx_access_lower].*access_lower_masked[nonmissing_idx_access_lower])/sum(pop_masked[nonmissing_idx_access_lower])
            
            # Calculate scaling constant required to adjust map estimates of NPC to match SNF estimates of access
            access_mean_scaling_k = access_snf_mean_estimate/access_spatial_mean_estimate
            access_upper_scaling_k = access_snf_upper_estimate/access_spatial_upper_estimate
            access_lower_scaling_k = access_snf_lower_estimate/access_spatial_lower_estimate

            # Calcualte adjusted maps/rasters using the scaling constant + appropriate truncation
            adj_access_mean_masked = max.(min.(access_mean_scaling_k.*access_mean_masked, 1), 0)
            adj_access_upper_masked = max.(min.(access_upper_scaling_k.*access_upper_masked, 1), 0)
            adj_access_lower_masked = max.(min.(access_lower_scaling_k.*access_lower_masked, 1), 0)

            ### Save to storage variables
            adj_npc_nat_mean_rasters[ISO_i] = adj_npc_mean_masked
            adj_npc_nat_upper_rasters[ISO_i] = adj_npc_upper_masked
            adj_npc_nat_lower_rasters[ISO_i] = adj_npc_lower_masked

            adj_access_nat_mean_rasters[ISO_i] = adj_access_mean_masked
            adj_access_nat_upper_rasters[ISO_i] = adj_access_upper_masked
            adj_access_nat_lower_rasters[ISO_i] = adj_access_lower_masked
        end
        
        # Combine national rasters into single raster and write (Multi threaded version - Hardcoded)
        println("Mosaic national rasters to Africa level, save and write rasters")
        for thread_i in 1:6
            if thread_i == 1
                println("Constructing adjusted national NPC mean raster...")
                combined_npc_adj_mean_raster = mosaic(first, adj_npc_nat_mean_rasters, atol = 0.01)
                write(output_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_mean.tif", combined_npc_adj_mean_raster, force = true)
                println("Constructed adjusted national NPC mean raster...")
            elseif thread_i == 2
                println("Constructing adjusted national NPC upper CI raster...")
                combined_npc_adj_upper_raster = mosaic(first, adj_npc_nat_upper_rasters, atol = 0.01)
                write(output_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_upper.tif", combined_npc_adj_upper_raster, force = true)
                println("Constructed adjusted national NPC upper CI raster...")
            elseif thread_i == 3
                println("Constructing adjusted national NPC lower CI raster...")
                combined_npc_adj_lower_raster = mosaic(first, adj_npc_nat_lower_rasters, atol = 0.01)
                write(output_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_lower.tif", combined_npc_adj_lower_raster, force = true)
                println("Constructed adjusted national NPC lower CI raster...")
            elseif thread_i == 4
                println("Constructing adjusted national access mean raster...")
                combined_access_adj_mean_raster = mosaic(first, adj_access_nat_mean_rasters, atol = 0.01)
                write(output_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_mean.tif", combined_access_adj_mean_raster, force = true)
                println("Constructed adjusted national access mean raster...")
            elseif thread_i == 5
                println("Constructing adjusted national access upper CI raster...")
                combined_access_adj_upper_raster = mosaic(first, adj_access_nat_upper_rasters, atol = 0.01)
                write(output_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_upper.tif", combined_access_adj_upper_raster, force = true)
                println("Constructed adjusted national access upper CI raster...")
            elseif thread_i == 6
                println("Constructing adjusted national access lower CI raster...")
                combined_access_adj_lower_raster = mosaic(first, adj_access_nat_lower_rasters, atol = 0.01)
                write(output_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_lower.tif", combined_access_adj_lower_raster, force = true)
                println("Constructed adjusted national access upper CI raster...")
            end
        end

        println("Raster construction complete.")
    end
end

# %% # Calculate Use Rasters and country level quantiles
for year in ProgressBar(YEAR_START:YEAR_END, leave = false)
    # Import population raster
    pop_year = min(max(year, 2000), 2020)
    population_raster = replace_missing(Raster(pop_dir*"WorldPop_UNAdj_v3_DRC_fix.$(pop_year).Annual.Data.5km.sum.tif"), missingval = NaN)

    for month in 1:12
        println("Processing raster year [$(year)/$(YEAR_END)], month [$(month)/12]")

        # Calculate reference monthidx to access data
        monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEAR_START)

        # Get correct year/month string for importing file
        year_str = "$(year)"
        month_str = "$(month)"
        if month < 10
            month_str = "0"*month_str
        end

        println("Importing rasters...")
        # Import calculated Access rasters 
        adj_access_mean_raster = Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_mean.tif")
        adj_access_upper_raster = Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_upper.tif")
        adj_access_lower_raster = Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_lower.tif")
        
        # Import INLA regression of use rasters
        # logis_use_mean_raster = resample(replace_missing(Raster(inla_dir*"inla_use_logis/USE_logismodel_$(year)_$(month)_mean.tif"), missingval = NaN), to = adj_access_mean_raster)
        # logis_use_mean_raster = resample(replace_missing(Raster(inla_dir*"use_logis/USE_logismodel_$(year)_$(month)_mean.tif"), missingval = NaN), to = adj_access_mean_raster)
        logis_use_mean_raster = resample(replace_missing(Raster("outputs/rasters/use_logis/USE_logismodel_$(year)_$(month)_mean.tif"), missingval = NaN), to = adj_access_mean_raster)
        
        logis_use_sample_rasters = Vector{Raster}(undef, n_samples)
        for sample_i in 1:n_samples
            # logis_use_sample_rasters[sample_i] = resample(replace_missing(Raster(inla_dir*"inla_use_logis/USE_logismodel_$(year)_$(month)_sample_$(sample_i).tif"), missingval = NaN), to = adj_access_mean_raster)
            # logis_use_sample_rasters[sample_i] = resample(replace_missing(Raster(inla_dir*"use_logis/USE_logismodel_$(year)_$(month)_sample_$(sample_i).tif"), missingval = NaN), to = adj_access_mean_raster)
            logis_use_sample_rasters[sample_i] = resample(replace_missing(Raster("outputs/rasters/use_logis/USE_logismodel_$(year)_$(month)_sample_$(sample_i).tif"), missingval = NaN), to = adj_access_mean_raster)
        end

        println("Calculating use rasters for each country (with CI)")
        # Pre-calculate continent level use rasters for mean and samples
        use_mean_raster = inv_p_transform.(logis_use_mean_raster, adj_access_mean_raster, n=2)

        use_sample_upper_rasters = Vector{Raster}(undef, n_samples)
        use_sample_lower_rasters = Vector{Raster}(undef, n_samples)
        for sample_i in 1:n_samples
            use_sample_upper_rasters[sample_i] = inv_p_transform.(logis_use_sample_rasters[sample_i], adj_access_upper_raster, n=2)
            use_sample_lower_rasters[sample_i] = inv_p_transform.(logis_use_sample_rasters[sample_i], adj_access_lower_raster, n=2)
        end

        # Storage variables for country rasters of ITN coverage with CI
        use_nat_upper_rasters = Vector{Raster}(undef, length(filt_ISOs))
        use_nat_lower_rasters = Vector{Raster}(undef, length(filt_ISOs))

        # Do extractions and calculations for each country
        for ISO_i in ProgressBar(1:length(filt_ISOs), leave = false)
            # Get ISO
            ISO = filt_ISOs[ISO_i]

            # Get required national geometry
            admin0_geometry = admin0_shapes_geoIO[admin0_shapes_geoIO.ISO .== ISO,:].geometry

            # Get trimmed subset of population raster for target country
            pop_nat_masked = Rasters.trim(mask(population_raster, with = admin0_geometry); pad=0)

            # Get trimmed subset of use rasters for target country
            # This line is not needed, but really just used to make code more readable in population weighting line
            use_mean_raster_masked = resample(Rasters.trim(mask(use_mean_raster, with = admin0_geometry); pad=0), to = pop_nat_masked)

            use_sample_upper_rasters_masked = Vector{Raster}(undef, n_samples)
            use_sample_lower_rasters_masked = Vector{Raster}(undef, n_samples)
            for sample_i in 1:n_samples
                use_sample_upper_rasters_masked[sample_i] = resample(Rasters.trim(mask(use_sample_upper_rasters[sample_i], with = admin0_geometry); pad=0), to = pop_nat_masked)
                use_sample_lower_rasters_masked[sample_i] = resample(Rasters.trim(mask(use_sample_lower_rasters[sample_i], with = admin0_geometry); pad=0), to = pop_nat_masked)
            end

            # Calculate the population weighted average the upper and lower use rasters to get joint CI

            # Find non-missing entries in raster
            nonmissing_idxs = intersect(findall(.!isnan.(pop_nat_masked)),findall(.!isnan.(use_mean_raster_masked)))

            # Calculated weighted means
            use_rates_samples_upper = Vector{Float64}(undef, n_samples)
            use_rates_samples_lower = Vector{Float64}(undef, n_samples)
            for sample_i in 1:n_samples
                use_rates_samples_upper[sample_i] = sum(pop_nat_masked[nonmissing_idxs].*use_sample_upper_rasters_masked[sample_i][nonmissing_idxs])/sum(pop_nat_masked[nonmissing_idxs])
                use_rates_samples_lower[sample_i] = sum(pop_nat_masked[nonmissing_idxs].*use_sample_lower_rasters_masked[sample_i][nonmissing_idxs])/sum(pop_nat_masked[nonmissing_idxs])
            end

            # Get the quantile rasters based on population weighted means on national level
            upper_rank_idx = round(Int,0.95*n_samples)
            lower_rank_idx = max(round(Int, 0.05*n_samples),1)

            use_nat_upper_rasters[ISO_i] = use_sample_upper_rasters_masked[sortperm(use_rates_samples_upper)][upper_rank_idx]
            use_nat_lower_rasters[ISO_i] = use_sample_lower_rasters_masked[sortperm(use_rates_samples_lower)][lower_rank_idx]
        end

        # Mosaic together the upper and lower national rasters
        use_upper_raster = resample(mosaic(first, use_nat_upper_rasters, atol = 0.01), to = use_mean_raster)
        use_lower_raster = resample(mosaic(first, use_nat_lower_rasters, atol = 0.01), to = use_mean_raster)
        
        # # Write rasters
        println("Saving rasters...")
        write(output_dir*"final_use/logis_use/use_$(year)_$(month_str)_mean.tif", use_mean_raster, force = true)
        write(output_dir*"final_use/logis_use/use_$(year)_$(month_str)_upper.tif", use_upper_raster, force = true)
        write(output_dir*"final_use/logis_use/use_$(year)_$(month_str)_lower.tif", use_lower_raster, force = true)

        println("Raster construction complete.")
    end
end