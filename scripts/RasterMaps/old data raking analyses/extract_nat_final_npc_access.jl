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
using CSV

using Plots
# Custom packages
using DateConversions

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

# Num of samples to import from INLA outputs for use
n_samples = 20

# Get base raster with required resolution to build from
# raster_base = replace_missing(Raster("outputs/rasters/inla_logmodel_npc/NPC_logmodel_$(2000)_mean.tif"), missingval = -NaN)

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
exclusion_ISOs = ["CPV","ZAF"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Time bounds
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

# %%

# %% Construct final NPC and Access rasters from final INLA estimates


# %% Loop to construct unadjusted NPC and Access rasters using INLA deviation regression outputs
# for year in ProgressBar(YEAR_START:YEAR_END)
    year = 2015
    

    # Import log model npc rasters
    logmodel_npc_mean_raster = replace_missing(Raster(inla_dir*"inla_logmodel_final_npc/final_NPC_logmodel_$(year)_mean.tif"), missingval = NaN)
    logmodel_npc_sample_rasters = Vector{Any}(undef, n_samples)
    for sample_i in 1:n_samples
        logmodel_npc_sample_rasters[sample_i] = replace_missing(Raster(inla_dir*"inla_logmodel_final_npc/final_NPC_pmodel_$(year)_sample_$(sample_i).tif"), missingval = NaN)
    end

    # Import pmodel access rasters
    pmodel_access_mean_raster = replace_missing(Raster(inla_dir*"inla_pmodel_final_access/final_ACCESS_pmodel_$(year)_mean.tif"), missingval = NaN)
    pmodel_access_sample_rasters = Vector{Any}(undef, n_samples)
    for sample_i in 1:n_samples
        pmodel_access_sample_rasters[sample_i] = replace_missing(Raster(inla_dir*"inla_pmodel_final_access/final_ACCESS_pmodel_$(year)_sample_$(sample_i).tif"), missingval = NaN)
    end

    # Import population raster
    pop_year = min(max(year, 2000), 2020)
    population_raster = replace_missing(Raster(pop_dir*"WorldPop_UNAdj_v3_DRC_fix.$(pop_year).Annual.Data.5km.sum.tif"), missingval = NaN)
    # Resample and align population raster to ITN Coverage rasters
    population_raster = resample(population_raster, to = logmodel_npc_mean_raster)

#     for month in 1:12

        month = 1
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
        npc_snf_rasters = [npc_snf_lower_raster, npc_snf_mean_raster, npc_snf_upper_raster]

        access_snf_mean_raster = replace_missing(Raster(input_dir*"final_access/snf_access/access_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        access_snf_upper_raster = replace_missing(Raster(input_dir*"final_access/snf_access/access_$(year)_$(month_str)_upper.tif"), missingval = NaN)
        access_snf_lower_raster = replace_missing(Raster(input_dir*"final_access/snf_access/access_$(year)_$(month_str)_lower.tif"), missingval = NaN)
        access_snf_rasters = [access_snf_lower_raster, access_snf_mean_raster, access_snf_upper_raster]
        
        logmodel_unadj_npc_mean_raster = replace_missing(Raster(input_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        unadj_mean_raster = exp.(logmodel_unadj_npc_mean_raster) .*npc_snf_mean_raster

        pmodel_unadj_access_mean_raster = replace_missing(Raster(input_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_mean.tif"), missingval = NaN)
        
        
        # Calculate NPC rasters
        # Apply logmodel ratios to get estimate of npc map
        println("Calculating NPC rasters...")

        # Mean rasters
        logmodel_npc_ratio_mean_raster = exp.(logmodel_npc_mean_raster)
        npc_mean_raster = logmodel_npc_ratio_mean_raster.* npc_snf_mean_raster

        # Sample rasters
        npc_sample_rasters = Array{Any}(undef, 3, n_samples)
        for quantile_j in 1:3 # lower, mean, upper
            for sample_i in 1:n_samples
                npc_sample_rasters[quantile_j, sample_i] = exp.(logmodel_npc_sample_rasters[sample_i]) .* npc_snf_rasters[quantile_j]
            end
        end

        # Calculate Access rasters
        println("Calculating Access rasters...")

        # Mean rasters
        access_map_mean_raster = inv_p_transform.(pmodel_access_mean_raster, access_snf_mean_raster, n=2)

        # Sample rasters
        access_sample_rasters = Array{Any}(undef, 3, n_samples)
        for quantile_j in 1:3 # lower, mean, upper
            for sample_i in 1:n_samples
                access_sample_rasters[quantile_j, sample_i] = inv_p_transform.(pmodel_access_sample_rasters[sample_i], access_snf_rasters[quantile_j], n=2)
            end
        end

        # Get upper and lower population weighted means for each country, and associated rasters
        npc_iso_lower_rasters = Vector{Any}(undef, length(filt_ISOs))
        npc_iso_mean_rasters = Vector{Any}(undef, length(filt_ISOs))
        npc_iso_upper_rasters = Vector{Any}(undef, length(filt_ISOs))

        access_iso_lower_rasters = Vector{Any}(undef, length(filt_ISOs))
        access_iso_mean_rasters = Vector{Any}(undef, length(filt_ISOs))
        access_iso_upper_rasters = Vector{Any}(undef, length(filt_ISOs))

        ISO_i = 10
        ISO = filt_ISOs[ISO_i]
        ISO = "TZA"
        admin0_geometry = admin0_shapes_geoIO[admin0_shapes_geoIO.ISO .== ISO,"geometry"]

        # Get trimmed subset of population raster for target country
        pop_nat_masked = Rasters.trim(mask(population_raster, with = admin0_geometry); pad=0)

        # Get trimmed subset of ITN coverage rasters and align with population raster
        npc_mean_masked = resample(Rasters.trim(mask(npc_mean_raster, with = admin0_geometry); pad=0), to = pop_nat_masked)
		npc_sample_rasters_masked = Array{Any}(undef, 3, n_samples)
        npc_sample_values = Array{Float64}(undef, 3, n_samples)
        for quantile_j in 1:3 # lower, mean, upper
            for sample_i in 1:n_samples
                npc_sample_rasters_masked[quantile_j, sample_i] = resample(Rasters.trim(mask(npc_sample_rasters[quantile_j, sample_i], with = admin0_geometry); pad=0), to = pop_nat_masked)
                nonmissing_idxs_npc = intersect(findall(.!isnan.(pop_nat_masked)), findall(.!isnan.(npc_sample_rasters_masked[quantile_j, sample_i])))
                npc_sample_values[quantile_j, sample_i] = sum((npc_sample_rasters_masked[quantile_j, sample_i].*pop_nat_masked)[nonmissing_idxs_npc])/sum(pop_nat_masked[nonmissing_idxs_npc])
            end
        end

        access_mean_masked = resample(Rasters.trim(mask(access_mean_raster, with = admin0_geometry); pad=0), to = pop_nat_masked)
		access_sample_rasters_masked = Array{Any}(undef, 3, n_samples)
        access_sample_values = Array{Float64}(undef, 3, n_samples)
        for quantile_j in 1:3 # lower, mean, upper
            for sample_i in 1:n_samples
                access_sample_rasters_masked[quantile_j, sample_i] = resample(Rasters.trim(mask(access_sample_rasters[quantile_j, sample_i], with = admin0_geometry); pad=0), to = pop_nat_masked)
                nonmissing_idxs_npc = intersect(findall(.!isnan.(pop_nat_masked)), findall(.!isnan.(access_sample_rasters_masked[quantile_j, sample_i])))
                access_sample_values[quantile_j, sample_i] = sum((access_sample_rasters_masked[quantile_j, sample_i].*pop_nat_masked)[nonmissing_idxs_npc])/sum(pop_nat_masked[nonmissing_idxs_npc])
            end
        end

        # Calculate NPC 95% CI rasters
        npc_lower_val, npc_upper_val = quantile(npc_sample_values[:], [0.025, 0.975])
        upper_rank_idx = round(Int,0.975*n_samples)
        lower_rank_idx = max(round(Int, 0.025*n_samples),1)

        npc_lower_raster = npc_sample_rasters_masked[lower_rank_idx]
        npc_upper_raster = npc_sample_rasters_masked[upper_rank_idx]

        # Save rasters to collection
        npc_iso_lower_rasters[ISO_i] = npc_lower_raster
        npc_iso_mean_rasters[ISO_i] = npc_mean_masked
        npc_iso_upper_rasters[ISO_i] = npc_upper_raster
        
        
        # npc_snf_mean_masked = resample(Rasters.trim(mask(npc_snf_mean_raster, with = admin0_geometry); pad=0), to = pop_nat_masked)
        # unadj_npc_mean_masked = resample(Rasters.trim(mask(unadj_mean_raster, with = admin0_geometry); pad=0), to = pop_nat_masked)
        
        # nonmissing_idxs_npc = intersect(findall(.!isnan.(pop_nat_masked)), findall(.!isnan.(npc_mean_masked)))
        
        # npc_snf = sum((npc_snf_mean_masked.*pop_nat_masked)[nonmissing_idxs_npc])/sum(pop_nat_masked[nonmissing_idxs_npc])
        # npc_unadj = sum((unadj_npc_mean_masked.*pop_nat_masked)[nonmissing_idxs_npc])/sum(pop_nat_masked[nonmissing_idxs_npc])
        # npc_adj = sum((npc_mean_masked.*pop_nat_masked)[nonmissing_idxs_npc])/sum(pop_nat_masked[nonmissing_idxs_npc])
        # plot(resample(Rasters.trim(mask(logmodel_npc_ratio_mean_raster, with = admin0_geometry); pad=0), to = pop_nat_masked),clims = (0,2))
        # plot(plot(npc_snf_mean_masked, clims = (0,1)), 
        #     plot(unadj_npc_mean_masked, clims = (0,1)), 
        #     plot(npc_mean_masked, clims = (0,1)), layout = (1,3), size = (800,300))
        
        
        # npc_snf/npc_unadj
        # npc_snf/npc_adj

        # npc_snf-npc_unadj
        # npc_snf-npc_adj

        
        
        # plot(npc_sample_rasters_masked[2,2], clims = (0,0.5))
        # plot(npc_mean_masked, clims = (0,0.5))





        # # Inverse transform inla estimates to get access map estimates
        # println("Calculating Access rasters...")
        # access_map_mean_raster = inv_p_transform.(pmodel_access_mean_raster, access_snf_mean_raster, n=2)
        # access_map_upper_raster = inv_p_transform.(pmodel_access_mean_raster, access_snf_upper_raster, n=2)
        # access_map_lower_raster = inv_p_transform.(pmodel_access_mean_raster, access_snf_lower_raster, n=2)

        # # Save and write rasters
        # println("Save and write rasters")
        # write(output_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_mean.tif", npc_map_mean_raster, force = true)
        # write(output_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_upper.tif", npc_map_upper_raster, force = true)
        # write(output_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_lower.tif", npc_map_lower_raster, force = true)

        # write(output_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_mean.tif", access_map_mean_raster, force = true)
        # write(output_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_upper.tif", access_map_upper_raster, force = true)
        # write(output_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_lower.tif", access_map_lower_raster, force = true)

        println("Raster construction complete.")
    end
end
