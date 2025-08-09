"""
Author: Eugene Tan
Date Created: 14/5/2025
Last Updated: 21/5/2025
Extracts time series from rasters and combines with Stock and flow estimates to yield a large dataframe summary of the time series.
Saves as a CSV. Done on a country level.
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
using Shapefile
using GeoInterface
using GeoIO
using JLD2
using StatsBase

# Custom packages
using DateConversions
using UsefulTransformations
using RasterLookup

#######################################
# %% Define File paths and parameters from TOML file.
#######################################
# MITN Posterior Estimates
snf_post_dir = OUTPUT_DRAWS_DIR*"subnational/"

# Region Admin 1 area id legend
admin1_legend_dir = RAW_DATASET_DIR*"subnational/"
admin1_legend_filename = ADMIN1_AREAID_LEGEND_FILENAME

# %% Location to save data to
output_dir = OUTPUT_DIR*"coverage_timeseries/master_extractions_parts/admin0/"

# Population Rasters directory
pop_dir = POPULATION_RASTER_DIR

# Input and output directory for ITN rasters
raster_dir = OUTPUT_RASTERS_DIR

# INLA Spatial Disaggregation survey data
survey_data_dir = OUTPUT_DATAPREP_DIR
survey_data_filename = INLA_REDUCED_DATAPREP_FILENAME

# Num of samples to import from INLA outputs for use
n_samples = INLA_N_SAMPLES # Number of samples that were saved in the INLA raster process

# %% Get list of countries to analyse
ISO_list = ISO_LIST
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

#######################################
# %% Load Datasets
#######################################

# Region boundaries
admin0_shapes_geoIO = GeoIO.load(ADMIN0_SHAPEFILE)
admin1_shapes_geoIO = GeoIO.load(ADMIN1_SHAPEFILE)


#######################################
# %% Take input of year from arg list
#######################################

# %% Run extraction based on input arguments
year = parse(Int, ARGS[1])
month = parse(Int, ARGS[2])

# # Loop extraction code for each year and month ############### UNCOMMENT THIS IF MANUALLY RUNNING SCRIPT
# for year in YEAR_START:YEAR_END

# Import population raster
pop_year = year
if pop_year <= 2020
    # pop_year = min(max(year, 2000), 2020)
    population_raster = replace_missing(Raster(pop_dir*"WorldPop_UNAdj_v3_DRC_fix.$(pop_year).Annual.Data.5km.sum.tif"), missingval = NaN)
else
    population_raster = replace_missing(Raster("/mnt/s3/mastergrids/Other_Global_Covariates/Population/WorldPop/v3/PopulationCounts_DRC_fixed_projected/5km/"*"WorldPop_UNAdj_v3_DRC_fix_projected.$(pop_year).Annual.Data.5km.sum.tif"), missingval = NaN)
end

# for month in 1:12  ############### UNCOMMENT THIS IF MANUALLY RUNNING SCRIPT

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

# Import ITN sample rasters
npc_sample_raster_filenames = "npc_$(year)_$(month_str)_sample_".*string.(1:n_samples).*".tif"
# npc_sample_rasters = replace_missing.(Raster.(raster_dir*"final_npc/posterior_samples/".*npc_sample_raster_filenames), missingval = NaN)

access_sample_raster_filenames = "access_$(year)_$(month_str)_sample_".*string.(1:n_samples).*".tif"
# access_sample_rasters = replace_missing.(Raster.(raster_dir*"final_access/posterior_samples/".*access_sample_raster_filenames), missingval = NaN)

use_sample_raster_filenames = "use_$(year)_$(month_str)_sample_".*string.(1:n_samples).*".tif"
# use_sample_rasters = replace_missing.(Raster.(raster_dir*"final_use/posterior_samples/".*use_sample_raster_filenames), missingval = NaN)

util_sample_raster_filenames = "utilisation_$(year)_$(month_str)_sample_".*string.(1:n_samples).*".tif"
# util_sample_rasters = replace_missing.(Raster.(raster_dir*"final_utilisation/posterior_samples/monthly/".*util_sample_raster_filenames), missingval = NaN)

# Pre define storage variables
# Admin0 data
admin0_npc_raster_vals = Matrix{Float32}(undef, n_samples, length(filt_ISOs))
admin0_access_raster_vals = Matrix{Float32}(undef, n_samples, length(filt_ISOs))
admin0_use_raster_vals = Matrix{Float32}(undef, n_samples, length(filt_ISOs))
admin0_util_raster_vals = Matrix{Float32}(undef, n_samples, length(filt_ISOs))
admin0_eff_raster_vals = Matrix{Float32}(undef, n_samples, length(filt_ISOs))
admin0_snf_npc_ci = Matrix{Float32}(undef, length(filt_ISOs),3)
admin0_snf_access_ci = Matrix{Float32}(undef, length(filt_ISOs),3)

admin0_pop_vals = Vector{Float32}(undef, length(filt_ISOs))
admin0_name = Vector{String}(undef, length(filt_ISOs))
admin0_areaid = Vector{Int64}(undef, length(filt_ISOs))
admin0_geometry = Vector{Any}(undef, length(filt_ISOs))

admin1_npc_collection = Vector{Matrix{Float32}}(undef, length(filt_ISOs))
admin1_access_collection = Vector{Matrix{Float32}}(undef, length(filt_ISOs))
admin1_use_collection = Vector{Matrix{Float32}}(undef, length(filt_ISOs))
admin1_util_collection = Vector{Matrix{Float32}}(undef, length(filt_ISOs))
admin1_eff_collection = Vector{Matrix{Float32}}(undef, length(filt_ISOs))
admin1_snf_npc_collection = Vector{Matrix{Float32}}(undef, length(filt_ISOs))
admin1_snf_access_collection = Vector{Matrix{Float32}}(undef, length(filt_ISOs))

admin1_pop_collection = Vector{Vector{Float32}}(undef, length(filt_ISOs))
admin1_name_collection = Vector{Vector{String}}(undef, length(filt_ISOs))
admin1_areaid_collection = Vector{Vector{Int64}}(undef, length(filt_ISOs))
admin1_geometry_collection = Vector{Vector{Any}}(undef, length(filt_ISOs))

###########################################################
# Scrape SNF and Metadata from SNF Draws
###########################################################

Threads.@threads for ISO_i in 1:length(filt_ISOs)
    ISO = filt_ISOs[ISO_i]
    println("SNF Extraction Country [$(ISO_i)/$(length(filt_ISOs))]: $(ISO)...")

    # Load Subnational Stock and Flow Data to get metadata
    subnat_snf_filename = "$(ISO)_SUBNAT_draws.jld2"
    subnat_snf_post_draws = JLD2.load(snf_post_dir*subnat_snf_filename)

    # Get n_admin1 for loop bounds
    n_admin1 = length(subnat_snf_post_draws["admin1_names"])

    # Define required storage variables
    admin1_snf_npc_ci = Matrix{Float32}(undef, n_admin1, 3)
    admin1_snf_access_ci = Matrix{Float32}(undef, n_admin1, 3)
    admin1_pop = Vector{Float32}(undef, n_admin1)
    admin1_names = Vector{String}(undef, n_admin1)
    admin1_areaid = Vector{Int32}(undef, n_admin1)
    admin1_geometry = Vector{Any}(undef, n_admin1)
    # admin1_popmask = Vector{Raster}(undef, n_admin1)

    # Scrape Admin1 data from SNF Draws
    for admin1_i in 1:n_admin1
        println("Extracting SNF [$(ISO_i)/$(length(filt_ISOs))]: $(ISO), region $(admin1_i) of $(n_admin1)")
        admin1_names[admin1_i] = subnat_snf_post_draws["admin1_names"][admin1_i]
        admin1_areaid[admin1_i] = subnat_snf_post_draws["merged_outputs"][admin1_i]["area_id"]
        admin1_geometry[admin1_i] = admin1_shapes_geoIO[findfirst(admin1_shapes_geoIO.area_id .== admin1_areaid[admin1_i]),"geometry"]
        pop_mask = Rasters.trim(mask(population_raster, with = admin1_geometry[admin1_i]); pad=0)
        admin1_pop[admin1_i] = sum(pop_mask[findall(.!isnan.(pop_mask))])

        # Stock and flow NPC values and CI
        admin1_snf_npc_ci[admin1_i,[1,3]] = quantile(subnat_snf_post_draws["merged_outputs"][admin1_i]["ADJ_NPC_MONTHLY_TOTAL_samples"][:,monthidx], [0.025, 0.975])
        admin1_snf_npc_ci[admin1_i,2] = mean(subnat_snf_post_draws["merged_outputs"][admin1_i]["ADJ_NPC_MONTHLY_TOTAL_samples"][:,monthidx])

        # Stock and flow Access values and CI
        admin1_snf_access_ci[admin1_i,[1,3]] = quantile(subnat_snf_post_draws["merged_outputs"][admin1_i]["ADJ_λ_ACCESS_samples"][:,monthidx], [0.025, 0.975])
        admin1_snf_access_ci[admin1_i,2] = mean(subnat_snf_post_draws["merged_outputs"][admin1_i]["ADJ_λ_ACCESS_samples"][:,monthidx])
    end

    # Store Admin1 data to storage variable
    admin1_name_collection[ISO_i] = admin1_names
    admin1_areaid_collection[ISO_i] = admin1_areaid
    admin1_geometry_collection[ISO_i] = admin1_geometry
    # admin1_popmask_collection[ISO_i] = admin1_popmask
    admin1_pop_collection[ISO_i] = admin1_pop
    admin1_snf_npc_collection[ISO_i] = admin1_snf_npc_ci
    admin1_snf_access_collection[ISO_i] = admin1_snf_access_ci

    # Scrape Admin0 data from SNF Draws
    admin0_name[ISO_i] = admin0_shapes_geoIO[findfirst(admin0_shapes_geoIO.ISO .== ISO),"Name_0"]
    admin0_areaid[ISO_i] = admin0_shapes_geoIO[findfirst(admin0_shapes_geoIO.ISO .== ISO),"area_id"]
    admin0_geometry[ISO_i] = admin0_shapes_geoIO[findfirst(admin0_shapes_geoIO.ISO .== ISO),"geometry"]
    pop_mask = Rasters.trim(mask(population_raster, with = admin0_geometry[ISO_i]); pad=0)
    admin0_pop_vals[ISO_i] = sum(pop_mask[findall(.!isnan.(pop_mask))])

    admin0_snf_npc_ci[ISO_i,1] = sum(admin1_pop.*admin1_snf_npc_ci[:,1])/sum(admin1_pop)
    admin0_snf_npc_ci[ISO_i,2] = sum(admin1_pop.*admin1_snf_npc_ci[:,2])/sum(admin1_pop)
    admin0_snf_npc_ci[ISO_i,3] = sum(admin1_pop.*admin1_snf_npc_ci[:,3])/sum(admin1_pop)

    admin0_snf_access_ci[ISO_i,1] = sum(admin1_pop.*admin1_snf_access_ci[:,1])/sum(admin1_pop)
    admin0_snf_access_ci[ISO_i,2] = sum(admin1_pop.*admin1_snf_access_ci[:,2])/sum(admin1_pop)
    admin0_snf_access_ci[ISO_i,3] = sum(admin1_pop.*admin1_snf_access_ci[:,3])/sum(admin1_pop)

    # Assign storage values for raster sample aggregates at admin1
    npc_raster_vals = Matrix{Float32}(undef, n_samples, n_admin1)
    access_raster_vals = Matrix{Float32}(undef, n_samples, n_admin1)
    use_raster_vals = Matrix{Float32}(undef, n_samples, n_admin1)
    util_raster_vals = Matrix{Float32}(undef, n_samples, n_admin1)
    eff_raster_vals = Matrix{Float32}(undef, n_samples, n_admin1)

    admin1_npc_collection[ISO_i] = npc_raster_vals
    admin1_access_collection[ISO_i] = access_raster_vals
    admin1_use_collection[ISO_i] = use_raster_vals
    admin1_util_collection[ISO_i] = util_raster_vals
    admin1_util_collection[ISO_i] = eff_raster_vals
end

###########################################################
# Aggregate ITN Metrics from sampled INLA rasters
###########################################################
for sample_i in 1:n_samples
    # Import raster
    npc_sample_raster = replace_missing(Raster(raster_dir*"final_npc/posterior_samples/"*npc_sample_raster_filenames[sample_i]), missingval = NaN)
    access_sample_raster = replace_missing(Raster(raster_dir*"final_access/posterior_samples/"*access_sample_raster_filenames[sample_i]), missingval = NaN)
    use_sample_raster = replace_missing(Raster(raster_dir*"final_use/posterior_samples/"*use_sample_raster_filenames[sample_i]), missingval = NaN)
    # util_sample_raster = replace_missing(Raster(raster_dir*"final_utilisation/posterior_samples/"*util_sample_raster_filenames[sample_i]), missingval = NaN)

    for ISO_i in 1:length(filt_ISOs)

        ISO = filt_ISOs[ISO_i]

        # Admin 1 Raster Extraction
        n_admin1 = length(admin1_name_collection[ISO_i])

        # for admin1_i in 1:n_admin1
        #     println("Sample: [$(sample_i)/$(n_samples)], Country: [$(ISO_i)/$(length(filt_ISOs))] - ($(ISO)), Admin 1 Region [$(admin1_i)/$(n_admin1)]")
        #     geometry = admin1_geometry_collection[ISO_i][admin1_i]

        #     pop_masked = Rasters.trim(mask(population_raster, with = geometry); pad=0)

        #     npc_sample_raster_masked = resample(Rasters.trim(mask(npc_sample_raster, with = geometry); pad = 0), to = pop_masked)
        #     access_sample_raster_masked = resample(Rasters.trim(mask(access_sample_raster, with = geometry); pad = 0), to = pop_masked)
        #     use_sample_raster_masked = resample(Rasters.trim(mask(use_sample_raster, with = geometry); pad = 0), to = pop_masked)
        #     util_sample_raster_masked = resample(Rasters.trim(mask(util_sample_raster, with = geometry); pad = 0), to = pop_masked)

        #     admin1_npc_collection[ISO_i][sample_i,admin1_i] = aggregate_raster_weighted_mean(npc_sample_raster_masked, pop_masked)
        #     admin1_access_collection[ISO_i][sample_i,admin1_i] = aggregate_raster_weighted_mean(access_sample_raster_masked, pop_masked)
        #     admin1_use_collection[ISO_i][sample_i,admin1_i] = aggregate_raster_weighted_mean(use_sample_raster_masked, pop_masked)
        #     admin1_util_collection[ISO_i][sample_i,admin1_i] = aggregate_raster_weighted_mean(util_sample_raster_masked, pop_masked)
        # end

        # Admin 0 Raster Extraction
        println("Sample: [$(sample_i)/$(n_samples)], Country: [$(ISO_i)/$(length(filt_ISOs))] - ($(ISO)), Admin 0 Extraction")
        geometry = admin0_geometry[ISO_i]
        pop_masked = Rasters.trim(mask(population_raster, with = geometry); pad=0)

        npc_sample_raster_masked = resample(Rasters.trim(mask(npc_sample_raster, with = geometry); pad = 0), to = pop_masked)
        access_sample_raster_masked = resample(Rasters.trim(mask(access_sample_raster, with = geometry); pad = 0), to = pop_masked)
        use_sample_raster_masked = resample(Rasters.trim(mask(use_sample_raster, with = geometry); pad = 0), to = pop_masked)
        # util_sample_raster_masked = resample(Rasters.trim(mask(util_sample_raster, with = geometry); pad = 0), to = pop_masked)

        admin0_npc_raster_vals[sample_i,ISO_i] = aggregate_raster_weighted_mean(npc_sample_raster_masked, pop_masked)
        admin0_access_raster_vals[sample_i,ISO_i] = aggregate_raster_weighted_mean(access_sample_raster_masked, pop_masked)
        admin0_use_raster_vals[sample_i,ISO_i] = aggregate_raster_weighted_mean(use_sample_raster_masked, pop_masked)
        admin0_util_raster_vals[sample_i,ISO_i] = (admin0_use_raster_vals[sample_i,ISO_i]/2)/(admin0_npc_raster_vals[sample_i,ISO_i]) #aggregate_raster_weighted_mean(util_sample_raster_masked, pop_masked)
        admin0_eff_raster_vals[sample_i,ISO_i] = (admin0_use_raster_vals[sample_i,ISO_i])/(admin0_access_raster_vals[sample_i,ISO_i])
    end
end

###########################################################
# Compute Confidence Intervals
###########################################################

for ISO_i in 1:length(filt_ISOs)
    ISO = filt_ISOs[ISO_i]

    # Admin 1 Raster Extraction
    n_admin1 = length(admin1_name_collection[ISO_i])

    # # Dataframe storage variable for data
    # df_admin1_collection = []

    # for admin1_i in 1:n_admin1

    #     println("$(ISO), $admin1_i")
    #     admin1_name = admin1_name_collection[ISO_i][admin1_i]
    #     admin1_id = admin1_areaid_collection[ISO_i][admin1_i]
    #     pop_total = Float64(admin1_pop_collection[ISO_i][admin1_i])

    #     raster_npc_mean = Float64(mean(admin1_npc_collection[ISO_i][:, admin1_i]))
    #     raster_access_mean = Float64(mean(admin1_access_collection[ISO_i][:, admin1_i]))
    #     raster_use_mean = Float64(mean(admin1_use_collection[ISO_i][:, admin1_i]))
    #     raster_util_mean = Float64(mean(admin1_npc_collection[ISO_i][:, admin1_i]))

    #     npc_quantiles = Float64.(quantile(admin1_npc_collection[ISO_i][:, admin1_i], [0.025, 0.975]))
    #     access_quantiles = Float64.(quantile(admin1_access_collection[ISO_i][:, admin1_i], [0.025, 0.975]))
    #     use_quantiles = Float64.(quantile(admin1_use_collection[ISO_i][:, admin1_i], [0.025, 0.975]))
    #     util_quantiles = Float64.(quantile(admin1_util_collection[ISO_i][:, admin1_i], [0.025, 0.975]))

    #     snf_npc_lower, snf_npc_mean, snf_npc_upper = Float64.(admin1_snf_npc_collection[ISO_i][admin1_i,:])
    #     snf_access_lower, snf_access_mean, snf_access_upper = Float64.(admin1_snf_npc_collection[ISO_i][admin1_i,:])

    #     # Compile into Dataframe entry
    #     df1_i = DataFrame(ISO = ISO, category = "Admin1", 
    #                 admin_name = admin1_name, 
    #                 area_id = admin1_id,
    #                 month = month, year = year,
    #                 population = pop_total,
    #                 snf_npc_95lower = snf_npc_lower,
    #                 snf_npc_mean = snf_npc_mean,
    #                 snf_npc_95upper = snf_npc_upper,
    #                 raster_npc_95lower = npc_quantiles[1],
    #                 raster_npc_mean = raster_npc_mean,
    #                 raster_npc_95upper = npc_quantiles[2],
    #                 snf_access_95lower = snf_access_lower,
    #                 snf_access_mean = snf_access_mean,
    #                 snf_access_95upper = snf_access_upper,
    #                 raster_access_95lower = access_quantiles[1],
    #                 raster_access_mean = raster_access_mean,
    #                 raster_access_95upper = access_quantiles[2],
    #                 raster_use_95lower = use_quantiles[1],
    #                 raster_use_mean = raster_use_mean,
    #                 raster_use_95upper = use_quantiles[2],
    #                 raster_util_95lower = util_quantiles[1],
    #                 raster_util_mean = raster_util_mean,
    #                 raster_util_95upper = util_quantiles[2]
    #                 )

    #     push!(df_admin1_collection, df1_i)
    # end

    # # Concatenate admin1 entries for current country
    # df1_i_entries = vcat(df_admin1_collection...)
    # push!(df1_full_collection, df1_i_entries)

    # Compile DataFrame entry for admin0
    pop_total = Float64(admin0_pop_vals[ISO_i])

    raster_npc_mean = Float64(mean(admin0_npc_raster_vals[:, ISO_i]))
    raster_access_mean = Float64(mean(admin0_access_raster_vals[:, ISO_i]))
    raster_use_mean = Float64(mean(admin0_use_raster_vals[:, ISO_i]))
    npc_quantiles = Float64.(quantile(admin0_npc_raster_vals[:, ISO_i], [0.025, 0.975]))
    access_quantiles = Float64.(quantile(admin0_access_raster_vals[:, ISO_i], [0.025, 0.975]))
    use_quantiles = Float64.(quantile(admin0_use_raster_vals[:, ISO_i], [0.025, 0.975]))
    
    # Default assume no util or eff available
    util_quantiles = [NaN, NaN]
    eff_quantiles = [NaN, NaN]
    raster_util_mean = NaN
    raster_eff_mean = NaN
    
    nonnanidx_util = findall(.!isnan.(admin0_util_raster_vals[:, ISO_i]))
    nonnanidx_eff = findall(.!isnan.(admin0_eff_raster_vals[:, ISO_i]))
    
    if length(nonnanidx_util) > 2
      raster_util_mean = Float64(mean(admin0_util_raster_vals[nonnanidx_util, ISO_i]))
      util_quantiles = Float64.(quantile(admin0_util_raster_vals[nonnanidx_util, ISO_i], [0.025, 0.975]))
    end
    
    if length(nonnanidx_eff) > 2
      raster_eff_mean = Float64(mean(admin0_eff_raster_vals[nonnanidx_eff, ISO_i]))
      eff_quantiles = Float64.(quantile(admin0_eff_raster_vals[nonnanidx_eff, ISO_i], [0.025, 0.975]))
    end

    snf_npc_lower, snf_npc_mean, snf_npc_upper = Float64.(admin0_snf_npc_ci[ISO_i,:])
    snf_access_lower, snf_access_mean, snf_access_upper = Float64.(admin0_snf_access_ci[ISO_i,:])

    df0_i = DataFrame(ISO = ISO, category = "Admin0", 
                    admin_name = admin0_name[ISO_i], 
                    area_id = admin0_areaid[ISO_i],
                    month = month, year = year,
                    population = pop_total,
                    snf_npc_95lower = snf_npc_lower,
                    snf_npc_mean = snf_npc_mean,
                    snf_npc_95upper = snf_npc_upper,
                    raster_npc_95lower = npc_quantiles[1],
                    raster_npc_mean = raster_npc_mean,
                    raster_npc_95upper = npc_quantiles[2],
                    snf_access_95lower = snf_access_lower,
                    snf_access_mean = snf_access_mean,
                    snf_access_95upper = snf_access_upper,
                    raster_access_95lower = access_quantiles[1],
                    raster_access_mean = raster_access_mean,
                    raster_access_95upper = access_quantiles[2],
                    raster_use_95lower = use_quantiles[1],
                    raster_use_mean = raster_use_mean,
                    raster_use_95upper = use_quantiles[2],
                    raster_util_95lower = util_quantiles[1],
                    raster_util_mean = raster_util_mean,
                    raster_util_95upper = util_quantiles[2],
                    raster_eff_95lower = eff_quantiles[1],
                    raster_eff_mean = raster_eff_mean,
                    raster_eff_95upper = eff_quantiles[2]
                    )

    # Add data frame to admin0 df collection
    push!(df0_full_collection, df0_i)
end

# %% Combine collections together and save
compiled_dataframe = vcat(df0_full_collection...,df1_full_collection...)
mkpath(output_dir)
rm(output_dir*"mitn_extraction_$(year)_$(month).csv"; force=true)
CSV.write(output_dir*"mitn_extraction_$(year)_$(month).csv", compiled_dataframe)


#     end
# end

