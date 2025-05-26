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
raster_base = replace_missing(Raster(OUTPUT_RASTERS_DIR*"inla_logmodel_npc/NPC_logmodel_$(2000)_sample_1.tif"), missingval = -NaN)

# Input and output directory for rasters
inla_dir = OUTPUT_RASTERS_DIR
raster_dir = OUTPUT_RASTERS_DIR
input_dir = OUTPUT_RASTERS_DIR
output_dir = OUTPUT_RASTERS_DIR

# Make paths tosave location if doesn't already exist
mkpath(output_dir*"final_netage/snf_netage/")
mkpath(output_dir*"final_npc/snf_npc/")
mkpath(output_dir*"final_npc/posterior_samples/")
mkpath(output_dir*"final_access/snf_access/")
mkpath(output_dir*"final_access/posterior_samples/")
mkpath(output_dir*"final_use/posterior_samples/")
mkpath(output_dir*"final_utilisation/posterior_samples/")

# %% Perform draws and save outputs. Filter out unwanted countries
ISO_list = ISO_LIST
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# # %% Time bounds
year = parse(Int, ARGS[1])
month = parse(Int, ARGS[2])

# YEAR_START = 2023
# YEAR_END = 2023

#####################################################
# %% Construct mean net age rasters from SNF
#####################################################
# Import net age data
mean_netage_data = CSV.read(OUTPUT_DATAPREP_DIR*"snf_mean_netage.csv", DataFrame)

# for year in ProgressBar(YEAR_START:YEAR_END)
#     for month in 1:12




monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEAR_NAT_START)

println("Constructing stock and flow rasters year $(year), month $(month)")

# %% Declare storage variables
netage_nat_snf_mean_rasters = Vector{Any}(undef, length(filt_ISOs))

println("Constructing rasters using subnational SNF values...$(year), month $(month)")
for ISO_i in 1:length(filt_ISOs)

    ISO = filt_ISOs[ISO_i]
    println("Selected country $(ISO)")

    snf_post_filename = ISO*"_SUBNAT_draws.jld2"
    snf_post_draws = load(snf_post_dir*snf_post_filename)
    println("Imported draws for $(ISO)")

    n_snf_posteriors = size(snf_post_draws["merged_outputs"][1]["ADJ_NPC_MONTHLY_TOTAL_samples"])[1] # Infer number of posteriors to sample from based on snf draws
    n_admin1 = length(snf_post_draws["merged_outputs"])
    println("Extracted loop bounds.")

    # Declare storage variables
    netage_subnat_snf_mean_rasters = []

    for subnat_i in 1:n_admin1
        println("Admin region $(subnat_i) of $(n_admin1)")
        admin1_id = snf_post_draws["merged_outputs"][subnat_i]["area_id"]

        netage_subnat_mean_estimate = mean_netage_data[mean_netage_data.area_id .== admin1_id .&&
                                        mean_netage_data.month .== month .&&
                                        mean_netage_data.year .== year,"mean_age_months_mean"][1]

        
        # Get required geometry
        admin1_geometry = admin1_shapes_geoIO[admin1_shapes_geoIO.area_id .== admin1_id,:].geometry

        # Mask raster base and trim to desired subregion and set values to SNF estimates of NPC and access
        subnat_masked = Rasters.trim(mask(raster_base, with = admin1_geometry); pad=0)
        
        # check if able to find region in default tiling. If region is too small, then skip in analysis
        if size(raster_base) == size(subnat_masked)
            continue
        end

        netage_subnat_mean_raster = copy(subnat_masked)
        netage_subnat_mean_raster[findall(.!isnan.(netage_subnat_mean_raster))] .= netage_subnat_mean_estimate

        push!(netage_subnat_snf_mean_rasters, netage_subnat_mean_raster)
    end

    netage_nat_snf_mean_rasters[ISO_i] = copy(netage_subnat_snf_mean_rasters)
end

# Combine Subnational rasters into national level rasters
Threads.@threads for ISO_i in ProgressBar(1:length(filt_ISOs), leave = false)
    netage_nat_snf_mean_rasters[ISO_i] = mosaic(first, netage_nat_snf_mean_rasters[ISO_i]..., atol = 0.01)
end

# Combine national level rasters into Africa raster
# Month string
month_str = "$(month)"
if month < 10
    month_str = "0$(month)"
end

# Mosaic and save to disk.
println("Combining all national rasters into Africa level raster and write to disk...")

println("Constructing Net Age SNF mean raster...")
combined_netage_snf_mean_raster = resample(mosaic(first, netage_nat_snf_mean_rasters..., atol = 0.01), to = raster_base)
write(output_dir*"final_netage/snf_netage/netage_$(year)_$(month_str)_mean.tif", combined_netage_snf_mean_raster, force = true)
println("Constructed Net Age SNF mean raster...")







#     end
# end




#####################################################
# %% Loop to first construct SNF block map upto subnational resolution
#####################################################


# for year in ProgressBar(YEAR_START:YEAR_END)
#     for month in 1:12







monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEAR_NAT_START)

println("Constructing stock and flow rasters year $(year), month $(month)")

# %% Import subnat draws for first country just to get n_sample posteriors
n_snf_posteriors = size((load(snf_post_dir*(filt_ISOs[1]*"_SUBNAT_draws.jld2")))["merged_outputs"][1]["ADJ_NPC_MONTHLY_TOTAL_samples"])[1] # Infer number of posteriors to sample from based on snf draws

for sample_i in 1:n_samples
    # %% Generate index to take SNF sample value from
    sample_idx = rand(1:n_snf_posteriors)

    # %% Declare storage variables
    npc_nat_snf_sample_rasters = Vector{Any}(undef, length(filt_ISOs))
    access_nat_snf_sample_rasters = Vector{Any}(undef, length(filt_ISOs))

    println("Constructing rasters using subnational SNF values...$(year), month $(month), sample $(sample_i)/$(n_samples)")
    for ISO_i in 1:length(filt_ISOs)
        ISO = filt_ISOs[ISO_i]
        println("Sample $(sample_i)/$(n_samples), Country $(ISO_i)/$(length(filt_ISOs)) - $(ISO)")

        snf_post_filename = ISO*"_SUBNAT_draws.jld2"
        snf_post_draws = load(snf_post_dir*snf_post_filename)
        
        n_snf_posteriors = size(snf_post_draws["merged_outputs"][1]["ADJ_NPC_MONTHLY_TOTAL_samples"])[1] # Infer number of posteriors to sample from based on snf draws
        n_admin1 = length(snf_post_draws["merged_outputs"])

        # Declare storage variables
        npc_subnat_snf_sample_rasters = []
        access_subnat_snf_sample_rasters = []

        for subnat_i in 1:n_admin1
            println("$(ISO) Admin1 Region $(subnat_i)/$(n_admin1)")
            admin1_id = snf_post_draws["merged_outputs"][subnat_i]["area_id"]

            # Get SNF estimates for NPC and Access
            npc_subnat_sample_estimate = snf_post_draws["merged_outputs"][subnat_i]["ADJ_NPC_MONTHLY_TOTAL_samples"][sample_idx,monthidx]

            ### Remove NaNs for this draw possibly due to transition period in problematic countries (MIGHT WANT TO CHECK)
            access_subnat_draws = snf_post_draws["merged_outputs"][subnat_i]["ADJ_λ_ACCESS_samples"][:,monthidx]
            access_subnat_draws = access_subnat_draws[findall(.!isnan.(access_subnat_draws))]
            access_subnat_sample_estimate = access_subnat_draws[sample_idx]
            
            if isnan(access_subnat_sample_estimate) # i.e. no data available, so no valid access. Assume 0 access
                access_subnat_sample_estimate = 0
            end
            
            # Get required geometry
            admin1_geometry = admin1_shapes_geoIO[admin1_shapes_geoIO.area_id .== admin1_id,:].geometry

            # Mask raster base and trim to desired subregion and set values to SNF estimates of NPC and access
            subnat_masked = Rasters.trim(mask(raster_base, with = admin1_geometry); pad=0)
            
            # check if able to find region in default tiling. If region is too small, then skip in analysis
            if size(raster_base) == size(subnat_masked)
                continue
            end

            npc_subnat_sample_raster = copy(subnat_masked)
            access_subnat_sample_raster = copy(subnat_masked)

            # Check if regions were valid to begin with
            if isempty(findall(.!isnan.(npc_subnat_sample_raster))) # For some reason if not accounted for, causes problems with mosaic()
                continue
            end

            npc_subnat_sample_raster[findall(.!isnan.(npc_subnat_sample_raster))] .= npc_subnat_sample_estimate
            access_subnat_sample_raster[findall(.!isnan.(access_subnat_sample_raster))] .= access_subnat_sample_estimate

            push!(npc_subnat_snf_sample_rasters, npc_subnat_sample_raster)
            push!(access_subnat_snf_sample_rasters, access_subnat_sample_raster)
        end

        npc_nat_snf_sample_rasters[ISO_i] = copy(npc_subnat_snf_sample_rasters)
        access_nat_snf_sample_rasters[ISO_i] = copy(access_subnat_snf_sample_rasters)
    end

    # Combine Subnational rasters into national level rasters
    println("Mosaicing subnational rasters")
    Threads.@threads for ISO_i in ProgressBar(1:length(filt_ISOs), leave = false)
        npc_nat_snf_sample_rasters[ISO_i] = mosaic(first, npc_nat_snf_sample_rasters[ISO_i]..., atol = 0.01)
        access_nat_snf_sample_rasters[ISO_i] = mosaic(first, access_nat_snf_sample_rasters[ISO_i]..., atol = 0.01)
    end
    
    # Combine national level rasters into Africa raster
    # Month string
    month_str = "$(month)"
    if month < 10
        month_str = "0$(month)"
    end

    # Mosaic and save to disk.
    println("Combining all national rasters into Africa level raster and write to disk...")
    for thread_i in 1:2
        if thread_i == 1
            println("Constructing NPC SNF sampled raster...")
            combined_npc_snf_sample_raster = resample(mosaic(first, npc_nat_snf_sample_rasters..., atol = 0.01), to = raster_base)
            write(output_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_sample_$(sample_i).tif", combined_npc_snf_sample_raster, force = true)
            println("Constructed NPC SNF sampled raster...")
        else thread_i == 2
            println("Constructing Access SNF sampled raster...")
            combined_access_snf_sample_raster = resample(mosaic(first, access_nat_snf_sample_rasters..., atol = 0.01), to = raster_base)
            write(output_dir*"final_access/snf_access/access_$(year)_$(month_str)_sample_$(sample_i).tif", combined_access_snf_sample_raster, force = true)
            println("Constructed Access SNF sampled raster...")
        end
    end

    println("Raster construction complete for $(year), month $(month), sample $(sample_i)/$(n_samples)")
end







#     end
# end





#####################################################
# %% Loop to construct unadjusted NPC, Access and Use rasters using INLA deviation regression outputs
#####################################################
# for year in YEAR_START:YEAR_END
#     for month in 1:12






        println("Constructing spatial disaggregated rasters year $(year), month $(month)")

        # Get month string for importing files
        month_str = "$(month)"
        if month < 10
            month_str = "0$(month)"
        end

        for sample_i in 1:n_samples
            println("Constructing rasters for year $(year), month $(month),  sample $(sample_i)/$(n_samples)")
            
            # Import stock and flow rasters
            println("Importing rasters...")

            npc_snf_sample_raster = replace_missing(Raster(input_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_sample_$(sample_i).tif"), missingval = NaN)
            access_snf_sample_raster = replace_missing(Raster(input_dir*"final_access/snf_access/access_$(year)_$(month_str)_sample_$(sample_i).tif"), missingval = NaN)

            #################################
            # Construct mean rasters
            #################################
            # Import log model npc rasters
            logmodel_npc_sample_raster = replace_missing(Raster(inla_dir*"inla_logmodel_npc/NPC_logmodel_$(year)_sample_$(sample_i).tif"), missingval = NaN)

            # Import pmodel access rasters
            pmodel_access_sample_raster =  replace_missing(Raster(inla_dir*"inla_pmodel_access/ACCESS_pmodel_$(year)_sample_$(sample_i).tif"), missingval = NaN)

            # Import logismodel use rasters
            logis_use_sample_raster = replace_missing(Raster(inla_dir*"inla_use_logis/USE_logismodel_$(year)_$(month)_sample_$(sample_i).tif"), missingval = NaN)

            # Consruct rasters
            println("Calculating sample NPC raster...")
            logmodel_npc_ratio_sample_raster = exp.(logmodel_npc_sample_raster)
            npc_map_sample_raster = logmodel_npc_ratio_sample_raster.* npc_snf_sample_raster
            println("Calculating sample Access raster...")
            access_map_sample_raster = inv_p_transform.(pmodel_access_sample_raster, access_snf_sample_raster, n=2)
            println("Calculating sample Use raster...")
            use_map_sample_raster = inv_p_transform.(logis_use_sample_raster, access_map_sample_raster, n=2)
            println("Calculating sample Utilisation raster...")
            util_map_sample_raster = use_map_sample_raster./access_map_sample_raster

            # Save mean rasters, overwrites previous rasters
            println("Saving sample rasters...")
            rm(output_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_sample_$(sample_i).tif"; force=true)
            rm(output_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_sample_$(sample_i).tif"; force=true)
            rm(output_dir*"final_use/logis_use/use_$(year)_$(month_str)_sample_$(sample_i).tif"; force=true)
            rm(output_dir*"final_utilisation/monthly/utilisation_$(year)_$(month_str)_sample_$(sample_i).tif"; force=true)
            
            write(output_dir*"final_npc/posterior_samples/npc_$(year)_$(month_str)_sample_$(sample_i).tif", npc_map_sample_raster, force = true)
            write(output_dir*"final_access/posterior_samples/access_$(year)_$(month_str)_sample_$(sample_i).tif", access_map_sample_raster, force = true)
            write(output_dir*"final_use/posterior_samples/use_$(year)_$(month_str)_sample_$(sample_i).tif", use_map_sample_raster, force = true)
            write(output_dir*"final_utilisation/posterior_samples/utilisation_$(year)_$(month_str)_sample_$(sample_i).tif", util_map_sample_raster, force = true)

            println("Raster construction complete.")
        end





        
#     end
# end

