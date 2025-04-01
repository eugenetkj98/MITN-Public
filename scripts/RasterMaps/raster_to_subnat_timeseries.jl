"""
Author: Eugene Tan
Date Created: 16/1/2025
Last Updated: 18/3/2025
Uses generated rasters to construct timeseries population weighted average estimates for ITN coverage at subnational resolution
Also extracts the model estimates from the previous BV model rasters and compiles into a master .jld2 dict file
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
using CSV

# Custom packages
using DateConversions

# %% File paths
# Population Rasters directory
pop_dir = "Z:/mastergrids/Other_Global_Covariates/Population/WorldPop/v3/PopulationCounts_DRC_fixed/5km/"

# Region boundaries
admin1_shapes_geoIO = GeoIO.load("Z:/master_geometries/Admin_Units/Global/MAP/2023/MG_5K/admin2023_1_MG_5K.shp")

# Directory of rasters
raster_dir = "Z:/eugene/Final Rasters/rasters/"#"outputs/rasters/"

# Output save dir
mkpath("outputs/coverage_timeseries/subnational/")
output_dir = "outputs/coverage_timeseries/subnational/"

# %% Perform draws and save outputs. Filter out unwanted countries
ISO_list = String.(CSV.read("datasets/ISO_list.csv", DataFrame)[:,1])
exclusion_ISOs = ["CPV","ZAF"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Time bounds
YEAR_START = 2000
YEAR_END = 2023

# %% Construct Timeseries for model predictions

# Year bounds
YEAR_LIST = YEAR_START:YEAR_END
n_years = length(YEAR_LIST)

# %% Define Storage variables for timeseries with confidence 95% intervals
mitn_subnat_npc = Array{Float64}(undef, length(filt_ISOs), n_years*12, 3)
mitn_subnat_access = Array{Float64}(undef, length(filt_ISOs), n_years*12, 3)
mitn_subnat_use = Array{Float64}(undef, length(filt_ISOs), n_years*12, 3)

# %% Define temporary DataFrame storage for compiling all the data
df_entries = []

# %% For each year, extract required population weighted estimates from rasters (Extraction for MITN)
for year in ProgressBar(YEAR_START:YEAR_END, leave = false)
	for month in ProgressBar(1:12, leave = false)
		# Calculate reference monthidx to access data
		monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEAR_START)

		# Get correct year/month string for importing file
		year_str = "$(year)"
		month_str = "$(month)"
		if month < 10
			month_str = "0"*month_str
		end

        # Import ITN coverage rasters (MITN)
        bv_npc_mean_raster = replace_missing(Raster("Z:/map-data-eng/airflow/itn_model/output/itn_cube/20241015/nets_per_capita/ITN_$(year)_percapita_nets_mean.tif"), missingval = NaN)
        
        bv_use_mean_raster = replace_missing(Raster("Z:/map-data-eng/airflow/itn_model/output/itn_cube/20241015/ITN_$(year)_use_mean.tif"), missingval = NaN)
        bv_use_upper_raster = replace_missing(Raster("Z:/map-data-eng/airflow/itn_model/output/itn_cube/20241015/ci/ITN_$(year)_use_upper.tif"), missingval = NaN)
        bv_use_lower_raster = replace_missing(Raster("Z:/map-data-eng/airflow/itn_model/output/itn_cube/20241015/ci/ITN_$(year)_use_lower.tif"), missingval = NaN)

		# Import ITN coverage rasters (MITN)
		mitn_npc_mean_raster = Raster(raster_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_mean.tif")
		mitn_npc_upper_raster = Raster(raster_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_upper.tif")
		mitn_npc_lower_raster = Raster(raster_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_lower.tif")

		mitn_access_mean_raster = Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_mean.tif")
		mitn_access_upper_raster = Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_upper.tif")
		mitn_access_lower_raster = Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_lower.tif")
		
		mitn_use_mean_raster = Raster(raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_mean.tif")
		mitn_use_upper_raster = Raster(raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_upper.tif")
		mitn_use_lower_raster = Raster(raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_lower.tif")

		# Import population raster
		pop_year = min(max(year, 2000), 2020)
		population_raster = replace_missing(Raster(pop_dir*"WorldPop_UNAdj_v3_DRC_fix.$(pop_year).Annual.Data.5km.sum.tif"), missingval = NaN)
		
        # Resample and align population raster to ITN Coverage rasters
		population_raster = resample(population_raster, to = bv_npc_mean_raster)

        df_entries_monthidx = []

        for ISO_i in ProgressBar(1:length(filt_ISOs), leave = false)
            # Get current country of interest
            ISO = filt_ISOs[ISO_i]

            # Retrieve metadata and geometries
            filt_admin1_metadata = admin1_shapes_geoIO[admin1_shapes_geoIO.ISO .== ISO,:]
            admin1_area_ids = unique(filt_admin1_metadata.area_id)
            
            n_admin1 = length(admin1_area_ids)

            admin1_names = []
            admin1_geometries = []
            for i in 1:n_admin1
                admin1_name = (filt_admin1_metadata[filt_admin1_metadata.area_id .== admin1_area_ids[i],:].Name_1)[1]
                admin1_geometry = (filt_admin1_metadata[filt_admin1_metadata.area_id .== admin1_area_ids[i],:].geometry)[1]

                push!(admin1_names, admin1_name)
                push!(admin1_geometries, admin1_geometry)
            end
            
            # Make Storage Variables for model estimates
            bv_npc_ests_mean = Vector{Float64}(undef, n_admin1)
            bv_npc_ests_upper = Vector{Float64}(undef, n_admin1)
            bv_npc_ests_lower = Vector{Float64}(undef, n_admin1)

            bv_use_ests_mean = Vector{Float64}(undef, n_admin1)
            bv_use_ests_upper = Vector{Float64}(undef, n_admin1)
            bv_use_ests_lower = Vector{Float64}(undef, n_admin1)

            mitn_npc_ests_mean = Vector{Float64}(undef, n_admin1)
            mitn_npc_ests_upper = Vector{Float64}(undef, n_admin1)
            mitn_npc_ests_lower = Vector{Float64}(undef, n_admin1)

            mitn_access_ests_mean = Vector{Float64}(undef, n_admin1)
            mitn_access_ests_upper = Vector{Float64}(undef, n_admin1)
            mitn_access_ests_lower = Vector{Float64}(undef, n_admin1)

            mitn_use_ests_mean = Vector{Float64}(undef, n_admin1)
            mitn_use_ests_upper = Vector{Float64}(undef, n_admin1)
            mitn_use_ests_lower = Vector{Float64}(undef, n_admin1)

            # Extract population weighted means of measures from rasters
            for i in ProgressBar(1:n_admin1, leave = false)
                # Get geometry
                admin1_geometry = admin1_geometries[i]

                # Get masked version of rasters
                pop_masked = Rasters.trim(mask(population_raster, with = admin1_geometry), pad = 0)

                bv_npc_mean_raster_masked = resample(Rasters.trim(mask(bv_npc_mean_raster, with = admin1_geometry), pad = 0), to = pop_masked)

                bv_use_mean_raster_masked = resample(Rasters.trim(mask(bv_use_mean_raster, with = admin1_geometry), pad = 0), to = pop_masked)
                bv_use_upper_raster_masked = resample(Rasters.trim(mask(bv_use_upper_raster, with = admin1_geometry), pad = 0), to = pop_masked)
                bv_use_lower_raster_masked = resample(Rasters.trim(mask(bv_use_lower_raster, with = admin1_geometry), pad = 0), to = pop_masked)

                mitn_npc_mean_raster_masked = resample(Rasters.trim(mask(mitn_npc_mean_raster, with = admin1_geometry), pad = 0), to = pop_masked)
                mitn_npc_upper_raster_masked = resample(Rasters.trim(mask(mitn_npc_upper_raster, with = admin1_geometry), pad = 0), to = pop_masked)
                mitn_npc_lower_raster_masked = resample(Rasters.trim(mask(mitn_npc_lower_raster, with = admin1_geometry), pad = 0), to = pop_masked)

                mitn_access_mean_raster_masked = resample(Rasters.trim(mask(mitn_access_mean_raster, with = admin1_geometry), pad = 0), to = pop_masked)
                mitn_access_upper_raster_masked = resample(Rasters.trim(mask(mitn_access_upper_raster, with = admin1_geometry), pad = 0), to = pop_masked)
                mitn_access_lower_raster_masked = resample(Rasters.trim(mask(mitn_access_lower_raster, with = admin1_geometry), pad = 0), to = pop_masked)

                mitn_use_mean_raster_masked = resample(Rasters.trim(mask(mitn_use_mean_raster, with = admin1_geometry), pad = 0), to = pop_masked)
                mitn_use_upper_raster_masked = resample(Rasters.trim(mask(mitn_use_upper_raster, with = admin1_geometry), pad = 0), to = pop_masked)
                mitn_use_lower_raster_masked = resample(Rasters.trim(mask(mitn_use_lower_raster, with = admin1_geometry), pad = 0), to = pop_masked)
                
                # Get index of pixels with non missing values
                nonmissing_idx_bv_npc_mean = findall(.!isnan.(pop_masked) .* .!isnan.(bv_npc_mean_raster_masked))

                nonmissing_idx_bv_use_mean = findall(.!isnan.(pop_masked) .* .!isnan.(bv_use_mean_raster_masked))
                nonmissing_idx_bv_use_upper = findall(.!isnan.(pop_masked) .* .!isnan.(bv_use_upper_raster_masked))
                nonmissing_idx_bv_use_lower = findall(.!isnan.(pop_masked) .* .!isnan.(bv_use_lower_raster_masked))

                nonmissing_idx_mitn_npc_mean = findall(.!isnan.(pop_masked) .* .!isnan.(mitn_npc_mean_raster_masked))
                nonmissing_idx_mitn_npc_upper = findall(.!isnan.(pop_masked) .* .!isnan.(mitn_npc_upper_raster_masked))
                nonmissing_idx_mitn_npc_lower = findall(.!isnan.(pop_masked) .* .!isnan.(mitn_npc_lower_raster_masked))

                nonmissing_idx_mitn_access_mean = findall(.!isnan.(pop_masked) .* .!isnan.(mitn_access_mean_raster_masked))
                nonmissing_idx_mitn_access_upper = findall(.!isnan.(pop_masked) .* .!isnan.(mitn_access_upper_raster_masked))
                nonmissing_idx_mitn_access_lower = findall(.!isnan.(pop_masked) .* .!isnan.(mitn_access_lower_raster_masked))

                nonmissing_idx_mitn_use_mean = findall(.!isnan.(pop_masked) .* .!isnan.(mitn_use_mean_raster_masked))
                nonmissing_idx_mitn_use_upper = findall(.!isnan.(pop_masked) .* .!isnan.(mitn_use_upper_raster_masked))
                nonmissing_idx_mitn_use_lower = findall(.!isnan.(pop_masked) .* .!isnan.(mitn_use_lower_raster_masked))

                # Get population weighted estimates
                bv_npc_ests_mean[i] = sum(pop_masked[nonmissing_idx_bv_npc_mean] .* bv_npc_mean_raster_masked[nonmissing_idx_bv_npc_mean])/sum(pop_masked[nonmissing_idx_bv_npc_mean])

                bv_use_ests_mean[i] = sum(pop_masked[nonmissing_idx_bv_use_mean] .* bv_use_mean_raster_masked[nonmissing_idx_bv_use_mean])/sum(pop_masked[nonmissing_idx_bv_use_mean])
                bv_use_ests_upper[i] = sum(pop_masked[nonmissing_idx_bv_use_upper] .* bv_use_upper_raster_masked[nonmissing_idx_bv_use_upper])/sum(pop_masked[nonmissing_idx_bv_use_upper])
                bv_use_ests_lower[i] = sum(pop_masked[nonmissing_idx_bv_use_lower] .* bv_use_lower_raster_masked[nonmissing_idx_bv_use_lower])/sum(pop_masked[nonmissing_idx_bv_use_lower])

                mitn_npc_ests_mean[i] = sum(pop_masked[nonmissing_idx_mitn_npc_mean] .* mitn_npc_mean_raster_masked[nonmissing_idx_mitn_npc_mean])/sum(pop_masked[nonmissing_idx_mitn_npc_mean])
                mitn_npc_ests_upper[i] = sum(pop_masked[nonmissing_idx_mitn_npc_upper] .* mitn_npc_upper_raster_masked[nonmissing_idx_mitn_npc_upper])/sum(pop_masked[nonmissing_idx_mitn_npc_upper])
                mitn_npc_ests_lower[i] = sum(pop_masked[nonmissing_idx_mitn_npc_lower] .* mitn_npc_lower_raster_masked[nonmissing_idx_mitn_npc_lower])/sum(pop_masked[nonmissing_idx_mitn_npc_lower])

                mitn_access_ests_mean[i] = sum(pop_masked[nonmissing_idx_mitn_access_mean] .* mitn_access_mean_raster_masked[nonmissing_idx_mitn_access_mean])/sum(pop_masked[nonmissing_idx_mitn_access_mean])
                mitn_access_ests_upper[i] = sum(pop_masked[nonmissing_idx_mitn_access_upper] .* mitn_access_upper_raster_masked[nonmissing_idx_mitn_access_upper])/sum(pop_masked[nonmissing_idx_mitn_access_upper])
                mitn_access_ests_lower[i] = sum(pop_masked[nonmissing_idx_mitn_access_lower] .* mitn_access_lower_raster_masked[nonmissing_idx_mitn_access_lower])/sum(pop_masked[nonmissing_idx_mitn_access_lower])

                mitn_use_ests_mean[i] = sum(pop_masked[nonmissing_idx_mitn_use_mean] .* mitn_use_mean_raster_masked[nonmissing_idx_mitn_use_mean])/sum(pop_masked[nonmissing_idx_mitn_use_mean])
                mitn_use_ests_upper[i] = sum(pop_masked[nonmissing_idx_mitn_use_upper] .* mitn_use_upper_raster_masked[nonmissing_idx_mitn_use_upper])/sum(pop_masked[nonmissing_idx_mitn_use_upper])
                mitn_use_ests_lower[i] = sum(pop_masked[nonmissing_idx_mitn_use_lower] .* mitn_use_lower_raster_masked[nonmissing_idx_mitn_use_lower])/sum(pop_masked[nonmissing_idx_mitn_use_lower])
            end

            df_entry = DataFrame(ISO = ISO,
                                    month = month,
                                    year = year,
                                    monthidx = monthidx,
                                    admin1_name = admin1_names,
                                    area_id = admin1_area_ids,
                                    bv_npc_lower = NaN,
                                    bv_npc_mean = bv_npc_ests_mean,
                                    bv_npc_upper = NaN,
                                    bv_access_lower = NaN,
                                    bv_access_mean = NaN,
                                    bv_access_upper = NaN,
                                    bv_use_lower = bv_use_ests_lower,
                                    bv_use_mean = bv_use_ests_mean,
                                    bv_use_upper = bv_use_ests_upper,
                                    mitn_npc_lower = mitn_npc_ests_lower,
                                    mitn_npc_mean = mitn_npc_ests_mean,
                                    mitn_npc_upper = mitn_npc_ests_upper,
                                    mitn_access_lower = mitn_access_ests_lower,
                                    mitn_access_mean = mitn_access_ests_mean,
                                    mitn_access_upper = mitn_access_ests_upper,
                                    mitn_use_lower = mitn_use_ests_lower,
                                    mitn_use_mean = mitn_use_ests_mean,
                                    mitn_use_upper = mitn_use_ests_upper)

            
            push!(df_entries_monthidx, df_entry)
        end

        # CSV.write("outputs/coverage_timeseries/model_prediction_comparisons_part_$(monthidx).csv", vcat(df_entries_monthidx...))
        push!(df_entries, vcat(df_entries_monthidx...))

	end
end

CSV.write(output_dir*"model_prediction_subnational_comparisons.csv", vcat(df_entries...))

# # %%
# df_temp = []
# for monthidx in ProgressBar(1:48)
#     push!(df_temp, CSV.read("outputs/coverage_timeseries/model_prediction_comparisons_part_$(monthidx).csv", DataFrame))
# end


# %% Reform DataFrame into time series dictionary to be saved as JLD2 file
full_dataset = CSV.read(output_dir*"model_prediction_subnational_comparisons.csv", DataFrame)
n_months = length(YEAR_LIST)*12
filt_ISOs = unique(full_dataset.ISO)

for ISO in filt_ISOs
    filt_data = full_dataset[full_dataset.ISO .== ISO,:]

    admin1_area_ids = unique(filt_data.area_id)
    n_admin1 = length(admin1_area_ids)

    # Make storage variables
    admin1_names = Array{String}(undef, n_admin1)

    bv_subnat_npc = Array{Float64}(undef, n_admin1, length(YEAR_LIST)*12, 3)
    bv_subnat_use = Array{Float64}(undef, n_admin1, length(YEAR_LIST)*12, 3)

    mitn_subnat_npc = Array{Float64}(undef, n_admin1, length(YEAR_LIST)*12, 3)
    mitn_subnat_access = Array{Float64}(undef, n_admin1, length(YEAR_LIST)*12, 3)
    mitn_subnat_use = Array{Float64}(undef, n_admin1, length(YEAR_LIST)*12, 3)
    filt_data
    # Extract data from DataFrame and put in storage variable
    for admin1_i in ProgressBar(1:n_admin1, leave = false)
        admin1_names[admin1_i] = (filt_data[filt_data.area_id .== admin1_area_ids[admin1_i],:].admin1_name)[1]
        area_id = admin1_area_ids[admin1_i]

        for monthidx in 1:n_months
            bv_subnat_npc[admin1_i,monthidx,1] = NaN
            bv_subnat_npc[admin1_i,monthidx,2] = (filt_data[(filt_data.area_id .== area_id) .& 
                                                            (filt_data.monthidx .== monthidx),:].bv_npc_mean)[1]
            bv_subnat_npc[admin1_i,monthidx,3] = NaN

            bv_subnat_use[admin1_i,monthidx,1] = (filt_data[(filt_data.area_id .== area_id) .& 
                                                            (filt_data.monthidx .== monthidx),:].bv_use_lower)[1]
            bv_subnat_use[admin1_i,monthidx,2] = (filt_data[(filt_data.area_id .== area_id) .& 
                                                            (filt_data.monthidx .== monthidx),:].bv_use_mean)[1]
            bv_subnat_use[admin1_i,monthidx,3] = (filt_data[(filt_data.area_id .== area_id) .& 
                                                            (filt_data.monthidx .== monthidx),:].bv_use_upper)[1]

            mitn_subnat_npc[admin1_i,monthidx,1] = (filt_data[(filt_data.area_id .== area_id) .& 
                                                            (filt_data.monthidx .== monthidx),:].mitn_npc_lower)[1]
            mitn_subnat_npc[admin1_i,monthidx,2] = (filt_data[(filt_data.area_id .== area_id) .& 
                                                            (filt_data.monthidx .== monthidx),:].mitn_npc_mean)[1]
            mitn_subnat_npc[admin1_i,monthidx,3] = (filt_data[(filt_data.area_id .== area_id) .& 
                                                            (filt_data.monthidx .== monthidx),:].mitn_npc_upper)[1]

            mitn_subnat_access[admin1_i,monthidx,1] = (filt_data[(filt_data.area_id .== area_id) .& 
                                                            (filt_data.monthidx .== monthidx),:].mitn_access_lower)[1]
            mitn_subnat_access[admin1_i,monthidx,2] = (filt_data[(filt_data.area_id .== area_id) .& 
                                                            (filt_data.monthidx .== monthidx),:].mitn_access_mean)[1]
            mitn_subnat_access[admin1_i,monthidx,3] = (filt_data[(filt_data.area_id .== area_id) .& 
                                                            (filt_data.monthidx .== monthidx),:].mitn_access_upper)[1]                                                

            mitn_subnat_use[admin1_i,monthidx,1] = (filt_data[(filt_data.area_id .== area_id) .& 
                                                            (filt_data.monthidx .== monthidx),:].mitn_use_lower)[1]
            mitn_subnat_use[admin1_i,monthidx,2] = (filt_data[(filt_data.area_id .== area_id) .& 
                                                            (filt_data.monthidx .== monthidx),:].mitn_use_mean)[1]
            mitn_subnat_use[admin1_i,monthidx,3] = (filt_data[(filt_data.area_id .== area_id) .& 
                                                            (filt_data.monthidx .== monthidx),:].mitn_use_upper)[1]
        end
    end
    
    # %% Save model predictions to .jld2 file
    jldsave(output_dir*"$(ISO)_subnat_model_coverage.jld2";
            filt_ISOs,
            YEAR_LIST,
            admin1_area_ids,
            admin1_names,
            mitn_subnat_npc,
            mitn_subnat_access,
            mitn_subnat_use,
            bv_subnat_npc,
            bv_subnat_use)
end
