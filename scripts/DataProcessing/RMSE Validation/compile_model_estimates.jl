"""
Author: Eugene Tan
Date Created: 3/2/2025
Last Updated: 16/3/2025
Compiles input INLA dataset of household surveys, and model estimates from BV and MITN models to get .csv file of point 
estimates, used to calculate RMSE
"""
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import packages
using ProgressBars
using CSV
using DataFrames
using GeoIO
using Rasters
using JLD2
using NetAccessModel

using DateConversions

# %% Define raster file locations
map_rasters_dir = "Z:/eugene/Final Rasters/rasters/"#"outputs/rasters/"
pop_rasters_dir = "Z:/mastergrids/Other_Global_Covariates/Population/WorldPop/v3/PopulationCounts_DRC_fixed/5km/"

# %% Define dataset to post-process
dataset_dir = "datasets/"
inla_dataset_filename = "INLA/inla_dataset.csv"

# %% Get region boundaries
# Region boundaries
admin0_shapes_geoIO = GeoIO.load("Z:/master_geometries/Admin_Units/Global/MAP/2023/MG_5K/admin2023_0_MG_5K.shp")
admin1_shapes_geoIO = GeoIO.load("Z:/master_geometries/Admin_Units/Global/MAP/2023/MG_5K/admin2023_1_MG_5K.shp")

# %% Extracted time series from models (for National level estimates)
model_nat_est_timeseries = load("outputs/coverage_timeseries/nat_model_coverage.jld2")

# %% National level survey measurements
nat_survey_data = CSV.read("datasets/npc_monthly_data.csv", DataFrame)

# %% Reference years
YEAR_START = 2000
YEAR_END = 2023

# %% Import data
data = CSV.read(dataset_dir*inla_dataset_filename, DataFrame)

# Storage variable to store summary entries for final dataset
admin1_df_entries = []

# %% Extract data procedurally for each monthidx
monthidxs = sort(unique(data.monthidx))


# %%
for i in ProgressBar(1:length(monthidxs), leave = false)
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
    bv_npc_raster = replace_missing(Raster("Z:/map-data-eng/airflow/itn_model/output/itn_cube/20241015/nets_per_capita/ITN_$(year)_percapita_nets_mean.tif"), missingval = NaN)
    bv_use_raster = replace_missing(Raster("Z:/map-data-eng/airflow/itn_model/output/itn_cube/20241015/ITN_$(year)_use_mean.tif"), missingval = NaN)
    
    # Import MITN model rasters
    mitn_npc_raster = Raster(map_rasters_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_mean.tif")
    mitn_access_raster = Raster(map_rasters_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_mean.tif")
    mitn_use_raster = Raster(map_rasters_dir*"final_use/logis_use/use_$(year)_$(month_str)_mean.tif")
    # mitn_use_raster = Raster("outputs/rasters/final_use/logis_use/use_$(year)_$(month_str)_mean.tif")

    # Find all entries in relevant monthidx
    monthidx_filt_data = data[data.monthidx .== monthidx,:]

    # Get list of countries in given monthidx slice
    ISO_list = unique(monthidx_filt_data.ISO)

    # %% Extract summary metrics for each country
    for ISO_i in ProgressBar(1:length(ISO_list), leave = false)
        ISO = ISO_list[ISO_i]
        filt_admin0_data = monthidx_filt_data[monthidx_filt_data.ISO .== ISO,:]

        # Summarise data for all admin1 entries
        admin1_names = unique(filt_admin0_data.admin1_name)
        admin1_ids = unique(filt_admin0_data.area_id)
        n_admin1 = length(admin1_names)

        populations = Vector{Float64}(undef, n_admin1)

        for admin1_i in ProgressBar(1:n_admin1, leave = false)
            admin1_id = admin1_ids[admin1_i]
            data_slice = filt_admin0_data[(filt_admin0_data.area_id .== admin1_id),:]

            # Get population weights based on latitude and longitude
            pop_vals = Vector{Float64}(undef, size(data_slice)[1])
            for j in 1:size(data_slice)[1]
                lat, lon = data_slice[j,["latitude","longitude"]]
                pop_vals[j] = population_raster[X=Near(lon), Y=Near(lat)]
            end
            
            # Get weighted average survey estimates (weighted by cluster sample weight)
            npc = sum(data_slice.npc .* data_slice.cluster_sample_wt)/sum(data_slice.cluster_sample_wt)
            access = sum(data_slice.access .* data_slice.cluster_sample_wt)/sum(data_slice.cluster_sample_wt)
            use = sum(data_slice.use .* data_slice.cluster_sample_wt)/sum(data_slice.cluster_sample_wt)

            # Get region geographical boundaries
            admin1_geometry = admin1_shapes_geoIO[admin1_shapes_geoIO.area_id .== admin1_id,:].geometry

            # Get masked population
            pop_masked = Rasters.trim(mask(population_raster, with = admin1_geometry), pad = 0)
            
            # Get Subnational estimates from BV rasters (for subnational only have npc and use, no access)
            bv_npc_raster_masked = resample(Rasters.trim(mask(bv_npc_raster, with = admin1_geometry), pad = 0), to = pop_masked)
            bv_use_raster_masked = resample(Rasters.trim(mask(bv_use_raster, with = admin1_geometry); pad = 0), to = pop_masked)
            
            nonmissing_bv_npc_idxs = findall(.!isnan.(pop_masked) .& .!isnan.(bv_npc_raster_masked))
            nonmissing_bv_use_idxs = findall(.!isnan.(pop_masked) .& .!isnan.(bv_use_raster_masked))
            
            bv_npc = sum(pop_masked[nonmissing_bv_npc_idxs].*bv_npc_raster_masked[nonmissing_bv_npc_idxs])/sum(pop_masked[nonmissing_bv_npc_idxs])
            bv_use = sum(pop_masked[nonmissing_bv_use_idxs].*bv_use_raster_masked[nonmissing_bv_use_idxs])/sum(pop_masked[nonmissing_bv_use_idxs])

            # Get subnational estimates from MITN rasters for npc, access and use
            mitn_npc_raster_masked = resample(Rasters.trim(mask(mitn_npc_raster, with = admin1_geometry); pad = 0), to = pop_masked)
            mitn_access_raster_masked = resample(Rasters.trim(mask(mitn_access_raster, with = admin1_geometry); pad = 0), to = pop_masked)
            mitn_use_raster_masked = resample(Rasters.trim(mask(mitn_use_raster, with = admin1_geometry); pad = 0), to = pop_masked)

            nonmissing_mitn_npc_idxs = findall(.!isnan.(pop_masked) .& .!isnan.(mitn_npc_raster_masked))
            nonmissing_mitn_access_idxs = findall(.!isnan.(pop_masked) .& .!isnan.(mitn_access_raster_masked))
            nonmissing_mitn_use_idxs = findall(.!isnan.(pop_masked) .& .!isnan.(mitn_use_raster_masked))

            mitn_npc = sum(pop_masked[nonmissing_mitn_npc_idxs].*mitn_npc_raster_masked[nonmissing_mitn_npc_idxs])/sum(pop_masked[nonmissing_mitn_npc_idxs])
            mitn_access = sum(pop_masked[nonmissing_mitn_access_idxs].*mitn_access_raster_masked[nonmissing_mitn_access_idxs])/sum(pop_masked[nonmissing_mitn_access_idxs])
            mitn_use = sum(pop_masked[nonmissing_mitn_use_idxs].*mitn_use_raster_masked[nonmissing_mitn_use_idxs])/sum(pop_masked[nonmissing_mitn_use_idxs])

            output_df_entry = DataFrame(ISO = ISO,
                                        admin_level = 1,
                                        admin_name = admin1_names[admin1_i],
                                        admin_id = admin1_ids[admin1_i],
                                        interview_month = month,
                                        interview_year = year,
                                        monthidx = monthidx,
                                        population = sum(pop_vals),
                                        npc = npc,
                                        access = access,
                                        use = use,
                                        bv_npc = bv_npc,
                                        bv_access = NaN,
                                        bv_use = bv_use,
                                        mitn_npc = mitn_npc,
                                        mitn_access = mitn_access,
                                        mitn_use = mitn_use)
            push!(admin1_df_entries, output_df_entry)
        end
    end
end


# %%
admin0_df_entries = []
ISO_list = String.(CSV.read(raw"C:\Users\ETan\Documents\Prototype Analyses\itn-updated\datasets\ISO_list.csv", DataFrame)[:,1])
exclusion_ISOs = ["CPV","ZAF"]
YEAR_START = 2000
YEAR_END = 2023
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

for ISO_i in ProgressBar(1:length(filt_ISOs))
    ISO = filt_ISOs[ISO_i]
    # Import access data and calculate household access
    net_access_input_dict = load("outputs/extractions/access/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")

    # Calculate reference values for National household Access from Surveys
    national_H_aggregated = net_access_input_dict["H_aggregated"]
    national_access_aggregated = zeros(size(national_H_aggregated)[1])
    for survey_i in 1:length(national_access_aggregated)
        national_access_aggregated[survey_i] = sum(H_to_access(national_H_aggregated[survey_i,:,:]))
    end
    

    # Get list of monthidx to extract values from (those present in surveys)
    nat_monthidxs = monthyear_to_monthidx.(nat_survey_data[(nat_survey_data.ISO .== ISO) .& 
                                                            (nat_survey_data.year .>= YEAR_START) .&
                                                            (nat_survey_data.year .<= YEAR_END),:].month, 
                                            nat_survey_data[(nat_survey_data.ISO .== ISO) .& 
                                                    (nat_survey_data.year .>= YEAR_START) .&
                                                    (nat_survey_data.year .<= YEAR_END),:].year, YEAR_START = YEAR_START)
    access_monthidxs = net_access_input_dict["survey_monthidx_aggregated"]
    use_monthidxs = monthyear_to_monthidx.(nat_survey_data[(nat_survey_data.ISO .== ISO) .& 
                                                            (nat_survey_data.year .>= YEAR_START) .&
                                                            (nat_survey_data.year .<= YEAR_END),:].month, 
                                                            nat_survey_data[(nat_survey_data.ISO .== ISO) .& 
                                                            (nat_survey_data.year .>= YEAR_START) .&
                                                            (nat_survey_data.year .<= YEAR_END),:].year, YEAR_START = YEAR_START)

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
        bv_npc_raster = replace_missing(Raster("Z:/map-data-eng/airflow/itn_model/output/itn_cube/20241015/nets_per_capita/ITN_$(year)_percapita_nets_mean.tif"), missingval = NaN)
        bv_use_raster = replace_missing(Raster("Z:/map-data-eng/airflow/itn_model/output/itn_cube/20241015/ITN_$(year)_use_mean.tif"), missingval = NaN)

        # Import MITN model rasters
        mitn_npc_raster = Raster(map_rasters_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_mean.tif")
        mitn_access_raster = Raster(map_rasters_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_mean.tif")
        mitn_use_raster = Raster(map_rasters_dir*"final_use/logis_use/use_$(year)_$(month_str)_mean.tif")
        # mitn_use_raster = Raster("outputs/rasters/final_use/logis_use/use_$(year)_$(month_str)_mean.tif")

        # Find all entries in relevant monthidx
        monthidx_filt_data = data[data.monthidx .== monthidx,:]

        # Extract summary for admin0 region based on population weights
        admin0_npc = NaN
        admin0_access = NaN
        admin0_use = NaN

        if monthidx ∈ nat_monthidxs
            admin0_npc = nat_survey_data[(nat_survey_data.ISO .== ISO) .& 
                                            (nat_survey_data.month .== month) .&
                                            (nat_survey_data.year .== year),"NPC_mean"][1]
        end

        if monthidx ∈ access_monthidxs
            admin0_access = national_access_aggregated[findfirst(access_monthidxs .== monthidx)]
        end

        if monthidx ∈ use_monthidxs
            admin0_use = nat_survey_data[(nat_survey_data.ISO .== ISO) .& 
                                            (nat_survey_data.month .== month) .&
                                            (nat_survey_data.year .== year),"use_mean"][1]
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

        # Combine all entries and add to collection)
        push!(admin0_df_entries, output_df_entry)
    end
end

# %%
output_csv = vcat(admin1_df_entries..., admin0_df_entries...)

# %%

CSV.write("outputs/coverage_timeseries/model_prediction_comparisons.csv", output_csv)


