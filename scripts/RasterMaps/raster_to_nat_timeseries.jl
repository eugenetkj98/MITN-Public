"""
Author: Eugene Tan
Date Created: 16/1/2025
Last Updated: 18/3/2025
Uses generated rasters to construct timeseries population weighted average estimates for ITN coverage.
Also extracts the model estimates from the previous BV model and compiles into a master .jld2 dict file
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
# Population Rasters directory
pop_dir = POPULATION_RASTER_DIR

# Directory of rasters
raster_dir = OUTPUT_RASTERS_DIR

# Output save dir
output_dir = OUTPUT_DIR*"coverage_timeseries/"
mkpath(output_dir)


# %% Perform draws and save outputs. Filter out unwanted countries
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = ["CPV","ZAF"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Time bounds
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

# %% Construct Timeseries for model predictions

# Raster Directory
raster_dir = OUTPUT_RASTERS_DIR

# Population Raster directory
pop_dir = POPULATION_RASTER_DIR

# Region boundaries
admin0_shapes_geoIO = GeoIO.load(ADMIN0_SHAPEFILE)

# Year bounds
YEAR_LIST = YEAR_START:YEAR_END
n_years = length(YEAR_LIST)

# %% Define Storage variables for timeseries with confidence 95% intervals
mitn_continent_npc = Array{Float64}(undef, n_years*12, 3)
mitn_continent_access = Array{Float64}(undef, n_years*12, 3)
mitn_continent_use = Array{Float64}(undef, n_years*12, 3)

mitn_nat_npc = Array{Float64}(undef, length(filt_ISOs), n_years*12, 3)
mitn_nat_access = Array{Float64}(undef, length(filt_ISOs), n_years*12, 3)
mitn_nat_use = Array{Float64}(undef, length(filt_ISOs), n_years*12, 3)

# %% For each year, extract required population weighted estimates from rasters (Extraction for MITN)
for year_i in ProgressBar(1:n_years, leave = false)
	for month in ProgressBar(1:12, leave = false)
		year = YEAR_LIST[year_i]

		# Calculate reference monthidx to access data
		monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEAR_START)

		# Get correct year/month string for importing file
		year_str = "$(year)"
		month_str = "$(month)"
		if month < 10
			month_str = "0"*month_str
		end

		# Import ITN coverage rasters (MITN)
		npc_mean_raster = replace_missing(Raster(raster_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_mean.tif"), missingval = NaN)
		npc_upper_raster = replace_missing(Raster(raster_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_upper.tif"), missingval = NaN)
		npc_lower_raster = replace_missing(Raster(raster_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_lower.tif"), missingval = NaN)

		access_mean_raster = replace_missing(Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_mean.tif"), missingval = NaN)
		access_upper_raster = replace_missing(Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_upper.tif"), missingval = NaN)
		access_lower_raster = replace_missing(Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_lower.tif"), missingval = NaN)
		
		use_mean_raster = replace_missing(Raster(raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_mean.tif"), missingval = NaN)
		use_upper_raster = replace_missing(Raster(raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_upper.tif"), missingval = NaN)
		use_lower_raster = replace_missing(Raster(raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_lower.tif"), missingval = NaN)
		# use_mean_raster = Raster("outputs/rasters/final_use/logis_use/use_$(year)_$(month_str)_mean.tif")
		# use_upper_raster = Raster("outputs/rasters/final_use/logis_use/use_$(year)_$(month_str)_upper.tif")
		# use_lower_raster = Raster("outputs/rasters/final_use/logis_use/use_$(year)_$(month_str)_lower.tif")

		# Import population raster
		pop_year = min(max(year, 2000), 2020)
		population_raster = replace_missing(Raster(pop_dir*"WorldPop_UNAdj_v3_DRC_fix.$(pop_year).Annual.Data.5km.sum.tif"), missingval = NaN)
		# Resample and align population raster to ITN Coverage rasters
		population_raster = resample(population_raster, to = npc_mean_raster)

		####
		# Calculate continent measures of ITN Coverage
		####
		# Get all locations with non-missing values and calculate population weighted averages
		# TEMP FIX: Calculate separate nonmissing_idxs for npc, access and use (INLA transformation function didn't deal with access = 0 and 1 very well)
		nonmissing_idxs_npc = intersect(findall(.!isnan.(population_raster)), findall(.!isnan.(npc_mean_raster)))
		mitn_continent_npc[monthidx,1] = sum(population_raster[nonmissing_idxs_npc].*npc_lower_raster[nonmissing_idxs_npc])/sum(population_raster[nonmissing_idxs_npc])
		mitn_continent_npc[monthidx,2] = sum(population_raster[nonmissing_idxs_npc].*npc_mean_raster[nonmissing_idxs_npc])/sum(population_raster[nonmissing_idxs_npc])
		mitn_continent_npc[monthidx,3] = sum(population_raster[nonmissing_idxs_npc].*npc_upper_raster[nonmissing_idxs_npc])/sum(population_raster[nonmissing_idxs_npc])

		nonmissing_idxs_access = intersect(findall(.!isnan.(population_raster)), findall(.!isnan.(access_mean_raster)))
		mitn_continent_access[monthidx,1] = sum(population_raster[nonmissing_idxs_access].*access_lower_raster[nonmissing_idxs_access])/sum(population_raster[nonmissing_idxs_access])
		mitn_continent_access[monthidx,2] = sum(population_raster[nonmissing_idxs_access].*access_mean_raster[nonmissing_idxs_access])/sum(population_raster[nonmissing_idxs_access])
		mitn_continent_access[monthidx,3] = sum(population_raster[nonmissing_idxs_access].*access_upper_raster[nonmissing_idxs_access])/sum(population_raster[nonmissing_idxs_access])

		nonmissing_idxs_use = intersect(findall(.!isnan.(population_raster)), findall(.!isnan.(use_mean_raster)))
		mitn_continent_use[monthidx,1] = sum(population_raster[nonmissing_idxs_use].*use_lower_raster[nonmissing_idxs_use])/sum(population_raster[nonmissing_idxs_use])
		mitn_continent_use[monthidx,2] = sum(population_raster[nonmissing_idxs_use].*use_mean_raster[nonmissing_idxs_use])/sum(population_raster[nonmissing_idxs_use])
		mitn_continent_use[monthidx,3] = sum(population_raster[nonmissing_idxs_use].*use_upper_raster[nonmissing_idxs_use])/sum(population_raster[nonmissing_idxs_use])

		####
		# Calculate national measures of ITN Coverage
		####
		for ISO_i in ProgressBar(1:length(filt_ISOs), leave = false)
			ISO = filt_ISOs[ISO_i]

			# Get required geometry
			admin0_geometry = admin0_shapes_geoIO[admin0_shapes_geoIO.ISO .== ISO,:].geometry

			# Get trimmed subset of population raster for target country
			pop_nat_masked = Rasters.trim(mask(population_raster, with = admin0_geometry); pad=0)

			# Get trimmed subset of ITN coverage rasters and align with population raster
			npc_mean_nat_masked = resample(Rasters.trim(mask(npc_mean_raster, with = admin0_geometry); pad=0), to = pop_nat_masked)
			npc_upper_nat_masked = resample(Rasters.trim(mask(npc_upper_raster, with = admin0_geometry); pad=0), to = pop_nat_masked)
			npc_lower_nat_masked = resample(Rasters.trim(mask(npc_lower_raster, with = admin0_geometry); pad=0), to = pop_nat_masked)

			access_mean_nat_masked = resample(Rasters.trim(mask(access_mean_raster, with = admin0_geometry); pad=0), to = pop_nat_masked)
			access_upper_nat_masked = resample(Rasters.trim(mask(access_upper_raster, with = admin0_geometry); pad=0), to = pop_nat_masked)
			access_lower_nat_masked = resample(Rasters.trim(mask(access_lower_raster, with = admin0_geometry); pad=0), to = pop_nat_masked)

			use_mean_nat_masked = resample(Rasters.trim(mask(use_mean_raster, with = admin0_geometry); pad=0), to = pop_nat_masked)
			use_upper_nat_masked = resample(Rasters.trim(mask(use_upper_raster, with = admin0_geometry); pad=0), to = pop_nat_masked)
			use_lower_nat_masked = resample(Rasters.trim(mask(use_lower_raster, with = admin0_geometry); pad=0), to = pop_nat_masked)

			# Get all locations with non-missing values and calculate population weighted averages
			# TEMP FIX: Calculate separate nonmissing_idxs for npc, access and use (INLA transformation function didn't deal with access = 0 and 1 very well)
			nonmissing_idxs_npc = intersect(findall(.!isnan.(pop_nat_masked)), findall(.!isnan.(npc_mean_nat_masked)))
			mitn_nat_npc[ISO_i,monthidx,1] = sum(pop_nat_masked[nonmissing_idxs_npc].*npc_lower_nat_masked[nonmissing_idxs_npc])/sum(pop_nat_masked[nonmissing_idxs_npc])
			mitn_nat_npc[ISO_i,monthidx,2] = sum(pop_nat_masked[nonmissing_idxs_npc].*npc_mean_nat_masked[nonmissing_idxs_npc])/sum(pop_nat_masked[nonmissing_idxs_npc])
			mitn_nat_npc[ISO_i,monthidx,3] = sum(pop_nat_masked[nonmissing_idxs_npc].*npc_upper_nat_masked[nonmissing_idxs_npc])/sum(pop_nat_masked[nonmissing_idxs_npc])

			nonmissing_idxs_access = intersect(findall(.!isnan.(pop_nat_masked)), findall(.!isnan.(access_mean_nat_masked)))
			mitn_nat_access[ISO_i,monthidx,1] = sum(pop_nat_masked[nonmissing_idxs_access].*access_lower_nat_masked[nonmissing_idxs_access])/sum(pop_nat_masked[nonmissing_idxs_access])
			mitn_nat_access[ISO_i,monthidx,2] = sum(pop_nat_masked[nonmissing_idxs_access].*access_mean_nat_masked[nonmissing_idxs_access])/sum(pop_nat_masked[nonmissing_idxs_access])
			mitn_nat_access[ISO_i,monthidx,3] = sum(pop_nat_masked[nonmissing_idxs_access].*access_upper_nat_masked[nonmissing_idxs_access])/sum(pop_nat_masked[nonmissing_idxs_access])

			nonmissing_idxs_use = intersect(findall(.!isnan.(pop_nat_masked)), findall(.!isnan.(use_mean_nat_masked)))
			mitn_nat_use[ISO_i,monthidx,1] = sum(pop_nat_masked[nonmissing_idxs_use].*use_lower_nat_masked[nonmissing_idxs_use])/sum(pop_nat_masked[nonmissing_idxs_use])
			mitn_nat_use[ISO_i,monthidx,2] = sum(pop_nat_masked[nonmissing_idxs_use].*use_mean_nat_masked[nonmissing_idxs_use])/sum(pop_nat_masked[nonmissing_idxs_use])
			mitn_nat_use[ISO_i,monthidx,3] = sum(pop_nat_masked[nonmissing_idxs_use].*use_upper_nat_masked[nonmissing_idxs_use])/sum(pop_nat_masked[nonmissing_idxs_use])
		end
	end
end

# %% BV Data extraction

# Import BV model data
bv_outputs_dir = BV_OUTPUT_DIR*"aggregated/"

# Declare storage variables (data stored in quarters)

bv_population = Array{Float64}(undef, length(filt_ISOs), n_years*12)
bv_nat_npc = Array{Float64}(undef, length(filt_ISOs), n_years*12, 3)
bv_nat_access = Array{Float64}(undef, length(filt_ISOs), n_years*12, 3)
bv_nat_use = Array{Float64}(undef, length(filt_ISOs), n_years*12, 3)

for year in ProgressBar(YEAR_START:YEAR_END, leave = false)
	# Import BV model output file for given year
	bv_filename = "aggregated_predictions_$(year).csv"
	bv_data = CSV.read(bv_outputs_dir*bv_filename, DataFrame)
	for month in 1:12
		monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEAR_START)
		Threads.@threads for ISO_i in 1:length(filt_ISOs)
			ISO = filt_ISOs[ISO_i]

			filt_data = bv_data[(bv_data.month .== "$(month)") .* (bv_data.iso3 .== ISO),:]
			if size(filt_data)[1] == 0# Country didn't have a prediction
				bv_nat_npc[ISO_i, monthidx,1] = NaN
				bv_nat_npc[ISO_i, monthidx,2] = NaN
				bv_nat_npc[ISO_i, monthidx,3] = NaN
				
				bv_nat_access[ISO_i, monthidx,1] = NaN
				bv_nat_access[ISO_i, monthidx,2] = NaN
				bv_nat_access[ISO_i, monthidx,3] = NaN

				bv_nat_use[ISO_i, monthidx,1] = NaN
				bv_nat_use[ISO_i, monthidx,2] = NaN
				bv_nat_use[ISO_i, monthidx,3] = NaN
			else
				bv_nat_npc[ISO_i, monthidx,1] = filt_data[filt_data.variable .== "percapita_nets","lower"][1]
				bv_nat_npc[ISO_i, monthidx,2] = filt_data[filt_data.variable .== "percapita_nets","mean"][1]
				bv_nat_npc[ISO_i, monthidx,3] = filt_data[filt_data.variable .== "percapita_nets","upper"][1]
				
				bv_nat_access[ISO_i, monthidx,1] = filt_data[filt_data.variable .== "access","lower"][1]
				bv_nat_access[ISO_i, monthidx,2] = filt_data[filt_data.variable .== "access","mean"][1]
				bv_nat_access[ISO_i, monthidx,3] = filt_data[filt_data.variable .== "access","upper"][1]

				bv_nat_use[ISO_i, monthidx,1] = filt_data[filt_data.variable .== "use","lower"][1]
				bv_nat_use[ISO_i, monthidx,2] = filt_data[filt_data.variable .== "use","mean"][1]
				bv_nat_use[ISO_i, monthidx,3] = filt_data[filt_data.variable .== "use","upper"][1]
			end

		end
	end
end

# %% Save model predictions to .jld2 file
jldsave(output_dir*"adj_nat_model_coverage.jld2";
        filt_ISOs,
        YEAR_LIST,
        mitn_continent_npc,
        mitn_continent_access,
        mitn_continent_use,
        mitn_nat_npc,
        mitn_nat_access,
        mitn_nat_use,
		bv_nat_npc,
		bv_nat_access,
		bv_nat_use)
