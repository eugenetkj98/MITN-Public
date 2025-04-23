"""
Author: Eugene Tan
Date Created: 4/11/2024
Last Updated: 4/11/2024
Code to post-process household survey data and add subnational area IDs
Saves output as a subnational version of itn_hh_surveydata_complete_subnat.csv
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %%
using GeoIO
using GeoStats
using GeoInterface
using ProgressBars
using CSV
using DataFrames
using Missings

# %% Import datasets

# Map geometries
dataset = GeoIO.load(ADMIN1_SHAPEFILE)

# Household Survey Entries
full_survey_data = CSV.read(RAW_DATASET_DIR*HOUSEHOLD_SURVEY_DATA_FILENAME, DataFrame)
nonmissing_idx = findall(.!ismissing.(full_survey_data[:,"latitude"]))

# %% Test for intersection and extract area_ids
area_ids = missings(Int, size(full_survey_data)[1])

Threads.@threads for i in ProgressBar(1:length(nonmissing_idx))

    idx = nonmissing_idx[i]

    # Get country ISO and relevant boundary data from shapefile
    ISO = full_survey_data[idx,:].ISO
    country_dataset = dataset[dataset.ISO.==ISO,:]
    n_admin1regions = length(country_dataset[:,1])
    
    # Get mercator projection of survey point (lat, lon)
    lat,lon = full_survey_data[idx,["latitude", "longitude"]]
    survey_point = Proj(Mercator).(Meshes.Point(LatLon(lat,lon)))

    # Begin search to test for intersection
    in_region_bool = false
    region_idx = -1

    for j in 1:n_admin1regions
        admin1_region = Proj(Mercator)(country_dataset[j,:].geometry)
        
        in_region_bool = intersects(survey_point,admin1_region)

        # Intersection found, end search
        if in_region_bool
            region_idx = j
            break
        end
    end

    # If there was a valid membership/intersection found, extract area_id
    if in_region_bool && (region_idx > 0)
        area_ids[idx] = country_dataset[region_idx,"area_id"]
    end
end


# %% Write and save subnat dataset
subnat_full_survey_data = hcat(full_survey_data, DataFrame(area_id = area_ids))
CSV.write(OUTPUT_SUBNAT_DATAPREP_DIR*HOUSEHOLD_SUBNAT_SURVEY_DATA_FILENAME, subnat_full_survey_data)



