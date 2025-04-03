"""
Author: Eugene Tan
Date Created: 17/3/2025
Last Updated: 17/3/2025
This script takes into the cleaned survey entries from the data engineering team and attaches area_ids
"""

# %% Declare required paths
push!(LOAD_PATH,"src/") # Path for source files and modules

# %% Activate required environment
# Package manager (This should install the required packages. If not, run the "instantiate" command via Pkg)
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Data Wrangling
using CSV
using DataFrames
using GeoIO
using GeoStats
using Missings
using ProgressBars

# %% Define Directories
datasets_dir = HOUSEHOLD_SURVEY_DIR
survey_data_filename = HOUSEHOLD_SURVEY_FILENAME

# %% Define output file name
output_dir = OUTPUT_DATAPREP_DIR
output_filename = HOUSEHOLD_SURVEY_DATA_FILENAME

# %% Geometry Data to lookup area id
# Region boundaries
admin1_shapes_geoIO = GeoIO.load(ADMIN1_SHAPEFILE)

#############################
# %% Data Eng Processing
#############################
# Load Data
dataeng_survey_data = CSV.read(datasets_dir*survey_data_filename, DataFrame)

# Get all entries that are nonmissing
nonmissing_idx = findall(.!ismissing.(dataeng_survey_data[:,"latitude"]))

# %% Test for intersection and extract area_ids
area_ids = missings(Int, size(dataeng_survey_data)[1])

Threads.@threads for i in ProgressBar(1:length(nonmissing_idx))
    idx = nonmissing_idx[i]

    # Get country ISO and relevant boundary data from shapefile
    ISO = dataeng_survey_data[idx,:].ISO
    country_dataset = admin1_shapes_geoIO[admin1_shapes_geoIO.ISO.==ISO,:]
    n_admin1regions = length(country_dataset[:,1])
    
    # Get mercator projection of survey point (lat, lon)
    lat,lon = dataeng_survey_data[idx,["latitude", "longitude"]]
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

# %% Write and save cleaned and standardised survey dataset
cleaned_full_survey_data = hcat(dataeng_survey_data[:,setdiff(names(dataeng_survey_data),["area_id"])], DataFrame(area_id = area_ids))
CSV.write(output_dir*output_filename, cleaned_full_survey_data)

println("JULIA EXTRACTION COMPLETED :D")
flush(stdout)