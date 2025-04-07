"""
Author: Eugene Tan
Date Created: 23/7/2024
Last Updated: 16/1/2025
Code to extract and prep input data to for SNF modelling
"""

module DataExtractions
export extract_data_netcrop
export SA_extract_data_netcrop
export extract_data_netaccess

# %% Import settings, filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Import required packages
# Data Wrangling
using CSV
using DataFrames
using JLD2
using Missings

# Maths packages
using LinearAlgebra
using StatsBase
using Random

# Parsing
using Unicode
using Dates

# Helper Functions
using DateConversions
using NetCropModel

########################################
# %% Define all filepaths.
# TO DO: CONSIDER INPUT AS A CONFIG TEXT FILE
########################################
# Filepaths
deliveries_data_dir = RAW_DATASET_DIR
distribution_data_dir = RAW_DATASET_DIR
population_data_dir = RAW_DATASET_DIR
hh_data_dir = RAW_DATASET_DIR
hh_dataprep_dir = OUTPUT_DATAPREP_DIR
countrycodes_dir = RAW_DATASET_DIR

# Filenames
deliveries_data_filename = DELIVERIES_DATA_FILENAME
distribution_data_filename = DISTRIBUTION_DATA_FILENAME
population_data_filename = POPULATION_DATA_FILENAME
hh_listing_filename = HOUSEHOLD_SURVEY_DATA_FILENAME
hh_npc_data = HOUSEHOLD_NAT_SUMMARY_DATA_FILENAME
countrycodes_filename = COUNTRY_CODES_FILENAME

# Output Directory
output_dir = OUTPUT_EXTRACTIONS_DIR

########################################
# %% Define extraction settings
########################################
# Window smoothing settings
R_window = SMOOTHING_WINDOW_WIDTH

########################################
# %% Data Extraction Code
########################################
"""
    extract_data_netcrop(ISO::String, YEAR_START::Int64, YEAR_END::Int64)
    
Read from pre-declared directories of CSV files and extract required raw data to perform model analysis on. Extraction code also handles the required upsampling, interpolation and averaging. Output is returned as a Dict. A copy of the extraction is also saved in the listed output directory as a JLD2 file with name format: e.g. `AGO_2010_2021_extract.jld2`
"""
function extract_data_netcrop(ISO::String, YEAR_START::Int64, YEAR_END::Int64;
                deliveries_data_dir=deliveries_data_dir, 
                distribution_data_dir=distribution_data_dir,
                population_data_dir=population_data_dir,
                hh_data_dir=hh_data_dir, 
                countrycodes_dir=countrycodes_dir,
                deliveries_data_filename=deliveries_data_filename,
                distribution_data_filename=distribution_data_filename,
                population_data_filename=population_data_filename,
                # hh_listing_filename=hh_listing_filename,
                hh_npc_data = hh_npc_data,
                countrycodes_filename=countrycodes_filename,
                R_window = R_window,
                save_output = true)

    ########################################
    # %% Import Country Code Conversion Table
    ########################################
    country_codes = CSV.read(countrycodes_dir*countrycodes_filename, DataFrame)
    DHS = country_codes.DHS_CODE[findfirst(country_codes.ISO3 .== ISO)]

    ########################################
    # %% Calculate time bounds for conversion from annual to monthly
    ########################################
    N_MONTHS = (YEAR_END-YEAR_START+1)*12
    
    ########################################
    # %% Declare storage variables
    ########################################
    YEARS_ANNUAL = Vector(YEAR_START:YEAR_END)
    DELIVERIES_ANNUAL = zeros(Int64, length(YEARS_ANNUAL))
    DISTRIBUTION_ANNUAL = missings(Int64, length(YEARS_ANNUAL)) # TEMPORARY SIZE, will be expanded to account for multiple nets later in the code
    # Needs to be extended by one year for interpolation puposes
    POPULATION_YEARS = Vector(YEAR_START:(YEAR_END+1))
    POPULATION_ANNUAL = zeros(Int64, length(POPULATION_YEARS))
    

    MONTHS_MONTHLY = Vector(1:N_MONTHS)
    HOUSEHOLD_cITNPC_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
    HOUSEHOLD_LLINPC_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
    HOUSEHOLD_NPC_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
    HOUSEHOLD_NPC_STD_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
    cITN_CROP_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
    LLIN_CROP_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
    NET_CROP_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
    NET_CROP_STD_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
    # NET_CROP_MONTHLY_SMOOTHED = missings(Float64, length(MONTHS_MONTHLY))
    POPULATION_MONTHLY = zeros(Int64, length(MONTHS_MONTHLY))

    ########################################
    # %% Deliveries Data
    ########################################

    # Import Data
    deliveries_data_RAW = CSV.read(deliveries_data_dir*deliveries_data_filename, DataFrame)
    deliveries_data_RAW_ISOSLICE = deliveries_data_RAW[findall(deliveries_data_RAW.iso3 .== ISO),:]

    # Parse Data into DataFrame and extract relevant values for annual delivery data
    col_names = names(deliveries_data_RAW_ISOSLICE)
    for i in 1:length(YEARS_ANNUAL)
        year = YEARS_ANNUAL[i]
        year_string = string(year)
        if year_string ∈ col_names
            col_idx = findfirst(names(deliveries_data_RAW_ISOSLICE) .== year_string)
            DELIVERIES_ANNUAL[i] = deliveries_data_RAW_ISOSLICE[1,col_idx]
        end
    end

    ########################################
    # %% Distributions Data
    ########################################
    # Import Data
    distribution_data_RAW = CSV.read(distribution_data_dir*distribution_data_filename, DataFrame)
    distribution_data_RAW_ISOSLICE = distribution_data_RAW[findall((distribution_data_RAW.iso .== ISO).*
                                                                (distribution_data_RAW.year .>= YEAR_START).*
                                                                (distribution_data_RAW.year .<= YEAR_END)),:]

    
    # Get list of Net Type col_names that have been previously distributed
    net_extract = Matrix(distribution_data_RAW_ISOSLICE[:,6:end])
    net_extract[findall(ismissing.(net_extract))] .= 0
    net_tally = sum(net_extract, dims = 1)[:]
    included_net_idxs = union([1,2],findall(net_tally.>0))
    NET_NAMES = names(distribution_data_RAW)[6:end][included_net_idxs]

    # Reshape storage size for annual distribution array to account for multiple nets
    DISTRIBUTION_ANNUAL = repeat(DISTRIBUTION_ANNUAL, 1,length(NET_NAMES)+1)

    
    
    # Parse Data into DataFrame and extract relevant values for annual delivery data
    for i in 1:length(YEARS_ANNUAL)
        year = YEARS_ANNUAL[i]
        year_vals = distribution_data_RAW_ISOSLICE.year
        if year ∈ year_vals
            row_idx = findfirst(year_vals .== year)
            entry = Vector(distribution_data_RAW_ISOSLICE[row_idx,vcat("Total Nets", NET_NAMES)])
            DISTRIBUTION_ANNUAL[i,:].= copy(entry)
        end
    end

    ########################################
    # %% Population Data
    ########################################

    # Import Data
    population_data_RAW = CSV.read(population_data_dir*population_data_filename, DataFrame)
    population_data_RAW_ISOSLICE = population_data_RAW[findall( 
                                                                (population_data_RAW.iso3 .== ISO).*
                                                                (population_data_RAW.admin_unit_level .== "ADMIN0").*
                                                                (population_data_RAW.age_bin .== "All_Ages")
                                                                ),:]

    # Parse Data into DataFrame and extract relevant values for annual delivery data
    for i in 1:(length(POPULATION_YEARS))
        year = POPULATION_YEARS[i]
        year_vals = population_data_RAW_ISOSLICE.year
        if year ∈ year_vals
            row_idx = findfirst(year_vals .== year)
            POPULATION_ANNUAL[i] = round(Int64,population_data_RAW_ISOSLICE[row_idx,"total_pop"])
        end
    end

    # Linearly Interpolate for monthly population data
    for i in 1:(length(POPULATION_ANNUAL)-1)
        POPULATION_MONTHLY[((i-1)*12+1):((i*12))] = round.(Int,LinRange(POPULATION_ANNUAL[i],POPULATION_ANNUAL[i+1],13)[1:12])
    end

    ########################################
    # %% Household Survey Data
    ########################################

    # Import NPC data
    household_npc_data = CSV.read(hh_dataprep_dir*hh_npc_data, DataFrame)

    # Filter by require country ISO
    household_npc_ISOSLICE = household_npc_data[findall((household_npc_data.ISO .== ISO).*
                                                        (household_npc_data.year .>= YEAR_START).*
                                                        (household_npc_data.year .<= YEAR_END)),:]

    
    # Filter by year range
    HOUSEHOLD_cITNPC_LISTS = Vector{Any}( undef, length(MONTHS_MONTHLY))
    HOUSEHOLD_LLINPC_LISTS = Vector{Any}( undef, length(MONTHS_MONTHLY))
    HOUSEHOLD_NPC_LISTS = Vector{Any}( undef, length(MONTHS_MONTHLY))
    HOUSEHOLD_NPC_STD_LISTS = Vector{Any}( undef, length(MONTHS_MONTHLY))
    for i in 1:length(HOUSEHOLD_NPC_LISTS)
        HOUSEHOLD_cITNPC_LISTS[i] = []
        HOUSEHOLD_LLINPC_LISTS[i] = []
        HOUSEHOLD_NPC_LISTS[i] = []
        HOUSEHOLD_NPC_STD_LISTS[i] = []
    end

    monthidxs = monthyear_to_monthidx.(household_npc_ISOSLICE.month,household_npc_ISOSLICE.year, YEAR_START = YEAR_START)

    for i in 1:size(household_npc_ISOSLICE)[1]
        month = monthidxs[i]
        normalised_cITNPC_mean = household_npc_ISOSLICE.NPC_mean[i]*(household_npc_ISOSLICE.cITNPC_mean[i]/(household_npc_ISOSLICE.cITNPC_mean[i]+household_npc_ISOSLICE.LLINPC_mean[i]))
        normalised_LLINPC_mean = household_npc_ISOSLICE.NPC_mean[i]*(household_npc_ISOSLICE.LLINPC_mean[i]/(household_npc_ISOSLICE.cITNPC_mean[i]+household_npc_ISOSLICE.LLINPC_mean[i]))
        push!(HOUSEHOLD_cITNPC_LISTS[month], household_npc_ISOSLICE.cITNPC_mean[i])
        push!(HOUSEHOLD_LLINPC_LISTS[month], household_npc_ISOSLICE.LLINPC_mean[i])
        push!(HOUSEHOLD_NPC_LISTS[month], household_npc_ISOSLICE.NPC_mean[i])
        push!(HOUSEHOLD_NPC_STD_LISTS[month], household_npc_ISOSLICE.NPC_adj_se[i])
    end

    valid_idx = findall(.!isempty.(HOUSEHOLD_NPC_LISTS))
    HOUSEHOLD_cITNPC_MONTHLY[valid_idx] = mean.(HOUSEHOLD_cITNPC_LISTS[valid_idx])
    HOUSEHOLD_LLINPC_MONTHLY[valid_idx] = mean.(HOUSEHOLD_LLINPC_LISTS[valid_idx])
    HOUSEHOLD_NPC_MONTHLY[valid_idx] = mean.(HOUSEHOLD_NPC_LISTS[valid_idx])
    HOUSEHOLD_NPC_STD_MONTHLY[valid_idx] = mean.(HOUSEHOLD_NPC_STD_LISTS[valid_idx])

    ########################################
    # %% Smooth out monthly NCP Data
    ########################################
    HOUSEHOLD_NPC_MONTHLY_SMOOTHED = missings(Float64, length(MONTHS_MONTHLY))
    for i in 1:length(MONTHS_MONTHLY)
        L_index = max(i-R_window,1)
        R_index = min(i,length(MONTHS_MONTHLY))
        window = HOUSEHOLD_NPC_MONTHLY[L_index:R_index]
        if !ismissing(HOUSEHOLD_NPC_MONTHLY[i])
            HOUSEHOLD_NPC_MONTHLY_SMOOTHED[i] = mean(window[findall(.!ismissing.(window))])
        end
    end

    ########################################
    # %% Estimate Monthly Net Crop based on population data NPC survey data
    ########################################
    # NET_CROP_MONTHLY_SMOOTHED = missings(Float64, length(MONTHS_MONTHLY))
    
    for i in 1:length(MONTHS_MONTHLY)
        if !ismissing(HOUSEHOLD_NPC_MONTHLY[i])
            cITN_CROP_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_cITNPC_MONTHLY[i]
            LLIN_CROP_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_LLINPC_MONTHLY[i]
            NET_CROP_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_NPC_MONTHLY[i]
            # NET_CROP_MONTHLY_SMOOTHED[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_NPC_MONTHLY_SMOOTHED[i]
            NET_CROP_STD_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_NPC_STD_MONTHLY[i]
        end
    end


    input_data_dict = Dict("ISO" => ISO,
                            "YEAR_START" => YEAR_START,
                            "YEAR_END" => YEAR_END,
                            "YEARS_ANNUAL" => YEARS_ANNUAL,
                            "DELIVERIES_ANNUAL" => DELIVERIES_ANNUAL,
                            "DISTRIBUTION_ANNUAL" => DISTRIBUTION_ANNUAL,
                            "POPULATION_ANNUAL" => POPULATION_ANNUAL,
                            "POPULATION_MONTHLY" => POPULATION_MONTHLY,
                            "MONTHS_MONTHLY" => MONTHS_MONTHLY,
                            "HOUSEHOLD_NPC_MONTHLY" => HOUSEHOLD_NPC_MONTHLY,
                            "HOUSEHOLD_NPC_STD_MONTHLY" => HOUSEHOLD_NPC_STD_MONTHLY,
                            "NET_CROP_MONTHLY" => NET_CROP_MONTHLY,
                            "NET_CROP_STD_MONTHLY" => NET_CROP_STD_MONTHLY,
                            "cITN_CROP_MONTHLY" => cITN_CROP_MONTHLY,
                            "LLIN_CROP_MONTHLY" => LLIN_CROP_MONTHLY,
                            # "NET_CROP_MONTHLY_SMOOTHED" => NET_CROP_MONTHLY_SMOOTHED,
                            "NET_NAMES" => NET_NAMES)
    
    ########################################
    # %% Save Extracted Data into JLD2 file
    ########################################

    # Make directory path to store data
    target_path = output_dir*"/crop/$(YEAR_START)_$(YEAR_END)/"
    mkpath(target_path)

    if save_output
        output_filename = (target_path)*(ISO*"_"*"$(YEAR_START)"*"_"*"$(YEAR_END)"*"_cropextract.jld2")
        save(output_filename, input_data_dict)
    end

    return input_data_dict
end

########################################
# %% Data Extraction Code
########################################
"""
    SA_extract_data_netcrop(ISO::String, YEAR_START::Int64, YEAR_END::Int64)
    
Extract data into appropriate jld2 files and saves them into the directory: outputs/extractions/crop/sesnsitivity_analysis/YEAR_START_YEAR_END/...
Extracted data is prepared as inputs for a sensitivity analysis with the following settings:
    - chronological: Procedurally include survey entries chronologically
    - random: Randomly include survey entries
    - cross: 1-fold cross validation
"""
function SA_extract_data_netcrop(ISO::String, YEAR_START::Int64, YEAR_END::Int64;
                deliveries_data_dir=deliveries_data_dir, 
                distribution_data_dir=distribution_data_dir,
                population_data_dir=population_data_dir,
                hh_data_dir=hh_data_dir, 
                countrycodes_dir=countrycodes_dir,
                deliveries_data_filename=deliveries_data_filename,
                distribution_data_filename=distribution_data_filename,
                population_data_filename=population_data_filename,
                # hh_listing_filename=hh_listing_filename,
                hh_npc_data = hh_npc_data,
                countrycodes_filename=countrycodes_filename,
                R_window = R_window,
                save_output = true)

    ########################################
    # %% Import Country Code Conversion Table
    ########################################
    country_codes = CSV.read(countrycodes_dir*countrycodes_filename, DataFrame)
    DHS = country_codes.DHS_CODE[findfirst(country_codes.ISO3 .== ISO)]

    ########################################
    # %% Calculate time bounds for conversion from annual to monthly
    ########################################
    N_MONTHS = (YEAR_END-YEAR_START+1)*12
    
    ########################################
    # %% Declare storage variables
    ########################################
    YEARS_ANNUAL = Vector(YEAR_START:YEAR_END)
    DELIVERIES_ANNUAL = zeros(Int64, length(YEARS_ANNUAL))
    DISTRIBUTION_ANNUAL = missings(Int64, length(YEARS_ANNUAL)) # TEMPORARY SIZE, will be expanded to account for multiple nets later in the code
    # Needs to be extended by one year for interpolation puposes
    POPULATION_YEARS = Vector(YEAR_START:(YEAR_END+1))
    POPULATION_ANNUAL = zeros(Int64, length(POPULATION_YEARS))
    

    MONTHS_MONTHLY = Vector(1:N_MONTHS)
    # NET_CROP_MONTHLY_SMOOTHED = missings(Float64, length(MONTHS_MONTHLY))
    POPULATION_MONTHLY = zeros(Int64, length(MONTHS_MONTHLY))

    ########################################
    # %% Deliveries Data
    ########################################

    # Import Data
    deliveries_data_RAW = CSV.read(deliveries_data_dir*deliveries_data_filename, DataFrame)
    deliveries_data_RAW_ISOSLICE = deliveries_data_RAW[findall(deliveries_data_RAW.iso3 .== ISO),:]

    # Parse Data into DataFrame and extract relevant values for annual delivery data
    col_names = names(deliveries_data_RAW_ISOSLICE)
    for i in 1:length(YEARS_ANNUAL)
        year = YEARS_ANNUAL[i]
        year_string = string(year)
        if year_string ∈ col_names
            col_idx = findfirst(names(deliveries_data_RAW_ISOSLICE) .== year_string)
            DELIVERIES_ANNUAL[i] = deliveries_data_RAW_ISOSLICE[1,col_idx]
        end
    end

    ########################################
    # %% Distributions Data
    ########################################
    # Import Data
    distribution_data_RAW = CSV.read(distribution_data_dir*distribution_data_filename, DataFrame)
    distribution_data_RAW_ISOSLICE = distribution_data_RAW[findall((distribution_data_RAW.iso .== ISO).*
                                                                (distribution_data_RAW.year .>= YEAR_START).*
                                                                (distribution_data_RAW.year .<= YEAR_END)),:]

    
    # Get list of Net Type col_names that have been previously distributed
    net_extract = Matrix(distribution_data_RAW_ISOSLICE[:,6:end])
    net_extract[findall(ismissing.(net_extract))] .= 0
    net_tally = sum(net_extract, dims = 1)[:]
    included_net_idxs = union([1,2],findall(net_tally.>0))
    NET_NAMES = names(distribution_data_RAW)[6:end][included_net_idxs]

    # Reshape storage size for annual distribution array to account for multiple nets
    DISTRIBUTION_ANNUAL = repeat(DISTRIBUTION_ANNUAL, 1,length(NET_NAMES)+1)

    
    
    # Parse Data into DataFrame and extract relevant values for annual delivery data
    for i in 1:length(YEARS_ANNUAL)
        year = YEARS_ANNUAL[i]
        year_vals = distribution_data_RAW_ISOSLICE.year
        if year ∈ year_vals
            row_idx = findfirst(year_vals .== year)
            entry = Vector(distribution_data_RAW_ISOSLICE[row_idx,vcat("Total Nets", NET_NAMES)])
            DISTRIBUTION_ANNUAL[i,:].= copy(entry)
        end
    end

    ########################################
    # %% Population Data
    ########################################

    # Import Data
    population_data_RAW = CSV.read(population_data_dir*population_data_filename, DataFrame)
    population_data_RAW_ISOSLICE = population_data_RAW[findall( 
                                                                (population_data_RAW.iso3 .== ISO).*
                                                                (population_data_RAW.admin_unit_level .== "ADMIN0").*
                                                                (population_data_RAW.age_bin .== "All_Ages")
                                                                ),:]

    # Parse Data into DataFrame and extract relevant values for annual delivery data
    for i in 1:(length(POPULATION_YEARS))
        year = POPULATION_YEARS[i]
        year_vals = population_data_RAW_ISOSLICE.year
        if year ∈ year_vals
            row_idx = findfirst(year_vals .== year)
            POPULATION_ANNUAL[i] = round(Int64,population_data_RAW_ISOSLICE[row_idx,"total_pop"])
        end
    end

    # Linearly Interpolate for monthly population data
    for i in 1:(length(POPULATION_ANNUAL)-1)
        POPULATION_MONTHLY[((i-1)*12+1):((i*12))] = round.(Int,LinRange(POPULATION_ANNUAL[i],POPULATION_ANNUAL[i+1],13)[1:12])
    end

    ########################################
    # %% Household Survey Data
    ########################################

    # Import NPC data
    household_npc_data = CSV.read(hh_dataprep_dir*hh_npc_data, DataFrame)

    # Filter by require country ISO
    household_npc_ISOSLICE = household_npc_data[findall((household_npc_data.ISO .== ISO).*
                                                        (household_npc_data.year .>= YEAR_START).*
                                                        (household_npc_data.year .<= YEAR_END)),:]

    # Get list of unique survey ids to filter for sensitivity analysis
    survey_names = unique(household_npc_ISOSLICE.SurveyId[sortperm(household_npc_ISOSLICE.year),:])
    n_surveys = length(survey_names)

    # Storage variables to survey extracts
    # Chronological inclusion
    CHRONO_SA_HOUSEHOLD_cITNPC_LISTS = Array{Any}( undef, n_surveys, length(MONTHS_MONTHLY))
    CHRONO_SA_HOUSEHOLD_LLINPC_LISTS = Array{Any}( undef, n_surveys, length(MONTHS_MONTHLY))
    CHRONO_SA_HOUSEHOLD_NPC_LISTS = Array{Any}( undef, n_surveys, length(MONTHS_MONTHLY))
    CHRONO_SA_HOUSEHOLD_NPC_STD_LISTS = Array{Any}( undef, n_surveys, length(MONTHS_MONTHLY))
    for n in 1:length(n_surveys)
        for i in 1:length(CHRONO_SA_HOUSEHOLD_NPC_LISTS)
            CHRONO_SA_HOUSEHOLD_cITNPC_LISTS[i] = []
            CHRONO_SA_HOUSEHOLD_LLINPC_LISTS[i] = []
            CHRONO_SA_HOUSEHOLD_NPC_LISTS[i] = []
            CHRONO_SA_HOUSEHOLD_NPC_STD_LISTS[i] = []
        end
    end

    # Random inclusion
    RANDOM_SA_HOUSEHOLD_cITNPC_LISTS = Array{Any}( undef, n_surveys, length(MONTHS_MONTHLY))
    RANDOM_SA_HOUSEHOLD_LLINPC_LISTS = Array{Any}( undef, n_surveys, length(MONTHS_MONTHLY))
    RANDOM_SA_HOUSEHOLD_NPC_LISTS = Array{Any}( undef, n_surveys, length(MONTHS_MONTHLY))
    RANDOM_SA_HOUSEHOLD_NPC_STD_LISTS = Array{Any}( undef, n_surveys, length(MONTHS_MONTHLY))
    for n in 1:length(n_surveys)
        for i in 1:length(RANDOM_SA_HOUSEHOLD_NPC_LISTS)
            RANDOM_SA_HOUSEHOLD_cITNPC_LISTS[i] = []
            RANDOM_SA_HOUSEHOLD_LLINPC_LISTS[i] = []
            RANDOM_SA_HOUSEHOLD_NPC_LISTS[i] = []
            RANDOM_SA_HOUSEHOLD_NPC_STD_LISTS[i] = []
        end
    end

    # 1-fold cross validation
    CROSS_SA_HOUSEHOLD_cITNPC_LISTS = Array{Any}( undef, n_surveys, length(MONTHS_MONTHLY))
    CROSS_SA_HOUSEHOLD_LLINPC_LISTS = Array{Any}( undef, n_surveys, length(MONTHS_MONTHLY))
    CROSS_SA_HOUSEHOLD_NPC_LISTS = Array{Any}( undef, n_surveys, length(MONTHS_MONTHLY))
    CROSS_SA_HOUSEHOLD_NPC_STD_LISTS = Array{Any}( undef, n_surveys, length(MONTHS_MONTHLY))
    for n in 1:length(n_surveys)
        for i in 1:length(CROSS_SA_HOUSEHOLD_NPC_LISTS)
            CROSS_SA_HOUSEHOLD_cITNPC_LISTS[i] = []
            CROSS_SA_HOUSEHOLD_LLINPC_LISTS[i] = []
            CROSS_SA_HOUSEHOLD_NPC_LISTS[i] = []
            CROSS_SA_HOUSEHOLD_NPC_STD_LISTS[i] = []
        end
    end

    # Scenario 1: Chronological Filtering (Procedurally include surveys chronologically i.e. add one by one)
    for n in 1:n_surveys
        # Filter out on a range of desired surveys
        incl_survey_names = survey_names[1:n]

        membership_bool = zeros(Bool, length(household_npc_ISOSLICE.SurveyId))
        for i in 1:length(membership_bool)
            if household_npc_ISOSLICE.SurveyId[i] ∈ incl_survey_names
                membership_bool[i] = 1
            end
        end

        CHRONO_SA_FILTER_household_npc_ISOSLICE = household_npc_ISOSLICE[findall(membership_bool),:]

        monthidxs = monthyear_to_monthidx.(CHRONO_SA_FILTER_household_npc_ISOSLICE.month,CHRONO_SA_FILTER_household_npc_ISOSLICE.year, YEAR_START = YEAR_START)

        for i in 1:size(CHRONO_SA_FILTER_household_npc_ISOSLICE)[1]
            month = monthidxs[i]
            normalised_cITNPC_mean = CHRONO_SA_FILTER_household_npc_ISOSLICE.NPC_mean[i]*(CHRONO_SA_FILTER_household_npc_ISOSLICE.cITNPC_mean[i]/(CHRONO_SA_FILTER_household_npc_ISOSLICE.cITNPC_mean[i]+CHRONO_SA_FILTER_household_npc_ISOSLICE.LLINPC_mean[i]))
            normalised_LLINPC_mean = CHRONO_SA_FILTER_household_npc_ISOSLICE.NPC_mean[i]*(CHRONO_SA_FILTER_household_npc_ISOSLICE.LLINPC_mean[i]/(CHRONO_SA_FILTER_household_npc_ISOSLICE.cITNPC_mean[i]+CHRONO_SA_FILTER_household_npc_ISOSLICE.LLINPC_mean[i]))
            push!(CHRONO_SA_HOUSEHOLD_cITNPC_LISTS[n,month], CHRONO_SA_FILTER_household_npc_ISOSLICE.cITNPC_mean[i])
            push!(CHRONO_SA_HOUSEHOLD_LLINPC_LISTS[n,month], CHRONO_SA_FILTER_household_npc_ISOSLICE.LLINPC_mean[i])
            push!(CHRONO_SA_HOUSEHOLD_NPC_LISTS[n,month], CHRONO_SA_FILTER_household_npc_ISOSLICE.NPC_mean[i])
            push!(CHRONO_SA_HOUSEHOLD_NPC_STD_LISTS[n,month], CHRONO_SA_FILTER_household_npc_ISOSLICE.NPC_adj_se[i])
        end
    end

    # Scenario 2: Random inclusion ordering
    rand_survey_names = shuffle(survey_names)
    for n in 1:n_surveys
        # Filter out on a range of desired surveys
        incl_survey_names = rand_survey_names[1:n]

        membership_bool = zeros(Bool, length(household_npc_ISOSLICE.SurveyId))
        for i in 1:length(membership_bool)
            if household_npc_ISOSLICE.SurveyId[i] ∈ incl_survey_names
                membership_bool[i] = 1
            end
        end

        RANDOM_SA_FILTER_household_npc_ISOSLICE = household_npc_ISOSLICE[findall(membership_bool),:]

        monthidxs = monthyear_to_monthidx.(RANDOM_SA_FILTER_household_npc_ISOSLICE.month,RANDOM_SA_FILTER_household_npc_ISOSLICE.year, YEAR_START = YEAR_START)

        for i in 1:size(RANDOM_SA_FILTER_household_npc_ISOSLICE)[1]
            month = monthidxs[i]
            normalised_cITNPC_mean = RANDOM_SA_FILTER_household_npc_ISOSLICE.NPC_mean[i]*(RANDOM_SA_FILTER_household_npc_ISOSLICE.cITNPC_mean[i]/(RANDOM_SA_FILTER_household_npc_ISOSLICE.cITNPC_mean[i]+RANDOM_SA_FILTER_household_npc_ISOSLICE.LLINPC_mean[i]))
            normalised_LLINPC_mean = RANDOM_SA_FILTER_household_npc_ISOSLICE.NPC_mean[i]*(RANDOM_SA_FILTER_household_npc_ISOSLICE.LLINPC_mean[i]/(RANDOM_SA_FILTER_household_npc_ISOSLICE.cITNPC_mean[i]+RANDOM_SA_FILTER_household_npc_ISOSLICE.LLINPC_mean[i]))
            push!(RANDOM_SA_HOUSEHOLD_cITNPC_LISTS[n,month], RANDOM_SA_FILTER_household_npc_ISOSLICE.cITNPC_mean[i])
            push!(RANDOM_SA_HOUSEHOLD_LLINPC_LISTS[n,month], RANDOM_SA_FILTER_household_npc_ISOSLICE.LLINPC_mean[i])
            push!(RANDOM_SA_HOUSEHOLD_NPC_LISTS[n,month], RANDOM_SA_FILTER_household_npc_ISOSLICE.NPC_mean[i])
            push!(RANDOM_SA_HOUSEHOLD_NPC_STD_LISTS[n,month], RANDOM_SA_FILTER_household_npc_ISOSLICE.NPC_adj_se[i])
        end
    end

    # Scenario 3: 1-fold cross validation
    for n in 1:n_surveys
        # Filter out on a range of desired surveys
        incl_survey_names = setdiff(survey_names,[survey_names[n]])

        membership_bool = zeros(Bool, length(household_npc_ISOSLICE.SurveyId))
        for i in 1:length(membership_bool)
            if household_npc_ISOSLICE.SurveyId[i] ∈ incl_survey_names
                membership_bool[i] = 1
            end
        end

        CROSS_SA_FILTER_household_npc_ISOSLICE = household_npc_ISOSLICE[findall(membership_bool),:]

        monthidxs = monthyear_to_monthidx.(CROSS_SA_FILTER_household_npc_ISOSLICE.month,CROSS_SA_FILTER_household_npc_ISOSLICE.year, YEAR_START = YEAR_START)

        for i in 1:size(CROSS_SA_FILTER_household_npc_ISOSLICE)[1]
            month = monthidxs[i]
            normalised_cITNPC_mean = CROSS_SA_FILTER_household_npc_ISOSLICE.NPC_mean[i]*(CROSS_SA_FILTER_household_npc_ISOSLICE.cITNPC_mean[i]/(CROSS_SA_FILTER_household_npc_ISOSLICE.cITNPC_mean[i]+CROSS_SA_FILTER_household_npc_ISOSLICE.LLINPC_mean[i]))
            normalised_LLINPC_mean = CROSS_SA_FILTER_household_npc_ISOSLICE.NPC_mean[i]*(CROSS_SA_FILTER_household_npc_ISOSLICE.LLINPC_mean[i]/(CROSS_SA_FILTER_household_npc_ISOSLICE.cITNPC_mean[i]+CROSS_SA_FILTER_household_npc_ISOSLICE.LLINPC_mean[i]))
            push!(CROSS_SA_HOUSEHOLD_cITNPC_LISTS[n,month], CROSS_SA_FILTER_household_npc_ISOSLICE.cITNPC_mean[i])
            push!(CROSS_SA_HOUSEHOLD_LLINPC_LISTS[n,month], CROSS_SA_FILTER_household_npc_ISOSLICE.LLINPC_mean[i])
            push!(CROSS_SA_HOUSEHOLD_NPC_LISTS[n,month], CROSS_SA_FILTER_household_npc_ISOSLICE.NPC_mean[i])
            push!(CROSS_SA_HOUSEHOLD_NPC_STD_LISTS[n,month], CROSS_SA_FILTER_household_npc_ISOSLICE.NPC_adj_se[i])
        end
    end

    # Calculate NPC for each SA filter and then save
    # Chronological Inclusion
    for n in 1:n_surveys
        ########################################
        # %% Declare storage variables
        ########################################
        HOUSEHOLD_cITNPC_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        HOUSEHOLD_LLINPC_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        HOUSEHOLD_NPC_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        HOUSEHOLD_NPC_STD_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        cITN_CROP_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        LLIN_CROP_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        NET_CROP_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        NET_CROP_STD_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))

        ########################################
        # %% Average out list entries and save NPC
        ########################################
        
        valid_idx = findall(.!isempty.(CHRONO_SA_HOUSEHOLD_NPC_LISTS[n,:]))
        HOUSEHOLD_cITNPC_MONTHLY[valid_idx] = mean.(CHRONO_SA_HOUSEHOLD_cITNPC_LISTS[n,:][valid_idx])
        HOUSEHOLD_LLINPC_MONTHLY[valid_idx] = mean.(CHRONO_SA_HOUSEHOLD_LLINPC_LISTS[n,:][valid_idx])
        HOUSEHOLD_NPC_MONTHLY[valid_idx] = mean.(CHRONO_SA_HOUSEHOLD_NPC_LISTS[n,:][valid_idx])
        HOUSEHOLD_NPC_STD_MONTHLY[valid_idx] = mean.(CHRONO_SA_HOUSEHOLD_NPC_STD_LISTS[n,:][valid_idx])

        ########################################
        # %% Smooth out monthly NCP Data
        ########################################
        HOUSEHOLD_NPC_MONTHLY_SMOOTHED = missings(Float64, length(MONTHS_MONTHLY))
        for i in 1:length(MONTHS_MONTHLY)
            L_index = max(i-R_window,1)
            R_index = min(i,length(MONTHS_MONTHLY))
            window = HOUSEHOLD_NPC_MONTHLY[L_index:R_index]
            if !ismissing(HOUSEHOLD_NPC_MONTHLY[i])
                HOUSEHOLD_NPC_MONTHLY_SMOOTHED[i] = mean(window[findall(.!ismissing.(window))])
            end
        end

        ########################################
        # %% Estimate Monthly Net Crop based on population data NPC survey data
        ########################################
        # NET_CROP_MONTHLY_SMOOTHED = missings(Float64, length(MONTHS_MONTHLY))
        
        for i in 1:length(MONTHS_MONTHLY)
            if !ismissing(HOUSEHOLD_NPC_MONTHLY[i])
                cITN_CROP_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_cITNPC_MONTHLY[i]
                LLIN_CROP_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_LLINPC_MONTHLY[i]
                NET_CROP_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_NPC_MONTHLY[i]
                # NET_CROP_MONTHLY_SMOOTHED[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_NPC_MONTHLY_SMOOTHED[i]
                NET_CROP_STD_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_NPC_STD_MONTHLY[i]
            end
        end

        input_data_dict = Dict("ISO" => ISO,
                                "YEAR_START" => YEAR_START,
                                "YEAR_END" => YEAR_END,
                                "YEARS_ANNUAL" => YEARS_ANNUAL,
                                "DELIVERIES_ANNUAL" => DELIVERIES_ANNUAL,
                                "DISTRIBUTION_ANNUAL" => DISTRIBUTION_ANNUAL,
                                "POPULATION_ANNUAL" => POPULATION_ANNUAL,
                                "POPULATION_MONTHLY" => POPULATION_MONTHLY,
                                "MONTHS_MONTHLY" => MONTHS_MONTHLY,
                                "HOUSEHOLD_NPC_MONTHLY" => HOUSEHOLD_NPC_MONTHLY,
                                "HOUSEHOLD_NPC_STD_MONTHLY" => HOUSEHOLD_NPC_STD_MONTHLY,
                                "NET_CROP_MONTHLY" => NET_CROP_MONTHLY,
                                "NET_CROP_STD_MONTHLY" => NET_CROP_STD_MONTHLY,
                                "cITN_CROP_MONTHLY" => cITN_CROP_MONTHLY,
                                "LLIN_CROP_MONTHLY" => LLIN_CROP_MONTHLY,
                                # "NET_CROP_MONTHLY_SMOOTHED" => NET_CROP_MONTHLY_SMOOTHED,
                                "NET_NAMES" => NET_NAMES)

        ########################################
        # %% Save Extracted Data into JLD2 file
        ########################################

        # Make directory path to store data
        target_path = output_dir*"/crop/sensitivity_analysis/$(YEAR_START)_$(YEAR_END)/chronological/$(ISO)/$(n)/"
        mkpath(target_path)
        
        if save_output
            output_filename = (target_path)*(ISO*"_"*"$(YEAR_START)"*"_"*"$(YEAR_END)"*"_cropextract_CSA.jld2")
            save(output_filename, input_data_dict)
        end
    end

    # Random Inclusion
    for n in 1:n_surveys
        ########################################
        # %% Declare storage variables
        ########################################
        HOUSEHOLD_cITNPC_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        HOUSEHOLD_LLINPC_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        HOUSEHOLD_NPC_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        HOUSEHOLD_NPC_STD_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        cITN_CROP_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        LLIN_CROP_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        NET_CROP_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        NET_CROP_STD_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))

        ########################################
        # %% Average out list entries and save NPC
        ########################################
        
        valid_idx = findall(.!isempty.(RANDOM_SA_HOUSEHOLD_NPC_LISTS[n,:]))
        HOUSEHOLD_cITNPC_MONTHLY[valid_idx] = mean.(RANDOM_SA_HOUSEHOLD_cITNPC_LISTS[n,:][valid_idx])
        HOUSEHOLD_LLINPC_MONTHLY[valid_idx] = mean.(RANDOM_SA_HOUSEHOLD_LLINPC_LISTS[n,:][valid_idx])
        HOUSEHOLD_NPC_MONTHLY[valid_idx] = mean.(RANDOM_SA_HOUSEHOLD_NPC_LISTS[n,:][valid_idx])
        HOUSEHOLD_NPC_STD_MONTHLY[valid_idx] = mean.(RANDOM_SA_HOUSEHOLD_NPC_STD_LISTS[n,:][valid_idx])

        ########################################
        # %% Smooth out monthly NCP Data
        ########################################
        HOUSEHOLD_NPC_MONTHLY_SMOOTHED = missings(Float64, length(MONTHS_MONTHLY))
        for i in 1:length(MONTHS_MONTHLY)
            L_index = max(i-R_window,1)
            R_index = min(i,length(MONTHS_MONTHLY))
            window = HOUSEHOLD_NPC_MONTHLY[L_index:R_index]
            if !ismissing(HOUSEHOLD_NPC_MONTHLY[i])
                HOUSEHOLD_NPC_MONTHLY_SMOOTHED[i] = mean(window[findall(.!ismissing.(window))])
            end
        end

        ########################################
        # %% Estimate Monthly Net Crop based on population data NPC survey data
        ########################################
        # NET_CROP_MONTHLY_SMOOTHED = missings(Float64, length(MONTHS_MONTHLY))
        
        for i in 1:length(MONTHS_MONTHLY)
            if !ismissing(HOUSEHOLD_NPC_MONTHLY[i])
                cITN_CROP_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_cITNPC_MONTHLY[i]
                LLIN_CROP_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_LLINPC_MONTHLY[i]
                NET_CROP_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_NPC_MONTHLY[i]
                # NET_CROP_MONTHLY_SMOOTHED[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_NPC_MONTHLY_SMOOTHED[i]
                NET_CROP_STD_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_NPC_STD_MONTHLY[i]
            end
        end

        input_data_dict = Dict("ISO" => ISO,
                                "YEAR_START" => YEAR_START,
                                "YEAR_END" => YEAR_END,
                                "YEARS_ANNUAL" => YEARS_ANNUAL,
                                "DELIVERIES_ANNUAL" => DELIVERIES_ANNUAL,
                                "DISTRIBUTION_ANNUAL" => DISTRIBUTION_ANNUAL,
                                "POPULATION_ANNUAL" => POPULATION_ANNUAL,
                                "POPULATION_MONTHLY" => POPULATION_MONTHLY,
                                "MONTHS_MONTHLY" => MONTHS_MONTHLY,
                                "HOUSEHOLD_NPC_MONTHLY" => HOUSEHOLD_NPC_MONTHLY,
                                "HOUSEHOLD_NPC_STD_MONTHLY" => HOUSEHOLD_NPC_STD_MONTHLY,
                                "NET_CROP_MONTHLY" => NET_CROP_MONTHLY,
                                "NET_CROP_STD_MONTHLY" => NET_CROP_STD_MONTHLY,
                                "cITN_CROP_MONTHLY" => cITN_CROP_MONTHLY,
                                "LLIN_CROP_MONTHLY" => LLIN_CROP_MONTHLY,
                                # "NET_CROP_MONTHLY_SMOOTHED" => NET_CROP_MONTHLY_SMOOTHED,
                                "NET_NAMES" => NET_NAMES)

        ########################################
        # %% Save Extracted Data into JLD2 file
        ########################################

        # Make directory path to store data
        target_path = output_dir*"/crop/sensitivity_analysis/$(YEAR_START)_$(YEAR_END)/random/$(ISO)/$(n)/"
        mkpath(target_path)
        
        if save_output
            output_filename = (target_path)*(ISO*"_"*"$(YEAR_START)"*"_"*"$(YEAR_END)"*"_cropextract_RSA.jld2")
            save(output_filename, input_data_dict)
        end
    end

    # Cross Validation
    for n in 1:n_surveys
        ########################################
        # %% Declare storage variables
        ########################################
        HOUSEHOLD_cITNPC_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        HOUSEHOLD_LLINPC_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        HOUSEHOLD_NPC_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        HOUSEHOLD_NPC_STD_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        cITN_CROP_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        LLIN_CROP_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        NET_CROP_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
        NET_CROP_STD_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))

        ########################################
        # %% Average out list entries and save NPC
        ########################################
        
        valid_idx = findall(.!isempty.(CROSS_SA_HOUSEHOLD_NPC_LISTS[n,:]))
        HOUSEHOLD_cITNPC_MONTHLY[valid_idx] = mean.(CROSS_SA_HOUSEHOLD_cITNPC_LISTS[n,:][valid_idx])
        HOUSEHOLD_LLINPC_MONTHLY[valid_idx] = mean.(CROSS_SA_HOUSEHOLD_LLINPC_LISTS[n,:][valid_idx])
        HOUSEHOLD_NPC_MONTHLY[valid_idx] = mean.(CROSS_SA_HOUSEHOLD_NPC_LISTS[n,:][valid_idx])
        HOUSEHOLD_NPC_STD_MONTHLY[valid_idx] = mean.(CROSS_SA_HOUSEHOLD_NPC_STD_LISTS[n,:][valid_idx])

        ########################################
        # %% Smooth out monthly NCP Data
        ########################################
        HOUSEHOLD_NPC_MONTHLY_SMOOTHED = missings(Float64, length(MONTHS_MONTHLY))
        for i in 1:length(MONTHS_MONTHLY)
            L_index = max(i-R_window,1)
            R_index = min(i,length(MONTHS_MONTHLY))
            window = HOUSEHOLD_NPC_MONTHLY[L_index:R_index]
            if !ismissing(HOUSEHOLD_NPC_MONTHLY[i])
                HOUSEHOLD_NPC_MONTHLY_SMOOTHED[i] = mean(window[findall(.!ismissing.(window))])
            end
        end

        ########################################
        # %% Estimate Monthly Net Crop based on population data NPC survey data
        ########################################
        # NET_CROP_MONTHLY_SMOOTHED = missings(Float64, length(MONTHS_MONTHLY))
        
        for i in 1:length(MONTHS_MONTHLY)
            if !ismissing(HOUSEHOLD_NPC_MONTHLY[i])
                cITN_CROP_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_cITNPC_MONTHLY[i]
                LLIN_CROP_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_LLINPC_MONTHLY[i]
                NET_CROP_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_NPC_MONTHLY[i]
                # NET_CROP_MONTHLY_SMOOTHED[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_NPC_MONTHLY_SMOOTHED[i]
                NET_CROP_STD_MONTHLY[i] = POPULATION_MONTHLY[i]*HOUSEHOLD_NPC_STD_MONTHLY[i]
            end
        end

        input_data_dict = Dict("ISO" => ISO,
                                "YEAR_START" => YEAR_START,
                                "YEAR_END" => YEAR_END,
                                "YEARS_ANNUAL" => YEARS_ANNUAL,
                                "DELIVERIES_ANNUAL" => DELIVERIES_ANNUAL,
                                "DISTRIBUTION_ANNUAL" => DISTRIBUTION_ANNUAL,
                                "POPULATION_ANNUAL" => POPULATION_ANNUAL,
                                "POPULATION_MONTHLY" => POPULATION_MONTHLY,
                                "MONTHS_MONTHLY" => MONTHS_MONTHLY,
                                "HOUSEHOLD_NPC_MONTHLY" => HOUSEHOLD_NPC_MONTHLY,
                                "HOUSEHOLD_NPC_STD_MONTHLY" => HOUSEHOLD_NPC_STD_MONTHLY,
                                "NET_CROP_MONTHLY" => NET_CROP_MONTHLY,
                                "NET_CROP_STD_MONTHLY" => NET_CROP_STD_MONTHLY,
                                "cITN_CROP_MONTHLY" => cITN_CROP_MONTHLY,
                                "LLIN_CROP_MONTHLY" => LLIN_CROP_MONTHLY,
                                # "NET_CROP_MONTHLY_SMOOTHED" => NET_CROP_MONTHLY_SMOOTHED,
                                "NET_NAMES" => NET_NAMES)

        ########################################
        # %% Save Extracted Data into JLD2 file
        ########################################

        # Make directory path to store data
        target_path = output_dir*"/crop/sensitivity_analysis/$(YEAR_START)_$(YEAR_END)/cross/$(ISO)/$(n)/"
        mkpath(target_path)
        
        if save_output
            output_filename = (target_path)*(ISO*"_"*"$(YEAR_START)"*"_"*"$(YEAR_END)"*"_cropextract_CVA.jld2")
            save(output_filename, input_data_dict)
        end
    end
end

"""
    extract_data_netaccess(ISO::String, YEAR_START::Int64, YEAR_END::Int64)
    
Read from pre-declared directories of CSV files and extract required raw data to perform model analysis on. Similar as `extract_data()` but specifically outputs data needed for calculating net access.
- n_max: Maximum number of nets that are considered in aggregated data
- h_max: Maximum household size that is considered in aggregated data
"""
function extract_data_netaccess(ISO::String, YEAR_START::Int64, YEAR_END::Int64,
                input_dict, 
                regressions_dict;
                n_max = 20, 
                h_max = 10,
                # population_data_dir=population_data_dir,
                hh_data_dir=hh_dataprep_dir, 
                countrycodes_dir=countrycodes_dir,
                population_data_filename=population_data_filename,
                hh_listing_filename=hh_listing_filename,
                countrycodes_filename=countrycodes_filename,
                save_output = true, reg_mode = false)

    ########################################
    # %% Calculate time bounds for conversion from annual to monthly
    ########################################
    N_MONTHS = (YEAR_END-YEAR_START+1)*12

    ########################################
    # %% Import Country Code Conversion Table (for survey data)
    ########################################
    country_codes = CSV.read(countrycodes_dir*countrycodes_filename, DataFrame)
    DHS = country_codes.DHS_CODE[findfirst(country_codes.ISO3 .== ISO)]
    
    ########################################
    # %% Declare storage variables
    ########################################
    # Define year bounds
    YEARS_ANNUAL = Vector(YEAR_START:YEAR_END)
    # Needs to be extended by one year for interpolation puposes
    POPULATION_YEARS = Vector(YEAR_START:(YEAR_END+1))
    POPULATION_ANNUAL = zeros(Int64, length(POPULATION_YEARS))

    MONTHS_MONTHLY = Vector(1:N_MONTHS)
    # POPULATION_MONTHLY = input_dict["POPULATION_MONTHLY"]

    # Import Household Data
    
    hh_surveydata = CSV.read(hh_data_dir*hh_listing_filename,DataFrame)
    hh_surveydata = hh_surveydata[findall(.!ismissing.(hh_surveydata.ISO)),:]
    hh_surveydata_ISOSLICE = hh_surveydata[findall((hh_surveydata.ISO .== ISO).* 
                                            (hh_surveydata.interview_year .>= YEAR_START).*
                                            (hh_surveydata.interview_year .<= YEAR_END)),:]

    # household_survey_listing = CSV.read(hh_data_dir*hh_listing_filename, DataFrame)
    # household_survey_listing_ISOSLICE = household_survey_listing[findall(   
    #                                                                     (household_survey_listing.DHS_CountryCode.==DHS).*
    #                                                                     (household_survey_listing.IndicatorData.==1)
    #                                                                     ),:]
    
    # Get all survey years
    # year_vals = intersect(YEAR_START:YEAR_END,household_survey_listing_ISOSLICE.SurveyYear)
    hh_surveyids_ISOSLICE = unique(hh_surveydata_ISOSLICE.SurveyId)
    n_surveys = length(hh_surveyids_ISOSLICE)#length(year_vals)

    # NEW ADDITION! - TESTING TO ACCOUNT FOR CASES WHERE THERE ARE LIMITED SURVEYS AVAILABLE (Added 16/1/2025)
    if !reg_mode # If data is to be used for regression of net access model, then exclude countries with no surveys
        if n_surveys < 3 # If below threshold of 5 surveys, just use the entire dataset
            hh_surveydata_ISOSLICE = hh_surveydata[findall((hh_surveydata.interview_year .>= YEAR_START).*
                                                (hh_surveydata.interview_year .<= YEAR_END)),:]
            hh_surveyids_ISOSLICE = unique(hh_surveydata_ISOSLICE.SurveyId)
            n_surveys = length(hh_surveyids_ISOSLICE)#length(year_vals)
        end
    end

    # Aggregated values for net access model parameters
    ρ_h_aggregated = missings(Float64, n_surveys, h_max)
    μ_h_aggregated = missings(Float64, n_surveys, h_max)
    p_h_aggregated = missings(Float64, n_surveys, h_max)
    H_aggregated = missings(Float64, n_surveys,h_max,n_max)
    γ_aggregated = missings(Float64, n_surveys)
    monthidx_aggregated = missings(Int64, n_surveys)

    ########################################
    # %% Population Data
    ########################################

    POPULATION_MONTHLY = input_dict["POPULATION_MONTHLY"]

    ########################################
    # %% Get Γ_MONTHLY simulated net crop data for calculating net crop per capita from surveys
    ########################################

    # Need to first simulate Γ_MONTHLY from regressed net crop model
    YEARS_ANNUAL = input_dict["YEARS_ANNUAL"]
    MONTHS_MONTHLY = input_dict["MONTHS_MONTHLY"]
    DELIVERIES_ANNUAL = input_dict["DELIVERIES_ANNUAL"]
    DISTRIBUTION_ANNUAL = input_dict["DISTRIBUTION_ANNUAL"]
    monthly_p = regressions_dict["monthly_p"]
    chain = regressions_dict["chain"]
    NET_NAMES = input_dict["NET_NAMES"]
    
    # Import number of nets
    n_net_types = length(NET_NAMES)

    # Extract posterior draws
    ϕ_est = mean(chain[:, :ϕ])
    α_init_est = mean(chain[:,:α_init])
    α_LLIN_est = mean(chain[:,:α_LLIN])
    b_net_est = mean(Matrix(DataFrame(chain)[:,5:2:5+2*(n_net_types-1)]), dims = 1)[:]
    k_net_est = mean(Matrix(DataFrame(chain)[:,6:2:6+2*(n_net_types-1)]), dims = 1)[:]

    n_missing_nets_vals = sum(ismissing.(DISTRIBUTION_ANNUAL[:,1]))
    if n_missing_nets_vals > 0
        missing_nets_vals = Matrix(DataFrame(chain)[:,((4+2*(n_net_types))+1):end])
        missing_nets_est = mean(missing_nets_vals, dims = 1)[:]
    else
        missing_nets_est = zeros(size(chain)[1],0)
    end

    Γ_epochs_BYNET = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                        DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                        ϕ_est, b_net_est, k_net_est,
                                        α_init_est, α_LLIN_est,
                                        missing_nets_est; 
                                        monthly_p = monthly_p)

    Γ_MONTHLY = sum(Γ_epochs_BYNET, dims = 2)[:]

    

    ########################################
    # %% Survey Data
    ########################################
    skipped_surveys_idx = [] # Need to track which survey files were problematic and need to remove from analysis
    # Select a year
    for i in 1:n_surveys

        surveyid = hh_surveyids_ISOSLICE[i]

        # year = year_vals[i]
        # row_idx = findfirst(household_survey_listing_ISOSLICE.SurveyYear .== year)

        # # Get the raw survey filename
        # survey_num = household_survey_listing_ISOSLICE.SurveyNum[row_idx]
        # survey_filename = "\\standard_tables\\Svy_$(survey_num)_ITN_HH_Res.csv"

        # # Import the survey
        # hh_data_RAW = nothing
        # try # Try to import result
        #     hh_data_RAW = CSV.read(hh_data_dir*survey_filename, DataFrame)
        # catch # If can't find file
        #     println("Can't find file: $(survey_filename). Skipping to next survey")
        #     push!(skipped_surveys_idx,i)
        #     continue
        # end

        current_hh_survey = hh_surveydata_ISOSLICE[findall(hh_surveydata_ISOSLICE.SurveyId .== surveyid),:]

        hh_sizes = current_hh_survey.hh_size
        n_itns = current_hh_survey.n_itn
        hh_weight = current_hh_survey.hh_sample_wt
        # # Get interview years
        # interview_years = unique(hh_data_RAW.interview_year)
        # interview_years = intersect(interview_years[findall(.!ismissing.(interview_years))], YEARS_ANNUAL)

        # # Need to check and ignore potential missing data
        # nonmissing_idxs = findall((.!ismissing.(hh_data_RAW.hh_size)).*
        #                             (.!ismissing.(hh_data_RAW.n_itn)))
        # nonmissing_hh_size = hh_data_RAW.hh_size[nonmissing_idxs]
        # nonmissing_n_itn = hh_data_RAW.n_itn[nonmissing_idxs]
        

        # %% Calculate H matrix
        # Matrix of weighted household counts
        HH_matrix = zeros(h_max,n_max) 

        for n in 0:size(HH_matrix)[2]-1
            for h in 1:size(HH_matrix)[1]
                HH_matrix[h,n+1] = sum((hh_sizes .== h).*
                                                (n_itns .== n) .* hh_weight)
            end
        end

        ρ_h = missings(Float64, h_max)
        μ_h = missings(Float64, h_max)

        for h in 1:length(ρ_h)

            # Calculate proportion of households of size h with no nets
            num_nonets_hh = sum(HH_matrix[h,:])
            if num_nonets_hh != 0
                ρ_h[h] = HH_matrix[h,1]/sum(HH_matrix[h,:])
            end

            # Get idxs for households of size h with at least one net
            idxs = findall((hh_sizes .== h).*(n_itns .>= 1))
            if !isempty(idxs)
                μ_h[h] = mean(n_itns[idxs])
            end
        end

        # Proportion of households with size h
        p_h = sum(HH_matrix, dims = 2)[:]./sum(HH_matrix)

        # Calculate H matrix for survey
        H = zeros(h_max,n_max)

        for h in 1:size(HH_matrix)[1]
            if (!ismissing(ρ_h[h])) && (!ismissing(μ_h[h]))
                for n in 0:size(HH_matrix)[2]-1
                    if n == 0
                        H[h,n+1] = p_h[h]*ρ_h[h]
                    else
                        H[h,n+1] = p_h[h]*(1-ρ_h[h])*((μ_h[h]^n)/(factorial(big(n))*(exp(μ_h[h])-1)))
                    end
                end
            end
        end

        # Calculate net crop per capita from model posterior draws for net crop
        monthidx_vals = []#zeros(Int,size(hh_data_RAW)[1])
        
        for j in 1:size(current_hh_survey)[1]
            idx_val = nothing 
            
            try # Try to see if can compute monthidx
                idx_val = monthyear_to_monthidx(current_hh_survey.interview_month[j],
                                                    current_hh_survey.interview_year[j], YEAR_START = YEAR_START)
            catch
                # There was an error. Skip this entry
                println("Skipping survey num $(survey_num), entry $(i)...")
                continue
            end
            if !ismissing(idx_val)
                push!(monthidx_vals, idx_val)
            end
        end

        γ = mean(Γ_MONTHLY[monthidx_vals]./POPULATION_MONTHLY[monthidx_vals])

        # Store results in aggregate variable
        ρ_h_aggregated[i,:] = copy(ρ_h)
        μ_h_aggregated[i,:] = copy(μ_h)
        p_h_aggregated[i,:] = copy(p_h)
        H_aggregated[i,:,:] = copy(H)
        γ_aggregated[i] = copy(γ)

        # Need to get monthidx to align with coincide with existing survey entry
        # for algorithm running purposes. NOTE: THIS IS ONLY A PROBLEM FOR MICS datasets
        average_monthidx = round(Int64, mean(monthidx_vals))
        if average_monthidx ∈ monthidx_vals
            # There is at least one entry coinciding, therefore will not cause issues in
            # running algorithm. 
            # Do saving as per normal.
            monthidx_aggregated[i] = round(Int64, mean(monthidx_vals))
        else
            # Average monthidx doesn't align with any of the survey monthidxs
            # This is the case if a given survey contains entries for months with large gaps 
            # in between. Use the closest entry instead.
            new_target_monthidx = monthidx_vals[argmin(abs.(monthidx_vals.-average_monthidx))]
            monthidx_aggregated[i] = new_target_monthidx
        end
    end

    # # Remove rows in aggregated matrix that contains problematic surveys
    # idxs = setdiff(1:n_surveys, skipped_surveys_idx)
    # ρ_h_aggregated = ρ_h_aggregated[idxs,:]
    # μ_h_aggregated = μ_h_aggregated[idxs,:]
    # p_h_aggregated = p_h_aggregated[idxs,:]
    # H_aggregated = H_aggregated[idxs,:,:]
    # γ_aggregated = γ_aggregated[idxs]

    # Store results as a dict
    net_access_input_dict = Dict("ISO" => ISO,
                            "YEAR_START" => YEAR_START,
                            "YEAR_END" => YEAR_END,
                            "YEARS_ANNUAL" => YEARS_ANNUAL,
                            "POPULATION_ANNUAL" => POPULATION_ANNUAL,
                            "POPULATION_MONTHLY" => POPULATION_MONTHLY,
                            "ρ_h_aggregated" => ρ_h_aggregated,
                            "μ_h_aggregated" => μ_h_aggregated,
                            "p_h_aggregated" => p_h_aggregated,
                            "H_aggregated" => H_aggregated,
                            "γ_aggregated" => γ_aggregated,
                            "survey_monthidx_aggregated" => monthidx_aggregated,
                            "n_max" => n_max,
                            "h_max" => h_max)

    
    ########################################
    # %% Save Extracted Data into JLD2 file
    ########################################
    target_path = ""

    if reg_mode
        target_path = output_dir*"/access/reg_data/$(YEAR_START)_$(YEAR_END)/"
        mkpath(target_path)
    else
        target_path = output_dir*"/access/pred_data/$(YEAR_START)_$(YEAR_END)/"
        mkpath(target_path)
    end
    
    if save_output
        output_filename = (target_path)*ISO*"_"*"$(YEAR_START)"*"_"*"$(YEAR_END)"*"_accessextract.jld2"
        save(output_filename, net_access_input_dict)
    end
    return net_access_input_dict
end

end