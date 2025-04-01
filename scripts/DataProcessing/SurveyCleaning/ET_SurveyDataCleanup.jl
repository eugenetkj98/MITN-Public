"""
Author: Eugene Tan
Date Created: 22/8/2024
Last Updated: 22/8/2024
This is a script to consolidate the DHS and MICS4-6 survey datasets into a single list.
Outputs 2 csv files:
(1) survey_listing.csv - Contains a metadata pertaining to each survey IDs
(2) itn_hh_surveydata_complete.csv - Compilation of all survey entries for all countries
"""

# %% Declare required paths
push!(LOAD_PATH,"src/") # Path for source files and modules

# %% Activate required environment
# Package manager (This should install the required packages. If not, run the "instantiate" command via Pkg)
include(pwd()*"/scripts/init_env.jl")

# %% Data Wrangling
using CSV
using DataFrames

# %% Define Directories
datasets_dir = "Z:/eugene/ITN Input Datasets/"

###########################################
# %% Get List of Survey IDs and Meta Data
###########################################

# %% Country codes
country_codes_key = CSV.read(datasets_dir*"country_codes.csv", DataFrame)

# %% DHS Data
dhs_survey_filename = "dhs_survey_listing.csv"
dhs_legend = CSV.read(datasets_dir*dhs_survey_filename, DataFrame)

# Filter rows to countries of interest
ind_completeness = (dhs_legend.IndicatorData.==1)
ind_dhsmembership = zeros(Bool, length(ind_completeness))

for i in 1:length(ind_dhsmembership)
    if dhs_legend.DHS_CountryCode[i] ∈ country_codes_key.DHS_CODE
        ind_dhsmembership[i] = 1
    end
end

dhs_filt_row_idxs = findall(ind_completeness.*ind_dhsmembership)
dhs_filt_legend = dhs_legend[dhs_filt_row_idxs,["SurveyId", "DHS_CountryCode", "SurveyYear","SurveyType","CountryName"]]

# Get ISO codes to append to legend
DHS_codes = dhs_filt_legend.DHS_CountryCode
ISO_codes = Array{String}(undef, length(DHS_codes))

for i in 1:length(DHS_codes)
    DHS_id = DHS_codes[i]
    ISO_id = country_codes_key.ISO3[findfirst(country_codes_key.DHS_CODE .== DHS_id)]
    ISO_codes[i] = ISO_id
end
dhs_filt_legend[:, :ISO] = ISO_codes

dhs_filt_legend

# %% Other HH Data
# File names

survey_filename = "other_hh.csv"
other_key_filename = "KEY_080817.csv"



# Import files
other_survey = CSV.read(datasets_dir*survey_filename, DataFrame)
other_legend = CSV.read(datasets_dir*other_key_filename, DataFrame)

# Get number of unique surveys
other_surveyids = unique(other_survey[:,"Survey.hh"])

#
other_dhscode = Array{String}(undef, size(other_surveyids)[1])
other_surveyyear = Array{Int}(undef, size(other_surveyids)[1])
other_surveytype = Array{String}(undef, size(other_surveyids)[1])
other_countrynames = Array{String}(undef, size(other_surveyids)[1])
other_countryisos = Array{String}(undef, size(other_surveyids)[1])

for i in 1:size(other_surveyids)[1]
    # Extract Metadata
    survey_id, country_name = other_legend[findfirst(other_legend[:,"Svy Name"] .== other_surveyids[i]),["Svy Name", "Name"]]
    surveyyear = other_survey.year[findfirst(other_survey[:,"Survey.hh"].== other_surveyids[1])]
    surveytype = "OTHER"
    dhscode = "N/A"
    country_iso = country_codes_key.ISO3[findfirst(country_codes_key.Country .== country_name)]

    # Save to DataFrame
    other_surveyids[i] = survey_id
    other_dhscode[i] = dhscode
    other_surveyyear[i] = surveyyear
    other_surveytype[i] = surveytype
    other_countrynames[i] = country_name
    other_countryisos[i] = country_iso
end

other_filt_legend = DataFrame("SurveyId" => other_surveyids, 
                            "DHS_CountryCode" => other_dhscode, 
                            "SurveyYear" => other_surveyyear,
                            "SurveyType" => other_surveytype,
                            "CountryName" => other_countrynames,
                            "ISO" => other_countryisos)

# %% MICS4 Data
# File names

mics4_filename = "mics4_hh_21_january.csv"
mics4_key_filename = "KEY_080817.csv"

# Import files
mics4_survey = CSV.read(datasets_dir*mics4_filename, DataFrame)
mics4_legend = CSV.read(datasets_dir*mics4_key_filename, DataFrame)

# Get number of unique surveys
mics4_surveyids = unique(mics4_survey[:,"Survey.hh"])

#
mics4_dhscode = Array{String}(undef, size(mics4_surveyids)[1])
mics4_surveyyear = Array{Int}(undef, size(mics4_surveyids)[1])
mics4_surveytype = Array{String}(undef, size(mics4_surveyids)[1])
mics4_countrynames = Array{String}(undef, size(mics4_surveyids)[1])
mics4_countryisos = Array{String}(undef, size(mics4_surveyids)[1])


for i in 1:size(mics4_surveyids)[1]
    # Extract Metadata
    survey_id, country_name = mics4_legend[findfirst(mics4_legend[:,"Svy Name"] .== mics4_surveyids[i]),["Svy Name", "Name"]]
    surveyyear = mics4_survey.year[findfirst(mics4_survey[:,"Survey.hh"].== mics4_surveyids[1])]
    surveytype = "MICS4"
    dhscode = "N/A"
    country_iso = country_codes_key.ISO3[findfirst(country_codes_key.Country .== country_name)]

    # Save to DataFrame
    mics4_surveyids[i] = survey_id
    mics4_dhscode[i] = dhscode
    mics4_surveyyear[i] = surveyyear
    mics4_surveytype[i] = surveytype
    mics4_countrynames[i] = country_name
    mics4_countryisos[i] = country_iso
end

mics4_filt_legend = DataFrame("SurveyId" => mics4_surveyids, 
                            "DHS_CountryCode" => mics4_dhscode, 
                            "SurveyYear" => mics4_surveyyear,
                            "SurveyType" => mics4_surveytype,
                            "CountryName" => mics4_countrynames,
                            "ISO" => mics4_countryisos)

# %% MICS5 Data

# File names

mics5_filename = "mics5_hh_18_june_2020.csv"

# Import files
mics5_survey = CSV.read(datasets_dir*mics5_filename, DataFrame)

# 

# Get number of unique surveys
mics5_surveyids = unique(mics5_survey.surveyid)

#
mics5_dhscode = Array{String}(undef, size(mics5_surveyids)[1])
mics5_surveyyear = Array{Int}(undef, size(mics5_surveyids)[1])
mics5_surveytype = Array{String}(undef, size(mics5_surveyids)[1])
mics5_countrynames = Array{String}(undef, size(mics5_surveyids)[1])
mics5_countryisos = Array{String}(undef, size(mics5_surveyids)[1])

for i in 1:size(mics5_surveyids)[1]
    # Extract Metadata
    rowidx = findfirst(mics5_survey.surveyid .== mics5_surveyids[i])
    surveyid = mics5_surveyids[i]
    country_name = mics5_survey.country[rowidx]
    surveyyear = mics5_survey.year[rowidx]
    surveytype = "MICS5"
    dhscode = "N/A"
    country_iso = mics5_survey.iso3[rowidx]

    # Save to DataFrame
    mics5_surveyids[i] = surveyid
    mics5_dhscode[i] = dhscode
    mics5_surveyyear[i] = surveyyear
    mics5_surveytype[i] = surveytype
    mics5_countrynames[i] = country_name
    mics5_countryisos[i] = country_iso
end

mics5_filt_legend = DataFrame("SurveyId" => mics5_surveyids, 
                            "DHS_CountryCode" => mics5_dhscode, 
                            "SurveyYear" => mics5_surveyyear,
                            "SurveyType" => mics5_surveytype,
                            "CountryName" => mics5_countrynames,
                            "ISO" => mics5_countryisos)

# %% MICS6 Data

# File names
mics6_filename = "mics6_hh_23_june_2020.csv"

# Import files
mics6_survey = CSV.read(datasets_dir*mics6_filename, DataFrame)

# Get number of unique surveys
mics6_surveyids = unique(mics6_survey.surveyid)

#
mics6_dhscode = Array{String}(undef, size(mics6_surveyids)[1])
mics6_surveyyear = Array{Int}(undef, size(mics6_surveyids)[1])
mics6_surveytype = Array{String}(undef, size(mics6_surveyids)[1])
mics6_countrynames = Array{String}(undef, size(mics6_surveyids)[1])
mics6_countryisos = Array{String}(undef, size(mics6_surveyids)[1])

for i in 1:size(mics6_surveyids)[1]
    # Extract Metadata
    rowidx = findfirst(mics6_survey.surveyid .== mics6_surveyids[i])
    surveyid = mics6_surveyids[i]
    country_name = mics6_survey.country[rowidx]
    surveyyear = mics6_survey.year[rowidx]
    surveytype = "MICS6"
    dhscode = "N/A"
    country_iso = mics6_survey.iso3[rowidx]

    # Save to DataFrame
    mics6_surveyids[i] = surveyid
    mics6_dhscode[i] = dhscode
    mics6_surveyyear[i] = surveyyear
    mics6_surveytype[i] = surveytype
    mics6_countrynames[i] = country_name
    mics6_countryisos[i] = country_iso
end

mics6_filt_legend = DataFrame("SurveyId" => mics6_surveyids, 
                            "DHS_CountryCode" => mics6_dhscode, 
                            "SurveyYear" => mics6_surveyyear,
                            "SurveyType" => mics6_surveytype,
                            "CountryName" => mics6_countrynames,
                            "ISO" => mics6_countryisos)

# %% Compile Legends together into single dataframe
master_survey_key = vcat(dhs_filt_legend,other_filt_legend,mics4_filt_legend,
                        mics5_filt_legend,mics6_filt_legend)
CSV.write(datasets_dir*"survey_listing.csv", master_survey_key)

###########################################
# %% Compile Master List of All Survey Entries
###########################################
master_survey_key = CSV.read(datasets_dir*"survey_listing.csv", DataFrame)

dhs_survey_filename = "dhs_survey_listing.csv"
dhs_legend = CSV.read(datasets_dir*dhs_survey_filename, DataFrame)
############
# %% DHS
############



# Function to scrape surveys
function scrape_DHS_survey(survey_data, dhs_legend, country_codes_key)
    # Get survey id to lookup in dhs legend
    survey_num = survey_data.surveyid[1]

    # Get number of survey entries
    n_entries = size(survey_data)[1]

    # Metadata for survey entries
    surveyid = String(dhs_legend.SurveyId[findfirst(dhs_legend.SurveyNum .== survey_num)])
    surveycountry = dhs_legend.CountryName[findfirst(dhs_legend.SurveyNum .== survey_num)]
    surveyiso = country_codes_key.ISO3[findfirst(country_codes_key.Country .== surveycountry)]

    # Scrape data
    surveyid_entries = Array{String}(undef,n_entries)
    surveycountry_entries = Array{String}(undef,n_entries)
    surveyiso_entries = Array{String}(undef,n_entries)
    surveytype_entries = Array{String}(undef,n_entries)
    surveylatitude_entries = Array{Union{Missing, Float64}}(undef,n_entries)
    surveylongitude_entries = Array{Union{Missing, Float64}}(undef,n_entries)
    n_citn_entries = Array{Union{Missing, Int64}}(undef,n_entries)
    n_llin_entries = Array{Union{Missing, Int64}}(undef,n_entries)
    surveyid_entries .= surveyid
    surveycountry_entries .= surveycountry
    surveyiso_entries .= surveyiso
    surveytype_entries .= "DHS"

    survey_clusterid_entries = survey_data[:,"clusterid"]
    for i in 1:n_entries
        lat = survey_data[i,"latitude"]
        long = survey_data[i,"longitude"] 
        if !ismissing(lat)
            dattype = typeof(lat)
            if (dattype == Float64)||(dattype == Int)
                surveylatitude_entries[i] = lat
                surveylongitude_entries[i] = long
            end
        end

        # Data Format Corrrection
        if typeof(survey_data.n_conv_itn[i]) == Int64
            n_citn_entries[i] = survey_data.n_conv_itn[i]
            n_llin_entries[i] = survey_data.n_llin[i]
        else
            n_citn_entries[i] = 0
            n_llin_entries[i] = survey_data.n_itn[i]
        end
    end

    surveymonth_entries = survey_data.interview_month
    surveyyear_entries = survey_data.interview_year
    hh_wt_entries = survey_data.hh_sample_wt
    hh_size_entries = survey_data.hh_size
    n_itn_entries = survey_data.n_itn
    n_itn_used_entries = survey_data.n_itn_used
    n_slept_under_itn_entries = survey_data.n_slept_under_itn
    
    # Concatenate data together

    survey_entries =    DataFrame("SurveyId" => surveyid_entries,
                                        "Country" => surveycountry_entries,
                                        "ISO" => surveyiso_entries,
                                        "SurveyType" => surveytype_entries,
                                        "clusterid" => survey_clusterid_entries,
                                        "latitude" => surveylatitude_entries,
                                        "longitude" => surveylongitude_entries,
                                        "interview_month" => surveymonth_entries,
                                        "interview_year" => surveyyear_entries,
                                        "hh_sample_wt" => hh_wt_entries,
                                        "hh_size" => hh_size_entries,
                                        "n_itn" => n_itn_entries,
                                        "n_itn_used" => n_itn_used_entries,
                                        "n_slept_under_itn" => n_slept_under_itn_entries,
                                        "n_citn" =>  n_citn_entries,
                                        "n_llin" => n_llin_entries)

    return survey_entries
end

# Define directory for DHS surveys and get filenames
dhs_survey_dir = datasets_dir*"standard_tables/"
dhs_survey_filenames = readdir(dhs_survey_dir)

# For each DHS survey, extract data
survey_dataframes = []

using ProgressBars
for i in ProgressBar(1:length(dhs_survey_filenames))
    survey_filename = dhs_survey_filenames[i]
    survey_data = CSV.read(dhs_survey_dir*survey_filename, DataFrame)
    survey_dataframe = scrape_DHS_survey(survey_data, dhs_legend, country_codes_key)
    push!(survey_dataframes, survey_dataframe)
end

dhs_compiled_survey_entries = vcat(survey_dataframes...)

############
# %% UNKNOWN Origin Other HH_data
############
other_filename = "other_hh.csv"
other_survey = CSV.read(datasets_dir*other_filename, DataFrame)

other_n_entries = size(other_survey)[1]
other_country_entries = Array{String}(undef, other_n_entries)
other_iso_entries = Array{String}(undef, other_n_entries)
other_surveytype_entries = Array{String}(undef, other_n_entries)
other_lat_entries = Array{Union{Missing,Float64}}(missing, other_n_entries)
other_long_entries = Array{Union{Missing,Float64}}(missing, other_n_entries)
other_n_itn_used_entries = Array{Union{Missing,Int}}(missing, other_n_entries)
other_n_slept_under_itn_entries = Array{Union{Missing,Int}}(missing, other_n_entries)

# %%

for i in 1:other_n_entries
    legend_lookup_idx = findfirst(master_survey_key.SurveyId .== other_survey[i,"Survey.hh"])
    other_country_entries[i] = master_survey_key[legend_lookup_idx,"CountryName"]
    other_iso_entries[i] = master_survey_key[legend_lookup_idx,"ISO"]
    other_surveytype_entries[i] = "OTHER"

    try
        other_n_itn_used_entries[i] = parse(Int64,other_survey[i,"n.ITN.used"])
        other_n_slept_under_itn_entries[i] = parse(Int64,other_survey[i,"n.individuals.that.slept.under.ITN"])
    catch
        continue
    end 
end
surveyid_entries = other_survey[:, "Survey.hh"]



other_compiled_survey_entries = DataFrame("SurveyId" => surveyid_entries,
                                            "Country" => other_country_entries,
                                            "ISO" => other_iso_entries,
                                            "SurveyType" => other_surveytype_entries,
                                            "clusterid" => other_survey[:,"Cluster.hh"],
                                            "latitude" => other_lat_entries,
                                            "longitude" => other_long_entries,
                                            "interview_month" => other_survey[:,"month"],
                                            "interview_year" => other_survey[:,"year"],
                                            "hh_sample_wt" => other_survey[:,"sample.w"],
                                            "hh_size" => other_survey[:,"hh.size"],
                                            "n_itn" => other_survey[:,"n.ITN.per.hh"],
                                            "n_itn_used" => other_n_itn_used_entries,
                                            "n_slept_under_itn" => other_n_slept_under_itn_entries,
                                            "n_citn" => other_survey[:,"n.conventional.ITNs"],
                                            "n_llin" => other_survey[:,"n.LLINs"])


############
# %% MICS4
############
mics4_survey = CSV.read(datasets_dir*mics4_filename, DataFrame)

mics4_n_entries = size(mics4_survey)[1]
mics4_country_entries = Array{String}(undef, mics4_n_entries)
mics4_iso_entries = Array{String}(undef, mics4_n_entries)
mics4_surveytype_entries = Array{String}(undef, mics4_n_entries)
mics4_lat_entries = Array{Union{Missing,Float64}}(missing, mics4_n_entries)
mics4_long_entries = Array{Union{Missing,Float64}}(missing, mics4_n_entries)

for i in 1:mics4_n_entries
    legend_lookup_idx = findfirst(master_survey_key.SurveyId .== mics4_survey[i,"Survey.hh"])
    mics4_country_entries[i] = master_survey_key[legend_lookup_idx,"CountryName"]
    mics4_iso_entries[i] = master_survey_key[legend_lookup_idx,"ISO"]
    mics4_surveytype_entries[i] = "MICS4"
end
surveyid_entries = mics4_survey[:, "Survey.hh"]

mics4_compiled_survey_entries = DataFrame("SurveyId" => surveyid_entries,
                                            "Country" => mics4_country_entries,
                                            "ISO" => mics4_iso_entries,
                                            "SurveyType" => mics4_surveytype_entries,
                                            "clusterid" => mics4_survey[:,"Cluster.hh"],
                                            "latitude" => mics4_lat_entries,
                                            "longitude" => mics4_long_entries,
                                            "interview_month" => mics4_survey[:,"month"],
                                            "interview_year" => mics4_survey[:,"year"],
                                            "hh_sample_wt" => mics4_survey[:,"sample.w"],
                                            "hh_size" => mics4_survey[:,"hh.size"],
                                            "n_itn" => mics4_survey[:,"n.ITN.per.hh"],
                                            "n_itn_used" => mics4_survey[:,"n.ITN.used"],
                                            "n_slept_under_itn" => mics4_survey[:,"n.individuals.that.slept.under.ITN"],
                                            "n_citn" => mics4_survey[:,"n.conventional.ITNs"],
                                            "n_llin" => mics4_survey[:,"n.LLIN"])

############
# %% MICS5
############

mics5_survey = CSV.read(datasets_dir*mics5_filename, DataFrame)
mics5_n_entries = size(mics5_survey)[1]

mics5_surveytype_entries = Array{String}(undef, mics5_n_entries)
mics5_surveytype_entries .= "MICS5"
mics5_lat_entries = Array{Union{Missing,Float64}}(missing, mics5_n_entries)
mics5_long_entries = Array{Union{Missing,Float64}}(missing, mics5_n_entries)
mics5_month_entries = Array{Union{Missing,Int}}(missing, mics5_n_entries)
mics5_n_itn_entries = Array{Union{Missing,Int}}(missing, mics5_n_entries)
mics5_n_itn_used_entries = Array{Union{Missing,Int}}(missing, mics5_n_entries)
mics5_n_slept_under_itn_entries = Array{Union{Missing,Int}}(missing, mics5_n_entries)
mics5_n_citn_entries = Array{Union{Missing,Int}}(missing, mics5_n_entries)
mics5_n_llin_entries = Array{Union{Missing,Int}}(missing, mics5_n_entries)

# Correct data formatting of entries

for i in 1:mics5_n_entries
    try
        mics5_n_itn_entries[i] = parse(Int64,mics5_survey.n_itn[i])
        mics5_n_itn_used_entries[i] = parse(Int64,mics5_survey.n_itn_used[i])
        mics5_n_slept_under_itn_entries[i] = parse(Int64,mics5_survey.n_slept_under_itn[i])
        mics5_month_entries[i] = mics5_survey.month[i]
        
        mics5_n_llin_entries[i] = parse(Int64,mics5_survey.n_llin[i])
    catch
        continue
    end 
    mics5_n_citn_entries[i] = mics5_survey.n_conv_itn[i]
end

mics5_compiled_survey_entries = DataFrame("SurveyId" => mics5_survey.surveyid,
                                            "Country" => mics5_survey.country,
                                            "ISO" => mics5_survey.iso3,
                                            "SurveyType" => mics5_surveytype_entries,
                                            "clusterid" => mics5_survey.clusterid,
                                            "latitude" => mics5_lat_entries,
                                            "longitude" => mics5_long_entries,
                                            "interview_month" => mics5_month_entries,
                                            "interview_year" => mics5_survey.year,
                                            "hh_sample_wt" => mics5_survey.hh_sample_wt,
                                            "hh_size" => mics5_survey.hh_size,
                                            "n_itn" => mics5_n_itn_entries,
                                            "n_itn_used" => mics5_n_itn_used_entries,
                                            "n_slept_under_itn" => mics5_n_slept_under_itn_entries,
                                            "n_citn" => mics5_n_citn_entries,
                                            "n_llin" => mics5_n_llin_entries)


############
# %% MICS6
############

mics6_survey = CSV.read(datasets_dir*mics6_filename, DataFrame)
mics6_n_entries = size(mics6_survey)[1]

mics6_surveytype_entries = Array{String}(undef, mics6_n_entries)
mics6_surveytype_entries .= "MICS6"
mics6_lat_entries = Array{Union{Missing,Float64}}(missing, mics6_n_entries)
mics6_long_entries = Array{Union{Missing,Float64}}(missing, mics6_n_entries)
mics6_month_entries = Array{Union{Missing,Int}}(missing, mics6_n_entries)
mics6_n_itn_entries = Array{Union{Missing,Int}}(missing, mics6_n_entries)
mics6_n_slept_under_itn_entries = Array{Union{Missing,Int}}(missing, mics6_n_entries)
mics6_n_itn_used_entries = Array{Union{Missing,Int}}(missing, mics6_n_entries)
mics6_n_citn_entries = Array{Union{Missing,Int}}(missing, mics6_n_entries)
mics6_n_llin_entries = Array{Union{Missing,Int}}(missing, mics6_n_entries)

# Correct data formatting of entries
for i in 1:mics6_n_entries
    try
        mics6_month_entries[i] = Int(parse(Float64,mics6_survey.month[i]))
        
    catch
        continue
    end
    mics6_n_itn_entries[i] = mics6_survey.n_itn[i]
    mics6_n_slept_under_itn_entries[i] = mics6_survey.n_slept_under_itn[i]
    mics6_n_itn_used_entries[i] = mics6_survey.n_itn_used[i]
    mics6_n_citn_entries[i] = mics6_survey.n_conv_itn[i]
    mics6_n_llin_entries[i] = mics6_survey.n_llin[i]
end

mics6_compiled_survey_entries = DataFrame("SurveyId" => mics6_survey.surveyid,
                                            "Country" => mics6_survey.country,
                                            "ISO" => mics6_survey.iso3,
                                            "SurveyType" => mics6_surveytype_entries,
                                            "clusterid" => mics6_survey.clusterid,
                                            "latitude" => mics6_lat_entries,
                                            "longitude" => mics6_long_entries,
                                            "interview_month" => mics6_month_entries,
                                            "interview_year" => mics6_survey.year,
                                            "hh_sample_wt" => mics6_survey.hh_sample_wt,
                                            "hh_size" => mics6_survey.hh_size,
                                            "n_itn" => mics6_n_itn_entries,
                                            "n_itn_used" => mics6_n_itn_used_entries,
                                            "n_slept_under_itn" => mics6_n_slept_under_itn_entries,
                                            "n_citn" => mics6_n_citn_entries,
                                            "n_llin" => mics6_n_llin_entries)

############
# %% Compile into single dataset
############

compiled_surveys = vcat(dhs_compiled_survey_entries, 
                        other_compiled_survey_entries,
                        mics4_compiled_survey_entries,
                        mics5_compiled_survey_entries,
                        mics6_compiled_survey_entries)

# %% Identify all data with missing
month_check = .!ismissing.(compiled_surveys.interview_month)
year_check = .!ismissing.(compiled_surveys.interview_year)
weight_check = .!ismissing.(compiled_surveys.hh_sample_wt)
hh_check = .!(ismissing.(compiled_surveys.hh_size))
itn_check = .!ismissing.(compiled_surveys.n_itn)
use_check = .!ismissing.(compiled_surveys.n_itn_used)
slept_check = .!ismissing.(compiled_surveys.n_slept_under_itn)

# Check for missing or corrupted data
valid_entry_idxs = intersect(1:size(compiled_surveys)[1], findall(month_check.*year_check.*weight_check.*hh_check.*itn_check.*use_check.*slept_check))
compiled_surveys = compiled_surveys[valid_entry_idxs,:]

# Filter our zero household size entries
valid_nonzero_idxs = findall((compiled_surveys.hh_size).!=0)
compiled_surveys = compiled_surveys[valid_nonzero_idxs,:]

# Save outputs
CSV.write("datasets/itn_hh_surveydata_complete.csv", compiled_surveys)
# CSV.write(datasets_dir*"itn_hh_surveydata_complete.csv", compiled_surveys)