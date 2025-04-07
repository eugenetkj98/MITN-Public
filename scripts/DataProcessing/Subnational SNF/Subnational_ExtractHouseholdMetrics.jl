"""
Author: Eugene Tan
Date Created: 11/11/2024
Last Updated: 11/11/2024
Extract observed net crop estimates for household survey data by area_id.
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Import Public Packages
using JLD2
using CSV
using DataFrames
using DateConversions
using ProgressBars
using LinearAlgebra
using StatsBase

# %%
function inflation_factor(n; saturation_size = DEFAULT_INFLATION_SAT_SIZE)
    return 1 ./sqrt((min.(1, n/saturation_size)))
end

# %% Define paths
dataset_dir = RAW_SUBNAT_DATASET_DIR
dataprep_dir = OUTPUT_DATAPREP_DIR

# Admin 1 ID Legend
id_legend_filename = ADMIN1_AREAID_LEGEND_FILENAME

# National aggregated survey data
nat_npc_monthly_data_filename = HOUSEHOLD_NAT_SUMMARY_DATA_FILENAME

# Full household survey data entries
hh_survey_data_filename = OUTPUT_DATAPREP_DIR*HOUSEHOLD_SURVEY_DATA_FILENAME #"datasets/subnational/itn_hh_surveydata_complete_subnat.csv"

# Paths to save output
output_path = OUTPUT_DATAPREP_DIR
output_filename = HOUSEHOLD_SUBNAT_SUMMARY_DATA_FILENAME

# %% Import datasets
master_id_legend = CSV.read(dataset_dir*id_legend_filename, DataFrame)
nat_npc_monthly_data = CSV.read(OUTPUT_DATAPREP_DIR*nat_npc_monthly_data_filename, DataFrame)
subnat_full_survey_data = CSV.read(hh_survey_data_filename, DataFrame)
# Exclude data entries with missing ISO
subnat_full_survey_data = subnat_full_survey_data[.!ismissing.(subnat_full_survey_data.ISO),:]

# %% Choose Country and Analysis Settings
sat_survey_size = DEFAULT_INFLATION_SAT_SIZE
YEAR_START = YEAR_NAT_START # Start year for data parsing
YEAR_END = YEAR_NAT_END # End year for data parsing

# %% Get list of countries
ISO_list = unique(subnat_full_survey_data.ISO)

# Storage variable for summarised data
dataframerows = []

# %% Select country
for ISO_idx in ProgressBar(1:length(ISO_list))
    ISO = ISO_list[ISO_idx]

    # %% Filter survey data to desired country
    
    findall(ismissing.(subnat_full_survey_data.ISO .== ISO))
    country_full_survey_data = subnat_full_survey_data[subnat_full_survey_data.ISO .== ISO,:]
    country_name = country_full_survey_data.Country[1]

    # Get list of unique surveys
    country_surveyids = unique(country_full_survey_data.SurveyId)

    # Focus on a specific survey
    for surveyid_idx in 1:length(country_surveyids)
        surveyid = country_surveyids[surveyid_idx]

        # Filter survey entries by id
        survey_entries = country_full_survey_data[(country_full_survey_data.SurveyId .== surveyid),:]

        # %% Get list of admin1 area_ids
        country_legend = master_id_legend[master_id_legend.ISO.==ISO,:]
        admin1_names = country_legend.Name_1
        n_admin1 = length(admin1_names)

        # %% For each admin1_region, find netcrop estimate
        for subnat_i in 1:n_admin1
            admin1_name = admin1_names[subnat_i]
            area_id = country_legend[country_legend.Name_1 .== admin1_name, "area_id"][1]

            # Check if there is admin1 level information
            nonmissing_idxs = findall(.!ismissing.(survey_entries.area_id))
            full_admin1_survey_entries = survey_entries[nonmissing_idxs,:]
            
            if area_id ∈ unique(full_admin1_survey_entries.area_id) # i.e. there is admin1 specific information
                filt_survey_entries = full_admin1_survey_entries[(full_admin1_survey_entries.area_id .== area_id) .&
                                                                   (full_admin1_survey_entries.interview_year .>= YEAR_START) .&
                                                                   (full_admin1_survey_entries.interview_year .<= YEAR_END),:]

                monthidx_vals = monthyear_to_monthidx.(filt_survey_entries.interview_month, filt_survey_entries.interview_year, YEAR_START = YEAR_START)
                unique_months = unique(monthidx_vals)

                # %% Calculate net crop estimates for each monthidx
                for monthidx_i in 1:length(unique_months)
                    monthidx = unique_months[monthidx_i]

                    month_val, ref_year_val = monthidx_to_monthyear(monthidx)
                    year_val = ref_year_val+ YEAR_START-1

                    # Calculate net crop estimate
                    monthyear_entries = filt_survey_entries[(filt_survey_entries.interview_month .== month_val) .& 
                                                            (filt_survey_entries.interview_year .== year_val),:]

                    NPC_entries = (monthyear_entries.n_itn)./(monthyear_entries.hh_size)
                    NPC_month_mean = sum(NPC_entries.*(monthyear_entries.hh_sample_wt))/sum(monthyear_entries.hh_sample_wt)

                    weighted_std = sqrt(sum(((NPC_entries .- NPC_month_mean).^2).*(monthyear_entries.hh_sample_wt))/sum(monthyear_entries.hh_sample_wt))
                    weighted_stderr = weighted_std/sqrt(length(NPC_entries))
                    hh_size_mean = sum(monthyear_entries.hh_size.*(monthyear_entries.hh_sample_wt))/sum(monthyear_entries.hh_sample_wt)
                    NPC_adj_se = weighted_stderr*inflation_factor(length(NPC_entries), saturation_size = sat_survey_size)

                    # Calculate use estimate
                    use_month_mean = dot(monthyear_entries.hh_sample_wt, monthyear_entries.n_slept_under_itn)/dot(monthyear_entries.hh_sample_wt, monthyear_entries.hh_size)

                    # Push to dataframerows but first check if sample size > 10 (FOR ERROR PROPAGATION REASONS AND AVOIDING NaNs)
                    if length(NPC_entries) > 10
                        push!(dataframerows, DataFrame(SurveyId = surveyid, Country = country_name, ISO = ISO,
                                                        admin1_name = admin1_name, area_id = area_id,
                                                        month = month_val, year = year_val,
                                                        source = "Subnational",
                                                        samplesize = length(NPC_entries),
                                                        hh_size_mean = hh_size_mean,
                                                        NPC_mean = NPC_month_mean,
                                                        NPC_adj_se = NPC_adj_se,
                                                        use_mean = use_month_mean))
                    end
                end

            else # There is no admin1_data available, so just copy from national extraction data
                filt_npc_monthly_data = nat_npc_monthly_data[nat_npc_monthly_data.SurveyId .== surveyid,:]
                for row_i in 1:size(filt_npc_monthly_data)[1]
                    datarow = filt_npc_monthly_data[row_i,:]

                    # Extract entry from national summary data and push to list of of dataframerows
                    push!(dataframerows, DataFrame(SurveyId = surveyid, Country = country_name, ISO = ISO,
                                            admin1_name = admin1_name, area_id = area_id,
                                            month = datarow.month, year = datarow.year,
                                            source = "National",
                                            samplesize = datarow.sample_size,
                                            hh_size_mean = datarow.hh_size_mean,
                                            NPC_mean = datarow.NPC_mean,
                                            NPC_adj_se = datarow.NPC_adj_se,
                                            use_mean = datarow.use_mean))
                end
            end
        end
    end
end

# %% Join dataframe rows and save data
subnat_npc_monthly_data = vcat(dataframerows...)
mkpath(output_path)
CSV.write(output_path*output_filename, subnat_npc_monthly_data)