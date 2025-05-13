"""
Author: Eugene Tan
Date Created: 22/8/2024
Last Updated: 24/8/2024
This script calculates household nets per capita values from survey metrics and 
applies adjustment to the standard errors. Outputs are saved in a csv file.
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Data Wrangling
using CSV
using DataFrames

# %% Mathematics packages
using LinearAlgebra
using StatsBase
using Distributions

# %% Custom modules
using DateConversions

# %% Define Directories
datasets_dir = RAW_DATASET_DIR
dataprep_dir = OUTPUT_DATAPREP_DIR

# %% Set start and end year bounds for reference
YEAR_START = YEAR_REF_START
YEAR_END = YEAR_NAT_END

# %% Import required datasets and legend/keys
country_codes_key = CSV.read(datasets_dir*COUNTRY_CODES_FILENAME, DataFrame)
master_survey_data = CSV.read(dataprep_dir*HOUSEHOLD_SURVEY_DATA_FILENAME, DataFrame)

# %% Get list of surveyIds that need to be extracted
surveyid_list = unique(master_survey_data.SurveyId)

# %% Inflation factor for Std Error
τ = STD_ERR_TAU

function inflation_factor(n; saturation_size = DEFAULT_INFLATION_SAT_SIZE)
    return 1 ./sqrt((min.(1, n/saturation_size)))
end


# %% Compute metrics for each survey id
survey_summaries = []

for i in 1:length(surveyid_list)
    # Extract metadata
    surveyid = String(surveyid_list[i])
    survey_data = master_survey_data[findall(master_survey_data.SurveyId .== surveyid),:]
    surveycountry = String(survey_data.Country[1])
    if ismissing(survey_data.ISO[1]) # missing data entry. Skip data
        continue
    end
    surveyiso = String(survey_data.ISO[1])

    # Calculate month indexing based on YEAR_START
    monthidxs = monthyear_to_monthidx.(survey_data.interview_month, survey_data.interview_year; YEAR_START = YEAR_START)

    # identify all unique month indexes to group survey entries into
    unique_months = unique(monthidxs)


    # Storage variables
    months = []
    HH_size_mean = []
    NPC_mean = []
    cITNPC_mean = []
    LLINPC_mean = []
    NPC_unadj_se = []
    use_mean = []
    sample_sizes = []


    # Calculate NPC means and unadj standard errors for each month in the survey
    for monthidx in unique_months
        # Extract data
        rowidx = findall(monthidxs .== monthidx)
        n_sample = length(rowidx)
        month_survey_data = survey_data[rowidx,:]
        month_weights = month_survey_data.hh_sample_wt

        # Check if no entry for survey sample weights (common in MICS5)
        if sum(month_weights.==0) == length(month_weights)
            # If no weight, just assume equally weighted
            month_weights .= 1
        end

        # Calculate mean HH_size
        HH_size = dot(month_weights, month_survey_data.hh_size)/sum(month_weights)
        
        # Calculate NPC metrics
        NPC_month_entries = (month_survey_data.n_itn)./(month_survey_data.hh_size)
        NPC_month_mean = sum(month_survey_data.n_itn.*month_weights)./sum(month_survey_data.hh_size.*month_weights) #dot(month_weights, NPC_month_entries)/sum(month_weights)
        NPC_month_unadjusted_SE = sqrt((dot(month_weights, (NPC_month_entries.-NPC_month_mean).^2)/sum(month_weights))/n_sample)

        cITNPC_month_entries = (month_survey_data.n_citn)./(month_survey_data.hh_size)
        cITNPC_month_mean = sum(month_survey_data.n_citn.*month_weights)./sum(month_survey_data.hh_size.*month_weights) #dot(month_weights, cITNPC_month_entries)/sum(month_weights)

        LLINPC_month_entries = (month_survey_data.n_llin)./(month_survey_data.hh_size)
        LLINPC_month_mean = sum(month_survey_data.n_llin.*month_weights)./sum(month_survey_data.hh_size.*month_weights) #dot(month_weights, LLINPC_month_entries)/sum(month_weights)

        # Calculate mean use rate
        use_month_mean = dot(month_weights, month_survey_data.n_slept_under_itn)/dot(month_weights, month_survey_data.hh_size)

        # Push to storage variable
        if NPC_month_unadjusted_SE > 0
            push!(months, monthidx)
            push!(HH_size_mean, HH_size)
            push!(NPC_mean, NPC_month_mean)
            push!(cITNPC_mean, cITNPC_month_mean)
            push!(LLINPC_mean, LLINPC_month_mean)
            push!(NPC_unadj_se, NPC_month_unadjusted_SE)
            push!(sample_sizes, n_sample)
            push!(use_mean, use_month_mean)
        end
    end

    if isempty(months) # no usable month data
        continue
    end

    # Calculate inflation factor for adjusted standard error
    NPC_adj_se = zeros(length(months))
    projected_samples = zeros(length(months), sum(sample_sizes))
    
    for j in 1:length(months)
        month_ref = months[j]

        month_sample_vals = Float64[]
        
        for k in 1:length(months)
            month_idx = months[k]
            attrition_factor = exp(log(0.5)*((month_ref - month_idx)/12)*(1/τ))
            month_sample_vals = vcat(month_sample_vals,rand(Normal(NPC_mean[k], NPC_unadj_se[k]), sample_sizes[k]).*attrition_factor)
        end

        projected_samples[j,:] = month_sample_vals
    end
    projected_samples = max.(projected_samples, 0)

    NPC_adj_se = inflation_factor.(sample_sizes).*std(projected_samples, dims = 2)[:]
    real_months = zeros(Int,length(months))
    real_years = zeros(Int,length(months))
    for k in 1:length(months)
        real_months[k], real_years[k] = monthidx_to_monthyear(months[k])
        real_years[k] = real_years[k] + YEAR_START - 1# Original function outputs years starting from 0
    end

    summarised_survey_data = DataFrame("SurveyId" => repeat([surveyid],length(months)), 
                                    "Country" => repeat([surveycountry],length(months)), 
                                    "ISO" => repeat([surveyiso],length(months)),
                                    "month" => real_months, 
                                    "year" => real_years,
                                    "sample_size" => sample_sizes,
                                    "hh_size_mean" => HH_size_mean, 
                                    "NPC_mean" => NPC_mean,
                                    "NPC_unadj_se" => NPC_unadj_se,
                                    "NPC_adj_se" => NPC_adj_se,
                                    "cITNPC_mean" => cITNPC_mean,
                                    "LLINPC_mean" => LLINPC_mean,
                                    "use_mean" => use_mean)
    
    push!(survey_summaries, summarised_survey_data)
end

# %% Compile data and save
summarised_data = vcat(survey_summaries...)

CSV.write(dataprep_dir*HOUSEHOLD_NAT_SUMMARY_DATA_FILENAME, summarised_data)

println("JULIA EXTRACTION COMPLETED :D")
flush(stdout)