"""
Author: Eugene Tan
Date Created: 22/10/2024
Last Updated: 22/10/2024
Makes reference plot to demonstrate how ExtractHouseholdMetrics.jl works
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Data Wrangling
using CSV
using DataFrames

# %% Mathematics packages
using LinearAlgebra
using StatsBase
using Distributions

# %% Custom modules
using Plots
using DateConversions

# %% Define Directories
datasets_dir = "datasets/"

# %% Set start and end year bounds for reference
YEAR_START = 1950
YEAR_END = 2023

# %% Import required datasets and legend/keys
country_codes_key = CSV.read(datasets_dir*"country_codes.csv", DataFrame)
master_survey_data = CSV.read(datasets_dir*"itn_hh_surveydata_complete.csv", DataFrame)
# country_codes_key = CSV.read("datasets/country_codes.csv", DataFrame)
# master_survey_data = CSV.read("datasets/itn_hh_surveydata_complete.csv", DataFrame)

# %% Get list of surveyIds that need to be extracted
surveyid_list = unique(master_survey_data.SurveyId)

# %%
τ = 1.5

function inflation_factor(n; saturation_size = 4000)
    return 1 ./sqrt((min.(1, n/saturation_size)))
end


# %% Compute metrics for each survey id
survey_summaries = []



# %%
i = 159

# Extract metadata
surveyid = String(surveyid_list[i])
survey_data = master_survey_data[findall(master_survey_data.SurveyId .== surveyid),:]
surveycountry = String(survey_data.Country[1])
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
    NPC_month_mean = dot(month_weights, NPC_month_entries)/sum(month_weights)
    NPC_month_unadjusted_SE = sqrt((dot(month_weights, (NPC_month_entries.-NPC_month_mean).^2)/sum(month_weights))/n_sample)

    cITNPC_month_entries = (month_survey_data.n_citn)./(month_survey_data.hh_size)
    cITNPC_month_mean = dot(month_weights, cITNPC_month_entries)/sum(month_weights)

    LLINPC_month_entries = (month_survey_data.n_llin)./(month_survey_data.hh_size)
    LLINPC_month_mean = dot(month_weights, LLINPC_month_entries)/sum(month_weights)
    if NPC_month_unadjusted_SE > 0
        push!(months, monthidx)
        push!(HH_size_mean, HH_size)
        push!(NPC_mean, NPC_month_mean)
        push!(cITNPC_mean, cITNPC_month_mean)
        push!(LLINPC_mean, LLINPC_month_mean)
        push!(NPC_unadj_se, NPC_month_unadjusted_SE)
        push!(sample_sizes, n_sample)
    end
end

NPC_adj_se = zeros(length(months))

projected_samples_single = []
projected_samples_full = []
    
for j in 1:length(months)
    month_ref = months[j]
    idxs = findall(monthidxs .== month_ref)

    projection_single = zeros(length(months), length(idxs))
    projection_full = zeros(length(months), length(idxs))

    single_sample_point = rand(Normal(NPC_mean[j], 0.0), sample_sizes[j])
    full_sample_points = rand(Normal(NPC_mean[j], NPC_unadj_se[j]), sample_sizes[j])
    for k in 1:length(months)
        month_idx = months[k]
        attrition_factor = exp(log(0.5)*((month_ref-month_idx)/12)*(1/τ))
        projection_single[k,:] = single_sample_point.*attrition_factor
        projection_full[k,:] = full_sample_points.*attrition_factor
    end

    push!(projected_samples_single, projection_single)
    push!(projected_samples_full, projection_full)
end
# %%
months[sortperm(months)]
sortperm(months)
plot(hcat(projected_samples_full[sortperm(months)]...)[:,1:10:end], legend = false)
# %%
projected_samples_full[3]



# %%
pythonplot()
theme(:vibrant)
ylims = (minimum(minimum.(projected_samples_full))-0.015, maximum(maximum.(projected_samples_full))+0.015)


fig1 = plot(ylims = ylims, xlabel = "Month Index", ylabel = "γ", title = "Projection of Point Estimates",
            legendfontsize = 10)

for j in 1:length(months)
    # scatter!(fig, months,projected_samples[j][:,1:10:end], legend = false, markersize = 2, color = j,
    #         alpha = 0.1)
    if j == 1
        label1 = "Projection"
        label2 = "Survey Estimate"
    else
        label1 = nothing
        label2 = nothing
    end
    plot!(fig1, months[sortperm(months)], projected_samples_single[j][sortperm(months),1], color = j,
            alpha = 1, linewidth = 1.5, label = label1, linestyle = :dash)
    scatter!(fig1, [months[j]], [NPC_mean[j]], color = j, markersize = 7, label = label2)
end
fig1

# %%

adjusted_se = inflation_factor.(sample_sizes).*std(hcat(projected_samples_full...), dims = 2)[:]

# %%
fig2 = plot(ylims = ylims, xlabel = "Month Index", ylabel = "γ", title = "Full Projection")

for j in 1:length(months)
    plot!(fig2, months[sortperm(months)], projected_samples_full[j][sortperm(months),1:20:end], color = j,
            alpha = 0.5, linewidth = 0.2, label = nothing, linestyle = :dash)
end

for j in 1:length(months)
    n_points = length(projected_samples_full[j][j,1:50:end])
    scatter!(fig2, repeat([months[j]], 1,n_points), projected_samples_full[j][j,1:50:end], color = j, markersize = 3, label = nothing)
end

for j in 1:length(months)
    plot!(fig2, repeat([months[j]], 1,2), [NPC_mean[j]-adjusted_se[j],NPC_mean[j]+adjusted_se[j]], color = j, 
            linewidth= 1, label = nothing)
end

fig2

# %%
fig = plot(fig1, fig2, layout = (1,2), size = (1000,400))
savefig(fig, "Variance_Inflation_Plot.pdf")




# %% Generate std samples and do calculations
i = 159
# Extract metadata
surveyid = String(surveyid_list[i])
survey_data = master_survey_data[findall(master_survey_data.SurveyId .== surveyid),:]
surveycountry = String(survey_data.Country[1])
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
    NPC_month_mean = dot(month_weights, NPC_month_entries)/sum(month_weights)
    NPC_month_unadjusted_SE = sqrt((dot(month_weights, (NPC_month_entries.-NPC_month_mean).^2)/sum(month_weights))/n_sample)

    cITNPC_month_entries = (month_survey_data.n_citn)./(month_survey_data.hh_size)
    cITNPC_month_mean = dot(month_weights, cITNPC_month_entries)/sum(month_weights)

    LLINPC_month_entries = (month_survey_data.n_llin)./(month_survey_data.hh_size)
    LLINPC_month_mean = dot(month_weights, LLINPC_month_entries)/sum(month_weights)
    if NPC_month_unadjusted_SE > 0
        push!(months, monthidx)
        push!(HH_size_mean, HH_size)
        push!(NPC_mean, NPC_month_mean)
        push!(cITNPC_mean, cITNPC_month_mean)
        push!(LLINPC_mean, LLINPC_month_mean)
        push!(NPC_unadj_se, NPC_month_unadjusted_SE)
        push!(sample_sizes, n_sample)
    end
end

NPC_adj_se = zeros(length(months))

projected_samples = zeros(length(months), sum(sample_sizes))

random_draws = []
for k in 1:length(months)
    push!(random_draws, rand(Normal(NPC_mean[k], NPC_unadj_se[k]), sample_sizes[k]))
end

for j in 1:length(months)
    month_ref = months[j]

    month_sample_vals = Float64[]
    
    for k in 1:length(months)
        month_idx = months[k]
        attrition_factor = exp(log(0.5)*((month_ref - month_idx)/12)*(1/τ))
        month_sample_vals = vcat(month_sample_vals,random_draws[k].*attrition_factor)
    end

    projected_samples[j,:] = month_sample_vals
end

projected_samples = max.(projected_samples, 0)

NPC_adj_se = inflation_factor.(sample_sizes).*std(projected_samples, dims = 2)[:]

# %% Get Index Intervals
month_cumul_idx = cumsum(sample_sizes)
sorted_project_samples = projected_samples[sortperm(months),:]

# %% Get plot ylims

ylims = (minimum(projected_samples)-0.015, maximum(projected_samples)+0.015)

# %%
fig1 = plot(ylims = ylims, xlabel = "Month Index", ylabel = "γ", title = "Projection of Point Estimates",
legendfontsize = 10)

for i in 1:length(months)
    if i == 1
        idx_1 = 1
        idx_2 = month_cumul_idx[i]
    else
        idx_1 = month_cumul_idx[i-1]+1
        idx_2 = month_cumul_idx[i]
    end

    middle_idx = round(Int,(idx_1+idx_2)/2)
    plot!(fig1, months[sortperm(months)], sorted_project_samples[:, middle_idx],
            legend = false,
            color = palette(:tab10)[i],
            linestyle = :dash, linewidth = 1)
end
for i in 1:length(months)
    if i == 1
        idx_1 = 1
        idx_2 = month_cumul_idx[i]
    else
        idx_1 = month_cumul_idx[i-1]+1
        idx_2 = month_cumul_idx[i]
    end

    middle_idx = round(Int,(idx_1+idx_2)/2)

    scatter!(fig1,[months[i]] ,[projected_samples[i, middle_idx]], markercolor = palette(:tab10)[i], markersize = 8)

end



fig1

# %% split into generated values into sample sample_sizes
fig2 = plot(ylims = ylims, xlabel = "Month Index", ylabel = "γ", title = "Full Projection")
for i in 1:length(months)
    if i == 1
        idx_1 = 1
        idx_2 = month_cumul_idx[i]
    else
        idx_1 = month_cumul_idx[i-1]+1
        idx_2 = month_cumul_idx[i]
    end

    plot!(fig2, months[sortperm(months)], sorted_project_samples[:, idx_1:30:idx_2], legend = false, 
            color = palette(:tab10)[i],
            alpha = 0.08, linewidth = 0.8)

    
end

for i in 1:length(months)
    scatter!(fig2, [months[i]],[NPC_mean[i]], markercolor = palette(:tab10)[i], markersize = 8)
end
for i in 1:length(months)
    plot!(fig2, repeat([months[i]], 1,2)[:], [NPC_mean[i]-NPC_adj_se[i]  , NPC_mean[i]+NPC_adj_se[i]],
            color= palette(:tab10)[i])
end
fig2


# %%
fig = plot(fig1, fig2, layout = (1,2), size = (1000,400))
savefig(fig, "Variance_Inflation_Plot.png")
