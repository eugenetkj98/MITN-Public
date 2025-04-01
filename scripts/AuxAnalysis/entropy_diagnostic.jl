"""
Author: Eugene Tan
Date Created: 3/10/2024
Last Updated: 3/10/2024
Need to write documentation
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import Public Packages
using JLD2
using CSV
using DataFrames

# Maths packages
using LinearAlgebra
using StatsBase

# Plotting Functions
using Plots

# Helper packages
using DateConversions

# %% Analysis Year ranges
YEAR_START = 2000
YEAR_END = 2023

# %%
surveydata = CSV.read("datasets/itn_hh_surveydata_complete.csv", DataFrame)

# %% Get list of countries to include
ISO_list = String.(CSV.read(raw"C:\Users\ETan\Documents\Prototype Analyses\itn-updated\datasets\ISO_list.csv", DataFrame)[:,1])
exclusion_ISOs = ["CPV","BWA","GNQ","DJI","ETH","SOM","ZAF","SSD"]
ISO_list_filt = setdiff(ISO_list, exclusion_ISOs)

# %% Calculate bool filter to only include countries in ISO_list_filt
country_bool = zeros(Bool,size(surveydata)[1])

for i in 1:length(country_bool)
    country_bool[i] = (surveydata[i,"ISO"] ∈ ISO_list_filt)
end

# %% Filter out survey data to only take DHS (i.e. those with latitude and longitude values)
nonmissing_entries = surveydata[findall(Vector(.!ismissing.(surveydata.latitude)).*
                                            (surveydata.interview_year.>=YEAR_START).*
                                            (surveydata.interview_year.<=YEAR_END).*
                                            (country_bool)),:]

# %% Get list of unique SurveyIds
surveyid_list = unique(nonmissing_entries.SurveyId)

# %%
# Gaussian Helper function
gaussian(x, μ, σ) = sqrt(1/(2*pi*σ^2)).*exp.(-((x.-μ).^2)./(2 .*σ.^2))
mv_gaussian(x, μ, Σ) = ((2pi)^(-length(μ)/2))*(det(Σ)^(-1/2)).*exp.((-1/2)*(x-μ)' * inv(Σ) * (x-μ))

function survey_entropy(x_centres, y_centres, pos; α = 3)
    # Calculate reference Gaussian kernel pdf
    gaussian_unit_image = zeros(length(x_centres), length(y_centres))
    μ_ref = [(x_centres[end]+x_centres[1])/2, (y_centres[end]+y_centres[1])/2]

    # Get widths of grid elements
    dx = x_centres[2]-x_centres[1]
    dy = y_centres[2]-y_centres[1]

    # Define Covariance matrix for 2D Gaussian, scaled by spread term α
    Σ = Matrix(Diagonal([α*dx, α*dy].^2))
    for i in 1:length(x_centres)
        for j in 1:length(y_centres)
            r = [x_centres[i], y_centres[j]]
            gaussian_unit_image[i,j] = mv_gaussian(r, μ_ref, Σ)
        end
    end

    # Calculate superposition of Gaussians based on locations in list of survey entries
    resolution = length(x_centres)
    z_pdf_components = zeros(size(pos)[1], resolution, resolution)

    Threads.@threads for n in ProgressBar(1:size(pos)[1], leave = false)
        μ = pos[n,:]
        Δμ = μ.-μ_ref

        Δi = Int(round(Δμ[1]./dx))
        Δj = Int(round(Δμ[2]./dy))

        shifted_image = zeros(resolution, resolution)

        if Δi<0 && Δj<0
            x_idxs = ((-Δi+1):resolution)
            y_idxs = ((-Δj+1):resolution)
            shifted_image[1:length(x_idxs), 1:length(y_idxs)] = gaussian_unit_image[x_idxs, y_idxs]
        elseif Δi>=0 && Δj<0
            x_idxs = (1:resolution-Δi)
            y_idxs = ((-Δj+1):resolution)
            shifted_image[end-length(x_idxs)+1:end, 1:length(y_idxs)] = gaussian_unit_image[x_idxs, y_idxs]
        elseif Δi<0 && Δj>=0
            x_idxs = ((-Δi+1):resolution)
            y_idxs = (1:resolution-Δj)
            shifted_image[1:length(x_idxs), end-length(y_idxs)+1:end] = gaussian_unit_image[x_idxs, y_idxs]
        else
            x_idxs = (1:resolution-Δi)
            y_idxs = (1:resolution-Δj)
            shifted_image[end-length(x_idxs)+1:end, end-length(y_idxs)+1:end] = gaussian_unit_image[x_idxs, y_idxs]
        end

        z_pdf_components[n,:,:] = shifted_image
    end

    # Sum and normalise components
    z_pdf = sum(z_pdf_components, dims = 1)[1,:,:]
    z_pdf = z_pdf./sum(z_pdf)

    # Calculate entropy of distribution
    H = .-z_pdf.*log.(z_pdf)

    # Correct for Infs and NaNs
    H[findall(isinf.(H))] .= 0
    H[findall(isnan.(H))] .= 0

    # Sum for total entropy
    H_survey = sum(H)

    return z_pdf, H_survey
end

# %% Calculate relative entropies for survey by month
function survey_homog_score(survey_filt; resolution = 100, YEAR_START = 2000, α = 3)

    # Extract position data from specific surveyid DataFrame
    pos = Matrix(survey_filt[:,["latitude","longitude"]])

    # Get integration bounds unique to country for entire survey
    x_bounds = quantile(pos[:,1], [0.01, 0.99])
    y_bounds = quantile(pos[:,2], [0.01, 0.99])

    # Define resolution ratio to divide integration region and calculate bins, grid size
    x_bins = LinRange(x_bounds[1], x_bounds[2], resolution+1)
    y_bins = LinRange(y_bounds[1], y_bounds[2], resolution+1)
    x_centres = (x_bins[2:end].+x_bins[1:end-1])./2
    y_centres = (y_bins[2:end].+y_bins[1:end-1])./2

    # Calculate entropy for survey
    z_survey, H_survey = survey_entropy(x_centres, y_centres, pos, α = α)

    # Get list of month x_idxs and number of unique months
    monthidxs = zeros(Int, size(survey_filt)[1])
    survey_filt[i,:interview_month]
    for i in 1:length(monthidxs)
        monthidxs[i] = monthyear_to_monthidx(survey_filt[i,"interview_month"], survey_filt[i,"interview_year"], YEAR_START = YEAR_START)
    end
    unique_months = unique(monthidxs)

    # Calculate entropy for each unique month Select
    H_months = zeros(length(unique_months))
    KL_months = zeros(length(unique_months))

    for i in 1:length(unique_months)
        monthidx = unique_months[i]
        pos_month = Matrix(survey_filt[findall(monthidxs .== monthidx),["latitude", "longitude"]])
        z_month, H_month = survey_entropy(x_centres, y_centres, pos_month, α = α)
        H_months[i] = H_month

        # # Calculate KL divergence
        # KL_matrix = z_survey.*log.(z_survey./z_month)
        # KL_matrix[findall(isnan.(KL_matrix))].=0
        # KL_matrix[findall(isinf.(KL_matrix))].= 10000
        # KL_months[i] = sum(KL_matrix)

        # Calculate JS divergence
        P = z_survey
        Q = z_month
        M = (z_survey .+ z_month)./2
        
        KL_PM = P.*log2.(P./M)
        KL_PM[findall(isnan.(KL_PM))].=0

        KL_QM = Q.*log2.(Q./M)
        KL_QM[findall(isnan.(KL_QM))].=0
        
        KL_months[i] = sum((KL_PM .+ KL_QM)./2)
    end

    # Return normalised entropy scores and monthidx as scores
    norm_H_months = H_months./H_survey

    return unique_months, norm_H_months, H_survey, KL_months
end

# %%

# BUFFER_YEARS = 5

# BUFFER_YEAR_START = YEAR_START-BUFFER_YEARS
# BUFFER_MONTHS = BUFFER_YEARS*12

# n_months =  (YEAR_END-BUFFER_YEAR_START+1)*12
n_months =  (YEAR_END-YEAR_START+1)*12

score_breakdown_matrix = missings(Float64, length(surveyid_list), n_months)
survey_entropies = zeros(length(surveyid_list))
survey_KLs = missings(Float64, length(surveyid_list), n_months)

# Calculate normalised entropy score for each survey
for i in ProgressBar(1:length(surveyid_list))

    # Get surveyid and filter out relevant data
    surveyid = surveyid_list[i]
    survey_filt = nonmissing_entries[findall(nonmissing_entries.SurveyId.==surveyid),:]

    # Calculate scores and store
    unique_months, norm_H_months, H_survey, KL_months = survey_homog_score(survey_filt; YEAR_START = YEAR_START)
    score_breakdown_matrix[i, unique_months] .= KL_months#norm_H_months
    survey_entropies[i] = H_survey
    # survey_KLs[i, unique_months] = KL_months
end


# %%
# Get Country based on surveyid_list
country_list = []

for surveyid in surveyid_list
    push!(country_list,nonmissing_entries[findfirst(nonmissing_entries[:,"SurveyId"].==surveyid),"ISO"])
end

# %%
country_score_matrix = missings(Float64, length(ISO_list_filt), n_months)

# %%

for i in 1:size(score_breakdown_matrix)[1]
    monthidxs = findall(.!ismissing.(score_breakdown_matrix[i,:]))
    country_idx = findfirst(ISO_list_filt .== country_list[i])
    country_score_matrix[country_idx, monthidxs] = score_breakdown_matrix[i,monthidxs]
end

# %% Make summary plot of condition
# fig = plot()
# for i in 1:size(country_score_matrix)[1]
#     nonzero_idxs = findall(country_score_matrix[i,:].>0)
#     scatter!(fig, nonzero_idxs, i.*ones(length(nonzero_idxs)), zcolor = country_score_matrix[i,nonzero_idxs], markerstrokewidth = 0,
#             label = nothing)
# end
# fig
pythonplot()
gr()
theme(:vibrant)

fig = plot(xticks = (1:12:n_months, YEAR_START:YEAR_END), xrotation = 90,
            yticks = (1:length(ISO_list_filt), ISO_list_filt),
            xlabel = "Years", colorbartitle = "Survey Homogeneity", clims = (0,1),
            size = (700,500))
heatmap!(fig, 1 .-country_score_matrix, cmap = cgrad(["#FC5549", "#34BD7D"]), yflip = true)



# %%
savefig(fig, "SurveyHomogeneity_JSD.pdf")

# %% Generate survey component figs for illustrative purposes
gr()
surveyid = "NG2013DHS"
# surveyid = "MD2016MIS"
# surveyid = "MZ2022DHS"
resolution = 100
α = 3

# Get full survey data
survey_filt = nonmissing_entries[findall(nonmissing_entries.SurveyId.==surveyid),:]

# Extract position data from specific surveyid DataFrame
pos = Matrix(survey_filt[:,["latitude","longitude"]])

# Get integration bounds unique to country for entire survey
x_bounds = quantile(pos[:,1], [0.01, 0.99])
y_bounds = quantile(pos[:,2], [0.01, 0.99])

# Define resolution ratio to divide integration region and calculate bins, grid size
x_bins = LinRange(x_bounds[1], x_bounds[2], resolution+1)
y_bins = LinRange(y_bounds[1], y_bounds[2], resolution+1)
x_centres = (x_bins[2:end].+x_bins[1:end-1])./2
y_centres = (y_bins[2:end].+y_bins[1:end-1])./2


# %% Calculate full survey heatmap
z_survey, H_survey = survey_entropy(x_centres, y_centres, pos, α = α)

# %% Calculate heatmap for each separate month and save
z_survey_months = []
H_survey_months = []
KL_survey_months = []

# Get list of month x_idxs and number of unique months
monthidxs = zeros(Int, size(survey_filt)[1])
for i in 1:length(monthidxs)
    monthidxs[i] = monthyear_to_monthidx(survey_filt[i,"interview_month"], survey_filt[i,"interview_year"], YEAR_START = YEAR_START)
end
unique_months = unique(monthidxs)

# %%
month_scores = survey_homog_score(survey_filt)

# Calculate entropy for each unique month Select
for i in 1:length(unique_months)
    monthidx = unique_months[i]
    pos_month = Matrix(survey_filt[findall(monthidxs .== monthidx),["latitude", "longitude"]])
    z_month, H_month = survey_entropy(x_centres, y_centres, pos_month, α = α)
    push!(z_survey_months, z_month)
    push!(H_survey_months, H_month)

    # # Calculate KL divergence
    # KL_matrix = z_survey.*log.(z_survey./z_survey_months[i])
    # KL_matrix[findall(isnan.(KL_matrix))].=0
    # KL_matrix[findall(isinf.(KL_matrix))].=0
    # KL_div = sum(KL_matrix)
    # push!(KL_survey_months, KL_div)

end

# %% Make Plots
heatmap_collections = []
clims = (0,0.002)

for i in 1:length(month_scores)
    fig = heatmap(x_centres, y_centres, z_survey_months[i], yflip = true,
                    clims = clims, title = "Month $(i), H_norm = $(round(H_survey_months[i]./H_survey, digits = 3)), \nJSD = $(round(month_scores[2][i], digits = 3))")
    push!(heatmap_collections, fig)
end

fig = heatmap(x_centres, y_centres, z_survey, yflip = true,
                clims = clims, title = "Full Survey $(surveyid)")

push!(heatmap_collections, fig)
# %%

# fig_scat = scatter(H_survey_months./H_survey)
# fig_scat = scatter!(fig_scat, H_survey_months./H_survey)

# %%
fig_all = plot(heatmap_collections[[1,2,5,3,4]]..., layout = (2,3), size = (1200,600), margins = 5mm)


# %%
savefig(fig_all, "$(surveyid)_entropies.pdf")
# %%


