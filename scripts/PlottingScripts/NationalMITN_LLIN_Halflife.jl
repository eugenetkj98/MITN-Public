"""
Author: Eugene Tan
Date Created: 27/3/2025
Last Updated: 27/3/2025
Script to plot a figure of the LLIN halflives across a list of analysed countries. Compares BV vs MITN
"""
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Import Public Packages
using DataFrames
using Missings
using JLD2
using CSV
# %% Maths packages
using LinearAlgebra
using StatsBase

# %% Plot packages
using LaTeXStrings
using Plots
using Measures

# %% MITN Model packages
using NetAccessModel
using DateConversions
using NetLoss



# %% Define save paths
output_path = OUTPUT_PLOTS_DIR*"snf/national/summaries/"
mkpath(output_path)

# %% Country Data
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

BV_halflife = CSV.read(RAW_DATASET_DIR*BV_OUTPUTS_FILENAME, DataFrame)
country_code = CSV.read(RAW_DATASET_DIR*COUNTRY_CODES_FILENAME, DataFrame)

# %% Raw Survey data to find out how much data was available
survey_data = CSV.read(HOUSEHOLD_SURVEY_DIR*HOUSEHOLD_SURVEY_FILENAME, DataFrame)
survey_data = survey_data[.!ismissing.(survey_data.ISO),:] # Remove entries with missing ISO
# country_name = country_code[findfirst(country_code.ISO3 .== ISO),"Country"]

# %% Create storage variable for DataFrame rows
df_collection = []

# %%
for ISO in filt_ISOs

    # %% Load Datasets
    input_dict = load(OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
    regression_dict = load(OUTPUT_REGRESSIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")

    # %%
    # Get Metadata and define bounds
    NET_NAMES = regression_dict["NET_NAMES"]
    LLIN_idx = findfirst(NET_NAMES .== "LLIN")

    # Extract chain
    τ_chain = regression_dict["chain"][:,5:2:5+(length(NET_NAMES)-1)*2]
    κ_chain = regression_dict["chain"][:,6:2:6+(length(NET_NAMES)-1)*2]

    # Calculate halflifes from posterior sampled chain
    LLIN_halflife_chain = zeros(size(τ_chain)[1])

    for i in 1:size(τ_chain)[1]
        LLIN_halflife_chain[i] = net_life_compact(0.5, τ_chain[i,LLIN_idx], κ_chain[i,LLIN_idx])
    end

    # Calculate quantiles of halflifes
    mitn_halflife_lower, mitn_halflife_median, mitn_halflife_upper = quantile(LLIN_halflife_chain, [0.025, 0.5, 0.975])

    # Check if ISO was analysed by BV
    if ISO ∈ BV_halflife.ISO
        df = DataFrame(ISO = ISO,
                n_surveys = length(unique(survey_data[survey_data.ISO .==ISO, "SurveyId"])),
                mitn_halflife_lower = mitn_halflife_lower,
                mitn_halflife_median = mitn_halflife_median,
                mitn_halflife_upper = mitn_halflife_upper,
                bv_halflife_lower = BV_halflife[BV_halflife.ISO .== ISO,"halflife_lower"][1],
                bv_halflife_median = BV_halflife[BV_halflife.ISO .== ISO,"halflife_median"][1],
                bv_halflife_upper = BV_halflife[BV_halflife.ISO .== ISO,"halflife_upper"][1])
        push!(df_collection, df)
    else
        df = DataFrame(ISO = ISO,
                n_surveys = length(unique(survey_data[survey_data.ISO .==ISO, "SurveyId"])),
                mitn_halflife_lower = mitn_halflife_lower,
                mitn_halflife_median = mitn_halflife_median,
                mitn_halflife_upper = mitn_halflife_upper,
                bv_halflife_lower = NaN,
                bv_halflife_median = NaN,
                bv_halflife_upper = NaN)
        push!(df_collection, df)
    end
end

# %% Compile halflife dataframe entries and sort according to MITN median
model_halflifes = vcat(df_collection...)
model_halflifes = model_halflifes[sortperm(model_halflifes[:, "mitn_halflife_median"]),:]

# %% Make plot
# Plot visual settings
pythonplot()
theme(:vibrant)
sep = 0.15
lw = 2.8
ms = 6
col_bv = colorant"#005684"
col_mitn = colorant"#E72A3D"

# String labels
xtickstrings = []
for i in 1:size(model_halflifes)[1]
    push!(xtickstrings, "[$(model_halflifes[i,"n_surveys"])] $(model_halflifes[i,"ISO"])")
end

# Plot outputs
fig = plot(title = "LLIN Retention Halflives",
            xticks = (1:size(model_halflifes)[1], xtickstrings),
            xtickfontrotation = 90, ylims = (-0.05,8.05),
            grid = true, minorgrid = false, framestyle = :box,
            size = (1200,350), legend = :topleft,
            ylabel = "Years")

for i in 1:size(model_halflifes)[1]
    plot!(fig, [i,i].-sep, Array(model_halflifes[i,["bv_halflife_lower", "bv_halflife_upper"]]),
            linewidth = lw, linecolor = col_bv, label = nothing)
    plot!(fig, [i,i].+sep, Array(model_halflifes[i,["mitn_halflife_lower", "mitn_halflife_upper"]]),
            linewidth = lw, linecolor = col_mitn, label = nothing)
end
scatter!(fig, (1:size(model_halflifes)[1]).-sep, model_halflifes.bv_halflife_median,
            color = col_bv, markersize = ms, label = "BV")
scatter!(fig, (1:size(model_halflifes)[1]).+sep, model_halflifes.mitn_halflife_median,
            color = col_mitn, markersize = ms, label = "MITN")

savefig(fig, output_path*"snf_nat_halflife.pdf")

fig

