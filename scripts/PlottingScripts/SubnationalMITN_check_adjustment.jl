"""
Author: Eugene Tan
Date Created: 27/3/2025
Last Updated: 27/3/2025
Script to plot adjusted vs non adjusted NPC time series for national estimates in subnational model.
Used to check if general properties of national is inherited by subnational.
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/read_toml.jl")

# %% Import Public Packages
using DataFrames
using Missings
using JLD2
using ProgressBars
using CSV
using StatsBase

# %% Plot packages
using LaTeXStrings
using Plots
using Measures

# %%
ISO_list = String.(CSV.read(RAW_DATASET_DIR*"/ISO_list.csv", DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)
fig_collection = []

# %% Load Data
for ISO_i in ProgressBar(1:length(filt_ISOs))
    ISO = filt_ISOs[ISO_i]
    nat_reg_data = JLD2.load(OUTPUT_DIR*"draws/national/crop_access/$(ISO)_2000_2023_post_crop_access.jld2")
    subnat_reg_data = JLD2.load(OUTPUT_DIR*"draws/subnational/$(ISO)_SUBNAT_draws.jld2")

    
    # %% Get meta_data and loop ranges
    n_regions = length(subnat_reg_data["merged_outputs"])
    YEAR_START = subnat_reg_data["YEAR_START_NAT"]
    YEAR_END = subnat_reg_data["YEAR_END"]
    n_months = 12*(YEAR_END-YEAR_START+1)

    # %% Declare storage variables for time series
    UNADJ_SUBNAT_Γ_MONTHLY_TOTAL_mean = zeros(n_regions, n_months)
    ADJ_SUBNAT_Γ_MONTHLY_TOTAL_mean = zeros(n_regions, n_months)

    # %% Extract national and subnational estimates for country
    # National
    NAT_Γ_MONTHLY_TOTAL_mean = nat_reg_data["Γ_MONTHLY_mean_TOTAL"]

    # Subnational
    for i in 1:n_regions
        UNADJ_SUBNAT_Γ_MONTHLY_TOTAL_mean[i,:] = mean(subnat_reg_data["merged_outputs"][i]["Γ_MONTHLY_TOTAL_samples"], dims = 1)[1,:]
        ADJ_SUBNAT_Γ_MONTHLY_TOTAL_mean[i,:] = mean(subnat_reg_data["merged_outputs"][i]["ADJ_Γ_MONTHLY_TOTAL_samples"], dims = 1)[1,:]
    end

    plot((ADJ_SUBNAT_Γ_MONTHLY_TOTAL_mean./(repeat(nat_reg_data["POPULATION_MONTHLY"],1,n_regions)'))')

    subnat_reg_data["merged_outputs"][1]["NPC_MONTHLY_TOTAL_mean"]
    

    UNADJ_Γ_MONTHLY_TOTAL_mean = sum(UNADJ_SUBNAT_Γ_MONTHLY_TOTAL_mean, dims = 1)[1,:]
    ADJ_Γ_MONTHLY_TOTAL_mean = sum(ADJ_SUBNAT_Γ_MONTHLY_TOTAL_mean, dims = 1)[1,:]

    # %% Plot Comparison
    # Plot visual settings
    pythonplot()
    theme(:vibrant)
    lw = 1.2

    # Define Net colors
    colors = [  colorant"#005684",
            colorant"#00976A",
            colorant"#E72A3D",
            colorant"#F7B801",
            colorant"#7018B3",
            ]

    fig = plot(title = "$(ISO) Net Crop", 
                    xlims = (-1,290),
                    xticks = (1:12:n_months,YEAR_START:YEAR_END), xtickfontrotation = 90,
                    xlabel = "Year", ylabel = "Net Crop (mil)", legend = :topleft)
    vline!(fig, 1:12:n_months, 
            linecolor = colorant"#45332C", linestyle = :dash, linewidth = 1, linealpha = 0.2,
            label = nothing)
    vline!(fig, [(2010-YEAR_START)*12+1], label = nothing, color = :red)

    plot!(fig, 1:n_months, NAT_Γ_MONTHLY_TOTAL_mean./nat_reg_data["POPULATION_MONTHLY"],
            linecolor = colors[1], linewidth = lw,
            label= "National Estimate")
    plot!(fig, 1:n_months, UNADJ_Γ_MONTHLY_TOTAL_mean./nat_reg_data["POPULATION_MONTHLY"],
            linecolor = colors[2], linewidth = lw,
            label= "Unadj. Subnat Estimate")
    plot!(fig, 1:n_months, ADJ_Γ_MONTHLY_TOTAL_mean./nat_reg_data["POPULATION_MONTHLY"],
            linecolor = colors[3], linewidth = lw,
            label= "Adj. Subnat Estimate")

    push!(fig_collection, fig)
end

# %%
fig = plot(fig_collection..., layout = (6,8), size = (2880,1620))

savefig(fig, OUTPUT_PLOTS_DIR*"nat_vs_subnat_netcrop.pdf")