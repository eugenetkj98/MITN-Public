"""
Author: Eugene Tan
Date Created: 15/4/2025
Last Updated: 15/4/2025
Script to make plots of the predicted netcrop breakdown by type for all malaria endemic countries
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Prep environment and subdirectories
include(pwd()*"/scripts/dir_configs.jl")

# %% Load packages
using JLD2
using CSV
using DataFrames
using Plots

# %% Data Directories
netcrop_dir = "outputs/predictions/"
survey_data_dir = "outputs/extractions/crop/$(YEAR_NAT_START)_$(YEAR_NAT_END)/"

# Save directories
output_dir = OUTPUT_PLOTS_DIR*"snf/national/multitype_predictions/"
mkpath(output_dir)
# %% Get ISO List
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = ["CPV","ZAF"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Define Metadata
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END
n_months = (YEAR_END-YEAR_START+1)*12

# %% Plot visual settings
pythonplot()
theme(:vibrant)
ms = 4
# Color settings
total_net_col = :black;
total_dist_net_col = :blue;
colors = [  colorant"#005684",
                colorant"#00976A",
                colorant"#E72A3D",
                colorant"#F7B801",
                colorant"#7018B3",
                ];
component_alpha_mult = 0.5 # Multiplier to reduce component nets time series

# %% 
uncertainty_crop_fig_collection = []
cumul_crop_fig_collection = []
cumul_prop_fig_collection = []

# %% Import data
for ISO in filt_ISOs
    # Survey Data
    survey_net_crop = JLD2.load(survey_data_dir*"$(ISO)_2000_2023_cropextract.jld2")["NET_CROP_MONTHLY"]
    survey_months = JLD2.load(survey_data_dir*"$(ISO)_2000_2023_cropextract.jld2")["MONTHS_MONTHLY"]

    # Predictions
    netcrop_preds = JLD2.load(netcrop_dir*"$(ISO)_netcrop_prediction.jld2")
    Γ_BYNET_pred_CI = netcrop_preds["Γ_BYNET_pred"]
    Γ_TOTAL_pred_CI = netcrop_preds["Γ_TOTAL_pred"]
    NET_NAMES = netcrop_preds["NET_NAMES"]
    n_net_types = length(NET_NAMES)

    # %% Make plot of volumes with uncertainty
    fig = plot(xticks = (1:12:n_months, YEAR_START:YEAR_END), 
                xlabel = "Year", ylabel = "Net Crop (mil)",
                xtickfontrotation = 90, legend = :topleft,
                title = "Predicted Net Crop ($ISO)")
    plot!(fig, Γ_TOTAL_pred_CI[:,1]./1e6, fillrange = Γ_TOTAL_pred_CI[:,3]./1e6, 
                fillcolor = total_net_col, fillalpha = 0.2,
                linealpha = 0, label = nothing)
    plot!(fig, Γ_TOTAL_pred_CI[:,2]./1e6,
                linecolor = total_net_col, linealpha = 1,
                label = "Total Nets", linewidth = 1.2)

    for j in 1:n_net_types
        plot!(fig, Γ_BYNET_pred_CI[:,1,j]./1e6, fillrange = Γ_BYNET_pred_CI[:,3,j]./1e6, 
                fillcolor = colors[j], fillalpha = 0.2*component_alpha_mult,
                linealpha = 0, label = nothing)
        plot!(fig, Γ_BYNET_pred_CI[:,2,j]./1e6,
                linecolor = colors[j], linealpha = 1*component_alpha_mult,
                label = NET_NAMES[j], linewidth = 1.2)
    end

    # %% Make cumulative volumes plot
    cumul_data = zeros(size(Γ_BYNET_pred_CI)[1], size(Γ_BYNET_pred_CI)[3]+1)
    for i in 1:size(Γ_BYNET_pred_CI)[3]
        cumul_data[:, i+1] = sum(Γ_BYNET_pred_CI[:,2,1:i], dims = 2)
    end


    cumul_fig = plot(xticks = (1:12:n_months, YEAR_START:YEAR_END), 
                xlabel = "Year", ylabel = "Net Crop (mil)",
                xtickfontrotation = 90, legend = :topleft,
                title = "Predicted Net Crop ($ISO)")

    for i in 1:length(NET_NAMES)
        plot!(cumul_fig, cumul_data[:,i]./1e6, fillrange = cumul_data[:,i+1]./1e6,
                fillcolor = colors[i], fillalpha = component_alpha_mult,
                linealpha = 0, label = NET_NAMES[i])
    end

    # %% Make makeup proportions plot
    Γ_BYNET_proportions = zeros(size(Γ_BYNET_pred_CI)[1], size(Γ_BYNET_pred_CI)[3])

    for i in 1:n_net_types
        Γ_BYNET_proportions[:,i] = Γ_BYNET_pred_CI[:,2,i]./sum(Γ_BYNET_pred_CI[:,2,:], dims = 2)
    end

    cumul_proportions = zeros(size(Γ_BYNET_pred_CI)[1], size(Γ_BYNET_pred_CI)[3]+1)
    for i in 1:size(Γ_BYNET_pred_CI)[3]
        cumul_proportions[:, i+1] = sum(Γ_BYNET_proportions[:,1:i], dims = 2)
    end

    cumul_prop_fig = plot(xticks = (1:12:n_months, YEAR_START:YEAR_END), 
                xlabel = "Year", ylabel = "Proportion of Active Nets",
                xtickfontrotation = 90, legend = :topleft,
                title = "Predicted Net Type Breakdown ($ISO)")

    for i in 1:length(NET_NAMES)
        plot!(cumul_prop_fig, cumul_proportions[:,i], fillrange = cumul_proportions[:,i+1],
                fillcolor = colors[i], fillalpha = component_alpha_mult,
                linealpha = 0, label = NET_NAMES[i])
    end
    
    # Store plots
    push!(uncertainty_crop_fig_collection, fig)
    push!(cumul_crop_fig_collection, cumul_fig)
    push!(cumul_prop_fig_collection, cumul_prop_fig)
end

# %% Make combined plots and save

# Plot layout settings
layout = (8,6)
figsize = (2560,1800)

# Make plots and save
fig_uncertainty_crop_comb = plot(uncertainty_crop_fig_collection..., layout = layout, size = figsize)
fig_cumul_crop_comb = plot(cumul_crop_fig_collection..., layout = layout, size = figsize)
fig_cumul_prop_comb = plot(cumul_prop_fig_collection..., layout = layout, size = figsize)

savefig(fig_uncertainty_crop_comb, output_dir*"netcrop_type_breakdown_uncertainty.pdf")
savefig(fig_cumul_crop_comb, output_dir*"netcrop_type_breakdown.pdf")
savefig(fig_cumul_prop_comb, output_dir*"netcrop_cumul_proportions.pdf")

# %% Country Level Plots

# Extract net crop by type for each country
country_BYNET_crop = zeros(length(filt_ISOs),n_months, 3,n_net_types)
for ISO_i in 1:length(filt_ISOs)
    ISO = filt_ISOs[ISO_i]

    # Survey Data
    survey_net_crop = JLD2.load(survey_data_dir*"$(ISO)_2000_2023_cropextract.jld2")["NET_CROP_MONTHLY"]
    survey_months = JLD2.load(survey_data_dir*"$(ISO)_2000_2023_cropextract.jld2")["MONTHS_MONTHLY"]

    # Predictions
    netcrop_preds = JLD2.load(netcrop_dir*"$(ISO)_netcrop_prediction.jld2")
    country_BYNET_crop[ISO_i,:,:,:] = netcrop_preds["Γ_BYNET_pred"]
end

# Tally up total net crop by type across all countries
africa_BYNET_crop = sum(country_BYNET_crop, dims = 1)[1,:,:,:]

# Make plots
cumul_data = zeros(size(africa_BYNET_crop)[1], size(africa_BYNET_crop)[3]+1)
for i in 1:size(africa_BYNET_crop)[3]
    cumul_data[:, i+1] = sum(africa_BYNET_crop[:,2,1:i], dims = 2)
end


cumul_fig = plot(xticks = (1:12:n_months, YEAR_NAT_START:YEAR_NAT_END), 
            xlabel = "Year", ylabel = "Net Crop (mil)",
            xtickfontrotation = 90, legend = :topleft,
            title = "Predicted Net Crop (Africa)")

for i in 1:length(NET_NAMES)
    plot!(cumul_fig, cumul_data[:,i]./1e6, fillrange = cumul_data[:,i+1]./1e6,
            fillcolor = colors[i], fillalpha = component_alpha_mult,
            linealpha = 0, label = NET_NAMES[i])
end

cumul_fig
# %% Make makeup proportions plot
Γ_BYNET_proportions = zeros(size(africa_BYNET_crop)[1], size(africa_BYNET_crop)[3])

for i in 1:n_net_types
    Γ_BYNET_proportions[:,i] = africa_BYNET_crop[:,2,i]./sum(africa_BYNET_crop[:,2,:], dims = 2)
end

cumul_proportions = zeros(size(africa_BYNET_crop)[1], size(africa_BYNET_crop)[3]+1)
for i in 1:size(Γ_BYNET_pred_CI)[3]
    cumul_proportions[:, i+1] = sum(Γ_BYNET_proportions[:,1:i], dims = 2)
end

cumul_prop_fig = plot(xticks = (1:12:n_months, YEAR_START:YEAR_END), 
            xlabel = "Year", ylabel = "Proportion of Active Nets",
            xtickfontrotation = 90, legend = :topleft,
            title = "Predicted Net Type Breakdown (Africa)")

for i in 1:length(NET_NAMES)
    plot!(cumul_prop_fig, cumul_proportions[:,i], fillrange = cumul_proportions[:,i+1],
            fillcolor = colors[i], fillalpha = component_alpha_mult,
            linealpha = 0, label = NET_NAMES[i])
end

cumul_prop_fig

# Save plots
savefig(cumul_fig, output_dir*"africa_netcrop_type.pdf")
savefig(cumul_prop_fig, output_dir*"africa_cumul_proportions.pdf")