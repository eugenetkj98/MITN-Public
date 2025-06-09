"""
Author: Eugene Tan
Date Created: 27/3/2025
Last Updated: 27/3/2025
Make continent level plots for paper
"""
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/read_toml.jl")

# %% Import Public Packages
using DataFrames
using JLD2
using CSV
using ProgressBars
using GeoIO
using Rasters

# %% Maths packages
using LinearAlgebra
using StatsBase

# %% General useful functions
using DateConversions
using TS_filters

# %% Plot packages
using LaTeXStrings
using CairoMakie

# %% Plotting Theme and general settings
set_theme!(theme_ggplot2())

# Color settings
colors = [  colorant"#0082C7", # NPC
            colorant"#E72A3D", # Access
            colorant"#538255", # Use
            colorant"#45332C", # Guidelines
            ];

# General settings
fillalpha = 0.2
la = 0.5
lw = 2
label_fontsize = 25

# Plot color and limit settings
netage_clims = (0,3)
clims = (0,1)
cbar_size = 20
netage_cmap = :roma
npc_cmap = :navia
access_cmap = :lajolla
use_cmap = :devon
util_cmap = :lipari

# %% Choose years to plot snapshots of
year_vals = [2010, 2015, 2020, 2023]

# Define raster filenames and import
netage_raster_filenames = ["netage_$(year)_mean.tif" for year in year_vals]
npc_raster_filenames = ["npc_$(year)_mean.tif" for year in year_vals]
access_raster_filenames = ["access_$(year)_mean.tif" for year in year_vals]
use_raster_filenames = ["use_$(year)_mean.tif" for year in year_vals]
util_raster_filenames = ["utilisation_$(year)_mean.tif" for year in year_vals]

netage_rasters = replace_missing.(Raster.(OUTPUT_RASTERS_DIR.*"/final_netage/annual/".*netage_raster_filenames), missingval = NaN)
npc_rasters = replace_missing.(Raster.(OUTPUT_RASTERS_DIR.*"final_npc/mean/annual/".*npc_raster_filenames), missingval = NaN)
access_rasters = replace_missing.(Raster.(OUTPUT_RASTERS_DIR.*"final_access/mean/annual/".*access_raster_filenames), missingval = NaN)
use_rasters = replace_missing.(Raster.(OUTPUT_RASTERS_DIR.*"final_use/mean/annual/".*use_raster_filenames), missingval = NaN)
util_rasters = replace_missing.(Raster.(OUTPUT_RASTERS_DIR.*"final_utilisation/mean/annual/".*util_raster_filenames), missingval = NaN)

# %% Construct snapshot plot
# Make Figure and define axes
layout_resolution = (1000,1100)
fig = Figure(size = layout_resolution)

plot_axs = Matrix{Any}(undef, 5, length(year_vals))
for i in 1:length(year_vals)
    # Raster axes
    for j in 1:5
        ax = Axis(fig[j,i])
        hidedecorations!(ax)
        plot_axs[j,i] = ax
    end
end
# Colorbar Axes: Four metrics: NPC, Access, Use, Utilisation
Colorbar(fig[1,5], limits = netage_clims, colormap = netage_cmap, size = cbar_size, 
        ticklabelsize = 15, ticks = 0:0.5:3)
Colorbar(fig[2,5], limits = clims, colormap = npc_cmap, size = cbar_size, 
        ticklabelsize = 15, ticks = 0:0.2:1)
Colorbar(fig[3,5], limits = clims, colormap = access_cmap, size = cbar_size, 
        ticklabelsize = 15, ticks = 0:0.2:1)
Colorbar(fig[4,5], limits = clims, colormap = use_cmap, size = cbar_size, 
        ticklabelsize = 15, ticks = 0:0.2:1)
Colorbar(fig[5,5], limits = clims, colormap = util_cmap, size = cbar_size, 
        ticklabelsize = 15, ticks = 0:0.2:1)

# Add Text labels
for i in 1:length(year_vals)
    Label(fig[0,i], "$(year_vals[i])", fontsize = label_fontsize,
            padding = (0,0,-10,-10,), tellwidth = false)
end
Label(fig[1,6], "Mean Net Age\n(Years)", fontsize = label_fontsize, rotation = pi/2,
        tellheight = false,)
Label(fig[2,6], "NPC", fontsize = label_fontsize, rotation = pi/2,
        tellheight = false,)
Label(fig[3,6], "Access", fontsize = label_fontsize, rotation = pi/2,
        tellheight = false,)
Label(fig[4,6], "Use", fontsize = label_fontsize, rotation = pi/2,
        tellheight = false,)
Label(fig[5,6], "Utilisation", fontsize = label_fontsize, rotation = pi/2,
        tellheight = false,)

# Plot NPC rasters
for i in 1:length(year_vals)
    plot!(plot_axs[1,i], netage_rasters[i]./12, colormap = netage_cmap, colorrange = netage_clims)
    plot!(plot_axs[2,i], npc_rasters[i], colormap = npc_cmap, colorrange = clims)
    plot!(plot_axs[3,i], access_rasters[i], colormap = access_cmap, colorrange = clims)
    plot!(plot_axs[4,i], use_rasters[i], colormap = use_cmap, colorrange = clims)
    plot!(plot_axs[5,i], use_rasters[i]./(2 .*npc_rasters[i]), colormap = util_cmap, colorrange = clims)
end
fig
# %% Save figure
# save(OUTPUT_PLOTS_DIR*"PaperFigures/Raster_snapshots.pdf", fig, pdf_version = "1.4")
save(OUTPUT_PLOTS_DIR*"PaperFigures/Raster_snapshots.png", fig)
fig