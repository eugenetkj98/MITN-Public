"""
Author: Eugene Tan
Date Created: 22/7/2025
Last Updated: 22/7/2025
Plot NPC maps for all years quick and dirty
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

# %% Just grab snf output
year = 2020

for year in ProgressBar(2000:2023)

    # %% Construct NPC raster
    test_snf = mean([Raster.(OUTPUT_RASTERS_DIR.*"final_npc/snf_npc/npc_$(year)_06_sample_$(i)".*".tif") for i in 1:10])
    test_logmodel_ratio = mean([Raster.(OUTPUT_RASTERS_DIR.*"inla_logmodel_npc/NPC_RATIO_logmodel_$(year)_sample_$(i)".*".tif") for i in 1:10])
    test_logmodel_residual = mean([Raster.(OUTPUT_RASTERS_DIR.*"inla_logmodel_npc/NPC_RESIDUAL_logmodel_$(year)_sample_$(i)".*".tif") for i in 1:10])

    mean_npc_raster = test_snf.*test_logmodel_ratio .+ test_logmodel_residual

    # %%
    fig = Figure(size = (950,900))

    ax11 = Axis(fig[1,1], title = "Deviation Ratio", titlesize = 22)
    ax12 = Axis(fig[1,3], title = "Residuals Field", titlesize = 22)
    ax21 = Axis(fig[2,1], title = "SNF NPC", titlesize = 22)
    ax22 = Axis(fig[2,3], title = "Raster NPC", titlesize = 22)

    plot!(ax11, test_logmodel_ratio, colorrange = (0,1.5))
    plot!(ax12, test_logmodel_residual, colorrange = (-0.5,0.5))
    plot!(ax21, test_snf, colorrange = (0,1))
    plot!(ax22, mean_npc_raster, colorrange = (0,1))

    Colorbar(fig[1,2], colorrange = (0,1.5))
    Colorbar(fig[1,4], colorrange = (-0.5,0.5))
    Colorbar(fig[2,2], colorrange = (0,1))
    Colorbar(fig[2,4], colorrange = (0,1))

    Label(fig[0,:], "YEAR $(year)", fontsize = 30)
    

    # %% Save figure
    save(OUTPUT_PLOTS_DIR*"Special_Request/"*"Raster_Checks/NPC_$(year).png", fig)
end