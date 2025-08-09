"""
Author: Eugene Tan
Date Created: 27/3/2025
Last Updated: 27/3/2025
Quick and dirty code to check and troubleshoot raster plots.
CURRENTLY IN USE TO TROUBLESHOOT RASTER OUTPUTS
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


monthly_mean = []

for i in 1:12
    if i < 10
        snf_rasters = replace_missing.(Raster.(OUTPUT_RASTERS_DIR.*"final_npc/snf_npc/npc_$(year)_0$(i)_sample_".*string.(1:10:100).*".tif"), missingval = NaN)
        mean_npc = mean(snf_rasters)
        push!(monthly_mean, mean_npc)
    else
        snf_rasters = replace_missing.(Raster.(OUTPUT_RASTERS_DIR.*"final_npc/snf_npc/npc_$(year)_$(i)_sample_".*string.(1:10:100).*".tif"), missingval = NaN)
        mean_npc = mean(snf_rasters)
        push!(monthly_mean, mean_npc)
    end
end

annual_mean = mean(monthly_mean)


# %%
test_snf = mean([Raster.(OUTPUT_RASTERS_DIR.*"final_npc/snf_npc/npc_$(year)_06_sample_$(i)".*".tif") for i in 1:10])
test_logmodel = mean([Raster.(OUTPUT_RASTERS_DIR.*"inla_logmodel_npc/NPC_logmodel_$(year)_sample_$(i)".*".tif") for i in 1:10])
test_logmodel_nospde = mean([Raster.(OUTPUT_RASTERS_DIR.*"inla_logmodel_npc_nospde/NPC_logmodel_$(year)_sample_$(i)".*".tif") for i in 1:10])


fig2 = Figure(size = (1600,800))
ax_11 = Axis(fig2[1,1], title = "SPDE Deviation Ratio", titlesize = 22)
ax_12 = Axis(fig2[1,2], title = "NO SPDE Deviation Ratio", titlesize = 22)
ax_21 = Axis(fig2[2,1], title = "SPDE NPC", titlesize = 22)
ax_22 = Axis(fig2[2,2], title = "NO SPDE NPC", titlesize = 22)

plot!(ax_11, exp.(test_logmodel), colorrange = (0,2))
plot!(ax_12, exp.(test_logmodel_nospde), colorrange = (0,2))
plot!(ax_21, exp.(test_logmodel).*test_snf, colorrange = (0,1))
plot!(ax_22, exp.(test_logmodel_nospde).*test_snf, colorrange = (0,1))

Colorbar(fig2[1,3], colorrange = (0,2), label = "Deviation Ratio", labelsize = 20)
Colorbar(fig2[2,3], colorrange = (0,1), label = "NPC", labelsize = 20)

ax_snf = Axis(fig2[1:2,4:5], title = "SNF Subnat", titlesize = 22)
plot!(ax_snf, test_snf, colorrange = (0,1))
Colorbar(fig2[1:2,6], colorrange = (0,1), label = "NPC", labelsize = 20)

fig2

save(OUTPUT_PLOTS_DIR*"Special_Request/2020_spatial_comparisons.png",fig2)
# %% Posteriors
year = 2023

ras_monthly_mean = []

for i in 1:12
    if i < 10
        snf_rasters = replace_missing.(Raster.(OUTPUT_RASTERS_DIR.*"final_npc/posterior_samples/npc_$(year)_0$(i)_sample_".*string.(1:10:100).*".tif"), missingval = NaN)
        mean_npc = mean(snf_rasters)
        push!(ras_monthly_mean, mean_npc)
    else
        snf_rasters = replace_missing.(Raster.(OUTPUT_RASTERS_DIR.*"final_npc/posterior_samples/npc_$(year)_$(i)_sample_".*string.(1:10:100).*".tif"), missingval = NaN)
        mean_npc = mean(snf_rasters)
        push!(ras_monthly_mean, mean_npc)
    end
end

ras_annual_mean = mean(ras_monthly_mean)

# %% INLA deviations

logmodel_rasters = replace_missing.(Raster.(OUTPUT_RASTERS_DIR.*"inla_logmodel_npc/NPC_logmodel_$(year)_sample_".*string.(5:5:100).*".tif"), missingval = NaN)
mean_logmodel = mean(logmodel_rasters)
fig = Figure()
plot(fig[1,1],mean_logmodel, colorrange = (0,2))
Colorbar(fig[1,2], limits = (0,2), size = cbar_size, 
        ticklabelsize = 15, ticks = 0:0.2:2)
fig

# %% Manually calculate posterior
year = 2023
n_batch = 20
n_sample_batch = 40
post_monthly = []

for i in ProgressBar(1:12, leave = false)
    if i < 10
        temp = []
        # for batch in ProgressBar(1:n_batch, leave = false)
        for batch in 1:n_batch
            logmodel_raster_samples = []
            j = rand(1:100)
            snf_raster = replace_missing.(Raster.(OUTPUT_RASTERS_DIR.*"final_npc/posterior_samples/npc_$(year)_0$(i)_sample_$(j).tif"), missingval = NaN)

            # for sample in ProgressBar(1:n_sample_batch, leave = false)
            for sample in 1:n_sample_batch
                println("month$(i), batch $(batch)/$(n_batch), sample $(sample)/$(n_sample_batch)")
                k = rand(1:100)
                logmodel_rasters = replace_missing.(Raster.(OUTPUT_RASTERS_DIR.*"inla_logmodel_npc/NPC_logmodel_$(year)_sample_$(k).tif"), missingval = NaN)
                # npc_raster = (exp.(logmodel_rasters)).*(snf_raster)
                push!(logmodel_raster_samples, logmodel_rasters)
            end
            logmodel_raster_sample_mean = mean(logmodel_raster_samples)
            
            npc_raster_sample_mean = snf_raster .* exp.(logmodel_raster_sample_mean)
            push!(temp, npc_raster_sample_mean)
        end
        push!(post_monthly, mean(temp))
    else
        temp = []
        # for batch in ProgressBar(1:n_batch, leave = false)
        for batch in 1:n_batch
            logmodel_raster_samples = []
            j = rand(1:100)
            snf_raster = replace_missing.(Raster.(OUTPUT_RASTERS_DIR.*"final_npc/posterior_samples/npc_$(year)_$(i)_sample_$(j).tif"), missingval = NaN)
            
            # for sample in ProgressBar(1:n_sample_batch, leave = false)
            for sample in 1:n_sample_batch
                println("month$(i), batch $(batch)/$(n_batch), sample $(sample)/$(n_sample_batch)")
                k = rand(1:100)
                logmodel_rasters = replace_missing.(Raster.(OUTPUT_RASTERS_DIR.*"inla_logmodel_npc/NPC_logmodel_$(year)_sample_$(k).tif"), missingval = NaN)
                # npc_raster = (exp.(logmodel_rasters)).*(snf_raster)
                push!(logmodel_raster_samples, logmodel_rasters)
            end
            logmodel_raster_sample_mean = mean(logmodel_raster_samples)

            npc_raster_sample_mean = snf_raster .* exp.(logmodel_raster_sample_mean)
            push!(temp, npc_raster_sample_mean)
        end
        push!(post_monthly, mean(temp))
    end
end

post_annual = mean(post_monthly)


# %%
fig = Figure()
plot(fig[1,1],annual_mean, colorrange = (0,1))
Colorbar(fig[1,2], limits = (0,1), size = cbar_size, 
        ticklabelsize = 15, ticks = 0:0.2:1)
fig

# %%
test = exp.(mean_logmodel) .* annual_mean
fig = Figure()
plot(fig[1,1],test, colorrange = (0,1))
Colorbar(fig[1,2], limits = (0,1), size = cbar_size, 
        ticklabelsize = 15, ticks = 0:0.2:1)
fig

# %%
fig = Figure()
plot(fig[1,1],ras_annual_mean, colorrange = (0,1))
Colorbar(fig[1,2], limits = (0,1), size = cbar_size, 
        ticklabelsize = 15, ticks = 0:0.2:1)
fig

# %%
fig = Figure()
plot(fig[1,1],post_annual, colorrange = (0,1))
Colorbar(fig[1,2], limits = (0,1), size = cbar_size, 
        ticklabelsize = 15, ticks = 0:0.2:1)
fig

# %%
fig = Figure()
plot(fig[1,1],post_annual, colorrange = (0,1))
Colorbar(fig[1,2], limits = (0,1), size = cbar_size, 
        ticklabelsize = 15, ticks = 0:0.2:1)
fig





# %%
fig = Figure()
plot(fig[1,1],exp.(mean_logmodel), colorrange = (0,2))
Colorbar(fig[1,2], limits = (0,2), size = cbar_size, 
        ticklabelsize = 15, ticks = 0:0.2:2)
fig

# %%
fig = Figure()
plot(fig[1,1],exp.(mean_logmodel).>1, colorrange = (0,1))
Colorbar(fig[1,2], limits = (0,1), size = cbar_size, 
        ticklabelsize = 15, ticks = 0:0.2:2)
fig



