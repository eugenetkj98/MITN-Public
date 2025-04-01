"""
Author: Eugene Tan
Date Created: 4/11/2024
Last Updated: 11/11/2024
Script for plotting subnational NPC and access gap across regions
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %%
using GLMakie
using GeoIO
using GeoStats
using GeoInterface
using ProgressBars
using CSV
using DataFrames
using Missings
using JLD2
using DateConversions

# %%
ISO = "GHA"
YEAR_START = 2000
YEAR_END = 2023

# %% Import datasets
# Map geometries
geometries_dataset = GeoIO.load(raw"Z:\master_geometries\Admin_Units\Global\MAP\2023\MG_5K\admin2023_1_MG_5K.shp")
# National Estimates
nat_dataset = JLD2.load("outputs/draws/national/crop_access/$(ISO)_$(YEAR_START)_$(YEAR_END)_post_crop_access.jld2")
# Subnational Estimates
subnat_dataset = JLD2.load("outputs/draws/subnational/$(ISO)_SUBNAT_draws.jld2")

# %% Get number of regions
n_admin1 = length(subnat_dataset["admin1_names"])

# %% Extract gap values
NAT_A_NPC_mean_BYNET = nat_dataset["A_NPC_mean_BYNET"]
NAT_NPC_mean = sum(NAT_A_NPC_mean_BYNET, dims = (2,3))[:]
NAT_λ_access_mean = nat_dataset["λ_access_mean"]

SUBNAT_NPC_mean = zeros(n_admin1, length(NAT_NPC_mean))
SUBNAT_λ_access_mean = zeros(n_admin1, length(NAT_λ_access_mean))
for i in 1:n_admin1
    # SUBNAT_NPC_mean[i,:] = mean(subnat_dataset["admin1_outputs"][i]["NPC_MONTHLY_TOTAL_samples"], dims = 1)[:]
    # SUBNAT_λ_access_mean[i,:] = mean(subnat_dataset["admin1_outputs"][i]["λ_ACCESS_samples"], dims = 1)[:]

    SUBNAT_NPC_mean[i,:] = mean(subnat_dataset["merged_outputs"][i]["ADJ_NPC_MONTHLY_TOTAL_samples"], dims = 1)[:]
    SUBNAT_λ_access_mean[i,:] = mean(subnat_dataset["merged_outputs"][i]["ADJ_λ_ACCESS_samples"], dims = 1)[:]
end

SUBNAT_NPC_GAP = SUBNAT_NPC_mean .- repeat(NAT_NPC_mean, 1,n_admin1)'
SUBNAT_ACCESS_GAP = SUBNAT_λ_access_mean .- repeat(NAT_λ_access_mean, 1,n_admin1)'


# %% Make MAP plots coloured by gap

# Update functions
update_color(idx, data) = data[idx]

function monthyear(monthidx; YEAR_START = YEAR_START)
    month_ref,year_ref =  monthidx_to_monthyear(monthidx)
    month_val = month_ref
    year_val = year_ref + YEAR_START - 1
    return "$(year_val) - $(month_val)"
end    

# Plot options
cmap = :bam
NPC_clims = (-0.5,0.5)
ACCESS_clims = (-0.5,0.5)
cbarsize = 30
cbarfontsize = 20

# Initialise plot
fig = Figure(resolution = (1400,800))

# Calculate Starting year when transitioning to subnational
# monthstart = (subnat_dataset["YEAR_START"]-YEAR_START)*12+1
monthstart = 1

# Create sliders for time index
sg = SliderGrid(
    fig[2,1:4],
    (label = "Months", 
    range = monthstart:size(NAT_NPC_mean)[1],
    startvalue = monthstart)
)
month_index = sg.sliders[1].value

# Add axis
ax1 = fig[1,1] = Axis(fig, 
                        title = lift(idx -> "NPC Gap\n"*monthyear(idx), month_index),
                        titlesize = 30)
Colorbar(fig[1,2], limits = NPC_clims, colormap = cmap,
                    size = cbarsize, ticklabelsize = cbarfontsize)
ax2 = fig[1,3] = Axis(fig,
                        title = lift(idx -> "Access Gap\n"*monthyear(idx), month_index),
                        titlesize = 30)
Colorbar(fig[1,4], limits = ACCESS_clims, colormap = cmap,
                        size = cbarsize, ticklabelsize = cbarfontsize)

# Plot maps with adaptive colour
for i in 1:n_admin1
    admin1_name = subnat_dataset["admin1_names"][i] # Need to use area_id in the future after code has been fixed!!!
    admin1_geometry = geometries_dataset[(geometries_dataset.ISO .== ISO) .& 
                                            (geometries_dataset.Name_1 .== admin1_name),:][1,:].geometry
    plot!(ax1, Proj(Mercator)(admin1_geometry), 
                color = lift(update_color, month_index, SUBNAT_NPC_GAP[i,:]), 
                colormap = :bam,
                colorrange = (-0.5,0.5))
    plot!(ax2, Proj(Mercator)(admin1_geometry), 
                color = lift(update_color, month_index, SUBNAT_ACCESS_GAP[i,:]), 
                colormap = :bam,
                colorrange = (-0.5,0.5))
end
hidedecorations!(ax1)
hidedecorations!(ax2)

# Add Supertitle for country
supertitle = Label(fig[0,:], ISO, fontsize = 40)
fig

 