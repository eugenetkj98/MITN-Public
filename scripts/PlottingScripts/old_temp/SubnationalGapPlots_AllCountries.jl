"""
Author: Eugene Tan
Date Created: 4/11/2024
Last Updated: 11/11/2024
Script for plotting subnational NPC and access gap across all African Regions
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %%
# using GLMakie
using GeoIO
using GeoStats
using GeoInterface
using ProgressBars
using CSV
using DataFrames
using Missings
using JLD2
using DateConversions
using ProgressBars
using GLMakie

# 
using Theil

# %% Import datasets
# Map geometries
nat_geometries_dataset = GeoIO.load(raw"Z:\master_geometries\Admin_Units\Global\MAP\2023\MG_5K\admin2023_0_MG_5K.shp")
subnat_geometries_dataset = GeoIO.load(raw"Z:\master_geometries\Admin_Units\Global\MAP\2023\MG_5K\admin2023_1_MG_5K.shp")

# BV NPC Estimates
BV_dataset = CSV.read("datasets/subnational/BV_inla_npc_estimates.csv", DataFrame)
# %% Get ISO List

# Perform draws and save outputs. Filter out unwanted countries
ISO_list = String.(CSV.read(raw"C:\Users\ETan\Documents\Prototype Analyses\itn-updated\datasets\ISO_list.csv", DataFrame)[:,1])
exclusion_ISOs = ["CPV","ZAF"]#["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Make Storage variables
NAT_NPC_collection = Vector{Any}(undef, length(filt_ISOs))
NAT_ACCESS_collection = Vector{Any}(undef, length(filt_ISOs))
SUBNAT_NPC_collection = Vector{Any}(undef, length(filt_ISOs))
SUBNAT_ACCESS_collection = Vector{Any}(undef, length(filt_ISOs))
NPC_GAP_collection = Vector{Any}(undef, length(filt_ISOs))
NPC_GAP_normalised_collection = Vector{Any}(undef, length(filt_ISOs))
ACCESS_GAP_collection = Vector{Any}(undef, length(filt_ISOs))
ACCESS_GAP_normalised_collection = Vector{Any}(undef, length(filt_ISOs))
# NPC_VAR_collection = Vector{Any}(undef, length(filt_ISOs))
ACCESS_STD_collection = Vector{Any}(undef, length(filt_ISOs))
NAT_geometries_collection = Vector{Any}(undef, length(filt_ISOs))
SUBNAT_geometries_collection = Vector{Any}(undef, length(filt_ISOs))
THEIL_collection =  Vector{Any}(undef, length(filt_ISOs))

BV_SUBNAT_NPC_collection = Vector{Any}(undef, length(filt_ISOs))
BV_MITN_SUBNAT_NPC_GAP_collection = Vector{Any}(undef, length(filt_ISOs))
BV_MITN_SUBNAT_NPC_GAP_normalised_collection = Vector{Any}(undef, length(filt_ISOs))

# %%
for ISO_i in ProgressBar(1:length(filt_ISOs))

    # %%
    ISO = filt_ISOs[ISO_i]
    YEAR_START = 2000
    YEAR_END = 2023

    # %% Load Draw data
    # National Estimates
    nat_dataset = JLD2.load("outputs/draws/national/crop_access/$(ISO)_$(YEAR_START)_$(YEAR_END)_post_crop_access.jld2")
    # Subnational Estimates
    subnat_dataset = JLD2.load("outputs/draws/subnational/$(ISO)_SUBNAT_draws.jld2")

    # %% Get number of regions
    n_admin1 = length(subnat_dataset["admin1_names"])

    # %% Get population data
    subnat_reg_dir = "outputs/regressions/subnational/"
    REG_YEAR_START = 2010 # YEAR CHOSEN FOR START OF SUBNATIONAL REGRESSION WHEN DOING MCMC DRAWS
    subnat_reg_filename = "$(ISO)_SUBNAT_NETCROP_$(YEAR_START)_$(REG_YEAR_START)_$(YEAR_END)_regression.jld2"
    subnat_reg_data = JLD2.load(subnat_reg_dir*subnat_reg_filename)

    

    # Extract population values
    # Get month index w.r.t to trained indices
    TRAINED_YEAR_START = subnat_dataset["YEAR_START_NAT"]
    year_idx_ref_1 = YEAR_START-TRAINED_YEAR_START+1
    year_idx_ref_2 = YEAR_END-TRAINED_YEAR_START+1

    monthidx_ref_1 = (year_idx_ref_1-1)*12 + 1
    monthidx_ref_2 = (year_idx_ref_2-1)*12 + 12

    subnat_monthidxs = monthidx_ref_1:monthidx_ref_2
    n_months = length(subnat_monthidxs)

    # Extract admin1 populations
    FULL_POPULATION_MONTHLY = zeros(n_admin1, n_months)

    for i in 1:n_admin1
        FULL_POPULATION_MONTHLY[i,:] = subnat_reg_data["admin1_outputs"][i]["FULL_POPULATION_MONTHLY"][subnat_monthidxs]
    end

    SUBNAT_POPULATIONS = zeros(n_admin1, n_months)
    for i in 1:n_admin1
        SUBNAT_POPULATIONS[i,:] = subnat_reg_data["admin1_outputs"][i]["FULL_POPULATION_MONTHLY"]
    end
    
    # Extracted adjusted subnational means and calculate gap

    SUBNAT_NPC_mean = zeros(n_admin1, n_months)
    SUBNAT_λ_access_mean = zeros(n_admin1, n_months)

    for i in 1:n_admin1
        # SUBNAT_NPC_mean[i,:] = mean(subnat_dataset["admin1_outputs"][i]["NPC_MONTHLY_TOTAL_samples"], dims = 1)[:]
        # SUBNAT_λ_access_mean[i,:] = mean(subnat_dataset["admin1_outputs"][i]["λ_ACCESS_samples"], dims = 1)[:]

        SUBNAT_NPC_mean[i,:] = mean(subnat_dataset["merged_outputs"][i]["ADJ_NPC_MONTHLY_TOTAL_samples"], dims = 1)[:]
        SUBNAT_λ_access_mean[i,:] = mean(subnat_dataset["merged_outputs"][i]["ADJ_λ_ACCESS_samples"], dims = 1)[:]
    end

    # %% Calculate National NPC using Subnational estimates
    NAT_NPC_mean = (sum(SUBNAT_NPC_mean .* FULL_POPULATION_MONTHLY, dims = 1)[:])./(sum(FULL_POPULATION_MONTHLY, dims = 1)[:])
    NAT_λ_access_mean = (sum(SUBNAT_λ_access_mean .* FULL_POPULATION_MONTHLY, dims = 1)[:])./(sum(FULL_POPULATION_MONTHLY, dims = 1)[:])
    
    # %% Extract gap values
    # NAT_A_NPC_mean_BYNET = nat_dataset["A_NPC_mean_BYNET"]
    # NAT_NPC_mean = sum(NAT_A_NPC_mean_BYNET, dims = (2,3))[:]
    # NAT_λ_access_mean = nat_dataset["λ_access_mean"]

    SUBNAT_NPC_GAP = SUBNAT_NPC_mean .- repeat(NAT_NPC_mean, 1,n_admin1)'
    SUBNAT_ACCESS_GAP = SUBNAT_λ_access_mean .- repeat(NAT_λ_access_mean, 1,n_admin1)'

    # Calculate normalised gap variation coefficient    
    NPC_GAP_mean = mean(SUBNAT_NPC_GAP, dims = 1)[:]
    ACCESS_GAP_mean = mean(SUBNAT_NPC_GAP, dims = 1)[:]

    # Calculate normalised Gap values
    SUBNAT_NPC_GAP_normalised = SUBNAT_NPC_GAP./repeat(NAT_NPC_mean', n_admin1, 1)
    SUBNAT_ACCESS_GAP_normalised = SUBNAT_ACCESS_GAP./repeat(NAT_λ_access_mean', n_admin1, 1)
    
   

    # Theil values
    THEIL_values = zeros(n_months)
    
    for month_i in 1:n_months
        THEIL_values[month_i] = NPC_to_theil(SUBNAT_NPC_mean[:,month_i], FULL_POPULATION_MONTHLY[:,month_i])
    end


    # NPC_VAR = (sqrt.(sum(((SUBNAT_NPC_GAP.- repeat(NPC_GAP_mean', n_admin1)).^2).*SUBNAT_POPULATIONS, dims = 1)./sum(SUBNAT_POPULATIONS, dims = 1))[:])
    ACCESS_STD = (sqrt.(sum(((SUBNAT_ACCESS_GAP.- repeat(ACCESS_GAP_mean', n_admin1)).^2).*SUBNAT_POPULATIONS, dims = 1)./sum(SUBNAT_POPULATIONS, dims = 1))[:])

    # %% Extract geometries
    # National
    NAT_geometries_collection[ISO_i] = nat_geometries_dataset[nat_geometries_dataset.ISO .== ISO,:][1,:].geometry

    # Subnational
    SUBNAT_geometries = Vector{Any}(undef, n_admin1)

    for i in 1:n_admin1
        area_id = subnat_dataset["merged_outputs"][i]["area_id"]
        admin1_geometry = subnat_geometries_dataset[subnat_geometries_dataset.area_id .== area_id,:][1,:].geometry
        SUBNAT_geometries[i] = admin1_geometry
    end

    # Get BV estimates for NPC

    
    # Storage variable
    BV_SUBNAT_NPC_mean = zeros(size(SUBNAT_NPC_mean))
    for i in 1:n_admin1
        area_id = subnat_dataset["merged_outputs"][i]["area_id"]
        # Import data from BV extraction
        BV_SUBNAT_NPC_ANNUAL = BV_dataset[BV_dataset.area_id .== area_id,:].NPC
        for year_i in 1:(size(BV_SUBNAT_NPC_ANNUAL)[1]-1)
            BV_SUBNAT_NPC_mean[i,12*(year_i-1)+1:12*year_i] = LinRange(BV_SUBNAT_NPC_ANNUAL[year_i],BV_SUBNAT_NPC_ANNUAL[year_i+1],13)[1:12]
        end
    end

    BV_MITN_SUBNAT_NPC_GAP = (SUBNAT_NPC_mean .- BV_SUBNAT_NPC_mean)
    BV_MITN_SUBNAT_NPC_GAP_normalised = BV_MITN_SUBNAT_NPC_GAP./BV_SUBNAT_NPC_mean

    # %% Store in variable
    NAT_NPC_collection[ISO_i] = NAT_NPC_mean
    NAT_ACCESS_collection[ISO_i] = NAT_λ_access_mean
    SUBNAT_NPC_collection[ISO_i] = SUBNAT_NPC_mean
    SUBNAT_ACCESS_collection[ISO_i] = SUBNAT_λ_access_mean
    NPC_GAP_collection[ISO_i] = SUBNAT_NPC_GAP
    NPC_GAP_normalised_collection[ISO_i] = SUBNAT_NPC_GAP_normalised
    ACCESS_GAP_collection[ISO_i] = SUBNAT_ACCESS_GAP
    ACCESS_GAP_normalised_collection[ISO_i] = SUBNAT_ACCESS_GAP_normalised
    THEIL_collection[ISO_i] = THEIL_values
    ACCESS_STD_collection[ISO_i] = ACCESS_STD
    SUBNAT_geometries_collection[ISO_i] = SUBNAT_geometries

    BV_SUBNAT_NPC_collection[ISO_i] = BV_SUBNAT_NPC_mean
    BV_MITN_SUBNAT_NPC_GAP_collection[ISO_i] = BV_MITN_SUBNAT_NPC_GAP
    BV_MITN_SUBNAT_NPC_GAP_normalised_collection[ISO_i] = BV_MITN_SUBNAT_NPC_GAP_normalised
end

# %% Make MAP plots coloured by gap

# Update functions
update_color(idx, data) = data[idx]
update_color_withalpha(idx, data, alpha) = (data[idx], )
YEAR_START = 2000
function monthyear(monthidx; YEAR_START = YEAR_START)
    month_ref,year_ref =  monthidx_to_monthyear(monthidx)
    month_val = month_ref
    year_val = year_ref + YEAR_START - 1
    return "$(year_val) - $(month_val)"
end   
# %%


# Plot options
npc_cmap = :devon
access_cmap = :bamako
theil_cmap = :amp
gap_cmap = :balance
norm_gap_cmap = :balance
NPC_clims = (0,1)
ACCESS_clims = (0,1)
THEIL_clims = (0.0,0.5)
ACCESS_STD_clims = (0,0.5)
NPC_GAP_clims = (-0.5,0.5)
ACCESS_GAP_clims = (-0.5,0.5)
NPC_GAP_normalised_clims = (-2,2)
ACCESS_GAP_normalised_clims = (-2,2)
cbarsize = 30
cbarfontsize = 20
alpha = 0.85
segmentsize = 0.5
# Initialise plot
fig = Figure(resolution = (1400,800))

# Calculate Starting year when transitioning to subnational
# monthstart = (subnat_dataset["YEAR_START"]-YEAR_START)*12+1
monthstart = 1

# Create sliders for time index
sg = SliderGrid(
    fig[4,1:8],
    (label = "Months", 
    range = monthstart:length(NAT_NPC_collection[1]),
    startvalue = monthstart)
)
month_index = sg.sliders[1].value

# Add Axis
ax1 = fig[1,1] = Axis(fig, 
                        title = lift(idx -> "National NPC\n"*monthyear(idx), month_index),
                        titlesize = 30)
Colorbar(fig[1,2], limits = NPC_clims, colormap = (npc_cmap, alpha),
                    size = cbarsize, ticklabelsize = cbarfontsize)
ax2 = fig[1,3] = Axis(fig, 
                    title = lift(idx -> "Subnational NPC\n"*monthyear(idx), month_index),
                    titlesize = 30)
Colorbar(fig[1,4], limits = NPC_clims, colormap = (npc_cmap, alpha),
                size = cbarsize, ticklabelsize = cbarfontsize)

ax3 = fig[2,1] = Axis(fig, 
                title = lift(idx -> "National Access\n"*monthyear(idx), month_index),
                titlesize = 30)
Colorbar(fig[2,2], limits = ACCESS_clims, colormap = (access_cmap, alpha),
            size = cbarsize, ticklabelsize = cbarfontsize)

ax4 = fig[2,3] = Axis(fig, 
            title = lift(idx -> "Subnational Access\n"*monthyear(idx), month_index),
            titlesize = 30)
Colorbar(fig[2,4], limits = ACCESS_clims, colormap = (access_cmap, alpha),
        size = cbarsize, ticklabelsize = cbarfontsize)

ax7 = fig[1,5] = Axis(fig,
        title = lift(idx -> "NPC Gap\n"*monthyear(idx), month_index),
        titlesize = 30)
Colorbar(fig[1,6], limits = NPC_GAP_clims, colormap = gap_cmap,
        size = cbarsize, ticklabelsize = cbarfontsize)
ax8 = fig[2,5] = Axis(fig, 
        title = lift(idx -> "Access Gap\n"*monthyear(idx), month_index),
        titlesize = 30)
Colorbar(fig[2,6], limits = ACCESS_GAP_clims, colormap = gap_cmap,
    size = cbarsize, ticklabelsize = cbarfontsize)

ax9 = fig[1,7] = Axis(fig,
    title = lift(idx -> "NPC Normalised Gap\n"*monthyear(idx), month_index),
    titlesize = 30)
Colorbar(fig[1,8], limits = NPC_GAP_normalised_clims, colormap = norm_gap_cmap,
    size = cbarsize, ticklabelsize = cbarfontsize)
ax10 = fig[2,7] = Axis(fig, 
    title = lift(idx -> "Access Normalised Gap\n"*monthyear(idx), month_index),
    titlesize = 30)
Colorbar(fig[2,8], limits = ACCESS_GAP_normalised_clims, colormap = norm_gap_cmap,
size = cbarsize, ticklabelsize = cbarfontsize)


        


ax5 = fig[3,1] = Axis(fig, 
        title = lift(idx -> "NPC Theil\n"*monthyear(idx), month_index),
        titlesize = 30)
Colorbar(fig[3,2], limits = THEIL_clims, colormap = (theil_cmap, alpha),
    size = cbarsize, ticklabelsize = cbarfontsize)
ax6 = fig[3,3] = Axis(fig, 
    title = lift(idx -> "Access Gap σ \n"*monthyear(idx), month_index),
    titlesize = 30)
Colorbar(fig[3,4], limits = ACCESS_STD_clims, colormap = (theil_cmap, alpha),
size = cbarsize, ticklabelsize = cbarfontsize)






ax11 = fig[3,5] = Axis(fig,
                title = lift(idx -> "BV-MITN NPC Gap\n"*monthyear(idx), month_index),
                titlesize = 30)
Colorbar(fig[3,6], limits = NPC_GAP_clims, colormap = gap_cmap,
                size = cbarsize, ticklabelsize = cbarfontsize)
ax12 = fig[3,7] = Axis(fig, 
                title = lift(idx -> "BV-MITN Normalised NPC Gap\n"*monthyear(idx), month_index),
                titlesize = 30)
Colorbar(fig[3,8], limits = NPC_GAP_normalised_clims, colormap = norm_gap_cmap,
            size = cbarsize, ticklabelsize = cbarfontsize)

# For each country
for ISO_i in 1:length(filt_ISOs)
    nat_geometry = NAT_geometries_collection[ISO_i]

    Makie.plot!(ax1, Proj(Mercator)(nat_geometry), 
                color = lift(update_color, month_index, NAT_NPC_collection[ISO_i]), 
                colormap = npc_cmap,
                colorrange = NPC_clims,
                alpha = alpha)
    Makie.plot!(ax3, Proj(Mercator)(nat_geometry), 
                color = lift(update_color, month_index, NAT_ACCESS_collection[ISO_i]), 
                colormap = access_cmap,
                colorrange = ACCESS_clims,
                alpha = alpha)
    Makie.plot!(ax5, Proj(Mercator)(nat_geometry), 
                color = lift(update_color, month_index, THEIL_collection[ISO_i]), 
                colormap = theil_cmap,
                colorrange = THEIL_clims,
                alpha = alpha)
    Makie.plot!(ax6, Proj(Mercator)(nat_geometry), 
                color = lift(update_color, month_index, ACCESS_STD_collection[ISO_i]), 
                colormap = theil_cmap,
                colorrange = ACCESS_clims,
                alpha = alpha)
                
    n_admin1 = size(SUBNAT_NPC_collection[ISO_i])[1]

    # Plot boundaries for country borders
    Makie.plot!(ax1, Proj(Mercator)(boundary(nat_geometry)), 
                        segmentsize = segmentsize, color = :black)
    Makie.plot!(ax3, Proj(Mercator)(boundary(nat_geometry)), 
                        segmentsize = segmentsize, color = :black)
    Makie.plot!(ax2, Proj(Mercator)(boundary(nat_geometry)), 
                        segmentsize = segmentsize, color = :black)

    Makie.plot!(ax4, Proj(Mercator)(boundary(nat_geometry)), 
                        segmentsize = segmentsize, color = :black)

    Makie.plot!(ax5, Proj(Mercator)(boundary(nat_geometry)), 
                        segmentsize = segmentsize, color = :black)

    Makie.plot!(ax6, Proj(Mercator)(boundary(nat_geometry)), 
                        segmentsize = segmentsize, color = :black)

    Makie.plot!(ax7, Proj(Mercator)(boundary(nat_geometry)), 
                        segmentsize = segmentsize, color = :black)

    Makie.plot!(ax8, Proj(Mercator)(boundary(nat_geometry)), 
                        segmentsize = segmentsize, color = :black)
                    
    Makie.plot!(ax9, Proj(Mercator)(boundary(nat_geometry)), 
                        segmentsize = segmentsize, color = :black)

    Makie.plot!(ax10, Proj(Mercator)(boundary(nat_geometry)), 
                        segmentsize = segmentsize, color = :black)

    Makie.plot!(ax11, Proj(Mercator)(boundary(nat_geometry)), 
                        segmentsize = segmentsize, color = :black)

    Makie.plot!(ax12, Proj(Mercator)(boundary(nat_geometry)), 
                        segmentsize = segmentsize, color = :black)

    

    for i in 1:n_admin1
        admin1_geometry = SUBNAT_geometries_collection[ISO_i][i]
        Makie.plot!(ax2, Proj(Mercator)(admin1_geometry), 
                        color = lift(update_color, month_index, SUBNAT_NPC_collection[ISO_i][i,:]), 
                        colormap = npc_cmap,
                        colorrange = NPC_clims,
                        alpha = alpha,
                        showsegment = false)
        Makie.plot!(ax4, Proj(Mercator)(admin1_geometry), 
                        color = lift(update_color, month_index, SUBNAT_ACCESS_collection[ISO_i][i,:]), 
                        colormap = access_cmap,
                        colorrange = ACCESS_clims,
                        alpha = alpha,
                        showsegment = false)

        Makie.plot!(ax7, Proj(Mercator)(admin1_geometry), 
                        color = lift(update_color, month_index, NPC_GAP_collection[ISO_i][i,:]), 
                        colormap = gap_cmap,
                        colorrange = NPC_GAP_clims,
                        showsegment = false)
        Makie.plot!(ax8, Proj(Mercator)(admin1_geometry), 
                        color = lift(update_color, month_index, ACCESS_GAP_collection[ISO_i][i,:]), 
                        colormap = gap_cmap,
                        colorrange = ACCESS_GAP_clims,
                        showsegment = false)

        Makie.plot!(ax9, Proj(Mercator)(admin1_geometry), 
                        color = lift(update_color, month_index, NPC_GAP_normalised_collection[ISO_i][i,:]), 
                        colormap = norm_gap_cmap,
                        colorrange = NPC_GAP_normalised_clims,
                        showsegment = false)
        Makie.plot!(ax10, Proj(Mercator)(admin1_geometry), 
                        color = lift(update_color, month_index, ACCESS_GAP_normalised_collection[ISO_i][i,:]), 
                        colormap = norm_gap_cmap,
                        colorrange = ACCESS_GAP_normalised_clims,
                        showsegment = false)
        
        Makie.plot!(ax11, Proj(Mercator)(admin1_geometry), 
                        color = lift(update_color, month_index, BV_MITN_SUBNAT_NPC_GAP_collection[ISO_i][i,:]), 
                        colormap = gap_cmap,
                        colorrange = NPC_GAP_clims,
                        showsegment = false)
        Makie.plot!(ax12, Proj(Mercator)(admin1_geometry), 
                        color = lift(update_color, month_index, BV_MITN_SUBNAT_NPC_GAP_normalised_collection[ISO_i][i,:]), 
                        colormap = norm_gap_cmap,
                        colorrange = NPC_GAP_normalised_clims,
                        showsegment = false)
    end
end
hidedecorations!(ax1)
hidedecorations!(ax2)
hidedecorations!(ax3)
hidedecorations!(ax4)
hidedecorations!(ax5)
hidedecorations!(ax6)
hidedecorations!(ax7)
hidedecorations!(ax8)
hidedecorations!(ax9)
hidedecorations!(ax10)
hidedecorations!(ax11)
hidedecorations!(ax12)
fig
