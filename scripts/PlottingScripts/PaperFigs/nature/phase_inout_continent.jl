"""
Author: Eugene Tan
Date Created: 12/5/2025
Last Updated: 12/5/2025
Make Net Crop plots to show estimated breakdown of nets by type
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

# %% Maths packages
using LinearAlgebra
using StatsBase

# %% General useful functions
using DateConversions

# %% Plot packages
using LaTeXStrings
using CairoMakie

# %% Plotting Theme and general settings
set_theme!(theme_ggplot2())

# Color settings
colors = [  colorant"#094D92", # cITN
            colorant"#EE4266", # LLIN
            colorant"#38B673", # PBO
            colorant"#EEC643", # G2
            colorant"#8775C9" # ROYAL
            ];

# General settings
fillalpha = 0.65
la = 0.5
lw = 2
label_fontsize = 25

# %%
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END
n_months = 12*(YEAR_END-YEAR_START+1)

# %% Import and pre-process raster output data. Used to adjust SNF estimates of type and rake against raster estimates.
raster_timeseries = CSV.read(OUTPUT_DIR*"coverage_timeseries/master_extraction.csv", DataFrame)

# %% Define list of countries to extract data from and plot
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %%
NET_NAMES = ["cITN", "LLIN", "PBO", "G2", "ROYAL"]

###################################################
# %% POSTPROCESSING: Calculate Raster Estimates of Net Crop disaggrated by type (According to SNF ratios)
###################################################
Γ_BYNET_RASTER_combined = zeros(length(filt_ISOs), 12*(YEAR_END-YEAR_START+1), length(NET_NAMES))

# %%
for ISO_i in 1:length(filt_ISOs)
    ISO = filt_ISOs[ISO_i]

    Γ_BYNET_SNF = mean(JLD2.load(OUTPUT_DIR*"predictions/$(ISO)_netcrop_prediction.jld2")["Γ_BYNET_pred_samples"], dims = 1)[1,:,:]
    Γ_TOTAL_SNF = sum(Γ_BYNET_SNF, dims = 2)[:]

    # Extract net type ratios from predictive SNF
    Γ_BYNET_RATIOS = Γ_BYNET_SNF./repeat(Γ_TOTAL_SNF, 1, length(NET_NAMES))
    Γ_BYNET_RATIOS[isnan.(Γ_BYNET_RATIOS)] .= 0

    # Adjust for periods early in the time series where predictive gives 0 nets, because forward posterior predictive
    # doesn't impute missing distribution data
    for idx in findall(sum(Γ_BYNET_RATIOS, dims = 2)[:,1] .== 0)
        year = YEAR_START + monthidx_to_monthyear(idx)[2] - 1
        if year < YEAR_SUBNAT_TRANS
            Γ_BYNET_RATIOS[idx,findfirst(NET_NAMES .== "cITN")] = 1
        end
    end

    # Extract total net crop time series from raster extraction
    raster_country_timeseries = raster_timeseries[raster_timeseries.category .== "Admin0" .&&
                        raster_timeseries.ISO .== ISO,:]
    Γ_TOTAL_RASTER = zeros(size(Γ_BYNET_SNF)[1])

    for year in YEAR_START:YEAR_END
        for month in 1:12
            monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEAR_START)
            population, npc = raster_country_timeseries[raster_country_timeseries.month .== month .&&
                                raster_country_timeseries.year .== year,["population", "snf_npc_mean"]][1,:]
            Γ_TOTAL_RASTER[monthidx] = population*npc
        end
    end

    # %% Disaggregate raster crop estimate according to ratios from SNF
    Γ_BYNET_RASTER_combined[ISO_i,:,:] .= repeat(Γ_TOTAL_RASTER, 1, size(Γ_BYNET_RATIOS)[2]).*Γ_BYNET_RATIOS
end

###################################################
# %% Figure 1: Post Net crop breakdown plot by type for each country
###################################################

for ISO_i in 1:length(filt_ISOs)
    ISO = filt_ISOs[ISO_i]
    fig = Figure(size = (800,500))
    ax = Axis(fig[1,1], title = "$(ISO) Net Breakdown",
                xlabel = "Years", 
                xticks = (1:12:n_months, string.(YEAR_START:YEAR_END)),
                xticklabelrotation = pi/2,
                xlabelsize = 20,
                titlesize = 25,
                ylabel = "Net Crop (mil)",
                ylabelsize = 20
                )
    legend_bands = []
    legend_net_bool = []

    for i in 1:length(NET_NAMES)
        if sum(Γ_BYNET_RASTER_combined[ISO_i,:,i]) > 0
            if i == 1
                band = band!(ax, 1:n_months, zeros(n_months), Γ_BYNET_RASTER_combined[ISO_i,:,i]./1e6,
                        color = (colors[i], fillalpha))
            else
                band = band!(ax, 1:n_months, sum(Γ_BYNET_RASTER_combined[ISO_i,:,1:i-1], dims = 2)[:]./1e6, 
                        sum(Γ_BYNET_RASTER_combined[ISO_i,:,1:i], dims = 2)[:]./1e6,
                        color = (colors[i], fillalpha))
            end
            push!(legend_bands, band)
            push!(legend_net_bool, i)
        end
    end
    Legend(fig[1, 1],
            legend_bands,
            NET_NAMES[legend_net_bool],
            tellheight = false,
            tellwidth = false,
            margin = (15, 10, 10, 10),
            halign = :left, valign = :top)
    # Save fig
    mkpath(OUTPUT_PLOTS_DIR*"PaperFigures/Country_ITN_Types/")
    save(OUTPUT_PLOTS_DIR*"PaperFigures/Country_ITN_Types/$(ISO)_nettype.pdf", fig, pdf_version = "1.4")
end

###################################################
# %% Figure 2: Post Net crop breakdown plot for entire continent
###################################################

# Calculate continent level totals and ratios for Net Crop
Γ_BYNET_RASTER_continent = sum(Γ_BYNET_RASTER_combined, dims = 1)[1,:,:]
RATIOS_BYNET_RASTER_continent = Γ_BYNET_RASTER_continent./repeat(sum(Γ_BYNET_RASTER_continent, dims = 2)[:,1], 1, length(NET_NAMES))

# Make plot
fig = Figure(size = (800,600))
ax_crop = Axis(fig[1,1], title = "Africa Net Breakdown",
            xlabel = "Years", 
            xticks = (1:12:n_months, string.(YEAR_START:YEAR_END)),
            xticklabelrotation = pi/2,
            xlabelsize = 20,
            titlesize = 25,
            ylabel = "Net Crop (mil)",
            ylabelsize = 20
            )
ax_ratio = Axis(fig[2,1],
            xlabel = "Years", 
            xticks = (1:12:n_months, string.(YEAR_START:YEAR_END)),
            xticklabelrotation = pi/2,
            xlabelsize = 20,
            titlesize = 25,
            yticks = 0:20:100,
            ylabel = "Crop Share (%)",
            ylabelsize = 20
            )

xlims!(ax_crop, -1,n_months+1)
xlims!(ax_ratio, -1,n_months+1)
ylims!(ax_ratio, -1,101)
hidexdecorations!(ax_crop, ticks = false, grid = false)

legend_bands = []
legend_net_bool = []

for i in 1:length(NET_NAMES)
    if sum(Γ_BYNET_RASTER_continent[:,i]) > 0
        if i == 1
            band = band!(ax_crop, 1:n_months, zeros(n_months), Γ_BYNET_RASTER_continent[:,i]./1e6,
                    color = (colors[i], fillalpha))
            band!(ax_ratio, 1:n_months, zeros(n_months), RATIOS_BYNET_RASTER_continent[:,i].*100,
                    color = (colors[i], fillalpha))
        else
            band = band!(ax_crop, 1:n_months, sum(Γ_BYNET_RASTER_continent[:,1:i-1], dims = 2)[:]./1e6, 
                    sum(Γ_BYNET_RASTER_continent[:,1:i], dims = 2)[:]./1e6,
                    color = (colors[i], fillalpha))
            band!(ax_ratio, 1:n_months, sum(RATIOS_BYNET_RASTER_continent[:,1:i-1], dims = 2)[:].*100, 
                    sum(RATIOS_BYNET_RASTER_continent[:,1:i], dims = 2)[:].*100,
                    color = (colors[i], fillalpha))
        end
        push!(legend_bands, band)
        push!(legend_net_bool, i)
    end
end
Legend(fig[1, 1],
        legend_bands,
        NET_NAMES[legend_net_bool],
        tellheight = false,
        tellwidth = false,
        margin = (15, 10, 10, 10),
        halign = :left, valign = :top)
fig

# %%
# Save fig
save(OUTPUT_PLOTS_DIR*"PaperFigures/Continent_ITN_nettype.pdf", fig, pdf_version = "1.4")

