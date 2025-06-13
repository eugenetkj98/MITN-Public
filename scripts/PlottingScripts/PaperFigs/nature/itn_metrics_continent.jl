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
using Missings
using JLD2
using CSV
using ProgressBars

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

# %% Dataset Directories
input_dir = OUTPUT_DIR*"coverage_timeseries/"
input_filename = "master_extraction.csv"

# %% Load Required Datasets
raster_timeseries = CSV.read(input_dir*input_filename, DataFrame)

# %% Population values have annual resolution. So just temporarily calculated the monthly interpolated values
Threads.@threads for i in ProgressBar(1:size(raster_timeseries)[1])
    area_id, month_val, year_val, population = raster_timeseries[i,["area_id","month","year","population"]]
    interp_pop = population
    if month_val != 1 && year_val < 2023
        # Get pop values to interpolate between
        pop_lower = raster_timeseries[raster_timeseries.area_id .== area_id .&&
                                        raster_timeseries.month .== 1 .&&
                                        raster_timeseries.year .== year_val,"population"][1]
        pop_upper = raster_timeseries[raster_timeseries.area_id .== area_id .&&
                                        raster_timeseries.month .== 1 .&&
                                        raster_timeseries.year .== year_val+1,"population"][1]
        interp_pop = (1-(month_val-1)/12)*pop_lower + ((month_val-1)/12)*pop_upper

        # replace raw population with interpolated value
        raster_timeseries[i,"population"] = interp_pop
    end
end

# %% Define year bounds
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END
n_years = YEAR_END-YEAR_START+1
n_months = n_years*12

# %% Construct continent time series to plot
filt_data = raster_timeseries[raster_timeseries.category .== "Continent",:]

continent_population = zeros(n_months)
continent_npc = zeros(n_months, 3)
continent_access = zeros(n_months, 3)
continent_use = zeros(n_months, 3)
continent_util = zeros(n_months, 3)
continent_eff = zeros(n_months, 3)

for monthidx in 1:n_months
    # Get timestamps
    month_val, year_idx = monthidx_to_monthyear(monthidx)
    year_val = (YEAR_START:YEAR_END)[year_idx]

    # Further filter data and extract data
    data_slice = filt_data[filt_data.month .== month_val .&& filt_data.year .== year_val,:]
    continent_population[monthidx] = sum(data_slice[1,:].population)

    
    continent_npc[monthidx,1] = data_slice[1,"raster_npc_95lower"]
    continent_npc[monthidx,2] = data_slice[1,"raster_npc_mean"]
    continent_npc[monthidx,3] = data_slice[1,"raster_npc_95upper"]

    continent_access[monthidx,1] = data_slice[1,"raster_access_95lower"]
    continent_access[monthidx,2] = data_slice[1,"raster_access_mean"]
    continent_access[monthidx,3] = data_slice[1,"raster_access_95upper"]

    continent_use[monthidx,1] = data_slice[1,"raster_use_95lower"]
    continent_use[monthidx,2] = data_slice[1,"raster_use_mean"]
    continent_use[monthidx,3] = data_slice[1,"raster_use_95upper"]

    continent_util[monthidx,1] = data_slice[1,"raster_util_95lower"]
    continent_util[monthidx,2] = data_slice[1,"raster_util_mean"]
    continent_util[monthidx,3] = data_slice[1,"raster_util_95upper"]

    continent_eff[monthidx,1] = data_slice[1,"raster_eff_95lower"]
    continent_eff[monthidx,2] = data_slice[1,"raster_eff_mean"]
    continent_eff[monthidx,3] = data_slice[1,"raster_eff_95upper"]
end

# %% Calculate annual downsampled


continent_npc_annual = zeros(n_years,3)
continent_access_annual = zeros(n_years,3)
continent_use_annual = zeros(n_years,3)
continent_util_annual = zeros(n_years,3)
continent_eff_annual = zeros(n_years,3)

continent_npc_annual[:,1] = [mean(continent_npc[((i-1)*12+1):((i-1)*12+12),1]) for i in 1:n_years]
continent_npc_annual[:,2] = [mean(continent_npc[((i-1)*12+1):((i-1)*12+12),2]) for i in 1:n_years]
continent_npc_annual[:,3] = [mean(continent_npc[((i-1)*12+1):((i-1)*12+12),3]) for i in 1:n_years]

continent_access_annual[:,1] = [mean(continent_access[((i-1)*12+1):((i-1)*12+12),1]) for i in 1:n_years]
continent_access_annual[:,2] = [mean(continent_access[((i-1)*12+1):((i-1)*12+12),2]) for i in 1:n_years]
continent_access_annual[:,3] = [mean(continent_access[((i-1)*12+1):((i-1)*12+12),3]) for i in 1:n_years]

continent_use_annual[:,1] = [mean(continent_use[((i-1)*12+1):((i-1)*12+12),1]) for i in 1:n_years]
continent_use_annual[:,2] = [mean(continent_use[((i-1)*12+1):((i-1)*12+12),2]) for i in 1:n_years]
continent_use_annual[:,3] = [mean(continent_use[((i-1)*12+1):((i-1)*12+12),3]) for i in 1:n_years]

continent_util_annual[:,1] = [mean(continent_util[((i-1)*12+1):((i-1)*12+12),1]) for i in 1:n_years]
continent_util_annual[:,2] = [mean(continent_util[((i-1)*12+1):((i-1)*12+12),2]) for i in 1:n_years]
continent_util_annual[:,3] = [mean(continent_util[((i-1)*12+1):((i-1)*12+12),3]) for i in 1:n_years]

continent_eff_annual[:,1] = [mean(continent_eff[((i-1)*12+1):((i-1)*12+12),1]) for i in 1:n_years]
continent_eff_annual[:,2] = [mean(continent_eff[((i-1)*12+1):((i-1)*12+12),2]) for i in 1:n_years]
continent_eff_annual[:,3] = [mean(continent_eff[((i-1)*12+1):((i-1)*12+12),3]) for i in 1:n_years]


####################################################################
# %% FIGURE 1: ITN Coverage plots for Africa continent
####################################################################
# %% Construct plot
# Base axes
fig = Figure(size = (800,500))
ax = Axis(fig[1,1],
        title = "Continent ITN Metrics",
        xlabel = "Years", 
        xticks = (1:12:n_months, string.(YEAR_START:YEAR_END)),
        xticklabelrotation = pi/2,
        xlabelsize = 20,
        titlesize = 25,
        ylabel = "Coverage Metric",
        yticks = (0:0.2:1),
        ylabelsize = 20
        )
xlims!(-0.5,n_months+0.5)
ylims!(-0.02, 1.02)

# Add lines

### NPC
band!(ax, 1:n_months, continent_npc[:,1], continent_npc[:,3],
        color = (colors[1], fillalpha)
        )
npc_line = lines!(ax, 1:n_months, continent_npc[:,2],
        color = colors[1], linewidth = lw)
lines!(ax, 1:n_months, continent_npc[:,1],
        color = (colors[1], la), 
        linewidth = lw, linestyle = :dash)
lines!(ax, 1:n_months, continent_npc[:,3],
        color = (colors[1], la), linewidth = lw, 
        linestyle = :dash)

### Access
band!(ax, 1:n_months, continent_access[:,1], continent_access[:,3],
color = (colors[2], fillalpha)
)
access_line = lines!(ax, 1:n_months, continent_access[:,2],
        color = colors[2], linewidth = lw)
lines!(ax, 1:n_months, continent_access[:,1],
        color = (colors[2], la), 
        linewidth = lw, linestyle = :dash)
lines!(ax, 1:n_months, continent_access[:,3],
        color = (colors[2], la), linewidth = lw, 
        linestyle = :dash)

### Use
band!(ax, 1:n_months, continent_use[:,1], continent_use[:,3],
color = (colors[3], fillalpha)
)
use_line = lines!(ax, 1:n_months, continent_use[:,2],
        color = colors[3], linewidth = lw)
lines!(ax, 1:n_months, continent_use[:,1],
        color = (colors[3], la), 
        linewidth = lw, linestyle = :dash)
lines!(ax, 1:n_months, continent_use[:,3],
        color = (colors[3], la), linewidth = lw, 
        linestyle = :dash)

# Legend
Legend(fig[1, 2],
    [npc_line, access_line, use_line],
    ["NPC", "Access", "Use"])

save(OUTPUT_PLOTS_DIR*"PaperFigures/Continent_ITN_Metrics_monthly.pdf", fig, pdf_version = "1.4")
fig

####################################################################
# %% FIGURE 1B: ITN Coverage plots for Africa continent (ANNUAL)
####################################################################
# %% Construct plot
# Base axes
fig = Figure(size = (800,500))
ax = Axis(fig[1,1],
        title = "Continent ITN Metrics",
        xlabel = "Years", 
        xticks = (1:n_years, string.(YEAR_START:YEAR_END)),
        xticklabelrotation = pi/2,
        xlabelsize = 20,
        titlesize = 25,
        ylabel = "Coverage Metric",
        yticks = (0:0.2:1),
        ylabelsize = 20
        )
xlims!(0.5,n_years+0.5)
ylims!(-0.02, 1.02)

# Add lines

### NPC
band!(ax, 1:n_years, continent_npc_annual[:,1], continent_npc_annual[:,3],
        color = (colors[1], fillalpha)
        )
npc_line = lines!(ax, 1:n_years, continent_npc_annual[:,2],
        color = colors[1], linewidth = lw)
lines!(ax, 1:n_years, continent_npc_annual[:,1],
        color = (colors[1], la), 
        linewidth = lw, linestyle = :dash)
lines!(ax, 1:n_years, continent_npc_annual[:,3],
        color = (colors[1], la), linewidth = lw, 
        linestyle = :dash)

### Access
band!(ax, 1:n_years, continent_access_annual[:,1], continent_access_annual[:,3],
color = (colors[2], fillalpha)
)
access_line = lines!(ax, 1:n_years, continent_access_annual[:,2],
        color = colors[2], linewidth = lw)
lines!(ax, 1:n_years, continent_access_annual[:,1],
        color = (colors[2], la), 
        linewidth = lw, linestyle = :dash)
lines!(ax, 1:n_years, continent_access_annual[:,3],
        color = (colors[2], la), linewidth = lw, 
        linestyle = :dash)

### Use
band!(ax, 1:n_years, continent_use_annual[:,1], continent_use_annual[:,3],
color = (colors[3], fillalpha)
)
use_line = lines!(ax, 1:n_years, continent_use_annual[:,2],
        color = colors[3], linewidth = lw)
lines!(ax, 1:n_years, continent_use_annual[:,1],
        color = (colors[3], la), 
        linewidth = lw, linestyle = :dash)
lines!(ax, 1:n_years, continent_use_annual[:,3],
        color = (colors[3], la), linewidth = lw, 
        linestyle = :dash)

# Legend
Legend(fig[1, 2],
    [npc_line, access_line, use_line],
    ["NPC", "Access", "Use"])

save(OUTPUT_PLOTS_DIR*"PaperFigures/Continent_ITN_Metrics_annual.pdf", fig, pdf_version = "1.4")
fig

####################################################################
# %% FIGURE 4: Continent Utilisation Plot (MONTHLY)
####################################################################
# Make figure
layout_res = (1500,1300)
fig = Figure(size = (800,500))

ax = Axis(fig[1,1],
        title = "Continent ITN Metrics",
        xlabel = "Years", 
        xticks = (1:12:n_months, string.(YEAR_START:YEAR_END)),
        xticklabelrotation = pi/2,
        xlabelsize = 20,
        titlesize = 25,
        ylabel = "Metric Rate",
        yticks = (0:0.2:2),
        ylabelsize = 20
        )
xlims!(-0.5,n_months+0.5)
ylims!(-0.02, 2.02)

# Add Line
band!(ax, 1:n_months, continent_util[:,1], continent_util[:,3],
                color = (colors[1], fillalpha)
                )
util_line = lines!(ax, 1:n_months, continent_util[:,2],
        color = colors[1], linewidth = lw)
hlines!(ax, [1], color = :black, linestyle = :dash, linewidth = lw)

# Add Line
band!(ax, 1:n_months, continent_eff[:,1], continent_eff[:,3],
                color = (colors[2], fillalpha)
                )
eff_line = lines!(ax, 1:n_months, continent_eff[:,2],
        color = colors[2], linewidth = lw)
hlines!(ax, [2], color = :black, linestyle = :dash, linewidth = lw)

# Legend
Legend(fig[1, 2],
    [util_line, eff_line],
    ["Utilisation (η)", "Efficiency (α)"])

save(OUTPUT_PLOTS_DIR*"PaperFigures/Continent_ITN_utilisation_monthly.pdf", fig, pdf_version = "1.4")
fig

####################################################################
# %% FIGURE 4B: Continent Utilisation Plot (ANNUAL)
####################################################################
# Make figure
layout_res = (1500,1300)
fig = Figure(size = (800,500))

ax = Axis(fig[1,1],
        title = "Continent ITN Metrics",
        xlabel = "Years", 
        xticks = (1:n_years, string.(YEAR_START:YEAR_END)),
        xticklabelrotation = pi/2,
        xlabelsize = 20,
        titlesize = 25,
        ylabel = "Metric Rate",
        yticks = (0:0.2:2),
        ylabelsize = 20
        )
xlims!(0.5,n_years+0.5)
ylims!(-0.02, 2.02)

# Add Line
band!(ax, 1:n_years, continent_util_annual[:,1], continent_util_annual[:,3],
                color = (colors[1], fillalpha)
                )
util_line = lines!(ax, 1:n_years, continent_util_annual[:,2],
        color = colors[1], linewidth = lw)
hlines!(ax, [1], color = :black, linestyle = :dash, linewidth = lw)

# Add Line
band!(ax, 1:n_years, continent_eff_annual[:,1], continent_eff_annual[:,3],
                color = (colors[2], fillalpha)
                )
eff_line = lines!(ax, 1:n_years, continent_eff_annual[:,2],
        color = colors[2], linewidth = lw)

# Legend
Legend(fig[1, 2],
    [util_line, eff_line],
    ["Utilisation (η)", "Efficiency (α)"])

save(OUTPUT_PLOTS_DIR*"PaperFigures/Continent_ITN_utilisation_annual.pdf", fig, pdf_version = "1.4")
fig

####################################################################
# %% FIGURE 2: Moving Window Gains plot
####################################################################
# Linearly extrapolate population for period of 2021-2023
est_pop_change = continent_population[20*12]-continent_population[20*12-1]

est_population = continent_population
for monthidx in 20*12+1:length(continent_population)
        est_population[monthidx] = est_population[monthidx-1] + est_pop_change
end

# Mean downsample to annual data
continent_crop_annual = [mean((continent_npc.*continent_population)[((i-1)*12+1):(12*i),2]) for i in 1:n_years]
continent_pop_annual = [mean(est_population[((i-1)*12+1):(12*i)]) for i in 1:n_years]

# Calculate Moving Averages
ma_crop = MA_filter(continent_crop_annual, window = 3)
ma_pop = MA_filter(continent_pop_annual, window = 3)

# Calculate percentage gains
delta_crop = ma_crop[2:end]./ma_crop[1:end-1]
delta_pop = ma_pop[2:end]./ma_pop[1:end-1]

# 
fig = Figure(size = (1000,400))

ax = Axis(fig[1,1], title = "Africa Continent Population & ITN Trends",
            xlabel = "Year",
            xticks = (1:2:n_years, (string.(YEAR_START:YEAR_END))[1:2:end]),
            xticklabelrotation = pi/2,
            ylabel = "3 Year Δ%",
            yticks = -20:10:100,
            xlabelsize = 20, ylabelsize = 20,
            titlesize = 23)
xlims!(ax, 0, n_years + 1)
ylims!(ax, -25, 110)

# Plot barplot
data_tbl = (year_cat = repeat(1:length(delta_crop),2) .+ 1,
                vals = vcat(delta_crop, delta_pop).-1,
                grp_dodge = vcat(repeat([1], length(delta_crop)), repeat([2], length(delta_pop))),
                grp_col = colors[vcat(repeat([1], length(delta_crop)), repeat([2], length(delta_pop)))],
                pop = ma_pop)
barplot!(ax, data_tbl.year_cat, data_tbl.vals.*100,
        dodge = data_tbl.grp_dodge,
        dodge_gap = 0.1,
        color = data_tbl.grp_col,
        width = 0.6)
        colors

# Add Legend
labels = ["Net Crop", "Population"]
elements = [PolyElement(polycolor = colors[i]) for i in 1:length(labels)]
Legend(
        fig[1, 1], elements, labels,
        tellheight = false,
        tellwidth = false,
        margin = (10, 10, 10, 10),
        halign = :right, valign = :top
    ) 

# Add numerical labels on surplus for each year
delta_diff = delta_crop - delta_pop
text_color = []
for i in 1:length(delta_diff)
    if delta_diff[i] > 0
        push!(text_color, colors[1])
    else
        push!(text_color, colors[2])
    end
end

text!([(i+1,-10) for i in 1:length(delta_crop)],
        text = string.(round.(Int,delta_diff*100)).*"%",
        align = (:center, :center),
        color = text_color, fontsize = 13)

# Save fig
save(OUTPUT_PLOTS_DIR*"PaperFigures/Continent_ITN_trend.pdf", fig, pdf_version = "1.4")

fig