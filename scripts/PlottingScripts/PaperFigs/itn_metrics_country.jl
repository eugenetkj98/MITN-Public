"""
Author: Eugene Tan
Date Created: 8/5/2025
Last Updated: 8/5/2025
Make country level plots for paper
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
africa_lw = 1.4

# %% Dataset Directories
input_dir = OUTPUT_DIR*"coverage_timeseries/"
input_filename = "master_extraction.csv"

# Net age timeseries directory
net_age_dir = OUTPUT_DATAPREP_DIR
net_age_filename = "snf_mean_netage.csv"

# %% Load Required Datasets
raster_timeseries = CSV.read(input_dir*input_filename, DataFrame)
netage_timeseries = CSV.read(net_age_dir*net_age_filename, DataFrame)

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

# %% Define list of countries to plot
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Storage variables
country_population = zeros(length(filt_ISOs), n_months)
country_crop = zeros(length(filt_ISOs), n_months, 3)
country_access_pop = zeros(length(filt_ISOs), n_months, 3)
country_use_pop = zeros(length(filt_ISOs), n_months, 3)
country_npc = zeros(length(filt_ISOs), n_months, 3)
country_access = zeros(length(filt_ISOs), n_months, 3)
country_use = zeros(length(filt_ISOs), n_months, 3)
country_utilisation = zeros(length(filt_ISOs), n_months, 3)


# %% Select ISO
for ISO_i in 1:length(filt_ISOs)
    ISO = filt_ISOs[ISO_i]
    # %% Construct continent time series to plot
    filt_data = raster_timeseries[raster_timeseries.category .== "Admin0" .&&
                                    raster_timeseries.ISO .== ISO,:]

    for monthidx in 1:n_months
        # Get timestamps
        month_val, year_idx = monthidx_to_monthyear(monthidx)
        year_val = (YEAR_START:YEAR_END)[year_idx]

        # Further filter data and extract data
        data_slice = filt_data[filt_data.month .== month_val .&& filt_data.year .== year_val,:]

        country_population[ISO_i,monthidx] = sum(data_slice.population)

        country_crop[ISO_i,monthidx,1] = sum(data_slice.raster_npc_95lower.*data_slice.population)
        country_crop[ISO_i,monthidx,2] = sum(data_slice.raster_npc_mean.*data_slice.population)
        country_crop[ISO_i,monthidx,3] = sum(data_slice.raster_npc_95upper.*data_slice.population)

        country_access_pop[ISO_i,monthidx,1] = sum(data_slice.raster_access_95lower.*data_slice.population)
        country_access_pop[ISO_i,monthidx,2] = sum(data_slice.raster_access_mean.*data_slice.population)
        country_access_pop[ISO_i,monthidx,3] = sum(data_slice.raster_access_95upper.*data_slice.population)

        country_use_pop[ISO_i,monthidx,1] = sum(data_slice.raster_use_95lower.*data_slice.population)
        country_use_pop[ISO_i,monthidx,2] = sum(data_slice.raster_use_mean.*data_slice.population)
        country_use_pop[ISO_i,monthidx,3] = sum(data_slice.raster_use_95upper.*data_slice.population)
    end

    country_npc[ISO_i,:,:] = country_crop[ISO_i,:,:]./repeat(country_population[ISO_i,:],1,3)
    country_access[ISO_i,:,:] = country_access_pop[ISO_i,:,:]./repeat(country_population[ISO_i,:],1,3)
    country_use[ISO_i,:,:] = country_use_pop[ISO_i,:,:]./repeat(country_population[ISO_i,:],1,3)
    
    # Calculate utilisation rates. CI are worst case
    country_utilisation[ISO_i,:,1] = country_use_pop[ISO_i,:,1]./country_access_pop[ISO_i,:,3]
    country_utilisation[ISO_i,:,2] = country_use_pop[ISO_i,:,2]./country_access_pop[ISO_i,:,2]
    country_utilisation[ISO_i,:,3] = country_use_pop[ISO_i,:,3]./country_access_pop[ISO_i,:,1]
end

####################################################################
# %% FIGURE 1: ITN Coverage plots for each individual country
####################################################################

for ISO_i in ProgressBar(1:length(filt_ISOs))
    # Select ISO
    ISO = filt_ISOs[ISO_i]
    # Base axes
    fig = Figure(size = (800,500))
    ax = Axis(fig[1,1],
            title = "$(ISO) ITN Metrics",
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
    band!(ax, 1:n_months, country_npc[ISO_i,:,1], country_npc[ISO_i,:,3],
            color = (colors[1], fillalpha)
            )
    npc_line = lines!(ax, 1:n_months, country_npc[ISO_i,:,2],
            color = colors[1], linewidth = lw)
    lines!(ax, 1:n_months, country_npc[ISO_i,:,1],
            color = (colors[1], la), 
            linewidth = lw, linestyle = :dash)
    lines!(ax, 1:n_months, country_npc[ISO_i,:,3],
            color = (colors[1], la), linewidth = lw, 
            linestyle = :dash)

    ### Access
    band!(ax, 1:n_months, country_access[ISO_i,:,1], country_access[ISO_i,:,3],
    color = (colors[2], fillalpha)
    )
    access_line = lines!(ax, 1:n_months, country_access[ISO_i,:,2],
            color = colors[2], linewidth = lw)
    lines!(ax, 1:n_months, country_access[ISO_i,:,1],
            color = (colors[2], la), 
            linewidth = lw, linestyle = :dash)
    lines!(ax, 1:n_months, country_access[ISO_i,:,3],
            color = (colors[2], la), linewidth = lw, 
            linestyle = :dash)

    ### Use
    band!(ax, 1:n_months, country_use[ISO_i,:,1], country_use[ISO_i,:,3],
    color = (colors[2], fillalpha)
    )
    use_line = lines!(ax, 1:n_months, country_use[ISO_i,:,2],
            color = colors[3], linewidth = lw)
    lines!(ax, 1:n_months, country_use[ISO_i,:,1],
            color = (colors[3], la), 
            linewidth = lw, linestyle = :dash)
    lines!(ax, 1:n_months, country_use[ISO_i,:,3],
            color = (colors[3], la), linewidth = lw, 
            linestyle = :dash)

    # Legend
    Legend(fig[1, 2],
        [npc_line, access_line, use_line],
        ["NPC", "Access", "Use"])

    save(OUTPUT_PLOTS_DIR*"PaperFigures/Country_ITN_Metrics/$(ISO)_ITN_Metrics.pdf", fig, pdf_version = "1.4")
end
fig

####################################################################
# %% FIGURE 2: Make summary figure tiled by geographical location in africa
####################################################################
# Define index locations in Africa
subplot_lookup = Dict(   "MRT" => (1,2,"Mauritania"),
                        "ERI" => (1,8,"Eritrea"),
                        "GMB" => (2,1,"Gambia"),
                        "SEN" => (2,2,"Senegal"),
                        "GNB" => (2,3,"Guinea-Bissau"),
                        "MLI" => (2,4,"Mali"),
                        "NER" => (2,5,"Niger"),
                        "TCD" => (2,6,"Chad"),
                        "SDN" => (2,7,"Sudan"),
                        "DJI" => (2,8,"Djibouti"),
                        "SLE" => (3,1,"Sierra Leone"),
                        "GIN" => (3,2,"Guinea"),
                        "GHA" => (3,3,"Ghana"),
                        "BFA" => (3,4,"Burkina Faso"),
                        "BEN" => (3,5,"Benin"),
                        "CAF" => (3,6,"Centr. Afr. Rep."),
                        "SSD" => (3,7,"South Sudan"),
                        "ETH" => (3,8,"Ethiopia"),
                        "SOM" => (3,9,"Somalia"),
                        "LBR" => (4,2,"Liberia"),
                        "CIV" => (4,3,"Cote d'Ivoire"),
                        "TGO" => (4,4,"Togo"),
                        "NGA" => (4,5,"Nigeria"),
                        "CMR" => (4,6,"Cameroon"),
                        "RWA" => (4,7,"Rwanda"),
                        "KEN" => (4,8,"Kenya"),
                        "STP" => (5,3,"Sao Tome & Prin."),
                        "GAB" => (5,4,"Gabon"),
                        "COG" => (5,5,"Rep. of Congo"),
                        "GNQ" => (5,6,"Equatorial Guinea"),
                        "BDI" => (5,7,"Burundi"),
                        "UGA" => (5,8,"Uganda"),
                        "AGO" => (6,4,"Angola"),
                        "COD" => (6,5,"Dem. Rep. Congo"),
                        "ZMB" => (6,6,"Zambia"),
                        "TZA" => (6,7,"Tanzania"),
                        "NAM" => (7,3,"Namibia"),
                        "BWA" => (7,4,"Botswana"),
                        "ZWE" => (7,5,"Zimbabwe"),
                        "MWI" => (7,6,"Malawi"),
                        "MOZ" => (7,7,"Mozambique"),
                        "SWZ" => (7,8,"Eswatini"),
                        "COM" => (7,9,"Comoros"),
                        "MDG" => (7,10,"Madagascar"))
show_xlabel_ISOs = ["SLE","LBR","STP","NAM","BWA",
                    "ZWE","MWI","MOZ","SWZ","COM","MDG",
                    "UGA","SOM"]
show_ylabel_ISOs = ["MRT","ERI","GMB","SLE","LBR","STP",
                    "AGO","NAM"]

# construct year strings for xticks
year_strings = Vector{String}(undef, length(YEAR_START:YEAR_END))
for year_i in 1:length(YEAR_START:YEAR_END)
    year = (YEAR_START:YEAR_END)[year_i]
    year_mod_val = year%2000
    if year_mod_val < 10
        year_strings[year_i] = "'0"*string(year_mod_val)
    else
        year_strings[year_i] = "'"*string(year_mod_val)
    end
end

# Make figure
layout_res = (1500,1300)
fig = Figure(size = layout_res)

plot_axs = []
title_axs = []
for ISO_i in 1:length(filt_ISOs)
    ISO = filt_ISOs[ISO_i]
    x_idx, y_idx, country_name = subplot_lookup[ISO]
    ax = Axis(fig[2*x_idx,y_idx], width = layout_res[1]/14,
                # xlabel = "Years", 
                xticks = ((1:12:n_months)[1:5:end], year_strings[1:5:end]),
                xticklabelrotation = pi/2,
                # ylabel = "Metric",
                yticks = (0:0.2:1),
                ylabelsize = 20)
    xlims!(ax,-0.5,n_months+0.5)
    ylims!(ax,-0.02, 1.02)

    lb = Label(fig[2*(x_idx-1)+1,y_idx], "$(country_name)")
    
    # Hide labels as required
    if !(ISO ∈ show_xlabel_ISOs)
        ax.xticklabelsvisible = false
    end

    if !(ISO ∈ show_ylabel_ISOs)
        ax.yticklabelsvisible = false
    end
    push!(plot_axs, ax)
    push!(title_axs, lb)
end

# Super Labels
Label(fig[:,0],"Coverage Metric", rotation = pi/2, fontsize = 40)
Label(fig[15,:],"Years", fontsize = 40)

# Add plot lines for metric
for ISO_i in ProgressBar(1:length(filt_ISOs), leave = false)
    # Select ISO
    ISO = filt_ISOs[ISO_i]

    # Add lines

    ### NPC
    npc_line = lines!(plot_axs[ISO_i], 1:n_months, country_npc[ISO_i,:,2],
            color = colors[1], linewidth = africa_lw)

    ### Access
    access_line = lines!(plot_axs[ISO_i], 1:n_months, country_access[ISO_i,:,2],
            color = colors[2], linewidth = africa_lw)

    ### Use
    use_line = lines!(plot_axs[ISO_i], 1:n_months, country_use[ISO_i,:,2],
            color = colors[3], linewidth = africa_lw)

    # Add legend
    if ISO_i == 1
        Legend(fig[10, 1],
            [npc_line, access_line, use_line],
            ["NPC", "Access", "Use"])
    end
end

# Save fig
save(OUTPUT_PLOTS_DIR*"PaperFigures/Africa_country_ITN_Metrics.pdf", fig, pdf_version = "1.4")
fig

####################################################################
# %% FIGURE 3: ITN Utilisation Plot for each individual country
####################################################################
for ISO_i in ProgressBar(1:length(filt_ISOs))
    # Select ISO
    ISO = filt_ISOs[ISO_i]
    # Base axes
    fig = Figure(size = (800,500))
    ax = Axis(fig[1,1],
            title = "$(ISO) ITN Utilisation",
            xlabel = "Years", 
            xticks = (1:12:n_months, string.(YEAR_START:YEAR_END)),
            xticklabelrotation = pi/2,
            xlabelsize = 20,
            titlesize = 25,
            ylabel = "Net Utilisation (η)",
            yticks = (0:0.2:1),
            ylabelsize = 20
            )
    xlims!(-0.5,n_months+0.5)
    ylims!(-0.02, 1.02)

    # Add lines

    ### NPC
    # band!(ax, 1:n_months, country_utilisation[ISO_i,:,1], country_utilisation[ISO_i,:,3],
    #         color = (colors[1], fillalpha)
    #         )
    lines!(ax, 1:n_months, country_utilisation[ISO_i,:,2],
            color = colors[1], linewidth = lw)
    # lines!(ax, 1:n_months, country_utilisation[ISO_i,:,1],
    #         color = (colors[1], la), 
    #         linewidth = lw, linestyle = :dash)
    # lines!(ax, 1:n_months, country_utilisation[ISO_i,:,3],
    #         color = (colors[1], la), linewidth = lw, 
    #         linestyle = :dash)
    
    save(OUTPUT_PLOTS_DIR*"PaperFigures/Country_ITN_Metrics/$(ISO)_ITN_Utilisation.pdf", fig, pdf_version = "1.4")
end
fig

####################################################################
# %% FIGURE 4: ITN Utilisation Plot arranged according to Africa country location
####################################################################
# Make figure
layout_res = (1500,1300)
fig = Figure(size = layout_res)

plot_axs = []
title_axs = []
for ISO_i in 1:length(filt_ISOs)
    ISO = filt_ISOs[ISO_i]
    x_idx, y_idx, country_name = subplot_lookup[ISO]
    ax = Axis(fig[2*x_idx,y_idx], width = layout_res[1]/14,
                # xlabel = "Years", 
                xticks = ((1:12:n_months)[1:5:end], year_strings[1:5:end]),
                xticklabelrotation = pi/2,
                # ylabel = "Metric",
                yticks = (0:0.2:1),
                ylabelsize = 20)
    xlims!(ax,-0.5,n_months+0.5)
    ylims!(ax,-0.02, 1.02)

    lb = Label(fig[2*(x_idx-1)+1,y_idx], "$(country_name)")
    
    # Hide labels as required
    if !(ISO ∈ show_xlabel_ISOs)
        ax.xticklabelsvisible = false
    end

    if !(ISO ∈ show_ylabel_ISOs)
        ax.yticklabelsvisible = false
    end
    push!(plot_axs, ax)
    push!(title_axs, lb)
end

# Super Labels
Label(fig[:,0],"Utilisation Rate (η)", rotation = pi/2, fontsize = 40)
Label(fig[15,:],"Years", fontsize = 40)

# Add plot lines for metric
for ISO_i in ProgressBar(1:length(filt_ISOs), leave = false)
    # Select ISO
    ISO = filt_ISOs[ISO_i]

    # Add lines

    ### NPC
    npc_line = lines!(plot_axs[ISO_i], 1:n_months, country_utilisation[ISO_i,:,2],
            color = colors[1], linewidth = africa_lw)
end

save(OUTPUT_PLOTS_DIR*"PaperFigures/Africa_country_ITN_utilisation.pdf", fig, pdf_version = "1.4")

####################################################################
# %% FIGURE 5: Moving Window Gains plot
####################################################################
# %%

for ISO_i in ProgressBar(1:length(filt_ISOs))
    ISO = filt_ISOs[ISO_i]

    # Mean downsample to annual data
    country_crop_annual = [mean(country_crop[ISO_i,((i-1)*12+1):(12*i),2]) for i in 1:n_years]
    country_pop_annual = [mean(country_population[ISO_i,((i-1)*12+1):(12*i)]) for i in 1:n_years]

    # Calculate Moving Averages
    ma_crop = MA_filter(country_crop_annual, window = 3)
    ma_pop = MA_filter(country_pop_annual, window = 3)

    # Calculate percentage gains
    delta_crop = ma_crop[2:end]./ma_crop[1:end-1]
    delta_pop = ma_pop[2:end]./ma_pop[1:end-1]
    # 
    fig = Figure(size = (1000,400))

    ax = Axis(fig[1,1], title = "$(ISO) Population & ITN Trends",
                xlabel = "Year",
                xticks = (1:2:n_years, (string.(YEAR_START:YEAR_END))[1:2:end]),
                xticklabelrotation = pi/2,
                ylabel = "3 Year Δ%",
                yticks = 80:10:200,
                xlabelsize = 20, ylabelsize = 20)
    xlims!(ax, 0, n_years + 1)
    ylims!(ax, 75, 210)

    # Plot barplot
    data_tbl = (year_cat = repeat(1:length(delta_crop),2),
                    vals = vcat(delta_crop, delta_pop),
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
    # Save fig
    save(OUTPUT_PLOTS_DIR*"PaperFigures/Country_ITN_Trends/$(ISO)_ITN_trend.pdf", fig, pdf_version = "1.4")
end

####################################################################
# %% FIGURE 6: Mean Net Age Plot for each country
####################################################################

for ISO_i in ProgressBar(1:length(filt_ISOs))
    ISO = filt_ISOs[ISO_i]

    # Get required data slice from net age extract
    data_slice = netage_timeseries[  netage_timeseries.ISO .== ISO .&&
                        netage_timeseries.category .== "Admin0",:]

    # Extract timeseries and arrange as array
    netage = Vector{Float64}(undef, n_months)
    for monthidx in 1:n_months
        month, year_idx = monthidx_to_monthyear(monthidx)
        year = year_idx + YEAR_START - 1

        netage[monthidx] = data_slice[data_slice.month .== month .&&
                                        data_slice.year .== year, "mean_age_months"][1]
    end
    
    # Base axes
    fig = Figure(size = (800,500))
    ax = Axis(fig[1,1],
            title = "$(ISO) Net Age",
            xlabel = "Years", 
            xticks = (1:12:n_months, string.(YEAR_START:YEAR_END)),
            xticklabelrotation = pi/2,
            xlabelsize = 20,
            titlesize = 25,
            ylabel = "Age (Years)",
            yticks = (0:0.5:4.5),
            ylabelsize = 20
            )
    xlims!(-0.5,n_months+0.5)
    ylims!(-0.05, 4.55)

    # Add lines
    lines!(ax, 1:n_months, netage./12,
            color = colors[1], linewidth = lw)
    save(OUTPUT_PLOTS_DIR*"PaperFigures/Country_ITN_Metrics/$(ISO)_ITN_NetAge.pdf", fig, pdf_version = "1.4")
end


####################################################################
# %% FIGURE 7: ITN Net Age Plot arranged according to Africa country location
####################################################################
# Make figure
layout_res = (1500,1300)
fig = Figure(size = layout_res)

plot_axs = []
title_axs = []
for ISO_i in 1:length(filt_ISOs)
    ISO = filt_ISOs[ISO_i]
    x_idx, y_idx, country_name = subplot_lookup[ISO]
    ax = Axis(fig[2*x_idx,y_idx], width = layout_res[1]/14,
                xticks = ((1:12:n_months)[1:5:end], year_strings[1:5:end]),
                xticklabelrotation = pi/2,
                yticks = (0:1:4.5),
                ylabelsize = 20)
    xlims!(ax,-0.5,n_months+0.5)
    ylims!(-0.05, 4.55)

    lb = Label(fig[2*(x_idx-1)+1,y_idx], "$(country_name)")
    
    # Hide labels as required
    if !(ISO ∈ show_xlabel_ISOs)
        ax.xticklabelsvisible = false
    end

    if !(ISO ∈ show_ylabel_ISOs)
        ax.yticklabelsvisible = false
    end
    push!(plot_axs, ax)
    push!(title_axs, lb)
end

# Super Labels
Label(fig[:,0],"Mean Net Age (Years)", rotation = pi/2, fontsize = 40,
        padding = (-20,20,0,0))
Label(fig[15,:],"Years", fontsize = 40)

# Add plot lines for metric
for ISO_i in ProgressBar(1:length(filt_ISOs), leave = false)
    # Select ISO
    ISO = filt_ISOs[ISO_i]

    # Get required data slice from net age extract
    data_slice = netage_timeseries[  netage_timeseries.ISO .== ISO .&&
                        netage_timeseries.category .== "Admin0",:]

    # Extract timeseries and arrange as array
    netage = Vector{Float64}(undef, n_months)
    for monthidx in 1:n_months
        month, year_idx = monthidx_to_monthyear(monthidx)
        year = year_idx + YEAR_START - 1

        netage[monthidx] = data_slice[data_slice.month .== month .&&
                                        data_slice.year .== year, "mean_age_months"][1]
    end

    # Add lines

    ### NPC
    npc_line = lines!(plot_axs[ISO_i], 1:n_months, netage./12,
            color = colors[1], linewidth = africa_lw)
end

fig
# %%
save(OUTPUT_PLOTS_DIR*"PaperFigures/Africa_country_ITN_netage.pdf", fig, pdf_version = "1.4")
