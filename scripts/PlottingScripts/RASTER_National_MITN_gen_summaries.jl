"""
Author: Eugene Tan
Date Created: 10/4/2025
Last Updated: 10/4/2025
Script to make time series plots from full ITN model outputs
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
using ProgressMeter
using DateConversions
using Plots

# %% Dataset Directories
input_dir = OUTPUT_DIR*"coverage_timeseries/"
input_filename = "master_extraction.csv"

# %% Plot save directories
output_dir = OUTPUT_PLOTS_DIR*"raster summaries/"
mkpath(output_dir*"countries/")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Get list of ISOs to analyse
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = ["CPV","ZAF"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Load Time Series data and raw survey data (for comparison)
data = CSV.read(input_dir*input_filename, DataFrame)
cluster_obs = CSV.read(OUTPUT_DATAPREP_DIR*"INLA/inla_dataset_reduced.csv", DataFrame)

# %% Plot visual settings
pythonplot()
theme(:ggplot2)
fillalpha = 0.15
la = 0.07
lw = 1.2
ms = 3
ma = 0.45
colors = [  colorant"#00976A", # NPC
            colorant"#E72A3D", # Access
            colorant"#0082C7", # Use
            colorant"#45332C", # Guidelines
            ];

# %% Storage variable for plots
NPC_deviation_plot_collection = []
Access_deviation_plot_collection = []
ITN_coverage_plot_collection = []
npc_access_collection = []
npc_use_collection = []
access_use_collection = []

# %% Construct plots for each country
@showprogress for ISO in filt_ISOs
    # %% Load Data
    # Filter data
    filt_data = data[(data.ISO .== ISO) .* (data.category .== "Admin0"),:]

    # Get month and year vals, to calculate monthidxs for sorting and alignment
    month_vals = filt_data.month
    year_vals = filt_data.year
    monthidxs = monthyear_to_monthidx.(month_vals, year_vals, YEAR_START = YEAR_NAT_START)
    YEAR_START, YEAR_END = minimum(year_vals), maximum(year_vals)

    # Sort filtered data according to monthidxs
    sorted_data = filt_data[sortperm(monthidxs),:]
    sorted_monthidxs = monthidxs[sortperm(monthidxs)]

    # Extract required metric values
    snf_npc = sorted_data[:,["snf_npc_95lower", "snf_npc_mean", "snf_npc_95upper"]]
    raster_npc = sorted_data[:,["raster_npc_95lower", "raster_npc_mean", "raster_npc_95upper"]]

    snf_access = sorted_data[:,["snf_access_95lower", "snf_access_mean", "snf_access_95upper"]]
    raster_access = sorted_data[:,["raster_access_95lower", "raster_access_mean", "raster_access_95upper"]]

    raster_use = sorted_data[:,["raster_use_95lower", "raster_use_mean", "raster_use_95upper"]]

    # %% Make plots
    # Construct base plots
    1:12:((YEAR_NAT_END-YEAR_NAT_START+1)*12)
    
    
    fig_npc = plot(title = "NPC Model Dev. \n$(ISO)",
                xticks = (1:12:(((YEAR_NAT_END-YEAR_NAT_START+1)*12)+1), YEAR_NAT_START:(YEAR_NAT_END+1)), 
                xtickfontrotation = 90,
                xlims = (0, 288),
                ylims = (-0.01, 1.01),
                xlabel = "Years",
                ylabel = "NPC",
                legend = :topleft)
    fig_access = plot(title = "Access Model Dev. \n$(ISO)",
                xticks = (1:12:(((YEAR_NAT_END-YEAR_NAT_START+1)*12)+1), YEAR_NAT_START:(YEAR_NAT_END+1)), 
                xtickfontrotation = 90,
                xlims = (0, 288),
                ylims = (-0.01, 1.01),
                xlabel = "Years",
                ylabel = "Access",
                legend = :topleft)
    fig_comb = plot(title = "ITN Coverage \n$(ISO)",
                xticks = (1:12:(((YEAR_NAT_END-YEAR_NAT_START+1)*12)+1), YEAR_NAT_START:(YEAR_NAT_END+1)), 
                xtickfontrotation = 90,
                xlims = (0, 288),
                ylims = (-0.01, 1.01),
                xlabel = "Years",
                ylabel = "Coverage Metric",
                legend = :topleft)

    # %% NPC plot
    plot!(fig_npc, 
            snf_npc[:,1], fillrange = snf_npc[:,3],
            linealpha = 0, fillalpha = fillalpha,
            color = colors[1], fillcolor = colors[1],
            label = nothing)
    plot!(fig_npc, 
            snf_npc[:,2],
            linealpha = 1,
            color = colors[1],
            linestyle = :solid,
            label = "SNF")
    plot!(fig_npc, 
            raster_npc[:,1], fillrange = raster_npc[:,3],
            linealpha = 0, fillalpha = fillalpha,
            color = colors[2], fillcolor = colors[2],
            label = nothing)
    plot!(fig_npc, 
            raster_npc[:,2],
            linealpha = 1,
            color = colors[2],
            linestyle = :solid,
            label = "Raster")

    # %% Access Plot
    plot!(fig_access, 
            snf_access[:,2],
            linealpha = 1,
            color = colors[1],
            linestyle = :solid,
            label = "SNF")
    plot!(fig_access, 
            raster_access[:,1], fillrange = raster_access[:,3],
            linealpha = 0, fillalpha = fillalpha,
            color = colors[2], fillcolor = colors[2],
            label = nothing)
    plot!(fig_access, 
            raster_access[:,2],
            linealpha = 1,
            color = colors[2],
            linestyle = :solid,
            label = "Raster")

    # %% Combined Plot
    plot!(fig_comb, 
            snf_npc[:,1], fillrange = snf_npc[:,3],
            linealpha = 0, fillalpha = fillalpha,
            color = colors[1], fillcolor = colors[1],
            label = nothing)
    plot!(fig_comb, 
            snf_npc[:,2],
            linealpha = 1,
            color = colors[1],
            linestyle = :solid,
            label = "NPC (SNF)")
    plot!(fig_comb, 
            raster_npc[:,1], fillrange = raster_npc[:,3],
            linealpha = 0, fillalpha = fillalpha,
            color = colors[1], fillcolor = colors[1],
            label = nothing)
    plot!(fig_comb, 
            raster_npc[:,2],
            linealpha = 1,
            color = colors[1],
            linestyle = :dash,
            label = "NPC (Raster)")

    plot!(fig_comb, 
            snf_access[:,2],
            linealpha = 1,
            color = colors[2],
            linestyle = :solid,
            label = "Access (SNF)")

    plot!(fig_comb, 
            raster_access[:,2],
            linealpha = 1,
            color = colors[2],
            linestyle = :dash,
            label = "Access (Raster)")

    plot!(fig_comb, 
            raster_use[:,1], fillrange = raster_use[:,3],
            linealpha = 0, fillalpha = fillalpha,
            color = colors[3], fillcolor = colors[3],
            label = nothing)
    plot!(fig_comb, 
            raster_use[:,2],
            linealpha = 1,
            color = colors[3],
            linestyle = :dash,
            label = "Use (Raster)")

    # %% Plot NPC, Access and Use Relationships
    fig_1 = plot(title = "NPC - Access\n$(ISO)",
                    xlabel = "NPC",
                    ylabel = "Access",
                    legend = :topleft,
                    xticks = (0:0.1:1),
                    xlims = (-0.01,1.01))

    scatter!(fig_1, snf_npc[:,2], snf_access[:,2],
                    markersize = ms,
                    color = colors[2],
                    label = "SNF")
    scatter!(fig_1, raster_npc[:,2], raster_access[:,2],
                    markersize = ms,
                    color = colors[3],
                    xlabel = "NPC",
                    ylabel = "Access",
                    label = "Raster")
    plot!(fig_1, 0:1, 0:1, color = :black, label = nothing, linestyle = :dash)

    fig_2 = plot(title = "NPC - Use\n$(ISO)",
                    xlabel = "NPC",
                    ylabel = "Use",
                    legend = :topleft,
                    xticks = (0:0.1:1),
                    xlims = (-0.01,1.01))

    scatter!(fig_2, snf_npc[:,2], raster_use[:,2],
                    markersize = ms,
                    color = colors[2],
                    label = "SNF")
    scatter!(fig_2, raster_npc[:,2], raster_use[:,2],
                    markersize = ms,
                    color = colors[3],
                    label = "Raster")
	plot!(fig_2, 0:1, 0:1, color = :black, label = nothing, linestyle = :dash)


    fig_3 = plot(title = "Access - Use\n$(ISO)",
                    xlabel = "Access",
                    ylabel = "Use",
                    legend = :topleft,
                    xticks = (0:0.1:1),
                    xlims = (-0.01,1.01))

    scatter!(fig_3, snf_access[:,2], raster_use[:,2],
                    markersize = ms,
                    color = colors[2],
                    label = "SNF")
    scatter!(fig_3, raster_access[:,2], raster_use[:,2],
                    markersize = ms,
                    color = colors[3],
                    label = "Raster")
	plot!(fig_3, 0:1, 0:1, color = :black, label = nothing, linestyle = :dash)

    # %% Make combined plot
    combined_iso_fig = plot(fig_npc, fig_access, fig_comb,
                                fig_1, fig_2, fig_3, layout = (2,3), size = (1000,600))

    # %% Save combined plot
    savefig(combined_iso_fig, output_dir*"countries/$(ISO)_summary_plot.pdf")

    # %% Add component plot to collection
    push!(NPC_deviation_plot_collection, fig_npc)
    push!(Access_deviation_plot_collection, fig_access)
    push!(ITN_coverage_plot_collection, fig_comb)
    push!(npc_access_collection, fig_1)
    push!(npc_use_collection, fig_2)
    push!(access_use_collection, fig_3)
end

# %% Save aggregated plots for all countries
# Plot layout settings
layout = (8,6)
figsize = (2560,2000)

# Make plots and save
fig_npc_comb = plot(NPC_deviation_plot_collection..., layout = layout, size = figsize)
fig_access_comb = plot(Access_deviation_plot_collection..., layout = layout, size = figsize)
fig_itn_coverage_comb = plot(ITN_coverage_plot_collection..., layout = layout, size = figsize)

fig_npc_access_comb = plot(npc_access_collection..., layout = layout, size = figsize)
fig_npc_use_comb = plot(npc_use_collection..., layout = layout, size = figsize)
fig_access_use_comb = plot(access_use_collection..., layout = layout, size = figsize)

savefig(fig_npc_comb, output_dir*"npc_model_diff.pdf")
savefig(fig_access_comb, output_dir*"access_model_diff.pdf")
savefig(fig_itn_coverage_comb, output_dir*"itn_coverage_summary.pdf")
savefig(fig_npc_access_comb, output_dir*"npc_access_country.pdf")
savefig(fig_npc_use_comb, output_dir*"npc_use_country.pdf")
savefig(fig_access_use_comb, output_dir*"access_use_country.pdf")

# %% Make combined NPC-Access-Use Relationship plot

# Pre process data again
# Filter data to just Admin0 (Country level)
filt_data = data[data.category .== "Admin0",:]

# Extract required metric values
snf_npc = filt_data[:,["snf_npc_95lower", "snf_npc_mean", "snf_npc_95upper"]]
raster_npc = filt_data[:,["raster_npc_95lower", "raster_npc_mean", "raster_npc_95upper"]]

snf_access = filt_data[:,["snf_access_95lower", "snf_access_mean", "snf_access_95upper"]]
raster_access = filt_data[:,["raster_access_95lower", "raster_access_mean", "raster_access_95upper"]]

raster_use = filt_data[:,["raster_use_95lower", "raster_use_mean", "raster_use_95upper"]]


# %% Plot NPC, Access and Use Relationships
cluster_obs = CSV.read("outputs/data_prep/subnational/subnat_npc_monthly_data.csv", DataFrame)
fig_1 = plot(title = "NPC - Access",
                xlabel = "NPC",
                ylabel = "Access",
                legend = :topleft,
                xticks = (0:0.1:1),
                xlims = (-0.05,1.05),
                ylims = (-0.05,1.05))
fig_2 = plot(title = "NPC - Use",
                xlabel = "NPC",
                ylabel = "Use",
                legend = :topleft,
                xlims = (-0.05,1.05),
                ylims = (-0.05,1.05))
fig_3 = plot(title = "Access - Use",
                xlabel = "Access",
                ylabel = "Use",
                legend = :topleft,
                xlims = (-0.05,1.05),
                ylims = (-0.05,1.05))

scatter!(fig_1, snf_npc[:,2], snf_access[:,2],
                markersize = ms*0.3, markeralpha = ma,
                color = colors[2],
                label = "SNF")
scatter!(fig_1, raster_npc[:,2], raster_access[:,2],
                markersize = ms*0.3, markeralpha = ma,
                color = colors[3],
                xlabel = "NPC",
                ylabel = "Access",
                label = "Raster")
plot!(fig_1, 0:1, 0:1, color = :black, label = nothing, linestyle = :dash)


scatter!(fig_2, snf_npc[:,2], raster_use[:,2],
                markersize = ms*0.3, markeralpha = ma,
                color = colors[2],
                label = "SNF")
scatter!(fig_2, raster_npc[:,2], raster_use[:,2],
                markersize = ms*0.3, markeralpha = ma,
                color = colors[3],
                label = "Raster")
plot!(fig_2, 0:1, 0:1, color = :black, label = nothing, linestyle = :dash)


scatter!(fig_3, snf_access[:,2], raster_use[:,2],
                markersize = ms*0.3, markeralpha = ma,
                color = colors[2],
                label = "SNF")
scatter!(fig_3, raster_access[:,2], raster_use[:,2],
                markersize = ms*0.3, markeralpha = ma,
                color = colors[3],
                label = "Raster")
plot!(fig_3, 0:1, 0:1, color = :black, label = nothing, linestyle = :dash)

combined_iso_fig = plot(fig_1, fig_2, fig_3, layout = (1,3), size = (1200,350));

# %% Save combined plot
savefig(combined_iso_fig, output_dir*"global_metric_relationships.pdf")