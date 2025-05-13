"""
Author: Eugene Tan
Date Created: 24/4/2025
Last Updated: 24/4/2025
Make animation of raster plots over time period
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Import Public Packages
using GeoIO
using Rasters
using DateConversions
using Plots
using ProgressBars
using Measures

# %% Input Raster Directory
input_dir = OUTPUT_RASTERS_DIR

# %% Year bounds
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

# %% Create storage variables for rasters
month_vals = []
year_vals = []
snf_npc_rasters = []
snf_access_rasters = []
mitn_npc_rasters = []
mitn_access_rasters = []
mitn_use_rasters = []

# %%
for year in ProgressBar(YEAR_START:YEAR_END)
    for month in 1:12
        month_str = "$(month)"
        if month < 10
            month_str = "0$(month)"
        end

        snf_npc = Raster(input_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_mean.tif")
        snf_access = Raster(input_dir*"final_access/snf_access/access_$(year)_$(month_str)_mean.tif")
        mitn_npc = Raster(input_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_mean.tif")
        mitn_access = Raster(input_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_mean.tif")
        mitn_use = Raster(input_dir*"final_use/logis_use/use_$(year)_$(month_str)_mean.tif")

        push!(month_vals, month)
        push!(year_vals, year)
        push!(snf_npc_rasters, snf_npc)
        push!(snf_access_rasters, snf_access)
        push!(mitn_npc_rasters, mitn_npc)
        push!(mitn_access_rasters, mitn_access)
        push!(mitn_use_rasters, mitn_use)
    end
end

# %% Construct composite plot

composite_plot_collection = []
for i in ProgressBar(1:length(month_vals))
    fig_snf_npc = plot(snf_npc_rasters[i], clims = (0,1),
                        xlabel = "Longitude", ylabel = "Latitude",
                        # colorbar_title = "NPC",
                        title = "SNF NPC\n$(year_vals[i])-$(month_vals[i])\n",
                        titlefontsize = 16,
                        cmap = :navia)
    fig_access_npc = plot(snf_access_rasters[i], clims = (0,1),
                        xlabel = "Longitude", ylabel = "Latitude",
                        # colorbar_title = "Access",
                        title = "SNF Access\n$(year_vals[i])-$(month_vals[i])\n",
                        titlefontsize = 16,
                        cmap = :lajolla)
    fig_mitn_npc = plot(mitn_npc_rasters[i], clims = (0,1),
                        xlabel = "Longitude", ylabel = "Latitude",
                        # colorbar_title = "NPC",
                        title = "MITN NPC\n$(year_vals[i])-$(month_vals[i])\n",
                        titlefontsize = 16,
                        cmap = :navia)
    fig_mitn_access = plot(mitn_access_rasters[i], clims = (0,1),
                        xlabel = "Longitude", ylabel = "Latitude",
                        # colorbar_title = "Access",
                        title = "MITN Access\n$(year_vals[i])-$(month_vals[i])\n",
                        titlefontsize = 16,
                        cmap = :lajolla)
    fig_mitn_use = plot(mitn_use_rasters[i], clims = (0,1),
                        xlabel = "Longitude", ylabel = "Latitude",
                        # colorbar_title = "Use",
                        title = "MITN Use\n$(year_vals[i])-$(month_vals[i])\n",
                        titlefontsize = 16,
                        cmap = :devon)

    fig = plot(fig_snf_npc, fig_mitn_npc, fig_mitn_use, fig_access_npc, fig_mitn_access,
            layout = (2,3), size = (1800,950), 
            leftmargin = 8mm, rightmargin = -8mm,
            bottommargin = 6mm, topmargin = -10mm)
    push!(composite_plot_collection, fig)
end

# %% Plot animation
anim = @animate for i in ProgressBar(1:length(composite_plot_collection))
    plot(composite_plot_collection[i])
end
gif(anim, "output_plots/ITN_map_animation.gif", fps = 10)