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
year_vals = []
mitn_npc_rasters = []
mitn_access_rasters = []
mitn_use_rasters = []

# %%
for year in ProgressBar(YEAR_START:YEAR_END)
    mitn_npc = Raster(input_dir*"final_npc/annual/npc_$(year)_mean.tif")
    mitn_access = Raster(input_dir*"final_access/annual/access_$(year)_mean.tif")
    mitn_use = Raster(input_dir*"final_use/annual/use_$(year)_mean.tif")

    push!(year_vals, year)
    push!(mitn_npc_rasters, mitn_npc)
    push!(mitn_access_rasters, mitn_access)
    push!(mitn_use_rasters, mitn_use)
end

# %% Construct composite plot
composite_plot_collection = []
for i in ProgressBar(1:length(year_vals))
    fig_mitn_npc = plot(mitn_npc_rasters[i], clims = (0,1),
                        xlabel = "Longitude", ylabel = "Latitude",
                        # colorbar_title = "NPC",
                        title = "MITN NPC\n$(year_vals[i])\n",
                        titlefontsize = 16,
                        cmap = :navia)
    fig_mitn_access = plot(mitn_access_rasters[i], clims = (0,1),
                        xlabel = "Longitude", ylabel = "Latitude",
                        # colorbar_title = "Access",
                        title = "MITN Access\n$(year_vals[i])\n",
                        titlefontsize = 16,
                        cmap = :lajolla)
    fig_mitn_use = plot(mitn_use_rasters[i], clims = (0,1),
                        xlabel = "Longitude", ylabel = "Latitude",
                        # colorbar_title = "Use",
                        title = "MITN Use\n$(year_vals[i])\n",
                        titlefontsize = 16,
                        cmap = :devon)

    fig = plot(fig_mitn_npc, fig_mitn_access, fig_mitn_use, 
            layout = (1,3), size = (1800,500), 
            leftmargin = 8mm, rightmargin = -8mm,
            bottommargin = 6mm, topmargin = -10mm)
    push!(composite_plot_collection, fig)
end

# %% Plot animation
anim = @animate for i in ProgressBar(1:length(composite_plot_collection))
    plot(composite_plot_collection[i])
end
gif(anim, "output_plots/ITN_map_annual_animation.gif", fps = 3)


