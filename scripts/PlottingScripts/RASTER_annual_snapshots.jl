
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Import Public Packages
using GeoIO
using Rasters
using DateConversions
using CairoMakie
using ProgressBars
using Measures

# %% Input Raster Directory
input_dir = OUTPUT_RASTERS_DIR

# %% Year bounds
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

# %% Create storage variables for rasters
year_vals = [2000,2004,2008,2012,2016,2020,2023]

# %% Create Figure and axes
fig = Figure(size = (3000,1200))
npc_axes = [Axis(fig[1,2*(i-1)+1], title = "NPC\n$(year_vals[i])", titlesize = 30) for i in 1:length(year_vals)]
access_axes = [Axis(fig[2,2*(i-1)+1], title = "Access\n$(year_vals[i])", titlesize = 30) for i in 1:length(year_vals)]
use_axes = [Axis(fig[3,2*(i-1)+1], title = "Use\n$(year_vals[i])", titlesize = 30) for i in 1:length(year_vals)]
for i in 1:length(year_vals)
    Colorbar(fig[1,2*(i-1)+2], limits = (0,1), colormap = :navia, size = 30)
    Colorbar(fig[2,2*(i-1)+2], limits = (0,1), colormap = :lajolla, size = 30)
    Colorbar(fig[3,2*(i-1)+2], limits = (0,1), colormap = :devon, size = 30)
end

# %%
for year_i in ProgressBar(1:length(year_vals), leave = false)
    year = year_vals[year_i]
    mitn_npc = Raster(input_dir*"final_npc/annual/npc_$(year)_mean.tif")
    mitn_access = Raster(input_dir*"final_access/annual/access_$(year)_mean.tif")
    mitn_use = Raster(input_dir*"final_use/annual/use_$(year)_mean.tif")

    plot!(npc_axes[year_i], mitn_npc, colormap = :navia, colorrange = (0,1))
    plot!(access_axes[year_i], mitn_access, colormap = :lajolla, colorrange = (0,1))
    plot!(use_axes[year_i], mitn_use, colormap = :devon, colorrange = (0,1))
end

save(OUTPUT_PLOTS_DIR*"raster summaries/annual_snapshots.png",fig)
