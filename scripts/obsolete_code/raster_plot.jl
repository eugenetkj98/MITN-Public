"""
This seems like code to plot an animation of the raster maps.
May be useful for quick figure prep, but not much more. 
Consider refactoring into a different code base.
"""


# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import relevant packages
using ProgressBars
using DataFrames
using CSV
using Rasters
using Shapefile
using GeoInterface
using GeoIO
using JLD2
using StatsBase

# Custom packages
using DateConversions

using Plots
using Measures

# %%

YEAR_START = 2000
YEAR_END = 2023
n_months = (YEAR_END-YEAR_START+1)*12

anim = @animate for i in ProgressBar(1:n_months, leave = false)
    month, year_ref = monthidx_to_monthyear(i)
    year = year_ref + YEAR_START -1

    npc_gap_raster = Raster("outputs/maps/NPC_gap_$(year)_$(month).tif")
    plot(npc_gap_raster, clims = (-0.3,0.3), title = "NPC Gap\n$(year)-$(month)",
        xlabel = "long", ylabel = "lat", colorbar_title = "NPC Gap", margins = 4mm)
end fps = 5 

gif(anim, "output_plots/npc_gap.mp4", fps = 15)

# %%

YEAR_START = 2000
YEAR_END = 2023

anim = @animate for year in ProgressBar(YEAR_START:YEAR_END, leave = false)

    npc_gap_raster = Raster("outputs/maps/INLA/maps/Access_gap_reducedmodel_$year.tif")
    plot(npc_gap_raster, clims = (-0.3,0.3), title = "Access Gap\n$(year)",
        xlabel = "long", ylabel = "lat", colorbar_title = "Access Gap", margins = 4mm)
end fps = 3

gif(anim, "output_plots/access_gap.gif", fps = 7)

# %%

YEAR_START = 2003
YEAR_END = 2010
n_months = (YEAR_END-YEAR_START+1)*12

@gif for i in ProgressBar(1:n_months, leave = false)
    month, year_ref = monthidx_to_monthyear(i)
    year = year_ref + YEAR_START -1

    year_str = "$(year)"
    if month < 10
        month_str = "0$(month)"
    else
        month_str = "$(month)"
    end

    cov_raster = Raster("Z:mastergrids/MODIS_Global/MCD43D6_v6_BRDF_Reflectance/TCB_v6/5km/Monthly/TCB_v6.$(year_str).$(month_str).max.5km.max.tif")
    plot(cov_raster, clims = (0,2), title = "NPC Gap\n$(year)-$(month)",
                xlims = (-20,60), ylims = (-35,35),
                xlabel = "long", ylabel = "lat", colorbar_title = "Cov", margins = 4mm)
end fps = 10

# %%



