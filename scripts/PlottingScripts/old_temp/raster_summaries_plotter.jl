"""
Author: Eugene Tan
Date Created: 28/1/2025
Last Updated: 28/1/2025
Makie script to generate summary figures of NPC, Access and Use maps from rasters
"""
# %% Load packages
using Rasters
using Plots

# %% Declare directories
raster_dir = "outputs/rasters/" # Input directory to read rasters from
plot_output_dir = "output_plots/maps/summaries/"
gif_output_dir = "output_plots/maps/renamed_summaries/"

# %% Time bounds for extractions
YEAR_START = 2000
YEAR_END = 2002

# %% Script to generate list of summaries
"""
Generates 2 copies of figures. The first copy is saved to plot_output_dir and is
named according to year-month. The second copy is saved to "output_plots/maps/renamed summaries/"
for the purposes of collating into a gif (i.e. renamed in the format 000001.png etc...)
"""
k = 0 # Declare counter for frames
for (month, year) in ProgressBar(collect(Iterators.product(1:12,YEAR_START:YEAR_END))[:])
    # Get correct year/month string for importing file
    year_str = "$(year)"
    month_str = "$(month)"
    if month < 10
        month_str = "0"*month_str
    end

    # Import Rasters
    combined_npc_snf_mean_raster = Raster(raster_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_mean.tif")
    combined_npc_snf_upper_raster = Raster(raster_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_upper.tif")
    combined_npc_snf_lower_raster = Raster(raster_dir*"final_npc/snf_npc/npc_$(year)_$(month_str)_lower.tif")

    combined_npc_map_mean_raster = Raster(raster_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_mean.tif")
    combined_npc_map_upper_raster = Raster(raster_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_upper.tif")
    combined_npc_map_lower_raster = Raster(raster_dir*"final_npc/logmodel_npc/npc_$(year)_$(month_str)_lower.tif")

    combined_npc_adj_map_mean_raster = Raster(raster_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_mean.tif")
    combined_npc_adj_map_upper_raster = Raster(raster_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_upper.tif")
    combined_npc_adj_map_lower_raster = Raster(raster_dir*"final_npc/logmodel_npc/adj_npc_$(year)_$(month_str)_lower.tif")

    combined_access_snf_mean_raster = Raster(raster_dir*"final_access/snf_access/access_$(year)_$(month_str)_mean.tif")
    combined_access_snf_upper_raster = Raster(raster_dir*"final_access/snf_access/access_$(year)_$(month_str)_upper.tif")
    combined_access_snf_lower_raster = Raster(raster_dir*"final_access/snf_access/access_$(year)_$(month_str)_lower.tif")

    combined_access_map_mean_raster = Raster(raster_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_mean.tif")
    combined_access_map_upper_raster = Raster(raster_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_upper.tif")
    combined_access_map_lower_raster = Raster(raster_dir*"final_access/pmodel_access/access_$(year)_$(month_str)_lower.tif")

    combined_access_adj_map_mean_raster = Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_mean.tif")
    combined_access_adj_map_upper_raster = Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_upper.tif")
    combined_access_adj_map_lower_raster = Raster(raster_dir*"final_access/pmodel_access/adj_access_$(year)_$(month_str)_lower.tif")

    combined_use_map_mean_raster = Raster(raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_mean.tif")
    combined_use_map_upper_raster = Raster(raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_upper.tif")
    combined_use_map_lower_raster = Raster(raster_dir*"final_use/logis_use/use_$(year)_$(month_str)_lower.tif")

    # combined_use_adj_map_mean_raster = Raster(raster_dir*"final_use/logis_use/adj_use_$(year)_$(month_str)_mean.tif")
    # combined_use_adj_map_upper_raster = Raster(raster_dir*"final_use/logis_use/adj_use_$(year)_$(month_str)_upper.tif")
    # combined_use_adj_map_lower_raster = Raster(raster_dir*"final_use/logis_use/adj_use_$(year)_$(month_str)_lower.tif")

    # Plot settings
    cbar_range = (0,1)
    npc_cmap = :thermal
    access_cmap = :viridis
    use_cmap = :devon
    npc_cbar_range = cbar_range
    access_cbar_range = cbar_range
    use_cbar_range = cbar_range
    suptitle_fs = 38
    sup_fs = 28
    subtitlesize = 20
    cbarsize = 30
    cbarfontsize = 25
    labelsize = 30
    
    # Construct figure layouts
    fig = Figure(size = (2000,800))
    Label(fig[0, 1:7], "ITN Coverage Maps ($(year)-$month)", valign = :bottom,
                        font = :bold, fontsize = suptitle_fs)
    grid_npc = fig[1,1:3] = GridLayout()
    grid_access = fig[1,4:6] = GridLayout()
    grid_use = fig[1,7] = GridLayout()
    

    # SNF NPC
    ax11 = grid_npc[1,1] = Axis(fig, 
                        title = "SNF", ylabel = "95% Lower",
                        titlesize = subtitlesize,
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax11)
    hideydecorations!(ax11, label= false)
    
    ax21 = grid_npc[2,1] = Axis(fig, 
                        ylabel = "Mean",
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax21)
    hideydecorations!(ax21, label= false)

    ax31 = grid_npc[3,1] = Axis(fig, 
                        ylabel = "95% Upper",
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax31)
    hideydecorations!(ax31, label= false)

    # Spatial NPC
    ax12 = grid_npc[1,2] = Axis(fig, 
                        title = "Spatial", ylabel = "95% Lower",
                        titlesize = subtitlesize,
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax12)
    hideydecorations!(ax12)
    
    ax22 = grid_npc[2,2] = Axis(fig, 
                         ylabel = "Mean",
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax22)
    hideydecorations!(ax22)

    ax32 = grid_npc[3,2] = Axis(fig, 
                        ylabel = "95% Upper",
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax32)
    hideydecorations!(ax32)

    # Adjusted Spatial NPC
    ax13 = grid_npc[1,3] = Axis(fig, 
                        title = "Adj. Spatial", ylabel = "95% Lower",
                        titlesize = subtitlesize,
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax13)
    hideydecorations!(ax13)
    
    ax23 = grid_npc[2,3] = Axis(fig, 
                         ylabel = "Mean",
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax23)
    hideydecorations!(ax23)

    ax33 = grid_npc[3,3] = Axis(fig, 
                        ylabel = "95% Upper",
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax33)
    hideydecorations!(ax33)

    Colorbar(grid_npc[1:3,4], limits = (0,1), colormap = npc_cmap,
                        size = cbarsize, ticklabelsize = cbarfontsize)

    Label(grid_npc[0, 1:4], "NPC", valign = :bottom,
                        font = :bold, fontsize = sup_fs)

    # SNF Access
    ax14 = grid_access[1,1] = Axis(fig, 
                        title = "SNF", ylabel = "95% Lower",
                        titlesize = subtitlesize,
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax14)
    hideydecorations!(ax14)
    
    ax24 = grid_access[2,1] = Axis(fig, 
                         ylabel = "Mean",
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax24)
    hideydecorations!(ax24)

    ax34 = grid_access[3,1] = Axis(fig, 
                        ylabel = "95% Upper",
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax34)
    hideydecorations!(ax34)
    
    # Spatial Access
    ax15 = grid_access[1,2] = Axis(fig, 
                        title = "Spatial", ylabel = "95% Lower",
                        titlesize = subtitlesize,
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax15)
    hideydecorations!(ax15)
    
    ax25 = grid_access[2,2] = Axis(fig, 
                         ylabel = "Mean",
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax25)
    hideydecorations!(ax25)

    ax35 = grid_access[3,2] = Axis(fig, 
                        ylabel = "95% Upper",
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax35)
    hideydecorations!(ax35)
    
    # Adjusted Spatial Access
    ax16 = grid_access[1,3] = Axis(fig, 
                        title = "Adj. Spatial", ylabel = "95% Lower",
                        titlesize = subtitlesize,
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax16)
    hideydecorations!(ax16)
    
    ax26 = grid_access[2,3] = Axis(fig, 
                         ylabel = "Mean",
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax26)
    hideydecorations!(ax26)

    ax36 = grid_access[3,3] = Axis(fig, 
                        ylabel = "95% Upper",
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax36)
    hideydecorations!(ax36)

    Colorbar(grid_access[1:3,4], limits = (0,1), colormap = access_cmap,
                        size = cbarsize, ticklabelsize = cbarfontsize)

    Label(grid_access[0, 1:4], "Access", valign = :bottom,
                        font = :bold, fontsize = sup_fs)

    # Spatial Use
    ax17 = grid_use[1,1] = Axis(fig, 
                        title = "Spatial", ylabel = "95% Lower",
                        titlesize = subtitlesize,
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax17)
    hideydecorations!(ax17)
    
    ax27 = grid_use[2,1] = Axis(fig, 
                         ylabel = "Mean",
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax27)
    hideydecorations!(ax27)

    ax37 = grid_use[3,1] = Axis(fig, 
                        ylabel = "95% Upper",
                        xlabelsize = labelsize, ylabelsize = labelsize)
    hidexdecorations!(ax37)
    hideydecorations!(ax37)
    
    # # Adjusted Spatial Access
    # ax18 = grid_use[1,2] = Axis(fig, 
    #                     title = "Adj. Spatial", ylabel = "95% Lower",
    #                     titlesize = subtitlesize,
    #                     xlabelsize = labelsize, ylabelsize = labelsize)
    # hidexdecorations!(ax18)
    # hideydecorations!(ax18)
    
    # ax28 = grid_use[2,2] = Axis(fig, 
    #                      ylabel = "Mean",
    #                     xlabelsize = labelsize, ylabelsize = labelsize)
    # hidexdecorations!(ax28)
    # hideydecorations!(ax28)

    # ax38 = grid_use[3,2] = Axis(fig, 
    #                     ylabel = "95% Upper",
    #                     xlabelsize = labelsize, ylabelsize = labelsize)
    # hidexdecorations!(ax38)
    # hideydecorations!(ax38)

    Colorbar(grid_use[1:3,2], limits = (0,1), colormap = use_cmap,
                        size = cbarsize, ticklabelsize = cbarfontsize)
    
    Label(grid_use[0, 1:2], "Use", valign = :bottom,
                        font = :bold, fontsize = sup_fs)

    # Plot components
    CairoMakie.plot!(ax11, combined_npc_snf_lower_raster, 
                        colormap = npc_cmap, colorrange = npc_cbar_range)
    CairoMakie.plot!(ax21, combined_npc_snf_mean_raster, 
                        colormap = npc_cmap, colorrange = npc_cbar_range)
    CairoMakie.plot!(ax31, combined_npc_snf_upper_raster, 
                        colormap = npc_cmap, colorrange = npc_cbar_range)

    CairoMakie.plot!(ax12, combined_npc_map_lower_raster, 
                        colormap = npc_cmap, colorrange = npc_cbar_range)
    CairoMakie.plot!(ax22, combined_npc_map_mean_raster, 
                        colormap = npc_cmap, colorrange = npc_cbar_range)
    CairoMakie.plot!(ax32, combined_npc_map_upper_raster, 
                        colormap = npc_cmap, colorrange = npc_cbar_range)

    CairoMakie.plot!(ax13, combined_npc_adj_map_lower_raster, 
                        colormap = npc_cmap, colorrange = npc_cbar_range)
    CairoMakie.plot!(ax23, combined_npc_adj_map_mean_raster, 
                        colormap = npc_cmap, colorrange = npc_cbar_range)
    CairoMakie.plot!(ax33, combined_npc_adj_map_upper_raster, 
                        colormap = npc_cmap, colorrange = npc_cbar_range)

    CairoMakie.plot!(ax14, combined_access_snf_lower_raster, 
                        colormap = access_cmap, colorrange = access_cbar_range)
    CairoMakie.plot!(ax24, combined_access_snf_mean_raster, 
                        colormap = access_cmap, colorrange = access_cbar_range)
    CairoMakie.plot!(ax34, combined_access_snf_upper_raster, 
                        colormap = access_cmap, colorrange = access_cbar_range)

    CairoMakie.plot!(ax15, combined_access_map_lower_raster, 
                        colormap = access_cmap, colorrange = access_cbar_range)
    CairoMakie.plot!(ax25, combined_access_map_mean_raster, 
                        colormap = access_cmap, colorrange = access_cbar_range)
    CairoMakie.plot!(ax35, combined_access_map_upper_raster, 
                        colormap = access_cmap, colorrange = access_cbar_range)

    CairoMakie.plot!(ax16, combined_access_adj_map_lower_raster, 
                        colormap = access_cmap, colorrange = access_cbar_range)
    CairoMakie.plot!(ax26, combined_access_adj_map_mean_raster, 
                        colormap = access_cmap, colorrange = access_cbar_range)
    CairoMakie.plot!(ax36, combined_access_adj_map_upper_raster, 
                        colormap = access_cmap, colorrange = access_cbar_range)

    CairoMakie.plot!(ax17, combined_use_map_lower_raster, 
                        colormap = use_cmap, colorrange = use_cbar_range)
    CairoMakie.plot!(ax27, combined_use_map_mean_raster, 
                        colormap = use_cmap, colorrange = use_cbar_range)
    CairoMakie.plot!(ax37, combined_use_map_upper_raster, 
                        colormap = use_cmap, colorrange = use_cbar_range)

    # CairoMakie.plot!(ax18, combined_use_adj_map_lower_raster, 
    #                     colormap = use_cmap, colorrange = use_cbar_range)
    # CairoMakie.plot!(ax28, combined_use_adj_map_mean_raster, 
    #                     colormap = use_cmap, colorrange = use_cbar_range)
    # CairoMakie.plot!(ax38, combined_use_adj_map_upper_raster, 
    #                     colormap = use_cmap, colorrange = use_cbar_range)

    # Save figure
    mkpath(plot_output_dir)
    CairoMakie.save(plot_output_dir*"ITN_map_$(year)_$(month).png", fig)

    # Save renamed version use to construct gif later
    mkpath(gif_output_dir)
    filename=lpad(k, 6, "0")*".png"
    CairoMakie.save(gif_output_dir*filename, fig)

    # Update counter
    k=k+1
end

# %% Construct gif animation using generated figures
import Plots: Animation, buildanimation
fnames = String[]
k = 0
for (month, year) in ProgressBar(collect(Iterators.product(1:12, YEAR_START:YEAR_END))[:])
    # Gif component names
    filename=lpad(k, 6, "0")*".png"
    push!(fnames, filename)
    k= k+1
end

# Build animation
anim = Animation("output_plots/maps/renamed_summaries", fnames)
buildanimation(anim, "output_plots/maps/renamed_summaries/map_summaries.gif", fps = 6, show_msg=false)