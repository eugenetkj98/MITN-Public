"""
Author: Eugene Tan
Date Created: 10/4/2025
Last Updated: 10/4/2025
Script to calculate and compare rmse of BV and MITN and make associated plots
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
using CairoMakie
using StatsBase
using DateConversions

# %% Fig save dir
output_dir = OUTPUT_PLOTS_DIR*"raster summaries/"

# %% Import data lookup dataset
data = CSV.read(OUTPUT_DIR*"coverage_timeseries/bv_mitn_validation_values.csv", DataFrame)

# %% Filter all NaNs
filt_data = data[findall(.!isnan.(data.bv_npc) .&& .!isnan.(data.bv_access) .&& .!isnan.(data.bv_use) .&& 
                    .!isnan.(data.mitn_npc) .&& .!isnan.(data.mitn_access) .&& .!isnan.(data.mitn_use)),:]

# %% Calculate RMSE
bv_npc_rmse = sqrt(mean((filt_data.bv_npc .- filt_data.npc).^2))
bv_access_rmse = sqrt(mean((filt_data.bv_access .- filt_data.access).^2))
bv_use_rmse = sqrt(mean((filt_data.bv_use .- filt_data.use).^2))

mitn_npc_rmse = sqrt(mean((filt_data.mitn_npc .- filt_data.npc).^2))
mitn_access_rmse = sqrt(mean((filt_data.mitn_access .- filt_data.access).^2))
mitn_use_rmse = sqrt(mean((filt_data.mitn_use .- filt_data.use).^2))

# %% Calculate Residuals
bv_npc_res = mean(filt_data.bv_npc .- filt_data.npc)
bv_access_res = mean(filt_data.bv_access .- filt_data.access)
bv_use_res = mean(filt_data.bv_use .- filt_data.use)

mitn_npc_res = mean(filt_data.mitn_npc .- filt_data.npc)
mitn_access_res = mean(filt_data.mitn_access .- filt_data.access)
mitn_use_res = mean(filt_data.mitn_use .- filt_data.use)

# %% Calculate monthidxs
monthidxs = monthyear_to_monthidx.(filt_data.interview_month, filt_data.interview_year, YEAR_START = YEAR_NAT_START)

# %% Plot visual settings
# pythonplot()
# theme(:vibrant)
fillalpha = 0.15
la = 0.07
lw = 1.2
ms = 2
ma = 0.1
colors = [  colorant"#00976A", # NPC
            colorant"#E72A3D", # Access
            colorant"#0082C7", # Use
            colorant"#45332C", # Guidelines
            ]
colorrange = (YEAR_NAT_START, YEAR_NAT_END)
            

##############################
# %% RMSE Fit Plots using CairoMakie.jl
##############################
set_theme!(theme_ggplot2())
fig = Figure(size = (800,520))
BV_suptitle = fig[1,1:4] = GridLayout()
BV_plots = fig[2,1:4] = GridLayout()
MITN_suptitle = fig[3,1:4] = GridLayout()
MITN_plots = fig[4,1:4] = GridLayout()

# BV Title
Label(BV_suptitle[1,1:3],"BV Model", fontsize = 20, padding = (0,0,-20,-10))

# BV Plots
ax_bv_npc = Axis(BV_plots[1,1],
                    title = "NPC\nRMSE = $(round(bv_npc_rmse, digits = 4))", 
                    xlabel = "Cluster Observation", ylabel = "Model Prediction")
xlims!(ax_bv_npc, -0.05, 1.05)
ylims!(ax_bv_npc, -0.05, 1.05)
scatter!(ax_bv_npc, 
        filt_data.npc, filt_data.bv_npc, 
        color = (monthidxs./12).+YEAR_NAT_START,
        colorrange = colorrange, colormap = cgrad([colorant"#2A2A2A", colors[3]]),
        markersize = ms, alpha = ma)
lines!(ax_bv_npc, [0,1],[0,1], color = :black)

ax_bv_access = Axis(BV_plots[1,2],
                    title = "Access\nRMSE = $(round(bv_access_rmse, digits = 4))", 
                    xlabel = "Cluster Observation", ylabel = "Model Prediction")
xlims!(ax_bv_access, -0.05, 1.05)
ylims!(ax_bv_access, -0.05, 1.05)
scatter!(ax_bv_access, 
        filt_data.access, filt_data.bv_access, 
        color = (monthidxs./12).+YEAR_NAT_START,
        colorrange = colorrange, colormap = cgrad([colorant"#2A2A2A", colors[3]]),
        markersize = ms, alpha = ma)
lines!(ax_bv_access, [0,1], [0,1], color = :black)

ax_bv_use = Axis(BV_plots[1,3],
                    title = "Use\nRMSE = $(round(bv_use_rmse, digits = 4))", 
                    xlabel = "Cluster Observation", ylabel = "Model Prediction")
xlims!(ax_bv_use, -0.05, 1.05)
ylims!(ax_bv_use, -0.05, 1.05)
scatter!(ax_bv_use, 
        filt_data.use, filt_data.bv_use, 
        color = (monthidxs./12).+YEAR_NAT_START,
        colorrange = colorrange, colormap = cgrad([colorant"#2A2A2A", colors[3]]),
        markersize = ms, alpha = ma)
lines!(ax_bv_use, [0,1], [0,1], color = :black)

Colorbar(BV_plots[1,4], limits = colorrange, label = "Year",
        colormap = cgrad([colorant"#2A2A2A", colors[3]]))

# MITN Title
Label(MITN_suptitle[1,1:3],"MITN Model", fontsize = 20, padding = (0,0,-20,-10))

# MITN Plots

ax_mitn_npc = Axis(MITN_plots[1,1],
                    title = "NPC\nRMSE = $(round(mitn_npc_rmse, digits = 4))", 
                    xlabel = "Cluster Observation", ylabel = "Model Prediction")
xlims!(ax_mitn_npc, -0.05, 1.05)
ylims!(ax_mitn_npc, -0.05, 1.05)
scatter!(ax_mitn_npc, 
        filt_data.npc, filt_data.mitn_npc, 
        color = (monthidxs./12).+YEAR_NAT_START,
        colorrange = colorrange, colormap = cgrad([colorant"#2A2A2A", colors[2]]),
        markersize = ms, alpha = ma)
lines!(ax_mitn_npc, [0,1],[0,1], color = :black)

ax_mitn_access = Axis(MITN_plots[1,2],
                    title = "Access\nRMSE = $(round(mitn_access_rmse, digits = 4))", 
                    xlabel = "Cluster Observation", ylabel = "Model Prediction")
xlims!(ax_mitn_access, -0.05, 1.05)
ylims!(ax_mitn_access, -0.05, 1.05)
scatter!(ax_mitn_access, 
        filt_data.access, filt_data.mitn_access, 
        color = (monthidxs./12).+YEAR_NAT_START,
        colorrange = colorrange, colormap = cgrad([colorant"#2A2A2A", colors[2]]),
        markersize = ms, alpha = ma)
lines!(ax_mitn_access, [0,1], [0,1], color = :black)

ax_mitn_use = Axis(MITN_plots[1,3],
                    title = "Use\nRMSE = $(round(mitn_use_rmse, digits = 4))", 
                    xlabel = "Cluster Observation", ylabel = "Model Prediction")
xlims!(ax_mitn_use, -0.05, 1.05)
ylims!(ax_mitn_use, -0.05, 1.05)
scatter!(ax_mitn_use, 
        filt_data.use, filt_data.mitn_use, 
        color = (monthidxs./12).+YEAR_NAT_START,
        colorrange = colorrange, colormap = cgrad([colorant"#2A2A2A", colors[2]]),
        markersize = ms, alpha = ma)
lines!(ax_mitn_use, [0,1], [0,1], color = :black)

Colorbar(MITN_plots[1,4], limits = colorrange, label = "Year",
        colormap = cgrad([colorant"#2A2A2A", colors[2]]))

# Adjust Gaps

colgap!(BV_plots, 10)
rowgap!(BV_plots, 10)

colgap!(MITN_plots, 10)
rowgap!(MITN_plots, 10)

# Visualise Figure


fig


# %% Save figure
save(output_dir*"bv_mitn_rmse.pdf", fig)

##############################
# %% Residuals Fit Plots using CairoMakie.jl
##############################
set_theme!(theme_ggplot2())
fig = Figure(size = (800,520))
BV_suptitle = fig[1,1:4] = GridLayout()
BV_plots = fig[2,1:4] = GridLayout()
MITN_suptitle = fig[3,1:4] = GridLayout()
MITN_plots = fig[4,1:4] = GridLayout()

# BV Title
Label(BV_suptitle[1,1:3],"BV Model", fontsize = 20, padding = (0,0,-20,-10))

# BV Plots
ax_bv_npc = Axis(BV_plots[1,1],
                    title = "NPC\nMean Residual = $(round(bv_npc_res, digits = 4))", 
                    xlabel = "Cluster Observation", ylabel = "Residual")
xlims!(ax_bv_npc, -0.05, 1.05)
ylims!(ax_bv_npc, -0.55, 0.55)
scatter!(ax_bv_npc, 
        filt_data.npc, filt_data.bv_npc .- filt_data.npc, 
        color = (monthidxs./12).+YEAR_NAT_START,
        colorrange = colorrange, colormap = cgrad([colorant"#2A2A2A", colors[3]]),
        markersize = ms, alpha = ma)
lines!(ax_bv_npc, [0,1],[0,0], color = :black)

ax_bv_access = Axis(BV_plots[1,2],
                    title = "Access\nMean Residual = $(round(bv_access_res, digits = 4))", 
                    xlabel = "Cluster Observation", ylabel = "Residual")
xlims!(ax_bv_access, -0.05, 1.05)
ylims!(ax_bv_access, -0.55, 0.55)
scatter!(ax_bv_access, 
        filt_data.access, filt_data.bv_access .- filt_data.access, 
        color = (monthidxs./12).+YEAR_NAT_START,
        colorrange = colorrange, colormap = cgrad([colorant"#2A2A2A", colors[3]]),
        markersize = ms, alpha = ma)
lines!(ax_bv_access, [0,1], [0,0], color = :black)

ax_bv_use = Axis(BV_plots[1,3],
                    title = "Use\nMean Residual = $(round(bv_use_res, digits = 4))", 
                    xlabel = "Cluster Observation", ylabel = "Residual")
xlims!(ax_bv_use, -0.05, 1.05)
ylims!(ax_bv_use, -0.55, 0.55)
scatter!(ax_bv_use, 
        filt_data.use, filt_data.bv_use .- filt_data.use, 
        color = (monthidxs./12).+YEAR_NAT_START,
        colorrange = colorrange, colormap = cgrad([colorant"#2A2A2A", colors[3]]),
        markersize = ms, alpha = ma)
lines!(ax_bv_use, [0,1], [0,0], color = :black)

Colorbar(BV_plots[1,4], limits = colorrange, label = "Year",
        colormap = cgrad([colorant"#2A2A2A", colors[3]]))

# MITN Title
Label(MITN_suptitle[1,1:3],"MITN Model", fontsize = 20, padding = (0,0,-20,-10))

# MITN Plots

ax_mitn_npc = Axis(MITN_plots[1,1],
                    title = "NPC\nMean Residual = $(round(mitn_npc_res, digits = 4))", 
                    xlabel = "Cluster Observation", ylabel = "Residual")
xlims!(ax_mitn_npc, -0.05, 1.05)
ylims!(ax_mitn_npc, -0.55, 0.55)
scatter!(ax_mitn_npc, 
        filt_data.npc, filt_data.mitn_npc .- filt_data.npc, 
        color = (monthidxs./12).+YEAR_NAT_START,
        colorrange = colorrange, colormap = cgrad([colorant"#2A2A2A", colors[2]]),
        markersize = ms, alpha = ma)
lines!(ax_mitn_npc, [0,1],[0,0], color = :black)

ax_mitn_access = Axis(MITN_plots[1,2],
                    title = "Access\nMean Residual = $(round(mitn_access_res, digits = 4))", 
                    xlabel = "Cluster Observation", ylabel = "Residual")
xlims!(ax_mitn_access, -0.05, 1.05)
ylims!(ax_mitn_access, -0.55, 0.55)
scatter!(ax_mitn_access, 
        filt_data.access, filt_data.mitn_access .- filt_data.access, 
        color = (monthidxs./12).+YEAR_NAT_START,
        colorrange = colorrange, colormap = cgrad([colorant"#2A2A2A", colors[2]]),
        markersize = ms, alpha = ma)
lines!(ax_mitn_access, [0,1], [0,0], color = :black)

ax_mitn_use = Axis(MITN_plots[1,3],
                    title = "Use\nMean Residual = $(round(mitn_use_res, digits = 4))", 
                    xlabel = "Cluster Observation", ylabel = "Residual")
xlims!(ax_mitn_use, -0.05, 1.05)
ylims!(ax_mitn_use, -0.55, 0.55)
scatter!(ax_mitn_use, 
        filt_data.use, filt_data.mitn_use .- filt_data.use, 
        color = (monthidxs./12).+YEAR_NAT_START,
        colorrange = colorrange, colormap = cgrad([colorant"#2A2A2A", colors[2]]),
        markersize = ms, alpha = ma)
lines!(ax_mitn_use, [0,1], [0,0], color = :black)

Colorbar(MITN_plots[1,4], limits = colorrange, label = "Year",
        colormap = cgrad([colorant"#2A2A2A", colors[2]]))

# Adjust Gaps

colgap!(BV_plots, 10)
rowgap!(BV_plots, 10)

colgap!(MITN_plots, 10)
rowgap!(MITN_plots, 10)

# Visualise Figure


fig


# %% Save figure
save(output_dir*"bv_mitn_residuals.pdf", fig)

# ##############################
# # %% Plots using Plots.jl
# ##############################

# # %% Make plots of RMSE and Error 

# fig_1 = scatter(filt_data.npc[1:100:end], filt_data.bv_npc[1:100:end], zcolor = (monthidxs[1:100:end]./12).+YEAR_NAT_START,
#                 markersize = ms, markeralpha = ma, cmap = cgrad([colorant"#2A2A2A", colors[3]]),
#                 label = nothing,
#                 xlabel = "Cluster Observation",
#                 ylabel = "Model Prediction",
#                 colorbar_title = "Year",
#                 colorbar_ticks = YEAR_NAT_START:2:YEAR_NAT_END,
#                 title = "BV NPC\nRMSE = $(round(bv_npc_rmse, digits = 4))", 
#                 xlims = (-0.05, 1.05), ylims = (-0.05, 1.05))
# plot!(fig_1, 0:1, 0:1, linecolor = :black, label = nothing)

# fig_2 = scatter(filt_data.access, filt_data.bv_access, zcolor = (monthidxs./12).+YEAR_NAT_START,
#                 markersize = ms, markeralpha = ma,  cmap = cgrad([colorant"#2A2A2A", colors[3]]),
#                 label = nothing,
#                 xlabel = "Cluster Observation",
#                 ylabel = "Model Prediction",
#                 colorbar_title = "Year",
#                 title = "BV Access\nRMSE = $(round(bv_access_rmse, digits = 4))", 
#                 xlims = (-0.05, 1.05), ylims = (-0.05, 1.05))
# plot!(fig_2, 0:1, 0:1, linecolor = :black, label = nothing)

# fig_3 = scatter(filt_data.use, filt_data.bv_use, zcolor = (monthidxs./12).+YEAR_NAT_START,
#                 markersize = ms, markeralpha = ma,  cmap = cgrad([colorant"#2A2A2A", colors[3]]),
#                 label = nothing,
#                 xlabel = "Cluster Observation",
#                 ylabel = "Model Prediction",
#                 colorbar_title = "Year",
#                 title = "BV Use\nRMSE = $(round(bv_use_rmse, digits = 4))", 
#                 xlims = (-0.05, 1.05), y = (-0.05, 1.05))
# plot!(fig_3, 0:1, 0:1, linecolor = :black, label = nothing)

# fig_4 = scatter(filt_data.npc, filt_data.mitn_npc, zcolor = (monthidxs./12).+YEAR_NAT_START,
#                 markersize = ms, markeralpha = ma,  cmap = cgrad([colorant"#2A2A2A", colors[2]]),
#                 label = nothing,
#                 xlabel = "Cluster Observation",
#                 ylabel = "Model Prediction",
#                 colorbar_title = "Year",
#                 title = "MITN NPC\nRMSE = $(round(mitn_npc_rmse, digits = 4))", 
#                 xlims = (-0.05, 1.05), ylims = (-0.05, 1.05))
# plot!(fig_4, 0:1, 0:1, linecolor = :black, label = nothing)

# fig_5 = scatter(filt_data.access, filt_data.mitn_access, zcolor = (monthidxs./12).+YEAR_NAT_START,
#                 markersize = ms, markeralpha = ma,  cmap = cgrad([colorant"#2A2A2A", colors[2]]),
#                 label = nothing,
#                 xlabel = "Cluster Observation",
#                 ylabel = "Model Prediction",
#                 colorbar_title = "Year",
#                 title = "MITN Access\nRMSE = $(round(mitn_access_rmse, digits = 4))", 
#                 xlims = (-0.05, 1.05), ylims = (-0.05, 1.05))
# plot!(fig_5, 0:1, 0:1, linecolor = :black, label = nothing)

# fig_6 = scatter(filt_data.use, filt_data.mitn_use, zcolor = (monthidxs./12).+YEAR_NAT_START,
#                 markersize = ms, markeralpha = ma,  cmap = cgrad([colorant"#2A2A2A", colors[2]]),
#                 label = nothing,
#                 xlabel = "Cluster Observation",
#                 ylabel = "Model Prediction",
#                 title = "MITN Use\nRMSE = $(round(mitn_use_rmse, digits = 4))", 
#                 xlims = (-0.05, 1.05), ylims = (-0.05, 1.05))
# plot!(fig_6, 0:1, 0:1, linecolor = :black, label = nothing)

# rmse_plot = plot(fig_1, fig_2, fig_3, fig_4, fig_5, fig_6,
#         layout = (2,3), size = (800,400))

# savefig(rmse_plot, output_dir*"bv_mitn_rmse.pdf")