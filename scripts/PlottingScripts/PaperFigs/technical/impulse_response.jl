"""
Author: Eugene Tan
Date Created: 6/5/2025
Last Updated: 6/5/2025
Mathematical analysis to demonstrate the primary effects of net attrition function on overall system response.
Analysis is for technical paper
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Math packages
using LinearAlgebra
using StatsBase

# %% Plotting Packages
using CairoMakie

# %% General Packages
using ProgressBars
using LaTeXStrings

# %% Custom Packages
using Convolutions
using NetLoss

# %% Basic Plot settings
set_theme!(theme_latexfonts())


######################################
# %% Scenario 1: Uniform Annual Distribution
######################################
# %% Define loss/net attrition function
loss_t_vals = 0:1/12:5
τ = 6
κ = 40
L_base = net_loss_compact.(loss_t_vals, τ, κ)

# Nets distribution behaviour
n_years = 10
t_vals = 0:(1/12):(n_years)

# Once a year distribution
d = zeros(length(t_vals)).= 1/12

# # Plots.jl plot
# fig_attrition = plot(loss_t_vals, L_base, 
#                     title = "Net Attrition Curve\n" * 
#                             L"""$\tau = %$(τ), \kappa = %$(κ)$""",
#                     xlabel = "Years", ylabel = "Survival Rate",
#                     legend = nothing,
#                     xticks = 0:1:10)
                    
# fig_crop = plot(t_vals, convolution(L_base,d)[1:length(t_vals)],
#                     xlabel = "Years", ylabel = "Net Crop "*L"(\Gamma)",
#                     title = "Reduced Model Net Crop",
#                     xticks = 0:n_years)

# plot(fig_attrition, fig_crop, 
#         layout = (1,2), size = (800,300), margins = 5mm)


# Combine into a single plot
fig = Figure(size = (800,300))
ax1 = Axis(fig[1,1], 
                title = L"$\text{\tau} = %$(τ), \, \text{\kappa} = %$(κ)$",
                titlesize = 20,
                xlabel = "Years", xlabelsize = 18,
                ylabel = "Net Survival Rate", ylabelsize = 18,
                xticks = 0:1:10)
ax2 = Axis(fig[1,2], 
                title = L"$\text{\tau} = %$(τ), \, \text{\kappa} = %$(κ)$",
                titlesize = 20,
                xlabel = "Years", xlabelsize = 18,
                ylabel = "Net Crop (Γ)",ylabelsize = 18,
                xticks = 0:1:10)
lines!(ax1, loss_t_vals, L_base)
lines!(ax2, t_vals, convolution(L_base,d)[1:length(t_vals)])
fig

# %% Do a parameter scan

# Time series bounds
loss_t_vals = 0:1/12:5
n_years = 10
t_vals = 0:(1/12):(n_years)

# Uniform distribution behaviour
d = zeros(length(t_vals)).= 1/12

τ_vals = 0.1:0.01:5
κ_vals = exp.(-4:0.1:4)

SS_val = zeros(length(τ_vals), length(κ_vals))

for i in 1:length(τ_vals)
    for j in 1:length(κ_vals)
        τ = τ_vals[i]
        κ = κ_vals[j]
        L = net_loss_compact.(loss_t_vals, τ, κ)
        SS_val[i,j] = maximum(convolution(L,d)[1:length(t_vals)])
    end
end

# %%
# heatmap(τ_vals, κ_vals, SS_val',
#         xlabel = L"\tau", ylabel = L"\kappa",
#         colorbartitle = L"\Gamma_{SS}",
#         labelfontsize = 15, colorbar_titlefontsize = 15,
#         title = "Uniform Annual Distribution")

fig = Figure(size = (500,400))
ax = Axis(fig[1,1],
        title = L"\text{Uniform Annual Distribution}",
        titlesize = 25,
        xlabel = L"\tau", xlabelsize = 20,
        ylabel = L"\kappa", ylabelsize = 20,
        yscale = log,
        yticks = round.(exp.(-4:1:4), digits = 2)
        )
clims = (0,4)
Colorbar(fig[1,2], 
        colorrange = clims, colormap = :viridis, 
        label = L"\Gamma_{SS}", labelsize = 20)
heatmap!(ax, τ_vals, κ_vals, SS_val, colorrange = clims)
fig

######################################
# %% Scenario 2a: Periodic Distribution, varying τ
######################################
# Nets distribution behaviour
n_years = 10
t_vals = 0:(1/12):(n_years)

# Smooth periodic distribution
d = zeros(length(t_vals))
ϕ = 0.2
ω = 1
k = 2
d_shape = abs.(sin.(π .* (ω .* t_vals .- ϕ))).^k
d = ((n_years)/sum(d_shape)) .* d_shape

# %% Define loss/net attrition function
loss_t_vals = 0:1/12:10
τ_vals = [3,8,12,20]
κ_vals = [18,18,18,18]

L_base_vals = []
conv_vals = []


for i in 1:length(τ_vals)
    push!(L_base_vals, net_loss_compact.(loss_t_vals, τ_vals[i], κ_vals[i]))
end

L_base_cat = hcat(L_base_vals...)



# Combine into a single plot
fig = Figure(size = (1200,300))
ax1 = Axis(fig[1,1], 
                title = L"\text{Net Attrition Curve}",
                titlesize = 22,
                xlabel = "Years",  xlabelsize = 18,
                ylabel = "Survival Rate", ylabelsize = 18,
                xticks = 0:1:10)
ax2 = Axis(fig[1,2], 
                title = L"\text{Distribution Time Series}",
                titlesize = 22,
                xlabel = "Years",  xlabelsize = 18,
                ylabel = "Normalised Distribution", ylabelsize = 18,
                xticks = 0:n_years)
ax3 = Axis(fig[1,3], 
                title = L"\text{Effective Net Crop}",
                titlesize = 22,
                xlabel = "Years",  xlabelsize = 18, 
                ylabel = L"\text{Net Crop }(\Gamma)", ylabelsize = 18,
                xticks = 0:n_years)

survival_lines = []
conv_lines = []
for i in 1:length(τ_vals)
    push!(survival_lines, lines!(ax1, loss_t_vals, L_base_cat[:,i]))
    push!(conv_lines, lines!(ax3, t_vals, convolution(L_base_cat[:,i],d)[1:length(t_vals)]))
end
Legend(fig[1,4], survival_lines,
        # tellwidth = false, tellheight = false,
        labelsize = 18,
        halign = :right, valign = :center,
        padding = (5,5,5,5),
        margin = (5,5,5,5),
        [L"$\tau = %$(τ_vals[i]), \,\kappa = %$(round(κ_vals[i], digits = 2))$" for i in 1:length(τ_vals)])
# Legend(fig[1,3], conv_lines,
#         tellwidth = false, tellheight = false,
#         labelsize = 18,
#         halign = :right, valign = :bottom,
#         padding = (5,5,5,5),
#         margin = (5,5,5,5),
#         [L"$\tau = %$(τ_vals[i]), \,\kappa = %$(round(κ_vals[i], digits = 2))$" for i in 1:length(τ_vals)])

lines!(ax2, t_vals, d)
fig

######################################
# %% Scenario 2b: Periodic Distribution, varying κ
######################################
# Nets distribution behaviour
n_years = 10
t_vals = 0:(1/12):(n_years)

# Smooth periodic distribution
d = zeros(length(t_vals))
ϕ = 0.2
ω = 0.5
k = 2
d_shape = abs.(sin.(π .* (ω .* t_vals .- ϕ))).^k
d = ((n_years)/sum(d_shape)) .* d_shape

# %% Define loss/net attrition function
loss_t_vals = 0:1/12:5
τ_vals = [10,10,10,10]
κ_vals = [5,10,20,25]

L_base_vals = []
conv_vals = []


for i in 1:length(τ_vals)
    push!(L_base_vals, net_loss_compact.(loss_t_vals, τ_vals[i], κ_vals[i]))
end

L_base_cat = hcat(L_base_vals...)



# Combine into a single plot
fig = Figure(size = (1200,300))
ax1 = Axis(fig[1,1], 
                title = L"\text{Net Attrition Curve}",
                titlesize = 22,
                xlabel = "Years",  xlabelsize = 18,
                ylabel = "Survival Rate", ylabelsize = 18,
                xticks = 0:1:10)
ax2 = Axis(fig[1,2], 
                title = L"\text{Distribution Time Series}",
                titlesize = 22,
                xlabel = "Years",  xlabelsize = 18,
                ylabel = "Normalised Distribution", ylabelsize = 18,
                xticks = 0:n_years)
ax3 = Axis(fig[1,3], 
                title = L"\text{Effective Net Crop}",
                titlesize = 22,
                xlabel = "Years",  xlabelsize = 18, 
                ylabel = L"\text{Net Crop }(\Gamma)", ylabelsize = 18,
                xticks = 0:n_years)

survival_lines = []
conv_lines = []
for i in 1:length(τ_vals)
    push!(survival_lines, lines!(ax1, loss_t_vals, L_base_cat[:,i]))
    push!(conv_lines, lines!(ax3, t_vals, convolution(L_base_cat[:,i],d)[1:length(t_vals)]))
end
Legend(fig[1,4], survival_lines,
        # tellwidth = false, tellheight = false,
        labelsize = 18,
        halign = :right, valign = :center,
        padding = (5,5,5,5),
        margin = (5,5,5,5),
        [L"$\tau = %$(τ_vals[i]), \,\kappa = %$(round(κ_vals[i], digits = 2))$" for i in 1:length(τ_vals)])
# Legend(fig[1,3], conv_lines,
#         tellwidth = false, tellheight = false,
#         labelsize = 18,
#         halign = :right, valign = :bottom,
#         padding = (5,5,5,5),
#         margin = (5,5,5,5),
#         [L"$\tau = %$(τ_vals[i]), \,\kappa = %$(round(κ_vals[i], digits = 2))$" for i in 1:length(τ_vals)])

lines!(ax2, t_vals, d)
fig


######################################
# %% Figure 3: Parameter scan for periodic distribution
######################################
# %% Do a parameter scan
# Time series bounds
subsample_ratio = 36
time_resolution = 1/subsample_ratio
loss_t_vals = 0:time_resolution:5
n_years = 40
t_vals = 0:(time_resolution):(n_years)
# Periodic distribution behaviour
ϕ = 0
k = 3
ω_vals = 0.6:0.2:2
τ_vals = 1:0.25:20
κ_vals = exp.(-4:0.1:4)

MAX_val = zeros(length(ω_vals), length(τ_vals), length(κ_vals))
STD_val = zeros(length(ω_vals), length(τ_vals), length(κ_vals))

for ω_i in ProgressBar(1:length(ω_vals), leave = false)
    Threads.@threads for τ_i in ProgressBar(1:length(τ_vals), leave = false)
        for κ_i in 1:length(κ_vals)
            ω = ω_vals[ω_i]
            τ = τ_vals[τ_i]
            κ = κ_vals[κ_i]

            # Construct distribution time series
            d = zeros(length(t_vals))
            d_shape = abs.(sin.(π .* (ω .* t_vals .- ϕ))).^k
            d = ((n_years)/sum(d_shape)) .* d_shape

            # Construct response
            L = net_loss_compact.(loss_t_vals, τ, κ)
            response = convolution(L,d)[1:length(t_vals)]

            # Store Variables
            MAX_val[ω_i, τ_i, κ_i] = maximum(response)
            STD_val[ω_i, τ_i, κ_i] = std(response[end-20*subsample_ratio:end])
        end
    end
end

# %% Maximum Net Crop
fig = Figure(size = (900,450))
layout = (2,4)
idx_matrix = collect(Iterators.product(1:layout[1],1:layout[2]))
idx_matrix_transpose = Matrix{Tuple{Int64, Int64}}(undef, size(idx_matrix)[2], size(idx_matrix)[1])
for i in 1:size(idx_matrix)[1]
    for j in 1:size(idx_matrix)[2]
        idx_matrix_transpose[j,i] = idx_matrix[i,j]
    end
end
fig_idxs = idx_matrix_transpose[:]

axs = [Axis(fig[fig_idxs[i][1], fig_idxs[i][2]],
            xlabel = L"\tau", xlabelsize = 22,
            ylabel = L"\kappa", ylabelsize = 22,
            title = L"\omega = %$(ω_vals[i])", titlesize = 25,
            yscale = log, yticks = round.(exp.(-4:1:4), digits = 2)) for i in 1:length(ω_vals)]
clims = (0,6)
for ω_i in 1:length(ω_vals)
    heatmap!(axs[ω_i], τ_vals, κ_vals, MAX_val[ω_i, :,:], colorrange = clims)
    contour!(axs[ω_i], τ_vals, κ_vals, MAX_val[ω_i, :,:], 
        levels = round.(LinRange(0,maximum(MAX_val[ω_i, :,:]),5), digits = 2), 
        labels = true, color = (:black, 0.6))
end
Colorbar(fig[:,layout[2]+1], colorrange = clims, colormap = :viridis, 
        label = L"\Gamma_{MAX}", labelsize = 20)
fig

save(OUTPUT_PLOTS_DIR*"PaperFigures/analytical_periodic_max_heatmap.pdf", fig)

# %% Net Crop Std
fig = Figure(size = (900,450))
layout = (2,4)
idx_matrix = collect(Iterators.product(1:layout[1],1:layout[2]))
idx_matrix_transpose = Matrix{Tuple{Int64, Int64}}(undef, size(idx_matrix)[2], size(idx_matrix)[1])
for i in 1:size(idx_matrix)[1]
    for j in 1:size(idx_matrix)[2]
        idx_matrix_transpose[j,i] = idx_matrix[i,j]
    end
end
fig_idxs = idx_matrix_transpose[:]

axs = [Axis(fig[fig_idxs[i][1], fig_idxs[i][2]],
            xlabel = L"\tau", xlabelsize = 22,
            ylabel = L"\kappa", ylabelsize = 22,
            title = L"\omega = %$(ω_vals[i])", titlesize = 25,
            yscale = log, yticks = round.(exp.(-4:1:4), digits = 2)) for i in 1:length(ω_vals)]
clims = (0,0.4)
for ω_i in 1:length(ω_vals)
    
    heatmap!(axs[ω_i], τ_vals, κ_vals, STD_val[ω_i, :,:], colorrange = clims)
    contour!(axs[ω_i], τ_vals, κ_vals, STD_val[ω_i, :,:], 
            levels = round.(LinRange(0,maximum(STD_val[ω_i, :,:]),3), digits = 2), 
            labels = true, color = (:black, 0.6))
end

Colorbar(fig[:,layout[2]+1], colorrange = clims, colormap = :viridis, 
        label = L"\Gamma_{STD}", labelsize = 20)
fig

save(OUTPUT_PLOTS_DIR*"PaperFigures/analytical_periodic_std_heatmap.pdf", fig)