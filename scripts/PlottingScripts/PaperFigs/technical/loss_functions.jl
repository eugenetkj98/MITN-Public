"""
Author: Eugene Tan
Date Created: 6/5/2025
Last Updated: 6/5/2025
Plots to show loss functions for compact exponential and weibull (for illustrative purposes)
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/read_toml.jl")

# %% Math packages
using LinearAlgebra
using StatsBase

# %% Plotting Packages
using CairoMakie

# %% General Packages
using ProgressBars
using LaTeXStrings

# %% Custom Packages
using NetLoss

# %% Generate Curves
t_max = 4
t_vals = 0:0.01:4

# Compact Exponential Parameters
τ_vals = [2,3,4]
κ_vals = [0.4,6,25]

# Weibull Parameters
b_vals = [2,2,2,1]
k_vals = [0.8,4,20,4]

# Storage variables
compact_life = Matrix{Float64}(undef, length(κ_vals), length(t_vals))
weibull_life = Matrix{Float64}(undef, length(b_vals), length(t_vals))

# Generate survival curves
for i in 1:length(κ_vals)
    compact_life[i,:] = net_loss_compact.(t_vals, Float64(τ_vals[i]), Float64(κ_vals[i]))
end

for i in 1:length(k_vals)
    weibull_life[i,:] = net_loss_weibull.(t_vals, Float64(b_vals[i]), Float64(k_vals[i]))
end

# %% Basic Plot settings
set_theme!(theme_latexfonts())
ls = 17
ts = 20
lw = 2

# %% Make Plot
fig = Figure(size = (800,280))
ax1 = Axis(fig[1,1],
            title = "Compact Exponential",
            xlabel = L"t",
            ylabel = L"f(t)",
            xlabelsize = ls, ylabelsize = ls,
            titlesize = ts)
ax2 = Axis(fig[1,2],
            title = "Weibull",
            xlabel = L"t",
            ylabel = L"f(t)",
            xlabelsize = ls, ylabelsize = ls,
            titlesize = ts)

compact_elems = []
weibull_elems = []

for i in 1:length(κ_vals)
    push!(compact_elems, lines!(ax1, t_vals, compact_life[i,:],
            linewidth = lw,
            label = L"\tau = %$(τ_vals[i]),\, \kappa = %$(κ_vals[i])"))
end

for i in 1:length(k_vals)
    push!(weibull_elems, lines!(ax2, t_vals, weibull_life[i,:],
            linewidth = lw,
            label = L"b = %$(b_vals[i]),\, k = %$(k_vals[i])"))
end

Legend(fig[1,1], ax1,
        labelsize = ls, halign = :right, valign = :top,
        tellwidth = false, tellheight = false,
        margin = (5,5,5,5), rowgap = -2)
Legend(fig[1,2], ax2,
        labelsize = ls, halign = :right, valign = :top,
        tellwidth = false, tellheight = false,
        margin = (5,5,5,5), rowgap = -2)

save(OUTPUT_PLOTS_DIR*"PaperFigures/TechnicalPaper/Loss_Functions.pdf", fig, pdf_version = "1.4")
fig
