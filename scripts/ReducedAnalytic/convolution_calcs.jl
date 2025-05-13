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

# %% Required packages
using Plots
using ProgressBars
using Measures
using LinearAlgebra
using Convolutions
using NetLoss
using LaTeXStrings
using StatsBase

# %% Basic Plot settings
# set_theme!(theme_latexfonts())
plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing, grid=false)

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

# Plots.jl plot
fig_attrition = plot(loss_t_vals, L_base, 
                    title = "Net Attrition Curve\n" * 
                            L"""$\tau = %$(τ), \kappa = %$(κ)$""",
                    xlabel = "Years", ylabel = "Survival Rate",
                    legend = nothing,
                    xticks = 0:1:10)
                    
fig_crop = plot(t_vals, convolution(L_base,d)[1:length(t_vals)],
                    xlabel = "Years", ylabel = "Net Crop "*L"(\Gamma)",
                    title = "Reduced Model Net Crop",
                    xticks = 0:n_years)

plot(fig_attrition, fig_crop, 
        layout = (1,2), size = (800,300), margins = 5mm)


# # Combine into a single plot
# fig = Figure(size = (800,300))
# ax1 = Axis(fig[1,1], 
#                 title = #L"Net Attrition Curve\n" * 
#                 L"$\tau = %$(τ), \kappa = %$(κ)$",
#                 xlabel = "Years", ylabel = "Survival Rate",
#                 xticks = 0:1:10)
# ax2 = Axis(fig[1,2], 
#                 title = #"Net Attrition Curve\n" * 
#                 L"""$\tau = %$(τ), \kappa = %$(κ)$""",
#                 xlabel = "Years", ylabel = "Survival Rate",
#                 xticks = 0:1:10)
# lines!(ax1, loss_t_vals, L_base)
# lines!(ax2, t_vals, convolution(L_base,d)[1:length(t_vals)])
# fig

# %% Do a parameter scan

# Time series bounds
loss_t_vals = 0:1/12:5
n_years = 10
t_vals = 0:(1/12):(n_years)

# Uniform distribution behaviour
d = zeros(length(t_vals)).= 1/12

τ_vals = 0.1:0.01:5
κ_vals = 11:1:40

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
heatmap(τ_vals, κ_vals, SS_val',
        xlabel = L"\tau", ylabel = L"\kappa",
        colorbartitle = L"\Gamma_{SS}",
        labelfontsize = 15, colorbar_titlefontsize = 15,
        title = "Uniform Annual Distribution")

# fig = Figure(size = (500,400))
# ax = Axis(fig[1,1],
#         xlabel = L"\tau", ylabel = L"\kappa",
#         xlabelsize = 20, ylabelsize = 20,
#         title = "Uniform Annual Distribution")
# heatmap!(ax, τ_vals, κ_vals, SS_val)
# fig 
######################################
# %% Scenario 2: Periodic Distribution
######################################
# %% Define loss/net attrition function
loss_t_vals = 0:1/12:5
τ = 1.4
κ = exp(-1)

L_base = net_loss_compact.(loss_t_vals, τ, κ)


# Nets distribution behaviour
n_years = 10
t_vals = 0:(1/12):(n_years)

# Smooth periodic distribution
d = zeros(length(t_vals))
ϕ = 0.2
ω = 1.5
k = 2
d_shape = abs.(sin.(π .* (ω .* t_vals .- ϕ))).^k
d = ((n_years)/sum(d_shape)) .* d_shape



# Calculate net crop response and plot

# Plots.jl Combine into a single plot
fig_attrition = plot(loss_t_vals, L_base, 
                    title = "Net Attrition Curve\n" * 
                            L"""$\tau = %$(τ), \kappa = %$(κ)$""",
                    xlabel = "Years", ylabel = "Survival Rate",
                    legend = nothing,
                    xticks = 0:1:10)
fig_dist = plot(t_vals, d,
                    xlabel = "Years", ylabel = "Normalised Distribution",
                    xticks = 0:n_years)
fig_crop = plot(t_vals, convolution(L_base,d)[1:length(t_vals)],
                    xlabel = "Years", ylabel = "Net Crop "*L"(\Gamma)",
                    title = "Reduced Model Net Crop",
                    xticks = 0:n_years)
plot(fig_attrition, fig_dist, fig_crop, 
        layout = (1,3), size = (1200,300), margins = 5mm)

# # Combine into a single plot
# fig = Figure(size = (1200,300))
# ax1 = Axis(fig[1,1], 
#                 title = #L"Net Attrition Curve\n" * 
#                 L"$\tau = %$(τ), \kappa = %$(κ)$",
#                 xlabel = "Years", ylabel = "Survival Rate",
#                 xticks = 0:1:10)
# ax2 = Axis(fig[1,2], 
#                 title = "Distribution Time Series",
#                 xlabel = "Years", ylabel = "Normalised Distribution",
#                 xticks = 0:n_years)
# ax3 = Axis(fig[1,3], 
#                 title = #"Net Attrition Curve\n" * 
#                 L"""$\tau = %$(τ), \kappa = %$(κ)$""",
#                 xlabel = "Years", ylabel = "Survival Rate",
#                 xticks = 0:n_years)
# lines!(ax1, loss_t_vals, L_base)
# lines!(ax2, t_vals, d)
# lines!(ax3, t_vals, convolution(L_base,d)[1:length(t_vals)])
# fig

# %% Do a parameter scan
# Time series bounds
time_resolution = 1/24
loss_t_vals = 0:time_resolution:5
n_years = 10
t_vals = 0:(time_resolution):(n_years)

# Periodic distribution behaviour
ϕ = 0
k = 3

ω_vals = 1:0.1:2
τ_vals = 0.1:0.05:5
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
            STD_val[ω_i, τ_i, κ_i] = std(response[end-48:end])
        end
    end
end

# %%


heatmap(τ_vals, log.(κ_vals), MAX_val[end,:,:]',
        xlabel = L"\tau", ylabel = L"\log \kappa",
        colorbartitle = L"\Gamma_{SS}",
        labelfontsize = 15, colorbar_titlefontsize = 15,
        title = "Uniform Annual Distribution")

heatmap(τ_vals, log.(κ_vals), STD_val[5,:,:]',
        xlabel = L"\tau", ylabel = L"\log \kappa",
        colorbartitle = L"\Gamma_{SS}",
        labelfontsize = 15, colorbar_titlefontsize = 15,
        title = "Uniform Annual Distribution")

# %% 
max_val_collection = []
std_val_collection = []

for ω_i in 1:length(ω_vals)
    ω = ω_vals[ω_i]
    max_fig = heatmap(τ_vals, log.(κ_vals), MAX_val[ω_i,:,:]',
            xlabel = L"\tau", ylabel = L"\log \kappa",
            colorbartitle = L"\Gamma_{max}",
            labelfontsize = 15, colorbar_titlefontsize = 15,
            title = L"\omega = %$(ω)", clims = (0,5))
            
    std_fig = heatmap(τ_vals, log.(κ_vals), STD_val[ω_i,:,:]',
            xlabel = L"\tau", ylabel = L"\log \kappa",
            colorbartitle = L"\Gamma_{\sigma}",
            labelfontsize = 15, colorbar_titlefontsize = 15,
            title = L"\omega = %$(ω)", clims = (0,0.3))
    push!(max_val_collection, max_fig)
    push!(std_val_collection, std_fig)
end

# %%
plot(max_val_collection..., layout = (3,4), 
        size = (1920,800), margins = 2mm)
plot(std_val_collection..., layout = (3,4), 
        size = (1920,800), margins = 2mm)