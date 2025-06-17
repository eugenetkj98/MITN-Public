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
include(pwd()*"/scripts/read_toml.jl")

# %% Math packages
using LinearAlgebra
using StatsBase

# %% Plotting Packages
using CairoMakie

# %% General Packages
using ProgressBars
using LaTeXStrings
using CSV
using DataFrames
using JLD2

# %% Custom Packages
using Convolutions
using NetLoss

# %% Basic Plot settings
set_theme!(theme_latexfonts())

# %% Load collection of attrition parameters from current regressions
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %%
τ_regs = Vector{Float64}(undef, length(filt_ISOs))
κ_regs = Vector{Float64}(undef, length(filt_ISOs))
for ISO_i in ProgressBar(1:length(filt_ISOs))
    ISO = filt_ISOs[ISO_i]
    data = JLD2.load(OUTPUT_REGRESSIONS_DIR*"crop/$(YEAR_NAT_START)_$(YEAR_NAT_END)/$(ISO)_$(YEAR_NAT_START)_$(YEAR_NAT_END)_cropchains.jld2")
    τ_regs[ISO_i] = mean(vcat(data["chain_epochs"]...)[:,7])
    κ_regs[ISO_i] = mean(vcat(data["chain_epochs"]...)[:,8])
end

######################################
# %% Scenario 1: Uniform Annual Distribution
######################################
# %% Define loss/net attrition function
loss_t_vals = 0:1/12:5

## Compact Exponential
τ = 6
κ = 10
L_base_compact = net_loss_compact.(loss_t_vals, τ, κ)

# Weibull
b = 1.0
k = 3.0
L_base_weibull = net_loss_weibull.(loss_t_vals, b, k)

# Nets distribution behaviour
n_years = 10
t_vals = 0:(1/12):(n_years)

# Once a year distribution
d = zeros(length(t_vals)).= 1/12

# Combine into a single plot
fig = Figure(size = (800,450))
ax1 = Axis(fig[1,1], 
                title = L"$\text{\tau} = %$(τ), \, \text{\kappa} = %$(κ)$",
                titlesize = 20,
                xlabel = "Years", xlabelsize = 18,
                ylabel = "Survival Rate", ylabelsize = 18,
                xticks = 0:1:10)
ax2 = Axis(fig[1,2], 
                title = L"$\text{\tau} = %$(τ), \, \text{\kappa} = %$(κ)$",
                titlesize = 20,
                xlabel = "Years", xlabelsize = 18,
                ylabel = "Net Crop (Γ)",ylabelsize = 18,
                xticks = 0:1:10)
ax3 = Axis(fig[2,1], 
                title = L"$b = %$(b), \, k = %$(k)$",
                titlesize = 20,
                xlabel = "Years", xlabelsize = 18,
                ylabel = "Survival Rate", ylabelsize = 18,
                xticks = 0:1:10)
ax4 = Axis(fig[2,2], 
                title = L"$b = %$(b), \, k = %$(k)$",
                titlesize = 20,
                xlabel = "Years", xlabelsize = 18,
                ylabel = "Net Crop (Γ)",ylabelsize = 18,
                xticks = 0:1:10)
ylims!(ax1, -0.05, 1.05)
ylims!(ax3, -0.05, 1.05)
lines!(ax1, loss_t_vals, L_base_compact)
lines!(ax2, t_vals, convolution(L_base_compact,d)[1:length(t_vals)])

lines!(ax3, loss_t_vals, L_base_weibull)
lines!(ax4, t_vals, convolution(L_base_weibull,d)[1:length(t_vals)])
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

b_vals = 0.5:0.1:5
k_vals = exp.(-2:0.1:4)

SS_val_compact = zeros(length(τ_vals), length(κ_vals))
SS_val_weibull = zeros(length(b_vals), length(k_vals))

for i in 1:length(τ_vals)
    for j in 1:length(κ_vals)
        τ = τ_vals[i]
        κ = κ_vals[j]
        L = net_loss_compact.(loss_t_vals, τ, κ)
        SS_val_compact[i,j] = maximum(convolution(L,d)[1:length(t_vals)])
    end
end

for i in 1:length(b_vals)
    for j in 1:length(k_vals)
        b = b_vals[i]
        k = k_vals[j]
        L = net_loss_weibull.(loss_t_vals, b, k)
        SS_val_weibull[i,j] = maximum(convolution(L,d)[1:length(t_vals)])
    end
end

fig = Figure(size = (800,400))
ax1 = Axis(fig[1,1],
        title = L"\text{Compact Exponential}",
        titlesize = 25,
        xlabel = L"\tau", xlabelsize = 20,
        ylabel = L"\kappa", ylabelsize = 20,
        yscale = log,
        yticks = round.(exp.(-4:1:4), digits = 2)
        )
ax2 = Axis(fig[1,2],
        title = L"\text{Weibull}",
        titlesize = 25,
        xlabel = L"b", xlabelsize = 20,
        ylabel = L"k", ylabelsize = 20,
        yscale = log,
        yticks = round.(exp.(-1:1:4), digits = 2)
        )
clims = (0,4)
Colorbar(fig[1,3], 
        colorrange = clims, colormap = :viridis, 
        label = L"\Gamma_{SS}", labelsize = 20)
heatmap!(ax1, τ_vals, κ_vals, SS_val_compact, colorrange = clims)
heatmap!(ax2, b_vals, k_vals, SS_val_weibull, colorrange = clims)


contour!(ax1, τ_vals, κ_vals, [net_life_compact(0.5, τ, κ) for τ in τ_vals, κ in κ_vals],
            levels = 1:0.5:3, labels = true, color = :black, labelsize = 20)
contour!(ax2, b_vals, k_vals, [net_life_weibull(0.5, b, k) for b in b_vals, k in k_vals],
            levels = 1:0.5:3, labels = true, color = :black, labelsize = 20)

save(OUTPUT_PLOTS_DIR*"PaperFigures/TechnicalPaper/Uniform_Heatmaps.pdf", fig, pdf_version = "1.4")
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
τ_vals = [2,4,10,20]
κ_vals = [0.5,1.5,6,12]

# b_vals = [1.0,2.0,3.0,4.0]
# k_vals = [1.0,4.0,8.0,20.0]

b_vals = [1.5, 2.0, 2.5, 3.5]
k_vals = [0.8,1.5,4,10]

L_base_vals_compact = []
conv_vals_compact = []

for i in 1:length(τ_vals)
    push!(L_base_vals_compact, net_loss_compact.(loss_t_vals, τ_vals[i], κ_vals[i]))
end

L_base_vals_weibull = []
conv_vals_weibull = []

for i in 1:length(b_vals)
    push!(L_base_vals_weibull, net_loss_weibull.(loss_t_vals, b_vals[i], k_vals[i]))
end

L_base_cat_compact = hcat(L_base_vals_compact...)
L_base_cat_weibull = hcat(L_base_vals_weibull...)



# Combine into a single plot
fig = Figure(size = (1200,500))
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
                ylabel = "Normalised\nDistribution", ylabelsize = 18,
                xticks = 0:n_years)
ax3 = Axis(fig[1,3], 
                title = L"\text{Effective Net Crop}",
                titlesize = 22,
                xlabel = "Years",  xlabelsize = 18, 
                ylabel = L"\text{Net Crop }(\Gamma)", ylabelsize = 18,
                xticks = 0:n_years)
ax4 = Axis(fig[2,1], 
                title = L"\text{Net Attrition Curve}",
                titlesize = 22,
                xlabel = "Years",  xlabelsize = 18,
                ylabel = "Survival Rate", ylabelsize = 18,
                xticks = 0:1:10)
ax5 = Axis(fig[2,2], 
                title = L"\text{Distribution Time Series}",
                titlesize = 22,
                xlabel = "Years",  xlabelsize = 18,
                ylabel = "Normalised\nDistribution", ylabelsize = 18,
                xticks = 0:n_years)
ax6 = Axis(fig[2,3], 
                title = L"\text{Effective Net Crop}",
                titlesize = 22,
                xlabel = "Years",  xlabelsize = 18, 
                ylabel = L"\text{Net Crop }(\Gamma)", ylabelsize = 18,
                xticks = 0:n_years)

survival_lines_compact = []
conv_lines_compact = []
for i in 1:length(τ_vals)
    push!(survival_lines_compact, lines!(ax1, loss_t_vals, L_base_cat_compact[:,i]))
    push!(conv_lines_compact, lines!(ax3, t_vals, convolution(L_base_cat_compact[:,i],d)[1:length(t_vals)]))
end

survival_lines_weibull = []
conv_lines_weibull = []
for i in 1:length(b_vals)
    push!(survival_lines_weibull, lines!(ax4, loss_t_vals, L_base_cat_weibull[:,i]))
    push!(conv_lines_weibull, lines!(ax6, t_vals, convolution(L_base_cat_weibull[:,i],d)[1:length(t_vals)]))
end
Legend(fig[1,4], survival_lines_compact,
        # tellwidth = false, tellheight = false,
        labelsize = 18,
        halign = :right, valign = :center,
        padding = (5,5,5,5),
        margin = (5,5,5,5),
        [L"$\tau = %$(τ_vals[i]), \,\kappa = %$(round(κ_vals[i], digits = 2))$" for i in 1:length(τ_vals)])
Legend(fig[2,4], survival_lines_weibull,
        # tellwidth = false, tellheight = false,
        labelsize = 18,
        halign = :right, valign = :center,
        padding = (5,5,5,5),
        margin = (5,5,5,5),
        [L"$b = %$(b_vals[i]), \,k = %$(round(k_vals[i], digits = 2))$" for i in 1:length(b_vals)])
# Legend(fig[1,3], conv_lines,
#         tellwidth = false, tellheight = false,
#         labelsize = 18,
#         halign = :right, valign = :bottom,
#         padding = (5,5,5,5),
#         margin = (5,5,5,5),
#         [L"$\tau = %$(τ_vals[i]), \,\kappa = %$(round(κ_vals[i], digits = 2))$" for i in 1:length(τ_vals)])

lines!(ax2, t_vals, d)
lines!(ax5, t_vals, d)

xlims!(ax1, -0.5, 5.5)
xlims!(ax2, -0.5, 5.5)
# xlims!(ax3, -0.5, 5.5)
xlims!(ax4, -0.5, 5.5)
xlims!(ax5, -0.5, 5.5)
# xlims!(ax6, -0.5, 5.5)
# save(OUTPUT_PLOTS_DIR*"PaperFigures/TechnicalPaper/Periodic_Dist_Timeseries.pdf", fig, pdf_version = "1.4")
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

τ_vals = [1,2,5,20]
κ_vals = [10,10,10,10]



L_base_vals_compact = []
conv_vals_compact = []


for i in 1:length(τ_vals)
    push!(L_base_vals_compact, net_loss_compact.(loss_t_vals, τ_vals[i], κ_vals[i]))
end

L_base_cat_compact = hcat(L_base_vals_compact...)



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
    push!(survival_lines, lines!(ax1, loss_t_vals, L_base_cat_compact[:,i]))
    push!(conv_lines, lines!(ax3, t_vals, convolution(L_base_cat_compact[:,i],d)[1:length(t_vals)]))
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
# %% Figure 3a: Parameter scan for periodic distribution (COMPACT)
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
ω_vals = [0.5,1,1.5,2]
τ_vals = 1:0.25:20
κ_vals_linear = 5:0.25:22
κ_vals_exp = exp.(-4:0.1:4)

MEAN_val_linear = zeros(length(ω_vals), length(τ_vals), length(κ_vals_linear))
STD_val_linear = zeros(length(ω_vals), length(τ_vals), length(κ_vals_linear))

MEAN_val_exp = zeros(length(ω_vals), length(τ_vals), length(κ_vals_exp))
STD_val_exp = zeros(length(ω_vals), length(τ_vals), length(κ_vals_exp))

for ω_i in ProgressBar(1:length(ω_vals), leave = false)
    Threads.@threads for τ_i in ProgressBar(1:length(τ_vals), leave = false)
        for κ_i in 1:length(κ_vals_linear)
            ω = ω_vals[ω_i]
            τ = τ_vals[τ_i]
            κ = κ_vals_linear[κ_i]

            # Construct distribution time series
            d = zeros(length(t_vals))
            d_shape = abs.(sin.(π .* (ω .* t_vals .- ϕ))).^k
            d = ((n_years)/sum(d_shape)) .* d_shape

            # Construct response
            L = net_loss_compact.(loss_t_vals, τ, κ)
            response = convolution(L,d)[1:length(t_vals)]

            # Store Variables
            MEAN_val_linear[ω_i, τ_i, κ_i] = mean(response)
            STD_val_linear[ω_i, τ_i, κ_i] = std(response[end-20*subsample_ratio:end])
        end

        for κ_i in 1:length(κ_vals_exp)
            ω = ω_vals[ω_i]
            τ = τ_vals[τ_i]
            κ = κ_vals_exp[κ_i]

            # Construct distribution time series
            d = zeros(length(t_vals))
            d_shape = abs.(sin.(π .* (ω .* t_vals .- ϕ))).^k
            d = ((n_years)/sum(d_shape)) .* d_shape

            # Construct response
            L = net_loss_compact.(loss_t_vals, τ, κ)
            response = convolution(L,d)[1:length(t_vals)]

            # Store Variables
            MEAN_val_exp[ω_i, τ_i, κ_i] = mean(response)
            STD_val_exp[ω_i, τ_i, κ_i] = std(response[end-20*subsample_ratio:end])
        end
    end
end

# %% Heatmap for steady state mean
fig = Figure(size = (900,450))

# axs1 = fig[1,1] = GridLayout()
# axs2 = fig[2,1] = GridLayout()

clims = (0,5)
for ω_i in 1:length(ω_vals)
    # Linear Case
    ax = Axis(fig[1, ω_i],xlabel = L"\tau", xlabelsize = 22,
            ylabel = L"\kappa", ylabelsize = 22,
            title = L"\omega = %$(ω_vals[ω_i])", titlesize = 25,
            yticks = 5:5:22)
    heatmap!(ax, τ_vals, κ_vals_linear, MEAN_val_linear[ω_i, :,:], colorrange = clims)
    scatter!(ax, τ_regs, κ_regs, markersize = 4, color = :red)
    contour!(ax, τ_vals, κ_vals_linear, [net_life_compact(0.5, τ, κ) for τ in τ_vals, κ in κ_vals_linear],
                levels = 1:1:5, labels = true, color = (:black, 0.7), labelsize = 10)
    xlims!(ax, 1,20)
    ylims!(ax, 5,22)

    # Exponential Case
    ax = Axis(fig[2, ω_i],xlabel = L"\tau", xlabelsize = 22,
            ylabel = L"\kappa", ylabelsize = 22,
            title = L"\omega = %$(ω_vals[ω_i])", titlesize = 25,
            yscale = log, yticks = round.(exp.(-4:1:4), digits = 2))
    heatmap!(ax, τ_vals, κ_vals_exp, MEAN_val_exp[ω_i, :,:], colorrange = clims)
    scatter!(ax, τ_regs, κ_regs, markersize = 4, color = :red)
    contour!(ax, τ_vals, κ_vals_exp, [net_life_compact(0.5, τ, κ) for τ in τ_vals, κ in κ_vals_exp],
                levels = 1:1:5, labels = true, color = (:black, 0.7), labelsize = 10)
    xlims!(ax, 1,20)
    ylims!(ax, minimum(κ_vals_exp), maximum(κ_vals_exp)+0.5)
end

Colorbar(fig[1:2,length(ω_vals)+1], colorrange = clims, colormap = :viridis, 
        label = L"\bar{\Gamma}", labelsize = 20)
fig

# %%
save(OUTPUT_PLOTS_DIR*"PaperFigures/TechnicalPaper/Periodic_MEAN_Heatmap_COMPACT.pdf", fig)

# %% Net Crop Std
fig = Figure(size = (900,450))

# axs1 = fig[1,1] = GridLayout()
# axs2 = fig[2,1] = GridLayout()

clims = (0,0.5)
for ω_i in 1:length(ω_vals)
    # Linear Case
    ax = Axis(fig[1, ω_i],xlabel = L"\tau", xlabelsize = 22,
            ylabel = L"\kappa", ylabelsize = 22,
            title = L"\omega = %$(ω_vals[ω_i])", titlesize = 25,
            yticks = 5:5:22)
    heatmap!(ax, τ_vals, κ_vals_linear, STD_val_linear[ω_i, :,:], colorrange = clims)
    scatter!(ax, τ_regs, κ_regs, markersize = 4, color = :red)
    contour!(ax, τ_vals, κ_vals_linear, [net_life_compact(0.5, τ, κ) for τ in τ_vals, κ in κ_vals_linear],
                levels = 1:1:5, labels = true, color = (:black, 0.7), labelsize = 10)
    xlims!(ax, 1,20)
    ylims!(ax, 5,22)

    # Exponential Case
    ax = Axis(fig[2, ω_i],xlabel = L"\tau", xlabelsize = 22,
            ylabel = L"\kappa", ylabelsize = 22,
            title = L"\omega = %$(ω_vals[ω_i])", titlesize = 25,
            yscale = log, yticks = round.(exp.(-4:1:4), digits = 2))
    heatmap!(ax, τ_vals, κ_vals_exp, STD_val_exp[ω_i, :,:], colorrange = clims)
    scatter!(ax, τ_regs, κ_regs, markersize = 4, color = :red)
    contour!(ax, τ_vals, κ_vals_exp, [net_life_compact(0.5, τ, κ) for τ in τ_vals, κ in κ_vals_exp],
                levels = 1:1:5, labels = true, color = (:black, 0.7), labelsize = 10)
    xlims!(ax, 1,20)
    ylims!(ax, minimum(κ_vals_exp), maximum(κ_vals_exp)+0.5)
end

Colorbar(fig[1:2,length(ω_vals)+1], colorrange = clims, colormap = :viridis, 
        label = L"\sigma_{\Gamma}", labelsize = 20)
fig

# %%
save(OUTPUT_PLOTS_DIR*"PaperFigures/TechnicalPaper/Periodic_STD_Heatmap_COMPACT.pdf", fig)

######################################
# %% Figure 3b: Parameter scan for periodic distribution (WEIBULL)
######################################
# %% Do a parameter scan
# Time series bounds
subsample_ratio = 36
time_resolution = 1/subsample_ratio
loss_t_vals = 0:time_resolution:5
n_years = 40
t_vals = 0:(time_resolution):(n_years)
# Periodic distribution behaviour
ω_vals = [0.5,1,1.5,2]
ϕ = 0
K = 3
b_vals = 0.5:0.1:10
k_vals = exp.(-1:0.1:4)

MEAN_val = zeros(length(ω_vals), length(b_vals), length(k_vals))
STD_val = zeros(length(ω_vals), length(b_vals), length(k_vals))

for ω_i in ProgressBar(1:length(ω_vals), leave = false)
    Threads.@threads for b_i in ProgressBar(1:length(b_vals), leave = false)
        for k_i in 1:length(k_vals)
            ω = ω_vals[ω_i]
            b = b_vals[b_i]
            k = k_vals[k_i]

            # Construct distribution time series
            d = zeros(length(t_vals))
            d_shape = abs.(sin.(π .* (ω .* t_vals .- ϕ))).^K
            d = ((n_years)/sum(d_shape)) .* d_shape

            # Construct response
            L = net_loss_weibull.(loss_t_vals, b,k)
            response = convolution(L,d)[1:length(t_vals)]

            # Store Variables
            MEAN_val[ω_i, b_i, k_i] = mean(response)
            STD_val[ω_i, b_i, k_i] = std(response[end-20*subsample_ratio:end])
        end
    end
end

# %% Heatmap for steady state mean
fig = Figure(size = (900,450))

mean_clims = (0,5)
std_clims = (0,0.5)
for ω_i in 1:length(ω_vals)
    ax = Axis(fig[1, ω_i],
                xlabel = L"b", xlabelsize = 22,
                ylabel = L"k", ylabelsize = 22,
                title = L"\omega = %$(ω_vals[ω_i])", titlesize = 25,
                yscale = log, yticks = round.(exp.(-1:1:4), digits = 2))
    heatmap!(ax, b_vals, k_vals, MEAN_val[ω_i, :,:], colorrange = mean_clims)
    contour!(ax, b_vals, k_vals, [net_life_weibull(0.5, b, k) for b in b_vals, k in k_vals],
        levels = 1:1:5, labels = true, color = (:black, 0.7), labelsize = 10)

    ax = Axis(fig[2, ω_i],
            xlabel = L"b", xlabelsize = 22,
            ylabel = L"k", ylabelsize = 22,
            title = L"\omega = %$(ω_vals[ω_i])", titlesize = 25,
            yscale = log, yticks = round.(exp.(-1:1:4), digits = 2))
    heatmap!(ax, b_vals, k_vals, STD_val[ω_i, :,:], colorrange = std_clims)
    contour!(ax, b_vals, k_vals, [net_life_weibull(0.5, b, k) for b in b_vals, k in k_vals],
    levels = 1:1:5, labels = true, color = (:black, 0.7), labelsize = 10)
end

Colorbar(fig[1,length(ω_vals)+1], colorrange = mean_clims, colormap = :viridis, 
        label = L"\bar{\Gamma}", labelsize = 20)
Colorbar(fig[2,length(ω_vals)+1], colorrange = std_clims, colormap = :viridis, 
        label = L"\sigma_{\Gamma}", labelsize = 20)

fig


# %%
save(OUTPUT_PLOTS_DIR*"PaperFigures/TechnicalPaper/Periodic_MEANSTD_Heatmap_WEIBULL.pdf", fig)
