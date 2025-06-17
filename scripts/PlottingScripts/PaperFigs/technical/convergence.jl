"""
Author: Eugene Tan
Date Created: 6/5/2025
Last Updated: 6/5/2025
Make plots for convergence behaviour of EM algorithm (Nigeria)
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/read_toml.jl")

# %% Import packages
using JLD2
using CairoMakie
using LaTeXStrings
using DataFrames
using StatsBase
using NetCropRegression


# %%
input_dir = OUTPUT_REGRESSIONS_DIR*"crop/convergence_test/"
n_samples = 20
ISO = "NGA"

# %% Define Storage variables
p_monthly_samples = Vector{Matrix{Float64}}(undef, n_samples)
loss_samples = Vector{Vector{Float64}}(undef, n_samples)
τ_ci_samples = Vector{Matrix{Float64}}(undef, n_samples)
κ_ci_samples = Vector{Matrix{Float64}}(undef, n_samples)

# %%
for sample_idx in 1:n_samples
    # Load data
    data = JLD2.load(input_dir*"$(ISO)_convergence_sample_$(sample_idx).jld2")

    # Extract MCMC draws and CI
    n_epochs = length(data["chain_epochs"])

    τ_ci = Matrix{Float64}(undef, n_epochs, 3)
    κ_ci = Matrix{Float64}(undef, n_epochs, 3)

    for i in 1:n_epochs
        τ_ci[i,[1,3]] = quantile(data["chain_epochs"][i][:,7], [0.025,0.975])
        τ_ci[i,2] = mean(data["chain_epochs"][i][:,7])
        κ_ci[i,[1,3]] = quantile(data["chain_epochs"][i][:,8], [0.025,0.975])
        κ_ci[i,2] = mean(data["chain_epochs"][i][:,8])
    end
    τ_ci_samples[sample_idx] = τ_ci
    κ_ci_samples[sample_idx] = κ_ci

    # Extract convergence behaviour
    loss_samples[sample_idx] = data["loss"]
    p_monthly_samples[sample_idx] = hcat(data["monthly_weights"]...)
end

##############################################################
# %% Fig1: P monthly Convergence Plots
##############################################################
# Plot settings
set_theme!(theme_latexfonts())
ls = 20
lalpha = 0.2
lw = 1.5


# Calculate end state monthly p_mean
p_monthly_mean = normalised_monthly_weights(mean([p_monthly_samples[i][:,end] for i in 1:n_samples]))

# Calculate normalised RMSE
ϵ_values = Vector{Vector{Float64}}(undef, n_samples)

# Helper function to calculate RMSE
RMS(x) = sqrt(mean(x.^2))

# Extract p_monthly RMSE values
for sample_idx in 1:n_samples
    ϵ_values[sample_idx] = [RMS(p_monthly_samples[sample_idx][:,i] .- p_monthly_mean)/RMS(p_monthly_mean) for i in 1:size(p_monthly_samples[sample_idx])[2]]
end

# Get month idx where surveys were available
input_dict = load(OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_NAT_START)_$(YEAR_NAT_END)/$(ISO)_$(YEAR_NAT_START)_$(YEAR_NAT_END)_cropextract.jld2")

survey_monthidxs = findall(.!ismissing.(input_dict["NET_CROP_MONTHLY"]))

# %% Make Plots
fig = Figure(size = (800,550))
ax1 = Axis(fig[1,1],
            xlabel = "SGD Steps",
            ylabel = L"L",
            xlabelsize = ls,
            ylabelsize = ls*1.2)
ax2 = Axis(fig[1,2],
            xlabel = "SGD Steps",
            ylabel = L"\epsilon",
            xlabelsize = ls,
            ylabelsize = ls*1.2)
ax3 = Axis(fig[2,1:2],
            xlabel = "Month Index (t)",
            ylabel = L"p_t",
            xlabelsize = ls,
            ylabelsize = ls*1.2)

elems = []
for sample_idx in 1:n_samples
    lines!(ax1, loss_samples[sample_idx], 
            color = (:blue, lalpha), linewidth = lw)
    lines!(ax2, ϵ_values[sample_idx], 
            color = (:blue, lalpha), linewidth = lw)
    if sample_idx == 1
        push!(elems,lines!(ax3, p_monthly_samples[sample_idx][:,end],
                color = (:blue, lalpha), linewidth = lw))
    else
        lines!(ax3, p_monthly_samples[sample_idx][:,end],
                color = (:blue, lalpha), linewidth = lw)
    end
end

vlines!(ax1,2:30:(maximum([length(loss_samples[i]) for i in 1:n_samples])),
        color = :red, linestyle = :dash)

push!(elems,lines!(ax3, p_monthly_mean, color = :red))

push!(elems,vlines!(ax3, survey_monthidx, 
        color = (:green, 0.4), linewidth = 2))

Legend(fig[2,1:2], elems,
            ["Samples","Mean","Survey Times"],
            rowgap = -5,
            halign = :right, valign = :top, 
            margin = (5,5,5,5), labelsize = 14)

save(OUTPUT_PLOTS_DIR*"PaperFigures/TechnicalPaper/Convergence.pdf", fig, pdf_version = "1.4")
fig