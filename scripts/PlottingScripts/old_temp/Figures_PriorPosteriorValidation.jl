"""
Author: Eugene Tan
Date Created: 23/7/2024
Last Updated: 26/8/2024
Code to generate plots of prior vs posterior predictions
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import Public Packages
using JLD2
using CSV
using DataFrames

# %% Import Custom Packages
using NetCropModel

# Maths packages
using LinearAlgebra
using StatsBase
using Distributions
using Plots

# %% Plot theme
pythonplot()
theme(:vibrant)

##########################################
# %% Prior Predictive Trajectories
##########################################

# %% Define countries and import data
ISO = "SEN"
YEAR_START = 2000
YEAR_END = 2023

# %%
input_dict = load("outputs/extractions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
regression_dict = load("outputs/regressions/crop/Compact Regressions/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")

# %%
YEARS_ANNUAL = input_dict["YEARS_ANNUAL"]
MONTHS_MONTHLY = input_dict["MONTHS_MONTHLY"]
DELIVERIES_ANNUAL = input_dict["DELIVERIES_ANNUAL"]
DISTRIBUTION_ANNUAL = input_dict["DISTRIBUTION_ANNUAL"]
NET_CROP_MONTHLY = input_dict["NET_CROP_MONTHLY"]

# %% Generate priors observations
n_sample = 1000

# %% Plot ylims
ylims = (-0.1, 33)

##########################################
# %% Scenario 1: No SGD, all prior
##########################################

# Declare all priors to inference
ϕ_prior = Vector{Real}(undef, n_sample)

β_prior = rand(Uniform(20,24), n_sample) # Hyperprior for ϕ
for i in 1:n_sample
    ϕ_prior[i] = rand(truncated(Beta(2, β_prior[i]),0, 0.25)) # Redistribution parameter
end


α_init_prior = rand(LogNormal(0,1), n_sample) # N_initial nets
α_LLIN_prior = rand(Uniform(0,1), n_sample) # N_initial nets

n_net_types = size(DISTRIBUTION_ANNUAL)[2]-1
b_nets_prior = Array{Real}(undef, (n_sample, n_net_types))
k_nets_prior = Array{Real}(undef, (n_sample, n_net_types))

for i in 1:n_net_types
    # b_nets_prior[:,i] = rand(LogNormal(0.5,0.5), n_sample) # Location parameter Favouring Sigmoid
    # k_nets_prior[:,i] = rand(LogNormal(1.5,0.2), n_sample) # Location parameter Favouring Sigmoid

    b_nets_prior[:,i] = rand(Uniform(0.1,20.7), n_sample) # Location parameter Favouring Sigmoid
    k_nets_prior[:,i] = rand(Uniform(1,30), n_sample) # Location parameter Favouring Sigmoid
end

missing_dist_idxs = findall(ismissing.(DISTRIBUTION_ANNUAL[:,1]))
missing_dist_vals_prior = Array{Real}(undef, n_sample, length(missing_dist_idxs))
for i in 1:size(missing_dist_vals_prior)[2]
    missing_dist_vals_prior[:,i] .= rand(LogNormal(0,1), n_sample)
end

monthly_p = regression_dict["monthly_weights"][1]

# Scenario 1 Prior predictions of net crop
Γ_prior_BYNET = zeros(n_sample, length(MONTHS_MONTHLY), n_net_types)

for i in 1:n_sample
    Γ_sample_prior_BYNET = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                                    ϕ_prior[i], b_nets_prior[i,:], k_nets_prior[i,:], 
                                                    α_init_prior[i], α_LLIN_prior[i],
                                                    missing_dist_vals_prior[i,:];
                                                    monthly_p = monthly_p, return_age = false)
    Γ_prior_BYNET[i,:,:] = Γ_sample_prior_BYNET
end

Γ_prior_TOTAL =  sum(Γ_prior_BYNET, dims = 3)[:,:]

# Calculate quantiles
Γ_prior_TOTAL_quantiles = zeros(3, size(Γ_prior_TOTAL)[2])

for i in 1:size(Γ_prior_TOTAL)[2]
    Γ_prior_TOTAL_quantiles[:,i] = quantile(Γ_prior_TOTAL[:,i], [0.025, 0.5, 0.975])
end

# Make plot for Scenario 1
fig1 = plot(title = "$ISO Prior Predictions\n No SGD, full priors",
            xlabel = "Years", ylabel = "Net Crop (mil)",
            xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90,
            ylims = ylims, xrotation = 90)
plot!(fig1, Γ_prior_TOTAL_quantiles[1,:]./1e6, linealpha = 0,
            fillrange = Γ_prior_TOTAL_quantiles[3,:]./1e6, fillalpha = 0.2,
            color = 2, label = nothing)
plot!(fig1, Γ_prior_TOTAL_quantiles[2,:]./1e6, linealpha = 1, linewidth = 1.6,
            color = 2, label = nothing)

scatter!(NET_CROP_MONTHLY./1e6, color = 4, markersize = 4, label = nothing)

fig1

##########################################
# %% Scenario 2: with SGD, remaining all prior
##########################################

# Declare all priors to inference
ϕ_prior = Vector{Real}(undef, n_sample)

β_prior = rand(Uniform(20,24), n_sample) # Hyperprior for ϕ
for i in 1:n_sample
    ϕ_prior[i] = rand(truncated(Beta(2, β_prior[i]),0, 0.25)) # Redistribution parameter
end


α_init_prior = rand(LogNormal(0,1), n_sample) # N_initial nets
α_LLIN_prior = rand(Uniform(0,1), n_sample) # N_initial nets

n_net_types = size(DISTRIBUTION_ANNUAL)[2]-1
b_nets_prior = Array{Real}(undef, (n_sample, n_net_types))
k_nets_prior = Array{Real}(undef, (n_sample, n_net_types))

for i in 1:n_net_types
    # b_nets_prior[:,i] = rand(LogNormal(0.5,0.5), n_sample) # Location parameter Favouring Sigmoid
    # k_nets_prior[:,i] = rand(LogNormal(1.5,0.2), n_sample) # Location parameter Favouring Sigmoid

    b_nets_prior[:,i] = rand(Uniform(0.1,20.7), n_sample) # Location parameter Favouring Sigmoid
    k_nets_prior[:,i] = rand(Uniform(1,30), n_sample) # Location parameter Favouring Sigmoid
end

missing_dist_idxs = findall(ismissing.(DISTRIBUTION_ANNUAL[:,1]))
missing_dist_vals_prior = Array{Real}(undef, n_sample, length(missing_dist_idxs))
for i in 1:size(missing_dist_vals_prior)[2]
    missing_dist_vals_prior[:,i] .= rand(LogNormal(0,1), n_sample)
end

monthly_p = regression_dict["monthly_weights"][end]

# Scenario 2 Prior predictions of net crop
Γ_prior_BYNET = zeros(n_sample, length(MONTHS_MONTHLY), n_net_types)

for i in 1:n_sample
    Γ_sample_prior_BYNET = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                                    ϕ_prior[i], b_nets_prior[i,:], k_nets_prior[i,:], 
                                                    α_init_prior[i], α_LLIN_prior[i],
                                                    missing_dist_vals_prior[i,:];
                                                    monthly_p = monthly_p, return_age = false)
    Γ_prior_BYNET[i,:,:] = Γ_sample_prior_BYNET
end

Γ_prior_TOTAL =  sum(Γ_prior_BYNET, dims = 3)[:,:]

# Calculate quantiles
Γ_prior_TOTAL_quantiles = zeros(3, size(Γ_prior_TOTAL)[2])

for i in 1:size(Γ_prior_TOTAL)[2]
    Γ_prior_TOTAL_quantiles[:,i] = quantile(Γ_prior_TOTAL[:,i], [0.025, 0.5, 0.975])
end

# Make plot for Scenario 2
fig2 = plot(title = "$ISO Prior Predictions\n SGD, full priors",
            xlabel = "Years", ylabel = "Net Crop (mil)",
            xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90,
            ylims = ylims, xrotation = 90)
plot!(fig2, Γ_prior_TOTAL_quantiles[1,:]./1e6, linealpha = 0,
            fillrange = Γ_prior_TOTAL_quantiles[3,:]./1e6, fillalpha = 0.2,
            color = 2, label = nothing)
plot!(fig2, Γ_prior_TOTAL_quantiles[2,:]./1e6, linealpha = 1, linewidth = 1.6,
            color = 2, label = nothing)

scatter!(fig2, NET_CROP_MONTHLY./1e6, color = 4, markersize = 4, label = nothing)

fig2

##########################################
# %% Scenario 3: with SGD, α post, remaining all prior
##########################################

# Declare all priors to inference
ϕ_prior = Vector{Real}(undef, n_sample)

β_prior = rand(Uniform(20,24), n_sample) # Hyperprior for ϕ
for i in 1:n_sample
    ϕ_prior[i] = rand(truncated(Beta(2, β_prior[i]),0, 0.25)) # Redistribution parameter
end

n_post_draws = length(regression_dict["chain"][:,:α_init])
draw_idxs = sample(1:n_post_draws, n_sample, replace = false)


α_init_post = regression_dict["chain"][draw_idxs,:α_init] # N_initial nets
α_LLIN_post = regression_dict["chain"][draw_idxs,:α_LLIN] # N_initial nets

n_net_types = size(DISTRIBUTION_ANNUAL)[2]-1
b_nets_prior = Array{Real}(undef, (n_sample, n_net_types))
k_nets_prior = Array{Real}(undef, (n_sample, n_net_types))

for i in 1:n_net_types
    # b_nets_prior[:,i] = rand(LogNormal(0.5,0.5), n_sample) # Location parameter Favouring Sigmoid
    # k_nets_prior[:,i] = rand(LogNormal(1.5,0.2), n_sample) # Location parameter Favouring Sigmoid

    b_nets_prior[:,i] = rand(Uniform(0.1,20.7), n_sample) # Location parameter Favouring Sigmoid
    k_nets_prior[:,i] = rand(Uniform(1,30), n_sample) # Location parameter Favouring Sigmoid
end

missing_dist_idxs = findall(ismissing.(DISTRIBUTION_ANNUAL[:,1]))
missing_dist_vals_post = regression_dict["chain"][draw_idxs, end-length(missing_dist_idxs)+1:end]
monthly_p = regression_dict["monthly_weights"][end]

# Scenario 3 Prior predictions of net crop
Γ_prior_BYNET = zeros(n_sample, length(MONTHS_MONTHLY), n_net_types)

for i in 1:n_sample
    Γ_sample_prior_BYNET = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                                    ϕ_prior[i], b_nets_prior[i,:], k_nets_prior[i,:], 
                                                    α_init_post[i], α_LLIN_post[i],
                                                    missing_dist_vals_post[i,:];
                                                    monthly_p = monthly_p, return_age = false)
    Γ_prior_BYNET[i,:,:] = Γ_sample_prior_BYNET
end

Γ_prior_TOTAL =  sum(Γ_prior_BYNET, dims = 3)[:,:]

# Calculate quantiles
Γ_prior_TOTAL_quantiles = zeros(3, size(Γ_prior_TOTAL)[2])

for i in 1:size(Γ_prior_TOTAL)[2]
    Γ_prior_TOTAL_quantiles[:,i] = quantile(Γ_prior_TOTAL[:,i], [0.025, 0.5, 0.975])
end

# Make plot for Scenario 3
fig3 = plot(title = "$ISO Prior Predictions\n SGD, α posterior",
            xlabel = "Years", ylabel = "Net Crop (mil)",
            xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90,
            ylims = ylims, xrotation = 90)
plot!(fig3, Γ_prior_TOTAL_quantiles[1,:]./1e6, linealpha = 0,
            fillrange = Γ_prior_TOTAL_quantiles[3,:]./1e6, fillalpha = 0.2,
            color = 2, label = nothing)
plot!(fig3, Γ_prior_TOTAL_quantiles[2,:]./1e6, linealpha = 1, linewidth = 1.6,
            color = 2, label = nothing)

scatter!(fig3, NET_CROP_MONTHLY./1e6, color = 4, markersize = 4, label = nothing)

fig3

##########################################
# %% Scenario 4: with SGD, α, ϕ post, remaining all prior
##########################################
# Draw idxs from posterior
n_post_draws = length(regression_dict["chain"][:,:α_init])
draw_idxs = sample(1:n_post_draws, n_sample, replace = false)


# Declare all priors to inference
ϕ_post = regression_dict["chain"][draw_idxs,:ϕ]

α_init_post = regression_dict["chain"][draw_idxs,:α_init] # N_initial nets
α_LLIN_post = regression_dict["chain"][draw_idxs,:α_LLIN] # N_initial nets

n_net_types = size(DISTRIBUTION_ANNUAL)[2]-1
b_nets_prior = Array{Real}(undef, (n_sample, n_net_types))
k_nets_prior = Array{Real}(undef, (n_sample, n_net_types))

for i in 1:n_net_types
    # b_nets_prior[:,i] = rand(LogNormal(0.5,0.5), n_sample) # Location parameter Favouring Sigmoid
    # k_nets_prior[:,i] = rand(LogNormal(1.5,0.2), n_sample) # Location parameter Favouring Sigmoid

    b_nets_prior[:,i] = rand(Uniform(0.1,20.7), n_sample) # Location parameter Favouring Sigmoid
    k_nets_prior[:,i] = rand(Uniform(1,30), n_sample) # Location parameter Favouring Sigmoid
end

missing_dist_idxs = findall(ismissing.(DISTRIBUTION_ANNUAL[:,1]))
missing_dist_vals_post = regression_dict["chain"][draw_idxs, end-length(missing_dist_idxs)+1:end]
monthly_p = regression_dict["monthly_weights"][end]

# Scenario 3 Prior predictions of net crop
Γ_prior_BYNET = zeros(n_sample, length(MONTHS_MONTHLY), n_net_types)

for i in 1:n_sample
    Γ_sample_prior_BYNET = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                                    ϕ_post[i], b_nets_prior[i,:], k_nets_prior[i,:], 
                                                    α_init_post[i], α_LLIN_post[i],
                                                    missing_dist_vals_post[i,:];
                                                    monthly_p = monthly_p, return_age = false)
    Γ_prior_BYNET[i,:,:] = Γ_sample_prior_BYNET
end

Γ_prior_TOTAL =  sum(Γ_prior_BYNET, dims = 3)[:,:]

# Calculate quantiles
Γ_prior_TOTAL_quantiles = zeros(3, size(Γ_prior_TOTAL)[2])

for i in 1:size(Γ_prior_TOTAL)[2]
    Γ_prior_TOTAL_quantiles[:,i] = quantile(Γ_prior_TOTAL[:,i], [0.025, 0.5, 0.975])
end

# Make plot for Scenario 3
fig4 = plot(title = "$ISO Prior Predictions\n SGD, α, ϕ posterior",
            xlabel = "Years", ylabel = "Net Crop (mil)",
            xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90,
            ylims = ylims, xrotation = 90)
plot!(fig4, Γ_prior_TOTAL_quantiles[1,:]./1e6, linealpha = 0,
            fillrange = Γ_prior_TOTAL_quantiles[3,:]./1e6, fillalpha = 0.2,
            color = 2, label = nothing)
plot!(fig4, Γ_prior_TOTAL_quantiles[2,:]./1e6, linealpha = 1, linewidth = 1.6,
            color = 2, label = nothing)

scatter!(fig4, NET_CROP_MONTHLY./1e6, color = 4, markersize = 4, label = nothing)

fig4

##########################################
# %% Scenario 4: with SGD, full posterior
##########################################
# Draw idxs from posterior
chain = regression_dict["chain"]
n_post_draws = length(regression_dict["chain"][:,:α_init])
draw_idxs = sample(1:n_post_draws, n_sample, replace = false)


# Declare all priors to inference
ϕ_post = regression_dict["chain"][draw_idxs,:ϕ]

α_init_post = regression_dict["chain"][draw_idxs,:α_init] # N_initial nets
α_LLIN_post = regression_dict["chain"][draw_idxs,:α_LLIN] # N_initial nets

n_net_types = size(DISTRIBUTION_ANNUAL)[2]-1
b_nets_post = Matrix(DataFrame(chain)[:,5:2:5+2*(n_net_types-1)])[draw_idxs, :]
k_nets_post = Matrix(DataFrame(chain)[:,6:2:6+2*(n_net_types-1)])[draw_idxs, :]

missing_dist_idxs = findall(ismissing.(DISTRIBUTION_ANNUAL[:,1]))
missing_dist_vals_post = regression_dict["chain"][draw_idxs, end-length(missing_dist_idxs)+1:end]
monthly_p = regression_dict["monthly_weights"][end]

# Scenario 3 Prior predictions of net crop
Γ_prior_BYNET = zeros(n_sample, length(MONTHS_MONTHLY), n_net_types)

for i in 1:n_sample
    Γ_sample_prior_BYNET = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                                    ϕ_prior[i], b_nets_post[i,:], k_nets_post[i,:], 
                                                    α_init_post[i], α_LLIN_post[i],
                                                    missing_dist_vals_post[i,:];
                                                    monthly_p = monthly_p, return_age = false)
    Γ_prior_BYNET[i,:,:] = Γ_sample_prior_BYNET
end

Γ_prior_TOTAL =  sum(Γ_prior_BYNET, dims = 3)[:,:]

# Calculate quantiles
Γ_prior_TOTAL_quantiles = zeros(3, size(Γ_prior_TOTAL)[2])

for i in 1:size(Γ_prior_TOTAL)[2]
    Γ_prior_TOTAL_quantiles[:,i] = quantile(Γ_prior_TOTAL[:,i], [0.025, 0.5, 0.975])
end

# Make plot for Scenario 3
fig5 = plot(title = "$ISO Prior Predictions\n SGD, full posterior",
            xlabel = "Years", ylabel = "Net Crop (mil)",
            xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90,
            ylims = ylims, xrotation = 90)
plot!(fig5, Γ_prior_TOTAL_quantiles[1,:]./1e6, linealpha = 0,
            fillrange = Γ_prior_TOTAL_quantiles[3,:]./1e6, fillalpha = 0.2,
            color = 2, label = nothing)
plot!(fig5, Γ_prior_TOTAL_quantiles[2,:]./1e6, linealpha = 1, linewidth = 1.6,
            color = 2, label = nothing)

scatter!(fig5, NET_CROP_MONTHLY./1e6, color = 4, markersize = 4, label = nothing)

fig5

# %% Combine all scenario plots together
fig_full = plot(fig1, fig2, fig3, fig4, fig5, layout = (1,5), size = (2400,400))
savefig(fig_full, "$(ISO)_Prior_Posterior_predictions.pdf")