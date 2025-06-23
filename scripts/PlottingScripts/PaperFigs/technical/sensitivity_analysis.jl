"""
Author: Eugene Tan
Date Created: 6/5/2025
Last Updated: 6/5/2025
Make plots for effect of data on model convergence (Nigeria)
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
using NetLoss


# %% Define filenames and directories
input_dir = OUTPUT_REGRESSIONS_DIR*"crop/convergence_test/sensitivity_analysis/"
ISO = "NGA"
n_data_vals = min.(2:2:40,39) # Need to manually get from dataset
rand_sample_size = 10
chronological_filenames = "NGA_chronological_".*string.(n_data_vals).*".jld2"
rev_chronological_filenames = "NGA_rev_chronological_".*string.(n_data_vals).*".jld2"

# %% Storage Variables for posterior estimates
# chronological_τ_ests = Array{Float64}(undef, length(n_data_vals), 2, 3)
# chronological_κ_ests = Array{Float64}(undef, length(n_data_vals), 2, 3)
chronological_halflife_ests = Array{Float64}(undef, length(n_data_vals), 2, 3)

# rev_chronological_τ_ests = Array{Float64}(undef, length(n_data_vals), 2, 3)
# rev_chronological_κ_ests = Array{Float64}(undef, length(n_data_vals), 2, 3)
rev_chronological_halflife_ests = Array{Float64}(undef, length(n_data_vals), 2, 3)

rand_halflife_ests = Array{Float64}(undef, length(n_data_vals), 2, 3)

# %% Extract Posterior estimates
for i in 1:length(n_data_vals)
    # Chronological Case
    data = JLD2.load(input_dir*"chronological/"*chronological_filenames[i])
        
    τ_vals_citn = data["chain"][:,5]
    κ_vals_citn = data["chain"][:,6]
    τ_vals_llin = data["chain"][:,7]
    κ_vals_llin = data["chain"][:,8]

    # chronological_τ_ests[i,[1,3]] .= quantile(τ_vals_llin, [0.025, 0.975])
    # chronological_τ_ests[i,2] = mean(data["chain"][:,7])
    # chronological_κ_ests[i,[1,3]] .= quantile(κ_vals_llin, [0.025, 0.975])
    # chronological_κ_ests[i,2] = mean(data["chain"][:,8])

    halflife_vals_citn = net_life_compact.(0.5, τ_vals_citn, κ_vals_citn)
    halflife_vals_llin = net_life_compact.(0.5, τ_vals_llin, κ_vals_llin)
    
    chronological_halflife_ests[i,1,[1,3]] .= quantile(halflife_vals_citn, [0.05, 0.95])
    chronological_halflife_ests[i,1,2] = mean(halflife_vals_citn)
    chronological_halflife_ests[i,2,[1,3]] .= quantile(halflife_vals_llin, [0.05, 0.95])
    chronological_halflife_ests[i,2,2] = mean(halflife_vals_llin)

    # Reverse Chronological Case
    data = JLD2.load(input_dir*"rev_chronological/"*rev_chronological_filenames[i])
        
    τ_vals_citn = data["chain"][:,5]
    κ_vals_citn = data["chain"][:,6]
    τ_vals_llin = data["chain"][:,7]
    κ_vals_llin = data["chain"][:,8]

    # chronological_τ_ests[i,[1,3]] .= quantile(τ_vals_llin, [0.025, 0.975])
    # chronological_τ_ests[i,2] = mean(data["chain"][:,7])
    # chronological_κ_ests[i,[1,3]] .= quantile(κ_vals_llin, [0.025, 0.975])
    # chronological_κ_ests[i,2] = mean(data["chain"][:,8])

    halflife_vals_citn = net_life_compact.(0.5, τ_vals_citn, κ_vals_citn)
    halflife_vals_llin = net_life_compact.(0.5, τ_vals_llin, κ_vals_llin)
    
    rev_chronological_halflife_ests[i,1,[1,3]] .= quantile(halflife_vals_citn, [0.05, 0.95])
    rev_chronological_halflife_ests[i,1,2] = mean(halflife_vals_citn)
    rev_chronological_halflife_ests[i,2,[1,3]] .= quantile(halflife_vals_llin, [0.05, 0.95])
    rev_chronological_halflife_ests[i,2,2] = mean(halflife_vals_llin)

    # Random Case
    rand_halflife_vals_citn = []
    rand_halflife_vals_llin = []

    for rand_i in 1:rand_sample_size
        data = JLD2.load(input_dir*"random/"*"$(ISO)_random_$(n_data_vals[i])_sample_$(rand_i).jld2")
            
        τ_vals_citn = data["chain"][:,5]
        κ_vals_citn = data["chain"][:,6]
        τ_vals_llin = data["chain"][:,7]
        κ_vals_llin = data["chain"][:,8]

        halflife_vals_citn = net_life_compact.(0.5, τ_vals_citn, κ_vals_citn)
        halflife_vals_llin = net_life_compact.(0.5, τ_vals_llin, κ_vals_llin)
        
        rand_halflife_vals_citn = vcat(rand_halflife_vals_citn, halflife_vals_citn)
        rand_halflife_vals_llin = vcat(rand_halflife_vals_llin, halflife_vals_llin)
    end


    rand_halflife_ests[i,1,[1,3]] .= quantile(rand_halflife_vals_citn, [0.05, 0.95])
    rand_halflife_ests[i,1,2] = mean(rand_halflife_vals_citn)
    rand_halflife_ests[i,2,[1,3]] .= quantile(rand_halflife_vals_llin, [0.05, 0.95])
    rand_halflife_ests[i,2,2] = mean(rand_halflife_vals_llin)
end

#######################################################
# %% Make Plots
#######################################################
set_theme!(theme_latexfonts())
ls = 20
alpha = 0.4
lw = 2

fig = Figure(size = (1000,350))
ax1 = Axis(fig[1,1],
            xlabel = L"n",
            ylabel = "Posterior Halflife (Years)",
            xlabelsize = ls,
            ylabelsize = ls,
            title = "Chronological",
            titlesize = ls*1.3)
ax2 = Axis(fig[1,2],
            xlabel = L"n",
            ylabel = "Posterior Halflife (Years)",
            xlabelsize = ls,
            ylabelsize = ls,
            title = "Reverse Chrono.",
            titlesize = ls*1.3)
ax3 = Axis(fig[1,3],
            xlabel = L"n",
            ylabel = "Posterior Halflife (Years)",
            xlabelsize = ls,
            ylabelsize = ls,
            title = "Random",
            titlesize = ls*1.3)

band!(ax1, n_data_vals, chronological_halflife_ests[:,1,1], 
            chronological_halflife_ests[:,1,3],
        color = Makie.wong_colors()[2], alpha = alpha)
        
elem1 = lines!(ax1, n_data_vals, chronological_halflife_ests[:,1,2],
        color = Makie.wong_colors()[2], linewidth = lw)

band!(ax1, n_data_vals, chronological_halflife_ests[:,2,1], 
        chronological_halflife_ests[:,2,3],
            color = Makie.wong_colors()[1], alpha = alpha)
elem2 = lines!(ax1, n_data_vals, chronological_halflife_ests[:,2,2],
                color = Makie.wong_colors()[1], linewidth = lw)


band!(ax2, n_data_vals, rev_chronological_halflife_ests[:,1,1], 
        rev_chronological_halflife_ests[:,1,3],
    color = Makie.wong_colors()[2], alpha = alpha)
    
elem1 = lines!(ax2, n_data_vals, rev_chronological_halflife_ests[:,1,2],
    color = Makie.wong_colors()[2], linewidth = lw)

band!(ax2, n_data_vals, rev_chronological_halflife_ests[:,2,1], 
    rev_chronological_halflife_ests[:,2,3],
        color = Makie.wong_colors()[1], alpha = alpha)

elem2 = lines!(ax2, n_data_vals, rev_chronological_halflife_ests[:,2,2],
    color = Makie.wong_colors()[1], linewidth = lw)
  
    


    
band!(ax3, n_data_vals, rand_halflife_ests[:,1,1], 
                rand_halflife_ests[:,1,3],
    color = Makie.wong_colors()[2], alpha = alpha)
    
elem1 = lines!(ax3, n_data_vals, rand_halflife_ests[:,1,2],
    color = Makie.wong_colors()[2], linewidth = lw)

band!(ax3, n_data_vals, rand_halflife_ests[:,2,1], 
                rand_halflife_ests[:,2,3],
    color = Makie.wong_colors()[1], alpha = alpha)
    
elem2 = lines!(ax3, n_data_vals, rand_halflife_ests[:,2,2],
    color = Makie.wong_colors()[1], linewidth = lw)


Legend(fig[1,1], [elem1, elem2],
        ["cITN","LLIN"],
        rowgap = 0,
        halign = :right, valign = :top, 
        tellwidth = false, tellheight = false,
        margin = (5,5,5,5), labelsize = ls*0.9)
Legend(fig[1,2], [elem1, elem2],
        ["cITN","LLIN"],
        rowgap = 0,
        halign = :right, valign = :top, 
        tellwidth = false, tellheight = false,
        margin = (5,5,5,5), labelsize = ls*0.9)
Legend(fig[1,3], [elem1, elem2],
        ["cITN","LLIN"],
        rowgap = 0,
        halign = :right, valign = :top, 
        tellwidth = false, tellheight = false,
        margin = (5,5,5,5), labelsize = ls*0.9)
        

xlims!.([ax1,ax2,ax3], -1, 41)
ylims!.([ax1,ax2,ax3], -0.1, 6.8)

save(OUTPUT_PLOTS_DIR*"PaperFigures/TechnicalPaper/Sensitivity_Analysis.pdf", fig, pdf_version = "1.4")
fig