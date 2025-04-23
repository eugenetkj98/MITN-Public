"""
Author: Eugene Tan
Date Created: 23/7/2024
Last Updated: 26/8/2024
Script to compile and tidy outputs from sensitivity analysis into a single file 
for each country
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import Public Packages
using JLD2
using CSV
using DataFrames
using ProgressBars

# %% Custom packages/modules
using PlottingFunctions

# Maths packages
using LinearAlgebra
using StatsBase

# %% Get ISO List
ISO_list = String.(CSV.read(raw"C:\Users\ETan\Documents\Prototype Analyses\itn-updated\datasets\ISO_list.csv", DataFrame)[:,1])
exclusion_ISOs = ["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD"]#["CPV","BWA","GNQ","DJI","ETH","SOM","ZAF","SSD", ""]
# ["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD"]
# GAB excluded
# %% Run Analysis configs
YEAR_START = 2000
YEAR_END = 2023 # Inclusive until end of 2021 December, excludes 2022

full_reg_dir = "outputs/regressions/crop/$(YEAR_START)_$(YEAR_END)/"
SA_reg_dir = "outputs/regressions/crop/sensitivity_analysis/$(YEAR_START)_$(YEAR_END)/compiled_results/" # Directory for sensitivity analysis regression files
mode_names = Dict("chronological" => "CSA",
                    "random" => "RSA",
                    "cross" => "CVA")
analysis_modes = ["chronological", "random", "cross"]

pythonplot()
theme(:vibrant)

# %% Storage lists for plots
plot_lists = Vector{Any}(undef, length(analysis_modes))
for sa_index in 1:length(analysis_modes)
    plot_lists[sa_index] = []
end

for i in 1:length(ISO_list)
    # Select ISO
    ISO = ISO_list[i]

    println("Plotting SA results for Country $(i) of $(length(ISO_list)) → $(ISO)...")
    
    if ISO ∈ exclusion_ISOs
        println("$(ISO) is on exclusion list. Moving to next country.")
        continue
    else
        # Define filenames
        full_reg_data_filename = "$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2"
        data_filename = "$(ISO)_$(YEAR_START)_$(YEAR_END)_sensitivity_analysis_results.jld2"
        
        for sa_index in 1:length(analysis_modes)
            sa_mode = analysis_modes[sa_index]
            
            # Get sensitivity analysis mode name
            mode_name = mode_names[sa_mode]
        
            # Import data
            # Sensitivity data analysis
            data_dict = load(SA_reg_dir*data_filename)
            results_matrix = data_dict[mode_name]
            n_entries = size(results_matrix)[1]
            n_net_types = size(results_matrix)[3]
        
            # Full regression data analysis results (For convergence comparison)
            full_reg_data_dict = load(full_reg_dir*full_reg_data_filename)
            full_reg_halflife_quantiles_BYNET = netlife_posterior_draws(full_reg_data_dict)[2]
        
            # Make plots
            lw = 1.8
            vlw = 0.8
            la = 0.6
            legend_labels = ["cITN", "LLIN"]
            fig = plot(xticks = 1:n_entries, xlabel = "Samples",
                        ylabel = "Net Halflife (years)",
                        title = "$(ISO) Sensitivity ($(sa_mode))", legend = :outerright)
            for j in 1:2#n_net_types
                plot!(fig, 1:n_entries, results_matrix[:,2,j], color = j, linewidth = lw, label = nothing)
                scatter!(fig, 1:n_entries, results_matrix[:,2,j], color = j, label = "$(legend_labels[j])")
                for k in 1:n_entries
                    plot!(fig, [k,k], results_matrix[k,[1,3],j], color = j, label = nothing,
                            linewidth = vlw, markershape = :rect, markersize = 4)
                end
        
                # Plot values for full regression version
                hline!(fig, [full_reg_halflife_quantiles_BYNET[2,j]], color = j,
                            label = nothing)
                hline!(fig, [full_reg_halflife_quantiles_BYNET[[1,3],j]], color = j, 
                        linestyle = :dash, label = nothing, alpha = la)
            end
            push!(plot_lists[sa_index], fig)
        end

        println("Plot collection complete for $(ISO). Ready for plotting.")
    end
end



# %%
fig = plot(plot_lists[1]..., layout = (5,7), size = (3840,2160), margin=3*Plots.mm)
savefig(fig, "SA_plots_chronological.pdf")

# %%
fig = plot(plot_lists[2]..., layout = (5,7), size = (3840,2160), margin=3*Plots.mm)
savefig(fig, "SA_plots_random.pdf")

# %%
fig = plot(plot_lists[3]..., layout = (5,7), size = (3840,2160), margin=3*Plots.mm)
savefig(fig, "SA_plots_cross.pdf")

