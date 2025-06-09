"""
Author: Eugene Tan
Date Created: 8/5/2025
Last Updated: 8/5/2025
Make country level plots for paper
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
using ProgressBars


# %% Maths packages
using LinearAlgebra
using StatsBase

# %% General useful functions
using DateConversions
using TS_filters
using NetLoss

# %% Plot packages
using LaTeXStrings
using CairoMakie

# %% Constants
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

# %% Plotting Theme and general settings
set_theme!(theme_ggplot2())

# Color settings
color_dict = Dict("cITN" => colorant"#E72A3D",
                    "LLIN" => colorant"#0082C7",
                    "PBO" => colorant"#538255",
                    "G2" => colorant"#538255",
                    "ROYAL" => colorant"#A9197B")


# General settings
fillalpha = 0.1
la = 0.5
lw = 2
africa_lw = 1.4

# %% Define list of countries to plot
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Calculate attrition curves

#Range of time vals to calculate attrition curve over
t_vals = 0:0.1:5

#########################################################################
# %% Make Africa shaped plot of attrition curves across all countries
#########################################################################

# Define index locations in Africa
subplot_lookup = Dict(   "MRT" => (1,2,"Mauritania"),
                        "ERI" => (1,8,"Eritrea"),
                        "GMB" => (2,1,"Gambia"),
                        "SEN" => (2,2,"Senegal"),
                        "GNB" => (2,3,"Guinea-Bissau"),
                        "MLI" => (2,4,"Mali"),
                        "NER" => (2,5,"Niger"),
                        "TCD" => (2,6,"Chad"),
                        "SDN" => (2,7,"Sudan"),
                        "DJI" => (2,8,"Djibouti"),
                        "SLE" => (3,1,"Sierra Leone"),
                        "GIN" => (3,2,"Guinea"),
                        "GHA" => (3,3,"Ghana"),
                        "BFA" => (3,4,"Burkina Faso"),
                        "BEN" => (3,5,"Benin"),
                        "CAF" => (3,6,"Centr. Afr. Rep."),
                        "SSD" => (3,7,"South Sudan"),
                        "ETH" => (3,8,"Ethiopia"),
                        "SOM" => (3,9,"Somalia"),
                        "LBR" => (4,2,"Liberia"),
                        "CIV" => (4,3,"Cote d'Ivoire"),
                        "TGO" => (4,4,"Togo"),
                        "NGA" => (4,5,"Nigeria"),
                        "CMR" => (4,6,"Cameroon"),
                        "RWA" => (4,7,"Rwanda"),
                        "KEN" => (4,8,"Kenya"),
                        "STP" => (5,3,"Sao Tome & Prin."),
                        "GAB" => (5,4,"Gabon"),
                        "COG" => (5,5,"Rep. of Congo"),
                        "GNQ" => (5,6,"Equatorial Guinea"),
                        "BDI" => (5,7,"Burundi"),
                        "UGA" => (5,8,"Uganda"),
                        "AGO" => (6,4,"Angola"),
                        "COD" => (6,5,"Dem. Rep. Congo"),
                        "ZMB" => (6,6,"Zambia"),
                        "TZA" => (6,7,"Tanzania"),
                        "NAM" => (7,3,"Namibia"),
                        "BWA" => (7,4,"Botswana"),
                        "ZWE" => (7,5,"Zimbabwe"),
                        "MWI" => (7,6,"Malawi"),
                        "MOZ" => (7,7,"Mozambique"),
                        "SWZ" => (7,8,"Eswatini"),
                        "COM" => (7,9,"Comoros"),
                        "MDG" => (7,10,"Madagascar"))
show_xlabel_ISOs = ["SLE","LBR","STP","NAM","BWA",
                    "ZWE","MWI","MOZ","SWZ","COM","MDG",
                    "UGA","SOM"]
show_ylabel_ISOs = ["MRT","ERI","GMB","SLE","LBR","STP",
                    "AGO","NAM"]

# Make figure
layout_res = (1500,1300)
fig = Figure(size = layout_res)

plot_axs = []
title_axs = []
for ISO_i in 1:length(filt_ISOs)
    ISO = filt_ISOs[ISO_i]
    x_idx, y_idx, country_name = subplot_lookup[ISO]
    ax = Axis(fig[2*x_idx,y_idx], width = layout_res[1]/14,
                # xlabel = "Years", 
                xticks = (0:1:5),
                # ylabel = "Metric",
                yticks = (0:0.2:1),
                ylabelsize = 20)
    xlims!(ax,-0.5,5.05)
    ylims!(ax,-0.02, 1.02)

    lb = Label(fig[2*(x_idx-1)+1,y_idx], "$(country_name)")
    
    # Hide labels as required
    if !(ISO ∈ show_xlabel_ISOs)
        ax.xticklabelsvisible = false
    end

    if !(ISO ∈ show_ylabel_ISOs)
        ax.yticklabelsvisible = false
    end
    push!(plot_axs, ax)
    push!(title_axs, lb)
end

# Super Labels
Label(fig[:,0],"Survival Rate", rotation = pi/2, fontsize = 40)
Label(fig[15,:],"Years", fontsize = 40)

fig
# Add plot lines for metric
for ISO_i in ProgressBar(1:length(filt_ISOs), leave = false)
    # Select ISO
    ISO = filt_ISOs[ISO_i]

    # %% Load Datasets
    input_dict = load(OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
    regression_dict = load(OUTPUT_REGRESSIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")

    # %%
    # Get Metadata and define bounds
    NET_NAMES = regression_dict["NET_NAMES"]
    n_nets = length(NET_NAMES)

    # Extract chain
    τ_chain = regression_dict["chain"][:,5:2:5+(length(NET_NAMES)-1)*2]
    κ_chain = regression_dict["chain"][:,6:2:6+(length(NET_NAMES)-1)*2]

    for net_type_i in 1:n_nets
        NET_NAME = NET_NAMES[net_type_i]
    
        # Calculate attrition curve samples
        attrition_curves = Matrix{Float64}(undef, size(τ_chain)[1], length(t_vals))
        attrition_curves_quantiles = Matrix{Float64}(undef, length(t_vals),3)

        for i in 1:size(τ_chain)[1]
            attrition_curves[i,:] = net_loss_compact.(t_vals, τ_chain[i,net_type_i], κ_chain[i,net_type_i])
        end

        for t in 1:length(t_vals)
            attrition_curves_quantiles[t,:] = quantile(attrition_curves[:,t], [0.025,0.5,0.975])
        end

        color_dict[NET_NAME]
        band!(plot_axs[ISO_i], t_vals, attrition_curves_quantiles[:,1], attrition_curves_quantiles[:,3],
                color = (color_dict[NET_NAME], fillalpha))
        lines!(plot_axs[ISO_i], t_vals, attrition_curves_quantiles[:,2],
                linewidth = lw, color = color_dict[NET_NAME])
    end

end

elem_1 = [LineElement(color = color_dict["cITN"], linewidth = lw)]
elem_2 = [LineElement(color = color_dict["LLIN"], linewidth = lw)]
# elem_3 = [LineElement(color = color_dict["PBO"], linewidth = lw)]
# elem_4 = [LineElement(color = color_dict["G2"], linewidth = lw)]
# elem_5 = [LineElement(color = color_dict["ROYAL"], linewidth = lw)]

Legend(fig[12, 1],
    # [elem_1, elem_2, elem_3, elem_4, elem_5],
    [elem_1, elem_2],
    # ["cITN", "LLIN", "PBO", "G2", "Royal"],
    ["cITN", "LLIN"],
    labelsize = 25)

# Save fig
save(OUTPUT_PLOTS_DIR*"PaperFigures/Africa_country_Attrition.pdf", fig, pdf_version = "1.4")
fig