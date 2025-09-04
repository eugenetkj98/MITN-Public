"""
Author: Eugene Tan
Date Created: 8/5/2025
Last Updated: 8/5/2025
Make country level plots for paper
"""
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/read_toml.jl")

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

# %% Plot packages
using LaTeXStrings
using CairoMakie

####################################################
# %% Plotting Theme and general settings
####################################################
set_theme!(theme_ggplot2())

# Color settings
colors = [  colorant"#3C6E96", # cITN
            colorant"#C1434E", # LLIN
            colorant"#71366A", # PBO
            # colorant"#5E6F46", # G2
            colorant"#FDBB5F", # ROYAL
            ];

# General settings
fillalpha = 0.2
la = 0.5
lw = 2
africa_lw = 1.4
titlesize = 23
labelsize = 18

####################################################
# %% Import Country Draws
####################################################
# Max quarters to show in demography plot
max_age_quarters = 12

# Define list of countries to plot
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# Storage Variables
cumulative_net_age_demographics = Vector{Matrix{Float64}}(undef, length(filt_ISOs))
net_type_demographics = Vector{Matrix{Float64}}(undef, length(filt_ISOs))
net_names = Vector{Vector{String}}(undef, length(filt_ISOs))
# %%
Threads.@threads for ISO_i in ProgressBar(1:length(filt_ISOs))
    # Import sampled data
    ISO = filt_ISOs[ISO_i]

    pred_data = JLD2.load(OUTPUT_DIR*"multitype_predictions/$(ISO)_multitype_netcrop_predictions.jld2")

    NET_NAMES = pred_data["NET_NAMES"]
    n_net_types = length(NET_NAMES)

    # Calculate average net age matrix from imported samples
    A_BYNET = pred_data["admin0_A_MONTHLY_BYNET"]
    A_BYNET_TOTAL = sum(A_BYNET, dims = 3)[:,:,1]
    
    ##### Calculate age strata matrix for nets
    mean_net_age_matrix_BYNET = zeros(size(A_BYNET)[1], max_age_quarters, n_net_types)
    for i in 1:size(A_BYNET)[1] #current time
        for j in 1:i # All possible birth times
            for n in 1:n_net_types
                age_months = i-j
                age_quarter = (i-j)÷4
                age_index = min(age_quarter + 1, max_age_quarters)

                mean_net_age_matrix_BYNET[i,age_index,n] += A_BYNET[i,j,n]
            end
        end
    end

    mean_net_age_matrix_TOTAL = sum(mean_net_age_matrix_BYNET, dims = 3)[:,:]
    
    # Calculate cumulative net age matrix
    cumulative_net_age_matrix_BYNET = zeros(size(mean_net_age_matrix_BYNET))

    for strata_i in 1:size(mean_net_age_matrix_BYNET)[2]
        for n in 1:n_net_types
            cumulative_net_age_matrix_BYNET[:,strata_i,n] = sum(mean_net_age_matrix_BYNET[:,end-(strata_i-1):end,n], dims = 2)
        end
    end

    cumulative_net_age_matrix_TOTAL = sum(cumulative_net_age_matrix_BYNET, dims = 3)[:,:]

    # Save to variable
    cumulative_net_age_demographics[ISO_i] = cumulative_net_age_matrix_TOTAL
    net_type_demographics[ISO_i] = pred_data["admin0_Γ_MONTHLY_BYNET"]
    net_names[ISO_i] = NET_NAMES
end

####################################################
# %% Make Country Level Plot for Demography
####################################################
n_months = size(cumulative_net_age_demographics[1])[1]
fig = nothing
for ISO_i in 1:length(filt_ISOs)
    ISO = filt_ISOs[ISO_i]

    fig = Figure(size = (600,400))
    ax = Axis(fig[1,1],
                title = "$(ISO) Net Age Demography",
                titlesize = titlesize,
                xlabel = "Years", xlabelsize = labelsize,
                ylabel = "Net Crop (millions)", ylabelsize = labelsize,
                xticks = (1:12:n_months, string.(YEAR_NAT_START:YEAR_NAT_END)),
                xticklabelrotation = pi/2)
    legend_elems = Vector{Any}(undef, max_age_quarters)
    for strata_i in 1:max_age_quarters
        if strata_i == 1
            # legend_elems[strata_i] = band!(ax, 1:n_months, zeros(length(cumulative_net_age_demographics[ISO_i][:, strata_i]))./1e6, 
            #     cumulative_net_age_demographics[ISO_i][:, strata_i]./1e6,
            #     color = strata_i, colorrange = (1,max_age_quarters), colormap = Makie.reverse(cgrad(:roma)), alpha = 0.8)
            legend_elems[strata_i] = band!(ax, (6*12+1):n_months, zeros(length(cumulative_net_age_demographics[ISO_i][(6*12+1):n_months, strata_i]))./1e6, 
                cumulative_net_age_demographics[ISO_i][(6*12+1):n_months, strata_i]./1e6,
                color = strata_i, colorrange = (1,max_age_quarters), colormap = Makie.reverse(cgrad(:roma)), alpha = 0.8)
        else
            # legend_elems[strata_i] = band!(ax, 1:n_months, cumulative_net_age_demographics[ISO_i][:, strata_i-1]./1e6, 
            #     cumulative_net_age_demographics[ISO_i][:, strata_i]./1e6,
            #     color = strata_i, colorrange = (1,max_age_quarters), colormap = Makie.reverse(cgrad(:roma)), alpha = 0.8)
            legend_elems[strata_i] = band!(ax,(6*12+1):n_months, cumulative_net_age_demographics[ISO_i][(6*12+1):n_months, strata_i-1]./1e6, 
                cumulative_net_age_demographics[ISO_i][(6*12+1):n_months, strata_i]./1e6,
                color = strata_i, colorrange = (1,max_age_quarters), colormap = Makie.reverse(cgrad(:roma)), alpha = 0.8)
        end
    end

    Legend(fig[1,1], legend_elems,
            string.(1:max_age_quarters),
            "Net Age\n(Quarters)",
            tellheight = false, tellwidth = false,
            halign = :left, 
            valign = :top,
            margin = (10,5,10,5),
            patchsize = (15,8),
            rowgap = -3.5)

    # Save figure
    mkpath(OUTPUT_PLOTS_DIR*"PaperFigures/Country_ITN_Age_Demography/")
    save(OUTPUT_PLOTS_DIR*"PaperFigures/Country_ITN_Age_Demography/$(ISO)_age_demography.pdf", fig, pdf_version = "1.4")
end

fig
####################################################
# %% Make Country Level Plot for Type
####################################################
n_months = size(net_type_demographics[1])[1]
fig = nothing
for ISO_i in 1:length(filt_ISOs)
    ISO = filt_ISOs[ISO_i]
    net_type_demographic = net_type_demographics[ISO_i]
    type_proportion = (net_type_demographic./repeat(sum(net_type_demographic, dims = 2)[:,1],1,size(net_type_demographic)[2])).*100
    
    fig = Figure(size = (700,400))
    ax = Axis(fig[1,1],
                title = "$(ISO) Net Type Demography",
                titlesize = titlesize,
                xlabel = "Years", xlabelsize = labelsize,
                ylabel = "Net Crop Composition (%)", ylabelsize = labelsize,
                xticks = (1:12:n_months, string.(YEAR_NAT_START:YEAR_NAT_END)),
                xticklabelrotation = pi/2,
                yticks = 0:10:100)

    n_net_types = length(net_names[ISO_i])

    legend_elems = Vector{Any}(undef, n_net_types)
    
    for strata_i in 1:n_net_types
        if strata_i == 1
            # legend_elems[strata_i] = band!(ax, 1:n_months, zeros(length(type_proportion[:, strata_i])), 
            #             sum(type_proportion[:,1:strata_i], dims = 2)[:,1],
            #             color = (colors[strata_i],0.8))
            legend_elems[strata_i] = band!(ax, (6*12+1):n_months, zeros(length(type_proportion[(6*12+1):n_months, strata_i])), 
                        sum(type_proportion[(6*12+1):n_months,1:strata_i], dims = 2)[:,1],
                        color = (colors[strata_i],0.8))
        else
            # legend_elems[strata_i] = band!(ax, 1:n_months, sum(type_proportion[:,1:strata_i-1], dims = 2)[:,1],
            #             sum(type_proportion[:,1:strata_i], dims = 2)[:,1],
            #             color = (colors[strata_i], 0.8))
            legend_elems[strata_i] = band!(ax, (6*12+1):n_months, sum(type_proportion[(6*12+1):n_months,1:strata_i-1], dims = 2)[:,1],
                        sum(type_proportion[(6*12+1):n_months,1:strata_i], dims = 2)[:,1],
                        color = (colors[strata_i], 0.8))
        end
    end

    Legend(fig[1,2], legend_elems,
                net_names[ISO_i],
                "Net Type",
                tellheight = false, tellwidth = true,
                margin = (5,0,5,0),
                patchsize = (15,8),
                rowgap = 1)

    # Save figure
    mkpath(OUTPUT_PLOTS_DIR*"PaperFigures/Country_ITN_Type_Demography/")
    save(OUTPUT_PLOTS_DIR*"PaperFigures/Country_ITN_Type_Demography/$(ISO)_type_demography.pdf", fig, pdf_version = "1.4")
end

fig

####################################################
# %% Africa Net Age Demography
####################################################
continent_net_age_demographic = sum(cumulative_net_age_demographics)

fig = Figure(size = (600,400))
ax = Axis(fig[1,1],
            title = "Africa Net Age Demography",
            titlesize = titlesize,
            xlabel = "Years", xlabelsize = labelsize,
            ylabel = "Net Crop (millions)", ylabelsize = labelsize,
            xticks = (1:12:n_months, string.(YEAR_NAT_START:YEAR_NAT_END)),
            xticklabelrotation = pi/2)
legend_elems = Vector{Any}(undef, max_age_quarters)
for strata_i in 1:max_age_quarters
    if strata_i == 1
        # legend_elems[strata_i] = band!(ax, 1:n_months, zeros(length(continent_net_age_demographic[:, strata_i]))./1e6, 
        # continent_net_age_demographic[:, strata_i]./1e6,
        #     color = strata_i, colorrange = (1,max_age_quarters), colormap = Makie.reverse(cgrad(:roma)), alpha = 0.8)
        legend_elems[strata_i] = band!(ax, (6*12+1):n_months, zeros(length(continent_net_age_demographic[(6*12+1):n_months, strata_i]))./1e6, 
        continent_net_age_demographic[(6*12+1):n_months, strata_i]./1e6,
            color = strata_i, colorrange = (1,max_age_quarters), colormap = Makie.reverse(cgrad(:roma)), alpha = 0.8)
    else
        # legend_elems[strata_i] = band!(ax, 1:n_months, continent_net_age_demographic[:, strata_i-1]./1e6, 
        # continent_net_age_demographic[:, strata_i]./1e6,
        #     color = strata_i, colorrange = (1,max_age_quarters), colormap = Makie.reverse(cgrad(:roma)), alpha = 0.8)
        legend_elems[strata_i] = band!(ax, (6*12+1):n_months, continent_net_age_demographic[(6*12+1):n_months, strata_i-1]./1e6, 
        continent_net_age_demographic[(6*12+1):n_months, strata_i]./1e6,
            color = strata_i, colorrange = (1,max_age_quarters), colormap = Makie.reverse(cgrad(:roma)), alpha = 0.8)
    end
end

Legend(fig[1,1], legend_elems,
        string.(1:max_age_quarters),
        "Net Age\n(Quarters)",
        tellheight = false, tellwidth = false,
        halign = :left, 
        valign = :top,
        margin = (10,5,10,5),
        patchsize = (15,8),
        rowgap = -3.5)

fig

# %% Save figure
save(OUTPUT_PLOTS_DIR*"PaperFigures/Africa_age_demography.pdf", fig, pdf_version = "1.4")

####################################################
# %% Africa Net Age Demography
####################################################
continent_net_type_demographic = sum(net_type_demographics)

net_type_demographics
continent_proportion = (continent_net_type_demographic./repeat(sum(continent_net_type_demographic, dims = 2)[:,1],1,size(continent_net_type_demographic)[2])).*100
fig = Figure(size = (700,400))
ax = Axis(fig[1,1],
            title = "Africa Net Type Demography",
            titlesize = titlesize,
            xlabel = "Years", xlabelsize = labelsize,
            ylabel = "Net Type Composition (%)", ylabelsize = labelsize,
            xticks = (1:12:n_months, string.(YEAR_NAT_START:YEAR_NAT_END)),
            xticklabelrotation = pi/2,
            yticks = 0:10:100)

n_net_types = 4#5#length(net_names[ISO_i])

legend_elems = Vector{Any}(undef, n_net_types)

for strata_i in 1:n_net_types
    println(strata_i)
    if strata_i == 1
        # legend_elems[strata_i] = band!(ax, 1:n_months, zeros(length(continent_proportion[:, strata_i])), 
        #             sum(continent_proportion[:,1:strata_i], dims = 2)[:,1],
        #             color = (colors[strata_i],0.8))
        legend_elems[strata_i] = band!(ax, (6*12+1):n_months, zeros(length(continent_proportion[(6*12+1):n_months, strata_i])), 
                    sum(continent_proportion[(6*12+1):n_months,1:strata_i], dims = 2)[:,1],
                    color = (colors[strata_i],0.8))
    else
        # legend_elems[strata_i] = band!(ax, 1:n_months, sum(continent_proportion[:,1:strata_i-1], dims = 2)[:,1],
        #             sum(continent_proportion[:,1:strata_i], dims = 2)[:,1],
        #             color = (colors[strata_i], 0.8))
        legend_elems[strata_i] = band!(ax, (6*12+1):n_months, sum(continent_proportion[(6*12+1):n_months,1:strata_i-1], dims = 2)[:,1],
                    sum(continent_proportion[(6*12+1):n_months,1:strata_i], dims = 2)[:,1],
                    color = (colors[strata_i], 0.8))
    end
end

Legend(fig[1,2], legend_elems,
            # ["cITN","LLIN","PBO","G2","ROYAL"],
            ["cITN","LLIN","PBO","DAI"],
            "Net Type",
            tellheight = false, tellwidth = true,
            margin = (5,0,5,0),
            patchsize = (15,8),
            rowgap = 1)
fig

# %% Save figure
save(OUTPUT_PLOTS_DIR*"PaperFigures/Africa_type_demography.pdf", fig, pdf_version = "1.4")

####################################################
# %% Africa Net Age Demography (Next Gen vs Old Gen)
####################################################
continent_net_type_demographic = sum(net_type_demographics)

net_type_demographics
continent_proportion = (continent_net_type_demographic./repeat(sum(continent_net_type_demographic, dims = 2)[:,1],1,size(continent_net_type_demographic)[2])).*100
fig = Figure(size = (700,400))
ax = Axis(fig[1,1],
            title = "Africa Net Type Demography",
            titlesize = titlesize,
            xlabel = "Years", xlabelsize = labelsize,
            ylabel = "Net Type Composition (%)", ylabelsize = labelsize,
            xticks = (1:12:n_months, string.(YEAR_NAT_START:YEAR_NAT_END)),
            xticklabelrotation = pi/2,
            yticks = 0:10:100)

n_net_types = 4#5#length(net_names[ISO_i])

legend_elems = Vector{Any}(undef, 2)

legend_elems[1] = band!(ax, 1:n_months, zeros(length(continent_proportion[:, 1])), 
                    sum(continent_proportion[:,1:2], dims = 2)[:,1],
                    color = (colors[1],0.8))
legend_elems[2] = band!(ax, 1:n_months, sum(continent_proportion[:,1:2], dims = 2)[:,1],
                    sum(continent_proportion[:,1:end], dims = 2)[:,1],
                    color = (colors[2], 0.8))
    


Legend(fig[1,2], legend_elems,
            # ["cITN","LLIN","PBO","G2","ROYAL"],
            ["cITN+LLIN","Next Gen"],
            "Net Type",
            tellheight = false, tellwidth = true,
            margin = (17,0,5,0),
            patchsize = (15,8),
            rowgap = 1)

fig

# %% Save figure
save(OUTPUT_PLOTS_DIR*"PaperFigures/Africa_type_demography_twotype.pdf", fig, pdf_version = "1.4")