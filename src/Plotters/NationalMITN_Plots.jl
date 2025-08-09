"""
Author: Eugene Tan
Date Created: 20/3/2025
Last Updated: 20/3/2025
DEPRECATED IN MOVE TO USE MAKIE PLOTTING BACKEND!
Helper Function to generate nice plots using posterior draws of MITN outputs
"""

module NationalMITN_Plots
export plot_nat_timeseries
export plot_netcrop_demography
export plot_npc_demography
export plot_netcrop_bytype
export plot_npc_bytype
export plot_attrition_curves

# %% Import Public Packages
using DataFrames
using Missings
using JLD2
using CSV

# %% Maths packages
using LinearAlgebra
using StatsBase

# %% Plot packages
using LaTeXStrings
using Plots
using Measures

# %% MITN Model packages
using NetAccessModel
using DateConversions
using NetLoss


# %% Function to plot national time series
function plot_nat_timeseries(input_dict, net_access_input_dict, post_snf)
    # Plot visual settings
    pythonplot()
    theme(:vibrant)
    fillalpha = 0.15
    lw = 1.2
    ms = 4
    colors = [  colorant"#0082C7", # Net Crop
                colorant"#00976A", # NPC
                colorant"#E72A3D", # Access
                colorant"#45332C", # Guidelines
                ]

    #############################
    # %% Get metadata
    #############################
    ISO = input_dict["ISO"]
    MONTHS_MONTHLY = input_dict["MONTHS_MONTHLY"]
    YEARS_ANNUAL = input_dict["YEARS_ANNUAL"]
    n_months = length(MONTHS_MONTHLY)

    NET_CROP_SURVEY_MONTHLY = input_dict["NET_CROP_MONTHLY"]
    NPC_SURVEY_MONTHLY = input_dict["HOUSEHOLD_NPC_MONTHLY"]
    #############################
    # %% NPC and Access Time series
    #############################
    # Extract raw outputs
    POPULATION_MONTHLY = post_snf["POPULATION_MONTHLY"]
    Γ_MONTHLY_samples_BYNET = post_snf["Γ_MONTHLY_samples_BYNET"]
    λ_access_samples = post_snf["λ_access_samples"]

    # Get survey estimates of access
    λ_SURVEY_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
    if sum(.!ismissing.(NET_CROP_SURVEY_MONTHLY)) > 0
        national_H_aggregated = net_access_input_dict["H_aggregated"]
        national_access_aggregated = zeros(size(national_H_aggregated)[1])
        for survey_i in 1:length(national_access_aggregated)
            national_access_aggregated[survey_i] = sum(H_to_access(national_H_aggregated[survey_i,:,:]))
        end
        access_survey_monthidx = net_access_input_dict["survey_monthidx_aggregated"]
        # access_survey_monthidx = access_survey_monthidx[findall(.!ismissing.(access_survey_monthidx))]
        
        λ_SURVEY_MONTHLY[access_survey_monthidx] .= national_access_aggregated
    end

    # Sum across all net types
    Γ_MONTHLY_samples_TOTAL = sum(Γ_MONTHLY_samples_BYNET, dims = 3)[:,:,1]
    NPC_MONTHLY_samples_TOTAL = Γ_MONTHLY_samples_TOTAL./repeat(POPULATION_MONTHLY, 1, size(Γ_MONTHLY_samples_TOTAL)[1])'

    # Calculate coverage metrics with lower and upper 95% CI bounds
    Γ_MONTHLY_TOTAL = zeros(n_months, 3)
    NPC_MONTHLY_TOTAL = zeros(n_months, 3)
    λ_MONTHLY_TOTAL = zeros(n_months, 3)
    for i in 1:n_months
        Γ_MONTHLY_TOTAL[i,[1,3]] = quantile(Γ_MONTHLY_samples_TOTAL[:,i], [0.025, 0.975])
        Γ_MONTHLY_TOTAL[i,2] = mean(Γ_MONTHLY_samples_TOTAL[:,i])
        
        NPC_MONTHLY_TOTAL[i,[1,3]] = quantile(NPC_MONTHLY_samples_TOTAL[:,i], [0.025, 0.975])
        NPC_MONTHLY_TOTAL[i,2] = mean(NPC_MONTHLY_samples_TOTAL[:,i])
        
        λ_MONTHLY_TOTAL[i,[1,3]] = quantile(λ_access_samples[findall(.!isnan.(λ_access_samples[:,i])),i], [0.025, 0.975])
        λ_MONTHLY_TOTAL[i,2] = mean(λ_access_samples[findall(.!isnan.(λ_access_samples[:,i])),i])
    end

    # %% Make Plot
    fig = plot(title = "$(ISO) Net Coverage", 
                xlims = (-1,290), ylims = (-100000, quantile(Γ_MONTHLY_TOTAL[:], 0.99)*1.2)./1e6,
                xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90,
                xlabel = "Year", ylabel = "Net Crop", legend = :topleft)
    fig_2 = twinx(fig)
    plot!(fig_2, ylabel = "NPC, Access", xlims = (-1,290))

    # Annual guidelines
    vline!(fig, MONTHS_MONTHLY[1:12:end], 
            linecolor = colors[4], linestyle = :dash, linewidth = 1, linealpha = 0.2,
            label = nothing)

    # Raw Survey Data
    scatter!(fig, MONTHS_MONTHLY, NET_CROP_SURVEY_MONTHLY./1e6,
                markersize = ms, markercolor = colors[1],
                label = nothing)
    scatter!(fig_2, MONTHS_MONTHLY, NPC_SURVEY_MONTHLY,
                markersize = ms, markercolor = colors[2], ylims = (-0.01, 1.05),
                label = nothing)
    scatter!(fig_2, MONTHS_MONTHLY, λ_SURVEY_MONTHLY,
                markersize = ms, markercolor = colors[3], ylims = (-0.01, 1.05),
                label = nothing)

    # Plot empty Time series for legend
    plot!(fig, NaN.*MONTHS_MONTHLY, Γ_MONTHLY_TOTAL[:,2]./1e6,
            linealpha = 1, label = "Net Crop",
            linecolor = colors[1], linewidth = lw)
    plot!(fig, NaN.*MONTHS_MONTHLY, NPC_MONTHLY_TOTAL[:,2],
            linealpha = 1, label = "NPC",
            linecolor = colors[2], linewidth = lw)
    plot!(fig, NaN.*MONTHS_MONTHLY, λ_MONTHLY_TOTAL[:,2],
            linealpha = 1, label = "Access",
            linecolor = colors[3], linewidth = lw)

    # Plot Time series
    plot!(fig, MONTHS_MONTHLY, Γ_MONTHLY_TOTAL[:,1]./1e6, fillrange = Γ_MONTHLY_TOTAL[:,3]./1e6,
            linealpha = 0, fillalpha = fillalpha, label = nothing,
            linecolor = colors[1], fillcolor = colors[1])
    plot!(fig, MONTHS_MONTHLY, Γ_MONTHLY_TOTAL[:,2]./1e6,
            linealpha = 1, label = nothing,
            linecolor = colors[1], linewidth = lw)

    plot!(fig_2, MONTHS_MONTHLY, NPC_MONTHLY_TOTAL[:,1], fillrange = NPC_MONTHLY_TOTAL[:,3],
            linealpha = 0, fillalpha = fillalpha, label = nothing,
            linecolor = colors[2], fillcolor = colors[2], ylims = (-0.01, 1.05))
    plot!(fig_2, MONTHS_MONTHLY, NPC_MONTHLY_TOTAL[:,2],
            linealpha = 1, label = nothing,
            linecolor = colors[2], linewidth = lw, ylims = (-0.01, 1.05))

    plot!(fig_2, MONTHS_MONTHLY, λ_MONTHLY_TOTAL[:,1], fillrange = λ_MONTHLY_TOTAL[:,3],
            linealpha = 0, fillalpha = fillalpha, label = nothing,
            linecolor = colors[3], fillcolor = colors[3], ylims = (-0.01, 1.05))
    plot!(fig_2, MONTHS_MONTHLY, λ_MONTHLY_TOTAL[:,2],
            linealpha = 1, label = nothing,
            linecolor = colors[3], linewidth = lw, ylims = (-0.01, 1.05))

    return fig
end

# %% Plot Cumulative Net demography
function plot_netcrop_demography(input_dict, post_net_demography_mean)
    # Plot visual settings
    pythonplot()
    theme(:vibrant)
    fillalpha = 0.5
    #############################
    # %% Get metadata
    #############################
    ISO = input_dict["ISO"]
    MONTHS_MONTHLY = input_dict["MONTHS_MONTHLY"]
    YEARS_ANNUAL = input_dict["YEARS_ANNUAL"]
    POPULATION_MONTHLY = input_dict["POPULATION_MONTHLY"]
    n_months = length(MONTHS_MONTHLY)
    age_bins = parse.(Int64, names(post_net_demography_mean)[6:end])

    NET_CROP_SURVEY_MONTHLY = input_dict["NET_CROP_MONTHLY"]

    # Calculate net age breakdown during given month
    net_crop_demography = zeros(n_months, length(age_bins))
    for year in YEARS_ANNUAL
        for month in 1:12
            monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEARS_ANNUAL[1])
            slice = post_net_demography_mean[(post_net_demography_mean.ISO .== ISO) .&
                                        (post_net_demography_mean.YEAR .== year) .&
                                        (post_net_demography_mean.MONTH .== month),:]
            net_crop_demography[monthidx,:] = sum(Matrix(slice[:,6:end]), dims = 1)[1,:].*POPULATION_MONTHLY[monthidx]
        end
    end

    # Calculate cumulative demography
    cumul_net_crop_demography = zeros(size(net_crop_demography))
    for i in 1:size(net_crop_demography)[2]
        cumul_net_crop_demography[:,i] = sum(net_crop_demography[:,1:i], dims = 2)
    end

    quarterly_cumul_net_crop_demography = cumul_net_crop_demography[:, 1:4:end]


    # Make demography plot
    fig = plot(title = "$(ISO) Net Demography", 
                xlims = (-1,290), ylims = (-100000, quantile(cumul_net_crop_demography[:], 0.99)*1.2)./1e6,
                xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90,
                xlabel = "Year", ylabel = "Net Crop (mil)", legend = :topleft)
    colors = palette(:roma, size(quarterly_cumul_net_crop_demography)[2])

    for i in 1:size(quarterly_cumul_net_crop_demography)[2]
        if i == 1
            plot!(fig, MONTHS_MONTHLY, zeros(length(MONTHS_MONTHLY))./1e6,
                fillrange = quarterly_cumul_net_crop_demography[:,i]./1e6,
                linealpha = 0, linewidth = 0, fillalpha = fillalpha,
                label = "$i",
                fillcolor = colors[i])
        else
            plot!(fig, MONTHS_MONTHLY, quarterly_cumul_net_crop_demography[:,i-1]./1e6,
                fillrange = quarterly_cumul_net_crop_demography[:,i]./1e6,
                linealpha = 0, linewidth = 0, fillalpha = fillalpha,
                label = "$i",
                fillcolor = colors[i])
        end
    end
    scatter!(fig, MONTHS_MONTHLY, NET_CROP_SURVEY_MONTHLY./1e6,
                label = "Survey", markersize = 6, markercolor = colorant"#3E3936",
                marker = :xcross)

    return fig
end

function plot_npc_demography(input_dict, post_net_demography_mean)
    # Plot visual settings
    pythonplot()
    theme(:vibrant)
    fillalpha = 0.5
    #############################
    # %% Get metadata
    #############################
    ISO = input_dict["ISO"]
    MONTHS_MONTHLY = input_dict["MONTHS_MONTHLY"]
    YEARS_ANNUAL = input_dict["YEARS_ANNUAL"]
    n_months = length(MONTHS_MONTHLY)
    age_bins = parse.(Int64, names(post_net_demography_mean)[6:end])

    NPC_SURVEY_MONTHLY = input_dict["HOUSEHOLD_NPC_MONTHLY"]

    # Calculate net age breakdown during given month
    npc_demography = zeros(n_months, length(age_bins))
    for year in YEARS_ANNUAL
        for month in 1:12
            monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEARS_ANNUAL[1])
            slice = post_net_demography_mean[(post_net_demography_mean.ISO .== ISO) .&
                                        (post_net_demography_mean.YEAR .== year) .&
                                        (post_net_demography_mean.MONTH .== month),:]
            npc_demography[monthidx,:] = sum(Matrix(slice[:,6:end]), dims = 1)[1,:]
        end
    end

    # Calculate cumulative demography
    cumul_npc_demography = zeros(size(npc_demography))
    for i in 1:size(npc_demography)[2]
        cumul_npc_demography[:,i] = sum(npc_demography[:,1:i], dims = 2)
    end

    quarterly_cumul_npc_demography = cumul_npc_demography[:, 1:4:end]


    # Make demography plot
    fig = plot(title = "$(ISO) NPC Demography", 
                xlims = (-1,290), ylims = (-0.01, 1.05),
                xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90,
                xlabel = "Year", ylabel = "NPC", legend = :topleft)
    colors = palette(:roma, size(quarterly_cumul_npc_demography)[2])

    for i in 1:size(quarterly_cumul_npc_demography)[2]
        if i == 1
            plot!(fig, MONTHS_MONTHLY, zeros(length(MONTHS_MONTHLY)),
                fillrange = quarterly_cumul_npc_demography[:,i],
                linealpha = 0, linewidth = 0, fillalpha = fillalpha,
                label = "$i",
                fillcolor = colors[i])
        else
            plot!(fig, MONTHS_MONTHLY, quarterly_cumul_npc_demography[:,i-1],
                fillrange = quarterly_cumul_npc_demography[:,i],
                linealpha = 0, linewidth = 0, fillalpha = fillalpha,
                label = "$i",
                fillcolor = colors[i])
        end
    end
    scatter!(fig, MONTHS_MONTHLY, NPC_SURVEY_MONTHLY,
                label = "Survey", markersize = 6, markercolor = colorant"#3E3936",
                marker = :xcross)

    return fig
end

# %% Net Demography by Net Type
function plot_netcrop_bytype(input_dict, post_snf)
    # Plot visual settings
    pythonplot()
    theme(:vibrant)
    fillalpha = 0.4
    lw = 1.2
    # Define Net colors
    colors = [  colorant"#005684",
                colorant"#00976A",
                colorant"#E72A3D",
                colorant"#F7B801",
                colorant"#7018B3",
                ]

    #############################
    # %% Get metadata
    #############################
    ISO = input_dict["ISO"]
    MONTHS_MONTHLY = input_dict["MONTHS_MONTHLY"]
    YEARS_ANNUAL = input_dict["YEARS_ANNUAL"]
    NET_NAMES = input_dict["NET_NAMES"]

    NET_CROP_SURVEY_MONTHLY = input_dict["NET_CROP_MONTHLY"]

    #############################
    # %% NPC and Access Time series
    #############################
    # Extract raw outputs
    Γ_MONTHLY_samples_BYNET = post_snf["Γ_MONTHLY_samples_BYNET"]



    # Get average time series
    Γ_MONTHLY_mean_BYNET = mean(Γ_MONTHLY_samples_BYNET, dims = 1)[1,:,:]
    cumul_Γ_MONTHLY_mean_BYNET = zeros(size(Γ_MONTHLY_mean_BYNET)[1], size(Γ_MONTHLY_mean_BYNET)[2])
    for i in 1:size(cumul_Γ_MONTHLY_mean_BYNET)[2]
        cumul_Γ_MONTHLY_mean_BYNET[:,i] = sum(Γ_MONTHLY_mean_BYNET[:, 1:i], dims = 2)
    end

    # Make plot
    fig = plot(title = "$(ISO) Net Crop", 
                xlims = (-1,290), ylims = (-100000, quantile(cumul_Γ_MONTHLY_mean_BYNET[:], 0.99)*1.2)./1e6,
                xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90,
                xlabel = "Year", ylabel = "Net Crop (mil)", legend = :topleft)
    vline!(fig, MONTHS_MONTHLY[1:12:end], 
                linecolor = colorant"#45332C", linestyle = :dash, linewidth = 1, linealpha = 0.2,
                label = nothing)
    for i in size(cumul_Γ_MONTHLY_mean_BYNET)[2]:-1:1
        if i == 1
            plot!(fig, MONTHS_MONTHLY, zeros(size(cumul_Γ_MONTHLY_mean_BYNET)[1])./1e6,
                    fillrange = cumul_Γ_MONTHLY_mean_BYNET[:,i]./1e6,
                    linewidth = 0, linealpha = 0,
                    fillcolor = colors[i], fillalpha = fillalpha,
                    label = NET_NAMES[i])
            plot!(fig, MONTHS_MONTHLY, cumul_Γ_MONTHLY_mean_BYNET[:,i]./1e6,
                    linewidth = lw, linecolor = colors[i], label = nothing)
        else
            plot!(fig, MONTHS_MONTHLY, cumul_Γ_MONTHLY_mean_BYNET[:,i-1]./1e6,
                    fillrange = cumul_Γ_MONTHLY_mean_BYNET[:,i]./1e6,
                    linewidth = 0, linealpha = 0,
                    fillcolor = colors[i], fillalpha = fillalpha,
                    label = NET_NAMES[i])
            plot!(fig, MONTHS_MONTHLY, cumul_Γ_MONTHLY_mean_BYNET[:,i]./1e6,
                    linewidth = lw, linecolor = colors[i], label = nothing)
        end
    end

    scatter!(fig, MONTHS_MONTHLY, NET_CROP_SURVEY_MONTHLY./1e6,
                    label = "Survey", markersize = 6, markercolor = colorant"#3E3936",
                    marker = :xcross)

    return fig
end

function plot_npc_bytype(input_dict, post_snf)
    # Plot visual settings
    pythonplot()
    theme(:vibrant)
    fillalpha = 0.4
    lw = 1.2
    # Define Net colors
    colors = [  colorant"#005684",
                colorant"#00976A",
                colorant"#E72A3D",
                colorant"#F7B801",
                colorant"#7018B3",
                ]

    #############################
    # %% Get metadata
    #############################
    ISO = input_dict["ISO"]
    MONTHS_MONTHLY = input_dict["MONTHS_MONTHLY"]
    YEARS_ANNUAL = input_dict["YEARS_ANNUAL"]
    NET_NAMES = input_dict["NET_NAMES"]

    NPC_SURVEY_MONTHLY = input_dict["HOUSEHOLD_NPC_MONTHLY"]

    #############################
    # %% NPC and Access Time series
    #############################
    # Extract raw outputs
    POPULATION_MONTHLY = post_snf["POPULATION_MONTHLY"]
    Γ_MONTHLY_samples_BYNET = post_snf["Γ_MONTHLY_samples_BYNET"]

    # Get average time series
    Γ_MONTHLY_mean_BYNET = mean(Γ_MONTHLY_samples_BYNET, dims = 1)[1,:,:]
    cumul_NPC_MONTHLY_mean_BYNET = zeros(size(Γ_MONTHLY_mean_BYNET)[1], size(Γ_MONTHLY_mean_BYNET)[2])
    for i in 1:size(cumul_NPC_MONTHLY_mean_BYNET)[2]
        cumul_NPC_MONTHLY_mean_BYNET[:,i] = sum(Γ_MONTHLY_mean_BYNET[:, 1:i], dims = 2)./POPULATION_MONTHLY
    end

    # Make plot
    fig = plot(title = "$(ISO) NPC", 
                xlims = (-1,290), ylims = (-0.01, 1.05),
                xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90,
                xlabel = "Year", ylabel = "NPC", legend = :topleft)
    vline!(fig, MONTHS_MONTHLY[1:12:end], 
                linecolor = colorant"#45332C", linestyle = :dash, linewidth = 1, linealpha = 0.2,
                label = nothing)
    for i in size(cumul_NPC_MONTHLY_mean_BYNET)[2]:-1:1
        if i == 1
            plot!(fig, MONTHS_MONTHLY, zeros(size(cumul_NPC_MONTHLY_mean_BYNET)[1]),
                    fillrange = cumul_NPC_MONTHLY_mean_BYNET[:,i],
                    linewidth = 0, linealpha = 0,
                    fillcolor = colors[i], fillalpha = fillalpha,
                    label = NET_NAMES[i])
            plot!(fig, MONTHS_MONTHLY, cumul_NPC_MONTHLY_mean_BYNET[:,i],
                    linewidth = lw, linecolor = colors[i], label = nothing)
        else
            plot!(fig, MONTHS_MONTHLY, cumul_NPC_MONTHLY_mean_BYNET[:,i-1],
                    fillrange = cumul_NPC_MONTHLY_mean_BYNET[:,i],
                    linewidth = 0, linealpha = 0,
                    fillcolor = colors[i], fillalpha = fillalpha,
                    label = NET_NAMES[i])
            plot!(fig, MONTHS_MONTHLY, cumul_NPC_MONTHLY_mean_BYNET[:,i],
                    linewidth = lw, linecolor = colors[i], label = nothing)
        end
    end

    scatter!(fig, MONTHS_MONTHLY, NPC_SURVEY_MONTHLY,
                    label = "Survey", markersize = 6, markercolor = colorant"#3E3936",
                    marker = :xcross)

    return fig
end

# %% Plotting function for mean attrition curves
function plot_attrition_curves(regression_dict)
    # Plot visual settings
    pythonplot()
    theme(:vibrant)
    fillalpha = 0.12
    lw1 = 1.4
    lw2 = 0.8

    # Define Net colors
    colors = [  colorant"#005684",
                colorant"#00976A",
                colorant"#E72A3D",
                colorant"#F7B801",
                colorant"#7018B3",
                ]

    # Get Metadata and define bounds
    NET_NAMES = regression_dict["NET_NAMES"]
    t_vals = 0:0.01:5

    # Extract chain
    τ_chain = regression_dict["chain"][:,5:2:5+(length(NET_NAMES)-1)*2]
    κ_chain = regression_dict["chain"][:,6:2:6+(length(NET_NAMES)-1)*2]

    # Calculate attrition curves from posterior sampled chain
    attrition_curves = zeros(length(t_vals), size(τ_chain)[1], length(NET_NAMES))
    for i in 1:size(τ_chain)[1]
        for j in 1:length(NET_NAMES)
            attrition_curves[:,i,j] = net_loss_compact.(t_vals, τ_chain[i,j], κ_chain[i,j])
        end
    end

    # Summarise sample attrition curves into upper and lower 95% bounds for each time value
    attrition_curves_bounds = zeros(length(t_vals), length(NET_NAMES), 3)
    for t in 1:length(t_vals)
        for j in 1:length(NET_NAMES)
            attrition_curves_bounds[t,j,:] = quantile(attrition_curves[t,:,j], [0.025, 0.5, 0.975])
        end
    end

    fig = plot(title = "$(regression_dict["ISO"]) Net Attrition", 
                xlims = (-0.1, 4.1), ylims = (-0.05, 1.05),
                xlabel = "Years", ylabels = "Survival Rate")

    for i in 1:length(NET_NAMES)
        plot!(fig, t_vals, attrition_curves_bounds[:,i,1], fillrange = attrition_curves_bounds[:,i,3],
                linewidth = 0, linealpha = 0, fillalpha = fillalpha,
                fillcolor = colors[i], label = nothing)
        plot!(fig, t_vals, attrition_curves_bounds[:,i,2],
                linewidth = lw1, linecolor = colors[i], label = NET_NAMES[i])
        plot!(fig, t_vals, attrition_curves_bounds[:,i,[1,3]],
                linewidth = lw2, linecolor = colors[i], linestyle = :dash,
                label = nothing)
    end
    
    return fig
end

end

