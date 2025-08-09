"""
Author: Eugene Tan
Date Created: 13/8/2024
Last Updated: 26/8/2024
DEPRECATED IN MOVE TO USE MAKIE PLOTTING BACKEND!
Old plotting functions for showing various SNF outputs. No longer relevant. 
"""
module PlottingFunctions
export regression_timeseries_plot
export netlife_posterior_draws
export lifecurve_plots

# %% Import Public Packages
using DataFrames
using Missings

# %% Maths packages
using LinearAlgebra
using StatsBase

# %% Import SNF modules
using NetCropModel
using NetAccessModel
using NetLoss

# %% Plot packages
using LaTeXStrings
using Plots
using Measures

# %% Plotting Functions
"""
    regression_timeseries_plot(input_dict::Dict, regression_dict::Dict, net_access_input_dict::Dict, net_access_chain::Dict, n_samples = 100, kwargs...)
Takes in data extraction and regression results dictionaries from pipeline functions and returns a tuple of figures.
- A plot of net crop, net access, and inferred distribution proportion `monthly_p` with confidence bounds
- A plot of net crop stratified by age measured in quarters, with a maximum age of `max_age_quarters` *(default = 16)*
Other notable arguments
- A scatter plot comparing the RMSE of the BV model with uniform distribution, and new MAP ITN model.
- `n_samples`: Number of samples to randomly draw from calculate time series. *(default = 100)*
"""
function regression_timeseries_plot(input_dict, regression_dict, 
                                    net_access_input_dict, net_access_chain;# BV_dict;
                                    n_samples = 100,
                                    lw = 2, fillalpha = 0.2, msize = 6,
                                    max_age_quarters = 12, country_name = nothing)
    pythonplot()
    theme(:vibrant)

    # Crop Regression input data
    ISO = input_dict["ISO"]
    YEARS_ANNUAL = input_dict["YEARS_ANNUAL"]
    MONTHS_MONTHLY = input_dict["MONTHS_MONTHLY"]
    DELIVERIES_ANNUAL = input_dict["DELIVERIES_ANNUAL"]
    DISTRIBUTION_ANNUAL = input_dict["DISTRIBUTION_ANNUAL"]
    NET_CROP_MONTHLY = input_dict["NET_CROP_MONTHLY"]
    POPULATION_MONTHLY = input_dict["POPULATION_MONTHLY"]
    NET_NAMES = input_dict["NET_NAMES"]
    n_net_types = length(NET_NAMES)

    # Get number of missing net values
    n_missing_nets_vals = sum(ismissing.(DISTRIBUTION_ANNUAL[:,1]))

    # Crop Regression output data
    chain = regression_dict["chain"]
    monthly_p = regression_dict["monthly_p"]
    chain_UNIF = regression_dict["chain_epochs"][1]
    monthly_p_UNIF = regression_dict["monthly_weights"][1]

    # Access Regression output data
    ρ_chain_df = net_access_chain["ρ_chain_df"]
    μ_chain_df = net_access_chain["μ_chain_df"]
    p_h = mean(net_access_input_dict["p_h_aggregated"], dims = 1)[:]

    # Calculate reference values for National Access from Surveys for scatter plot
    national_H_aggregated = net_access_input_dict["H_aggregated"]
    national_access_aggregated = zeros(size(national_H_aggregated)[1])
    for survey_i in 1:length(national_access_aggregated)
        national_access_aggregated[survey_i] = sum(H_to_access(national_H_aggregated[survey_i,:,:]))
    end

    access_survey_monthidx = net_access_input_dict["survey_monthidx_aggregated"]
    access_survey_monthidx = access_survey_monthidx[findall(.!ismissing.(access_survey_monthidx))]
    NET_ACCESS_SURVEY_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
    NET_ACCESS_SURVEY_MONTHLY[access_survey_monthidx] .= national_access_aggregated
    
    # # Import BV model fits
    # ISO_BV = BV_dict["ISO_list"]
    # bv_iso_idx = findfirst(ISO_BV.==ISO)
    iso_in_bv = false
    # if !isnothing(bv_iso_idx)
    #     iso_in_bv = true
    # end

    # NPC_MONTHLY_BV_samples = []
    # ACCESS_MONTHLY_BV_samples = []
    # if iso_in_bv
    #     NPC_MONTHLY_BV_samples = BV_dict["NPC_MONTHLY_DRAWS"][bv_iso_idx,:,:]
    #     ACCESS_MONTHLY_BV_samples = BV_dict["ACCESS_MONTHLY_DRAWS"][bv_iso_idx,:,:]
    # end

    ##### Extract Net Crop MCMC parameters

    ### Part 1: Final Regressed Values
    # Randomly generate sample indexes to sample from chain
    chain_length = size(chain)[1]
    sample_idxs = sample(1:chain_length, n_samples, replace = false)

    # Extract MCMC draws for parameters
    ϕ_posterior_draws = chain[sample_idxs, :ϕ]
    α_init_posterior_draws = chain[sample_idxs,:α_init]
    α_LLIN_posterior_draws = chain[sample_idxs,:α_LLIN]
    
    b_net_posterior_draws = Matrix(DataFrame(chain)[:,5:2:5+2*(n_net_types-1)])[sample_idxs, :]
    k_net_posterior_draws = Matrix(DataFrame(chain)[:,6:2:6+2*(n_net_types-1)])[sample_idxs, :]
    if n_missing_nets_vals > 0
        n_missing_nets_posterior_draws = Matrix(DataFrame(chain)[:,(4+2*(n_net_types)+1):end])[sample_idxs, :]
    else
        # No missing data, so just feed in a zero matrix with no columns
        n_missing_nets_posterior_draws = zeros(size(chain)[1], 0)
    end
    

    ##### Generate Net Crop Trajectories
    Γ_MONTHLY_samples_BYNET = zeros(n_samples, length(MONTHS_MONTHLY), n_net_types)
    A_samples_BYNET = zeros(n_samples, length(MONTHS_MONTHLY), length(MONTHS_MONTHLY), n_net_types)

    for i in 1:n_samples
        # Select MCMC posterior draw parameters
        ϕ = ϕ_posterior_draws[i]
        α_init = α_init_posterior_draws[i]
        α_LLIN = α_LLIN_posterior_draws[i]
        b_nets = b_net_posterior_draws[i,:]
        k_nets = k_net_posterior_draws[i,:]
        missing_nets = n_missing_nets_posterior_draws[i,:]

        Γ_MONTHLY_samples_BYNET[i,:,:,:], A_samples_BYNET[i,:,:,:] = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                                                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                                                                    ϕ, b_nets, k_nets,
                                                                                    α_init, α_LLIN,
                                                                                    missing_nets; 
                                                                                    monthly_p = monthly_p,
                                                                                    return_age = true)
    end

    # Scale by to units of millions
    Γ_MONTHLY_samples_BYNET = Γ_MONTHLY_samples_BYNET./1e6
    A_samples_BYNET = A_samples_BYNET./1e6
    # Get Totals
    Γ_MONTHLY_samples_TOTAL = sum(Γ_MONTHLY_samples_BYNET, dims = 3)[:,:,1]
    A_samples_TOTAL = sum(A_samples_BYNET, dims = 4)[:,:,:,1]

    ### Part 2: UNIF Distribution Comparison Values (Not the real BV values, but just MAP ITN with no uniform distribution assumption)
    # Randomly generate sample indexes to sample from chain
    chain_length_UNIF = size(chain_UNIF)[1]
    sample_idxs_UNIF = sample(1:chain_length_UNIF, n_samples, replace = false)

    # Extract MCMC draws for parameters
    ϕ_posterior_draws_UNIF = chain_UNIF[sample_idxs_UNIF, :ϕ]
    α_LLIN_posterior_draws_UNIF = chain_UNIF[sample_idxs_UNIF,:α_LLIN]
    α_init_posterior_draws_UNIF = chain_UNIF[sample_idxs_UNIF,:α_init]

    b_net_posterior_draws_UNIF = Matrix(DataFrame(chain_UNIF)[:,5:2:5+2*(n_net_types-1)])[sample_idxs_UNIF, :]
    k_net_posterior_draws_UNIF = Matrix(DataFrame(chain_UNIF)[:,6:2:6+2*(n_net_types-1)])[sample_idxs_UNIF, :]
    
    if n_missing_nets_vals > 0
        n_missing_nets_posterior_draws_UNIF = Matrix(DataFrame(chain_UNIF)[:,(4+2*(n_net_types)+1):end])[sample_idxs_UNIF, :]
    else
        # No missing data, so just feed in a zero matrix with no columns
        n_missing_nets_posterior_draws_UNIF = zeros(size(chain_UNIF)[1], 0)
    end

    ##### Generate Net Crop Trajectories
    Γ_MONTHLY_samples_BYNET_UNIF = zeros(n_samples, length(MONTHS_MONTHLY), n_net_types)
    A_samples_BYNET_UNIF = zeros(n_samples, length(MONTHS_MONTHLY), length(MONTHS_MONTHLY), n_net_types)

    for i in 1:n_samples
        # Select MCMC posterior draw parameters
        ϕ = ϕ_posterior_draws_UNIF[i]
        α_init = α_init_posterior_draws_UNIF[i]
        α_LLIN = α_LLIN_posterior_draws_UNIF[i]
        b_nets = b_net_posterior_draws_UNIF[i,:]
        k_nets = k_net_posterior_draws_UNIF[i,:]
        missing_nets = n_missing_nets_posterior_draws_UNIF[i,:]

        Γ_MONTHLY_samples_BYNET_UNIF[i,:,:,:], A_samples_BYNET_UNIF[i,:,:,:] = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                                                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                                                                    ϕ, b_nets, k_nets,
                                                                                    α_init, α_LLIN,
                                                                                    missing_nets; 
                                                                                    monthly_p = monthly_p_UNIF,
                                                                                    return_age = true)
    end

    # Scale by to units of millions
    Γ_MONTHLY_samples_BYNET_UNIF = Γ_MONTHLY_samples_BYNET_UNIF./1e6
    A_samples_BYNET_UNIF = A_samples_BYNET_UNIF./1e6

    # Get Totals
    Γ_MONTHLY_samples_TOTAL_UNIF = sum(Γ_MONTHLY_samples_BYNET_UNIF, dims = 3)[:,:,1]
    A_samples_TOTAL_UNIF = sum(A_samples_BYNET_UNIF, dims = 4)[:,:,:,1]

    ### Part 3: Extract BV_model NET CROP and ACCESS draws
    Γ_MONTHLY_samples_TOTAL_BV = []
    ACCESS_samples_TOTAL_BV = []
    # if iso_in_bv
    #     Γ_MONTHLY_samples_TOTAL_BV = zeros(size(NPC_MONTHLY_BV_samples))
    #     ACCESS_samples_TOTAL_BV = zeros(size(ACCESS_MONTHLY_BV_samples))
    #     for i in 1:size(NPC_MONTHLY_BV_samples)[1]
    #         Γ_MONTHLY_samples_TOTAL_BV[i,:] = NPC_MONTHLY_BV_samples[i,:].*POPULATION_MONTHLY./1e6
    #         ACCESS_samples_TOTAL_BV[i,:] = ACCESS_MONTHLY_BV_samples[i,:]
    #     end
    # end

    # %%

    ##### Calculate Net Access
    λ_access_samples = sample_net_access(ρ_chain_df, μ_chain_df, p_h,
                                        POPULATION_MONTHLY, Γ_MONTHLY_samples_TOTAL.*1e6) # Need to temporarily unscale input by mil for calculating access

    ##### Calculate quantiles
    Γ_quantiles = zeros(3, length(MONTHS_MONTHLY))
    λ_access_quantiles = zeros(3, length(MONTHS_MONTHLY))

    for monthidx in 1:length(MONTHS_MONTHLY)
        Γ_quantiles[:,monthidx] = quantile(Γ_MONTHLY_samples_TOTAL[:,monthidx], [0.025,0.5,0.975])
        λ_access_quantiles[:,monthidx] = quantile(λ_access_samples[:,monthidx], [0.025,0.5,0.975])
    end

    ##### Get the true observed values for NET CROP MONTHLY and ACCESS for error plots later
    netcrop_monthidx = findall(.!ismissing.(NET_CROP_MONTHLY))

    netcrop_rmse = missings(Float64, length(netcrop_monthidx))
    netcrop_rmse_UNIF = missings(Float64, length(netcrop_monthidx))
    # netcrop_rmse_BV = missings(Float64, length(netcrop_monthidx))

    for i in 1:length(netcrop_monthidx)
        month_i = netcrop_monthidx[i]
        netcrop_rmse[i] = sqrt.(mean((Γ_MONTHLY_samples_TOTAL[:,month_i].-NET_CROP_MONTHLY[month_i]./1e6).^2))
        netcrop_rmse_UNIF[i] = sqrt.(mean((Γ_MONTHLY_samples_TOTAL_UNIF[:,month_i].-NET_CROP_MONTHLY[month_i]./1e6).^2))
        if iso_in_bv
            netcrop_rmse_BV[i] = sqrt.(mean((Γ_MONTHLY_samples_TOTAL_BV[:,month_i].-NET_CROP_MONTHLY[month_i]./1e6).^2))
        end
    end

    ##### Calculate age strata matrix for nets
    net_age_matrix_BYNET = zeros(n_samples, length(MONTHS_MONTHLY),max_age_quarters, n_net_types)
    for idx in 1:n_samples
        for i in 1:length(MONTHS_MONTHLY) #current time
            for j in 1:i # All possible birth times
                for n in 1:n_net_types
                    age_months = i-j
                    age_quarter = (i-j)÷4
                    age_index = min(age_quarter + 1, max_age_quarters)

                    net_age_matrix_BYNET[idx,i,age_index,n] += A_samples_BYNET[idx,i,j,n]
                end
            end
        end
    end


    ####### Calculate net age demographics

    # Calculate average across all samples and normalise to millions
    mean_net_age_matrix_BYNET = (mean(net_age_matrix_BYNET, dims = 1))[1,:,:,:]
    mean_net_age_matrix_TOTAL = sum(mean_net_age_matrix_BYNET, dims = 3)[:,:]

    # Calculate cumulative net age matrix
    cumulative_net_age_matrix_BYNET = zeros(size(mean_net_age_matrix_BYNET))

    for strata_i in 1:size(mean_net_age_matrix_BYNET)[2]
        for n in 1:n_net_types
            cumulative_net_age_matrix_BYNET[:,strata_i,n] = sum(mean_net_age_matrix_BYNET[:,end-(strata_i-1):end,n], dims = 2)
        end
    end
    cumulative_net_age_matrix_BYNET = cumulative_net_age_matrix_BYNET
    cumulative_net_age_matrix_TOTAL = sum(cumulative_net_age_matrix_BYNET, dims = 3)[:,:]

    # Calculate cumulative net age matrix of Total, but ordered by net types
    cumulative_net_age_matrix_STRATNET = zeros(size(mean_net_age_matrix_BYNET)[1],size(mean_net_age_matrix_BYNET)[2]*n_net_types)
    for n in 1:n_net_types
        for strata_i in 1:size(mean_net_age_matrix_BYNET)[2]
            index = (n-1)*max_age_quarters + strata_i
            if index == 1
                cumulative_net_age_matrix_STRATNET[:,index] = mean_net_age_matrix_BYNET[:,end-(strata_i-1),n]
            else
                cumulative_net_age_matrix_STRATNET[:,index] = cumulative_net_age_matrix_STRATNET[:,index-1] + mean_net_age_matrix_BYNET[:,end-(strata_i-1),n]
            end
        end
    end

    ##### Calculate net distribution ratios for plots
    monthly_p_scaled = zeros(size(monthly_p))
    for i in 1:(length(YEARS_ANNUAL))
        monthly_p_scaled[((i-1)*12+1):(i*12)] .= monthly_p[((i-1)*12+1):(i*12)]./sum(monthly_p[((i-1)*12+1):(i*12)])
    end

    ##### Make trajectory plots
    title = ISO
    if !isnothing(country_name)
        title = country_name
    end

    fig = Plots.plot(title = title, xlabel = "Year", ylabel = "Net Crop (mil)", margins = 2mm,
                ylims = (-0.05,maximum(Γ_MONTHLY_samples_TOTAL)*1.1))
    
    vline!(fig, MONTHS_MONTHLY[1:12:end], linecolor = 6, linestyle = :dash, alpha = 0.5,
            label = nothing,
            xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90)

    plot!(fig, MONTHS_MONTHLY, Γ_quantiles[1,:], fillrange = Γ_quantiles[3,:],
            color = 2, fillalpha = fillalpha, linealpha = 0, linewidth = lw, label = nothing)
    plot!(fig, MONTHS_MONTHLY, Γ_quantiles[2,:], linealpha = 1, linewidth = lw, color = 2,
            label = "Median Crop")
    scatter!(fig, MONTHS_MONTHLY, NET_CROP_MONTHLY./1e6,
                markersize = msize, markercolor = 2, label = "Survey Net Crop", legend = :topleft)

    plot!(fig, NaN.*MONTHS_MONTHLY, linealpha = 1, linewidth = lw,
                label = "Median Access", color = 4)    
    scatter!(fig, NaN.*MONTHS_MONTHLY, markersize = msize, markercolor = 4,
                label = "Survey Access")
    plot!(fig, NaN.*MONTHS_MONTHLY, linewidth = 1.5, label = "Dist. Proportion",
                linecolor = 1, linestyle = :dash)
    
    
    plot!(twinx(), monthly_p_scaled, linewidth = 1.5, ylabel = "Dist Proportion\n Net Access", 
                linecolor = 1, linestyle = :dash, linealpha = 0.8, legend = false,
                ylims = (-0.02, 1.05))
                
    plot!(twinx(), MONTHS_MONTHLY, λ_access_quantiles[2,:], linealpha = 1, linewidth = lw,
                label = "Median Access", legend = false, color = 4, yaxis = false,
                ylims = (-0.02, 1.05))
    plot!(twinx(), MONTHS_MONTHLY, λ_access_quantiles[1,:], fillrange = λ_access_quantiles[3,:],
            fillalpha = fillalpha, linealpha = 0, linewidth = lw, label = nothing,
            yaxis = false, color = 4,
            ylims = (-0.02, 1.05))
    scatter!(twinx(), MONTHS_MONTHLY, NET_ACCESS_SURVEY_MONTHLY,
                markersize = msize, markercolor = 4, label = nothing,
                yaxis = false, ylims = (-0.02, 1.05))

    ##### Make age strata plot (All nets combined)
    colors = palette(:roma, max_age_quarters)

    fig2_TOTAL = plot(title = "$title\nNet Age Stratification (Quarters)",
                xlabel = "Year", ylabel = "Net Crop (mil)", legend = :outerright,
                ylims = (-0.05,maximum(cumulative_net_age_matrix_TOTAL)*1.1))
    vline!(fig2_TOTAL, MONTHS_MONTHLY[1:12:end], linecolor = 6, linestyle = :dash, alpha = 0.5,
                label = nothing,
                xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90)
    for strata_i in 1:max_age_quarters
        if strata_i == 1
            plot!(fig2_TOTAL, zeros(size(mean_net_age_matrix_TOTAL)[1]),
                    fillrange = cumulative_net_age_matrix_TOTAL[:,strata_i],
                    linealpha = 0, fillalpha = 0.5, label = "$(max_age_quarters-(strata_i-1))",
                    color = colors[strata_i])
        else
            plot!(fig2_TOTAL, cumulative_net_age_matrix_TOTAL[:,strata_i-1],
                    fillrange = cumulative_net_age_matrix_TOTAL[:,strata_i],
                    linealpha = 0, fillalpha = 0.5, label = "$(max_age_quarters-(strata_i-1))",
                    color = colors[strata_i])
        end
    end

    ##### Make age strata plot separation_by_nets
    # Select base colors to represent each net types
    palette_offset = 3
    base_colors = palette(:Set1_9)[1:n_net_types]
    net_palettes = []
    for n in 1:n_net_types
        push!(net_palettes, palette([:black, base_colors[n],:white], max_age_quarters+palette_offset)[1:end-palette_offset])
    end
    
    fig2_BYNET = plot(title = "$title\nNet Type Age Stratification (Quarters)",
                xlabel = "Year", ylabel = "Net Crop (mil)", legend = :outerright,
                ylims = (-0.05,maximum(cumulative_net_age_matrix_TOTAL)*1.1))
    vline!(fig2_BYNET, MONTHS_MONTHLY[1:12:end], linecolor = 6, linestyle = :dash, alpha = 0.5,
                label = nothing,
                xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90)
    for n in 1:n_net_types
        for strata_i in 1:size(mean_net_age_matrix_BYNET)[2]
            if strata_i == 8
                label = NET_NAMES[n]
            else
                label = nothing
            end
            index = (n-1)*max_age_quarters + strata_i
            if index == 1
                plot!(fig2_BYNET, zeros(size(cumulative_net_age_matrix_STRATNET)[1]),
                        fillrange = cumulative_net_age_matrix_STRATNET[:,index],
                        linealpha = 0, fillalpha = 0.5, label = label,
                        color = net_palettes[n][strata_i])
            else
                plot!(fig2_BYNET, cumulative_net_age_matrix_STRATNET[:,index-1],
                        fillrange = cumulative_net_age_matrix_STRATNET[:,index],
                        linealpha = 0, fillalpha = 0.5, label = label,
                        color = net_palettes[n][strata_i])
            end
        end
    end

    ##### Make age strata plot only LLINs
    colors = net_palettes[2] # LLIN color palette
    
    fig2_LLIN = plot(title = "$title\nLLIN Age Stratification (Quarters)",
                xlabel = "Year", ylabel = "Net Crop (mil)", legend = :outerright,
                ylims = (-0.05,maximum(cumulative_net_age_matrix_TOTAL)*1.1))
    vline!(fig2_LLIN, MONTHS_MONTHLY[1:12:end], linecolor = 6, linestyle = :dash, alpha = 0.5,
                label = nothing,
                xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xtickfontrotation = 90)
    for strata_i in 1:max_age_quarters
        if strata_i == 1
            plot!(fig2_LLIN, zeros(size(mean_net_age_matrix_BYNET)[1]),
                    fillrange = cumulative_net_age_matrix_BYNET[:,strata_i,2],
                    linealpha = 0, fillalpha = 0.5, label = "$(max_age_quarters-(strata_i-1))",
                    color = colors[strata_i])
        else
            plot!(fig2_LLIN, cumulative_net_age_matrix_BYNET[:,strata_i-1,2],
                    fillrange = cumulative_net_age_matrix_BYNET[:,strata_i,2],
                    linealpha = 0, fillalpha = 0.5, label = "$(max_age_quarters-(strata_i-1))",
                    color = colors[strata_i])
        end
    end
    
            
    ##### Make RMSE Scatterplots (against UNIF)
    upper_lim = ceil(maximum(Γ_MONTHLY_samples_TOTAL[:,findall(.!ismissing.(NET_CROP_MONTHLY))]))
    netcrop_scatter_MAP = plot(xlabel = "Survey Estimates (mil)", ylabel = "Model Estimates (mil)",
                            title = "$(ISO) Net Crop MAP ITN\n RMSE = $(round(mean(netcrop_rmse), digits = 3))", legend = :bottomright,
                            xlims = (0,upper_lim), ylims = (0,upper_lim))
    netcrop_scatter_UNIF = plot(xlabel = "Survey Estimates (mil)", ylabel = "Model Estimates (mil)",
                            title = "$(ISO) Net Crop Uniform\n RMSE = $(round(mean(netcrop_rmse_UNIF), digits = 3))", legend = :bottomright,
                            xlims = (0,upper_lim), ylims = (0,upper_lim))
    plot!(netcrop_scatter_MAP,0:0.01:upper_lim,
            0:0.01:upper_lim, linealpha = 0.5,
            linestyle = :dash, color = :black, label = nothing)
    plot!(netcrop_scatter_UNIF,0:0.01:upper_lim,
            0:0.01:upper_lim,
            linestyle = :dash, color = :black, linealpha = 0.5,
            label = nothing)
    for i in 1:size(Γ_MONTHLY_samples_TOTAL)[1]
        if i == 1
            label1 = "MAP ITN Model Estimates"
            label2 = "Uniform Dist Model Estimates"
        else
            label1 = nothing
            label2 = nothing
        end
        scatter!(netcrop_scatter_MAP, NET_CROP_MONTHLY./1e6, Γ_MONTHLY_samples_TOTAL[i,:],
                markerstrokewidth = 0, markercolor = 3, markeralpha = 0.1, markersize = 2.5,
                label = label1)
        scatter!(netcrop_scatter_UNIF, NET_CROP_MONTHLY./1e6, Γ_MONTHLY_samples_TOTAL_UNIF[i,:],
                markerstrokewidth = 0, markercolor = 4, markeralpha = 0.1, markersize = 2.5,
                label = label2) 
    end
    scatter!(netcrop_scatter_MAP, NET_CROP_MONTHLY./1e6, Γ_quantiles[2,:],
                markerstrokewidth = 0.2, markercolor = 3, markersize = 8, markershape = :xcross,
                label = "MAP ITN Median") 
    scatter!(netcrop_scatter_UNIF, NET_CROP_MONTHLY./1e6, median(Γ_MONTHLY_samples_TOTAL_UNIF, dims = 1)[:],
                markerstrokewidth = 0.2, markercolor = 4, markersize = 8, markershape = :xcross,
                label = "Uniform Dist. Median")
    
    fig3 = plot(netcrop_scatter_MAP, netcrop_scatter_UNIF)

    ##### Make RMSE Scatterplots (against BV)
    fig3a = missing
    if iso_in_bv
        upper_lim = ceil(maximum(Γ_MONTHLY_samples_TOTAL[:,findall(.!ismissing.(NET_CROP_MONTHLY))]))
        netcrop_scatter_MAP = plot(xlabel = "Survey Estimates (mil)", ylabel = "Model Estimates (mil)",
                                title = "$(ISO) Net Crop MAP ITN\n RMSE = $(round(mean(netcrop_rmse), digits = 3))", legend = :bottomright,
                                xlims = (0,upper_lim), ylims = (0,upper_lim))
        netcrop_scatter_BV = plot(xlabel = "Survey Estimates (mil)", ylabel = "Model Estimates (mil)",
                                title = "$(ISO) Net Crop BV\n RMSE = $(round(mean(netcrop_rmse_BV), digits = 3))", legend = :bottomright,
                                xlims = (0,upper_lim), ylims = (0,upper_lim))
        plot!(netcrop_scatter_MAP,0:0.01:upper_lim,
                0:0.01:upper_lim, linealpha = 0.5,
                linestyle = :dash, color = :black, label = nothing)
        plot!(netcrop_scatter_BV,0:0.01:upper_lim,
                0:0.01:upper_lim,
                linestyle = :dash, color = :black, linealpha = 0.5,
                label = nothing)
        for i in 1:size(Γ_MONTHLY_samples_TOTAL)[1]
            if i == 1
                label1 = "MAP ITN Model Estimates"
                label2 = "BV Model Estimates"
            else
                label1 = nothing
                label2 = nothing
            end
            scatter!(netcrop_scatter_MAP, NET_CROP_MONTHLY./1e6, Γ_MONTHLY_samples_TOTAL[i,:],
                    markerstrokewidth = 0, markercolor = 3, markeralpha = 0.1, markersize = 2.5,
                    label = label1)
            scatter!(netcrop_scatter_BV, NET_CROP_MONTHLY./1e6, Γ_MONTHLY_samples_TOTAL_BV[i,:],
                    markerstrokewidth = 0, markercolor = 4, markeralpha = 0.1, markersize = 2.5,
                    label = label2) 
        end
        scatter!(netcrop_scatter_MAP, NET_CROP_MONTHLY, Γ_quantiles[2,:],
                    markerstrokewidth = 0.2, markercolor = 3, markersize = 8, markershape = :xcross,
                    label = "MAP ITN Median") 
        scatter!(netcrop_scatter_BV, NET_CROP_MONTHLY./1e6, median(Γ_MONTHLY_samples_TOTAL_BV, dims = 1)[:],
                    markerstrokewidth = 0.2, markercolor = 4, markersize = 8, markershape = :xcross,
                    label = "BV Median")
        
        fig3a = plot(netcrop_scatter_MAP, netcrop_scatter_BV)
    end

    return fig, fig2_TOTAL, fig2_BYNET, fig2_LLIN, fig3, fig3a, (NET_CROP_MONTHLY, Γ_MONTHLY_samples_TOTAL, Γ_MONTHLY_samples_TOTAL_UNIF, Γ_MONTHLY_samples_TOTAL_BV, NET_ACCESS_SURVEY_MONTHLY, λ_access_samples, ACCESS_samples_TOTAL_BV)
end

# %%
"""
    netlife_posterior_draws(regression_dict::Dict; t_vals = 0:0.01:5)
Calculate posterior draws of the survival curve based on the MCMC chain of the net crop model. Returns two outputs
- `lifecurve_quantiles`: 5%, 50% and 95% lifecurve profiles given in a (3, length(t_vals)) matrix.
- `halflife_quantiles`: 5%, 50% and 95% quantiles of the net half life.
"""
function netlife_posterior_draws(regression_dict; t_vals = 0:0.01:5)
    
    # Get number of net types
    NET_NAMES = regression_dict["NET_NAMES"]
    n_net_types = length(NET_NAMES)
    # n_net_types = (length(regression_dict["resampling_variances"])-3)÷2 #TEMPORARY CODE!!

    # Crop Regression output data
    chain = regression_dict["chain"]# Temporary
    n_samples = size(chain)[1]

    ##### Extract MCMC parameters
    # Extract MCMC draws for all parameters
    b_net_posterior_draws = Matrix(DataFrame(chain)[:,5:2:5+2*(n_net_types-1)])
    k_net_posterior_draws = Matrix(DataFrame(chain)[:,6:2:6+2*(n_net_types-1)])

    ##### Calculate lifecurves
    lifecurve_samples_BYNET = zeros(n_samples, length(t_vals), n_net_types)
    halflife_samples_BYNET = zeros(n_samples, n_net_types)

    for n in 1:n_net_types
        for i in 1:n_samples
            b = b_net_posterior_draws[i,n]
            k = k_net_posterior_draws[i,n]

            # Weibull Loss
            lifecurve_samples_BYNET[i,:,n] = net_loss_weibull.(t_vals, b, k)
            halflife_samples_BYNET[i,n] = net_life_weibull(0.5, b,k)

            # Compact Smooth Function
            lifecurve_samples_BYNET[i,:,n] = net_loss_compact.(t_vals, b, k)
            halflife_samples_BYNET[i,n] = net_life_compact(0.5, b,k)
        end
    end

    ##### Calculate quantiles
    lifecurve_quantiles_BYNET = zeros(3, length(t_vals), n_net_types)
    halflife_quantiles_BYNET = zeros(3, n_net_types)
    for n in 1:n_net_types
        for monthidx in 1:size(lifecurve_quantiles_BYNET)[2]
            lifecurve_quantiles_BYNET[:,monthidx,n] = quantile(lifecurve_samples_BYNET[:,monthidx,n], [0.025,0.5,0.975])
        end
        halflife_quantiles_BYNET[:,n] = quantile(halflife_samples_BYNET[:,n], [0.025,0.50,0.975])
    end

    return lifecurve_quantiles_BYNET, halflife_quantiles_BYNET
end

# %%
"""
    lifecurve_plots(ISO::String, lifecurve_quantiles::Matrix, halflife_quantiles::Vector; t_vals = 0:0.01:5, kwargs...)
Takes in the outputs of netlife_posterior draws and returns a plot of the attrition curve.
"""
function lifecurve_plots(ISO, NET_NAMES, lifecurve_quantiles_BYNET, halflife_quantiles_BYNET;
                        t_vals = 0:0.01:5,lw = 2, lw_minor = 1.5,
                        fillalpha = 0.2, country_name = nothing)

    pythonplot()
    theme(:vibrant)
    
    # Select base colors
    n_net_types = length(NET_NAMES)
    base_colors = palette(:Set1_9)[1:n_net_types]

    halflife_quantiles_BYNET
    # Make plot
    title = ISO
    if !isnothing(country_name)
        title = country_name
    end

    fig = plot(title = title, xlabel = "Years", ylabel = "Survival Rate",
                ylims = (-0.02, 1.05), xlims = (-0.1,5.1))
    hline!([0.5], linestyle = :dash, alpha = 0.8, linewidth = lw_minor, color = :black,
            label = nothing)
    for n in 1:n_net_types
        color = base_colors[n]
        vline!(halflife_quantiles_BYNET[:,n], linestyle = :dash,
                alpha = 0.8, linewidth = lw_minor, color = color, label = nothing)
        plot!(fig, t_vals, lifecurve_quantiles_BYNET[1,:,n], fillrange = lifecurve_quantiles_BYNET[3,:,n],
            xlabel = "Years", color = color, fillalpha = fillalpha,
            linealpha = 0, linewidth = lw, label = nothing)
        plot!(fig, t_vals, lifecurve_quantiles_BYNET[2,:,n], linealpha = 1, linewidth = lw,
            xlabel = "Years", label = NET_NAMES[n], color = color)
    end
    fig
 
    return fig

end

# %%
"""
    calc_scaling_region(input_dict, regression_dict,
                            net_access_input_dict, net_access_chain;
                            n_samples = 100)
Calculates linear scaling region given dicts of model inputs, and outputs. Returns a plot
of the data with a linear fit, and 95% confidence bounds of scaling, and width or scaling region
"""

function calc_scaling_region(input_dict, regression_dict,
                            net_access_input_dict, net_access_chain;
                            n_samples = 100)
    # n_samples = 100

    # ISO = "COG"
    # # Import Data
    # input_dict = load("outputs/extractions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
    # regression_dict = load("outputs/regressions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")
    # net_access_input_dict = load("outputs/extractions/access/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")
    # net_access_chain = load("outputs/regressions/access/netaccesschains.jld2")


    pythonplot()
    theme(:vibrant)

    # Crop Regression input data
    ISO = input_dict["ISO"]
    YEARS_ANNUAL = input_dict["YEARS_ANNUAL"]
    MONTHS_MONTHLY = input_dict["MONTHS_MONTHLY"]
    DELIVERIES_ANNUAL = input_dict["DELIVERIES_ANNUAL"]
    DISTRIBUTION_ANNUAL = input_dict["DISTRIBUTION_ANNUAL"]
    NET_CROP_MONTHLY = input_dict["NET_CROP_MONTHLY"]
    POPULATION_MONTHLY = input_dict["POPULATION_MONTHLY"]
    NET_NAMES = input_dict["NET_NAMES"]
    n_net_types = length(NET_NAMES)

    # Get number of missing net values
    n_missing_nets_vals = sum(ismissing.(DISTRIBUTION_ANNUAL[:,1]))

    # Crop Regression output data
    chain = regression_dict["chain"]# Temporary!
    monthly_p = regression_dict["monthly_p"]
    chain_UNIF = regression_dict["chain_epochs"][1]
    monthly_p_UNIF = regression_dict["monthly_weights"][1]

    # Access Regression output data
    ρ_chain_df = net_access_chain["ρ_chain_df"]
    μ_chain_df = net_access_chain["μ_chain_df"]
    p_h = mean(net_access_input_dict["p_h_aggregated"], dims = 1)[:]

    # Randomly generate sample indexes to sample from chain
    chain_length = size(chain)[1]
    sample_idxs = sample(1:chain_length, n_samples, replace = false)

    # Extract MCMC draws for parameters
    ϕ_posterior_draws = chain[sample_idxs, :ϕ]
    α_init_posterior_draws = chain[sample_idxs,:α_init]
    α_LLIN_posterior_draws = chain[sample_idxs,:α_LLIN]

    b_net_posterior_draws = Matrix(DataFrame(chain)[:,5:2:5+2*(n_net_types-1)])[sample_idxs, :]
    k_net_posterior_draws = Matrix(DataFrame(chain)[:,6:2:6+2*(n_net_types-1)])[sample_idxs, :]
    if n_missing_nets_vals > 0
        n_missing_nets_posterior_draws = Matrix(DataFrame(chain)[:,(4+2*(n_net_types)+1):end])[sample_idxs, :]
    else
        # No missing data, so just feed in a zero matrix with no columns
        n_missing_nets_posterior_draws = zeros(size(chain)[1], 0)
    end


    ##### Generate Net Crop Trajectories
    Γ_MONTHLY_samples_BYNET = zeros(n_samples, length(MONTHS_MONTHLY), n_net_types)
    A_samples_BYNET = zeros(n_samples, length(MONTHS_MONTHLY), length(MONTHS_MONTHLY), n_net_types)

    for i in 1:n_samples
        # Select MCMC posterior draw parameters
        ϕ = ϕ_posterior_draws[i]
        α_init = α_init_posterior_draws[i]
        α_LLIN = α_LLIN_posterior_draws[i]
        b_nets = b_net_posterior_draws[i,:]
        k_nets = k_net_posterior_draws[i,:]
        missing_nets = n_missing_nets_posterior_draws[i,:]

        Γ_MONTHLY_samples_BYNET[i,:,:,:], A_samples_BYNET[i,:,:,:] = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                                                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                                                                    ϕ, b_nets, k_nets,
                                                                                    α_init, α_LLIN,
                                                                                    missing_nets; 
                                                                                    monthly_p = monthly_p,
                                                                                    return_age = true)
    end

    # Get MONTHLY net crop totals
    Γ_MONTHLY_samples_BYNET = Γ_MONTHLY_samples_BYNET
    Γ_MONTHLY_samples_TOTAL = sum(Γ_MONTHLY_samples_BYNET, dims = 3)[:,:,1]

    λ_access_samples = sample_net_access(ρ_chain_df, μ_chain_df, p_h,
                                            POPULATION_MONTHLY, Γ_MONTHLY_samples_TOTAL) # Need to temporarily unscale input by mil for calculating access

    # %% Calculate mean Monthly Net crop
    NPC_MONTHLY_mean = mean(Γ_MONTHLY_samples_TOTAL, dims = 1)[:]./POPULATION_MONTHLY
    λ_ACCESS_mean = mean(λ_access_samples, dims = 1)[:]

    # %% Sort Net crop and Access by increasing Net Crop
    sortidx = sortperm(NPC_MONTHLY_mean)
    NPC_MONTHLY_mean = NPC_MONTHLY_mean[sortidx]
    λ_ACCESS_mean = λ_ACCESS_mean[sortidx]

    # %% Do Linear Interpolation to fill in gaps to get uniform spacing
    # Limit considered NPC scaling region to 2
    γ_vals = 0:0.005:(min(maximum(NPC_MONTHLY_mean), 0.4))
    λ_vals = zeros(length(γ_vals))

    for i in 1:length(γ_vals)
        γ = γ_vals[i]

        # Initialise values
        γ_left = 0
        λ_left = 0
        γ_right = NPC_MONTHLY_mean[end]
        λ_right = λ_ACCESS_mean[end]

        γ_idxs_left = findall((NPC_MONTHLY_mean.-γ).<0)
        if !isempty(γ_idxs_left)
            γ_left = NPC_MONTHLY_mean[γ_idxs_left[argmin(abs.(NPC_MONTHLY_mean[γ_idxs_left].-γ))]]
            λ_left = λ_ACCESS_mean[γ_idxs_left[argmin(abs.(NPC_MONTHLY_mean[γ_idxs_left].-γ))]]
        else
            γ_left = 0
            λ_left = 0
        end

        γ_idxs_right = findall((NPC_MONTHLY_mean.-γ).>0)

        if !isempty(γ_idxs_right)
            γ_right = NPC_MONTHLY_mean[γ_idxs_right[argmin(abs.(NPC_MONTHLY_mean[γ_idxs_right].-γ))]]
            λ_right = λ_ACCESS_mean[γ_idxs_right[argmin(abs.(NPC_MONTHLY_mean[γ_idxs_right].-γ))]]
        else
            γ_left = NPC_MONTHLY_mean[end-3]
            λ_left = λ_ACCESS_mean[end-3]
            γ_right = NPC_MONTHLY_mean[end]
            λ_right = λ_ACCESS_mean[end]
        end

        λ_vals[i] = λ_left + ((γ-γ_left)/(γ_right-γ_left))*(λ_right-λ_left)
    end


    # %% Calculate largest linear region
    # Find all candidate regions
    dxi = 10
    tol = 0.25
    lin_regions, slopes = linear_regions(γ_vals,λ_vals, dxi = dxi, tol = tol)
    llr_idx = argmax(length.(lin_regions)) # Index of the largest linear region
    llr_fit = LargestLinearRegion(dxi = dxi, tol = tol)
    
    m,m05,m95 = slopefit(Array(γ_vals[1:lin_regions[llr_idx][end]]), λ_vals[1:lin_regions[llr_idx][end]],
                    llr_fit, sat_threshold = 0.01)
    llr_idx_adjusted = argmin(abs.(slopes.-m))
    scaling_limit = γ_vals[lin_regions[llr_idx_adjusted][end]]
    # %% Plot
    y_offset = 0.05
    fig = scatter(NPC_MONTHLY_mean, λ_ACCESS_mean, markersize = 3, label = nothing, legend = false,
            labelfontsize = 13, title = "$(ISO)\nLSR = ($(round(m05, digits = 3)), $(round(m95, digits = 3)))",
            xlabel = "Nets Per Capita "*L"(\gamma)", ylabel = "Access "*L"(\lambda)",
            xlims = (-0.01, 1.05),
            ylims = (-0.01, 1.05))
    # scatter(γ_vals, λ_vals, markersize = 3, label = nothing, legend = false,
    #         labelfontsize = 13, title = "$(ISO)\nLSR = ($(round(m05, digits = 3)), $(round(m95, digits = 3)))",
    #         xlabel = "Nets Per Capita "*L"(\gamma)", ylabel = "Access "*L"(\lambda)",
    #         xlims = (-0.01, round(maximum(NPC_MONTHLY_mean), digits =1)),
    #         ylims = (-0.01, round(maximum(λ_ACCESS_mean), digits =1)+0.05))
    plot!(fig, 0:0.01:1, m.*(0:0.01:1).+ y_offset, color = 2, width = 2)
    plot!(fig, 0:0.01:1, m05.*(0:0.01:1).+ y_offset, color = 3, width = 0.8, linestyle = :dash)
    plot!(fig, 0:0.01:1, m95.*(0:0.01:1).+ y_offset, color = 3, width = 0.8, linestyle = :dash)

    return fig, (m, m05, m95), scaling_limit
end

end
