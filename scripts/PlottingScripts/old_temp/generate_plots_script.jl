# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")


# %% Import Public Packages
using JLD2
using CSV
using DataFrames
using ProgressBars
using Measures
using Plots
using StatsBase
pythonplot()

# %% Import helper packages
using PlottingFunctions
using LaTeXStrings
using DynamicalSystems

# %%
output_dir = "output_plots/"

# %% Get ISO List
ISO_list = String.(CSV.read("/datasets/ISO_list.csv", DataFrame)[:,1])
reg_results_list = Array{Any, 1}(undef, length(ISO_list))
exclusion_ISOs = ["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD"]
#["CPV","BWA","CAF","COM","GNQ","DJI","ERI","ETH","GAB","GNB","STP","SOM","SDN","SWZ","ZAF","SSD"]

YEAR_START = 2000
YEAR_END = 2023

country_codes_key = CSV.read("datasets/country_codes.csv", DataFrame)
# %%

# ETH Fits seem to fail for  net access
# COM Seems like NPC might be too low
# %% Storage list for figures
fig_timeseries_collection = []
fig_agestrata_TOTAL_collection = []
fig_agestrata_BYNET_collection = []
fig_rmsescatter_collection = []
fig_rmsescatter_BV_collection = []
fig_lifecurves_collection = []

# %% Storage list for halflife quantiles
halflife_quantiles_collection = []
RMSE_results_collection = []

for i in ProgressBar(1:length(ISO_list))
    # Select ISO
    ISO = ISO_list[i]

    if ISO ∈ exclusion_ISOs
        continue
    else
        println(ISO)

        # Import Data
        input_dict = load("outputs/extractions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
        # regression_dict = load("outputs/regressions/crop/Compact Regressions/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")
        regression_dict = load("outputs/regressions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")
        net_access_input_dict = load("outputs/extractions/access/pred_data/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")
        net_access_chain = load("outputs/regressions/access/netaccesschains.jld2")
        BV_dict = load("datasets/BV_draws.jld2")
        country_name = country_codes_key[findfirst(country_codes_key.ISO3 .== ISO), "Country"]

        # Calculate net life values
        lifecurve_quantiles, halflife_quantiles = netlife_posterior_draws(regression_dict)

        # Generate Pre-baked Plots
        fig_timeseries, fig_agestrata_TOTAL, fig_agestrata_BYNET, fig_agestrata_LLIN, fig_rmsescatter, fig_rmsescatter_BV, results_tuple = regression_timeseries_plot(input_dict, regression_dict,
                                                                    net_access_input_dict, net_access_chain;# BV_dict; 
                                                                    country_name = country_name)
        
        NET_NAMES = input_dict["NET_NAMES"]
        
        fig_lifecurves = lifecurve_plots(ISO,NET_NAMES,lifecurve_quantiles, halflife_quantiles; 
                                            country_name = country_name)

        push!(fig_timeseries_collection, fig_timeseries)
        push!(fig_agestrata_TOTAL_collection, fig_agestrata_TOTAL)
        push!(fig_agestrata_BYNET_collection, fig_agestrata_BYNET)
        push!(fig_rmsescatter_collection, fig_rmsescatter)
        push!(fig_rmsescatter_BV_collection, fig_rmsescatter_BV)
        push!(fig_lifecurves_collection, fig_lifecurves)
        push!(halflife_quantiles_collection, halflife_quantiles)
        push!(RMSE_results_collection, results_tuple)

        # Save individual component plots in output directory
        fig1 = plot(fig_timeseries_collection[end], size = (750,500),
        legendfontsize = 10, labelfontsize = 17, titlefontsize = 20,
        tickfontsize = 11)
        savefig(fig1, output_dir*"individual/trajectories/$(ISO)_trajectory.pdf")

        fig2 = plot(fig_agestrata_TOTAL_collection[end], size = (750,500),
                legendfontsize = 10, labelfontsize = 17, titlefontsize = 20,
                tickfontsize = 11)
        savefig(fig2, output_dir*"individual/net_age/$(ISO)_totalstrata.pdf")

        fig2a = plot(fig_agestrata_BYNET_collection[end], size = (750,500),
                legendfontsize = 10, labelfontsize = 17, titlefontsize = 20,
                tickfontsize = 11)
        savefig(fig2a, output_dir*"individual/net_age/$(ISO)_netstrata.pdf")

        fig3 = plot(fig_rmsescatter_collection[end], size = (750,500),
                legendfontsize = 10, labelfontsize = 13, titlefontsize = 16,
                tickfontsize = 10, legend = :topleft)
        savefig(fig3, output_dir*"individual/rmse/$(ISO)_rmse.pdf")

        fig4 = plot(fig_lifecurves_collection[end], size = (750,500),
                legendfontsize = 12, labelfontsize = 17, titlefontsize = 20,
                tickfontsize = 11)
        savefig(fig4, output_dir*"individual/attrition/$(ISO)_lifecurves.pdf")
    end
end


# %% Net crop and access trajectory plots
fig1 = plot(fig_timeseries_collection..., layout = (5,7), size = (3840,2160), margin=15*Plots.mm)
# fig1 = plot(fig_timeseries_collection..., layout = (3,4), size = (3840,2160), margin=15*Plots.mm)
savefig(fig1, output_dir*"Crop_Access_All.pdf")

# %% Age strata plots
fig2 = plot(fig_agestrata_TOTAL_collection..., layout = (5,7), size = (3840,2160), margin=10*Plots.mm)
# fig2 = plot(fig_agestrata_TOTAL_collection..., layout = (3,4), size = (3840,2160), margin=10*Plots.mm)
savefig(fig2, output_dir*"Net_AgeStrata_TOTAL_All.pdf")

# %% Age strata plots
fig2a = plot(fig_agestrata_BYNET_collection..., layout = (5,7), size = (3840,2160), margin=10*Plots.mm)
# fig2a = plot(fig_agestrata_BYNET_collection..., layout = (3,4), size = (3840,2160), margin=10*Plots.mm)
savefig(fig2a, output_dir*"Net_AgeStrata_BYNET_All.pdf")

# %% RMSE plots
fig3 = plot(fig_rmsescatter_collection..., layout = (5,7), size = (3840,2160), margin=10*Plots.mm)
# fig3 = plot(fig_rmsescatter_collection..., layout = (3,4), size = (3840,2160), margin=10*Plots.mm)
savefig(fig3, output_dir*"Net_RMSE_All.pdf")

# %% Lifecurve plots
fig4 = plot(fig_lifecurves_collection..., layout = (5,7), size = (3840,2160), margin=10*Plots.mm)
# fig4 = plot(fig_lifecurves_collection..., layout = (3,4), size = (3840,2160), margin=10*Plots.mm)
savefig(fig4, output_dir*"Lifecurves_All.pdf")


# %% Calculate RMSE per Country for MAP ITN vs Uniform Dist
n_samples = size(RMSE_results_collection[1][2])[1] # Get number of samples

NETCROP_survey_reference_data = zeros(0)
NETCROP_MAP_ITN_model_pred = zeros(n_samples,0)
NETCROP_UNIF_model_pred = zeros(n_samples,0)
NETCROP_BV_model_pred = zeros(n_samples,0)
ACCESS_survey_reference_data = zeros(0)
ACCESS_MAP_ITN_model_pred = zeros(n_samples,0)
ACCESS_BV_model_pred = zeros(n_samples,0)

for i in 1:length(RMSE_results_collection)
    println(i)
    if i ∈ [20,24,27] # Cases where there is no BV results
        continue
    end
    # Extract prediction trajectories
    NET_CROP_MONTHLY, Γ_MONTHLY_samples_TOTAL, Γ_MONTHLY_samples_TOTAL_UNIF, Γ_MONTHLY_samples_TOTAL_BV, NET_ACCESS_SURVEY_MONTHLY, λ_access_samples, ACCESS_samples_TOTAL_BV = RMSE_results_collection[i]
    
    # Find index of all non-missing data
    NETCROP_nonmissing_idx = findall(.!ismissing.(NET_CROP_MONTHLY))
    ACCESS_nonmissing_idx = findall(.!ismissing.(NET_ACCESS_SURVEY_MONTHLY))

    # Append to aggregate
    # Net Crop
    NETCROP_survey_reference_data = vcat(NETCROP_survey_reference_data, NET_CROP_MONTHLY[NETCROP_nonmissing_idx]./1e6)
    NETCROP_MAP_ITN_model_pred = hcat(NETCROP_MAP_ITN_model_pred, Γ_MONTHLY_samples_TOTAL[:, NETCROP_nonmissing_idx])
    NETCROP_BV_model_pred = hcat(NETCROP_BV_model_pred, Γ_MONTHLY_samples_TOTAL_BV[1:n_samples, NETCROP_nonmissing_idx])
    NETCROP_UNIF_model_pred = hcat(NETCROP_UNIF_model_pred, Γ_MONTHLY_samples_TOTAL_UNIF[:, NETCROP_nonmissing_idx])

    # Net Access
    ACCESS_survey_reference_data = vcat(ACCESS_survey_reference_data, NET_ACCESS_SURVEY_MONTHLY[ACCESS_nonmissing_idx])
    ACCESS_MAP_ITN_model_pred = hcat(ACCESS_MAP_ITN_model_pred, λ_access_samples[:, ACCESS_nonmissing_idx])
    ACCESS_BV_model_pred = hcat(ACCESS_BV_model_pred, ACCESS_samples_TOTAL_BV[1:n_samples, ACCESS_nonmissing_idx])
end


# %% NETCROP RMSE Plots
NETCROP_rmse = sqrt.(mean((NETCROP_MAP_ITN_model_pred .- repeat(NETCROP_survey_reference_data, 1,n_samples)').^2, dims = 1)[:])
NETCROP_rmse_UNIF = sqrt.(mean((NETCROP_UNIF_model_pred .- repeat(NETCROP_survey_reference_data, 1,n_samples)').^2, dims = 1)[:])
NETCROP_rmse_BV = sqrt.(mean((NETCROP_BV_model_pred .- repeat(NETCROP_survey_reference_data, 1,n_samples)').^2, dims = 1)[:])


NETCROP_scatter_MAP = plot(xlabel = "Log Survey Estimates (mil)", ylabel = "Log Model Estimates (log mil)",
                            title = "All Countries\nNet Crop MAP ITN\n RMSE = $(round(mean(NETCROP_rmse), digits = 3))", 
                            legend = :topleft, xlims = log.((1, (ceil(maximum(NETCROP_survey_reference_data))+8))),
                            ylims = log.((1, (ceil(maximum(NETCROP_survey_reference_data))+5))))
NETCROP_scatter_UNIF = plot(xlabel = "Log Survey Estimates (mil)", ylabel = "Log Model Estimates (log mil)",
                            title = "All Countries\nNet Crop MAP UNIF\n RMSE = $(round(mean(NETCROP_rmse_UNIF), digits = 3))", 
                            legend = :topleft, xlims = log.((1, (ceil(maximum(NETCROP_survey_reference_data))+8))),
                            ylims = log.((1, (ceil(maximum(NETCROP_survey_reference_data))+5))))
NETCROP_scatter_BV = plot(xlabel = "Log Survey Estimates (mil)", ylabel = "Log Model Estimates (log mil)",
                            title = "All Countries\nNet Crop BV\n RMSE = $(round(mean(NETCROP_rmse_BV), digits = 3))", 
                            legend = :topleft, xlims = log.((1, (ceil(maximum(NETCROP_survey_reference_data))+8))),
                            ylims = log.((1, (ceil(maximum(NETCROP_survey_reference_data))+5))))


plot!(NETCROP_scatter_MAP, 0:0.1:(ceil(maximum(NETCROP_survey_reference_data))+5), 0:0.1:(ceil(maximum(NETCROP_survey_reference_data))+5),
    linealpha = 0.5, linestyle = :dash, color = :black, label = nothing)
plot!(NETCROP_scatter_UNIF, 0:0.1:(ceil(maximum(NETCROP_survey_reference_data))+5), 0:0.1:(ceil(maximum(NETCROP_survey_reference_data))+5),
    linealpha = 0.5, linestyle = :dash, color = :black, label = nothing)
plot!(NETCROP_scatter_BV, 0:0.1:(ceil(maximum(NETCROP_survey_reference_data))+5), 0:0.1:(ceil(maximum(NETCROP_survey_reference_data))+5),
    linealpha = 0.5, linestyle = :dash, color = :black, label = nothing)


for i in 1:n_samples
    if i == 1
        label1 = "MAP ITN Model Estimates"
        label2 = "MAP Uniform Dist Estimates"
        label3 = "BV Estimates"
    else
        label1 = nothing
        label2 = nothing
        label3 = nothing
    end

    scatter!(NETCROP_scatter_MAP, log.(NETCROP_survey_reference_data), log.(NETCROP_MAP_ITN_model_pred[i,:]),
                markerstrokewidth = 0, markercolor = 3, markeralpha = 0.1, markersize = 2.5,
                label = label1)
    scatter!(NETCROP_scatter_UNIF, log.(NETCROP_survey_reference_data), log.(NETCROP_UNIF_model_pred[i,:]),
            markerstrokewidth = 0, markercolor = 6, markeralpha = 0.1, markersize = 2.5,
            label = label2) 
    scatter!(NETCROP_scatter_BV, log.(NETCROP_survey_reference_data), log.(NETCROP_BV_model_pred[i,:]),
            markerstrokewidth = 0, markercolor = 4, markeralpha = 0.1, markersize = 2.5,
            label = label3) 
end

scatter!(NETCROP_scatter_MAP, log.(NETCROP_survey_reference_data), log.(median(NETCROP_MAP_ITN_model_pred, dims = 1)[:]),
                markerstrokewidth = 0.2, markercolor = 3, markersize = 8, markershape = :xcross,
                label = "MAP ITN Median") 
scatter!(NETCROP_scatter_UNIF, log.(NETCROP_survey_reference_data), log.(median(NETCROP_UNIF_model_pred, dims = 1)[:]),
            markerstrokewidth = 0.2, markercolor = 6, markersize = 8, markershape = :xcross,
            label = "Uniform Median")
scatter!(NETCROP_scatter_BV, log.(NETCROP_survey_reference_data), log.(median(NETCROP_BV_model_pred, dims = 1)[:]),
            markerstrokewidth = 0.2, markercolor = 4, markersize = 8, markershape = :xcross,
            label = "BV Median")

aggregated_RMSE_fig = plot(NETCROP_scatter_MAP, NETCROP_scatter_UNIF, NETCROP_scatter_BV, layout = (1,3), size = (1000,400))

savefig(aggregated_RMSE_fig, "NETCROP_RMSE_Aggregated.pdf")

# %% NETCROP RMSE Plots
ACCESS_rmse = sqrt.(mean((ACCESS_MAP_ITN_model_pred .- repeat(ACCESS_survey_reference_data, 1,n_samples)').^2, dims = 1)[:])
ACCESS_rmse_BV = sqrt.(mean((ACCESS_BV_model_pred .- repeat(ACCESS_survey_reference_data, 1,n_samples)').^2, dims = 1)[:])


ACCESS_scatter_MAP = plot(xlabel = "Survey Estimates", ylabel = "Model Estimates",
                            title = "All Countries\nNet Access MAP ITN\n RMSE = $(round(mean(ACCESS_rmse), digits = 3))", 
                            legend = :topleft, xlims = (-0.05, 1.05),
                            ylims = (-0.05, 1.05))
ACCESS_scatter_BV = plot(xlabel = "Survey Estimates", ylabel = "Model Estimates",
                            title = "All Countries\nNet Access BV\n RMSE = $(round(mean(ACCESS_rmse_BV), digits = 3))", 
                            legend = :topleft, xlims = (-0.05, 1.05),
                            ylims = (-0.05, 1.05))


plot!(ACCESS_scatter_MAP, 0:0.1:(ceil(maximum(ACCESS_survey_reference_data))+5), 0:0.1:(ceil(maximum(ACCESS_survey_reference_data))+5),
    linealpha = 0.5, linestyle = :dash, color = :black, label = nothing)
plot!(ACCESS_scatter_BV, 0:0.1:(ceil(maximum(ACCESS_survey_reference_data))+5), 0:0.1:(ceil(maximum(ACCESS_survey_reference_data))+5),
    linealpha = 0.5, linestyle = :dash, color = :black, label = nothing)


for i in 1:n_samples
    if i == 1
        label1 = "MAP ITN Model Estimates"
        label3 = "BV Estimates"
    else
        label1 = nothing
        label2 = nothing
        label3 = nothing
    end

    scatter!(ACCESS_scatter_MAP, ACCESS_survey_reference_data, ACCESS_MAP_ITN_model_pred[i,:],
                markerstrokewidth = 0, markercolor = 3, markeralpha = 0.1, markersize = 2.5,
                label = label1)
    scatter!(ACCESS_scatter_BV, ACCESS_survey_reference_data, ACCESS_BV_model_pred[i,:],
            markerstrokewidth = 0, markercolor = 4, markeralpha = 0.1, markersize = 2.5,
            label = label3) 
end

scatter!(ACCESS_scatter_MAP, ACCESS_survey_reference_data, median(ACCESS_MAP_ITN_model_pred, dims = 1)[:],
                markerstrokewidth = 0.2, markercolor = 3, markersize = 8, markershape = :xcross,
                label = "MAP ITN Median")
scatter!(ACCESS_scatter_BV, ACCESS_survey_reference_data, median(ACCESS_BV_model_pred, dims = 1)[:],
            markerstrokewidth = 0.2, markercolor = 4, markersize = 8, markershape = :xcross,
            label = "BV Median")

aggregated_RMSE_fig = plot(ACCESS_scatter_MAP, ACCESS_scatter_BV, layout = (1,2), size = (700,400))

savefig(aggregated_RMSE_fig, "ACCESS_RMSE_Aggregated.pdf")


# %% Make Halflife plots
a = Float64[]
b = Float64[]
c = Float64[]
for i in 1:length(halflife_quantiles_collection)
    a = vcat(a,halflife_quantiles_collection[i][1,2])
    b = vcat(b,halflife_quantiles_collection[i][2,2])
    c = vcat(c,halflife_quantiles_collection[i][3,2])
end

halflife_quantiles_dataframe = DataFrame(ISO3 = setdiff(ISO_list, exclusion_ISOs),
                            a = a,
                            b = b,
                            c = c)
rename!(halflife_quantiles_dataframe, Dict(:a=> "NEW5%", :b=> "NEW50%", :c=> "NEW95%"))

# Import data from Bertozzi-Villa Results
LLIN_retention_ref = CSV.read(raw"itn-julia\database\bertozzi_villa_LLIN_retention_times.csv", DataFrame)

# Fill in missing values
for i in 1:size(LLIN_retention_ref)[1]
    if ismissing(LLIN_retention_ref[i,2])
        LLIN_retention_ref[i,2:end] .= missing
    end
end

rename!(LLIN_retention_ref, Dict("5%" => "BV5%", "50%" => "BV50%", "95%" => "BV95%"))

# %% Find list of countries that are shared between both datasets
ISO_list_intersect = intersect(halflife_quantiles_dataframe.ISO3,LLIN_retention_ref.ISO3)

function check_membership(list_A, list_B)
    membership_bool = zeros(Bool, length(list_A))
    for i in 1:length(list_A)
        if list_A[i] ∈ list_B
            membership_bool[i] = true
        else
            membership_bool[i] = false
        end
    end
    return membership_bool
end

# %% Merge quantiles and sort
filt_halflife_quantiles = halflife_quantiles_dataframe[check_membership(halflife_quantiles_dataframe.ISO3, ISO_list_intersect),:]
filt_halflife_quantiles_BV = LLIN_retention_ref[check_membership(LLIN_retention_ref.ISO3, ISO_list_intersect),:]


for name in names(filt_halflife_quantiles_BV)[2:end]
    filt_halflife_quantiles[!,name] = filt_halflife_quantiles_BV[:,name]
end

filt_halflife_quantiles = filt_halflife_quantiles[sortperm(filt_halflife_quantiles[:,"NEW50%"]),:]

# %%
gr()
lw = 4
la = 0.7
ms = 6
sep = 0.15
col1 = :steelblue;
col2 = :firebrick;
fig = plot(xticks = (1:size(filt_halflife_quantiles)[1], filt_halflife_quantiles[:,1]),
        minorticks = false, tickfontsize = 10, legendfontsize = 10, labelfontsize = 15, 
        xrotation = 90, size = (1200,400), framestyle = :box, legend  = :topleft,
        xlabel = "Country", ylabel = "LLIN Half Life (years)", gridalpha = 0.05, margins = 9mm,
        ylims = (-0.1, 7.05))

for i in 1:size(filt_halflife_quantiles)[1]
    if i == 1
        label1 = "MAP ITN"
        label2 = "BV Model"
    else
        label1 = nothing
        label2 = nothing
    end
    plot!(fig, [i,i].-sep, Vector(filt_halflife_quantiles[i,["NEW5%","NEW95%"]]),
            color = col1, linewidth = lw, label = label1, alpha = la)
    scatter!(fig, [i].-sep, [filt_halflife_quantiles[i,"NEW50%"]],
            markercolor = col1, markersize = ms, markerstrokewidth = 0, label = nothing)
    
    plot!(fig, [i,i].+sep, Vector(filt_halflife_quantiles[i,["BV5%","BV95%"]]),
            color = col2, linewidth = lw, label = label2, alpha = la)
    scatter!(fig, [i].+sep, [filt_halflife_quantiles[i,"BV50%"]],
            markercolor = col2, markersize = ms, markerstrokewidth = 0, label = nothing)

    # ISO = LLIN_retention_red[i,1]
    # ref_idx = findall(LLIN_retention_ref.ISO3.==ISO)[1]
    # plot!(fig, [i,i].+sep, Vector(LLIN_retention_ref[ref_idx,[2,4]]),
    #         color = col3, linewidth = lw, label = label3, alpha = la)
    # scatter!(fig, [i].+sep, [LLIN_retention_ref[ref_idx,3]],
    #         markercolor = col3, markersize = ms, markerstrokewidth = 0, label = nothing)
    
end
savefig(fig, "MAP_SNF_Halflife.pdf")
fig

# # %% Make Linear Scaling Region Slope Plots
# # Different set of exclusion ISOs, due to difficulty in estimating scaling region
# exclusion_ISOs = ["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD","STP"]
# ISOs = []
# scaling_plots_collection = []
# scaling_limit_values = []
# slope_values = zeros(3,0)

# for i in ProgressBar(1:length(ISO_list))
#     # Select ISO
#     ISO = ISO_list[i]

#     if ISO ∈ exclusion_ISOs
#         continue
#     else
#         println(ISO)

#         # Import Data
#         input_dict = load("outputs/extractions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
#         regression_dict = load("outputs/regressions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")
#         net_access_input_dict = load("outputs/extractions/access/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")
#         net_access_chain = load("outputs/regressions/access/netaccesschains.jld2")


#         fig, (m, m05, m95), scaling_limit = calc_scaling_region(input_dict, regression_dict,
#                                     net_access_input_dict, net_access_chain)

#         push!(ISOs, ISO)
#         push!(scaling_plots_collection, fig)
#         push!(scaling_limit_values, scaling_limit)
#         slope_values = hcat(slope_values, [m, m05,m95])
#     end
# end

# ISO = "BEN"
# regression_dict = load("outputs/regressions/crop/Compact Regressions/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")

# regression_dict["chain"]
# # %%
# y_offset = 0.01
# x_offset = 0.01
# fig = plot(xlabel = "Scaling Region", ylabel = "Slope", legend = false, ylims = (0.95, 3.05))
# for i in 1:length(ISOs)
#     plot!(fig, scaling_limit_values[i].*ones(2), slope_values[2:3,i], width = 2.5, linealpha = 0.3,
#             color = 1)
# end
# scatter!(scaling_limit_values, slope_values[1,:], color = 3)
# [annotate!(x+x_offset*rand([1,-1]), y+y_offset*rand([1,-1]), 
#             Plots.text(ISOs[i], :black, 8)) for (i,x,y) in zip(1:length(ISOs),scaling_limit_values,slope_values[1,:])];
# fig
# savefig(fig, "Scaling_Regions_Summary.pdf")

# # %%
# fig = plot(scaling_plots_collection..., layout = (5,7), size = (3840,2160), margin=10*Plots.mm)
# savefig(fig, "Scaling_Regions_All.pdf")
