# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import Native Packages
using ProgressBars
using JLD2
using CSV
using DataFrames

# %% Import National MITN Plotting functions
using NationalMITN_Plots

# %% Define save paths
output_path = "output_plots/snf/national/summaries/"
mkpath(output_path)

##############################################
# %% Get list of countries to plot
##############################################

ISO_list = String.(CSV.read("datasets/ISO_list.csv", DataFrame)[:,1])
exclusion_ISOs = ["CPV","ZAF"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

##############################################
# %% Make storage variables for plots
##############################################
n_countries = length(filt_ISOs)

nat_timeseries_collection = Array{Any}(undef, n_countries)
netcrop_demography_collection = Array{Any}(undef, n_countries)
npc_demography_collection = Array{Any}(undef, n_countries)
netcrop_bytype_collection = Array{Any}(undef, n_countries)
npc_bytype_collection = Array{Any}(undef, n_countries)
attrition_curves_collection = Array{Any}(undef, n_countries)

##############################################
# %% Generate Plots for each country on list
##############################################
# %% Analysed Time Bounds
YEAR_START = 2000
YEAR_END = 2023

for ISO_i in ProgressBar(1:length(filt_ISOs))

    # Select ISO
    ISO = filt_ISOs[ISO_i]

    # Load Datasets
    input_dict = load("outputs/extractions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
    regression_dict = load("outputs/regressions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")
    net_access_input_dict = load("outputs/extractions/access/reg_data/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")
    net_access_chain = load("outputs/regressions/access/netaccesschains.jld2")

    # Crop and Access draws
    post_snf = load("outputs/draws/national/crop_access/$(ISO)_$(YEAR_START)_$(YEAR_END)_post_crop_access.jld2")
    post_net_attrition = CSV.read("outputs/net_attrition_posteriors.csv", DataFrame)
    post_net_demography_mean = CSV.read("outputs/draws/national/demography/$(ISO)_net_age_demography_mean.csv", DataFrame)
    post_net_demography_samples = CSV.read("outputs/draws/national/demography/$(ISO)_net_age_demography_samples.csv", DataFrame)

    # Add to plot collections
    nat_timeseries_collection[ISO_i] = plot_nat_timeseries(input_dict, net_access_input_dict, post_snf)
    netcrop_demography_collection[ISO_i] = plot_netcrop_demography(input_dict, post_net_demography_mean)
    npc_demography_collection[ISO_i] = plot_npc_demography(input_dict, post_net_demography_mean)
    netcrop_bytype_collection[ISO_i] = plot_netcrop_bytype(input_dict, post_snf)
    npc_bytype_collection[ISO_i] = plot_npc_bytype(input_dict, post_snf)
    attrition_curves_collection[ISO_i] = plot_attrition_curves(regression_dict)

end

##############################################
# %% Construct country level aggregate plots
##############################################

for ISO_i in ProgressBar(1:length(filt_ISOs))
    ISO = filt_ISOs[ISO_i]

    country_fig = plot(   nat_timeseries_collection[ISO_i],
                            netcrop_bytype_collection[ISO_i],
                            npc_bytype_collection[ISO_i],
                            attrition_curves_collection[ISO_i],
                            netcrop_demography_collection[ISO_i],
                            npc_demography_collection[ISO_i],
                            layout = (2,3), size = (1600,700))

    savefig(country_fig, output_path*"country/$(ISO)_summary.pdf")
end

##############################################
# %% Construct country level aggregate plots
##############################################

# Plot layout settings
layout = (8,6)
figsize = (2560,1800)

# Make plots
nat_timeseries_fig = plot(nat_timeseries_collection..., layout = layout, size = figsize)
netcrop_demography_fig = plot(netcrop_demography_collection..., layout = layout, size = figsize, legendfontsize = 4.5, legend = :outerright)
npc_demography_fig = plot(npc_demography_collection..., layout = layout, size = figsize, legendfontsize = 4.5, legend = :outerright)
netcrop_bytype_fig = plot(netcrop_bytype_collection..., layout = layout, size = figsize)
npc_bytype_fig = plot(npc_bytype_collection..., layout = layout, size = figsize)
attrition_curves_fig = plot(attrition_curves_collection..., layout = layout, size = figsize)

# Save plots
savefig(nat_timeseries_fig, output_path*"snf_nat_timeseries.pdf")
savefig(netcrop_demography_fig, output_path*"snf_nat_crop_demography.pdf")
savefig(npc_demography_fig, output_path*"snf_nat_npc_demography.pdf")
savefig(netcrop_bytype_fig, output_path*"snf_nat_crop_bytype.pdf")
savefig(npc_bytype_fig, output_path*"snf_nat_npc_bytype.pdf")
savefig(attrition_curves_fig, output_path*"snf_nat_attrition_curves.pdf")

