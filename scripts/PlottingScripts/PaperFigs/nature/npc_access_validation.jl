"""
Author: Eugene Tan
Date Created: 12/5/2025
Last Updated: 12/5/2025
Make validation plots to verify NPC-Access saturation relationship
"""
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Import Public Packages
using DataFrames
using JLD2
using CSV
using ProgressBars

# %% Math packages
using StatsBase

# %% Plotting Packages
using CairoMakie

# %% Custom packages
using CropAccessFromSurvey
using NetAccessModel
using NetAccessPrediction


# %% Plotting Theme and general settings
set_theme!(theme_ggplot2())

# Color settings
colors = [  colorant"#0082C7",
            colorant"#E72A3D",
            colorant"#538255",
            colorant"#45332C",
            ];

# %%
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

# %% Define list of countries to extract data from and plot
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)
###################################################
# %% POSTPROCESSING: Calculate Raster Estimates of Net Crop disaggrated by type (According to SNF ratios)
###################################################

# %% Extract numerical/empirical estiamates of NPC-Access values from household survey data
NPC_SURVEY_aggregated = []
ACCESS_SURVEY_aggregated = []

for ISO_i in 1:length(filt_ISOs)
    # Select ISO
    ISO = filt_ISOs[ISO_i]

    # Import required data
    input_dict = load(OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
    regression_dict = load(OUTPUT_REGRESSIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2") # TEMP! JUST TO RENEW NET ACCESS EXTRACTIONS

    net_access_input_dict = load(OUTPUT_EXTRACTIONS_DIR*"access/reg_data/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")

    # Calculate reference values for National Access from Surveys for scatter plot
    national_H_aggregated = net_access_input_dict["H_aggregated"]

    # Extract survey estimates
    NET_CROP_SURVEY_MONTHLY, NPC_SURVEY_MONTHLY, NET_ACCESS_SURVEY_MONTHLY = calc_survey_estimates(input_dict, net_access_input_dict)

    # Find all non-missing entries and add store in aggregation variable
    nonzero_idxs = intersect(findall(.!ismissing.(NET_ACCESS_SURVEY_MONTHLY)), findall(.!ismissing.(NPC_SURVEY_MONTHLY)))#findall(.!ismissing.(NET_ACCESS_SURVEY_MONTHLY))

    NPC_SURVEY_aggregated = vcat(NPC_SURVEY_aggregated, Float64.(NPC_SURVEY_MONTHLY[nonzero_idxs]))
    ACCESS_SURVEY_aggregated = vcat(ACCESS_SURVEY_aggregated, Float64.(NET_ACCESS_SURVEY_MONTHLY[nonzero_idxs]))
end


# %%
npc_val = []
acc_val = []


for ISO in filt_ISOs

    data = JLD2.load("outputs/draws/national/crop_access/$(ISO)_$(YEAR_START)_$(YEAR_END)_post_crop_access.jld2")

    npc_val = vcat(npc_val,data["Γ_MONTHLY_mean_TOTAL"]./data["POPULATION_MONTHLY"])
    acc_val = vcat(acc_val,data["λ_access_mean"])
end

fig, ax = scatter(npc_val, acc_val, markersize = 5)
scatter!(ax, NPC_SURVEY_aggregated, ACCESS_SURVEY_aggregated, markersize = 5)
xlims!(ax, -0.05, 1.05)
ylims!(ax, -0.05, 1.05)

fig
# %% Extract SNF and RASTER Relationships for NPC-Access at national level
raster_timeseries = CSV.read(OUTPUT_DIR*"coverage_timeseries/master_extraction.csv", DataFrame)
raster_timeseries_country = raster_timeseries[raster_timeseries.category .== "Admin0",:]
NPC_SNF_aggregated = raster_timeseries_country.snf_npc_mean
ACCESS_SNF_aggregated = raster_timeseries_country.snf_access_mean
NPC_RASTER_aggregated = raster_timeseries_country.raster_npc_mean
ACCESS_RASTER_aggregated = raster_timeseries_country.raster_access_mean

# # %% Import Survey Raster Level observations
# inla_dataset = CSV.read(OUTPUT_DATAPREP_DIR*"INLA/unfiltered_inla_dataset.csv", DataFrame)
# NPC_SURVEY_CLUSTER = inla_dataset.npc
# ACCESS_SURVEY_CLUSTER = inla_dataset.access

# %% Cluster lookup data
lookup_data = CSV.read(OUTPUT_DIR*"coverage_timeseries/bv_mitn_validation_values.csv", DataFrame)
lookup_data = lookup_data[lookup_data.interview_year .>=2010,:]

NPC_SURVEY_CLUSTER = lookup_data.npc
ACCESS_SURVEY_CLUSTER = lookup_data.access
USE_SURVEY_CLUSTER = lookup_data.use
NPC_RASTER_CLUSTER = lookup_data.mitn_npc
ACCESS_RASTER_CLUSTER = lookup_data.mitn_access
USE_RASTER_CLUSTER = lookup_data.mitn_use

# %% ANALYTICAL MODEL CALCULATIONS - Import model
# Analysis domain/ranges
n_access_samples = 100 # number of access parameter posterior samples to draw
h_max = 10 # Maximum household size
γ_vals = 0:0.01:1 # Range of NPC values
POPULATION_MONTHLY = zeros(length(γ_vals)) .= 1e6 # Dummy population just to plug into model. Will end up being normalised anyway
Γ_MONTHLY_samples_TOTAL = repeat(γ_vals.*POPULATION_MONTHLY,1,n_access_samples)'

# Import MCMC model samples
net_access_chain = load(OUTPUT_REGRESSIONS_DIR*"access/netaccesschains.jld2")
ρ_chain_df = net_access_chain["ρ_chain_df"]
μ_chain_df = net_access_chain["μ_chain_df"]

# %%


# %% Analytical Model 1: Calculate analytical NPC-Access Relationship assuming uniform Household Size Distribution
# Access Regression output data
p_h_uniform = (ones(h_max+1) .= 1 ./(h_max+1))
acc_uniform_samples = sample_net_access(ρ_chain_df, μ_chain_df, p_h_uniform,
                                        POPULATION_MONTHLY, Γ_MONTHLY_samples_TOTAL) # Need to temporarily unscale input by mil for calculating access

# %% Analytical Model 2: Calculate analytical NPC-Access Relationship assuming Poisson Household Size Distribution
h_poisson_mean = 5
p_h_poisson = (h_poisson_mean.*(0:h_max)).*exp.(-h_poisson_mean)./factorial.(0:h_max)
acc_poisson_samples = sample_net_access(ρ_chain_df, μ_chain_df, p_h_poisson,
                                        POPULATION_MONTHLY, Γ_MONTHLY_samples_TOTAL) # Need to temporarily unscale input by mil for calculating access
access_poisson_mean = mean(acc_poisson_samples, dims = 1)[:]

# %% Analytical model 3: Household distribution based on Africa continent level household distribution average across all surveys
p_h_globaldata = JLD2.load(OUTPUT_EXTRACTIONS_DIR*"access/reg_data/aggregated_inputs/netaccess_allsurvey_inputs.jld2")["p_h_globaldata"]
p_h_global = mean(p_h_globaldata, dims = 1)[1,:]./sum(mean(p_h_globaldata, dims = 1)[1,:])
acc_global_samples = sample_net_access(ρ_chain_df, μ_chain_df, p_h_global,
                                        POPULATION_MONTHLY, Γ_MONTHLY_samples_TOTAL) # Need to temporarily unscale input by mil for calculating access
access_global_mean = mean(acc_global_samples, dims = 1)[:]

# %% Analytical model 4: Household distribution based on Africa continent level household distribution separated by Survey
p_h_globaldata = JLD2.load(OUTPUT_EXTRACTIONS_DIR*"access/reg_data/aggregated_inputs/netaccess_allsurvey_inputs.jld2")["p_h_globaldata"]
n_h_profiles = size(p_h_globaldata)[1]
access_profile_samples = Vector{Matrix}(undef, n_h_profiles)
access_profiles = Vector{Vector}(undef, n_h_profiles)

for profile_i in ProgressBar(1:n_h_profiles, leave = false)
    p_h_survey = p_h_globaldata[profile_i,:]
    acc_survey_samples = sample_net_access(ρ_chain_df, μ_chain_df, p_h_survey,
                                        POPULATION_MONTHLY, Γ_MONTHLY_samples_TOTAL) # Need to temporarily unscale input by mil for calculating access
    access_survey_mean = mean(acc_survey_samples, dims = 1)[:]

    access_profile_samples[profile_i] = acc_survey_samples
    access_profiles[profile_i] = access_survey_mean
end

access_profiles_samples_cat = vcat(access_profile_samples...)
access_profiles_mean_cat = hcat(access_profiles...)'

# %% Tally up and calculate confidence bounds
acc_uniform_ci = zeros(length(γ_vals),3)
acc_poisson_ci = zeros(length(γ_vals),3)
acc_global_ci = zeros(length(γ_vals),3)
acc_survey_ci = zeros(length(γ_vals),3)

for i in 1:length(γ_vals)
    acc_uniform_ci[i,[1,3]] = quantile(acc_uniform_samples[:,i], [0.025, 0.975])
    acc_uniform_ci[i,2] = mean(acc_uniform_samples[:,i])

    acc_poisson_ci[i,[1,3]] = quantile(acc_poisson_samples[:,i], [0.025, 0.975])
    acc_poisson_ci[i,2] = mean(acc_poisson_samples[:,i])

    acc_global_ci[i,[1,3]] = quantile(acc_global_samples[:,i], [0.025, 0.975])
    acc_global_ci[i,2] = mean(acc_global_samples[:,i])

    acc_survey_ci[i,[1,3]] = quantile(access_profiles_samples_cat[:,i], [0.025, 0.975])
    acc_survey_ci[i,2] = mean(access_profiles_samples_cat[:,i])
end


###################################################
# %% FIGURE 1: Measures and Metrics Relationships plot
###################################################
fig = Figure(size = (1200,700))
ax1 = Axis(fig[1,1],
            xlabel = "NPC (γ)",
            ylabel = "Access (λ)",
            xlabelsize = 16,
            ylabelsize = 16,
            xticks = 0:0.1:1,
            yticks = 0:0.1:1)
ax2 = Axis(fig[1,2],
            xlabel = "NPC (γ)",
            ylabel = "Access (λ)",
            xlabelsize = 16,
            ylabelsize = 16,
            xticks = 0:0.1:1,
            yticks = 0:0.1:1)
ax3 = Axis(fig[2,1],
            xlabel = "NPC (γ)",
            ylabel = "Access (λ)",
            xlabelsize = 16,
            ylabelsize = 16,
            xticks = 0:0.1:1,
            yticks = 0:0.1:1)
ax4 = Axis(fig[2,2],
            xlabel = "NPC (γ)",
            ylabel = "Access (λ)",
            xlabelsize = 16,
            ylabelsize = 16,
            xticks = 0:0.1:1,
            yticks = 0:0.1:1)
ax5 = Axis(fig[1,3],
                xlabel = "NPC (γ)",
                ylabel = "Utilisation",
                xlabelsize = 16,
                ylabelsize = 16,
                xticks = 0:0.2:2,
                yticks = 0:0.2:2)
ax6 = Axis(fig[2,3],
                xlabel = "Access (λ)",
                ylabel = "Efficiency",
                xlabelsize = 16,
                ylabelsize = 16,
                xticks = 0:0.2:1,
                yticks = 0:0.2:1)
                
Label(fig[0,:], "NPC-Access Relationships", font = :bold,
        halign = :center, tellwidth = false, fontsize = 30)

xlims!(ax1, -0.05,1.05)
xlims!(ax2, -0.05,1.05)
xlims!(ax3, -0.05,1.05)
xlims!(ax4, -0.05,1.05)
xlims!(ax5, -0.05,2.05)
xlims!(ax6, -0.05,1.05)

ylims!(ax1, -0.05,1.05)
ylims!(ax2, -0.05,1.05)
ylims!(ax3, -0.05,1.05)
ylims!(ax4, -0.05,1.05)
ylims!(ax5, -0.05,2.05)
ylims!(ax6, -0.05,1.05)

hidexdecorations!(ax1, ticks = false, grid = false)
hidexdecorations!(ax2, ticks = false, grid = false)
hideydecorations!(ax2, ticks = false, grid = false)
hideydecorations!(ax4, ticks = false, grid = false)

# %%
npc_snf = scatter!(ax1, NPC_SNF_aggregated[1], ACCESS_SNF_aggregated[1],
        markersize = 3, color = (colors[1],1))
scatter!(ax1, NPC_SNF_aggregated, ACCESS_SNF_aggregated,
        markersize = 3, color = (colors[1],0.2))
npc_survey = scatter!(ax1, NPC_SURVEY_aggregated, ACCESS_SURVEY_aggregated,
        markersize = 5, color = (colors[2],1))
Legend(fig[1,1], [npc_snf, npc_survey],
        ["SNF", "Survey"], 
        tellwidth = false, tellheight = false,
        halign = :left, valign = :top,
        padding = (10,10,10,10))

# %%
# scatter!(ax2, NPC_RASTER_aggregated, ACCESS_RASTER_aggregated,
#         markersize = 3, color = (colors[1],0.2))
npc_cluster_survey = scatter!(ax2, NPC_SURVEY_CLUSTER[1], ACCESS_SURVEY_CLUSTER[1],
        markersize = 3, color = (colors[2],1))
scatter!(ax2, NPC_SURVEY_CLUSTER, ACCESS_SURVEY_CLUSTER,
        markersize = 3, color = (colors[2],0.05))
npc_cluster_raster = scatter!(ax2, NPC_RASTER_CLUSTER[1], ACCESS_RASTER_CLUSTER[1],
        markersize = 5, color = (colors[1],1))
scatter!(ax2, NPC_RASTER_CLUSTER, ACCESS_RASTER_CLUSTER,
        markersize = 5, color = (colors[1],0.05))

Legend(fig[1,2], [npc_cluster_raster, npc_cluster_survey],
        ["Raster", "Survey"], 
        tellwidth = false, tellheight = false,
        halign = :left, valign = :top,
        padding = (10,10,10,10))

# %%
unif_line = lines!(ax3, γ_vals, acc_uniform_ci[:,2], color = colors[1])
band!(ax3, γ_vals, acc_uniform_ci[:,1], acc_uniform_ci[:,3], 
        color = (colors[1], 0.2))
poisson_line = lines!(ax3, γ_vals, acc_poisson_ci[:,2], color = colors[2])
band!(ax3, γ_vals, acc_poisson_ci[:,1], acc_poisson_ci[:,3], 
        color = (colors[2], 0.2))
Legend(fig[2,1], [unif_line, poisson_line],
        ["Uniform", "Poisson"], 
        tellwidth = false, tellheight = false,
        halign = :left, valign = :top,
        padding = (10,10,10,10))

# %%
global_line = lines!(ax4, γ_vals, acc_global_ci[:,2], color = colors[1])
band!(ax4, γ_vals, acc_global_ci[:,1], acc_global_ci[:,3], 
        color = (colors[1], 0.2))
scatter!(ax4, NPC_SURVEY_aggregated, ACCESS_SURVEY_aggregated,
        markersize = 5, color = (colors[2],1))
Legend(fig[2,2], [global_line],
        ["Africa HH Mean Demography"], 
        tellwidth = false, tellheight = false,
        halign = :left, valign = :top,
        padding = (10,10,10,10))

# %%

util_cluster_survey = scatter!(ax5, NPC_SURVEY_CLUSTER[1], (USE_SURVEY_CLUSTER[1]/2)/NPC_SURVEY_CLUSTER[1],
        markersize = 3, color = (colors[2],1))
scatter!(ax5, NPC_SURVEY_CLUSTER, (USE_SURVEY_CLUSTER./2)./NPC_SURVEY_CLUSTER,
        markersize = 3, color = (colors[2],0.1))
util_cluster_raster = scatter!(ax5, NPC_RASTER_CLUSTER[1], (USE_RASTER_CLUSTER[1]/2)/NPC_RASTER_CLUSTER[1],
        markersize = 5, color = (colors[1],1))
scatter!(ax5, NPC_RASTER_CLUSTER, (USE_RASTER_CLUSTER./2)./NPC_RASTER_CLUSTER,
        markersize = 5, color = (colors[1],0.05))

Legend(fig[1,3], [util_cluster_raster, util_cluster_survey],
        ["Raster", "Survey"], 
        tellwidth = false, tellheight = false,
        halign = :right, valign = :top,
        padding = (10,10,10,10))

# %%


eff_cluster_survey = scatter!(ax6, ACCESS_SURVEY_CLUSTER[1], USE_SURVEY_CLUSTER[1]/ACCESS_SURVEY_CLUSTER[1],
        markersize = 3, color = (colors[2],1))
scatter!(ax6, ACCESS_SURVEY_CLUSTER, USE_SURVEY_CLUSTER./ACCESS_SURVEY_CLUSTER,
        markersize = 3, color = (colors[2],0.1))
eff_cluster_raster = scatter!(ax6, ACCESS_RASTER_CLUSTER[1], USE_RASTER_CLUSTER[1]/USE_RASTER_CLUSTER[1],
        markersize = 5, color = (colors[1],1))
scatter!(ax6, ACCESS_RASTER_CLUSTER, USE_RASTER_CLUSTER./ACCESS_RASTER_CLUSTER,
        markersize = 5, color = (colors[1],0.05))

Legend(fig[2,3], [eff_cluster_raster, eff_cluster_survey],
        ["Raster", "Survey"], 
        tellwidth = false, tellheight = false,
        halign = :right, valign = :bottom,
        padding = (10,10,10,10))
 
fig

# %%
mkpath(OUTPUT_PLOTS_DIR*"PaperFigures/")
save(OUTPUT_PLOTS_DIR*"PaperFigures/NPC_Access_Relationship.pdf", fig)