"""
Author: Eugene Tan
Date Created: 6/6/2025
Last Updated: 6/6/2025
Calculate prediction error of the SNF outputs vs national aggregates
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from TOML file
include(pwd()*"/scripts/read_toml.jl")

# %% Import relevant packages
using ProgressBars
using DataFrames
using CSV
using JLD2
using DateConversions
using NetCropModel
using NetAccessModel
using NetAccessPrediction
using StatsBase
using CairoMakie

# %% Filenames
country_summaries_filename = HOUSEHOLD_NAT_SUMMARY_DATA_FILENAME

# %% Load Datasets
country_summaries = CSV.read(OUTPUT_DATAPREP_DIR*country_summaries_filename, DataFrame)


# %% Define list of countries to plot
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %%
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END
n_months = 12*(YEAR_END - YEAR_START + 1)
n_samples = 150

# %% Storage variable for posterior mean estimates of Γ, NPC and Access
Γ_posterior_mean = Matrix{Float64}(undef, length(filt_ISOs), n_months)
NPC_posterior_mean = Matrix{Float64}(undef, length(filt_ISOs), n_months)
λ_posterior_mean = Matrix{Float64}(undef, length(filt_ISOs), n_months)
λ_survey = Matrix{Union{Missing,Float64}}(undef, length(filt_ISOs), n_months)

#############################################################
# %% Extract mean estimates
#############################################################
for ISO_i in ProgressBar(1:length(filt_ISOs))

    ISO = filt_ISOs[ISO_i]

    input_dict = load(OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
    regression_dict = load(OUTPUT_REGRESSIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")
    net_access_input_dict = load(OUTPUT_EXTRACTIONS_DIR*"access/reg_data/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")
    net_access_chain = load(OUTPUT_REGRESSIONS_DIR*"access/netaccesschains.jld2")

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

    # Access Regression output data
    ρ_chain_df = net_access_chain["ρ_chain_df"]
    μ_chain_df = net_access_chain["μ_chain_df"]
    p_h = mean(net_access_input_dict["p_h_aggregated"], dims = 1)[:]

    # Calculate reference values for National Access from Surveys
    national_H_aggregated = net_access_input_dict["H_aggregated"]
    national_access_aggregated = zeros(size(national_H_aggregated)[1])
    for survey_i in 1:length(national_access_aggregated)
        national_access_aggregated[survey_i] = sum(H_to_access(national_H_aggregated[survey_i,:,:]))
    end

    access_survey_monthidx = net_access_input_dict["survey_monthidx_aggregated"]
    access_survey_monthidx = access_survey_monthidx[findall(.!ismissing.(access_survey_monthidx))]
    NET_ACCESS_SURVEY_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
    NET_ACCESS_SURVEY_MONTHLY[access_survey_monthidx] .= national_access_aggregated

    ##### Extract Net Crop MCMC parameters

    ### Part 1: Final Regressed Values
    # Randomly generate sample indexes to sample from chain
    chain_length = size(chain)[1]
    sample_idxs = sample(1:chain_length, n_samples, replace = false)

    # Extract MCMC draws for parameters
    ϕ_posterior_mean = mean(chain[:, :ϕ])

    α_init_posterior_mean = mean(chain[:,:α_init])
    α_LLIN_posterior_mean = mean(chain[:,:α_LLIN])

    b_net_posterior_mean = mean(Matrix(DataFrame(chain)[:,5:2:5+2*(n_net_types-1)]), dims = 1)[:]
    k_net_posterior_mean = mean(Matrix(DataFrame(chain)[:,6:2:6+2*(n_net_types-1)]), dims = 1)[:]

    if n_missing_nets_vals > 0
        n_missing_nets_posterior_mean = mean(Matrix(DataFrame(chain)[:,(4+2*(n_net_types)+1):end]), dims = 1)[:]
    else
        # No missing data, so just feed in a zero matrix with no columns
        n_missing_nets_posterior_mean = zeros(size(chain)[1], 0)
    end

    ##### Generate Mean Net Crop Trajectory
    Γ_MONTHLY_BYNET, A_BYNET = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                                    ϕ_posterior_mean, b_net_posterior_mean, k_net_posterior_mean,
                                                    α_init_posterior_mean, α_LLIN_posterior_mean,
                                                    n_missing_nets_posterior_mean; 
                                                    monthly_p = monthly_p,
                                                    return_age = true)
    # Get Totals
    Γ_MONTHLY_TOTAL = sum(Γ_MONTHLY_BYNET, dims = 2)[:]
    NPC_MONTHLY_TOTAL = Γ_MONTHLY_TOTAL./POPULATION_MONTHLY
    A_TOTAL = sum(A_BYNET, dims = 3)[:,:,1]

    ##### Calculate Net Access
    λ_access_samples = sample_net_access(ρ_chain_df, μ_chain_df, p_h,
                                        POPULATION_MONTHLY, repeat(Γ_MONTHLY_TOTAL', n_samples, 1)) # Need to temporarily unscale input by mil for calculating access

    λ_access_mean = mean(λ_access_samples, dims = 1)[:]

    ##### Store values into variable
    Γ_posterior_mean[ISO_i,:] = Γ_MONTHLY_TOTAL
    NPC_posterior_mean[ISO_i,:] = NPC_MONTHLY_TOTAL
    λ_posterior_mean[ISO_i,:] = λ_access_mean
    λ_survey[ISO_i,:] = NET_ACCESS_SURVEY_MONTHLY
end

#############################################################
# %% Construct data frame to calculate NPC and Access Aggregate RMSE
#############################################################
npc_df_row_collections = []

for i in 1:size(country_summaries)[1]
    if country_summaries[i,"ISO"] ∈ filt_ISOs
        ISO = country_summaries[i,"ISO"]
        ISO_i = findfirst(filt_ISOs .== ISO)

        month = country_summaries[i,"month"]
        year = country_summaries[i,"year"]
        monthidx = monthyear_to_monthidx(month, year, YEAR_START = YEAR_START)

        df = DataFrame(SurveyId = country_summaries[i,"SurveyId"],
                        Country = country_summaries[i,"Country"],
                        ISO = ISO, month = month, year = year, monthidx = monthidx,
                        sample_size = country_summaries[i,"sample_size"],
                        npc = country_summaries[i,"NPC_mean"],
                        snf_npc = NPC_posterior_mean[ISO_i,monthidx])
        push!(npc_df_row_collections, df)
    end
end

npc_vals = vcat(npc_df_row_collections...)

λ_df_row_collection = []
for idx in findall(.!ismissing.(λ_survey))
    ISO_i, monthidx = idx[1], idx[2]
    month, year_idx = monthidx_to_monthyear(monthidx)
    year = year_idx + YEAR_START - 1

    df = DataFrame(ISO = filt_ISOs[ISO_i], month = month, year = year, monthidx = monthidx,
                    access = λ_survey[idx],
                    snf_access = λ_posterior_mean[ISO_i, monthidx])
    push!(λ_df_row_collection, df)
end

access_vals = vcat(λ_df_row_collection...)

# Calculate metrics
npc_rmse = sqrt(mean((npc_vals.npc .- npc_vals.snf_npc).^2))
npc_cor = cor(npc_vals.npc, npc_vals.snf_npc)
access_rmse = sqrt(mean((access_vals.access .- access_vals.snf_access).^2))
access_cor = cor(access_vals.access, access_vals.snf_access)


#############################################################
# %% Make RMSE Plots
#############################################################
# %% Plotting Theme and general settings
set_theme!(theme_ggplot2())

# Set plot settings
labelsize = 20
titlesize = 25
ms = 5
ma = 0.5
lw = 1
colors = [  colorant"#0082C7", # NPC
            colorant"#E72A3D" # Access
            ]

# Make Plot
fig = Figure(size = (800,400))
ax1 = Axis(fig[1,1], title = "SNF NPC (γ) Fit\n"*
                    "RMSE = $(round(npc_rmse, digits = 3)), "* 
                        "ρ = $(round(npc_cor, digits = 3))",
            xlabel = "Survey NPC (γ)",
            ylabel = "Fitted NPC (γ)",
            titlesize = titlesize, 
            xlabelsize = labelsize,
            ylabelsize = labelsize,
            xticks = 0:0.2:1,
            yticks = 0:0.2:1)
ax2 = Axis(fig[1,2], title = "Access (λ) Fit\n"*
                    "RMSE = $(round(access_rmse, digits = 3)), "*
                    "ρ = $(round(access_cor, digits = 3))",
            xlabel = "Survey Access (λ)",
            ylabel = "Fitted Access (λ)",
            titlesize = titlesize, 
            xlabelsize = labelsize,
            ylabelsize = labelsize,
            xticks = 0:0.2:1,
            yticks = 0:0.2:1)
xlims!(ax1, -0.05, 1.05)
ylims!(ax1, -0.05, 1.05)
xlims!(ax2, -0.05, 1.05)
ylims!(ax2, -0.05, 1.05)

scatter!(ax1, npc_vals.npc, npc_vals.snf_npc, 
            markersize = ms, color = (colors[1], ma))
lines!(ax1, [-0.5,2], [-0.5,2],
            linewidth = lw, linestyle = :dash, color = :black)
scatter!(ax2, access_vals.access, access_vals.snf_access, 
            markersize = ms, color = (colors[2], ma))
lines!(ax2, [-0.5,2], [-0.5,2],
            linewidth = lw, linestyle = :dash, color = :black)

save(OUTPUT_PLOTS_DIR*"PaperFigures/SNF_RMSE.pdf", fig, pdf_version = "1.4")
fig 