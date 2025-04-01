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
theme(:vibrant)

# %% Import helper packages
using PlottingFunctions
using LaTeXStrings

# %% Import modules 
using CropAccessFromSurvey
using DataExtractions
using NetAccessModel

# %%
NPC_SURVEY_aggregated = []
NET_ACCESS_SURVEY_aggregated = []

# %% Compile NPC and Access observation estimates from survey data
ISO_list = String.(CSV.read(raw"C:\Users\ETan\Documents\Prototype Analyses\itn-updated\datasets\ISO_list.csv", DataFrame)[:,1])
reg_results_list = Array{Any, 1}(undef, length(ISO_list))
exclusion_ISOs = ["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD"]
#["CPV","BWA","CAF","COM","GNQ","DJI","ERI","ETH","GAB","GNB","STP","SOM","SDN","SWZ","ZAF","SSD"]

YEAR_START = 2000
YEAR_END = 2023

for i in ProgressBar(1:length(ISO_list), leave = false)
    # Select ISO
    ISO = ISO_list[i]

    if ISO ∈ exclusion_ISOs
        continue
    else
        println(ISO)
        
        # Import required data
        input_dict = load("outputs/extractions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
        regression_dict = load("outputs/regressions/crop/Compact Regressions/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2") # TEMP! JUST TO RENEW NET ACCESS EXTRACTIONS

        net_access_input_dict = load("outputs/extractions/access/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")

        # Extract survey estimates
        NET_CROP_SURVEY_MONTHLY, NPC_SURVEY_MONTHLY, NET_ACCESS_SURVEY_MONTHLY = calc_survey_estimates(input_dict, net_access_input_dict)

        # Find all non-missing entries and add store in aggregation variable
        nonzero_idxs = findall(.!ismissing.(NET_ACCESS_SURVEY_MONTHLY))

        NPC_SURVEY_aggregated = vcat(NPC_SURVEY_aggregated, Float64.(NPC_SURVEY_MONTHLY[nonzero_idxs]))
        NET_ACCESS_SURVEY_aggregated = vcat(NET_ACCESS_SURVEY_aggregated, Float64.(NET_ACCESS_SURVEY_MONTHLY[nonzero_idxs]))
    end
end





###############################################
# %% Do Analytical Poisson NPC vs Access Curve
###############################################

# Extract all net and household distribution data
ISO_list = String.(CSV.read(raw"C:\Users\ETan\Documents\Prototype Analyses\itn-updated\datasets\ISO_list.csv", DataFrame)[:,1])
exclusion_ISOs = ["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD"]
# ["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD"]
# GAB excluded

init_variables = false

p_h_globaldata = zeros(0,10)
ρ_h_globaldata = zeros(0,10)
μ_h_globaldata = zeros(0,10)
γ_globaldata = zeros(0)

for i in 1:length(ISO_list)
    # Select ISO
    ISO = ISO_list[i]

    if ISO ∈ exclusion_ISOs
        continue
    else
        # Import extracted data
        net_access_input_dict = load("outputs/extractions/access/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")
        
        γ_surveys = net_access_input_dict["γ_aggregated"]
        p_h = net_access_input_dict["p_h_aggregated"]
        ρ_h = net_access_input_dict["ρ_h_aggregated"]
        μ_h = net_access_input_dict["μ_h_aggregated"]
        
        h_max = size(p_h)[2]
        # Check if storage variable for aggregated global data has been made.
        if init_variables == false
            p_h_globaldata = zeros(0,h_max)
            ρ_h_globaldata = zeros(0,h_max)
            μ_h_globaldata = zeros(0,h_max)
            γ_globaldata = zeros(0)

            init_variables = true
        end

        # Add to aggregate storage variable
        p_h_globaldata = vcat(p_h_globaldata, p_h)
        ρ_h_globaldata = vcat(ρ_h_globaldata, ρ_h)
        μ_h_globaldata = vcat(μ_h_globaldata, μ_h)
        γ_globaldata = vcat(γ_globaldata, γ_surveys)
    end
end

# %% Parameter ranges to calculate for analytical case
λ = 4.5 # Poisson distribution parameter
γ_vals = Array(0:0.05:1) # Range for NPC to use
h_vals = Array(1:10) # Range of household sizes

# Poisson PDF to use
PDF = (λ.^(h_vals)).*(exp.(-λ))./factorial.(h_vals)
PDF = PDF./sum(PDF)

# %% Construct plot to demonstrate analytical case
fig1 = scatter(p_h_globaldata[1,:], markercolor = 4, markersize = 2, legend = :topright,
                xticks = 1:10, xlabel = "Household Size", ylabel = L"p_h", labelfontsize = 13,
                title = "Household Size Distribution", titlefontsize = 17,
                ylims = (-0.01, 0.25), label = "Survey", legendfontsize = 10)

scatter!(fig1, p_h_globaldata[2:end,:]', markercolor = 4, markersize = 2, legend = :topright,
                xticks = 1:10, xlabel = "Household Size", ylabel = L"p_h", labelfontsize = 13,
                title = "Household Size Distribution", titlefontsize = 17,
                ylims = (-0.01, 0.25), label = nothing, legendfontsize = 10)


plot!(fig1, h_vals, PDF, color = 2, label = "Poisson "*L"\lambda = 4.5", linewidth = 2)

savefig(fig1, "Poisson_HH_distribution.pdf")

# %% Calculate NPC - Access Curve for Analytical Poisson case

# Import regression parameters for Access model
access_chain = load("outputs/regressions/access/netaccesschains.jld2")
ρ_chain_df = access_chain["ρ_chain_df"]
μ_chain_df = access_chain["μ_chain_df"]
p_h = PDF

# Assume a population size of 1e7
POPULATION_MONTHLY = 10000000 .*ones(length(γ_vals))
Γ_MONTHLY_samples = repeat((γ_vals.*POPULATION_MONTHLY)',100)

# Calculate access
access = sample_net_access(ρ_chain_df, μ_chain_df, p_h,
                            POPULATION_MONTHLY, Γ_MONTHLY_samples;
                            n_max = 20)

access_quantiles = zeros(3, size(access)[2])

for i in 1:length(γ_vals)
    access_quantiles[:,i] = quantile(access[:,i], [0.05,0.5,0.95])
end

# %%
fig2 = plot(title = "NPC vs. Access", legendfontsize = 10,
            xlabel = "NPC", ylabel = "Access")

plot!(fig2, γ_vals, access_quantiles[1,:], fillrange = access_quantiles[3,:],
        linealpha = 0, label = nothing, color = 2,
        fillalpha = 0.2, legend = :topleft)

scatter!(fig2, NPC_SURVEY_aggregated, NET_ACCESS_SURVEY_aggregated, 
            markersize = 4, color = 4, label = "Survey")
plot!(fig2, γ_vals, access_quantiles[2,:], color = 2, linewidth = 1.6,
        label = "Model Fit (Poisson)")

savefig(fig2, "Poisson_NPC_ACCESS_fit.pdf")