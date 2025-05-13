# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")


# %% Import Public Packages
using JLD2
using CSV
using DataFrames
using ProgressBars
using Measures
using Plots
using StatsBase
using Distributions
using NetLoss
using JLD2
using DateConversions

# %% Directories
output_plots_dir = OUTPUT_PLOTS_DIR

# %% Time Parameters
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END

# %% Set Plotting Backend
pythonplot()
theme(:vibrant)



# %% Import helper packages
# using PlottingFunctions
# using NetCropModel
# using NetAccessModel
# using LaTeXStrings
# using DynamicalSystems

############################
# %% Construct Net Efficacy Model
############################

# %% First do Bayesian Inference (old school way) for getting decay behaviour of ITNs
# Import Compiled data
data = CSV.read("datasets/BioefficacyData/Net_BioEfficacy_Dataset.csv", DataFrame)

# Filter out all new-gen LLINs
filt_data = data[findall(data[:, "Net Type"].!="LLIN-2"),:]

t_data = filt_data[:, "Net Age"]
mortality_data = filt_data[:,"Mean Mortality"]
sd_data = filt_data[:,"SD Mortality"]

# Priors for net efficacy model parameters (assuming two parameter decay Weibull function)
α = 9/2
β = 3/2

# Do posterior Bayesian estimate
b_vals = 0:0.05:10
k_vals = 0:0.05:20
joint_post = zeros(length(b_vals), length(k_vals))
for i in ProgressBar(1:length(b_vals), leave = false)
    for j in 1:length(k_vals)
        b = b_vals[i]
        k = k_vals[j]

        d_prior = net_loss_weibull.(t_data, b, k)

        # Get Log components for bayesian inference
        log_p_b_prior = logpdf(Gamma(α, β), b)
        log_p_k_prior = log(pdf(Gamma(α, β), k))
        log_likelihood = zeros(length(t_data))
        for n in 1:length(t_data)
            log_likelihood[n] = log(pdf(Normal(d_prior[n], sd_data[n]),mortality_data[n]))
        end

        log_likelihood

        # Correct Infs because they will be low valued
        log_likelihood[findall(isinf.(log_likelihood))] .= -100

        log_posterior = mean(log_likelihood) + log_p_b_prior + log_p_k_prior

        joint_post[i,j] = exp(log_posterior)
    end
end

# Fix NaNs and Normalise to probability dist
joint_post[findall(isnan.(joint_post))] .= 0
joint_post = joint_post./sum(joint_post)
# heatmap(k_vals, b_vals , joint_post) #Plo just to check

# Get MAP estimate for b and k parameters of net efficacy model
b_idx, k_idx = argmax(joint_post)[1], argmax(joint_post)[2]
b_est = b_vals[b_idx]
k_est = k_vals[k_idx]

# %%
t_vals = 0:0.01:5
fig = plot(xlabel = "Years", ylabel = "24h Mortality (η)",
            title = "Waning Net Efficacy")
plot!(fig, t_vals,net_loss_weibull.(t_vals, b_est, k_est), label = "MAP")
scatter!(fig, t_data, mortality_data, label = "Raw Data (Paper)")

savefig(fig, output_plots_dir*"bioefficacy_fit.pdf")

############################
# %% Insecticide Resistance Model
############################

# Insecticide resistance time series
n_months = (YEAR_END-YEAR_START+1)*12
year_vals = (1:n_months)./12
β_t = 1 .- net_loss_weibull.(year_vals, 18.0, 6.0)
ξ = 0.4 # Maximum Kill rate reduction
ρ_t = 1 .- ξ.*β_t # Resistance effects
fig = plot(1:n_months, zeros(length(ρ_t)), fillrange = ρ_t, 
            ylims = (-0.05, 1.05), 
            legend = false,
            xlabel = "Year", ylabel = "Insecticide Resistance (ρ)",
            xticks = (1:12:n_months, YEAR_START:1:YEAR_END),
            xtickfontrotation = 90,
            color = colorant"#5BA58E", lw = 2.5,
            label = nothing,
            fillalpha = 0.18, linealpha = 0)
plot!(fig, 1:n_months, ρ_t,
        linewidth = 1.2,
        color = colorant"#5BA58E",
        label = "IR Penalty Rate")
plot!(fig, 1:n_months, ρ_t, fillrange = ones(length(ρ_t)),
        fillcolor = colorant"#AC0E26", fillalpha = 1,
        linealpha = 0, fillstyle =  :/)
plot!(fig, 1:n_months, ρ_t, fillrange = ones(length(ρ_t)),
        fillcolor = colorant"#AC0E26", fillalpha = 0.1,
        linealpha = 0, )

# %%
savefig(fig, output_plots_dir*"IR_weibull_model.pdf")

############################
# %% Simulate Insecticide Resistance fo each country and save plot
############################
fig_collection = []

ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

for ISO_i in ProgressBar(1:length(filt_ISOs))
    # Select country
    ISO = filt_ISOs[ISO_i]

    # %% Calculate age based penalty
    # Construct net age matrix given month bounds
    M_age = zeros(n_months, n_months)
    for i in 1:n_months # Current Time
        for j in 1:n_months # Birth time
            if j>i
                M_age[i,j] = 0
            else
                M_age[i,j] = (i-j)/12
            end
        end
    end

    # Use net age matrix to construct theoretical net loss matrix due to net wear
    η_netwear = net_loss_weibull.(M_age, b_est, k_est)

    # Calculate penalised net demography matrix and penalty timeseries
    snf_post_dir = OUTPUT_DRAWS_DIR
    snf_post_draws = JLD2.load(snf_post_dir*"subnational/$(ISO)_SUBNAT_draws.jld2")
    subnat_A_mean = snf_post_draws["merged_outputs"][1]["COMBINED_A_TOTAL_mean"]
    age_penalty_timeseries = (sum(subnat_A_mean.*η_netwear, dims = 2)./sum(subnat_A_mean, dims = 2))[:]

    # %% Calculate penalised use time series
    raster_mean_timeseries = CSV.read(OUTPUT_DIR*"coverage_timeseries/master_extraction_updated.csv", DataFrame)
    n_months = (YEAR_END - YEAR_START + 1)*12
    use_timeseries = (zeros(n_months, 3) .= NaN)
    for i in 1:n_months
        month, year_idx = monthidx_to_monthyear(i)
        year = YEAR_START + (year_idx - 1)
        data_slice = raster_mean_timeseries[ raster_mean_timeseries.ISO .== ISO .&&
                                raster_mean_timeseries.month .== month .&&
                                raster_mean_timeseries.year .== year,:]
        use_mean = sum(data_slice.population .* data_slice.raster_use_mean)./sum(data_slice.population)
        use_lower = sum(data_slice.population .* data_slice.raster_use_95lower)./sum(data_slice.population)
        use_upper = sum(data_slice.population .* data_slice.raster_use_95upper)./sum(data_slice.population)
        use_timeseries[i,:] .= use_lower, use_mean, use_upper
    end

    penalised_use_timeseries_age = zeros(size(use_timeseries))
    penalised_use_timeseries_IR = zeros(size(use_timeseries))
    penalised_use_timeseries_age_IR = zeros(size(use_timeseries))
    for j in 1:size(use_timeseries)[2]
        penalised_use_timeseries_age[:,j] = use_timeseries[:,j] .* age_penalty_timeseries
        penalised_use_timeseries_IR[:,j] = use_timeseries[:,j] .* ρ_t
        penalised_use_timeseries_age_IR[:,j] = use_timeseries[:,j] .* age_penalty_timeseries .* ρ_t
    end


    ############################
    # %% Construct plot of effective coverage
    ############################
    lw = 1.2
    la = 0.8
    fig = plot(title = "$(ISO)",
                xticks = (1:12:n_months, YEAR_START:YEAR_END),
                xlabel = "Year", xtickfontrotation = 90,
                ylabel = "Effective Net Coverage", legend = :outerright,
                ylims = (-0.05, 1.05),
                minorgrid = false, size = (680,400))
    plot!(fig, use_timeseries[:,2], color = :black, label = "Unpenalised",
            linewidth = lw)
    plot!(fig, penalised_use_timeseries_age[:,2], color = colorant"#2E72C6", label = "Age Penalty",
                linewidth = lw, linestyle = :dash, linealpha = la)
    plot!(fig, penalised_use_timeseries_IR[:,2], color = colorant"#AC0E26", label = "IR Penalty",
                linewidth = lw, linestyle = :dash, linealpha = la)
    plot!(fig, penalised_use_timeseries_age_IR[:,2], color = colorant"#A430B1", label = "Age + IR",
                linewidth = lw)

    plot!(zeros(length(ρ_t)), fillrange = ρ_t, 
            fillalpha = 0.15, linealpha = 0,
            color = colorant"#5BA58E", lw = 2,
            label = nothing)
    plot!(ρ_t, linewidth = lw, 
            color = colorant"#5BA58E", lw = 2, 
            label = "IR Penalty Rate (ρ)")

    savefig(fig, OUTPUT_PLOTS_DIR*"raster summaries/countries/$(ISO)_EffectiveCoverage.pdf")
    
    # Add plots to collection
    push!(fig_collection, fig)
end

# %% Combined all country plot

# Plot layout settings
layout = (9,5)
figsize = (2560,1800)

# Make plots and save
fig_comb = plot(fig_collection..., layout = layout, size = figsize)
savefig(fig_comb, OUTPUT_PLOTS_DIR*"raster summaries/coverage_efficacy.pdf")

























