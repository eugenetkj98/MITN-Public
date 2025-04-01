"""
Author: Eugene Tan
Date Created: 24/3/2025
Last Updated: 24/3/2025
Takes input timeseries and a reference country and produces randomly sampled timeseries for ITN nets
"""
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Load packages
using JLD2
using CSV
using DataFrames

# %% Other useful packages
using Dates
using ProgressBars
using StatsBase

# %% Custom modules
using DateConversions
using NetCropModel
using NetLoss
using NetCropPrediction

# %% Plotting Tools
using Plots

########################################################
# %% Define required directories and prediction settings
########################################################

ISO = "NGA" # Reference country to get parameter data from
input_timeseries_file = "datasets/forward_prediction/input_timeseries.csv"
mitn_posterior_chain = "outputs/regressions/crop/2000_2023/$(ISO)_2000_2023_cropchains.jld2"
n_samples = 5 # Number of posterior draws of time series
output_dir = "outputs/predictions/"
counter = 1
output_filename = "$(string(Dates.today()))_run_$(counter).csv" # NEED TO THINK OF A WAY TO AUTOMATE ADDING FILES


########################################################
# %% Load Input Time Series and Scrape required meta-data
########################################################
# Meta Data
input_timeseries = CSV.read(input_timeseries_file, DataFrame)
NET_NAMES = names(input_timeseries)[4:end]
n_net_types = length(NET_NAMES)
YEAR_START, YEAR_END = input_timeseries.year[[1,end]]
YEARS_ANNUAL = YEAR_START:YEAR_END

# Construct Delivery and Distribution Time Series
DELIVERIES_ANNUAL = Array(input_timeseries.LLIN_delivery)
DISTRIBUTION_ANNUAL = Array(input_timeseries[:,NET_NAMES])

########################################################
# %% Load regression chain and construct lookup datafile for posterior draws
########################################################
# Regressed YEAR RANGES
reg_YEAR_START = 2000
reg_YEAR_END = 2023

# Load required base data
mitn_model = JLD2.load(mitn_posterior_chain)
reg_monthly_p = mitn_model["monthly_p"]
reg_NET_NAMES = mitn_model["NET_NAMES"]
reg_ϕ_chain = mitn_model["chain"].ϕ
reg_α_init_chain = mitn_model["chain"].α_init
reg_α_LLIN_chain = mitn_model["chain"].α_LLIN
reg_τ_chain = rename!(mitn_model["chain"][:,5:2:5+2*(length(reg_NET_NAMES).-1)], reg_NET_NAMES)
reg_κ_chain = rename!(mitn_model["chain"][:,6:2:6+2*(length(reg_NET_NAMES).-1)], reg_NET_NAMES)

# Extrapolate monthly_p to required annual ranges
# Make storage variable
n_months = (YEAR_END-YEAR_START+1)*12
monthly_p = zeros(n_months)

# Fill in known range from data
reg_idx_start = monthyear_to_monthidx(1, reg_YEAR_START, YEAR_START = YEAR_START)
reg_idx_end = monthyear_to_monthidx(12, reg_YEAR_END, YEAR_START = YEAR_START)
monthly_p[reg_idx_start:reg_idx_end] = mitn_model["monthly_p"]

# For unknown range, copy from the last 12 months of the regression which is based on the average annual trend
annual_distribution_pattern = reg_monthly_p[end-11:end] # Get monthly pattern from last 12 months of regression
fill_in_years = setdiff(YEAR_START:YEAR_END, reg_YEAR_START:reg_YEAR_END)
for year in fill_in_years
    monthidx_start, monthidx_end = monthyear_to_monthidx.([1,12], year, YEAR_START = YEAR_START)
    monthly_p[monthidx_start:monthidx_end] = annual_distribution_pattern
end

# Construct posterior chains for parameter samples to align with desired NET_NAMES
# If prediction net type is not found in existing list, then just pretend it behaves similarly to an LLIN

τ_chain = DataFrame()
κ_chain = DataFrame()

for net_type in NET_NAMES
    if net_type ∈ reg_NET_NAMES
        τ_chain = hcat(τ_chain, reg_τ_chain[:, net_type], makeunique = true)
        κ_chain = hcat(κ_chain, reg_κ_chain[:, net_type], makeunique = true)
    else
        net_type_idx_LLIN = findfirst(reg_NET_NAMES .== "LLIN")
        τ_chain = hcat(τ_chain, reg_τ_chain[:, net_type_idx_LLIN], makeunique = true)
        κ_chain = hcat(κ_chain, reg_κ_chain[:, net_type_idx_LLIN], makeunique = true)
    end
end

rename!(τ_chain, NET_NAMES)
rename!(κ_chain, NET_NAMES)

########################################################
# %% Do forward simulation for time series
########################################################

# %% Define forward function
function mitn_national_predict(YEARS_ANNUAL,
                                DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                ϕ_est, τ_net_est, κ_net_est, 
                                α_LLIN_est;
                                monthly_p = nothing)

    # Get Number of net types from parameter MCMC chain
    n_net_types = length(NET_NAMES)
    
    # Calculate total net Distributions
    DISTRIBUTION_ANNUAL_TOTAL = sum(Array(DISTRIBUTION_ANNUAL), dims = 2)[:]

    # Declare intermediate variables
    COUNTRY_LLIN_STOCK_ANNUAL_SOY = missings(Float64,length(YEARS_ANNUAL)) # Available stock at start of year following annual delivery
    COUNTRY_LLIN_STOCK_ANNUAL_EOY = missings(Float64,length(YEARS_ANNUAL)) # Available stock at end of year preceding annual delivery

    UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL = zeros(Float64,length(YEARS_ANNUAL)) # Unadjusted number of distributed nets, totalled
    ADJUSTED_DISTRIBUTION_ANNUAL_TOTAL = zeros(Float64,length(YEARS_ANNUAL)) # Adjusted number of distributed nets, totalled

    UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL = zeros(Float64,length(YEARS_ANNUAL)) # Unadjusted number of distributed LLINs
    ADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL = zeros(Float64,length(YEARS_ANNUAL)) # Adjusted number of distributed LLINs

    UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET = zeros(Float64,length(YEARS_ANNUAL), n_net_types) # Unadjusted number of distributed nets, split by type
    ADJUSTED_DISTRIBUTION_ANNUAL_BYNET = zeros(Float64,length(YEARS_ANNUAL), n_net_types) # Adjusted number of distributed nets, split by type
    
    # Annual stock and flow calculation for delivery and distribution
    for i in 1:length(YEARS_ANNUAL)
        # Initialise LLIN STOCK delivery
        if i == 1
            COUNTRY_LLIN_STOCK_ANNUAL_SOY[i] = DELIVERIES_ANNUAL[i]
        else
            COUNTRY_LLIN_STOCK_ANNUAL_SOY[i] = COUNTRY_LLIN_STOCK_ANNUAL_EOY[i-1] + DELIVERIES_ANNUAL[i]
        end

        # Distributions are always defined i.e. no missing data in input
        UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] = DISTRIBUTION_ANNUAL_TOTAL[i]

        if UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] == 0 # No nets distributed in the year
            UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,:] .= 0 # i.e. no distributions
        else
            if YEARS_ANNUAL[i] < 2010 # Need to account for cITNs
                # Proportion of distributions attributed to cITNs
                UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = ((1/α_LLIN_est)/(1+1/α_LLIN_est))*UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]
                
                # Set remaining nets distributed to LLINs
                if n_net_types > 1 # i.e there are LLINs
                    if (UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]-UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1])>COUNTRY_LLIN_STOCK_ANNUAL_SOY[i]

                        # Number of actually distributed LLINs were based on current stock
                        UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2] = COUNTRY_LLIN_STOCK_ANNUAL_SOY[i]

                        # The rest of the distributions must have been cITNs
                        UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]-UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2]
                    else
                        UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2] = UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]-UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1]
                    end

                    # Fill in remainder of vector
                    if n_net_types > 2
                        UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,3:end] .= 0
                    end

                    # Tally up LLINs
                    UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i] = sum(UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2:end])
                end
            else # all are LLINs

                # Tally up LLINs Need to limit by available stock
                UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i] = min(COUNTRY_LLIN_STOCK_ANNUAL_SOY[i],sum(DISTRIBUTION_ANNUAL[i,:]))
                UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,:] = (UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i]/sum(DISTRIBUTION_ANNUAL[i,:])).*DISTRIBUTION_ANNUAL[i,:]
            end
        end

        # Calculate excess and adjust total distributions
        POTENTIAL_EXCESS_STOCK = max(COUNTRY_LLIN_STOCK_ANNUAL_SOY[i] - UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i],0)
        ADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i] = UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i] + ϕ_est*POTENTIAL_EXCESS_STOCK

        # Interpolate total distributions to individual net types
        if UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] == 0 # Deal with NaNs
            ADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,:] .= 0
        else 
            # cITN's don't get adjusted
            ADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1]

            # LLINs get adjusted by ϕ
            if ADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i] > 0
                if UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i] == 0
                    ADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2:end] .= 0
                else
                    ADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2:end] .= (ADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i]/UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i]).*UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2:end]
                end
            else
                ADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2:end] .= 0
            end

            ADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] = sum(ADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,:])
        end

        COUNTRY_LLIN_STOCK_ANNUAL_EOY[i] = COUNTRY_LLIN_STOCK_ANNUAL_SOY[i]-ADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i]
    end

    # Normalise monthly proportions on annual scale
    if isnothing(monthly_p)
        monthly_p = ones(length(YEARS_ANNUAL)*12) # Base case assumes regular uniform monthly equal distributions
    end

    effective_monthly_proportions = zeros(length(YEARS_ANNUAL)*12)
    
    for i in 1:length(YEARS_ANNUAL)
        effective_monthly_proportions[((i-1)*12+1):(i*12)] = monthly_p[((i-1)*12+1):(i*12)]./sum(monthly_p[((i-1)*12+1):(i*12)])
    end

    # Disaggregate annual distributions into monthly distributions with potential seasonality
    DISTRIBUTION_MONTHLY_BYNET = missings(Int64, length(YEARS_ANNUAL)*12, n_net_types) # Storage variable for distributions
    for i in 1:(length(YEARS_ANNUAL)*12)
        month_i, year_i = monthidx_to_monthyear(i)
        DISTRIBUTION_MONTHLY_BYNET[i,:] = round.(ADJUSTED_DISTRIBUTION_ANNUAL_BYNET[year_i,:].*effective_monthly_proportions[i])
    end

    # Survival net crop matrix calculation
    # A_BYNET[i,j] = Number of nets born in time j that has survived up to time i
    A_BYNET = zeros(length(YEARS_ANNUAL)*12, length(YEARS_ANNUAL)*12, n_net_types)

    for n in 1:n_net_types
        for j in 1:(length(YEARS_ANNUAL)*12)
            for i in j:(length(YEARS_ANNUAL)*12)
                if j == i # Monthly distribution
                    A_BYNET[i,j,n] = DISTRIBUTION_MONTHLY_BYNET[i,n]
                else # Calculate from Survival
                    # Probability of a net to survive most recent time step
                    net_survival_probability = net_loss_compact((i-j)/12,  τ_net_est[n], κ_net_est[n])
                    if isnan(net_survival_probability)
                        net_survival_probability = 0
                    end
                    A_BYNET[i,j,n] = A_BYNET[j,j,n]*net_survival_probability # Number from Previous
                end
            end
        end
    end

    Γ_MONTHLY_BYNET = sum(A_BYNET, dims = 2)[:,1,:]

    return Γ_MONTHLY_BYNET, A_BYNET
end


function mitn_national_dist_predict(YEARS_ANNUAL,
                                DISTRIBUTION_ANNUAL,
                                τ_net_est, κ_net_est;
                                monthly_p = monthly_p)

    # Get Number of net types from parameter MCMC chain
    n_net_types = length(NET_NAMES)
    
    # Calculate total net Distributions
    DISTRIBUTION_ANNUAL_TOTAL = sum(Array(DISTRIBUTION_ANNUAL), dims = 2)[:]

    # Normalise monthly proportions on annual scale
    if isnothing(monthly_p)
        monthly_p = ones(length(YEARS_ANNUAL)*12) # Base case assumes regular uniform monthly equal distributions
    end

    effective_monthly_proportions = zeros(length(YEARS_ANNUAL)*12)
    
    for i in 1:length(YEARS_ANNUAL)
        effective_monthly_proportions[((i-1)*12+1):(i*12)] = monthly_p[((i-1)*12+1):(i*12)]./sum(monthly_p[((i-1)*12+1):(i*12)])
    end

    # Disaggregate annual distributions into monthly distributions with potential seasonality
    DISTRIBUTION_MONTHLY_BYNET = missings(Int64, length(YEARS_ANNUAL)*12, n_net_types) # Storage variable for distributions
    for i in 1:(length(YEARS_ANNUAL)*12)
        month_i, year_i = monthidx_to_monthyear(i)
        DISTRIBUTION_MONTHLY_BYNET[i,:] = round.(DISTRIBUTION_ANNUAL[year_i,:].*effective_monthly_proportions[i])
    end
    
    # Survival net crop matrix calculation
    # A_BYNET[i,j] = Number of nets born in time j that has survived up to time i
    A_BYNET = zeros(length(YEARS_ANNUAL)*12, length(YEARS_ANNUAL)*12, n_net_types)

    for n in 1:n_net_types
        for j in 1:(length(YEARS_ANNUAL)*12)
            for i in j:(length(YEARS_ANNUAL)*12)
                if j == i # Monthly distribution
                    A_BYNET[i,j,n] = DISTRIBUTION_MONTHLY_BYNET[i,n]
                else # Calculate from Survival
                    # Probability of a net to survive most recent time step
                    net_survival_probability = net_loss_compact((i-j)/12,  τ_net_est[n], κ_net_est[n])
                    if isnan(net_survival_probability)
                        net_survival_probability = 0
                    end
                    A_BYNET[i,j,n] = A_BYNET[j,j,n]*net_survival_probability # Number from Previous
                end
            end
        end
    end

    Γ_MONTHLY_BYNET = sum(A_BYNET, dims = 2)[:,1,:]

    return Γ_MONTHLY_BYNET, A_BYNET
end


# %% Do prediction for multiple samples
n_samples = 100
# Storage variable for sampled posterior prediction
Γ_BYNET_pred = zeros(n_samples, n_months, n_net_types)
A_BYNET_pred = zeros(n_samples, n_months, n_months, n_net_types)

Γ_BYNET_dist_pred = zeros(n_samples, n_months, n_net_types)
A_BYNET_dist_pred = zeros(n_samples, n_months, n_months, n_net_types)

# Select idx values to sample from
sampled_idxs = rand(1:length(reg_ϕ_chain),n_samples)

# Perform predictive sampling for each selected idx in MCMC chain 
for i in ProgressBar(1:n_samples, leave = false)
    # Select sampled index
    idx = sampled_idxs[i]

    # Extract relevant chain values
    ϕ_est = reg_ϕ_chain[idx]
    α_LLIN_est = reg_α_LLIN_chain[idx]
    τ_net_est = τ_chain[idx,:]
    κ_net_est = κ_chain[idx,:]

    # Perform prediction
    Γ_MONTHLY_BYNET_sample, A_BYNET_sample = mitn_national_predict(YEARS_ANNUAL,
                                                        DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                                        ϕ_est, τ_net_est, κ_net_est, 
                                                        α_LLIN_est;
                                                        monthly_p = monthly_p)
                    
    Γ_MONTHLY_BYNET_dist_sample, A_BYNET_dist_sample = mitn_national_dist_predict(YEARS_ANNUAL,
                                                                        DISTRIBUTION_ANNUAL,
                                                                        τ_net_est, κ_net_est;
                                                                        monthly_p = monthly_p)

    # Save prediction to storage variable
    Γ_BYNET_pred[i,:,:] = Γ_MONTHLY_BYNET_sample
    A_BYNET_pred[i,:,:,:] = A_BYNET_sample

    Γ_BYNET_dist_pred[i,:,:] = Γ_MONTHLY_BYNET_dist_sample
    A_BYNET_dist_pred[i,:,:,:] = A_BYNET_dist_sample
end

# %% Calculate 95% CI
Γ_TOTAL_pred_sample = sum(Γ_BYNET_pred, dims = 3)[:,:,1]
Γ_TOTAL_dist_pred_sample = sum(Γ_BYNET_dist_pred, dims = 3)[:,:,1]

# Storag variable
Γ_BYNET_pred_CI = zeros(size(Γ_TOTAL_pred_sample)[2],3, n_net_types)
Γ_TOTAL_pred_CI = zeros(size(Γ_TOTAL_pred_sample)[2],3)

Γ_BYNET_dist_pred_CI = zeros(size(Γ_TOTAL_dist_pred_sample)[2],3, n_net_types)
Γ_TOTAL_dist_pred_CI = zeros(size(Γ_TOTAL_dist_pred_sample)[2],3)

for i in 1:(size(Γ_TOTAL_pred_CI)[1])
    Γ_TOTAL_pred_CI[i,:] = quantile(Γ_TOTAL_pred_sample[:,i],[0.025, 0.5, 0.975])
    Γ_TOTAL_dist_pred_CI[i,:] = quantile(Γ_TOTAL_dist_pred_sample[:,i],[0.025, 0.5, 0.975])
    for j in 1:n_net_types
        Γ_BYNET_pred_CI[i,:,j] = quantile(Γ_BYNET_pred[:,i,j],[0.025, 0.5, 0.975])
        Γ_BYNET_dist_pred_CI[i,:,j] = quantile(Γ_BYNET_dist_pred[:,i,j],[0.025, 0.5, 0.975])
    end
end

########################################################
# %% Save output as file
########################################################
# NEED TO DECIDE ON IDEAL DATA STRUCTURE FOR OUTPUTS


########################################################
# %% Make Plot
########################################################
# Import saved data
# Posterior Draws


# Reference Survey Data
survey_net_crop = JLD2.load("outputs/extractions/crop/2000_2023/$(ISO)_2000_2023_cropextract.jld2")["NET_CROP_MONTHLY"]
survey_months = JLD2.load("outputs/extractions/crop/2000_2023/$(ISO)_2000_2023_cropextract.jld2")["MONTHS_MONTHLY"]

# Plot visual settings
pythonplot()
theme(:vibrant)
ms = 4
# Color settings
total_net_col = :black
total_dist_net_col = :blue
colors = [  colorant"#005684",
                colorant"#00976A",
                colorant"#E72A3D",
                colorant"#F7B801",
                colorant"#7018B3",
                ]
component_alpha_mult = 0.4 # Multiplier to reduce component nets time series

# %% Validation/sanity check against survey Data

fig = plot(xticks = (1:12:n_months, YEAR_START:YEAR_END), 
            xlabel = "Year", ylabel = "Net Crop (mil)",
            xtickfontrotation = 90, legend = :topleft,
            title = "Predicted Net Crop ($ISO)")
plot!(fig, Γ_TOTAL_pred_CI[:,1]./1e6, fillrange = Γ_TOTAL_pred_CI[:,3]./1e6, 
            fillcolor = total_net_col, fillalpha = 0.2,
            linealpha = 0, label = nothing)
plot!(fig, Γ_TOTAL_pred_CI[:,2]./1e6,
            linecolor = total_net_col, linealpha = 1,
            label = "Total Nets", linewidth = 1.2)

plot!(fig, Γ_TOTAL_dist_pred_CI[:,1]./1e6, fillrange = Γ_TOTAL_dist_pred_CI[:,3]./1e6, 
            fillcolor = total_dist_net_col, fillalpha = 0.2,
            linealpha = 0, label = nothing)
plot!(fig, Γ_TOTAL_dist_pred_CI[:,2]./1e6,
            linecolor = total_dist_net_col, linealpha = 1,
            label = "Total Nets (Dist only)", linewidth = 1.2)
scatter!(fig, survey_months, survey_net_crop./1e6,
            label = "Survey", markersize = ms, markercolor = total_net_col)

for j in 1:n_net_types
    plot!(fig, Γ_BYNET_pred_CI[:,1,j]./1e6, fillrange = Γ_BYNET_pred_CI[:,3,j]./1e6, 
            fillcolor = colors[j], fillalpha = 0.2*component_alpha_mult,
            linealpha = 0, label = nothing)
    plot!(fig, Γ_BYNET_pred_CI[:,2,j]./1e6,
            linecolor = colors[j], linealpha = 1*component_alpha_mult,
            label = NET_NAMES[j], linewidth = 1.2)
end
            
vline!(fig, [(reg_YEAR_END - YEAR_START)*12], label = nothing, 
            color = :red, linestyle = :dash)