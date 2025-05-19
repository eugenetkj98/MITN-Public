"""
Author: Eugene Tan
Date Created: 13/8/2024
Last Updated: 30/9/2024
Code that defined the PPL model for the estimating the national net crop time series from extracted data. Also contains a forward evaluation mode to perform posterior predictions. This versions accounts for multiple net types.
"""

module NetCropModel
export monthly_reduced_itn
export model_evolve_forward

# %% Import filenames and directories from TOML file
include(pwd()*"/scripts/read_toml.jl")

using Distributions
using StatsBase
using LinearAlgebra
using Missings

using Turing
using AdvancedMH

using DateConversions
using NetLoss



########################################
# %% Forward evolution of monthly reduced ITN model
########################################
"""
    monthly_reduced_itn(YEARS_ANNUAL, MONTHS_MONTHLY,
                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                    NET_CROP_MONTHLY; 
                                    monthly_p = nothing)

Function of statistical evolution model for net crop component of the stock and flow model. This function uses the `@model` macro to pass the function into the Sampler construction in Turing.jl. The variable `monthly_p` contains real values for the relative proportions of each annual distribution quota that is distributed in each given month.
"""
@model function monthly_reduced_itn(YEARS_ANNUAL, MONTHS_MONTHLY,
                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                    cITN_CROP_MONTHLY,
                                    NET_CROP_MONTHLY, NET_CROP_STD_MONTHLY; 
                                    monthly_p = nothing,
                                    MISSING_NETS_SCALE = MISSING_NETS_SCALE)


    # Declare all priors to inference
    β ~ Uniform(20,24) # Hyperprior for ϕ
    ϕ ~ truncated(Beta(2, β), 0, 0.25) # Redistribution parameter
    # β = rand(Uniform(20,24)) # Hyperprior for ϕ
    # ϕ = rand(Beta(2, β)) # Redistribution parameter

    α_init ~ LogNormal(0,1) # N_initial nets
    α_LLIN ~ Uniform(0,1) # N_initial nets
    # α_init = rand(LogNormal(0,1))
    # α_LLIN = rand(Uniform(0,1)) # N_initial nets

    n_net_types = size(DISTRIBUTION_ANNUAL)[2]-1
    b_nets = Vector{Real}(undef, n_net_types)
    k_nets = Vector{Real}(undef, n_net_types)

    for i in 1:n_net_types
        # b_nets[i] ~ LogNormal(0.5,0.5) # Location parameter Favouring Sigmoid
        # k_nets[i] ~ LogNormal(1.5,0.2) # Location parameter Favouring Sigmoid

        b_nets[i] ~ Uniform(0.1,20.7) # Location parameter Favouring Sigmoid
        k_nets[i] ~ Uniform(1,30) # Location parameter Favouring Sigmoid
    end

    missing_dist_idxs = findall(ismissing.(DISTRIBUTION_ANNUAL[:,1]))
    missing_dist_vals = Vector{Real}(undef, length(missing_dist_idxs))
    for i in 1:length(missing_dist_vals)
        missing_dist_vals[i] ~ LogNormal(0,1)
        # missing_dist_vals[i] = rand(Uniform(0,1))
    end

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
            COUNTRY_LLIN_STOCK_ANNUAL_SOY[i] = DELIVERIES_ANNUAL[i] #α_init*MISSING_NETS_SCALE + 
        else
            COUNTRY_LLIN_STOCK_ANNUAL_SOY[i] = COUNTRY_LLIN_STOCK_ANNUAL_EOY[i-1] + DELIVERIES_ANNUAL[i]
        end

        if ismissing(DISTRIBUTION_ANNUAL[i,1])
            # Retrieve random variable representing missing distribution
            dist_idx = findfirst(missing_dist_idxs.==i)

            # Assume all distributed nets were LLINs or cITNs depending on starting year (Most conservative assumption)
            if YEARS_ANNUAL[i] < 2010
                # Needed to avoid missings in forward evolution

                if COUNTRY_LLIN_STOCK_ANNUAL_SOY[i] == 0 # No LLIN stock reference (especially in early years)
                    UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] = α_init*MISSING_NETS_SCALE
                    # Assign all to cITNs since no LLIN stock
                    UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]
                else
                    UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] = min((1+1/α_LLIN)*COUNTRY_LLIN_STOCK_ANNUAL_SOY[i], missing_dist_vals[dist_idx]*MISSING_NETS_SCALE)

                    # Assign a fraction of distribution towards cITN
                    UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = ((1/α_LLIN)/(1+1/α_LLIN))*UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]
                end
                
                
                # Set remaining nets distributed to LLINs
                if n_net_types > 1
                    # UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2:end] .= 0
                    UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2] = UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]-UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1]
                    if n_net_types > 2
                        UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,3:end] .= 0
                    end
                end
            else
                # If missing distribution data for stock and flow, assume dist is equal to delivery
                UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] = min(COUNTRY_LLIN_STOCK_ANNUAL_SOY[i], missing_dist_vals[dist_idx]*MISSING_NETS_SCALE)

                # Assign to LLINs (second column)
                UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2] = UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]
                
                # Set remaining nets distributed to 0
                UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = 0
                if n_net_types > 2
                    UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,3:end] .= 0
                end
            end

            # Update LLIN distribution total
            UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i] = sum(UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2:end])
        else
            UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] = DISTRIBUTION_ANNUAL[i,1]#min(COUNTRY_LLIN_STOCK_ANNUAL_SOY[i], DISTRIBUTION_ANNUAL[i,1])
            if UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] == 0 # No nets distributed in the year
                UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,:] .= 0 # i.e. no distributions
            else
                if YEARS_ANNUAL[i] < 2010 # Need to account for cITNs
                    # Proportion of distributions attributed to cITNs
                    UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = ((1/α_LLIN)/(1+1/α_LLIN))*UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]
                    # ADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1]
                    
                    # Set remaining nets distributed to LLINs
                    if n_net_types > 1 # i.e there are LLINs
                        # UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2:end] .= 0
                        if (UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]-UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1])>COUNTRY_LLIN_STOCK_ANNUAL_SOY[i]
                            # Number of actually distributed LLINs were based on current stock
                            UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2] = COUNTRY_LLIN_STOCK_ANNUAL_SOY[i]

                            # The rest of the distributions must have been c
                            UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]-UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2]
                            # ADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1]
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
                    UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i] = min(COUNTRY_LLIN_STOCK_ANNUAL_SOY[i],sum(DISTRIBUTION_ANNUAL[i,2:end]))

                    UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,:] = (UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i]/sum(DISTRIBUTION_ANNUAL[i,2:end])).*DISTRIBUTION_ANNUAL[i,2:end]
                end
            end
        end
        
        # Calculate excess and adjust total distributions
        POTENTIAL_EXCESS_STOCK = COUNTRY_LLIN_STOCK_ANNUAL_SOY[i] - UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i]#UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]
        ADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i] = UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i] + ϕ*POTENTIAL_EXCESS_STOCK
        # ADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] = UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] + ϕ*POTENTIAL_EXCESS_STOCK
        # ADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2] = UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2] + ϕ*POTENTIAL_EXCESS_STOCK

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

        # COUNTRY_STOCK_ANNUAL_EOY[i] = COUNTRY_STOCK_ANNUAL_SOY[i]-ADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]
        COUNTRY_LLIN_STOCK_ANNUAL_EOY[i] = COUNTRY_LLIN_STOCK_ANNUAL_SOY[i]-ADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i]
    end

    # Normalise monthly proportions on annual scale
    if isnothing(monthly_p)
        monthly_p = ones(length(MONTHS_MONTHLY)) # Base case assumes regular monthly equal distributions
    end

    effective_monthly_proportions = zeros(length(MONTHS_MONTHLY))
    
    for i in 1:length(YEARS_ANNUAL)
        effective_monthly_proportions[((i-1)*12+1):(i*12)] = monthly_p[((i-1)*12+1):(i*12)]./sum(monthly_p[((i-1)*12+1):(i*12)])
    end

    # Disaggregate annual distributions into monthly distributions with potential seasonality
    DISTRIBUTION_MONTHLY_BYNET = missings(Int64, length(MONTHS_MONTHLY), n_net_types) # Storage variable for distributions
    for i in 1:length(MONTHS_MONTHLY)
        month_i, year_i = monthidx_to_monthyear(i)
        DISTRIBUTION_MONTHLY_BYNET[i,:] = round.(ADJUSTED_DISTRIBUTION_ANNUAL_BYNET[year_i,:].*effective_monthly_proportions[i])
    end

    # Survival net crop matrix calculation
    # A_BYNET[i,j] = Number of nets born in time j that has survived up to time i
    A_BYNET = zeros(length(MONTHS_MONTHLY), length(MONTHS_MONTHLY), n_net_types)

    for n in 1:n_net_types
        for j in 1:length(MONTHS_MONTHLY)
            for i in j:length(MONTHS_MONTHLY)
                if j == i # Monthly distribution
                    A_BYNET[i,j,n] = DISTRIBUTION_MONTHLY_BYNET[i,n]
                else # Calculate from Survival
                    # Probability of a net to survive most recent time step
                    # net_survival_probability = net_loss_weibull((i-j)/12, b_nets[n], k_nets[n])

                    net_survival_probability = net_loss_compact((i-j)/12, b_nets[n], k_nets[n])
                    # net_survival_probability = net_loss_compact((i-j)/12, b_nets[n], κ = k_nets[n])
                    if isnan(net_survival_probability)
                        net_survival_probability = 0
                    end
                    # A_BYNET[i,j,n] = A_BYNET[i-1,j,n]*net_survival_probability # Number from Previous
                    A_BYNET[i,j,n] = A_BYNET[j,j,n]*net_survival_probability # Number from Previous
                end
            end
        end
    end

    Γ_MONTHLY_BYNET = sum(A_BYNET, dims = 2)[:,1,:]
    Γ_MONTHLY_TOTAL = sum(Γ_MONTHLY_BYNET, dims = 2)[:]

    # Regress all on missing values
    regression_idx = findall(.!ismissing.(NET_CROP_MONTHLY))
    for idx in regression_idx
        if (NET_CROP_STD_MONTHLY[idx] == 0)||(isnan(NET_CROP_STD_MONTHLY[idx]))||(isinf(NET_CROP_STD_MONTHLY[idx]))
            continue
        else
            cITN_CROP_MONTHLY[idx] ~ Normal(Γ_MONTHLY_BYNET[idx,1], NET_CROP_STD_MONTHLY[idx])
            NET_CROP_MONTHLY[idx] ~ Normal(Γ_MONTHLY_TOTAL[idx], NET_CROP_STD_MONTHLY[idx])
        end
    end
end


########################################
# %% Forward Simulations
########################################
"""
    model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL; 
                                monthly_p = nothing, 
                                ϕ_est = 0.5, b_est = 2, k_est = 4)
                                
Forward simulation of the net crop model as above. Requires parameter values for ``\\phi``, ``b`` and ``k`` to run. By default, these are set to 0.5, 2 and 4 respectively. If no 'monthly_p' vector is given, a uniform monthly distribution is assumed. Returns a time series as a vector of real numbers of the estimated national net crop.
"""
function model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                ϕ_est, b_net_est, k_net_est, 
                                α_init_est, α_LLIN_est,
                                missing_nets_est;
                                MISSING_NETS_SCALE = MISSING_NETS_SCALE, 
                                monthly_p = nothing, return_age = false,
                                return_stock = false)

    # Get Number of net types from parameter MCMC chain
    n_net_types = length(b_net_est)

    # Get indexes of where distribution data is missing and needs to be imputed from posterior
    missing_dist_idxs = findall(ismissing.(DISTRIBUTION_ANNUAL[:,1]))

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

        if ismissing(DISTRIBUTION_ANNUAL[i,1])
            # Retrieve random variable representing missing distribution
            dist_idx = findfirst(missing_dist_idxs.==i)

            # Assume all distributed nets were LLINs or cITNs depending on starting year (Most conservative assumption)
            if YEARS_ANNUAL[i] < 2010
                # Needed to avoid missings in forward evolution

                if COUNTRY_LLIN_STOCK_ANNUAL_SOY[i] == 0 # No LLIN stock reference (especially in early years)
                    UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] = α_init_est*MISSING_NETS_SCALE
                    # Assign all to cITNs since no LLIN stock
                    UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]
                else
                    UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] = min((1+1/α_LLIN_est)*COUNTRY_LLIN_STOCK_ANNUAL_SOY[i], missing_nets_est[dist_idx]*MISSING_NETS_SCALE)

                    # Assign a fraction of distribution towards cITN
                    UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = ((1/α_LLIN_est)/(1+1/α_LLIN_est))*UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]
                end
                
                
                # Set remaining nets distributed to LLINs
                if n_net_types > 1
                    # UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2:end] .= 0
                    UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2] = UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]-UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1]
                    if n_net_types > 2
                        UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,3:end] .= 0
                    end
                end
            else
                # If missing distribution data for stock and flow, assume dist is equal to delivery
                UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] = min(COUNTRY_LLIN_STOCK_ANNUAL_SOY[i], missing_nets_est[dist_idx]*MISSING_NETS_SCALE)

                # Assign to LLINs (second column)
                UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2] = UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]
                
                # Set remaining nets distributed to 0
                UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = 0
                if n_net_types > 2
                    UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,3:end] .= 0
                end
            end

            # Update LLIN distribution total
            UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i] = sum(UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2:end])
        else
            UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] = DISTRIBUTION_ANNUAL[i,1]#min(COUNTRY_LLIN_STOCK_ANNUAL_SOY[i], DISTRIBUTION_ANNUAL[i,1])
            if UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i] == 0 # No nets distributed in the year
                UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,:] .= 0 # i.e. no distributions
            else
                if YEARS_ANNUAL[i] < 2010 # Need to account for cITNs
                    # Proportion of distributions attributed to cITNs
                    UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = ((1/α_LLIN_est)/(1+1/α_LLIN_est))*UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]
                    # ADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1]
                    
                    # Set remaining nets distributed to LLINs
                    if n_net_types > 1 # i.e there are LLINs
                        # UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2:end] .= 0
                        if (UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]-UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1])>COUNTRY_LLIN_STOCK_ANNUAL_SOY[i]
                            # Number of actually distributed LLINs were based on current stock
                            UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2] = COUNTRY_LLIN_STOCK_ANNUAL_SOY[i]

                            # The rest of the distributions must have been c
                            UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]-UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,2]
                            # ADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1] = UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,1]
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
                    UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i] = min(COUNTRY_LLIN_STOCK_ANNUAL_SOY[i],sum(DISTRIBUTION_ANNUAL[i,2:end]))

                    UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET[i,:] = (UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i]/sum(DISTRIBUTION_ANNUAL[i,2:end])).*DISTRIBUTION_ANNUAL[i,2:end]
                end
            end
        end
        # Calculate excess and adjust total distributions
        POTENTIAL_EXCESS_STOCK = max(COUNTRY_LLIN_STOCK_ANNUAL_SOY[i] - UNADJUSTED_LLIN_DISTRIBUTION_ANNUAL_TOTAL[i],0)#UNADJUSTED_DISTRIBUTION_ANNUAL_TOTAL[i]
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
        monthly_p = ones(length(MONTHS_MONTHLY)) # Base case assumes regular monthly equal distributions
    end

    effective_monthly_proportions = zeros(length(MONTHS_MONTHLY))
    
    for i in 1:length(YEARS_ANNUAL)
        effective_monthly_proportions[((i-1)*12+1):(i*12)] = monthly_p[((i-1)*12+1):(i*12)]./sum(monthly_p[((i-1)*12+1):(i*12)])
    end

    # Disaggregate annual distributions into monthly distributions with potential seasonality
    DISTRIBUTION_MONTHLY_BYNET = missings(Int64, length(MONTHS_MONTHLY), n_net_types) # Storage variable for distributions
    for i in 1:length(MONTHS_MONTHLY)
        month_i, year_i = monthidx_to_monthyear(i)
        DISTRIBUTION_MONTHLY_BYNET[i,:] = round.(ADJUSTED_DISTRIBUTION_ANNUAL_BYNET[year_i,:].*effective_monthly_proportions[i])
    end

    # Survival net crop matrix calculation
    # A_BYNET[i,j] = Number of nets born in time j that has survived up to time i
    A_BYNET = zeros(length(MONTHS_MONTHLY), length(MONTHS_MONTHLY), n_net_types)

    for n in 1:n_net_types
        for j in 1:length(MONTHS_MONTHLY)
            for i in j:length(MONTHS_MONTHLY)
                if j == i # Monthly distribution
                    A_BYNET[i,j,n] = DISTRIBUTION_MONTHLY_BYNET[i,n]
                else # Calculate from Survival
                    # Probability of a net to survive most recent time step
                    # net_survival_probability = net_loss_weibull((i-j)/12, b_net_est[n], k_net_est[n])

                    net_survival_probability = net_loss_compact((i-j)/12, b_net_est[n], k_net_est[n])
                    if isnan(net_survival_probability)
                        net_survival_probability = 0
                    end
                    # A_BYNET[i,j,n] = A_BYNET[i-1,j,n]*net_survival_probability # Number from Previous
                    A_BYNET[i,j,n] = A_BYNET[j,j,n]*net_survival_probability # Number from Previous
                end
            end
        end
    end

    Γ_MONTHLY_BYNET = sum(A_BYNET, dims = 2)[:,1,:]
    
    if return_age
        if return_stock
            return Γ_MONTHLY_BYNET, A_BYNET, COUNTRY_LLIN_STOCK_ANNUAL_EOY, UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET, ADJUSTED_DISTRIBUTION_ANNUAL_BYNET
        else
            return Γ_MONTHLY_BYNET, A_BYNET
        end
    else
        if return_stock
            return Γ_MONTHLY_BYNET, COUNTRY_LLIN_STOCK_ANNUAL_EOY, UNADJUSTED_DISTRIBUTION_ANNUAL_BYNET, ADJUSTED_DISTRIBUTION_ANNUAL_BYNET
        else
            return Γ_MONTHLY_BYNET
        end
    end
end

end