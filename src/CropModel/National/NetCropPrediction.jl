"""
Author: Eugene Tan
Date Created: 31/3/2025
Last Updated: 31/3/2025
Helper function to generate posterior predictions/draws of NPC given model outputs
"""

module NetCropPrediction
export mitn_national_predict
export mitn_national_noredist_predict

using Missings
using DateConversions
using NetLoss
using NetCropModel


########################################
# %% Prediction using the full National MITN model with stock piling and redistributions
########################################
"""
    mitn_national_predict(YEARS_ANNUAL,
                                DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                ϕ_est, τ_net_est, κ_net_est, 
                                α_LLIN_est;
                                monthly_p = nothing)

Forward posterior prediction of net crop using the MITN model, given parameter values for ϕ, τ, κ, α etc.
"""
function mitn_national_predict(YEARS_ANNUAL,
                                DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                ϕ_est, τ_net_est, κ_net_est, 
                                α_LLIN_est;
                                monthly_p = nothing, return_distributions = false)

    # Get Number of net types from parameter MCMC chain
    n_net_types = size(DISTRIBUTION_ANNUAL)[2]
    
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

    if return_distributions
        return Γ_MONTHLY_BYNET, A_BYNET, DISTRIBUTION_MONTHLY_BYNET
    else
        return Γ_MONTHLY_BYNET, A_BYNET
    end
end

########################################
# %% Prediction a simplified National MITN model with no stockpiling or distributions
########################################
"""
    mitn_national_noredist_predict(YEARS_ANNUAL,
                                DISTRIBUTION_ANNUAL,
                                τ_net_est, κ_net_est;
                                monthly_p = monthly_p)

Predicts the estimated netcrop given some distribution time series. Does not use the redistribution parameter ϕ, or the imputed missing distribution values.
"""
# i.e. only use the raw annual distribution numbers as the driving signal
function mitn_national_noredist_predict(YEARS_ANNUAL,
                                DISTRIBUTION_ANNUAL,
                                τ_net_est, κ_net_est;
                                monthly_p = monthly_p)

    # Get Number of net types from parameter MCMC chain
    n_net_types = size(DISTRIBUTION_ANNUAL)[2]
    
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

end