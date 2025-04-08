"""
Author: Eugene Tan
Date Created: 26/7/2024
Last Updated: 5/8/2024
Code that defined the PPL model for the estimating parameters for Net Access
"""

module NetAccessPrediction
export sample_net_access

# Package requirements
using DataFrames
using LinearAlgebra
using StatsBase
using Distributions
using ProgressBars

# Custom packages
using Logit
using NetAccessModel

########################################
# %% Forward Simulations for net access
########################################

"""
    sample_net_access(ρ_chain_df::DataFrame, μ_chain_df::DataFrame, p_h::Vector,
                            POPULATION_MONTHLY::Vector, Γ_MONTHLY_samples::Matrix;
                            n_max = 20)

This function is used to generate posterior draws of net access time series conditioned on a set of posterior samples of monthly net crop. The inputs are defined as:
- `ρ_chain_df`: A DataFrame of MCMC chain draws for the parameters to calculate ``\\rho_h``
- `μ_chain_df`: A DataFrame of MCMC chain draws for the parameters to calculate ``\\mu_h``
- `p_h`: A vector of probabilities of length ``h_{max}`` for the proportion of people with household size ``h``
- `POPULATION_MONTHLY`: Monthly estimated time series of population. 
- `Γ_MONTHLY_samples`: ``N \\times m`` matrix of net crop posterior time series draws where ``M_{i,j}`` is the net crop of the ``i^{th}`` sample time series in month ``m``.
- `n_max`: Maximum number of a nets one can expect in a give household. Default limit is 20.
"""
function sample_net_access(ρ_chain_df, μ_chain_df, p_h,
                            POPULATION_MONTHLY, Γ_MONTHLY_samples;
                            n_max = 20)
    # Infer number of samples based on number of net crop samples
    n_samples = size(Γ_MONTHLY_samples)[1]

    # Generate randomly indexes to sample posterior at
    sample_idxs = sample(1:size(ρ_chain_df)[1], n_samples)

    # Storage variable for posterior samples of ρ_h and μ_h
    h_max = length(p_h)
    ρ_h_samples = zeros(n_samples, size(Γ_MONTHLY_samples)[2], h_max)
    μ_h_samples = zeros(n_samples, size(Γ_MONTHLY_samples)[2], h_max)

    # Get posterior samples from dataframe chain and arrange in neat format
    β_ρ_samples = ρ_chain_df[sample_idxs,1:6]
    τ_ρ_samples = ρ_chain_df[sample_idxs,7]

    β_μ_samples = μ_chain_df[sample_idxs,1:6]
    τ_μ_samples = μ_chain_df[sample_idxs,7]

    Threads.@threads for i in ProgressBar(1:n_samples, leave = false)
        γ_MONTHLY = Γ_MONTHLY_samples[i,:]./POPULATION_MONTHLY
        for month_idx in 1:size(γ_MONTHLY)[1]
            γ = γ_MONTHLY[month_idx]
            # Sample ρ_h from posterior
            β_vec = Vector(β_ρ_samples[i,:])
            for h in 1:h_max

                # # OLD BV INSPIRED MODEL
                # input_vec = [1/γ, h, h^2, γ, γ^2, γ^3]
                # nu_h = β_vec'*input_vec
                # ρ_h_samples[i,month_idx,h] = inv_emplogit(rand(Normal(nu_h,τ_ρ_samples[i])))

                # NEW MODEL ADJUSTED FOR NPC = 0
                input_vec = [asinh(((1+0.001)/(γ+0.001))-1), h, h^2, γ*h, γ^2, γ^3]
                nu_h = (β_vec'*input_vec) - (τ_ρ_samples[i]^2)/2

                ρ_h_samples[i,month_idx,h] = 2*inv_emplogit(rand(LogNormal(nu_h, τ_ρ_samples[i])))-1
            end
    
            # Sample μ_h from posterior
            for h in 1:h_max
                β_0, β_1, β_2, β_3, β_4, β_5 = β_μ_samples[i,:]
                τ_μ = τ_μ_samples[i]
                μ_mean = β_0 + β_1*(h) + β_2*sqrt(h) + β_3*γ + β_4*(γ^2)+ β_5*(h*γ)
                μ_h_samples[i,month_idx,h] = rand(Normal(μ_mean, τ_μ))
            end
        end
    end

    # Calculate H matrix for survey
    H_samples_MONTHLY = zeros(Float64, n_samples, size(Γ_MONTHLY_samples)[2],h_max,n_max);

    Threads.@threads for i in ProgressBar(1:n_samples, leave = false)
        ρ_h = ρ_h_samples[i,:,:]
        μ_h = μ_h_samples[i,:,:]
        for monthidx in 1:size(Γ_MONTHLY_samples)[2]
            for h in 1:h_max
                for n in 0:n_max-1
                    if n == 0
                        H_samples_MONTHLY[i,monthidx,h,n+1] = p_h[h]*ρ_h[monthidx,h]
                    else
                        H_samples_MONTHLY[i,monthidx,h,n+1] = p_h[h]*(1-ρ_h[monthidx,h])*((μ_h[monthidx,h]^n)/(factorial(big(n))*(exp(μ_h[monthidx,h])-1)))
                        # TEMPORARY FIX FOR NANS
                        if H_samples_MONTHLY[i,monthidx,h,n+1] == NaN
                            H_samples_MONTHLY[i,monthidx,h,n+1] = 0
                        end
                    end
                end
            end
        end
    end

    λ_access_samples_MONTHLY  = zeros(Float64, n_samples,size(Γ_MONTHLY_samples)[2]);
    for i in ProgressBar(1:n_samples, leave = false)
        for monthidx in 1:size(Γ_MONTHLY_samples)[2]
            λ_access_samples_MONTHLY[i,monthidx] = sum(H_to_access(H_samples_MONTHLY[i,monthidx,:,:]))
        end
    end

    return λ_access_samples_MONTHLY
end

end