"""
Author: Eugene Tan
Date Created: 26/7/2024
Last Updated: 5/8/2024
Code that defined the PPL model for the estimating parameters for Net Access
"""

module NetAccessModel

export H_to_access
export model_prop_h_nonets
export model_prop_h_meannets

# Package requirements
using DataFrames
using LinearAlgebra
using StatsBase
using Distributions
using Turing
using ProgressBars

# Custom packages
using Logit

########################################
# %% Helper Functions
########################################

"""
    H_to_access(H::Matrix)

This function takes in a matrix of household counts with ``h_{max}`` rows, and ``n_{max}-1`` columns. Normalises the entries of ``H`` to produce a row stochastic matrix ``\\hat{H}`` of probabilities on a per individual basis, and subsequently uses it to calulate the net access ``\\lambda_{\\text{access}}``.
"""
function H_to_access(H)
    # Get loop bounds
    h_max = size(H)[1]
    n_max = size(H)[2]

    # Calculate H matrix normalised by inidividuals in population
    H_hat = zeros(size(H))
    for h in 1:h_max
        for n in 1:n_max
            H_hat[h,n] = Float64(h*H[h,n]/(h*sum(H)))
        end
    end

    # Calculate net access
    λ_access = zeros(size(H))
    for h in 1:h_max
        for n in 1:n_max
            λ_access[h,n] = min(2*(n-1)/h,1)*H_hat[h,n]
        end
    end

    return λ_access
end

########################################
# %% Statistical model for inferring parameters for distribution of net across demographics
########################################

"""
    model_prop_h_nonets(emplogit_ρ_h_aggregated::Matrix, γ_aggregated::Matrix)

Turing PPL model for modelling the ``\\rho_h``, the proportion of households of size ``h`` that have no nets available.  ``emplogit_ρ_h_aggregated`` is a matrix of ``h_{max}`` where each row corresponds to a single DHS survey. `γ_aggregated` is a vector of net per capita values ∈ (0,1) corresponding to the survey of each row in `emplogit_ρ_h_aggregated`. The regression is modelled with a hierarchical model based on the Bertozzi-Villa and Bhatts model given as follows:
```math
\\begin{aligned}
    \\beta_0 &\\sim \\text{Uniform}(-50,50)\\\\
    \\beta_1 &\\sim \\text{Uniform}(-3,3)\\\\
    \\beta_2 &\\sim \\text{Uniform}(-1,1)\\\\
    \\beta_3 &\\sim \\text{Uniform}(-100,100)\\\\
    \\beta_4 &\\sim \\text{Uniform}(-300,300)\\\\
    \\beta_5 &\\sim \\text{Uniform}(-300,300)\\\\
    \\tau_\\rho &\\sim \\text{Gamma}(0.1,0.1)\\\\
    \\nu_h &= \\beta_0 + \\beta_1 h + \\beta_2 h^2 + \\beta_3 \\gamma + \\beta_4 \\gamma^2 + \\beta_5 \\gamma^3\\\\
    \\text{emplogit}(\\rho_h) &\\sim \\text{Normal}(\\nu_h, \\tau_\\rho)
\\end{aligned}
```
"""
@model function model_prop_h_nonets(emplogit_ρ_h_aggregated, γ_aggregated)
    # Get index bounds from data
    n_surveys = size(emplogit_ρ_h_aggregated)[1]
    h_max = size(emplogit_ρ_h_aggregated)[2]
    
    # # OLD BV INSPIRED MODEL Priors
    # β_0 ~ Uniform(-50,50)
    # β_1 ~ Uniform(-3,3)
    # β_2 ~ Uniform(-1,1)
    # β_3 ~ Uniform(-100,100)
    # β_4 ~ Uniform(-300,300)
    # β_5 ~ Uniform(-300,300)
    # τ_ρ ~ Uniform(0,0.5)

    # Priors # NEW MODEL ADJUSTED FOR NPC = 0
    β_0 ~ Uniform(-50,50)
    β_1 ~ Uniform(-5,5)
    β_2 ~ Uniform(-1,1)
    β_3 ~ Uniform(-100,100)
    β_4 ~ Uniform(-300,300)
    β_5 ~ Uniform(-500,500)
    τ_ρ ~ Gamma(0.1,0.1)
    
    # Intermediate variables
    for h in 1:h_max
        for i in 1:n_surveys
            γ = γ_aggregated[i]
            # nu_h = β_0*(1/γ) + β_1*h + β_2*(h^2) + β_3*γ + β_4*(γ^2) + β_5*(γ^3)
            # nu_h = (β_0*(γ-0.5) + β_1*h*(γ) + β_2*(h^2)*(γ) + β_3*((γ)^2) + β_4*((γ)^3))
            # emplogit_ρ_h_aggregated[i,h] ~ Normal(nu_h, τ_ρ)

            # NEW MODEL ADJUSTED FOR NPC = 0
            nu_h = β_0*asinh(((1+0.001)/(γ+0.001))-1) + β_1*h + β_2*(h^2) + β_3*γ*h + β_4*(γ^2) + β_5*(γ^3)-((τ_ρ)^2)/2
            emplogit_ρ_h_aggregated[i,h] ~ LogNormal(nu_h, τ_ρ)
        end
    end
end

"""
    model_prop_h_meannets(μ_h_aggregated::Matrix, γ_aggregated::Vector)

Turing PPL model for modelling the ``\\mu_h``, the mean number of net in households of size ``h``. Similar to the Bhatt paper, the mean is subsequently used to model entries of the probability matrix ``H``  as a zero truncated Poisson distribution. Inputs are matrix of aggregated empirical survey values.  `μ_h_aggregated` is a matrix of ``h_{max}`` where each row corresponds to a single DHS survey. `γ_aggregated` is a vector of net per capita values ∈ (0,1) corresponding to the survey of each row in `μ_h_aggregated`.
```math
\\begin{aligned}
    \\beta_0 &\\sim \\text{Uniform}(-50,50)\\\\
    \\beta_1 &\\sim \\text{Uniform}(-3,3)\\\\
    \\beta_2 &\\sim \\text{Uniform}(-10,10)\\\\
    \\beta_3 &\\sim \\text{Uniform}(-50,50)\\\\
    \\beta_4 &\\sim \\text{Uniform}(-50,50)\\\\
    \\beta_5 &\\sim \\text{Uniform}(-50,50)\\\\
    \\tau_\\mu &\\sim \\text{Gamma}(0.1,0.1)\\\\
    \\bar{\\mu} &= \\beta_0 + \\beta_1 h + \\beta_2 \\sqrt{h} + \\beta_3 \\gamma + \\beta_4 \\gamma^2 + \\beta_5 h \\gamma\\\\
    \\mu_h &\\sim \\text{Normal}(\\bar{\\mu}, \\tau_\\mu)
\\end{aligned}
```
"""
@model function model_prop_h_meannets(μ_h_aggregated, γ_aggregated)
    # Get index bounds from data
    n_surveys = size(μ_h_aggregated)[1]
    h_max = size(μ_h_aggregated)[2]
    
    # Priors
    β_0 ~ Uniform(-50,50)
    β_1 ~ Uniform(-3,3)
    β_2 ~ Uniform(-10,10)
    β_3 ~ Uniform(-50,50)
    β_4 ~ Uniform(-50,50)
    β_5 ~ Uniform(-50,50)
     
    τ_μ ~ Gamma(0.1,0.1)

    # Intermediate variables
    for h in 1:h_max
        for i in 1:n_surveys
            γ = γ_aggregated[i]
            # μ_mean = γ*(β_0 + β_1*(h) + β_2*(h^2) + β_3*sqrt(h))
            μ_mean = β_0 + β_1*(h) + β_2*sqrt(h) + β_3*γ + β_4*(γ^2)+ β_5*(h*γ)
            μ_h_aggregated[i,h] ~ Normal(μ_mean, τ_μ)
        end
    end
end

end