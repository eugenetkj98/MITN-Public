"""
Author: Eugene Tan
Date Created: 26/7/2024
Last Updated: 26/7/2024
Helper function for calculating country level normalised Theil coeffecient to
    measure inequality in net distirbution
"""

module Theil
export NPC_to_theil

using LinearAlgebra

"""
    NPC_to_theil(γ_i, N_i)
Calculates the Theil coefficient (measure of inequality in distribution) of nets,
given a collection of local subnational nets per capita γ_i, and population size N_i
"""
function NPC_to_theil(γ_i, N_i)
    # Normalise populations into proportions
    N_i_norm = N_i./sum(N_i)

    # Calculate total net crop
    Γ = dot(γ_i,N_i_norm)

    # Calculate total population
    N = 1#sum(N_i)

    # National NPC
    γ = Γ/N

    # Normalised Fractions NPC
    γ_hat = γ/Γ # National
    γ_hat_i = γ_i./Γ # Subnational

    # # Calculate maximum entropy (most equal distribution)
    # H_0 = -N*γ_hat*log(γ_hat)

    # Calculate actual entropy
    H_i = .-N_i_norm.*γ_hat_i.*log.(γ_hat_i)
    H_i[findall(isnan.(H_i))] .= 0 # Fix regions that have 0 NPC
    H = sum(H_i)

    # Return normalised Theil coefficient
    return -H #1-H/log(N)   
end


end