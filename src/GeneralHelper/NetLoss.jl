
"""
Author: Eugene Tan
Date: 8/7/2024
Updated: 23/7/2024
This module contains the net loss functions (and inverse) used for the stock and flow model.
"""

module NetLoss

export net_loss_compact,
        net_life_compact,
        net_loss_exp,
        net_life_exp,
        net_loss_weibull,
        net_life_weibull

########################################
# Old loss function used by Bertozzi-Villa
########################################


"""
    net_loss_compact(t::Float64, τ::Float64; κ = 20, precis_threshold = 1e-6)

Compute the survival rate α based on the compact ITN loss function based on Bertozzi-Villa, 2021.

```math
\\mathcal{L} = e^{\\kappa \\left( 1 - \\frac{1}{1-\\left(\\frac{t}{\\tau}\\right)^2} \\right)}
```

Default value of scale parameter is ``\\kappa = 20``, default rounding precision is 1e-6.

# Arguments
- t: the time elapsed in years
- τ: characteristic time scale of the sigmoidal function
"""
function net_loss_compact(t::Union{Int64,Float64}, τ::Union{Int64,Float64}, κ::Union{Int64,Float64}; precis_threshold = 1e-6)
  output = exp(κ.-(κ./(1 .-(t./τ).^2)))
  if (output < precis_threshold) || (t > τ) # Reals overflow/finite preci error check
    output = 0
  end
  return output
end;

"""
    net_life_compact(α::Float64, τ::Float64; κ = 20, precis_threshold = 1e-6)

Inverse function of net_loss used to compute the α survival lifetime for ITN nets.

Default value of scale parameter is κ = 20, default rounding precision is 1e-6.

# Arguments
- α: proportion of nets survived
- τ: characteristic time scale of the sigmoidal function
"""
function net_life_compact(α::Float64, τ::Float64, κ::Float64; precis_threshold = 1e-6)
  output = τ.*sqrt.((-log.(α))./(κ.-log.(α)))
  return output
end;

########################################
# 1 Parameter Exponential Loss FUnction
########################################

"""
    net_loss_exp(t::Float64, τ_max::Float64; κ_squared = 20, precis_threshold = 1e-6)

Compute the survival rate α based on the a modified loss function based on Bertozzi-Villa

```math
\\mathcal{L} = e^{\\kappa^2 \\left( 1 - \\frac{1}{1-\\left(\\frac{t}{\\tau_{\\text{max}}}\\right)^2} \\right)}
```

where ``\\kappa`` is given by the function kappa_coeff and is defined as

```math
\\kappa = \\frac{\\log(\\alpha)}{ 1 - \\frac{1}{ 1 - \\left( \\frac{ \\tau_{\\text{min}} }{ \\tau_{\\text{max}} } \\right)^2 } }
```
"""
function net_loss_exp(t::Float64, τ_max::Float64; κ_squared = 20, precis_threshold = 1e-6)
  output = exp(κ_squared.*(1 .- 1 ./(1-(t ./τ_max) .^2)))
  if (output < precis_threshold) || (output > 10) # Reals overflow/finite preci error check
    output = 0
  end
  return output
end;

"""
    net_life_exp(α::Float64, τ_max::Float64; κ_squared = 20, precis_threshold = 1e-6)

Inverse function for the modified exponential model
"""
function net_life_exp(α::Float64, τ_max::Float64; κ_squared = 20, precis_threshold = 1e-6)
  output = τ_max.*sqrt.((-log.(α))./(κ_squared.-log.(α)))
  return output
end;

"""
    kappa_coeff(α::Float64, τ_min::Float64, τ_max::Float64)

Takes in the shape parameter α∈[0,1] from regressed model and returns the
corresponding exponent coefficient ``\\kappa^2`` to be used in the net loss calculations  
"""
function kappa_coeff(α, τ_min, τ_max)
  return (log.(α)./(1-(1/(1-(τ_min/τ_max)^2)))).^2
end

########################################
# Weibull Loss function
########################################

"""
    net_loss_weibull(t::Float64, b::Float64, k::Float64; precis_threshold = 1e-6)
Two parameter Weibull loss model for net attrition given by
```math
\\mathcal{L} = e^{-\\left( \\frac{t}{b} \\right)^k }
```
"""
function net_loss_weibull(t::Float64, b::Float64, k::Float64; precis_threshold = 1e-6)
  output = exp(.-((t./b).^k))
  if (output < precis_threshold) || (output > 10) # Reals overflow/finite preci error check
    output = 0
  end
  return output
end;

"""
    net_life_weibull(α::Float64, b::Float64, k::Float64)
Inverse function for the Weibull model
"""
function net_life_weibull(α::Float64, b::Float64, k::Float64)
  output = b.*((-log.(α)).^(1 ./k))
  return output
end;

end