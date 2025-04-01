"""
Author: Eugene Tan
Date Created: 28/1/2024
Last Updated: 28/1/2024
Useful data transformatio functions used when processing INLA outputs.
"""

module UsefulTransformations

export piecewise_transform, inv_piecewise_transform
export p_transform, inv_p_transform

"""
    piecewise_transform(x; d = 0.05)

Piecewise transformation functions that is linear in the region x ∈ [0,d], and becomes
    exponential after such that f(x) → 1, x → ∞. Used to model saturation if needed.
    Example use case is for ensuring access is in the correct bounds. 
"""

function piecewise_transform(x; d = 0.05)
    if (0 < x) & (x < (1-d))
        return x
    else
        return (1-d) + d*(1 - exp(-(x-(1-d))/d))
    end
end

"""
    inv_piecewise_transform(y; d = 0.05)
Inverts ```piecewise_transform()```
"""

function inv_piecewise_transform(y; d = 0.05)
    if (0 < y) & (y < (1-d))
        return y
    else
        return -d * log(1 - (y-(1-d))/d) + (1-d)
    end
end

"""
    p_transform(x, μ; ϵ = 1e-3, n = 0.5)

Custom transformation function for an asymmetric logistic function centred at μ.
    Used in INLA transformaiton to fix model. Transformation is given by:

    INSERT EQUATION
"""
function p_transform(x, μ; ϵ = 1e-3, n = 0.5)
    if isnan(x)
        return NaN
    else
        if x < μ
            return -((abs(x-μ)/μ)^n)
        else
            return (abs(x-μ)/(1-μ))^n
        end
    end
end
"""
    inv_p_transform(x, μ; ϵ = 1e-3, n = 0.5)

Inverse of ```inv_p_transform()```
"""
function inv_p_transform(p, μ; n = 0.5)
    if isnan(p)
        return NaN
    else
        if p < 0
            return μ - μ*(-p)^(1/n)
        else
            return μ + (1 - μ)*(p^(1/n))
        end
    end
end

end