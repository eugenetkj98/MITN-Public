"""
Author: Eugene Tan
Date Created: 26/7/2024
Last Updated: 26/7/2024
Helper function for empirical logit and inverse empirical logit
"""

module Logit
export emplogit
export inv_emplogit

"""
    emplogit(x; ϵ = 0.001)
Function that calculates the empirical logit transform of x.
"""
function emplogit(x; ϵ = 0.001)
    return log.((x+ϵ)/(1-x+ϵ))
end

"""
    inv_emplogit(x; ϵ = 0.001)
Inverse of `emplogit()`
"""
function inv_emplogit(x; ϵ = 0.001)
    return (exp.(x).*(1-ϵ).-ϵ)./(1 .+ exp.(x))
end

end