"""
Author: Eugene Tan
Date Created: 6/5/2025
Last Updated: 6/5/2025
Numerical implementation of time series convolution. Used to do quick and dirty predictions of net crop.
"""

module Convolutions

export convolution

"""
    convolution(X,Y)

Convolution of two time series X and Y. X is stationary, and Y is sliding.
"""

function convolution(X,Y)
    output_length = length(X) + length(Y)
    output_vals = zeros(output_length)
    Y_slide = zeros(output_length)
    Y_slide[1:length(Y)] = Y
    for i in 1:output_length
        X_slice = zeros(output_length)
        X_slice[1:length(X)] = X

        Y_slice = zeros(output_length)
        Y_slice[1:i] = reverse(Y_slide[1:i])
        output_vals[i] = sum(X_slice.*Y_slice)
    end

    return output_vals
end


end