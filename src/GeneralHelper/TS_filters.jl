"""
Author: Eugene Tan
Date Created: 8/5/2025
Last Updated: 8/5/2025
Moving average filter implementation
"""
module TS_filters

export MA_filter

using StatsBase

function MA_filter(data; window = 5)
    MA_output = zeros(length(data))
    for i in 1:length(data)
        if i < window
            MA_output[i] = mean(data[1:i])
        else
            MA_output[i] = mean(data[i-window+1:i])
        end
    end

    return MA_output
end


end