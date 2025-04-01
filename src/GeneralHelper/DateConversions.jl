"""
Author: Eugene Tan
Date Created: 23/7/2024
Last Updated: 23/7/2024
Helper functions to convert from month/year format to month indexes.
"""

module DateConversions

export monthyear_to_monthidx, monthidx_to_monthyear

# Default Constants

YEAR_START = 2010

"""
    monthyear_to_monthidx(month::Int, year::Int)

Function to quickly convert month year format to generalised indexed month format.
Start bounds are taken from YEAR_START
"""
function monthyear_to_monthidx(month::Int, year::Int; YEAR_START = YEAR_START)
    if year < YEAR_START
        throw("Invalid year value")
    end
    if month > 12
        throw("Invalid month value")
    end
    return (year % YEAR_START)*12 + month
end

"""
    monthidx_to_monthyear(monthidx::Int)

Inverse function of monthyear_to_monthidx() but returns year as an Int starting
from 1.
"""
function monthidx_to_monthyear(monthidx::Int)
    month = monthidx % 12
    if month == 0
        month = 12
    end
    return month, (((monthidx-1)÷12)+1)
end

end