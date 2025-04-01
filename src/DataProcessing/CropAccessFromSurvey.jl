"""
Author: Eugene Tan
Date: 1/10/2024
Updated: 1/10/2024
This module contains a helper function to extract estimates of net crop and
    access purely based on survey observations with no dependence on regression. 

"""
module CropAccessFromSurvey
export calc_survey_estimates


# Import Package dependencies
using NetAccessModel
using Missings

"""
    calc_survey_estimates(input_dict, net_access_input_dict)
- Extracts net crop and net access estimates purely based on survey data only.
- Used for plotting validation.
"""
function calc_survey_estimates(input_dict, net_access_input_dict)
    # Crop Regression input data
    MONTHS_MONTHLY = input_dict["MONTHS_MONTHLY"]
    POPULATION_MONTHLY = input_dict["POPULATION_MONTHLY"]
    NET_CROP_SURVEY_MONTHLY = input_dict["NET_CROP_MONTHLY"]

    NPC_SURVEY_MONTHLY = NET_CROP_SURVEY_MONTHLY./POPULATION_MONTHLY

    # Calculate reference values for National Access from Surveys for scatter plot
    national_H_aggregated = net_access_input_dict["H_aggregated"]
    national_access_aggregated = zeros(size(national_H_aggregated)[1])
    for survey_i in 1:length(national_access_aggregated)
        national_access_aggregated[survey_i] = sum(H_to_access(national_H_aggregated[survey_i,:,:]))
    end

    access_survey_monthidx = net_access_input_dict["survey_monthidx_aggregated"]
    access_survey_monthidx = access_survey_monthidx[findall(.!ismissing.(access_survey_monthidx))]
    NET_ACCESS_SURVEY_MONTHLY = missings(Float64, length(MONTHS_MONTHLY))
    NET_ACCESS_SURVEY_MONTHLY[access_survey_monthidx] .= national_access_aggregated
    
    return NET_CROP_SURVEY_MONTHLY, NPC_SURVEY_MONTHLY, NET_ACCESS_SURVEY_MONTHLY
end

end