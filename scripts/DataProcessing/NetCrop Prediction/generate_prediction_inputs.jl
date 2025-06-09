"""
Author: Eugene Tan
Date Created: 15/4/2025
Last Updated: 15/4/2025
Construct prediction inputs for running forward/backcasting estimates with multiple net types. Dataset is reconstructed from existing distribution data.
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/read_toml.jl")

# %% Load packages
using JLD2
using CSV
using DataFrames

# %% Define filepaths
deliveries_filepath = RAW_DATASET_DIR*DELIVERIES_DATA_FILENAME
distributions_filepath = RAW_DATASET_DIR*MULTITYPE_DISTRIBUTION_DATA_FILENAME

# %% Output filepath
output_dir = OUTPUT_FORWARD_PRED_DIR
mkpath(output_dir)


# %% Import CSV files
deliveries_data = CSV.read(deliveries_filepath, DataFrame)
distribution_data = CSV.read(distributions_filepath, DataFrame)

# Replace all missing distribution data entries with 0 nets
function replace_df_missings(x)
    if (ismissing(x) || isnan(x))
        return 0
    else
        return x
    end
end
distribution_data[:, 6:end] .= replace_df_missings.(distribution_data[:, 6:end])

# %% Get ISO List
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = ["CPV","ZAF"]
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Year range
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END
YEAR_VALS = YEAR_START:YEAR_END

# %% Get Net Type Metadata
NET_NAMES = names(distribution_data)[6:end]

# %% Construct dataset for each country
for ISO in filt_ISOs
    net_distributions = zeros(length(YEAR_VALS), length(NET_NAMES))
    llin_deliveries = zeros(length(YEAR_VALS))

    for year_i in 1:length(YEAR_VALS)
        year = YEAR_VALS[year_i]
        llin_deliveries[year_i] = Int64(deliveries_data[deliveries_data.iso3 .== ISO,"$year"][1])
        for net_i in 1:length(NET_NAMES)
            net_distributions[year_i, net_i] = Int64(distribution_data[(distribution_data.iso .== ISO) .&&
                                                            (distribution_data.year .== year),NET_NAMES[net_i]][1])
        end
    end

    years_df = DataFrame(:year => Array(collect(YEAR_VALS)))
    llin_deliveries_df = DataFrame(:LLIN_delivery => llin_deliveries)
    net_distributions_df = DataFrame(net_distributions, NET_NAMES)
    total_distributions_df = DataFrame(sum(Array(net_distributions_df), dims = 2), ["total_distribution"])
    output_df = hcat(years_df, llin_deliveries_df, total_distributions_df, net_distributions_df)

    # Save dataset
    CSV.write(output_dir*"$(ISO)_prediction_input.csv", output_df)
end