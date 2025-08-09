"""
Author: Eugene Tan
Date Created: 8/8/2025
Last Updated: 8/8/2025
Uses extract posterior distribution timeseries to estimate multitype net crop at subnational level
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Prep environment and subdirectories
include(pwd()*"/scripts/read_toml.jl")

###########################################
# %% Load packages
###########################################
# %% Load packages
using JLD2
using CSV
using DataFrames

# %% Other useful packages
using Dates
using ProgressBars
using StatsBase

# %% Custom modules
using DateConversions
using NetCropPrediction

###########################################
# %% Define Directories, paths and parameters
###########################################
# %% Get ISO List
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Define countries where cITNs are forced to have the same decay curve as LLINs due to lack of surveys pre 2010
excl_citn_ISOs = EXCLUSION_CITN_ISOS

# %% Load distribution data
subnat_dist_data = CSV.read(OUTPUT_DATAPREP_DIR*"posterior_distribution_timeseries.csv", DataFrame)

###########################################
# %% Storage variable for Net Crop DataFrame entries
###########################################
netcrop_df_entries = []
age_breakdown_allnets_df_entries = []


###########################################
# %% Simulate for each ISO and subnational region
###########################################

for ISO in filt_ISOs
    println(ISO)
    ###########################################
    # %% Simulate SNF for subnational posterior estimate
    ###########################################

    # Load required JLD2 data files
    nat_extractions = JLD2.load(OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
    nat_posterior_chain = JLD2.load(OUTPUT_REGRESSIONS_DIR*"crop/$(YEAR_NAT_START)_$(YEAR_NAT_END)/$(ISO)_$(YEAR_NAT_START)_$(YEAR_NAT_END)_cropchains.jld2")
    subnat_posterior_chain = JLD2.load(OUTPUT_REGRESSIONS_DIR*"subnational/"*"$(ISO)_SUBNAT_NETCROP_$(YEAR_NAT_START)_$(YEAR_SUBNAT_TRANS)_$(YEAR_NAT_END)_regression.jld2")
    NET_NAMES = names(subnat_raw_dist_data)[findfirst(names(subnat_raw_dist_data) .== "cITN"):end]
    YEARS_ANNUAL = nat_extractions["YEARS_ANNUAL"]
    n_admin1 = length(subnat_posterior_chain["admin1_names"])

    for admin1_i in 1:n_admin1
        # Prepare required distribution input time series and metadata
        DISTRIBUTION_ANNUAL = zeros(length(YEARS_ANNUAL), length(NET_NAMES))
        for year_i in 1:length(YEARS_ANNUAL)
            year = YEARS_ANNUAL[year_i]

            # Metadata
            admin1_entry_metadata[year_i] = subnat_dist_data[(subnat_dist_data.year .== year) .&&
                                                                (subnat_dist_data.admin1 .== subnat_posterior_chain["admin1_names"][admin1_i]) .&&
                                                                (subnat_dist_data.iso .== ISO),["WHO_region","country","iso","admin1","admin1_id","area_id","year"]]

            # Distribution data
            DISTRIBUTION_ANNUAL[year_i,:] .= Matrix(subnat_dist_data[(subnat_dist_data.year .== year) .&&
                                                                (subnat_dist_data.admin1 .== subnat_posterior_chain["admin1_names"][admin1_i]) .&&
                                                                (subnat_dist_data.iso .== ISO),NET_NAMES])[:]
        end

        # Extract posterior estimates of net atrition parameters
        τ_citn_est, τ_llin_est = subnat_posterior_chain["admin1_outputs"][1]["τ_est"]
        κ_citn_est, κ_llin_est = subnat_posterior_chain["admin1_outputs"][1]["κ_est"]

        τ_net_est = Vector{Float64}(undef,length(NET_NAMES))
        κ_net_est = Vector{Float64}(undef,length(NET_NAMES))

        for net_type_i in 1:length(NET_NAMES)
            net_type = NET_NAMES[net_type_i]
            if net_type == "cITN" # cITN Attrition Parameters
                τ_net_est[net_type_i] = τ_citn_est
                κ_net_est[net_type_i] = κ_citn_est
            elseif net_type == "LLIN" # LLIN Attrition Parameters
                τ_net_est[net_type_i] = τ_llin_est
                κ_net_est[net_type_i] = κ_llin_est
            else # Any other net is assumed to follow LLIN attrition
                τ_net_est[net_type_i] = τ_llin_est
                κ_net_est[net_type_i] = κ_llin_est
            end
        end

        # Simulate SNF
        Γ_MONTHLY_BYNET, A_MONTHLY_BYNET = mitn_national_noredist_predict(YEARS_ANNUAL,
                                                                            DISTRIBUTION_ANNUAL,
                                                                            τ_net_est, κ_net_est;
                                                                            monthly_p = nat_posterior_chain["monthly_p"])

        # Prepare admin1 meta_data
        admin1_entry_metadata = Vector{DataFrame}(undef, length(YEARS_ANNUAL))

        for year_i in 1:length(YEARS_ANNUAL)
            year = YEARS_ANNUAL[year_i]

            df_left = subnat_dist_data[(subnat_dist_data.year .== year) .&&
                                                                (subnat_dist_data.admin1 .== subnat_posterior_chain["admin1_names"][admin1_i]) .&&
                                                                (subnat_dist_data.iso .== ISO),["WHO_region","country","iso","admin1","admin1_id","area_id","year"]]
            df_right = DataFrame(month = 1:12)
            admin1_entry_metadata[year_i] = hcat(repeat(df_left,12), df_right)
        end


        ###########################################
        # %% Construct DataFrame Entry for Net Crop
        ###########################################
        netcrop_df_entry = hcat(vcat(admin1_entry_metadata...),DataFrame(sum(Γ_MONTHLY_BYNET, dims = 2), ["Total Nets"]),DataFrame(Γ_MONTHLY_BYNET, NET_NAMES))

        # %% Save entry
        push!(netcrop_df_entries, netcrop_df_entry)

        ###########################################
        # %% Construct DataFrame Entry for Net Age Breakdown
        ###########################################
        # Calculate Age Weight Matrix
        max_age_months = 60
        n_months = size(A_MONTHLY_BYNET)[1]
        n_net_types = length(NET_NAMES)
        net_age_matrix_BYNET = zeros(n_months,max_age_months, n_net_types)

        # Calculate net age demography at monthly resolution
        for n in 1:n_net_types
            for month_i in 1:n_months #current time
                for month_j in 1:month_i # All possible birth times
                    age_months = month_i-month_j
                    age_index = min(age_months + 1, max_age_months)

                    net_age_matrix_BYNET[month_i,age_index,n] += A_MONTHLY_BYNET[month_i,month_j,n]
                end
                # Convert values to proportions
                total_crop = sum(net_age_matrix_BYNET[month_i,:,n])
                if total_crop > 0 
                    net_age_matrix_BYNET[month_i,:,n] .= net_age_matrix_BYNET[month_i,:,n]./total_crop
                end
            end
        end

        # Compile dataframe
        age_breakdown_allnets_df_entry = DataFrame()

        for net_type_i in 1:n_net_types
            age_breakdown_singlenet = DataFrame(net_age_matrix_BYNET[:,:,net_type_i],:auto)
            for j in 1:length(names(age_breakdown_singlenet))
                rename!(age_breakdown_singlenet, names(age_breakdown_singlenet)[j] => string.(0:1:max_age_months)[j])
            end

            age_breakdown_singlenet = hcat(DataFrame(NET_TYPE = repeat([NET_NAMES[net_type_i]], n_months)), age_breakdown_singlenet)

            age_breakdown_allnets_df_entry = vcat(age_breakdown_allnets_df_entry, hcat(vcat(admin1_entry_metadata...), age_breakdown_singlenet))
        end

        # %% Save entry
        push!(age_breakdown_allnets_df_entries, age_breakdown_allnets_df_entry)
    end
end

###########################################
# %% Save Outputs
###########################################
CSV.write(OUTPUT_DIR*"coverage_timeseries/netcrop_multitype_timeseries.csv", vcat(netcrop_df_entries...))
CSV.write(OUTPUT_DIR*"coverage_timeseries/agebreakdown_multitype_timeseries.csv", vcat(age_breakdown_allnets_df_entries...))