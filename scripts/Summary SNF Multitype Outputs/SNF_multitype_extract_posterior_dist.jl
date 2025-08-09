"""
Author: Eugene Tan
Date Created: 8/8/2025
Last Updated: 8/8/2025
Script to extract posterior adjusted estimates of net distributions at subnational level, accounts for multiple net type ratio from raw distribution data.
Saves as a .csv file and outputs are used later for posterior estimates on net crop of each each at subnational level
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
using NetCropModel
using NetLoss
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

###########################################
# %% Storage variable for distribution DataFrame entries
###########################################
dist_entries = []

for ISO in filt_ISOs
    println(ISO)
    ###########################################
    # %% Calculate Marginal Posterior Distribution Time Series from National Model (only cITN and LLIN for now)
    ###########################################

    # Load required JLD2 data files
    nat_posterior_chain = JLD2.load(OUTPUT_REGRESSIONS_DIR*"crop/$(YEAR_NAT_START)_$(YEAR_NAT_END)/$(ISO)_$(YEAR_NAT_START)_$(YEAR_NAT_END)_cropchains.jld2")
    nat_extractions = JLD2.load(OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
    YEARS_ANNUAL = nat_extractions["YEARS_ANNUAL"]

    # Generate posterior annual distributions given imputed missing data

    ## Find number of missing distributions entries
    nat_raw_distribution = nat_extractions["DISTRIBUTION_ANNUAL"]
    missingdist_idx = findall(ismissing.(nat_raw_distribution[:,1]))
    n_missing = length(missingdist_idx)

    ## Use imputed missing data posterior data to infer value at missing data points
    missingdist_vals = round.(Int,mean(Matrix(nat_posterior_chain["chain"][:, (end-n_missing+1):end]), dims = 1)[:].*MISSING_NETS_SCALE)

    ## Create imputed (but not delivery adjusted) distribution time series
    nat_imputed_distribution = copy(nat_raw_distribution)

    for i in 1:n_missing
        idx = missingdist_idx[i]

        if nat_extractions["YEARS_ANNUAL"][idx] < 2010 # if pre 2010, then give cITNs, other wise is LLINs
            nat_imputed_distribution[idx,:] .= [missingdist_vals[i], missingdist_vals[i], 0]
        else
            nat_imputed_distribution[idx,:] .= [missingdist_vals[i], 0, missingdist_vals[i]]
        end
    end

    # Get Attrition and conversion posterior parameters
    ϕ_est = mean(nat_posterior_chain["chain"][:,:ϕ])
    α_LLIN_est = mean(nat_posterior_chain["chain"][:,:α_LLIN])
    τ_net_est = mean(Matrix(nat_posterior_chain["chain"][:,[5,7]]), dims = 1)[1,:]
    κ_net_est = mean(Matrix(nat_posterior_chain["chain"][:,[6,8]]), dims = 1)[1,:]

    # Infer final net distributions from national level model
    nat_final_distribution_singletype = mitn_national_predict(nat_extractions["YEARS_ANNUAL"],
                                                                    nat_extractions["DELIVERIES_ANNUAL"], 
                                                                    nat_imputed_distribution[:,2:end],
                                                                    ϕ_est, τ_net_est, κ_net_est, 
                                                                    α_LLIN_est;
                                                                    monthly_p = nat_posterior_chain["monthly_p"], return_distributions = true)[3]
    nat_final_distribution_singletype_annual = vcat([sum(nat_final_distribution_singletype[((i-1)*12+1):i*12,:], dims = 1)[:] for i in 1:size(nat_raw_distribution)[1]]'...)

    ###########################################
    # %% Reconcile national and subnational distribution data and fill in missing values where required and account for multiple net types
    ###########################################
    # Load required data
    subnat_posterior_chain = JLD2.load(OUTPUT_REGRESSIONS_DIR*"subnational/"*"$(ISO)_SUBNAT_NETCROP_$(YEAR_NAT_START)_$(YEAR_SUBNAT_TRANS)_$(YEAR_NAT_END)_regression.jld2")
    subnat_raw_dist_data = CSV.read(RAW_SUBNAT_DATASET_DIR*SUBNAT_MULTITYPE_DISTRIBUTION_DATA_FILENAME, DataFrame)
    NET_NAMES = names(subnat_raw_dist_data)[findfirst(names(subnat_raw_dist_data) .== "cITN"):end]

    # Make temporary storage variable for distribution data
    n_admin1 = length(subnat_posterior_chain["admin1_names"])
    subnat_dist_annual_admin1 = Array{Float64}(undef, n_admin1, size(nat_final_distribution_singletype_annual)[1], size(nat_final_distribution_singletype_annual)[2])

    for admin1_i in ProgressBar(1:n_admin1, leave = false)

        # Calculate population disaggregation ratios for years < YEAR_SUBNAT_TRANS (i.e. pre subnational model)
        admin1_population = subnat_posterior_chain["admin1_outputs"][admin1_i]["FULL_POPULATION_MONTHLY"][1:12:end]
        nat_population = nat_extractions["POPULATION_ANNUAL"][1:end-1]
        pop_disagg_ratios = admin1_population./nat_population

        # Calculate distribution values for <YEAR_SUBNAT_TRANS
        pre_subnat_dist_annual_singletype = nat_final_distribution_singletype_annual[1:(YEAR_SUBNAT_TRANS-YEAR_NAT_START),:].*repeat(pop_disagg_ratios[1:(YEAR_SUBNAT_TRANS-YEAR_NAT_START)],1,size(nat_final_distribution_singletype_annual)[2])

        # Extract distribution values for >=YEAR_SUBNAT_TRANS from subnational adjusted model
        post_subnat_dist_annual_singletype = subnat_posterior_chain["admin1_outputs"][admin1_i]["ADJ_DISTRIBUTION_ANNUAL_BYNET"]

        # Concatenate distribution values into single DataFrame
        subnat_dist_annual_singletype = vcat(pre_subnat_dist_annual_singletype, post_subnat_dist_annual_singletype)

        subnat_dist_annual_multitype = zeros(length(YEARS_ANNUAL), length(NET_NAMES))
        subnat_metadata_annual = Vector{DataFrame}(undef, length(YEARS_ANNUAL))
        for year_i in 1:length(YEARS_ANNUAL)
            year = YEARS_ANNUAL[year_i]

            # Save metadata for year entry
            subnat_metadata_annual[year_i] = subnat_raw_dist_data[(subnat_raw_dist_data.year .== year) .&&
                                                                    (subnat_raw_dist_data.iso .== ISO) .&&
                                                                    (subnat_raw_dist_data.admin1 .== subnat_posterior_chain["admin1_names"][admin1_i]),["WHO_region","country","iso","admin1","admin1_id","area_id","year","type"]]
            # Get data slice from raw distribution values
            subnat_raw_dist_year_entry_multitype = Matrix(subnat_raw_dist_data[(subnat_raw_dist_data.year .== year) .&&
                                                                    (subnat_raw_dist_data.iso .== ISO) .&&
                                                                    (subnat_raw_dist_data.admin1 .== subnat_posterior_chain["admin1_names"][admin1_i]),NET_NAMES])[:]

            # Compare and reconcile distribution values

            if sum(ismissing.(subnat_raw_dist_year_entry_multitype))>0
                subnat_dist_annual_multitype[year_i,1:size(subnat_dist_annual_singletype)[2]] .= subnat_dist_annual_singletype[year_i,:]
            else
                if sum(subnat_raw_dist_year_entry_multitype) == 0 # Raw reported value is apparently 0. Then just take as real distribution value as 0.
                    subnat_dist_annual_multitype[year_i,:] .= (sum(subnat_raw_dist_year_entry_multitype)).*subnat_raw_dist_year_entry_multitype
                else
                    subnat_dist_annual_multitype[year_i,:] .= (sum(subnat_dist_annual_singletype[year_i,:])/sum(subnat_raw_dist_year_entry_multitype)).*subnat_raw_dist_year_entry_multitype
                end
            end
            
        end

        subnat_final_dist_annual_entry = hcat(vcat(subnat_metadata_annual...),DataFrame(sum(subnat_dist_annual_multitype, dims = 2), ["Total Nets"]),DataFrame(subnat_dist_annual_multitype, NET_NAMES))
        push!(dist_entries, subnat_final_dist_annual_entry)
    end
end

# %% Concatenate all entries
subnat_final_dist_annual = vcat(dist_entries...)
CSV.write(OUTPUT_DATAPREP_DIR*"posterior_distribution_timeseries.csv", subnat_final_dist_annual)
CSV.write(OUTPUT_DIR*"coverage_timeseries/posterior_distribution_timeseries.csv", subnat_final_dist_annual) # Put a copy here as well