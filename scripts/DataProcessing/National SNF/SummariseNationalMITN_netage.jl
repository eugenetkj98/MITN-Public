"""
Author: Eugene Tan
Date Created: 4/11/2024
Last Updated: 12/11/2024
Extract posterior net crop and attrition estimates from regressed national crop model to interface with subnational model.
"""
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Import Public Packages
using JLD2
using CSV
using DataFrames
using DateConversions
using ProgressBars
using LinearAlgebra
using StatsBase
using KernelDensity
using NetCropModel
using DateConversions

# %% Extraction Settings
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END
n_samples = NAT_CROPAGE_N_DRAWS
max_age_months = 4*12

# %% Define Posterior Models data location
netcrop_extractions_dir = OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/"
netcrop_regressions_dir = OUTPUT_REGRESSIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/"

# %% Define save directory
save_dir = OUTPUT_DRAWS_DIR*"national/demography/"
mkpath(save_dir)
# %%
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Extract net attrition parameters for all countries
for ISO in filt_ISOs

    println("Generating net age draws for $ISO...")
    # Get Net Attrition Posteriors MCMC chain
    input_filename = ISO*"_$(YEAR_START)_$(YEAR_END)_cropextract.jld2"
    regression_filename = ISO*"_$(YEAR_START)_$(YEAR_END)_cropchains.jld2"
    input_dict = load(netcrop_extractions_dir*input_filename)
    regression_dict = load(netcrop_regressions_dir*regression_filename)

    # Crop Regression input data
    ISO = input_dict["ISO"]
    YEARS_ANNUAL = input_dict["YEARS_ANNUAL"]
    MONTHS_MONTHLY = input_dict["MONTHS_MONTHLY"]
    DELIVERIES_ANNUAL = input_dict["DELIVERIES_ANNUAL"]
    DISTRIBUTION_ANNUAL = input_dict["DISTRIBUTION_ANNUAL"]
    NET_CROP_MONTHLY = input_dict["NET_CROP_MONTHLY"]
    POPULATION_MONTHLY = input_dict["POPULATION_MONTHLY"]
    NET_NAMES = input_dict["NET_NAMES"]
    n_net_types = length(NET_NAMES)

    # Get number of missing net values
    n_missing_nets_vals = sum(ismissing.(DISTRIBUTION_ANNUAL[:,1]))

    # Crop Regression output data
    chain = regression_dict["chain"]
    monthly_p = regression_dict["monthly_p"]
    chain_UNIF = regression_dict["chain_epochs"][1]
    monthly_p_UNIF = regression_dict["monthly_weights"][1]

    ##### Extract Net Crop MCMC parameters

    ### Part 1: Final Regressed Values
    # Randomly generate sample indexes to sample from chain
    chain_length = size(chain)[1]
    sample_idxs = sample(1:chain_length, n_samples, replace = false)

    # Extract MCMC draws for parameters
    ϕ_posterior_draws = chain[sample_idxs, :ϕ]
    α_init_posterior_draws = chain[sample_idxs,:α_init]
    α_LLIN_posterior_draws = chain[sample_idxs,:α_LLIN]

    b_net_posterior_draws = Matrix(DataFrame(chain)[:,5:2:5+2*(n_net_types-1)])[sample_idxs, :]
    k_net_posterior_draws = Matrix(DataFrame(chain)[:,6:2:6+2*(n_net_types-1)])[sample_idxs, :]
    if n_missing_nets_vals > 0
        n_missing_nets_posterior_draws = Matrix(DataFrame(chain)[:,(4+2*(n_net_types)+1):end])[sample_idxs, :]
    else
        # No missing data, so just feed in a zero matrix with no columns
        n_missing_nets_posterior_draws = zeros(size(chain)[1], 0)
    end

    # Calculate month, year value indices
    # Do age conversions from monthidx to month year format
    year_values = zeros(Int64, length(MONTHS_MONTHLY))
    month_values = zeros(Int64, length(MONTHS_MONTHLY))
    for i in 1:length(MONTHS_MONTHLY)
        month, relative_year = monthidx_to_monthyear(i)
        year = YEAR_START + relative_year - 1
        year_values[i] = year
        month_values[i] = month
    end

    ##### Generate Net Crop Trajectories
    Γ_MONTHLY_samples_BYNET = zeros(n_samples, length(MONTHS_MONTHLY), n_net_types)
    A_samples_BYNET = zeros(n_samples, length(MONTHS_MONTHLY), length(MONTHS_MONTHLY), n_net_types)
    net_age_matrix_samples_BYNET = zeros(n_samples, length(MONTHS_MONTHLY),max_age_months, n_net_types)
    age_breakdown_allnets_dataframerows = Vector{Any}(undef, n_samples)
    

    Threads.@threads for i in ProgressBar(1:n_samples, leave = false)
        # Select MCMC posterior draw parameters
        ϕ = ϕ_posterior_draws[i]
        α_init = α_init_posterior_draws[i]
        α_LLIN = α_LLIN_posterior_draws[i]
        b_nets = b_net_posterior_draws[i,:]
        k_nets = k_net_posterior_draws[i,:]
        missing_nets = n_missing_nets_posterior_draws[i,:]

        Γ_MONTHLY_samples_BYNET[i,:,:,:], A_samples_BYNET[i,:,:,:] = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                                                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                                                                    ϕ, b_nets, k_nets,
                                                                                    α_init, α_LLIN,
                                                                                    missing_nets; 
                                                                                    monthly_p = monthly_p,
                                                                                    return_age = true)
        

        monthyear = DataFrame(ISO = ISO, YEAR = year_values, MONTH = month_values, sample_idx = i)

        # Calculate net age demography at monthly resolution
        for month_i in 1:length(MONTHS_MONTHLY) #current time
            for month_j in 1:month_i # All possible birth times
                for n in 1:n_net_types
                    age_months = month_i-month_j
                    age_index = min(age_months + 1, max_age_months)

                    net_age_matrix_samples_BYNET[i,month_i,age_index,n] += A_samples_BYNET[i,month_i,month_j,n]
                end
            end

            net_age_matrix_samples_BYNET[i,month_i,:,:] = net_age_matrix_samples_BYNET[i,month_i,:,:]./POPULATION_MONTHLY[month_i]
        end

        # Compile dataframe
        age_breakdown_allnets_sample = DataFrame()

        for net_type_i in 1:n_net_types
            age_breakdown_singlenet = DataFrame(net_age_matrix_samples_BYNET[i,:,:,net_type_i],:auto)
            for j in 1:length(names(age_breakdown_singlenet))
                rename!(age_breakdown_singlenet, names(age_breakdown_singlenet)[j] => string.(0:1:max_age_months)[j])
            end

            age_breakdown_singlenet = hcat(DataFrame(NET_TYPE = repeat([NET_NAMES[net_type_i]], length(MONTHS_MONTHLY))), age_breakdown_singlenet)

            age_breakdown_allnets_sample = vcat(age_breakdown_allnets_sample, hcat(monthyear, age_breakdown_singlenet))
        end

        age_breakdown_allnets_dataframerows[i] = age_breakdown_allnets_sample
        
    end

    # Construct compiled DataFrame
    age_breakdown_allnets = vcat(age_breakdown_allnets_dataframerows...)

    # Calculate average age demography and write into DataFrame
    net_age_matrix_mean_BYNET = mean(net_age_matrix_samples_BYNET, dims = 1)[1,:,:,:]
    age_breakdown_allnets_means = DataFrame()

    monthyear = DataFrame(ISO = ISO, YEAR = year_values, MONTH = month_values, sample_idx = "MEANS")
    for net_type_i in 1:n_net_types
        age_breakdown_singlenet = DataFrame(net_age_matrix_mean_BYNET[:,:,net_type_i],:auto)
        for j in 1:length(names(age_breakdown_singlenet))
            rename!(age_breakdown_singlenet, names(age_breakdown_singlenet)[j] => string.(0:1:max_age_months)[j])
        end

        age_breakdown_singlenet = hcat(DataFrame(NET_TYPE = repeat([NET_NAMES[net_type_i]], length(MONTHS_MONTHLY))), age_breakdown_singlenet)

        age_breakdown_allnets_means = vcat(age_breakdown_allnets_means, hcat(monthyear, age_breakdown_singlenet))
    end
    
    # Save all dataframes as CSVs
    CSV.write(save_dir*"$(ISO)_net_age_demography_samples.csv", age_breakdown_allnets)
    CSV.write(save_dir*"$(ISO)_net_age_demography_mean.csv", age_breakdown_allnets_means)
end



