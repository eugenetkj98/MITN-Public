"""
Author: Eugene Tan
Date Created: 26/7/2024
Last Updated: 5/78/2024
Code that defined the PPL model for the estimating parameters for Net Access
"""

module NetAccessRegression
export bayes_access

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Default parameters for net access regression
# Bayesian MCMC hyperparameters
iterations = NAT_ACCESS_MCMC_ITERATIONS; # Number of MCMC iterations
n_chains = min(NAT_ACCESS_N_CHAINS, Threads.nthreads()); # Number of chains
# n_params = 5
burn_in = NAT_ACCESS_MCMC_BURNIN; # Number of iterations to discard for burn-in

using Logit # Import Logit helper function
using NetAccessModel
using Turing
using DataFrames
using JLD2

# Save directory for final chains
chain_output_dir = OUTPUT_REGRESSIONS_DIR*"access/"

########################################
# %% Function to perform model fit for Net Crop on survey data
########################################

# """
#     bayes_access(input_data_dict::Dict)
# Need to write documentation
# """
# function bayes_access(ISO::String, extractions_dict,
#                         net_access_input_dict;
#                         n_max = 20, h_max = 10,
#                         iterations = 10000, n_params = 2000,
#                         save_output = true, chain_output_dir = chain_output_dir)

#     # Extract data from dicts
#     YEAR_START = extractions_dict["YEAR_START"]
#     YEAR_END = extractions_dict["YEAR_END"]

#     POPULATION_MONTHLY = net_access_input_dict["POPULATION_MONTHLY"]
#     ρ_h_aggregated = net_access_input_dict["ρ_h_aggregated"]
#     μ_h_aggregated = net_access_input_dict["μ_h_aggregated"]
#     γ_aggregated = net_access_input_dict["γ_aggregated"]

#     H_aggregated = net_access_input_dict["H_aggregated"]
#     p_h_aggregated = net_access_input_dict["p_h_aggregated"]

#     ### Part 1: Regressing for ρ_h                         
#     # Calculate empirical logit of ρ_h to prep for regression as per model
#     emplogit_ρ_h_aggregated = emplogit.(ρ_h_aggregated)

#     # Bayesian inference for model parameters with Turing.jl
#     model_ρ = model_prop_h_nonets(emplogit_ρ_h_aggregated, γ_aggregated)
#     ρ_h_chain = sample(model_ρ, NUTS(), iterations; discard_initial = burn_in)

#     ### Part 2: Regressing for μ_h
#     model_μ = model_prop_h_meannets(μ_h_aggregated, γ_aggregated)
#     μ_h_chain = sample(model_μ, NUTS(), iterations; discard_initial = burn_in)

#     ### Part 3: Package results to dict and return as function output
#     ρ_chain_df = DataFrame(ρ_h_chain[:,1:7,:])[:,3:end]
#     μ_chain_df = DataFrame(μ_h_chain[:,:,:])[:,3:end]

#     netaccess_dict = Dict("ρ_chain_df" => ρ_chain_df,
#                         "μ_chain_df" => μ_chain_df,
#                         "ISO" => ISO,
#                         "YEAR_START" => YEAR_START,
#                         "YEAR_END" => YEAR_END,
#                         "POPULATION_MONTHLY" => POPULATION_MONTHLY,
#                         "H_aggregated" => H_aggregated,
#                         "p_h_aggregated" => p_h_aggregated)

#     # Save last epoch's chain result in directory
#     if save_output
#         # Define filename
#         chain_output_filepath = chain_output_dir*ISO*"_"*"$(YEAR_START)"*"_"*"$(YEAR_END)"*"_accesschains.jld2"
    
#         # Save .jld2 file
#         save(chain_output_filepath, netaccess_dict)

#         println("File saved at: "*chain_output_filepath)
#     end

#     return netaccess_dict
# end

"""
    bayes_access(access_survey_globaldata::Dict; iterations = 10000, burn_in = 2000, chain_output_dir = chain_output_di)
Function to perform MCMC draws of access model parameters `ρ_h` and `μ_h` using empirical survey data. Input variable `access_survey_globaldata` is a `Dict` with the following key values:
- `p_h_globaldata`: Matrix of aggregated data for distribution of households of size ``h`` for each survey.
- `ρ_h_globaldata`: Matrix of aggregated data of proportion of households of size ``h`` with no nets for each survey.
- `μ_h_globaldata`: Matrix of aggregated data of mean number of nets in households of size ``h`` given at least one net.

Does not return a result, but instead saves the posterior MCMC chain as DataFrames with the filename `netaccesschains.jld2` in the directory indicated by `chain_output_dir`.
"""
function bayes_access(access_survey_globaldata;
                            iterations = 10000, burn_in = 2000,
                            chain_output_dir = chain_output_dir, n_chains = n_chains)
    # Read saved aggregate data from surveys for regressing net access model
    p_h_globaldata = access_survey_globaldata["p_h_globaldata"]
    ρ_h_globaldata = access_survey_globaldata["ρ_h_globaldata"]
    μ_h_globaldata = access_survey_globaldata["μ_h_globaldata"]
    γ_globaldata = access_survey_globaldata["γ_globaldata"]

    # Pre-process data to be in required form for model
    # emplogit_ρ_h_globaldata = emplogit.(ρ_h_globaldata) # FOR OLD BV INSPIRED MODEL
    emplogit_ρ_h_globaldata = emplogit.((ρ_h_globaldata .+ 1)./2) # FOR NEW ACCESS MODEL

    # Construct PPL models and do MCMC regression
    model_ρ = model_prop_h_nonets(emplogit_ρ_h_globaldata, γ_globaldata)
    model_μ = model_prop_h_meannets(μ_h_globaldata, γ_globaldata)

    ρ_h_chain = []
    μ_h_chain = []

    if n_chains > 1
        ρ_h_chain = sample(model_ρ, NUTS(), MCMCThreads(), iterations, n_chains; progress = true, discard_initial = burn_in)   
        μ_h_chain = sample(model_μ, NUTS(), MCMCThreads(), iterations, n_chains; progress = true, discard_initial = burn_in)
    else
        ρ_h_chain = sample(model_ρ, NUTS(), iterations; discard_initial = burn_in)   
        μ_h_chain = sample(model_μ, NUTS(), iterations; discard_initial = burn_in)
    end
    

    # Convert MCMC Chain to DataFrame for storage
    ρ_chain_df = DataFrame(ρ_h_chain[:,1:7,:])[:,3:end]
    μ_chain_df = DataFrame(μ_h_chain[:,1:7,:])[:,3:end]

    # Save data

    # Make filepath if it doesn't already exist
    mkpath(chain_output_dir)
    
    chain_output_filepath = chain_output_dir*"netaccesschains.jld2"
    jldsave(chain_output_filepath;
                ρ_chain_df, μ_chain_df)
    println("File saved at: "*chain_output_filepath)
end

end