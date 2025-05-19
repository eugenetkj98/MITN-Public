"""
Author: Eugene Tan
Date Created: 13/8/2024
Last Updated: 5/9/2024
Code to execute iterative algorithm for fitting the national net crop model to survey data. This version accounts for multiple net types and must be used alongside NetCropModel.jl
"""

module NetCropRegression
export bayes_GD
export normalised_monthly_weights

# %% Import filenames and directories from TOML file
include(pwd()*"/scripts/read_toml.jl")

using ProgressBars
using LinearAlgebra
using Distributions
using Turing
using AdvancedMH
using JLD2
using DataFrames
using Missings

using NetCropModel
using NetLoss

# %% Default parameters for regression
# Bayesian MCMC hyperparameters
iterations = NAT_CROP_MCMC_ITERATIONS; # Number of MCMC iterations
n_chains = min(NAT_CROP_N_CHAINS, Threads.nthreads()); # Number of chains
burn_in = NAT_CROP_MCMC_BURNIN#10000; # Number of iterations to discard for burn-in
# proposal_resampling_variances = 0.004 .*[30,15,20,20,4,5] # Weibull sampling, last two entries are for b and k
proposal_resampling_variances = NAT_CROP_PROPOSAL_SAMPLING_VAR # Compact support sampling

# SGD hyperparameters
N_EPOCHS = NAT_CROP_SGD_EPOCHS # Maximum number of SGD alternating iterations
EPOCH_LEN = NAT_CROP_SGD_STEPS # Number of SGD steps in each iteration
NORMALISE_PERIOD = 1
α = NAT_CROP_SGD_ALPHA
ϵ = NAT_CROP_SGD_EPSILON 

# DIC Loss function
DIC_SAMPLE_SIZE = NAT_CROP_SGD_DIC_SAMPLESIZE

# Save directory for final chains
chain_output_dir = OUTPUT_REGRESSIONS_DIR*"crop/"

########################################
# %% Helper function to normalise monthly distribution weights
########################################

"""
    normalised_monthly_weights(monthly_p)
This function just normalises the monthly_p weights to represent a proportion out of 1, of the annual allocated distributions, that are distributed in a given month.
"""
function normalised_monthly_weights(monthly_p)
    normalised_monthly_p = zeros(length(monthly_p))
    
    for i in 1:(length(monthly_p)÷12)
            normalised_monthly_p[((i-1)*12+1):(i*12)] = monthly_p[((i-1)*12+1):(i*12)]./sum(monthly_p[((i-1)*12+1):(i*12)])
    end

    return normalised_monthly_p
end

########################################
# %% Helper function to calculate DIC
########################################

function calc_deviance(Y, Y_STD, θ_vals)#, σ)
    D_i = missings(Float64, length(Y))
    for i in 1:length(Y)
        # if θ_vals[i] == 0
        if Y_STD[i] == 0
            # D_i[i] = -2*logpdf(Normal(θ_vals[i], σ*abs(minimum(θ_vals[findall(θ_vals.>0)]))), Y[i])
            D_i[i] = -2*logpdf(Normal(θ_vals[i], mean(Y_STD)), Y[i])
        else
            # D_i[i] = -2*logpdf(Normal(θ_vals[i], σ*abs(θ_vals[i])), Y[i])
            D_i[i] = -2*logpdf(Normal(θ_vals[i], Y_STD[i]), Y[i])
        end
    end
    D = mean(D_i)
    return D
end

function DIC_netcrop(YEARS_ANNUAL, MONTHS_MONTHLY,
                        DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL, 
                        NET_CROP_MONTHLY, NET_CROP_STD_MONTHLY,
                        target_idx, monthly_p,
                        ϕ_vals, b_net_vals, k_net_vals, 
                        α_init_vals, α_LLIN_vals,
                        missing_nets_vals; λ = 1,
                        full_mode = false, printout = false)
    # Calculate Laplace regularisation

    # Construct Toeplitz matrix for finite difference
    A = zeros(length(monthly_p), length(monthly_p))
    for i in 1:length(monthly_p)
        if i == 1
            A[i,1:3] .= [1,-2,1]
        elseif i == length(monthly_p)
            A[i,(i-2):i] .= [1,-2,1]
        else
            A[i,(i-1):(i+1)] .= [1,-2,1]
        end
    end


    # Calculate Laplace penalty term
    L_penalty = λ*sum((A*monthly_p).^2)

    if full_mode # Calculate the full DIC across the entire joint posterior
        D_iter = missings(Float64, size(ϕ_vals)[1])
        
        for iter in 1:size(ϕ_vals)[1]
            Γ_posterior_draw_BYNET = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                            DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                            ϕ_vals[iter],
                                            b_net_vals[iter,:],
                                            k_net_vals[iter,:],
                                            α_init_vals[iter],
                                            α_LLIN_vals[iter],
                                            missing_nets_vals[iter,:]; 
                                            monthly_p = monthly_p)
            Γ_posterior_draw_TOTAL = sum(Γ_posterior_draw_BYNET, dims = 2)[:]

            D_iter[iter] = calc_deviance(NET_CROP_MONTHLY[target_idx], NET_CROP_STD_MONTHLY[target_idx],
                                        Γ_posterior_draw_TOTAL[target_idx])
        end

        DIC = mean(D_iter)

        if printout
            println("Total Loss = $(round(DIC + L_penalty, digits = 3)) DIC = $(round(DIC, digits = 3)), L_penalty = $(round(L_penalty, digits = 3))")
        end
        return DIC + L_penalty
    else # Calculate only using the mean of the joint posterior
        ϕ_est = mean(ϕ_vals)
        α_init = mean(α_init_vals)
        α_LLIN = mean(α_LLIN_vals)
        b_net_est = mean(b_net_vals, dims = 1)[:]
        k_net_est = mean(k_net_vals, dims = 1)[:]
        missing_nets_est = mean(missing_nets_vals, dims = 1)[:]
        

        Γ_posterior_draw_BYNET = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                            DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                            ϕ_est,
                                            b_net_est,
                                            k_net_est,
                                            α_init,
                                            α_LLIN,
                                            missing_nets_est; 
                                            monthly_p = monthly_p)

        Γ_posterior_draw_TOTAL = sum(Γ_posterior_draw_BYNET, dims = 2)[:]

        DIC = calc_deviance(NET_CROP_MONTHLY[target_idx], NET_CROP_STD_MONTHLY[target_idx],
                            Γ_posterior_draw_TOTAL[target_idx])

        if printout
            println("Total Loss = $(round(DIC + L_penalty, digits = 3)) DIC = $(round(DIC, digits = 3)), L_penalty = $(round(L_penalty, digits = 3))")
        end
        return DIC + L_penalty
    end
end

########################################
# %% Function to perform model fit for Net Crop on survey data
########################################
"""
    bayes_GD(input_data_dict::Dict)
Function that takes in the a dictionary of extracted data (see `DataExtraction.extract_data`) and performs an custom iterative fitting algorithm to identify model parameters `monthly_p`, ``\\phi,b,`` and ``k``. See default constant hyperparameter values for algorithm in source code. MCMC proposal mechanism is a random walk with manually tuned multivariate normal proposal distribution. Sampler is run using a combination of Turing and AdvancedMH. Returns a dictionary of the final sampled MCMC chains and fits in each iteration.

*Note: Does trim off the burn-in period and needs to be done in post. But if `save_output=true` (by default), then a copy of the final chain with burn-in trimmed, and the metadata is saved to the path specified by chain_output_dir. (Refer to source code for path)*
"""
function bayes_GD(input_dict;
                N_EPOCHS = N_EPOCHS, EPOCH_LEN = EPOCH_LEN,
                α = α, ϵ = ϵ, NORMALISE_PERIOD = NORMALISE_PERIOD,
                proposal_resampling_variances = proposal_resampling_variances,
                iterations = iterations, burn_in = burn_in, n_chains = n_chains,
                save_output = true, chain_output_dir = chain_output_dir,
                filename = nothing)

    println("Starting EM regression for net crop for $(n_chains) chains with $(Threads.nthreads()) threads")

    # Extract data from dict
    YEARS_ANNUAL = input_dict["YEARS_ANNUAL"]
    MONTHS_MONTHLY = input_dict["MONTHS_MONTHLY"]
    DELIVERIES_ANNUAL = input_dict["DELIVERIES_ANNUAL"]
    DISTRIBUTION_ANNUAL = input_dict["DISTRIBUTION_ANNUAL"]
    cITN_CROP_MONTHLY = input_dict["cITN_CROP_MONTHLY"]
    NET_CROP_MONTHLY = input_dict["NET_CROP_MONTHLY"]
    # NET_CROP_MONTHLY_SMOOTHED = input_dict["NET_CROP_MONTHLY_SMOOTHED"]
    NET_CROP_STD_MONTHLY = input_dict["NET_CROP_STD_MONTHLY"]
    NET_NAMES = input_dict["NET_NAMES"]
    n_net_types = length(NET_NAMES)

    # Storage variable for output
    Γ_epochs_BYNET = zeros(N_EPOCHS+1, length(MONTHS_MONTHLY), n_net_types)
    Γ_epochs_TOTAL = zeros(N_EPOCHS+1, length(MONTHS_MONTHLY))

    # Adjust RWMH proposal sampling variances to account for additional net types
    resampling_variances = vcat(proposal_resampling_variances[1:end-2], repeat(proposal_resampling_variances[end-1:end], n_net_types))

    # Add additional varainces to RWMH proposal for random variable accounting for initial number of nets
    resampling_variances = vcat(resampling_variances, 0.8)

    # Add additional variances to RWMH proposal to account for missing data
    n_missing = sum(ismissing.(DISTRIBUTION_ANNUAL[:,1]))
    resampling_variances = vcat(resampling_variances, ones(n_missing).*0.15)

    # Define Random Walk Metropolis Hastings external sampler. Taken from AdvancedMH package.
    rwmh = externalsampler(RWMH(MvNormal(zeros(length(resampling_variances)), Diagonal(resampling_variances))))

    # Create storage variables to track progress of optimisation
    monthly_p_weights = [] # Weights for each month interval
    outputs = [] # Collection of sampled MCMC chains for each epoch
    loss_vals = [] # Loss values across entire optimisation for each descent step
    epoch_loss_vals = [] # Loss values after each epoch

    # Temporary initialisation variables (will be replaced in loop)
    monthly_p = []
    Γ_MONTHLY = []

    α_eff = α # learning rate

    # Identify time/month index in data to include in SGD calculations (excludes nans and infs)
    target_idx = findall(.!ismissing.(NET_CROP_MONTHLY))
    target_idx = target_idx[findall((.!isnan.(NET_CROP_STD_MONTHLY[target_idx])).*(.!isinf.(NET_CROP_STD_MONTHLY[target_idx])))]

    # Bool for early end
    early_regression_end = false

    # Iterate for each epoch
    for n in 1:(N_EPOCHS+1)
        # Initialisation - Assume that distribution is uniform across each month in the year
        if n == 1
            monthly_p = normalised_monthly_weights(rand(length(MONTHS_MONTHLY))) # initialise iterations
        end

        ### Part 1: Bayesian Inference

        # Sample the MCMC chain, given the selected monthly_p
        if n_chains > 1 # Sample in multithreaded parallel
            output = sample(monthly_reduced_itn(YEARS_ANNUAL, MONTHS_MONTHLY,
                            DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                            cITN_CROP_MONTHLY,
                            NET_CROP_MONTHLY, NET_CROP_STD_MONTHLY; 
                            monthly_p = monthly_p), rwmh, MCMCThreads(), iterations, n_chains, 
                            progress = true, discard_initial = burn_in)
        else # Sample single thread and chain
            output = sample(monthly_reduced_itn(YEARS_ANNUAL, MONTHS_MONTHLY,
                            DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                            cITN_CROP_MONTHLY,
                            NET_CROP_MONTHLY, NET_CROP_STD_MONTHLY; 
                            monthly_p = monthly_p), rwmh, iterations,
                            progress = true, discard_initial = burn_in)
        end
        

        # %%
        # Save a copy of the output chain into collection
        push!(outputs, deepcopy(output))

        
        # Extract estimates of parameter fits from MCMC chain, excluding a burn_in
        ϕ_vals = output[:ϕ].data[1:end]
        α_init_vals = output[:α_init].data[1:end]
        α_LLIN_vals = output[:α_LLIN].data[1:end]
        b_net_vals = Matrix(DataFrame(output)[:,7:2:7+2*(n_net_types-1)])
        k_net_vals = Matrix(DataFrame(output)[:,8:2:8+2*(n_net_types-1)])

        n_missing_nets_vals = sum(ismissing.(DISTRIBUTION_ANNUAL[:,1]))
        missing_nets_vals = Matrix(DataFrame(output)[:,((7+2*(n_net_types-1))+2):(7+2*(n_net_types-1))+n_missing_nets_vals+1])
        
        ϕ_est = mean(ϕ_vals)
        α_init_est = mean(α_init_vals)
        α_LLIN_est = mean(α_LLIN_vals)
        b_net_est = mean(b_net_vals, dims = 1)[:]
        k_net_est = mean(k_net_vals, dims = 1)[:]
        missing_nets_est = mean(missing_nets_vals, dims = 1)[:]

        # Extract specific parameters values related to LLINs (only used for printout)
        llin_entry_idx = findfirst(NET_NAMES .== "LLIN")
        # LLIN_halflife = net_life_weibull(0.5, b_net_est[llin_entry_idx], k_net_est[llin_entry_idx])
        LLIN_halflife = net_life_compact(0.5, b_net_est[llin_entry_idx], k_net_est[llin_entry_idx])

        # println("B.I. Phase Iter [$(n)/$(N_EPOCHS)]...ϕ = $(round(ϕ_est, digits = 4)), LLIN_b = $(round(b_net_est[llin_entry_idx], digits = 4)), LLIN_k = $(round(k_net_est[llin_entry_idx], digits = 4)), σ = $(round(σ_est, digits = 4)), LLIN halflife = $(round(LLIN_halflife, digits = 4)) years")
        println("B.I. Phase Iter [$(n)/$(N_EPOCHS)]...ϕ = $(round(ϕ_est, digits = 4)), LLIN_b = $(round(b_net_est[llin_entry_idx], digits = 4)), LLIN_k = $(round(k_net_est[llin_entry_idx], digits = 4)), LLIN halflife = $(round(LLIN_halflife, digits = 4)) years")

        ### Part 2: Stochastic Gradient Descent Tweaking
        # Calculate the net crop estimated time series, given mean regressed parameters values
        Γ_epochs_BYNET[n,:,:] = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                    ϕ_est, b_net_est, k_net_est, 
                                    α_init_est, α_LLIN_est,
                                    missing_nets_est; 
                                    monthly_p = monthly_p)
        
        Γ_epochs_TOTAL[n,:] = sum(Γ_epochs_BYNET[n,:,:], dims = 2)[:]

        # %%
        
        Γ_MONTHLY = Γ_epochs_TOTAL[n,:]

        # Select sample indexes over which to optimise DIC
        DIC_sample_idx = sample(1:length(ϕ_vals), DIC_SAMPLE_SIZE, replace = false)
        ϕ_sample_vals = ϕ_vals[DIC_sample_idx]
        b_sample_vals = b_net_vals[DIC_sample_idx,:]
        k_sample_vals = k_net_vals[DIC_sample_idx,:]
        α_init_sample_vals = α_init_vals[DIC_sample_idx]
        α_LLIN_sample_vals = α_LLIN_vals[DIC_sample_idx]
        missing_nets_sample_vals = zeros(size(missing_nets_vals))
        if size(missing_nets_sample_vals)[1]>0 #i.e. there was missing data that was regressed
            missing_nets_sample_vals = missing_nets_vals[DIC_sample_idx,:]
        else
            missing_nets_sample_vals = zeros(length(ϕ_sample_vals),0)
        end

        # Calculate Loss function as the DIC
        loss = DIC_netcrop(YEARS_ANNUAL, MONTHS_MONTHLY,
                        DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL, 
                        NET_CROP_MONTHLY, NET_CROP_STD_MONTHLY,
                        target_idx, monthly_p,
                        ϕ_sample_vals, 
                        b_sample_vals, 
                        k_sample_vals, 
                        α_init_sample_vals, 
                        α_LLIN_sample_vals, 
                        missing_nets_sample_vals; 
                        full_mode = true, printout = true)

        # Save the loss values and monthly weights
        if (n == 1) || (n == (N_EPOCHS+1)) || (early_regression_end == true)
            push!(loss_vals, loss)
            push!(epoch_loss_vals, loss)
            push!(monthly_p_weights, deepcopy(monthly_p))
            if n==1
                println("Initial Pre GD DIC Loss: $(epoch_loss_vals[end])")
            else
                println("Final DIC Loss: $(epoch_loss_vals[end])")

                if early_regression_end
                    Γ_epochs_BYNET = Γ_epochs_BYNET[n,:,:]
                    Γ_epochs_TOTAL = Γ_epochs_TOTAL[n,:]
                    break
                end
            end
        end

        # No points to regress against
        if length(target_idx) == 0
            break
        end

        if (n <= N_EPOCHS)&&(early_regression_end == false)
            # Gradient descent optimisation
            for n_iter in 1:EPOCH_LEN
                # Stochastic Gradient Descent - Resample indexes again
                DIC_sample_idx = sample(1:length(ϕ_vals), DIC_SAMPLE_SIZE, replace = false)
                ϕ_sample_vals = ϕ_vals[DIC_sample_idx]
                b_sample_vals = b_net_vals[DIC_sample_idx,:]
                k_sample_vals = k_net_vals[DIC_sample_idx,:]
                α_init_sample_vals = α_init_vals[DIC_sample_idx]
                α_LLIN_sample_vals = α_LLIN_vals[DIC_sample_idx]
                missing_nets_sample_vals = zeros(size(missing_nets_vals))
                if size(missing_nets_sample_vals)[1]>0 #i.e. there was missing data that was regressed
                    missing_nets_sample_vals = missing_nets_vals[DIC_sample_idx,:]
                else
                    missing_nets_sample_vals = zeros(length(ϕ_sample_vals),0)
                end

                # Using finite difference to calculate gradient of loss w.r.t perturbation in each month weight
                # dLdδ = zeros(length(target_idx), length(monthly_p)) # For MSE GD algorithm
                dLdδ = zeros(length(monthly_p)) # For DIC GD algorithm
                Threads.@threads for dim_i in 1:maximum(target_idx) #length(monthly_p)
                    # Central finite difference to numerically estimate derivative
                    δ = zeros(length(monthly_p))

                    ##### DIC Loss function
                    δ[dim_i] = min(monthly_p[dim_i]*(ϵ/2), ϵ*0.1) # Put a lower bound on interval size
                    # ϕ_vals_thinned = ϕ_vals[1:DIC_thin:end]
                    # # σ_vals_thinned = σ_vals[1:DIC_thin:end]
                    # α_init_vals_thinned = α_init_vals[1:DIC_thin:end]

                    # b_net_vals_thinned = b_net_vals[1:DIC_thin:end,:]
                    # k_net_vals_thinned = k_net_vals[1:DIC_thin:end,:]
                    # missing_nets_vals_thinned = missing_nets_vals[1:DIC_thin:end,:]

                    monthly_p_plus = normalised_monthly_weights(monthly_p.+δ)
                    monthly_p_minus = normalised_monthly_weights(max.(monthly_p.-δ,0)) # Need to prevent negative distribution weights
                    
                    DIC_plus = DIC_netcrop(YEARS_ANNUAL, MONTHS_MONTHLY,
                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL, 
                                    NET_CROP_MONTHLY, NET_CROP_STD_MONTHLY,
                                    target_idx, monthly_p_plus,
                                    ϕ_sample_vals, 
                                    b_sample_vals, 
                                    k_sample_vals, 
                                    α_init_sample_vals, 
                                    α_LLIN_sample_vals, 
                                    missing_nets_sample_vals; 
                                    full_mode = true)
                    DIC_minus = DIC_netcrop(YEARS_ANNUAL, MONTHS_MONTHLY,
                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL, 
                                    NET_CROP_MONTHLY, NET_CROP_STD_MONTHLY,
                                    target_idx, monthly_p_minus,
                                    ϕ_sample_vals, 
                                    b_sample_vals, 
                                    k_sample_vals,  
                                    α_init_sample_vals, 
                                    α_LLIN_sample_vals, 
                                    missing_nets_sample_vals; 
                                    full_mode = true)
                    dLdδ[dim_i] = (DIC_plus.-DIC_minus)./(2*δ[dim_i])
                end

                # Calculate resulting gradient averaged across all observations + normalise
                # ∇ = sum(dLdδ, dims = 1)[:] # MSE Version
                # ∇ = -∇./norm(∇) # MSE Version
                ∇ = dLdδ # DIC Version
                ∇ = -∇./norm(∇) # DIC Version

                # Update month weights
                # monthly_p = max.(monthly_p .+ ∇.*(α_eff*sqrt(length(monthly_p))),1) # MSE Version
                # monthly_p = max.(monthly_p .+ ∇.*(α_eff*sqrt(maximum(target_idx))),1) # MSE Version adjusted for latest survey month
                # monthly_p = max.(monthly_p .+ ∇.*(α_eff*sqrt(maximum(target_idx))),1e-4) # MSE Version adjusted for latest survey month
                monthly_p = monthly_p.*(1 .+ ∇.*(α_eff*sqrt(maximum(target_idx)))) # DIC Version
                
                # Generate time series with updated monthly_p
                Γ_MONTHLY_BYNET = model_evolve_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                                DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL,
                                                ϕ_est, b_net_est, k_net_est, 
                                                α_init_est, α_LLIN_est,
                                                missing_nets_est; 
                                                monthly_p = monthly_p)
                Γ_MONTHLY = sum(Γ_MONTHLY_BYNET, dims = 2)

                loss = DIC_netcrop(YEARS_ANNUAL, MONTHS_MONTHLY,
                                    DELIVERIES_ANNUAL, DISTRIBUTION_ANNUAL, 
                                    NET_CROP_MONTHLY, NET_CROP_STD_MONTHLY,
                                    target_idx, monthly_p,
                                    ϕ_sample_vals, 
                                    b_sample_vals, 
                                    k_sample_vals, 
                                    α_init_sample_vals, 
                                    α_LLIN_sample_vals, 
                                    missing_nets_sample_vals; 
                                    full_mode = true,
                                    printout = true)

                # Store results in variables
                push!(loss_vals, loss)
                
                if ((n_iter % 5) == 0)||(n_iter == 1)
                    println("Epoch [$n/$N_EPOCHS], Iteration [$n_iter/$EPOCH_LEN] DIC = $(loss-loss_vals[1])")
                    flush(stdout)
                end

                # Normalise monthly_p after NORMALISE_PERIOD EPOCH STEPS
                if (n_iter % NORMALISE_PERIOD) == 0
                    monthly_p = normalised_monthly_weights(monthly_p)
                end
            end

            # Check for early regression end (if no improvement is gained)
            if loss_vals[end] > epoch_loss_vals[end]
                # Predictions end regression here
                early_regression_end = true
                monthly_p = monthly_p_weights[end]

                println("G.D. Phase Current [$(n)/$(N_EPOCHS)] ΔDIC = $(loss_vals[end]-epoch_loss_vals[1]).")
                println("No improvement gained. Backtracking and ending regression.")
            else
                # There was improvement in the regressions. Carry on...

                #### Fill in remaining monthly_p
                # Fill in monthly_p for later years after most recent survey based on regressed historical monthly_p
                n_ref_years = (maximum(target_idx.-1)÷12)+1 # Number of years in reference
                n_ref_months = 12*n_ref_years # Number of reference years

                # Get monthly p_ref and normalise to sum to 1
                p_ref_matrix = Matrix(reshape(monthly_p[1:n_ref_months], 12,:)')
                normalised_p_ref_matrix = zeros(size(p_ref_matrix))
                for i in 1:size(p_ref_matrix)[1]
                    normalised_p_ref_matrix[i,:] = p_ref_matrix[i,:]./sum(p_ref_matrix[i,:])
                end

                # Average monthly p
                monthly_p_mean = mean(p_ref_matrix, dims = 1)[:]
                monthly_p_mean = monthly_p_mean./sum(monthly_p_mean) # Normalise to distribution within a year.

                # Enforce historical average of monthly_p to the remainder of the time series
                for i in (n_ref_years+1):length(YEARS_ANNUAL)
                    monthly_p[(12*(i-1)+1):(12*i)] .= copy(monthly_p_mean)
                end

                # Store epoch results variables
                push!(epoch_loss_vals, loss_vals[end])
                push!(monthly_p_weights, deepcopy(monthly_p))

                println("G.D. Phase Current [$(n)/$(N_EPOCHS)] ΔDIC = $(epoch_loss_vals[end]-epoch_loss_vals[1])")
            end

            # # Scheduler for loss adjustment
            # if epoch_loss_vals[end] > epoch_loss_vals[end-1]
            #     α_eff = α_eff*0.5
            #     monthly_p = monthly_p_weights[argmin(epoch_loss_vals)]
            #     println("Learning Rate Decreased: α_eff = $(α_eff). Backstepping to previous best...")
            # end
            
        end
    end

    # Package results to dict and return as function output
    regression_dict = Dict("NET_CROP_EPOCHS" => Γ_epochs_TOTAL,
                "MONTHLY_WEIGHTS" => monthly_p_weights,
                "CHAINS_EPOCHS" => outputs,
                "LOSS" => loss_vals,
                "LOSS_EPOCHS" => epoch_loss_vals,
                "burn_in" => burn_in,
                "resampling_variances" => resampling_variances)

    # Save last epoch's chain result in directory
    if save_output
        ISO = input_dict["ISO"]
        YEAR_START = input_dict["YEAR_START"]
        YEAR_END = input_dict["YEAR_END"]
        burn_in = regression_dict["burn_in"]
        resampling_variances = regression_dict["resampling_variances"]
        loss = regression_dict["LOSS"]
        loss_epochs = regression_dict["LOSS_EPOCHS"]
        monthly_weights = regression_dict["MONTHLY_WEIGHTS"]

        # Get Fit results
        monthly_p = regression_dict["MONTHLY_WEIGHTS"][end]
        chain = DataFrame(regression_dict["CHAINS_EPOCHS"][end])[:,3:end-1]

        # Also save historical fits
        chain_epochs = []
        for i in 1:length(outputs)
            push!(chain_epochs, DataFrame(regression_dict["CHAINS_EPOCHS"][i])[:,3:end-1])
        end

        # Define filename
        # Default filepath
        if isnothing(filename)
            chain_output_dir = chain_output_dir*"$(YEAR_START)_$(YEAR_END)/"
            filename = ISO*"_"*"$(YEAR_START)"*"_"*"$(YEAR_END)"*"_cropchains.jld2"
        end
        # Make filepath if it doesn't already exist
        mkpath(chain_output_dir)

        chain_output_filepath = chain_output_dir*filename

        # Save .jld2 file
        jldsave(chain_output_filepath; ISO, YEAR_START, YEAR_END,
                    chain, chain_epochs, 
                    monthly_p, burn_in, resampling_variances,
                    loss, loss_epochs, monthly_weights,
                    NET_NAMES)

        println("File saved at: "*chain_output_filepath)
    end

    return regression_dict

end



end
