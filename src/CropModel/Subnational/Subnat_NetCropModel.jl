"""
Author: Eugene Tan
Date Created: 4/11/2024
Last Updated: 11/11/2024
MITN model for subnational
"""

module Subnat_NetCropModel
export SUBNAT_monthly_itn
export SUBNAT_monthly_itn_forward
export subnat_weighted_monthly_itn
export nat_subnat_calibration

# Import Public Packages
using DataFrames
using Missings
using LinearAlgebra
using StatsBase
using Distributions
using NetLoss
using Turing

# %% Define required functions
"""
    SUBNAT_monthly_itn(YEARS_ANNUAL, MONTHS_MONTHLY,
                                    POPULATION_ANNUAL,
                                    DISTRIBUTION_ANNUAL_BYNET,
                                    SUBNAT_NET_CROP_MONTHLY,
                                    SUBNAT_NET_CROP_STD_MONTHLY,
                                    NET_NAMES,
                                    initial_demography_condition,
                                    net_attrition_params;
                                    monthly_p = nothing,
                                    MISSING_NETS_SCALE = MISSING_NETS_SCALE)

Function of statistical evolution model for net crop component of the stock and flow model. 
This function uses the `@model` macro to pass the function into the Sampler construction in Turing.jl. 
The variable `monthly_p` contains real values for the relative proportions of each annual distribution quota that is distributed in each given month.
"""
MISSING_NETS_SCALE = 1e6

@model function SUBNAT_monthly_itn(YEARS_ANNUAL, MONTHS_MONTHLY,
                                    POPULATION_ANNUAL,
                                    DISTRIBUTION_ANNUAL_BYNET,
                                    SUBNAT_NET_CROP_MONTHLY,
                                    SUBNAT_NET_CROP_STD_MONTHLY,
                                    NET_NAMES,
                                    initial_demography_condition,
                                    net_attrition_params;
                                    monthly_p = nothing,
                                    MISSING_NETS_SCALE = MISSING_NETS_SCALE)

    # Get basic constants and indexing bounds
    n_net_types = length(NET_NAMES)

    # For linking with net age demography bounds of national model
    YEAR_START = YEARS_ANNUAL[1]
    YEAR_END = YEARS_ANNUAL[end]
    max_age_months = size(initial_demography_condition[:,6:end])[2]
    YEAR_HISTORY_START = YEAR_START-Int(max_age_months/12)

    # Declare all priors to inference
    #### Net Attrition Parameters

    τ_nets = Vector{Real}(undef, n_net_types)
    κ_nets = Vector{Real}(undef, n_net_types)

    for net_type_i in 1:n_net_types
        # Get net name
        NET_NAME = NET_NAMES[net_type_i]

        # Check if we have any prior knowledge of net characteristics
        if NET_NAME ∈ net_attrition_params.NET_TYPE # there is prior knowledge
            τ_μ, τ_σ = net_attrition_params[net_attrition_params.NET_TYPE .== NET_NAME, ["kde_gaussian_tau_mu", "kde_gaussian_tau_sigma"]][1,:]
            κ_μ, κ_σ = net_attrition_params[net_attrition_params.NET_TYPE .== NET_NAME, ["kde_gaussian_kappa_mu", "kde_gaussian_kappa_sigma"]][1,:]

            # Update initial distributions
            τ_nets[net_type_i] ~ truncated(Normal(τ_μ, τ_σ),0.1,40)
            κ_nets[net_type_i] ~ truncated(Normal(κ_μ, κ_σ),1,30)
        else # Just take default uniform priors
            τ_nets[net_type_i] ~ Uniform(0.1,40)
            κ_nets[net_type_i] ~ Uniform(1,30)
        end
    end

    #### Missing Distribution Values
    missing_dist_idxs = findall(ismissing.(DISTRIBUTION_ANNUAL_BYNET[:,1]))
    missing_dist_vals = Vector{Real}(undef, length(missing_dist_idxs)*n_net_types)
    for i in 1:length(missing_dist_idxs)
        for net_type_i in 1:n_net_types
            missing_dist_vals[(i-1)*n_net_types + net_type_i] ~ LogNormal(0,1)
        end
    end
    

    # Subnational SNF Forward evolution

    #### Normalise monthly proportions on annual scale
    if isnothing(monthly_p)
        monthly_p = ones(length(MONTHS_MONTHLY)) # Base case assumes regular monthly equal distributions
    end

    effective_monthly_proportions = zeros(length(MONTHS_MONTHLY))
    DISTRIBUTION_MONTHLY_BYNET = zeros(length(MONTHS_MONTHLY), n_net_types)

    # Disaggregate each year into monthly resolution
    for i in 1:length(YEARS_ANNUAL)
        effective_monthly_proportions[(12*(i-1)+1):(12*i)] = monthly_p[(12*(i-1)+1):(12*i)]./sum(monthly_p[(12*(i-1)+1):(12*i)])

        if i ∈ missing_dist_idxs # Missing annual net distribution data, need to use imputed random variable
            idx = findfirst(missing_dist_idxs .== i)
            for net_type_i in 1:n_net_types
                DISTRIBUTION_MONTHLY_BYNET[(12*(i-1)+1):(12*i),net_type_i] = effective_monthly_proportions[(12*(i-1)+1):(12*i)].*(missing_dist_vals[(idx-1)*n_net_types + net_type_i]*MISSING_NETS_SCALE)
            end
        else # Data is available, just disaggregate as per normal
            for net_type_i in 1:n_net_types
                DISTRIBUTION_MONTHLY_BYNET[(12*(i-1)+1):(12*i),net_type_i] = effective_monthly_proportions[(12*(i-1)+1):(12*i)].*DISTRIBUTION_ANNUAL_BYNET[i,net_type_i]
            end
        end
    end

    # Forward simulation
    A_BYNET = zeros((YEAR_END - YEAR_START + 1)*12,(YEAR_END - YEAR_HISTORY_START + 1)*12, n_net_types)

    # Initialise net demography for first time step

    for net_type_i in 1:n_net_types
        NET_NAME = NET_NAMES[net_type_i]
        
        # Use posterior mean final net demography condition at YEAR_START from national model to initialise
        average_initial_demography = mean(Matrix(initial_demography_condition[initial_demography_condition.NET_TYPE .== NET_NAME,6:end]), dims = 1)[:]
        A_BYNET[1,1:max_age_months,net_type_i] .= reverse(average_initial_demography.*POPULATION_ANNUAL[1])
    end

    # Start forward simulation for net type i
    for net_type_i in 1:n_net_types
        # Define net attrition parameters from prior estimates
        τ_est = τ_nets[net_type_i]
        κ_est = κ_nets[net_type_i]

        # Propagate distribution forward to calculate estimate age demography based on selected attrition parameters
        for i in 1:length(MONTHS_MONTHLY) # Current time
            if i >= 2 #Simulate time passing
                for j in 1:i+max_age_months-1 # Birth Times
                    birth_month = j
                    current_month = i+max_age_months
                    age_months = current_month-birth_month
                    net_survival_probability = net_loss_compact(age_months/12, τ_est, κ_est)/net_loss_compact((age_months-1)/12, τ_est, κ_est) #p(survival to next step, given current age)
                    if isnan(net_survival_probability) # probabilities are too small, just set survival to 0
                        net_survival_probability = 0
                    end
                    A_BYNET[i,j,net_type_i] = A_BYNET[i-1,j,net_type_i]*net_survival_probability
                end
            end
            A_BYNET[i,i+max_age_months,net_type_i] += DISTRIBUTION_MONTHLY_BYNET[i,net_type_i]
        end
    end

    # Tally net totals across age groups
    Γ_MONTHLY_TOTAL_BYNET = sum(A_BYNET, dims = 2)[:,1,:]

    # Tally net totals across age groups and net type
    Γ_MONTHLY_TOTAL = sum(Γ_MONTHLY_TOTAL_BYNET, dims = 2)[:,1]

    # Regress all on missing values
    regression_idx = findall(.!ismissing.(SUBNAT_NET_CROP_MONTHLY))
    for idx in regression_idx
        # Ignore problematic/error data
        if (SUBNAT_NET_CROP_STD_MONTHLY[idx] == 0)||(isnan(SUBNAT_NET_CROP_STD_MONTHLY[idx]))||(isinf(SUBNAT_NET_CROP_STD_MONTHLY[idx]))
            continue
        else
            SUBNAT_NET_CROP_MONTHLY[idx] ~ Normal(Γ_MONTHLY_TOTAL[idx], SUBNAT_NET_CROP_STD_MONTHLY[idx])
        end
    end
end
##########################################
# Posterior model prediction
##########################################
# %%

"""
    SUBNAT_monthly_itn_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                            POPULATION_ANNUAL,
                            DISTRIBUTION_ANNUAL_BYNET,
                            NET_NAMES,
                            initial_demography_condition,
                            τ_est, κ_est,
                            missing_nets_est;
                            monthly_p = nothing,
                            MISSING_NETS_SCALE = MISSING_NETS_SCALE)
                                
Forward simulation of the net crop model as above. Requires parameter values for ``\\phi``, ``b`` and ``k`` to run. By default, these are set to 0.5, 2 and 4 respectively. If no 'monthly_p' vector is given, a uniform monthly distribution is assumed. Returns a time series as a vector of real numbers of the estimated national net crop.
NEED TO WRITE PROPER DOCUMENTATION!!
"""
function SUBNAT_monthly_itn_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                            POPULATION_ANNUAL,
                            DISTRIBUTION_ANNUAL_BYNET,
                            NET_NAMES,
                            initial_demography_condition,
                            τ_est, κ_est,
                            missing_nets_est;
                            monthly_p = nothing,
                            MISSING_NETS_SCALE = MISSING_NETS_SCALE)

    # Get basic constants and indexing bounds
    n_net_types = length(NET_NAMES)

    # For linking with net age demography bounds of national model
    YEAR_START = YEARS_ANNUAL[1]
    YEAR_END = YEARS_ANNUAL[end]
    max_age_months = size(initial_demography_condition[:,6:end])[2]
    YEAR_HISTORY_START = YEAR_START-Int(max_age_months/12)

    # Declare all priors to inference

    #### Missing Distribution Values
    missing_dist_idxs = findall(ismissing.(DISTRIBUTION_ANNUAL_BYNET[:,1]))
    
    # Subnational SNF Forward evolution

    #### Normalise monthly proportions on annual scale
    if isnothing(monthly_p)
        monthly_p = ones(length(MONTHS_MONTHLY)) # Base case assumes regular monthly equal distributions
    end

    effective_monthly_proportions = zeros(length(MONTHS_MONTHLY))
    DISTRIBUTION_MONTHLY_BYNET = zeros(length(MONTHS_MONTHLY), n_net_types)

    # Disaggregate each year into monthly resolution
    for i in 1:length(YEARS_ANNUAL)
        effective_monthly_proportions[(12*(i-1)+1):(12*i)] = monthly_p[(12*(i-1)+1):(12*i)]./sum(monthly_p[(12*(i-1)+1):(12*i)])

        if i ∈ missing_dist_idxs # Missing annual net distribution data, need to use imputed random variable
            idx = findfirst(missing_dist_idxs .== i)
            for net_type_i in 1:n_net_types
                DISTRIBUTION_MONTHLY_BYNET[(12*(i-1)+1):(12*i),net_type_i] = effective_monthly_proportions[(12*(i-1)+1):(12*i)].*(missing_nets_est[(idx-1)*n_net_types + net_type_i]*MISSING_NETS_SCALE)
            end
        else # Data is available, just disaggregate as per normal
            for net_type_i in 1:n_net_types
                DISTRIBUTION_MONTHLY_BYNET[(12*(i-1)+1):(12*i),net_type_i] = effective_monthly_proportions[(12*(i-1)+1):(12*i)].*DISTRIBUTION_ANNUAL_BYNET[i,net_type_i]
            end
        end
    end

    # Forward simulation
    A_BYNET = zeros((YEAR_END - YEAR_START + 1)*12,(YEAR_END - YEAR_HISTORY_START + 1)*12, n_net_types)

    # Initialise net demography for first time step

    for net_type_i in 1:n_net_types
        NET_NAME = NET_NAMES[net_type_i]
        
        # Use posterior mean final net demography condition at YEAR_START from national model to initialise
        average_initial_demography = mean(Matrix(initial_demography_condition[initial_demography_condition.NET_TYPE .== NET_NAME,6:end]), dims = 1)[:]
        if isempty(Matrix(initial_demography_condition[initial_demography_condition.NET_TYPE .== NET_NAME,6:end])) # i.e. net not present at the time
            average_initial_demography .= 0
        end
        A_BYNET[1,1:max_age_months,net_type_i] .= reverse(average_initial_demography.*POPULATION_ANNUAL[1])
    end

    # Start forward simulation for net type i
    for net_type_i in 1:n_net_types

        # Propagate distribution forward to calculate estimate age demography based on selected attrition parameters
        for i in 1:length(MONTHS_MONTHLY) # Current time
            if i >= 2 #Simulate time passing
                for j in 1:i+max_age_months-1 # Birth Times
                    birth_month = j
                    current_month = i+max_age_months
                    age_months = current_month-birth_month
                    #p(survival to next step, given current age)
                    net_survival_probability = net_loss_compact(age_months/12, τ_est[net_type_i], κ_est[net_type_i])/net_loss_compact((age_months-1)/12, τ_est[net_type_i], κ_est[net_type_i])
                    if isnan(net_survival_probability) # probabilities are too small, just set survival to 0
                        net_survival_probability = 0
                    end
                    A_BYNET[i,j,net_type_i] = A_BYNET[i-1,j,net_type_i]*net_survival_probability
                end
            end
            A_BYNET[i,i+max_age_months,net_type_i] += DISTRIBUTION_MONTHLY_BYNET[i,net_type_i]
        end
    end

    # Tally net totals across age groups    
    Γ_MONTHLY_TOTAL_BYNET = sum(A_BYNET, dims = 2)[:,1,:]

    return Γ_MONTHLY_TOTAL_BYNET, A_BYNET
end

"""
    subnat_weighted_monthly_itn(YEARS_ANNUAL, MONTHS_MONTHLY,
                                        ALLREGIONS_POPULATION_ANNUAL_EXTRACT,
                                        ALLREGIONS_DISTRIBUTIONS_ARRAY,
                                        NAT_DISTRIBUTIONS_TOTAL,
                                        ALLREGIONS_NET_NAMES,
                                        net_attrition_params,
                                        initial_demography_condition,
                                        ALLREGIONS_NET_CROP_MONTHLY,
                                        ALLREGIONS_NET_CROP_STD_MONTHLY; 
                                        monthly_p = nothing,
                                        ϵ = 0.003)
                                
Function used for finding adjustment weights of the net distributions such that national netcrop and household net crop 
    match well.
    !!! NEED TO WRITE PROPER DOCUMENTATION
"""

@model function subnat_weighted_monthly_itn(YEARS_ANNUAL, MONTHS_MONTHLY,
                                        ALLREGIONS_POPULATION_ANNUAL_EXTRACT,
                                        ALLREGIONS_DISTRIBUTIONS_ARRAY,
                                        NAT_DISTRIBUTIONS_TOTAL,
                                        ALLREGIONS_NET_NAMES,
                                        net_attrition_params,
                                        initial_demography_condition,
                                        ALLREGIONS_NET_CROP_MONTHLY,
                                        ALLREGIONS_NET_CROP_STD_MONTHLY; 
                                        monthly_p = nothing,
                                        ϵ = 0.003)
        # Net crop conservation tolerance
        # ϵ = 0.003

        n_admin1 = length(ALLREGIONS_POPULATION_ANNUAL_EXTRACT)

        # Define weighting constants across n_admin1 regions
        ω = Vector{Float64}(undef, n_admin1)

        for i in 1:n_admin1
            ω[i] ~ LogNormal(0,1)
        end

        # Carry over posterior net_attrition parameters from national model
        τ_est = net_attrition_params.mean_tau
        κ_est = net_attrition_params.mean_kappa
        missing_nets_est = [] # Should have no missing nets because data will be taken from the national model

        # Storage variable
        ALLREGIONS_SUBNAT_Γ_MONTHLY_TOTAL = Matrix{Float64}(undef, n_admin1, length(MONTHS_MONTHLY))

        # Foward simulation of subnational model
        for admin1_name_i in 1:n_admin1
            # Select required input data
            POPULATION_ANNUAL = ALLREGIONS_POPULATION_ANNUAL_EXTRACT[admin1_name_i]
            NET_NAMES = ALLREGIONS_NET_NAMES[admin1_name_i]

            # Adjust annual net distribution for subnational region by weight
            ADJ_DISTRIBUTION_ANNUAL_BYNET = ALLREGIONS_DISTRIBUTIONS_ARRAY[admin1_name_i,:,:].*ω[admin1_name_i]

            # Subnational model forward calculation
            SUBNAT_Γ_MONTHLY_TOTAL_BYNET = SUBNAT_monthly_itn_forward(YEARS_ANNUAL, MONTHS_MONTHLY,
                                                                    POPULATION_ANNUAL,
                                                                    ADJ_DISTRIBUTION_ANNUAL_BYNET,
                                                                    NET_NAMES,
                                                                    initial_demography_condition,
                                                                    τ_est, κ_est,
                                                                    missing_nets_est;
                                                                    monthly_p = monthly_p)[1]
            # Save estimate of simulation to variable
            ALLREGIONS_SUBNAT_Γ_MONTHLY_TOTAL[admin1_name_i,:] = sum(SUBNAT_Γ_MONTHLY_TOTAL_BYNET, dims = 2)[:,1]
        end
        

        # Regression section
        # Make sure individual subnational data matches household surveys
        for admin1_name_i in 1:n_admin1
            # Find index of a all nonmissing data points to regress against
            nonmissing_idxs = findall(.!ismissing.(ALLREGIONS_NET_CROP_MONTHLY[admin1_name_i,:]))
            
            # Extract required data for regression
            SUBNAT_Γ_MONTHLY_TOTAL = ALLREGIONS_SUBNAT_Γ_MONTHLY_TOTAL[admin1_name_i,nonmissing_idxs][:]
            
            # Gaussian Regression
            for j in 1:length(nonmissing_idxs)
                ALLREGIONS_NET_CROP_MONTHLY[admin1_name_i,nonmissing_idxs[j]] ~ truncated(Normal(SUBNAT_Γ_MONTHLY_TOTAL[j], ALLREGIONS_NET_CROP_STD_MONTHLY[admin1_name_i,nonmissing_idxs[j]]), 0, Inf)
            end
        end

        # Make sure distribution totals do not differ too much from national estimates of distribution
        ALLREGIONS_DISTRIBUTIONS_TOTAL = sum(ALLREGIONS_DISTRIBUTIONS_ARRAY, dims = 3)[:,:,1]
        ADJ_NAT_DISTRIBUTIONS_TOTAL = ((ω')*(ALLREGIONS_DISTRIBUTIONS_TOTAL))[:]

        for i in 1:length(NAT_DISTRIBUTIONS_TOTAL)
            NAT_DISTRIBUTIONS_TOTAL[i] ~ Normal(ADJ_NAT_DISTRIBUTIONS_TOTAL[i], NAT_DISTRIBUTIONS_TOTAL[i].*ϵ)
        end

        # NAT_DISTRIBUTIONS_TOTAL ~ MvNormal(ADJ_NAT_DISTRIBUTIONS_TOTAL, Diagonal((NAT_DISTRIBUTIONS_TOTAL.*ϵ).^2))
    end

"""
    nat_subnat_calibration(NAT_NET_CROP_MONTHLY,
                            NAT_NET_CROP_STD_MONTHLY,
                            SUBNAT_NET_CROP_MONTHLY,
                            SUBNAT_NET_CROP_STD_MONTHLY,
                            ALLREGIONS_Γ_MONTHLY_TOTAL_mean)
                                
Function to regress and calculate weights for each admin region to reconcile against national net crop estimates.
    NEED TO WRITE PROPER DOCUMENTATION!!!!
"""

@model function nat_subnat_calibration(NAT_NET_CROP_MONTHLY,
    NAT_NET_CROP_STD_MONTHLY,
    SUBNAT_NET_CROP_MONTHLY,
    SUBNAT_NET_CROP_STD_MONTHLY,
    ALLREGIONS_Γ_MONTHLY_TOTAL_mean)

    # Get loop constants
    n_admin1 = size(ALLREGIONS_Γ_MONTHLY_TOTAL_mean)[1]

    # Priors for relative weights
    ω = Vector{Real}(undef, n_admin1)
    # α = Vector{Real}(undef, n_admin1)

    for i in 1:n_admin1
        ω[i] ~ truncated(Normal(1,1/3), 0, Inf)
        # ω[i] = rand(truncated(Normal(1,1/3), 0, Inf))
    end

    # Make such that total sum of weights are conserved
    # α = (ω./sum(ω))*n_admin1

    # # Calculated adjusted values of net crop
    # ADJ_ALLREGIONS_Γ_MONTHLY_TOTAL_mean = copy(ALLREGIONS_Γ_MONTHLY_TOTAL_mean)
    # ADJ_ALLREGIONS_Γ_MONTHLY_TOTAL_mean = ALLREGIONS_Γ_MONTHLY_TOTAL_mean.*repeat(α, 1,size(NAT_NET_CROP_MONTHLY)[1])
    # ADJ_NAT_Γ_MONTHLY_TOTAL_mean = sum(ADJ_ALLREGIONS_Γ_MONTHLY_TOTAL_mean, dims = 1)[:]

    # Do fitting for each admin region according to subnational surveys
    for i in 1:n_admin1
        admin1_reg_idxs = findall(.!ismissing.(SUBNAT_NET_CROP_MONTHLY[i,:]))
        for monthidx in admin1_reg_idxs
            # SUBNAT_NET_CROP_MONTHLY[i,monthidx] ~ truncated(Normal(ADJ_ALLREGIONS_Γ_MONTHLY_TOTAL_mean[i,monthidx], SUBNAT_NET_CROP_STD_MONTHLY[i,monthidx]), 0, Inf)
            # SUBNAT_NET_CROP_MONTHLY[i,monthidx] ~ truncated(Normal(ALLREGIONS_Γ_MONTHLY_TOTAL_mean[i,monthidx]*α[i], SUBNAT_NET_CROP_STD_MONTHLY[i,monthidx]), 0, Inf)
            SUBNAT_NET_CROP_MONTHLY[i,monthidx] ~ truncated(Normal(ALLREGIONS_Γ_MONTHLY_TOTAL_mean[i,monthidx]*ω[i], SUBNAT_NET_CROP_STD_MONTHLY[i,monthidx]), 0, Inf)
        end
    end

    # Do fitting for national household surveys
    reg_idxs = findall(.!ismissing.(NAT_NET_CROP_MONTHLY))
    for monthidx in reg_idxs

        # NAT_NET_CROP_MONTHLY[monthidx] ~ truncated(Normal(ADJ_NAT_Γ_MONTHLY_TOTAL_mean[monthidx], NAT_NET_CROP_STD_MONTHLY[monthidx]), 0, Inf)
        # NAT_NET_CROP_MONTHLY[monthidx] ~ truncated(Normal((ALLREGIONS_Γ_MONTHLY_TOTAL_mean[:,monthidx]')*α, NAT_NET_CROP_STD_MONTHLY[monthidx]), 0, Inf)
        NAT_NET_CROP_MONTHLY[monthidx] ~ truncated(Normal((ALLREGIONS_Γ_MONTHLY_TOTAL_mean[:,monthidx]')*ω, NAT_NET_CROP_STD_MONTHLY[monthidx]), 0, Inf)
    end
end

end


