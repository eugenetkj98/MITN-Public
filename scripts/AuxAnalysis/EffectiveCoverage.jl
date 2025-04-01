# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import Public Packages
using JLD2
using CSV
using DataFrames
using ProgressBars
using Measures
using Plots
using StatsBase
using Distributions
using NetLoss
pythonplot()
theme(:vibrant)

# %% Import helper packages
using PlottingFunctions
using NetCropModel
using NetAccessModel
using LaTeXStrings
using DynamicalSystems

# %% First do Bayesian Inference (old school way) for getting decay behaviour of ITNs
# Import Compiled data
data = CSV.read("datasets/Net_BioEfficacy_Dataset.csv", DataFrame)

# Filter out all new-gen LLINs
filt_data = data[findall(data[:, "Net Type"].!="LLIN-2"),:]

t_data = filt_data[:, "Net Age"]
mortality_data = filt_data[:,"Mean Mortality"]
sd_data = filt_data[:,"SD Mortality"]
# Priors
α = 9/2
β = 3/2

b_vals = 0:0.05:10
k_vals = 0:0.05:20

joint_post = zeros(length(b_vals), length(k_vals))
for i in ProgressBar(1:length(b_vals), leave = false)
    for j in 1:length(k_vals)
        b = b_vals[i]
        k = k_vals[j]

        d_prior = net_loss_weibull.(t_data, b, k)

        # Get Log components for bayesian inference
        log_p_b_prior = logpdf(Gamma(α, β), b)
        log_p_k_prior = log(pdf(Gamma(α, β), k))
        log_likelihood = zeros(length(t_data))
        for n in 1:length(t_data)
            log_likelihood[n] = log(pdf(Normal(d_prior[n], sd_data[n]),mortality_data[n]))
        end

        log_likelihood

        # Correct Infs because they will be low valued
        log_likelihood[findall(isinf.(log_likelihood))] .= -100

        log_posterior = mean(log_likelihood) + log_p_b_prior + log_p_k_prior

        

        joint_post[i,j] = exp(log_posterior)
    end
end

# Fix NaNs and Normalise to probability dist
joint_post[findall(isnan.(joint_post))] .= 0
joint_post = joint_post./sum(joint_post)

# %% Plot
heatmap(k_vals, b_vals , joint_post)

# %%
b_idx, k_idx = argmax(joint_post)[1], argmax(joint_post)[2]
b_est = b_vals[b_idx]
k_est = k_vals[k_idx]

# %%
t_vals = 0:0.01:5
fig = plot(xlabel = "Years", ylabel = "24h Mortality")
plot!(fig, t_vals,net_loss_weibull.(t_vals, b_est, k_est), label = "MLE")
scatter!(fig, t_data, mortality_data, label = "Raw Data")

savefig(fig, "bioefficacy_fit.pdf")

# %%
output_dir = "output_plots/"

# %% Get ISO List
ISO_list = String.(CSV.read(raw"C:\Users\ETan\Documents\Prototype Analyses\itn-updated\datasets\ISO_list.csv", DataFrame)[:,1])
reg_results_list = Array{Any, 1}(undef, length(ISO_list))
exclusion_ISOs = ["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD"]
#["CPV","BWA","CAF","COM","GNQ","DJI","ERI","ETH","GAB","GNB","STP","SOM","SDN","SWZ","ZAF","SSD"]

YEAR_START = 2000
YEAR_END = 2023

country_codes_key = CSV.read("datasets/country_codes.csv", DataFrame)
# %%
ISO = "KEN"
n_samples = 500

# Import Data
input_dict = load("outputs/extractions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
regression_dict = load("outputs/regressions/crop/Compact Regressions/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")
net_access_input_dict = load("outputs/extractions/access/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")
net_access_chain = load("outputs/regressions/access/netaccesschains.jld2")
BV_dict = load("datasets/BV_draws.jld2")
country_name = country_codes_key[findfirst(country_codes_key.ISO3 .== ISO), "Country"]

# Crop Regression input data
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
    

##### Generate Net Crop Trajectories
Γ_MONTHLY_samples_BYNET = zeros(n_samples, length(MONTHS_MONTHLY), n_net_types)
A_samples_BYNET = zeros(n_samples, length(MONTHS_MONTHLY), length(MONTHS_MONTHLY), n_net_types)

for i in 1:n_samples
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
end

# Scale by to units of millions
Γ_MONTHLY_samples_BYNET = Γ_MONTHLY_samples_BYNET./1e6
A_samples_BYNET = A_samples_BYNET./1e6

# Get Totals
Γ_MONTHLY_samples_TOTAL = sum(Γ_MONTHLY_samples_BYNET, dims = 3)[:,:,1]
A_samples_TOTAL = sum(A_samples_BYNET, dims = 4)[:,:,:,1]

# %%
heatmap(A_samples_TOTAL[1,:,:], yflip = true)

# %% Calculate Access
ρ_chain_df = net_access_chain["ρ_chain_df"]
μ_chain_df = net_access_chain["μ_chain_df"]
p_h = mean(net_access_input_dict["p_h_aggregated"], dims = 1)[:]
λ_access_samples = sample_net_access(ρ_chain_df, μ_chain_df, p_h,
                                        POPULATION_MONTHLY, Γ_MONTHLY_samples_TOTAL.*1e6) # Need to temporarily unscale input by mil for calculating access


# %%
##### Calculate Age Matrix
n_months = size(A_samples_TOTAL)[2]
M_age = zeros(n_months, n_months)

for i in 1:n_months # Current Time
    for j in 1:n_months # Birth time
        if j>i
            M_age[i,j] = 0
        else
            M_age[i,j] = (i-j)/12
        end
    end
end

# %%
# heatmap(M_age, yflip = true)

# %% Define Net Wear model
η_netwear = net_loss_weibull.(M_age, b_est,k_est)


# # %%
# heatmap(η_netwear, yflip = true)

# %%
t_vals = (1:n_months)./12

# Insecticide resistance time series
β_t = 1 .-net_loss_weibull.(t_vals, 18.0, 6.0)
ξ = 0.4 # Maximum Kill rate reduction
ρ_t = 1 .- ξ.*β_t # Resistance effects
fig = plot(t_vals.+YEAR_START, ρ_t, ylims = (-0.05, 1.05), legend = false)
plot!(fig, title = "b = $(18), k = $(6)", xlabel = "Year", ylabel = "ρ(t)")


# %% 

impact_unadjusted = zeros(size(λ_access_samples))
impact_IR_adjusted = zeros(size(λ_access_samples))
impact_age_adjusted = zeros(size(λ_access_samples))

impact_full_adjusted = zeros(size(λ_access_samples))

for i in 1:n_samples
    impact_unadjusted[i,:] = sum(A_samples_TOTAL[i,:,:].*λ_access_samples[i,:,:].*1e6, dims = 2)[:]./POPULATION_MONTHLY
    impact_IR_adjusted[i,:] = sum(A_samples_TOTAL[i,:,:].*λ_access_samples[i,:,:].*1e6, dims = 2)[:].*ρ_t./POPULATION_MONTHLY 
    impact_age_adjusted[i,:] = sum(A_samples_TOTAL[i,:,:].*λ_access_samples[i,:,:].*1e6.*η_netwear, dims = 2)[:]./POPULATION_MONTHLY 
    impact_full_adjusted[i,:] = sum(A_samples_TOTAL[i,:,:].*λ_access_samples[i,:,:].*1e6.*η_netwear, dims = 2)[:].*ρ_t./POPULATION_MONTHLY 
end

# %%
fig = plot(title = "$ISO", xlabel = "Year", ylabel = "Effective Efficacy per capita",
            xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), 
            xtickfontrotation = 90, legend = :left)
plot!(mean(impact_unadjusted, dims = 1)[:], label = "Unadjusted", color = 2)
# plot!(mean(impact_IR_adjusted, dims = 1)[:], label = "IR Adjusted", color = 2)
plot!(mean(impact_age_adjusted, dims = 1)[:], label = "Age Adjusted", color = 5)
plot!(mean(impact_full_adjusted, dims = 1)[:], label = "IR + Age Adjusted", color = 6)
# plot!(mean(impact_unadjusted, dims = 1)[:] .- mean(impact_IR_adjusted, dims = 1)[:], label = "IR Deficit", color = 2,
#         linestyle = :dash)
plot!(mean(impact_unadjusted, dims = 1)[:] .- mean(impact_age_adjusted, dims = 1)[:], label = "Age Deficit", color = 5,
        linestyle = :dash)
plot!(mean(impact_unadjusted, dims = 1)[:] .- mean(impact_full_adjusted, dims = 1)[:], label = "IR + Age Deficit", color = 6,
    linestyle = :dash)

plot!(twinx(), ρ_t, linewidth = 1.5, ylabel = "IR Efficacy Penalty", 
        linecolor = 4, linestyle = :dash, linealpha = 0.8, legend = false,
        ylims = (-0.02, 1.05))


savefig(fig, "$(ISO)_EffectiveCoverage.pdf")
 
#####################################
# %% Batch Analysis to make plots
#####################################

YEAR_START = 2000
YEAR_END = 2023
n_samples = 500
country_codes_key = CSV.read("datasets/country_codes.csv", DataFrame)
IR_b = 18.0
IR_k = 6.0
ξ = 0.4 # Maximum kill rate reduction


# %% Get ISO List
ISO_list = String.(CSV.read(raw"C:\Users\ETan\Documents\Prototype Analyses\itn-updated\datasets\ISO_list.csv", DataFrame)[:,1])
reg_results_list = Array{Any, 1}(undef, length(ISO_list))
exclusion_ISOs = ["CPV","BWA","CAF","GNQ","DJI","GAB","GNB","ERI","ETH","SOM","SDN","ZAF","SSD"]
filt_ISO_list = setdiff(ISO_list, exclusion_ISOs)

# %% Calculate penalty effects for each country
age_penalties = zeros(length(filt_ISO_list),(YEAR_END-YEAR_START+1)*12)
full_penalties = zeros(length(filt_ISO_list),(YEAR_END-YEAR_START+1)*12)

for i in ProgressBar(1:length(filt_ISO_list))
    # Select ISO
    ISO = filt_ISO_list[i]

    # Import Data
    input_dict = load("outputs/extractions/crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")
    regression_dict = load("outputs/regressions/crop/Compact Regressions/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropchains.jld2")
    net_access_input_dict = load("outputs/extractions/access/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_accessextract.jld2")
    net_access_chain = load("outputs/regressions/access/netaccesschains.jld2")
    BV_dict = load("datasets/BV_draws.jld2")
    country_name = country_codes_key[findfirst(country_codes_key.ISO3 .== ISO), "Country"]

    # Crop Regression input data
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
        

    ##### Generate Net Crop Trajectories
    Γ_MONTHLY_samples_BYNET = zeros(n_samples, length(MONTHS_MONTHLY), n_net_types)
    A_samples_BYNET = zeros(n_samples, length(MONTHS_MONTHLY), length(MONTHS_MONTHLY), n_net_types)

    for i in 1:n_samples
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
    end

    # Scale by to units of millions
    Γ_MONTHLY_samples_BYNET = Γ_MONTHLY_samples_BYNET./1e6
    A_samples_BYNET = A_samples_BYNET./1e6

    # Get Totals
    Γ_MONTHLY_samples_TOTAL = sum(Γ_MONTHLY_samples_BYNET, dims = 3)[:,:,1]
    A_samples_TOTAL = sum(A_samples_BYNET, dims = 4)[:,:,:,1]

    # %% Calculate Access
    ρ_chain_df = net_access_chain["ρ_chain_df"]
    μ_chain_df = net_access_chain["μ_chain_df"]
    p_h = mean(net_access_input_dict["p_h_aggregated"], dims = 1)[:]
    λ_access_samples = sample_net_access(ρ_chain_df, μ_chain_df, p_h,
                                            POPULATION_MONTHLY, Γ_MONTHLY_samples_TOTAL.*1e6) # Need to temporarily unscale input by mil for calculating access


    # %%
    ##### Calculate Age Matrix
    n_months = size(A_samples_TOTAL)[2]
    M_age = zeros(n_months, n_months)

    for i in 1:n_months # Current Time
        for j in 1:n_months # Birth time
            if j>i
                M_age[i,j] = 0
            else
                M_age[i,j] = (i-j)/12
            end
        end
    end

    # %% Define Net Wear model
    η_netwear = net_loss_weibull.(M_age, b_est,k_est)

    # %% Calculate IR effects
    t_vals = (1:n_months)./12

    # Insecticide resistance time series
    β_t = 1 .-net_loss_weibull.(t_vals, IR_b, IR_k)
    ρ_t = 1 .- ξ.*β_t # Resistance effects

    # %% Calculate impacts

    impact_unadjusted = zeros(size(λ_access_samples))
    impact_age_adjusted = zeros(size(λ_access_samples))
    impact_full_adjusted = zeros(size(λ_access_samples))

    for i in 1:n_samples
        impact_unadjusted[i,:] = sum(A_samples_TOTAL[i,:,:].*λ_access_samples[i,:,:].*1e6, dims = 2)[:]./POPULATION_MONTHLY
        impact_age_adjusted[i,:] = sum(A_samples_TOTAL[i,:,:].*λ_access_samples[i,:,:].*1e6.*η_netwear, dims = 2)[:]./POPULATION_MONTHLY 
        impact_full_adjusted[i,:] = sum(A_samples_TOTAL[i,:,:].*λ_access_samples[i,:,:].*1e6.*η_netwear, dims = 2)[:].*ρ_t./POPULATION_MONTHLY 
    end

    penalty_age = mean((impact_unadjusted.-impact_age_adjusted)./impact_unadjusted, dims = 1)[:]
    penalty_full = mean((impact_unadjusted.-impact_full_adjusted)./impact_unadjusted, dims = 1)[:]

    # Save penalty values
    age_penalties[i,:] = penalty_age
    full_penalties[i,:] = penalty_full
end

# %% Fix NaNs
age_penalties[findall(isnan.(age_penalties))] .= 0
full_penalties[findall(isnan.(full_penalties))] .= 0
# %%
MONTHS_MONTHLY = Vector(1:((YEAR_END-YEAR_START+1)*12))
YEARS_ANNUAL = Vector(YEAR_START:YEAR_END)

# %% Make Age Adjusted Plot
gr()
fig = plot()
heatmap!(fig, MONTHS_MONTHLY, 1:length(filt_ISO_list), age_penalties,
        yticks = (1:length(filt_ISO_list),filt_ISO_list),
        xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xrotation =90,
        xlabel = "Year", colorbartitle = "Age Adjusted Penalty", colormap = cgrad(["#1E86C2","#CE3648"]),
        clims = (0,1))
savefig(fig, "age_adjusted_penalty.pdf")
# %% Make Full Adjusted Plot
fig = plot()
heatmap!(fig,MONTHS_MONTHLY, 1:length(filt_ISO_list), full_penalties,
        yticks = (1:length(filt_ISO_list),filt_ISO_list),
        xticks = (MONTHS_MONTHLY[1:12:end],YEARS_ANNUAL[1]:YEARS_ANNUAL[end]), xrotation =90,
        xlabel = "Year", colorbartitle = "Full Adjusted Penalty", colormap = cgrad(["#1E86C2","#CE3648"]),
        clims = (0,1))
savefig(fig, "full_adjusted_penalty.pdf")

























