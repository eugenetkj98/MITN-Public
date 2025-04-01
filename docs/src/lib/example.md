# Example Code

Some example code to run net crop and access predictions of the SNF model for a single country.

Extract relevant data for net crop regression for a single country.

```julia
# Define analysis parameters
ISO = "SEN"
YEAR_START = 2010
YEAR_END = 2021

# Extract data
extract_data_netcrop(ISO, YEAR_START, YEAR_END)

# By default, output is saved into location "outputs/extractions/crop/"
```

Perform MCMC + gradient descent regression to get required net crop parameters

```julia
# Load saved data extraction for net crop
input_dict = load("outputs/extractions/crop/SEN_2010_2021_cropextract.jld2")

# Perform Bayesian inference + gradient descent algorithm and save the output
# Default save location is "outputs/regressions/crop/"
bayes_GD(input_dict, save_output = true, N_EPOCHS = 20)
```

Extract survey data for calculating net access in later steps.

```julia
# Load saved data for net crop
input_dict = load("outputs/extractions/crop/SEN_2010_2021_cropextract.jld2")
regression_dict = load("outputs/regressions/crop/SEN_2010_2021_cropchains.jld2")

# Extract required survey data for net access. 
# Default save location is "outputs/extractions/access/"
extract_data_netaccess(ISO, YEAR_START, YEAR_END, input_dict, regression_dict)
```

Extract and aggregate all the survey estimates for net access parameters.

**Note: this only needs to be done once and requires access to the 'regression_dict' files for all countries due to sparsity of data.**

```julia
# Aggregate survey data for regressing net access parameters and save .jld2
init_variables = false

# Pre-define storage variables. Dimensions are not important here and will be defined later.
p_h_globaldata = zeros(0,10)
ρ_h_globaldata = zeros(0,10)
μ_h_globaldata = zeros(0,10)
γ_globaldata = zeros(0)

# List of countries to exclude due to errors (e.g. NaNs, data format etc.)
exclusion_ISOs = ["CPV","BWA","CAF","COM","GNQ","DJI","ERI","ETH","GAB","GNB",
                    "STP","SOM","SDN","SWZ","ZAF","SSD"]

for i in 1:length(ISO_list)
    # Select ISO
    ISO = ISO_list[i]

    if ISO ∈ exclusion_ISOs # Choose which countries to exclude in extraction.
        continue
    else
        # Import extracted data
        net_access_input_dict = load("outputs/extractions/access/$(ISO)_2010_2021_accessextract.jld2")
        
        γ_surveys = net_access_input_dict["γ_aggregated"]
        p_h = net_access_input_dict["p_h_aggregated"]
        ρ_h = net_access_input_dict["ρ_h_aggregated"]
        μ_h = net_access_input_dict["μ_h_aggregated"]
        
        h_max = size(p_h)[2]
        # Check if storage variable for aggregated global data has been made.
        if init_variables == false
            p_h_globaldata = zeros(0,h_max)
            ρ_h_globaldata = zeros(0,h_max)
            μ_h_globaldata = zeros(0,h_max)
            γ_globaldata = zeros(0)

            init_variables = true
        end

        # Add to aggregate storage variable
        p_h_globaldata = vcat(p_h_globaldata, p_h)
        ρ_h_globaldata = vcat(ρ_h_globaldata, ρ_h)
        μ_h_globaldata = vcat(μ_h_globaldata, μ_h)
        γ_globaldata = vcat(γ_globaldata, γ_surveys)
    end
end

# Save aggregated data
jldsave("outputs/extractions/access/aggregated_inputs/netaccess_allsurvey_inputs.jld2";
            p_h_globaldata, ρ_h_globaldata, μ_h_globaldata, γ_globaldata)
```

Perform MCMC regression for net access model

```julia
# MCMC Regression for Net Access model
# Default save location is "output/extractions/access/aggregated_inputs/netaccesschains.jld2"
bayes_access(access_survey_globaldata;
                iterations = 10000, burn_in = 2000,
                chain_output_dir = chain_output_dir)
println("Analysis Complete. Check directory for outputs.")
```

Feed results into plotting functions for visualisation if required. 