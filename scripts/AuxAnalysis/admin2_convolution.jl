"""
Author: Eugene Tan
Date Created: 10/6/2025
Last Updated: 10/6/2025
Shortcut way for doing NPC predictions at Admin2 level, based on distributions time series and learned SNF attrition parameters.
    Shortcut is based on convolutions
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/read_toml.jl")

# %% Import Public Packages
using DataFrames
using Missings
using JLD2
using CSV
using ProgressBars

# %% Maths packages
using LinearAlgebra
using StatsBase

# %% General useful functions
using DateConversions
using NetLoss
using NetAccessPrediction
using Convolutions

# %% Load Admin2 Distribution Dataset
admin2_filepath = "/mnt/efs/model_dev/meerkat/data/gfatm_request_monthly_volume_02072025.csv"
admin2_data = CSV.read(admin2_filepath, DataFrame)

# %% Define save file name
input_filename = split(admin2_filepath, "/")[end]
save_dir = "/mnt/efs/model_dev/meerkat/data/admin2_predictions/"
save_filename = split(input_filename,".")[1]*"_pred.csv"


# %% Extract Required Country data
# List of ISOs in input dataset
data_ISOs = unique(admin2_data.ISO3)

# List of ISOs in our library
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# Find overlap of ISOs bewteen input dataset and our library
ISOs = intersect(data_ISOs, filt_ISOs)

# Get reference year to start from
YEAR_START = minimum(admin2_data.year)

# Subnational data file location
subnat_reg_dir = OUTPUT_REGRESSIONS_DIR*"subnational/"
REG_YEAR_START_NAT = YEAR_NAT_START # Start year for national model
REG_YEAR_START = YEAR_SUBNAT_TRANS #2011 # Start year for subnational model
REG_YEAR_END = YEAR_NAT_END

# Enumerate number of scenarios
scenario_names = names(admin2_data)[11:end]
n_scenario = length(scenario_names)

# Length of time (in years) to propagate impulse response
loss_t_vals = 0:1/12:30

# %% Storage variable for outputs
df_collection = Vector{DataFrame}(undef, length(ISOs))
# %% Do prediction for each ISO

for ISO_i in ProgressBar(1:length(ISOs))
    ISO = ISOs[ISO_i]

    # Import metadata regression outputs for parameters
    nat_reg = JLD2.load(OUTPUT_REGRESSIONS_DIR*"crop/$(REG_YEAR_START_NAT)_$(REG_YEAR_END)/$(ISO)_$(REG_YEAR_START_NAT)_$(REG_YEAR_END)_cropchains.jld2")
    reg_admin0_τ = mean(nat_reg["chain"][:,5])
    reg_admin0_κ = mean(nat_reg["chain"][:,6])

    net_access_input_dict = load(OUTPUT_EXTRACTIONS_DIR*"access/pred_data/$(REG_YEAR_START_NAT)_$(REG_YEAR_END)/$(ISO)_$(REG_YEAR_START_NAT)_$(REG_YEAR_END)_accessextract.jld2")
    net_access_chain = load(OUTPUT_REGRESSIONS_DIR*"access/netaccesschains.jld2")
    μ_chain_df = net_access_chain["μ_chain_df"]
    ρ_chain_df = net_access_chain["ρ_chain_df"]
    p_h_mean = mean(net_access_input_dict["p_h_aggregated"], dims = 1)[:]
    
    subnat_reg_filename = "$(ISO)_SUBNAT_NETCROP_$(REG_YEAR_START_NAT)_$(REG_YEAR_START)_$(REG_YEAR_END)_regression.jld2"
    subnat_reg = JLD2.load(subnat_reg_dir*subnat_reg_filename)

    reg_admin1_names = subnat_reg["admin1_names"]

    reg_admin1_τ = [subnat_reg["admin1_outputs"][i]["τ_est"][2] for i in 1:length(subnat_reg["admin1_outputs"])]
    reg_admin1_κ = [subnat_reg["admin1_outputs"][i]["κ_est"][2] for i in 1:length(subnat_reg["admin1_outputs"])]

    # Get admin1 ids - Need to do in layers because some countries have bad data
    admin1_ids = unique(admin2_data[admin2_data.ISO3 .== ISO,"ID_1"])
    admin1_df = Vector{DataFrame}(undef, length(admin1_ids))
    
    Threads.@threads for admin1_id_i in ProgressBar(1:length(admin1_ids), leave = false)
        admin1_id = admin1_ids[admin1_id_i]
        # Get admin2 ids
        admin2_ids = unique(admin2_data[(admin2_data.ISO3 .== ISO) .* (admin2_data.ID_1 .== admin1_ids[admin1_id_i]),"ID_2"])
        admin2_df = Vector{DataFrame}(undef, length(admin2_ids))

        for admin2_id_i in 1:length(admin2_ids)
            # Select admin2 id and filter data
            admin2_id = admin2_ids[admin2_id_i]
            filt_data = admin2_data[(admin2_data.ID_1 .== admin1_id).*(admin2_data.ID_2 .== admin2_id),:]
            sort!(filt_data , [order(:year), order(:month)])

            # Get time indexes
            monthidxs = monthyear_to_monthidx.(filt_data.month, filt_data.year, YEAR_START = YEAR_START)

            # Create storage variables
            admin2_dist = Matrix{Float64}(undef, n_scenario, maximum(monthidxs))
            admin2_population = Vector{Float64}(undef, maximum(monthidxs))
            admin2_crop = Matrix{Float64}(undef, n_scenario, maximum(monthidxs))
            admin2_access = Matrix{Float64}(undef, n_scenario, maximum(monthidxs))

            # Extract distribution data
            for i in 1:length(monthidxs)
                # Check for NA
                na_idxs = findall(Array(filt_data[i,scenario_names]) .== "NA")

                if !isempty(na_idxs)
                    filt_data[i,scenario_names][na_idxs] .= "0"
                end
                
                # Parse data
                admin2_dist[:,monthidxs[i]] = Vector(filt_data[i,scenario_names])
                admin2_population[monthidxs[i]] = filt_data[i,"population"]
            end

            # Determine which attrition parameter to use. If it is one of the regressed admin1 use that, otherwise use national SNF
            τ = reg_admin0_τ
            κ = reg_admin0_κ
            if filt_data[1,"Name_1"] ∈ reg_admin1_names
                admin1_idx = findfirst(reg_admin1_names .== filt_data[1,"Name_1"])
                τ = reg_admin1_τ[admin1_idx]
                κ = reg_admin1_κ[admin1_idx]
            end

            # Convolve to get net crop
            L = net_loss_compact.(loss_t_vals, τ, κ)
            for scenario_i in 1:n_scenario
                admin2_crop[scenario_i,:] = convolution(admin2_dist[scenario_i,:],L)[1:size(admin2_dist)[2]]
            end

            # Calculate npc and access
            admin2_npc = admin2_crop./repeat(admin2_population', n_scenario,1)
            for scenario_i in 1:n_scenario
                admin2_access[scenario_i,:] = mean_net_access(Matrix(ρ_chain_df), Matrix(μ_chain_df), p_h_mean,
                                                                admin2_population, admin2_crop[scenario_i,:])
            end            

            df1 = DataFrame(ISO3 = ISO,
                        Name_0 = filt_data.Name_0,
                        ID_0 = filt_data.ID_0,
                        Name_1 = filt_data.Name_1,
                        ID_1 = filt_data.ID_1,
                        Name_2 = filt_data.Name_2,
                        ID_2 = filt_data.ID_2,
                        year = filt_data.year,
                        month = filt_data.month,
                        population = filt_data.population)

            df2 = filt_data[:,11:end]

            df3 = DataFrame(admin2_crop', "crop_pred_".*scenario_names)

            df4 = DataFrame(admin2_npc', "npc_pred_".*scenario_names)

            df5 = DataFrame(admin2_access', "access_pred_".*scenario_names)

            admin2_df[admin2_id_i] = hcat(df1, df2, df3, df4, df5)
        end

        admin1_df[admin1_id_i] = vcat(admin2_df...)
    end

    # Concatenate data frame for each admin2 region and save to storage variable
    df_collection[ISO_i] = vcat(admin1_df...)

end

# %%
# Concatenate data frames and sort + format
master_df = vcat(df_collection...)
sort!(master_df, [order(:ISO3), order(:Name_1), order(:Name_2), order(:year), order(:month)])

# Define save dir and remove old output to overwrite
mkpath(save_dir)
# rm(save_dir*save_filename)

# Write file
CSV.write(save_dir*save_filename, master_df)
println("Saved time series admin2 predictions at: $(save_dir*save_filename)")

