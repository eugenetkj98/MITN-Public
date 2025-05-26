"""
Author: Eugene Tan
Date Created: 23/4/2025
Last Updated: 23/4/2025
Use raster level estimates to reconstruct survey aggregates that would have been used to
calibrate SNF model
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from TOML file
include(pwd()*"/scripts/read_toml.jl")

# %% Import relevant packages
using ProgressBars
using DataFrames
using CSV

# %% Import data lookup dataset
data = CSV.read(OUTPUT_DIR*"coverage_timeseries/bv_mitn_validation_values.csv", DataFrame)

# %% Get list of all unique combinations of ISO, interview month and year
metadata_combinations = unique(data[:,["ISO", "interview_month", "interview_year"]])
sort!(metadata_combinations, [order(:interview_year), order(:ISO), order(:interview_month)])

# %% Construct observation aggregates for each ISO region
df_entries = []

for i in 1:size(metadata_combinations)[1]
    # Extract metadata
    ISO, month, year = metadata_combinations[i,:]

    # Filter to required slice
    data_slice = data[data.ISO .== ISO .&&
                        data.interview_month .== month .&&
                        data.interview_year .== year,:]

    # Calculate NPC, Access and Use Values
    weights = data_slice.cluster_sample_wt
    npc_obs = data_slice.npc
    access_obs = data_slice.access
    use_obs = data_slice.use
    bv_npc_obs = data_slice.bv_npc
    bv_access_obs = data_slice.bv_access
    bv_use_obs = data_slice.bv_use
    mitn_npc_obs = data_slice.mitn_npc
    mitn_access_obs = data_slice.mitn_access
    mitn_use_obs = data_slice.mitn_use

    npc_val = sum(weights.*npc_obs)/sum(weights)
    access_val = sum(weights.*access_obs)/sum(weights)
    use_val = sum(weights.*use_obs)/sum(weights)

    bv_npc_val = sum(weights.*bv_npc_obs)/sum(weights)
    bv_access_val = sum(weights.*bv_access_obs)/sum(weights)
    bv_use_val = sum(weights.*bv_use_obs)/sum(weights)

    mitn_npc_val = sum(weights.*mitn_npc_obs)/sum(weights)
    mitn_access_val = sum(weights.*mitn_access_obs)/sum(weights)
    mitn_use_val = sum(weights.*mitn_use_obs)/sum(weights)

    # Construct data frame slice
    df_entry = DataFrame(ISO = ISO,
                            interview_month = month,
                            interview_year = year,
                            npc = npc_val,
                            bv_npc = bv_npc_val,
                            mitn_npc = mitn_npc_val,
                            access = access_val,
                            bv_access = bv_access_val,
                            mitn_access = mitn_access_val,
                            use = use_val,
                            bv_use = bv_use_val,
                            mitn_use = mitn_use_val)
    push!(df_entries, df_entry)
end

# %%
reconstructed_surveys = vcat(df_entries...)

# %% Save into outputs
CSV.write(OUTPUT_DIR*"coverage_timeseries/bv_mitn_survey_reconstructions.csv", reconstructed_surveys)