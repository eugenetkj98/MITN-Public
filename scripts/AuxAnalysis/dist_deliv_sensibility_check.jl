"""
Author: Eugene Tan
Date Created: 12/8/2025
Last Updated: 12/8/2025
Script to check raw national delivery and distribution data on per capita basis to check that it's reasonable
- Produces a CSV of compiled data with flags
- Makes plots
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Prep environment and subdirectories
include(pwd()*"/scripts/read_toml.jl")

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
# %% Storage variable for DataFrame entries
###########################################
df_entries = []

###########################################
# %% Extract data and calculate per capita values for distribution. Saves to CSV file
###########################################

for ISO in filt_ISOs
    ###########################################
    # %% Simulate SNF for subnational posterior estimate
    ###########################################

    # Load required JLD2 data files
    nat_extractions = JLD2.load(OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_NAT_START)_$(YEAR_NAT_END)/$(ISO)_$(YEAR_NAT_START)_$(YEAR_NAT_END)_cropextract.jld2")

    NET_NAMES = nat_extractions["NET_NAMES"]
    YEARS_ANNUAL = nat_extractions["YEARS_ANNUAL"]
    POPULATION_ANNUAL = nat_extractions["POPULATION_ANNUAL"][1:end-1]
    DELIVERY_ANNUAL = nat_extractions["DELIVERIES_ANNUAL"]
    DISTRIBUTION_ANNUAL_TOTAL =nat_extractions["DISTRIBUTION_ANNUAL"][:,1]
    DISTRIBUTION_ANNUAL_BYNET = nat_extractions["DISTRIBUTION_ANNUAL"][:,2:end]

    DELIVERY_NPC_ANNUAL = DELIVERY_ANNUAL./POPULATION_ANNUAL
    DISTRIBUTION_NPC_ANNUAL_BYNET = DISTRIBUTION_ANNUAL_BYNET./repeat(POPULATION_ANNUAL, 1, size(DISTRIBUTION_ANNUAL_BYNET)[2])
    DISTRIBUTION_NPC_ANNUAL_TOTAL = sum(DISTRIBUTION_NPC_ANNUAL_BYNET, dims = 2)[:]

    DELIVERY_CHECK = Vector{Union{Missing, String}}(undef, length(DELIVERY_NPC_ANNUAL))
    DISTRIBUTION_CHECK = Vector{Union{Missing, String}}(undef, length(DISTRIBUTION_NPC_ANNUAL_TOTAL))

    for i in 1:length(DELIVERY_NPC_ANNUAL)
        if ismissing(DELIVERY_NPC_ANNUAL[i])
            DELIVERY_CHECK[i] = "MISSING"
        elseif DELIVERY_NPC_ANNUAL[i] > 0.5
            DELIVERY_CHECK[i] = ">0.5"
        else
            DELIVERY_CHECK[i] = "<=0.5"
        end
    end

    for i in 1:length(DISTRIBUTION_NPC_ANNUAL_TOTAL)
        if ismissing(DISTRIBUTION_NPC_ANNUAL_TOTAL[i])
            DISTRIBUTION_CHECK[i] = "MISSING"
        elseif DISTRIBUTION_NPC_ANNUAL_TOTAL[i] > 0.5
            DISTRIBUTION_CHECK[i] = ">0.5"
        else
            DISTRIBUTION_CHECK[i] = "<=0.5"
        end
    end

    df_entry = hcat(   DataFrame(ISO = "$(ISO)", year = YEARS_ANNUAL,
                    LLIN_DELIVERY = DELIVERY_ANNUAL,
                    LLIN_DELIVERY_NPC = DELIVERY_NPC_ANNUAL,
                    DELIVERY_FLAG = DELIVERY_CHECK,
                    ITN_DISTRIBUTION = DISTRIBUTION_ANNUAL_TOTAL,
                    ITN_DISTRIBUTION_NPC = DISTRIBUTION_NPC_ANNUAL_TOTAL,
                    DISTRIBUTION_FLAG = DISTRIBUTION_CHECK),
            DataFrame(DISTRIBUTION_ANNUAL_BYNET, NET_NAMES.*"_DISTRIBUTION"),
            DataFrame(DISTRIBUTION_NPC_ANNUAL_BYNET, NET_NAMES.*"_NPC")
            )

    push!(df_entries, df_entry)
end

# Concatenate dataframe
combined_df = vcat(df_entries...)

# Save
CSV.write(OUTPUT_DATAPREP_DIR*"distribution_datacheck.csv", combined_df)
