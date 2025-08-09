"""
Author: Eugene Tan
Date Created: 6/8/2025
Last Updated: 6/8/2025
Pre-processing of the distribution data to get a standardised dataset. Combines the following
- WHO Reported distributions (Gold standard, although it's kind bad but oh well)
- Tas' cleaned up dataset for the year 2024

Additional notes
- Aggregated up to Annual resolution for the purposes of MITN regression
- Aggregated up to Admin0 and Admin1 respectively
- Mauricio was in charge of type disaggregation. Disaggregation in to Admin1 is by population

Limited ITN Types. Current final roster
- cITN
- LLIN
- PBO
- DAI (ROYAL + G2 + Others)
"""

####################################
# %% Prep Environment and Load Packages
####################################

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from TOML file
include(pwd()*"/scripts/read_toml.jl")

# %% Data Wrangling
using CSV
using DataFrames
using ProgressBars

# %% Mathematics packages
using LinearAlgebra
using StatsBase
using Distributions

# %% Custom modules
using DateConversions

# %% Shapefile Processing
using GeoIO

####################################
# %% Load Datasets
####################################
# National Distribution 2000-2023
nat_dist_multitype_WHO_raw = CSV.read(RAW_DATASET_DIR*"net_distributions_cITN_adjusted_multitype.csv", DataFrame)

# Subnational Distribution 2000-2023
subnat_dist_multitype_WHO_raw = CSV.read(RAW_SUBNAT_DATASET_DIR*"net_distributions_admin1_dummy_combined_amp_multitype.csv", DataFrame)

# Admin2 Tas Distribution Data
admin2_dist_multitype_TAS_raw = CSV.read("E:/userdata/tas/scratch/ITN_volumes_20002024.csv", DataFrame)

####################################
# %% Clean WHO datasets
####################################
# Aggregate DAI = G2 + ROYAL in WHO Datasets
nat_dist_multitype_WHO = hcat(nat_dist_multitype_WHO_raw[:,["WHO_region","country","iso","year","Total Nets","cITN","LLIN","PBO"]],DataFrame(sum(Matrix(nat_dist_multitype_WHO_raw[:,["G2","ROYAL"]]), dims = 2),["DAI"]))
subnat_dist_multitype_WHO = hcat(subnat_dist_multitype_WHO_raw[:,["WHO_region","country","iso","admin1","admin1_id","area_id","year","type","Total Nets","cITN","LLIN","PBO"]],DataFrame(sum(Matrix(subnat_dist_multitype_WHO_raw[:,["G2","ROYAL"]]), dims = 2),["DAI"]))

####################################
# %% Get Unique Admin0 and Admin1 Metadata
####################################
# Scrape admin0-admin1 labels
admin0_metadata_rows = [nat_dist_multitype_WHO[i,["WHO_region", "country", "iso"]] for i in 1:size(nat_dist_multitype_WHO)[1]]
admin0_metadata = vcat(DataFrame.(unique(admin0_metadata_rows))...)
admin0_isos = unique(admin0_metadata.iso)

admin1_metadata_rows = [subnat_dist_multitype_WHO[i,["WHO_region", "country", "iso", "admin1", "admin1_id", "area_id"]] for i in 1:size(subnat_dist_multitype_WHO)[1]]
admin1_metadata = vcat(DataFrame.(unique(admin1_metadata_rows))...)
admin1_ids = unique(admin1_metadata.admin1_id)

####################################
# %% Clean Tas Datasets
####################################
# Rename dual AI and LLIN headings in dataset
admin2_dist_multitype_TAS_raw.type[findall(admin2_dist_multitype_TAS_raw.year .< 2010)] .= "cITN"
admin2_dist_multitype_TAS_raw.type[findall(admin2_dist_multitype_TAS_raw.type .== "dual_AI")] .= "DAI"

# Get list of unique names
TAS_net_names = unique(admin2_dist_multitype_TAS_raw.type)
TAS_years = sort(unique(admin2_dist_multitype_TAS_raw.year))

# Storage variables
nat_df_rows = []
subnat_df_rows = []

# Do joins for national 
for admin0_iso_i in ProgressBar(1:length(admin0_isos), leave = false)
    
    iso = admin0_isos[admin0_iso_i]

    for year in TAS_years

        volumes = zeros(Int,length(TAS_net_names))

        # Aggregate and extract volumes for each net type
        for net_type_i in 1:length(TAS_net_names)
            type = TAS_net_names[net_type_i]
            admin2_filt_data = admin2_dist_multitype_TAS_raw[(admin2_dist_multitype_TAS_raw.year .== year) .&& 
                                                                (admin2_dist_multitype_TAS_raw.ISO3 .== iso) .&& 
                                                                (admin2_dist_multitype_TAS_raw.type .== type),:]
            volumes[net_type_i] = round(Int,sum(admin2_filt_data.volume))
        end

        # Construct Dataframe Entry
        df_entry = hcat(DataFrame(admin0_metadata[admin0_iso_i,:]),DataFrame(year = year),DataFrame(Matrix(vcat(sum(Int,volumes), volumes)'), vcat("Total Nets",TAS_net_names)))

        # Save to storage variable
        push!(nat_df_rows, df_entry)
    end
end

# Do joins for subnational 
for admin1_id_i in ProgressBar(1:length(admin1_ids), leave = false)
    
    id = admin1_ids[admin1_id_i]

    for year in TAS_years

        volumes = zeros(Int,length(TAS_net_names))

        # Aggregate and extract volumes for each net type
        for net_type_i in 1:length(TAS_net_names)
            type = TAS_net_names[net_type_i]
            admin2_filt_data = admin2_dist_multitype_TAS_raw[(admin2_dist_multitype_TAS_raw.year .== year) .&& 
                                                                (admin2_dist_multitype_TAS_raw.ID_1 .== id) .&& 
                                                                (admin2_dist_multitype_TAS_raw.type .== type),:]
            volumes[net_type_i] = round(Int,sum(admin2_filt_data.volume))
        end

        # Construct Dataframe Entry
        df_entry = hcat(DataFrame(admin1_metadata[admin1_id_i,:]),DataFrame(year = year, type = "TAS"),DataFrame(Matrix(vcat(sum(Int,volumes), volumes)'), vcat("Total Nets",TAS_net_names)))

        # Save to storage variable
        push!(subnat_df_rows, df_entry)
    end
end

# Join all rows together
nat_dist_multitype_TAS = vcat(nat_df_rows...)
subnat_dist_multitype_TAS = vcat(subnat_df_rows...)

# Sort datasets
sort!(nat_dist_multitype_TAS, [order(:iso), order(:year)])
sort!(subnat_dist_multitype_TAS, [order(:iso), order(:admin1), order(:year)])

####################################
# %% Join WHO and TAS datasets
####################################

# Combine datasets such that WHO for <2024, TAS for 2024
nat_dist_multitype_FINAL = vcat(nat_dist_multitype_WHO,nat_dist_multitype_TAS[findall(nat_dist_multitype_TAS.year .== 2024),:])
subnat_dist_multitype_FINAL = vcat(subnat_dist_multitype_WHO,subnat_dist_multitype_TAS[findall(subnat_dist_multitype_TAS.year .== 2024),:])

# Sort datasets
sort!(nat_dist_multitype_FINAL, [order(:iso), order(:year)])
sort!(subnat_dist_multitype_FINAL, [order(:iso), order(:admin1), order(:year)])

# Create LLIN only version of dataset
nat_dist_singletype_FINAL = hcat(nat_dist_multitype_FINAL[:,["WHO_region","country","iso","year","Total Nets","cITN"]],DataFrame(sum(Matrix(nat_dist_multitype_FINAL[:,["LLIN","PBO","DAI"]]), dims = 2), ["LLIN"]))
subnat_dist_singletype_FINAL = hcat(subnat_dist_multitype_FINAL[:,["WHO_region","country","iso","admin1","admin1_id","area_id","year","type","Total Nets","cITN"]],DataFrame(sum(Matrix(subnat_dist_multitype_FINAL[:,["LLIN","PBO","DAI"]]), dims = 2), ["LLIN"]))

####################################
# %% Save CSV datasets
####################################
CSV.write("datasets/NAT_net_distributions_final_multitype.csv", nat_dist_multitype_FINAL)
CSV.write("datasets/subnational/SUBNAT_net_distributions_final_multitype.csv", subnat_dist_multitype_FINAL)

CSV.write("datasets/NAT_net_distributions_final_singletype.csv", nat_dist_singletype_FINAL)
CSV.write("datasets/subnational/SUBNAT_net_distributions_final_singletype.csv", subnat_dist_singletype_FINAL)

####################################
# %% Append distribution values to deliveries dataset (Temporary measure)
####################################
delivery_data = CSV.read("datasets/base_manufacturer_deliveries.csv", DataFrame)

ISOs = delivery_data.ISO3
deliv_values = zeros(Int,length(ISOs))

for i in 1:length(ISOs)
    println(i)

    ISO = ISOs[i]

    value = nat_dist_singletype_FINAL[(nat_dist_singletype_FINAL.year .== 2024) .&& 
                                (nat_dist_multitype_FINAL.iso .== ISO),"Total Nets"]
    if !isempty(value)
        deliv_values[i] = value[1]
    end
end

CSV.write("datasets/net_deliveries_final.csv",hcat(delivery_data, DataFrame("2024" => deliv_values)))
