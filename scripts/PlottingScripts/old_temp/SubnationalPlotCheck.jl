"""
Author: Eugene Tan
Date Created: 4/11/2024
Last Updated: 11/11/2024
Script for plotting subnational NPC and access gap across regions
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %%
using Plots
using ProgressBars
using StatsBase
using CSV
using DataFrames
using Missings
using JLD2
using DateConversions

# %%
ISO = "MOZ"
YEAR_START = 2000
YEAR_END = 2023

# %% Import datasets
# National Estimates
nat_dataset = JLD2.load("outputs/draws/national/$(ISO)_$(YEAR_START)_$(YEAR_END)_post_crop_access.jld2")
# Subnational Estimates
subnat_dataset = JLD2.load("outputs/draws/subnational/$(ISO)_SUBNAT_draws.jld2")

# %% Get number of regions
n_admin1 = length(subnat_dataset["admin1_names"])

# %% Extract gap values
NAT_A_NPC_mean_BYNET = nat_dataset["A_NPC_mean_BYNET"]
NAT_NPC_mean = sum(NAT_A_NPC_mean_BYNET, dims = (2,3))[:]
NAT_λ_access_mean = nat_dataset["λ_access_mean"]

SUBNAT_NPC_mean = zeros(n_admin1, length(NAT_NPC_mean))
SUBNAT_λ_access_mean = zeros(n_admin1, length(NAT_λ_access_mean))
subnat_dataset["admin1_outputs"][i]
for i in 1:n_admin1
    SUBNAT_NPC_mean[i,:] = mean(subnat_dataset["admin1_outputs"][i]["NPC_MONTHLY_TOTAL_samples"], dims = 1)[:]
    SUBNAT_λ_access_mean[i,:] = mean(subnat_dataset["admin1_outputs"][i]["λ_ACCESS_samples"], dims = 1)[:]
end

# %%
dataset_dir = "datasets/subnational/"
# NPC monthly dataset for regression
subnat_npc_monthly_data_filename = "subnat_npc_monthly_data.csv"
subnat_npc_monthly_data = CSV.read(dataset_dir*subnat_npc_monthly_data_filename, DataFrame)

# %%
id_legend_filename = "admin2023_1_MG_5K_config.csv"
master_id_legend = CSV.read(dataset_dir*id_legend_filename, DataFrame)

# %%
population_filename = "map_admin1_zonal_stats_pops.csv"
master_populations = CSV.read(dataset_dir*population_filename, DataFrame)


# %%
YEAR_START_NAT = subnat_dataset["YEAR_START_NAT"]
YEARS_ANNUAL = Vector(YEAR_START:1:(YEAR_END))
FULL_YEARS_ANNUAL = Vector(YEAR_START_NAT:1:(YEAR_END))
N_MONTHS = (YEAR_END-YEAR_START+1)*12
MONTHS_MONTHLY = Vector(1:(YEAR_END-YEAR_START+1)*12)
FULL_MONTHS_MONTHLY = Vector(1:(YEAR_END-YEAR_START_NAT+1)*12)

# %% Population Data
admin1_master_populations = master_populations[(master_populations.area_id .== area_id) .&
                                                (master_populations.year .>= YEAR_START) .&
                                                (master_populations.year .<= YEAR_END+1),:]
FULL_admin1_master_populations = master_populations[(master_populations.area_id .== area_id) .&
                                                (master_populations.year .>= YEAR_START_NAT) .&
                                                (master_populations.year .<= YEAR_END+1),:]

POPULATION_ANNUAL = admin1_master_populations[sortperm(admin1_master_populations.year),"SUM"]
FULL_POPULATION_ANNUAL = FULL_admin1_master_populations[sortperm(FULL_admin1_master_populations.year),"SUM"]                                    

POPULATION_MONTHLY = zeros(length(MONTHS_MONTHLY))
FULL_POPULATION_MONTHLY = zeros(length(FULL_MONTHS_MONTHLY))

for i in 1:(length(YEARS_ANNUAL))
    POPULATION_MONTHLY[(12*(i-1)+1):(12*i)] = LinRange(POPULATION_ANNUAL[i], POPULATION_ANNUAL[i+1],13)[1:12]
end

for i in 1:(length(FULL_YEARS_ANNUAL))
    FULL_POPULATION_MONTHLY[(12*(i-1)+1):(12*i)] = LinRange(FULL_POPULATION_ANNUAL[i], FULL_POPULATION_ANNUAL[i+1],13)[1:12]
end
# %%

# Extract observed netcrop estimates for survey aggregate (subnat_npc_monthly_data.csv)
FULL_SUBNAT_NET_CROP_MONTHLY = missings(Float64,length(FULL_MONTHS_MONTHLY))
FULL_SUBNAT_NET_CROP_STD_MONTHLY = missings(Float64,length(FULL_MONTHS_MONTHLY))

admin1_i = 4
admin1_name = subnat_dataset["admin1_names"][admin1_i]
area_id = master_id_legend[master_id_legend.Name_1 .== admin1_name, "area_id"]

subnat_npc_entries = subnat_npc_monthly_data[subnat_npc_monthly_data.area_id .== area_id,:]

for row_i in 1:size(subnat_npc_entries)[1]
    month_val = subnat_npc_entries[row_i,"month"]
    year_val = subnat_npc_entries[row_i,"year"]
    monthidx = monthyear_to_monthidx(month_val, year_val, YEAR_START = YEAR_START_NAT)
    FULL_SUBNAT_NET_CROP_MONTHLY[monthidx] = subnat_npc_entries[row_i,"NPC_mean"]*FULL_POPULATION_MONTHLY[monthidx]
    FULL_SUBNAT_NET_CROP_STD_MONTHLY[monthidx] = subnat_npc_entries[row_i,"NPC_adj_se"]*FULL_POPULATION_MONTHLY[monthidx]
end

# %%
fig = plot()

plot!(fig, SUBNAT_NPC_mean[admin1_i,:].*FULL_POPULATION_MONTHLY)
scatter!(fig, FULL_SUBNAT_NET_CROP_MONTHLY)
plot!(fig, NAT_NPC_mean.*FULL_POPULATION_MONTHLY)