"""
Author: Eugene Tan
Date Created: 27/3/2025
Last Updated: 27/3/2025
Script to explore and check relationship between adjusted SE values original SE.
Used to provide metric estimates for datapoints taken from report summaries
"""

# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Data Wrangling
using CSV
using DataFrames

# %% Mathematics packages
using LinearAlgebra
using StatsBase
using Distributions

# %% Plotters
using Plots

# %%
npc_monthly_data = CSV.read("datasets/npc_monthly_data.csv", DataFrame)

# %%
scatter(npc_monthly_data.NPC_unadj_se, npc_monthly_data.NPC_adj_se./npc_monthly_data.NPC_unadj_se, ylims = (0,10))

scatter(npc_monthly_data.NPC_mean, npc_monthly_data.NPC_mean./npc_monthly_data.NPC_adj_se, ylims = (0,10))

histogram(npc_monthly_data.NPC_adj_se./npc_monthly_data.NPC_unadj_se, xlims = (0,20), bins = 0:20)

quantile(npc_monthly_data.NPC_adj_se./npc_monthly_data.NPC_unadj_se, [0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99])

quantile(npc_monthly_data.NPC_adj_se, [0, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99])

