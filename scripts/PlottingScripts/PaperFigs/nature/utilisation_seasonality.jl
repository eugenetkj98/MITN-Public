"""
Author: Eugene Tan
Date Created: 8/5/2025
Last Updated: 8/5/2025
Analyse seasonality behaviour of utilisation
"""
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from TOML file
include(pwd()*"/scripts/read_toml.jl")

# %% Import Public Packages
using CSV
using DataFrames

# %% Import Custom Packages
using DateConversions

# Maths packages
using LinearAlgebra
using StatsBase
using ProgressBars

# %% Plot packages
using LaTeXStrings
using CairoMakie

# %% Data Directories
input_data_filepath = OUTPUT_DIR*"coverage_timeseries/master_extraction.csv"

# %%
timeseries_data = CSV.read(input_data_filepath, DataFrame)


# %% Country List
ISO_list = ISO_LIST
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %%
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END


# %% Extract efficiency data
n_months = 12*(YEAR_END - YEAR_START + 1)
n_years = YEAR_END-YEAR_START+1
eff_data = Matrix{Float64}(undef, length(filt_ISOs), n_months)
eff_stack = Array{Float64}(undef, length(filt_ISOs), n_years, 12)
for ISO_i in 1:length(filt_ISOs)
    ISO = filt_ISOs[ISO_i]

    # Filter data
    filt_data = timeseries_data[findall(timeseries_data.ISO .== ISO),:]

    # Calculate monthidxs
    month_vals = filt_data.month
    year_vals = filt_data.year
    monthidxs = monthyear_to_monthidx.(month_vals, year_vals, YEAR_START = YEAR_START)

    for i in 1:size(filt_data)[1]
        eff_data[ISO_i, monthidxs[i]] = filt_data[i,"raster_eff_mean"]
        eff_stack[ISO_i, year_vals[i]-YEAR_START+1, month_vals[i]] = filt_data[i,"raster_eff_mean"]
    end
end

# %%
ISO_i = 1
ISO = filt_ISOs[ISO_i]

# Filter data
filt_data = timeseries_data[findall(timeseries_data.ISO .== ISO),:]

# %% Calculate autocorrelation integral
ACI = Vector{Float64}(undef, length(filt_ISOs))
τ_max = 12# Max num of months to calculate correlation against

for ISO_i in 1:length(filt_ISOs)
    ACI[ISO_i] = mean(autocor(eff_data[ISO_i,:],0:1:τ_max, demean = true))
end
scatter(ACI)

# %%
ACI[argmax(ACI)]
# %% Calculate annual mean modes
mean_modes = mean(eff_stack, dims = 2)[:,1,:]

ISO_i = 39
fig = Figure(size = (600,400))
ax = Axis(fig[1,1])
ylims!(ax, 0.5, 1.05)

for i in 1:size(eff_stack)[2]
    lines!(ax, 1:12, eff_stack[ISO_i,i,:], alpha = 0.2, color = :blue)
    lines!(ax, 1:12, mean_modes[ISO_i,:], color = :red)
end

fig
findfirst(filt_ISOs .== "TZA")
# %%
lines(util_data[12,:])