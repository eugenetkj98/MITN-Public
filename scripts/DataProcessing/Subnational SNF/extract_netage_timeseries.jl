"""
Author: Eugene Tan
Date Created: 8/5/2025
Last Updated: 8/5/2025
Postprocess SNF outputs to get data for mean net age by country
"""
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Import Public Packages
using DataFrames
using JLD2
using CSV
using ProgressBars
using GeoIO
using Rasters

# %% Maths packages
using LinearAlgebra
using StatsBase

# %% Useful functions
using DateConversions

# %% Plotting Packages
using CairoMakie

# %% Years
YEAR_START = YEAR_NAT_START
YEAR_END = YEAR_NAT_END
n_months = 12*(YEAR_END-YEAR_START+1)

# %%
input_dir = OUTPUT_DRAWS_DIR*"subnational/"
output_dir = OUTPUT_DATAPREP_DIR
output_filename = "snf_mean_netage.csv"

# %% Shapefile metadata (just to get Admin0 names)
admin0_shapes_geoIO = GeoIO.load(ADMIN0_SHAPEFILE)

# %% Define list of countries to plot
ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
exclusion_ISOs = EXCLUSION_ISOS
filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

# %% Extract age timeseries for each country and save as DataFrane
df1_collection = []
df0_collection = []
for ISO_i in ProgressBar(1:length(filt_ISOs))
    # Select country
    ISO = filt_ISOs[ISO_i]

    # Get Admin0 Name
    admin0_name = admin0_shapes_geoIO[findfirst(admin0_shapes_geoIO.ISO .== ISO),"Name_0"]

    # Import SNF data
    snf_posterior = JLD2.load(input_dir*"$(ISO)_SUBNAT_draws.jld2")

    # Get admin1 metadata
    admin1_names = snf_posterior["admin1_names"]
    n_admin1 = length(admin1_names)

    # Demography matrix collection to get national estimate
    A_collection = Vector{Matrix}(undef, n_admin1)

    # Extract age data for each admin1 region
    for admin1_i in 1:n_admin1
        # Get area id
        area_id = snf_posterior["merged_outputs"][admin1_i]["area_id"]

        # Get Demography Matrix
        A = snf_posterior["merged_outputs"][admin1_i]["ADJ_COMBINED_A_TOTAL_mean"]
        A_collection[admin1_i] = A

        # Calculate Age Weight Matrix
        n_months = size(A)[2]
        M = zeros(n_months, n_months)
        for i in 1:n_months
            for j in 1:i
                M[i,j] = i-j
            end
        end

        # Calculate average age in months
        mean_age = (sum(M.*A, dims = 2)[:,1])./(sum(A, dims = 2)[:,1])
        mean_age[findall(isnan.(mean_age))] .= 0

        # Calculate month and year values
        month_vals = [monthidx_to_monthyear(i)[1] for i in 1:n_months]
        year_vals = [monthidx_to_monthyear(i)[2]-1+YEAR_START for i in 1:n_months]

        # Construct and store dataframe
        df_entry = DataFrame(ISO = ISO,
                        name = admin1_names[admin1_i],
                        category = "Admin1",
                        area_id = area_id,
                        month = month_vals, year = year_vals,
                        mean_age_months = mean_age)

        push!(df1_collection, df_entry)
    end

    # Calculation National mean net age estimate
    nat_A = sum(A_collection)

    # Calculate Age Weight Matrix
    n_months = size(nat_A)[2]
    M = zeros(n_months, n_months)
    for i in 1:n_months
        for j in 1:i
            M[i,j] = i-j
        end
    end

    # Calculate average age in months
    mean_age = (sum(M.*A, dims = 2)[:,1])./(sum(A, dims = 2)[:,1])
    mean_age[findall(isnan.(mean_age))] .= 0

    # Calculate month and year values
    month_vals = [monthidx_to_monthyear(i)[1] for i in 1:n_months]
    year_vals = [monthidx_to_monthyear(i)[2]-1+YEAR_START for i in 1:n_months]

    # Construct and store dataframe
    df_entry = DataFrame(ISO = ISO,
                    name = admin0_name,
                    category = "Admin0",
                    area_id = area_id,
                    month = month_vals, year = year_vals,
                    mean_age_months = mean_age)

    push!(df0_collection, df_entry)
end

# %% Construct data set and save
output_df = vcat(df0_collection..., df1_collection...)
CSV.write(output_dir*output_filename, output_df)