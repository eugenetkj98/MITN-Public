"""
Author: Eugene Tan
Date Created: 2/4/2025
Last Updated: 2/4/2025
Script to run regression for net crop
"""
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

using Distributed

# %% Create worker processes
n_workers = N_WORKERS
addprocs(n_workers)
println("Starting National MITN Regression with $(nprocs()) workers, each with $(Threads.nthreads()) threads")
flush(stdout)


@everywhere begin
    # %% Prep environment and subdirectories
    include(pwd()*"/scripts/init_env.jl")

    # %% Import filenames and directories from config file
    include(pwd()*"/scripts/dir_configs.jl")

    # %% Multiprocessing Packages
    using Distributed

    # %% Import Public Packages
    using JLD2
    using CSV
    using DataFrames

    # %% Import Custom Packages
    using DataExtractions
    using DateConversions
    using NetCropModel
    using NetCropRegression

    # Maths packages
    using LinearAlgebra
    using StatsBase

    # %% Get ISO List
    ISO_list = String.(CSV.read(RAW_DATASET_DIR*ISO_LIST_FILENAME, DataFrame)[:,1])
    exclusion_ISOs = EXCLUSION_ISOS
    filt_ISOs = setdiff(ISO_list, exclusion_ISOs)

    # %% Run Analysis
    YEAR_START = YEAR_NAT_START
    YEAR_END = YEAR_NAT_END
end

# %%
for i in 1:length(filt_ISOs)
    # Select ISO
    ISO = filt_ISOs[i]

    println("Extracting Data for Country $(i) of $(length(filt_ISOs)) → $(ISO).")
    flush(stdout)
    
    # Net crop data extraction
    extract_data_netcrop(ISO, YEAR_START, YEAR_END)

    # %%
    println("Net Crop Data Extraction complete for $(ISO). Data saved")
    flush(stdout)
end

@sync @distributed for i in 1:length(filt_ISOs)
    # Select ISO
    ISO = filt_ISOs[i]

    println("Fitting model for Country $(i) of $(length(ISO_list)) → $(ISO), on worker $(Distributed.myid()) with $(Threads.nthreads())")
    flush(stdout)
    
    # if ISO ∈ exclusion_ISOs
    #     println("$(ISO) is on exclusion list. Moving to next country.")
    #     flush(stdout)
    #     continue
    # else

    # Load extracted data
    input_dict = load(OUTPUT_EXTRACTIONS_DIR*"crop/$(YEAR_START)_$(YEAR_END)/$(ISO)_$(YEAR_START)_$(YEAR_END)_cropextract.jld2")

    # Net crop regression
    bayes_GD(input_dict, save_output = true, N_EPOCHS = 5, n_chains = Threads.nthreads())

    # %%
    println("Net Crop Regression complete for $(ISO). Data saved")
    flush(stdout)
    # end
end


# %% Close all worker processes
println("Completed net crop regression for countries in ISO list. Closign worker processes.")
flush(stdout)
rmprocs(workers())

