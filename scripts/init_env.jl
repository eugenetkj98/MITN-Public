"""
Author: Eugene Tan
Date Created: 4/11/2024
Last Updated: 4/11/2024
Script to initialise environment and load all required sub_directories to check for modules
"""

# Function to recursively find all sub_directories and load all child directories
function recurse_loaddir(src_dir)
    push!(LOAD_PATH, src_dir)
    sub_dirs = src_dir.*readdir(src_dir)
    isdir_bool = isdir.(sub_dirs)
    isdir_idx = findall(isdir_bool)
    for idx in isdir_idx
        push!(LOAD_PATH,sub_dirs[idx]*"/")
        sum(isdir.(sub_dirs[idx]*"/".*readdir(sub_dirs[idx]*"/"))) > 0
        if sum(isdir.(sub_dirs[idx]*"/".*readdir(sub_dirs[idx]*"/"))) > 0
            subsub_dirs = sub_dirs[idx]*"/".*readdir(sub_dirs[idx]*"/")

            for subsub_dir in subsub_dirs
                recurse_loaddir(subsub_dir)
            end
        end
    end
end

# Load subdirectories from src folder
src_dir = pwd()*"/src/"
recurse_loaddir(src_dir)

# Activate required environment
# Package manager (This should install the required packages. If not, run the "instantiate" command via Pkg)
using Pkg
Pkg.activate("env_files")

# %% Temp install packages
Pkg.add(["JLD2",
            "CSV",
            "DataFrames",
            "Missings",
            "LinearAlgebra",
            "StatsBase",
            "StatsPlots",
            "Random",
            "Unicode",
            "Dates",
            "ProgressBars"])

Pkg.add(["GeoIO",
            "GeoInterface",
            "Shapefile",
            "Rasters",
            "MultivariateStats",
            "SparseArrays",
            "Documenter"])