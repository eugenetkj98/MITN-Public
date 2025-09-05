# MAP ITN Model

## Overview
This repository contains the codebase to construct newly developed MITN (Multitype ITN) models using raw national and subnational ITN delivery data, household survey data and various geographical covariates. The model construction pipeline consists of three main phases:
1. Data extraction
2. Stock and flow modelling
3. Spatial disaggregation
4. Raster and map construction

In order to provide better scalability with available computational resources, phases 1,2, and 4 are written in Julia and are multithreaded where possible. Phase 3 is implemented in R to utilise the existing mature INLA package for spatial statistics.

## Repository Structure (For illustrative reference)

NOTE: The following directories are recently deprecated due to implementation with cloud computing APIs and on-prem storage. This repository is used purely for the inspection of code.

### Non-Model Code
#### docs (WORK IN PROGRESS)
Sub-project directory containing Documenter script for producing documentation pages for the code base. Currently WIP, but will eventually contain information on the following:
- How-to-use guide
- Mathematical justification and code details
- Required input data formats and standards
- Outstanding code issues
- References (if relevant)
#### datasets (Note raw datasets are currently stored on-prem for storage capacity reasons. Please contact if looking that have access for running code)
Folder containing all the required raw numerical datasets (anything contained in a .csv file) such as 
- Cleaned household surveys
- Delivery data
- Distribution data
- Country code lookups
- Population numbers
#### outputs
Contains output files from the model construction pipeline. Files may either be intermediate components in the pipeline (e.g. data extractions, posterior draws of net parameters), or finalised model components (e.g. INLA draws, MCMC models for stock and flow mechanics, posterior draws of net crop and access). Files are generally stored in .jld2 format, with the exception of simple table data (.csv), raster files (.tiff) and INLA model outputs (.RData). The .jld2 format is based of HDF5 and should be compatible with other programming languages that support the HDF5 standard. Detail formats of each output file data type is provided in the documentation wiki.
#### output_plots
General storage folder for outputs of various plotting scripts. Contains diagnostic and visualisation outputs of different parts of the ITN model.

### Model Code
#### env_files
Contains the Julia Manifest.toml and Project.toml files listing all required packages in the environment for model construction. Will be instantiated with Pkg.activate() at the start of each pipeline run. Note: The Documenter code utilises its own set of .toml environment files.
#### src
Source files containing bulk of model code. Written as modules containing functions and wrappers to be called from pipeline scripts. (Do not edit unless absolutely necessary.) More details are provided in the documentation wiki.
#### scripts
Contains scripts files that procedurally run components of the pipeline (from src) used to construct the MITN model. General overview of the subfolders
##### DataProcessing
Scripts used in the pre-processing or raw or intermediate datasets into the required format to be input into the model. (e.g. scraping household data, aligning country boundaries, geo-labelling entries, extracting previous BV model estimates for comparison)
##### AuxAnalysis
Non-model critical analyses (e.g. homogeneity tests) and other miscellaneous experiments and exploratory analysis for the ITN model. Not required for outputs. 
##### MITN Folders: NationalMITN and SubnationalMITN
These two folders contain scripts that must be run in order to construct stock and flow model estimates on the national and subnational level with the MITN model. They primarily use functions declared in **src** and output to the **outputs** folder
##### INLA
Contains R code for running the INLA spatial disaggregation step for NPC, Access and Use. Runs via Posit and outputs .tiff component rasters needed to construct final ITN coverage maps. 
##### RasterMaps
Script that combines stock and flow outputs, and INLA rasters to provide final ITN coverage maps. Heavily focused around the Julia Rasters.jl package and mosaic() function. 


## List of Changes from BV-ITN
### Mathematical/Structural Updates
- Added support for multiple simultaneous net types in stock and flow
- Added support for processing and outputting subnational level data (currently prototyping with dummy data)
- Overhauled temporal interpolation approach to use gradient descent rather than random uniform.
- Increased temporal resolution from quarterly to monthly estimates
- Streamlined linear regression model for net life parameters (see documentation for details)
- MCMC adjustment parameters updates to hew subnational stock and flow estimates to national level estimates
- Added support for multiple net attrition models
### Coding/Technical Changes
- Partial multithreading support in stock and flow model
- Usage of standard HDF5 data storage types for better future inter-compatibility with other programming languages
- Porting of stock and flow deterministic model to Julia for better performance and memory management
- Pipeline with standardised input and output formats
### Miscellaneous
- Improved documentation and code organisation

## Documentation wiki
Coming soon once I figure out how to get Github-CI to play nice with Documenter.jl :(
