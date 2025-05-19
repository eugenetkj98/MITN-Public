"""
Author: Eugene Tan
Date Created: 15/5/2025
Last Updated: 15/5/2025
Script to parse TOML config file into nice Julia formatting that is consistent with existing codebase structure
"""
# %% Load packages
using TOML

# %% Load TOML file
dir_config = TOML.parsefile("scripts/awsbatch/configs/dir_config.toml")
model_config = TOML.parsefile("scripts/awsbatch/configs/model_config.toml")

########################################
# DIRECTORIES
########################################
# Input Data Directories
BASE_DIR = dir_config["base"]["BASE_DIR"];
HOUSEHOLD_SURVEY_DIR = dir_config["input_raw_dirs"]["HOUSEHOLD_SURVEY_DIR"];
RAW_DATASET_DIR = BASE_DIR*dir_config["input_rel_dirs"]["BASEDIR_REL_RAW_DATASET_DIR"];
RAW_SUBNAT_DATASET_DIR = RAW_DATASET_DIR*dir_config["input_rel_dirs"]["BASEDIR_REL_RAW_SUBNAT_DATASET_DIR"];
POPULATION_RASTER_DIR = dir_config["input_pop_raw_dir"]["POPULATION_RASTER_DIR"];

# Output Data Directories
OUTPUT_DIR = BASE_DIR*dir_config["output_rel_dirs"]["BASEDIR_REL_OUTPUT_DIR"];
OUTPUT_DATAPREP_DIR = OUTPUT_DIR*dir_config["output_rel_dirs"]["OUTPUTDIR_REL_DATAPREP_DIR"];
OUTPUT_SUBNAT_DATAPREP_DIR = OUTPUT_DATAPREP_DIR*dir_config["output_rel_dirs"]["OUTPUTDIR_REL_SUBNAT_DATAPREP_DIR"];
OUTPUT_EXTRACTIONS_DIR = OUTPUT_DIR*dir_config["output_rel_dirs"]["OUTPUTDIR_REL_EXTRACTIONS_DIR"];
OUTPUT_REGRESSIONS_DIR = OUTPUT_DIR*dir_config["output_rel_dirs"]["OUTPUTDIR_REL_REGRESSIONS_DIR"];
OUTPUT_DRAWS_DIR = OUTPUT_DIR*dir_config["output_rel_dirs"]["OUTPUTDIR_REL_DRAWS_DIR"];
OUTPUT_RASTERS_DIR = OUTPUT_DIR*dir_config["output_rel_dirs"]["OUTPUTDIR_REL_RASTERS_DIR"];
OUTPUT_PLOTS_DIR = BASE_DIR*dir_config["output_rel_dirs"]["BASEDIR_REL_OUTPUT_PLOTS_DIR"];
OUTPUT_FORWARD_PRED_DIR = OUTPUT_DIR*dir_config["output_rel_dirs"]["OUTPUTDIR_REL_FORWARD_PRED_DIR"];

# BV Output Cube Directory
BV_OUTPUT_DIR = dir_config["input_raw_dirs"]["BV_OUTPUT_DIR"];

# Covariates directories and names
## Static
ACCESSIBILITY_RASTER_DIR = dir_config["input_static_cov_raw_dirs"]["ACCESSIBILITY_RASTER_DIR"];
ACCESSIBILITY_COV_FILENAME = dir_config["input_static_cov_raw_dirs"]["ACCESSIBILITY_COV_FILENAME"];

COV_PET_DIR = dir_config["input_static_cov_raw_dirs"]["COV_PET_DIR"];
COV_PET_FILENAME = dir_config["input_static_cov_raw_dirs"]["COV_PET_FILENAME"];

COV_ARID_DIR = dir_config["input_static_cov_raw_dirs"]["COV_ARID_DIR"];
COV_ARID_FILENAME = dir_config["input_static_cov_raw_dirs"]["COV_ARID_FILENAME"];

COV_NTL_DIR = dir_config["input_static_cov_raw_dirs"]["COV_NTL_DIR"];
COV_NTL_FILENAME = dir_config["input_static_cov_raw_dirs"]["COV_NTL_FILENAME"];

COV_ELEV_DIR = dir_config["input_static_cov_raw_dirs"]["COV_ELEV_DIR"];
COV_ELEV_FILENAME = dir_config["input_static_cov_raw_dirs"]["COV_ELEV_FILENAME"];

COV_SLP_DIR = dir_config["input_static_cov_raw_dirs"]["COV_SLP_DIR"];
COV_SLP_FILENAME = dir_config["input_static_cov_raw_dirs"]["COV_SLP_FILENAME"];

## Monthly
COV_EVI_DIR = dir_config["input_monthly_cov_raw_dirs"]["COV_EVI_DIR"];
COV_LSTD_DIR = dir_config["input_monthly_cov_raw_dirs"]["COV_LSTD_DIR"];
COV_LSTN_DIR = dir_config["input_monthly_cov_raw_dirs"]["COV_LSTN_DIR"];
COV_LSTDELTA_DIR = dir_config["input_monthly_cov_raw_dirs"]["COV_LSTDELTA_DIR"];
COV_TCW_DIR = dir_config["input_monthly_cov_raw_dirs"]["COV_TCW_DIR"];
COV_TSI_DIR = dir_config["input_monthly_cov_raw_dirs"]["COV_TSI_DIR"];
COV_TCB_DIR = dir_config["input_monthly_cov_raw_dirs"]["COV_TCB_DIR"];

## Annual
COV_LAND00_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND00_DIR"];
COV_LAND01_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND01_DIR"];
COV_LAND02_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND02_DIR"];
COV_LAND03_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND03_DIR"];
COV_LAND04_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND04_DIR"];
COV_LAND05_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND05_DIR"];
COV_LAND06_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND06_DIR"];
COV_LAND07_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND07_DIR"];
COV_LAND08_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND08_DIR"];
COV_LAND09_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND09_DIR"];
COV_LAND10_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND10_DIR"];
COV_LAND11_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND11_DIR"];
COV_LAND12_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND12_DIR"];
COV_LAND13_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND13_DIR"];
COV_LAND14_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND14_DIR"];
COV_LAND15_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND15_DIR"];
COV_LAND16_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND16_DIR"];
COV_LAND17_DIR = dir_config["input_annual_cov_raw_dirs"]["COV_LAND17_DIR"];

########################################
# FILENAMES
########################################
# Raw/Input Data

## Shapefiles
ADMIN0_SHAPEFILE = dir_config["input_admin_shapefiles"]["ADMIN0_SHAPEFILE"];
ADMIN1_SHAPEFILE = dir_config["input_admin_shapefiles"]["ADMIN1_SHAPEFILE"];

## BC Comparison Data
BV_OUTPUTS_FILENAME = dir_config["filenames"]["BV_OUTPUTS_FILENAME"];

## National
HOUSEHOLD_SURVEY_FILENAME = dir_config["filenames"]["HOUSEHOLD_SURVEY_FILENAME"];
COUNTRY_CODES_FILENAME = dir_config["filenames"]["COUNTRY_CODES_FILENAME"];
ISO_LIST_FILENAME = dir_config["filenames"]["ISO_LIST_FILENAME"];
DELIVERIES_DATA_FILENAME = dir_config["filenames"]["DELIVERIES_DATA_FILENAME"];
DISTRIBUTION_DATA_FILENAME = dir_config["filenames"]["DISTRIBUTION_DATA_FILENAME"];
MULTITYPE_DISTRIBUTION_DATA_FILENAME = dir_config["filenames"]["MULTITYPE_DISTRIBUTION_DATA_FILENAME"];
POPULATION_DATA_FILENAME = dir_config["filenames"]["POPULATION_DATA_FILENAME"];

## Subnational
ADMIN1_AREAID_LEGEND_FILENAME = dir_config["filenames"]["ADMIN1_AREAID_LEGEND_FILENAME"];
SUBNAT_POPULATION_FILENAME = dir_config["filenames"]["SUBNAT_POPULATION_FILENAME"];
SUBNAT_DISTRIBUTION_DATA_FILENAME = dir_config["filenames"]["SUBNAT_DISTRIBUTION_DATA_FILENAME"];

# Data Prep Outputs
HOUSEHOLD_SURVEY_DATA_FILENAME = dir_config["filenames"]["HOUSEHOLD_SURVEY_DATA_FILENAME"];
HOUSEHOLD_SUBNAT_SURVEY_DATA_FILENAME = dir_config["filenames"]["HOUSEHOLD_SUBNAT_SURVEY_DATA_FILENAME"];
HOUSEHOLD_NAT_SUMMARY_DATA_FILENAME = dir_config["filenames"]["HOUSEHOLD_NAT_SUMMARY_DATA_FILENAME"];
HOUSEHOLD_SUBNAT_SUMMARY_DATA_FILENAME = dir_config["filenames"]["HOUSEHOLD_SUBNAT_SUMMARY_DATA_FILENAME"];
INLA_DATAPREP_FILENAME = dir_config["filenames"]["INLA_DATAPREP_FILENAME"];
INLA_REDUCED_DATAPREP_FILENAME = dir_config["filenames"]["INLA_REDUCED_DATAPREP_FILENAME"];
INLA_USE_REDUCED_DATAPREP_FILENAME = dir_config["filenames"]["INLA_USE_REDUCED_DATAPREP_FILENAME"];
COV_LEGEND_FILENAME = dir_config["filenames"]["COV_LEGEND_FILENAME"];


########################################
# Model Algorithm Settings
########################################
# List of ISOs
ISO_LIST = model_config["ISO_LIST"]

# Special rules for country
EXCLUSION_ISOS = model_config["EXCLUSION_ISOS"];
EXCLUSION_CITN_ISOS = model_config["EXCLUSION_CITN_ISOS"]; # Countries to lump cITN and LLINs together when doing calculation

# Regression Time Period
YEAR_REF_START = model_config["YEAR_REF_START"]; # Earliest possible year for referencing survey extracts (PLEASE DO NOT TOUCH!)
YEAR_NAT_START = model_config["YEAR_NAT_START"]; # Start of National SNF Regression
YEAR_NAT_END = model_config["YEAR_NAT_END"]; # End of National SNF Regression (inclusive)
YEAR_SUBNAT_TRANS = model_config["YEAR_SUBNAT_TRANS"]; # Year to transition into subnat regression form national regression

# Survey Standard Error Inflation
STD_ERR_TAU = model_config["STD_ERR_TAU"]; # Standard error decay time scale adjustment (years)
DEFAULT_INFLATION_SAT_SIZE = model_config["DEFAULT_INFLATION_SAT_SIZE"]; # Default number of survey points needed to assume accurate std err estimate

# Data Extraction Settings
SMOOTHING_WINDOW_WIDTH = model_config["SMOOTHING_WINDOW_WIDTH"]; # Window size overwhich to smooth survey data

# Missing data handling
MISSING_NETS_SCALE = model_config["MISSING_NETS_SCALE"];

# National MITN Regression
NAT_CROP_MCMC_ITERATIONS = model_config["NAT_CROP_MCMC_ITERATIONS"]; #50000; # Number of MCMC iterations
NAT_CROP_N_CHAINS = model_config["NAT_CROP_N_CHAINS"]; # Number of chains
NAT_CROP_MCMC_BURNIN = model_config["NAT_CROP_MCMC_BURNIN"]; #10000; # Number of iterations to discard for burn-in
NAT_CROP_PROPOSAL_SAMPLING_VAR = model_config["NAT_CROP_PROPOSAL_SAMPLING_VAR"]; # Compact support sampling
NAT_CROP_SGD_EPOCHS = model_config["NAT_CROP_SGD_EPOCHS"]; # Maximum number of SGD alternating iterations
NAT_CROP_SGD_STEPS = model_config["NAT_CROP_SGD_STEPS"]; # Number of SGD steps in each iteration
NAT_CROP_SGD_ALPHA = model_config["NAT_CROP_SGD_ALPHA"]; #0.05 SGD learning rate
NAT_CROP_SGD_EPSILON = model_config["NAT_CROP_SGD_EPSILON"]; #0.2 step size over which to calculate numerical gradient
NAT_CROP_SGD_DIC_SAMPLESIZE = model_config["NAT_CROP_SGD_DIC_SAMPLESIZE"]; # Number of random points overwhich to calculate average DIC in SGD phase

NAT_ACCESS_MCMC_ITERATIONS = model_config["NAT_ACCESS_MCMC_ITERATIONS"]; # Number of MCMC iterations
NAT_ACCESS_N_CHAINS = model_config["NAT_ACCESS_N_CHAINS"]; # Number of chains
NAT_ACCESS_MCMC_BURNIN = model_config["NAT_ACCESS_MCMC_BURNIN"]; #10000; # Number of iterations to discard for burn-in
NAT_ACCESS_PROPOSAL_SAMPLING_VAR = model_config["NAT_ACCESS_PROPOSAL_SAMPLING_VAR"]; # Compact support sampling

# National MITN Posterior Draws
NAT_CROP_ACCESS_N_DRAWS = model_config["NAT_CROP_ACCESS_N_DRAWS"];
NAT_CROPAGE_N_DRAWS = model_config["NAT_CROPAGE_N_DRAWS"];
SUBNAT_CROP_ACCESS_N_DRAWS = model_config["SUBNAT_CROP_ACCESS_N_DRAWS"];
NAT_SUBNAT_ADJ_MCMC_ITERATIONS = model_config["NAT_SUBNAT_ADJ_MCMC_ITERATIONS"];
NAT_SUBNAT_ADJ_MCMC_BURNIN = model_config["NAT_SUBNAT_ADJ_MCMC_BURNIN"];
NAT_SUBNAT_ADJ_MCMC_SAMPLING_VAR = model_config["NAT_SUBNAT_ADJ_MCMC_SAMPLING_VAR"];

# Subnational MITN Regresion
## Adjustment Weights
SUBNAT_CROP_ADJ_MCMC_ITERATIONS = model_config["SUBNAT_CROP_ADJ_MCMC_ITERATIONS"]; #50000; # Number of MCMC iterations
SUBNAT_CROP_ADJ_N_CHAINS = model_config["SUBNAT_CROP_ADJ_N_CHAINS"]; # Number of chains
SUBNAT_CROP_ADJ_MCMC_BURNIN = model_config["SUBNAT_CROP_ADJ_MCMC_BURNIN"];# Number of iterations to discard for burn-in

## Attrition Parameters
SUBNAT_CROP_ATR_MCMC_ITERATIONS = model_config["SUBNAT_CROP_ATR_MCMC_ITERATIONS"]; # Number of MCMC iterations
SUBNAT_CROP_ATR_N_CHAINS = model_config["SUBNAT_CROP_ATR_N_CHAINS"]; # Number of chains
SUBNAT_CROP_ATR_MCMC_BURNIN = model_config["SUBNAT_CROP_ATR_MCMC_BURNIN"]; # Number of iterations to discard for burn-in
SUBNAT_CROP_ATR_SCALING_CONSTANT = model_config["SUBNAT_CROP_ATR_SCALING_CONSTANT"];
SUBNAT_CROP_ATR_PROPOSAL_SAMPLING_VAR = model_config["SUBNAT_CROP_ATR_PROPOSAL_SAMPLING_VAR"]; # Compact support sampling

# INLA Raster Sample
INLA_N_SAMPLES = model_config["INLA_N_SAMPLES"];
INLA_USE_N_SAMPLES = model_config["INLA_USE_N_SAMPLES"]; # Number of samples to draw for INLA use deviation
INLA_UNCERTAINTY_N_SAMPLES = model_config["INLA_UNCERTAINTY_N_SAMPLES"]; # Number of samples to draw for calculating INLA use deviation uncertainty