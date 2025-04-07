########################################
# DISTRIBUTED COMPUTING SETTINGS
########################################
N_WORKERS = 5

########################################
# DIRECTORIES
########################################

# Input Data Directories
HOUSEHOLD_SURVEY_DIR = "/mnt/s3/malaria_survey_extraction/mse_v10/data/custom/20250326/"
RAW_DATASET_DIR = "datasets/"
RAW_SUBNAT_DATASET_DIR = RAW_DATASET_DIR*"subnational/"
POPULATION_RASTER_DIR = "/mnt/s3/mastergrids/Other_Global_Covariates/Population/WorldPop/v3/PopulationCounts_DRC_fixed/5km/"
ACCESSIBILITY_RASTER_DIR = "/mnt/s3/mastergrids/Other_Global_Covariates/Accessibility/Weiss/2015/"

# Output Data Directories
OUTPUT_DIR = "outputs/"
OUTPUT_DATAPREP_DIR = OUTPUT_DIR*"data_prep/"
OUTPUT_SUBNAT_DATAPREP_DIR = OUTPUT_DATAPREP_DIR*"subnational/"
OUTPUT_EXTRACTIONS_DIR = OUTPUT_DIR*"extractions/"
OUTPUT_REGRESSIONS_DIR = OUTPUT_DIR*"regressions/"
OUTPUT_DRAWS_DIR = OUTPUT_DIR*"draws/"
OUTPUT_RASTERS_DIR = OUTPUT_DIR*"INLA/rasters/"
OUTPUT_PLOTS_DIR = "output_plots/"

# BV Output Cube Directory
BV_OUTPUT_DIR = "/mnt/efs/userdata/mvandenberg/data/itn_cube/20241015/"

# Covariates directories and names

## Static
COV_PET_DIR = "/mnt/s3/mastergrids/Other_Global_Covariates/PET/5km/Synoptic/"
COV_PET_FILENAME = "PET_v3.Synoptic.Overall.Data.5km.mean.tif"

COV_ARID_DIR = "/mnt/s3/mastergrids/Other_Global_Covariates/Aridity/5km/Synoptic/"
COV_ARID_FILENAME = "Aridity_Index_v3.Synoptic.Overall.Data.5km.mean.tif"

COV_NTL_DIR = "/mnt/s3/mastergrids/Other_Global_Covariates/NightTimeLights/VIIRS_DNB_Composites/5km/Monthly/"
COV_NTL_FILENAME = "VIIRS-SLC.2014.01.Data.5km.mean.tif"

COV_ELEV_DIR = "/mnt/s3/mastergrids/Other_Global_Covariates/Elevation/SRTM-Elevation/5km/Synoptic/"
COV_ELEV_FILENAME = "SRTM_elevation.Synoptic.Overall.Data.5km.mean.tif"

COV_SLP_DIR = "/mnt/s3/mastergrids/Other_Global_Covariates/Elevation/SRTM-Slope/5km/Synoptic/"
COV_SLP_FILENAME = "SRTM_SlopePCT_Corrected.Synoptic.Overall.Data.5km.mean.tif"

## Monthly
COV_EVI_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD43D6_v061_BRDF_Reflectance/EVI_v061/5km/Monthly/"
COV_LSTD_DIR = "/mnt/s3/mastergrids/MODIS_Global/MOD11A2_v061_LST/LST_Day_v061/5km/Monthly/"
COV_LSTN_DIR = "/mnt/s3/mastergrids/MODIS_Global/MOD11A2_v061_LST/LST_Night_v061/5km/Monthly/"
COV_LSTDELTA_DIR = "/mnt/s3/mastergrids/MODIS_Global/MOD11A2_v061_LST/LST_DiurnalDifference_v061/5km/Monthly/"
COV_TCW_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD43D6_v061_BRDF_Reflectance/TCW_v061/5km/Monthly/"
COV_TSI_DIR = "/mnt/s3/mastergrids/Other_Global_Covariates/TemperatureSuitability/TSI_Pf_Dynamic/5km/Monthly/"
COV_TCB_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD43D6_v061_BRDF_Reflectance/TCB_v061/5km/Monthly/"

## Annual
COV_LAND00_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class00/5km/Annual/"
COV_LAND01_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class01/5km/Annual/"
COV_LAND02_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class02/5km/Annual/"
COV_LAND03_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class03/5km/Annual/"
COV_LAND04_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class04/5km/Annual/"
COV_LAND05_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class05/5km/Annual/"
COV_LAND06_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class06/5km/Annual/"
COV_LAND07_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class07/5km/Annual/"
COV_LAND08_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class08/5km/Annual/"
COV_LAND09_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class09/5km/Annual/"
COV_LAND10_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class10/5km/Annual/"
COV_LAND11_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class11/5km/Annual/"
COV_LAND12_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class12/5km/Annual/"
COV_LAND13_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class13/5km/Annual/"
COV_LAND14_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class14/5km/Annual/"
COV_LAND15_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class15/5km/Annual/"
COV_LAND16_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class16/5km/Annual/"
COV_LAND17_DIR = "/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class17/5km/Annual/"




########################################
# FILENAMES
########################################
# Raw/Input Data

## Accessibility Covariate
ACCESSIBILITY_COV_FILENAME = "accessibility_to_cities_2015_v1.0.tif"

## Shapefiles
ADMIN0_SHAPEFILE = "/mnt/s3/master_geometries/Admin_Units/Global/MAP/2023/MG_5K/admin2023_0_MG_5K.shp"
ADMIN1_SHAPEFILE = "/mnt/s3/master_geometries/Admin_Units/Global/MAP/2023/MG_5K/admin2023_1_MG_5K.shp"

## BC Comparison Data
BV_OUTPUTS_FILENAME = "bv_retention_times.csv"

## National
HOUSEHOLD_SURVEY_FILENAME = "Svy_ALL_ITN_HH_Res.csv"
COUNTRY_CODES_FILENAME = "country_codes.csv"
ISO_LIST_FILENAME = "ISO_list.csv"
DELIVERIES_DATA_FILENAME = "base_manufacturer_deliveries.csv"
DISTRIBUTION_DATA_FILENAME = "net_distributions_cITN_adjusted.csv"
POPULATION_DATA_FILENAME = "ihme_populations.csv"

## Subnational
ADMIN1_AREAID_LEGEND_FILENAME = "admin2023_1_MG_5K_config.csv"
SUBNAT_POPULATION_FILENAME = "map_admin1_zonal_stats_pops.csv"
SUBNAT_DISTRIBUTION_DATA_FILENAME = "net_distributions_admin1_dummy_combined_amp.csv"

# Data Prep Outputs
HOUSEHOLD_SURVEY_DATA_FILENAME = "itn_hh_surveydata_complete_dataeng.csv"
HOUSEHOLD_NAT_SUMMARY_DATA_FILENAME = "npc_monthly_data.csv"
HOUSEHOLD_SUBNAT_SUMMARY_DATA_FILENAME = "subnat_npc_monthly_data.csv"
INLA_DATAPREP_FILENAME = "INLA/inla_dataset.csv"
INLA_REDUCED_DATAPREP_FILENAME = "INLA/inla_dataset_reduced.csv"
INLA_USE_REDUCED_DATAPREP_FILENAME = "INLA/inla_dataset_reduced_use.csv"
COV_LEGEND_FILENAME = "inla_raw_covariates_legend.csv"


########################################
# Model Algorithm Settings
########################################

# Regression Time Period
YEAR_REF_START = 1950 # Earliest possible year for referencing survey extracts (PLEASE DO NOT TOUCH!)
YEAR_NAT_START = 2000 # Start of National SNF Regression
YEAR_NAT_END = 2023 # End of National SNF Regression (inclusive)
YEAR_SUBNAT_TRANS = 2010 # Year to transition into subnat regression form national regression

# Survey Standard Error Inflation
STD_ERR_TAU = 1.5 # Standard error decay time scale adjustment (years)
DEFAULT_INFLATION_SAT_SIZE = 4000 # Default number of survey points needed to assume accurate std err estimate

# Data Extraction Settings
SMOOTHING_WINDOW_WIDTH = 2 # Window size overwhich to smooth survey data

# Missing data handling
MISSING_NETS_SCALE = 1e6

# National MITN Regression
NAT_CROP_MCMC_ITERATIONS = 20000; #50000; # Number of MCMC iterations
NAT_CROP_N_CHAINS = 4; # Number of chains
NAT_CROP_MCMC_BURNIN = 2000 #10000; # Number of iterations to discard for burn-in
NAT_CROP_PROPOSAL_SAMPLING_VAR = 0.004 .*[30,15,20,20,0.25,10] # Compact support sampling
NAT_CROP_SGD_EPOCHS = 5 # Maximum number of SGD alternating iterations
NAT_CROP_SGD_STEPS = 30 # Number of SGD steps in each iteration
NAT_CROP_SGD_ALPHA = 0.03 #0.05 SGD learning rate
NAT_CROP_SGD_EPSILON = 0.03 #0.2 step size over which to calculate numerical gradient
NAT_CROP_SGD_DIC_SAMPLESIZE = 15 # Number of random points overwhich to calculate average DIC in SGD phase

NAT_ACCESS_MCMC_ITERATIONS = 10000; # Number of MCMC iterations
NAT_ACCESS_N_CHAINS = 4; # Number of chains
NAT_ACCESS_MCMC_BURNIN = 2000 #10000; # Number of iterations to discard for burn-in
NAT_ACCESS_PROPOSAL_SAMPLING_VAR = 0.004 .*[30,15,20,20,0.25,10] # Compact support sampling

# National MITN Posterior Draws
NAT_CROP_ACCESS_N_DRAWS = 100
SUBNAT_CROP_ACCESS_N_DRAWS = 100
NAT_SUBNAT_ADJ_MCMC_ITERATIONS = 20000
NAT_SUBNAT_ADJ_MCMC_BURNIN = 5000
NAT_SUBNAT_ADJ_MCMC_SAMPLING_VAR = 0.0004

# Subnational MITN Regresion
## Adjustment Weights
SUBNAT_CROP_ADJ_MCMC_ITERATIONS = 5000; #50000; # Number of MCMC iterations
SUBNAT_CROP_ADJ_N_CHAINS = 4; # Number of chains
SUBNAT_CROP_ADJ_MCMC_BURNIN = 2000 #10000; # Number of iterations to discard for burn-in

## Attrition Parameters
SUBNAT_CROP_ATR_MCMC_ITERATIONS = 15000; #50000; # Number of MCMC iterations
SUBNAT_CROP_ATR_N_CHAINS = 4; # Number of chains
SUBNAT_CROP_ATR_MCMC_BURNIN = 2000 #10000; # Number of iterations to discard for burn-in
SUBNAT_CROP_ATR_SCALING_CONSTANT = 0.1
SUBNAT_CROP_ATR_PROPOSAL_SAMPLING_VAR = SUBNAT_CROP_ATR_SCALING_CONSTANT .*[0.15,8,0.15] # Compact support sampling

# INLA Raster Sample
INLA_USE_N_SAMPLES = 10 # Number of samples to draw for INLA use deviation
INLA_USE_UNCERTAINTY_N_SAMPLES = 10 # Number of samples to draw for calculating INLA use deviation uncertainty