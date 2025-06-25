# Set working directory
setwd("/mnt/efs/userdata/etan/map-itn")

# Set timeout to allow enough time to install INLA
options(timeout=3600)

# Install required packages in case not in DockerImage by default
install.packages("tidyverse")
install.packages("raster")
install.packages("sf")
install.packages("lattice")
install.packages("grideExtra")
install.packages("tomledit")
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

# Load all packages
library(INLA)
library(raster)
library(tidyverse)
library(sf)
library(terra)
library(matrixStats)
library(tomledit)

# Load custom transformation functions
source("scripts/INLA/transforms.R")

# Import Models
load("/mnt/efs/userdata/etan/mitn_outputs/outputs/INLA/model1_npc_complete_logmodel.RData")

# Check summary
summary(m1)

# Load provided arguments in script
args <- commandArgs(trailingOnly = TRUE)
year <- strtoi(args[1])

# Load TOML Config file
model_config = from_toml(read_toml("/mnt/efs/userdata/etan/map-itn/scripts/awsbatch/configs/model_config.toml"))

##############################
# Select desired year and month to do predictions on
##############################

# Define number of posterior samples
n_samples <- 100
n_samples_saved <- 100

# Get raster image for coordinates to predict on (i.e. only restrict predictions to Africa)
reference_image <- raster('datasets/INLA/MAP_Regions_Unclipped_5k.tif')

# Filter and only extract coordinates from African region
in_country <- which( (!is.na(getValues(reference_image))) & (getValues(reference_image) == 1) )
reference_coordinates <- coordinates(reference_image)[in_country,]

# Get Africa Geometry
sf_use_s2(FALSE)
global_shp <- read_sf("/mnt/s3/master_geometries/Admin_Units/Global/MAP/2023/MG_5K/admin2023_0_MG_5K.shp")

# Filter for required countries
ISO_list <- model_config$ISO_LIST
exclusion_ISOs <- model_config$EXCLUSION_ISOS
filt_ISOs <- setdiff(ISO_list, exclusion_ISOs)

test <- global_shp[global_shp$ISO %in% filt_ISOs,]

africa_geometry <- st_union(test$geometry)

# Convert prediction coordinates to points
pred.points <- SpatialPoints(reference_coordinates, proj4string = crs(africa_geometry))

# Type 1 Covariates: Static/Non Time Varying

## Accessibility needs to be aggregated to 5km by 5km grid
access_raster <- raster('/mnt/s3/mastergrids/Other_Global_Covariates/Accessibility/Weiss/2015/accessibility_to_cities_2015_v1.0.tif')
access_agg_raster <- aggregate(access_raster, 5, fun = sum)
cov_ACCESS <- raster::extract(access_agg_raster, pred.points, df = TRUE)
cov_ACCESS[,2] <- log(cov_ACCESS[,2] + 1e-5)

cov_PET <- raster::extract(raster('/mnt/s3/mastergrids/Other_Global_Covariates/PET/5km/Synoptic/PET_v3.Synoptic.Overall.Data.5km.mean.tif'), pred.points, df = TRUE)
cov_ARID <- raster::extract(raster('/mnt/s3/mastergrids/Other_Global_Covariates/Aridity/5km/Synoptic/Aridity_Index_v3.Synoptic.Overall.Data.5km.mean.tif'), pred.points, df = TRUE)
cov_ARID[,2] <- sqrt(cov_ARID[,2])

cov_NTL <- raster::extract(raster('/mnt/s3/mastergrids/Other_Global_Covariates/NightTimeLights/VIIRS_DNB_Composites/5km/Monthly/VIIRS-SLC.2014.01.Data.5km.mean.tif'), pred.points, df = TRUE)
cov_ELEV <- raster::extract(raster('/mnt/s3/mastergrids/Other_Global_Covariates/Elevation/SRTM-Elevation/5km/Synoptic/SRTM_elevation.Synoptic.Overall.Data.5km.mean.tif'), pred.points, df = TRUE)
cov_SLP <- raster::extract(raster('/mnt/s3/mastergrids/Other_Global_Covariates/Elevation/SRTM-Slope/5km/Synoptic/SRTM_SlopePCT_Corrected.Synoptic.Overall.Data.5km.mean.tif'), pred.points, df = TRUE)

start_year <- model_config$YEAR_NAT_START
# end_year <- 2023

# for (year in start_year:end_year) {


  # print(str_glue("Analysing for y-{year} out of y-{end_year}"))
  
  # Adjust year string as needed
  adj_year <- max(min(year, 2021),2002)
  
  # Get string interpolate for year
  year_str <- str_glue("{adj_year}")
  
  ##############################
  # Extract required covariates (Annual)
  ##############################
  print("Extracting annual covariates...")
  # Type 2 Covariates: Annual varying
  # cov_LAND00 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class00/5km/Annual/IGBP_Landcover_Class-00_Unclassified.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  cov_LAND01 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class01/5km/Annual/IGBP_Landcover_Class-01_Evergreen_Needleleaf_Forest.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  cov_LAND02 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class02/5km/Annual/IGBP_Landcover_Class-02_Evergreen_Broadleaf_Forest.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  # cov_LAND03 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class03/5km/Annual/IGBP_Landcover_Class-03_Deciduous_Needleleaf_Forest.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  cov_LAND04 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class04/5km/Annual/IGBP_Landcover_Class-04_Deciduous_Broadleaf_Forest.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  cov_LAND05 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class05/5km/Annual/IGBP_Landcover_Class-05_Mixed_Forest.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  cov_LAND06 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class06/5km/Annual/IGBP_Landcover_Class-06_Closed_Shrublands.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  cov_LAND07 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class07/5km/Annual/IGBP_Landcover_Class-07_Open_Shrublands.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  cov_LAND08 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class08/5km/Annual/IGBP_Landcover_Class-08_Woody_Savannas.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  cov_LAND09 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class09/5km/Annual/IGBP_Landcover_Class-09_Savannas.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  cov_LAND10 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class10/5km/Annual/IGBP_Landcover_Class-10_Grasslands.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  cov_LAND11 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class11/5km/Annual/IGBP_Landcover_Class-11_Permanent_Wetlands.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  cov_LAND12 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class12/5km/Annual/IGBP_Landcover_Class-12_Croplands.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  cov_LAND13 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class13/5km/Annual/IGBP_Landcover_Class-13_Urban_And_Built_Up.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  cov_LAND14 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class14/5km/Annual/IGBP_Landcover_Class-14_Cropland_Natural_Vegetation_Mosaic.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  # cov_LAND15 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class15/5km/Annual/IGBP_Landcover_Class-15_Snow_And_Ice.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  cov_LAND16 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class16/5km/Annual/IGBP_Landcover_Class-16_Barren_Or_Sparsely_Populated.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  cov_LAND17 <- raster::extract(raster(str_glue('/mnt/s3/mastergrids/MODIS_Global/MCD12Q1_v061_Annual_Landcover/IGBP_Landcover_Class17/5km/Annual/IGBP_Landcover_Class-17_Water.{year_str}.Annual.Data.5km.fraction.tif')), pred.points, df = TRUE)
  
  ##############################
  # Preprocess Covariate Values
  ##############################
  # Normalise all covariates and write as new variables
  cov_norm_constants <- read.csv("/mnt/efs/userdata/etan/mitn_outputs/outputs/data_prep/INLA/covariate_normalisation_constants.csv")
  norm_cov_var_names <- c()
  
  for (i in 1:(length(cov_norm_constants$cov)-6)){
    cov_name <- cov_norm_constants$cov[i]
    norm_var_string <- paste("norm_", cov_name, sep = "")
    mu <- cov_norm_constants$mean[i]
    sigma <- cov_norm_constants$std[i]
    # assign(norm_var_string, (get(cov_name)[,2]  - mu) / sigma)
    assign(norm_var_string, (get(cov_name)[,2]) / sigma)
    norm_cov_var_names[i] <- norm_var_string
    print(norm_var_string)
  }
  
  # Combine all normalised covariates into a single matrix
  norm_cov_matrix <- as.matrix(do.call("rbind", as.list(lapply(norm_cov_var_names, get))))
  
  # Import projection matrix to covariates and calculated values for reduced proj_matrix
  M_proj_raw <- read.csv("/mnt/efs/userdata/etan/mitn_outputs/outputs/data_prep/INLA/proj_matrix.csv")
  M_proj <- as.matrix(M_proj_raw[1:(dim(M_proj_raw)[1]-2),2:(dim(M_proj_raw)[2]-6)])
  
  proj_cov_datavalues <- M_proj %*% norm_cov_matrix
  
  # Write projected data into a list dataset
  proj_cov_names <- M_proj_raw$raw_cov[1:dim(M_proj)[1]]
  
  proj_cov_dataset <- as.data.frame(t(proj_cov_datavalues))
  colnames(proj_cov_dataset) <- proj_cov_names
  
  ###########################################
  # Start Extracting required coefficients
  ###########################################
  
  print("Prediction...")
  
  yearidx = year-start_year+1
  
  # Define current groupid slice
  groupid <- replace(rep(0,dim(reference_coordinates)[1]), c(1:dim(reference_coordinates)[1]), yearidx)
  
  # Get interpolation tensor matrix A
  Aprediction <- inla.spde.make.A(mesh = africa_spde, loc = reference_coordinates,
                                  group = groupid,
                                  group.mesh = temporal_mesh_annual)
  
  ##########################################
  ####### Direct formula method to get mean raster
  ##########################################
  # Calculate spatial structure
  sfield_nodes <- m1$summary.random$field['mean']
  field <- (Aprediction %*% as.data.frame(sfield_nodes)[, 1])
  summary(m1)
  # Calculate Predicted values using regression formula
  pred_mean <- #m1$summary.fixed['Intercept', 'mean'] +
                m1$summary.fixed['static_1', 'mean'] * proj_cov_dataset$static_1 +
                m1$summary.fixed['static_2', 'mean'] * proj_cov_dataset$static_2 +
                m1$summary.fixed['static_3', 'mean'] * proj_cov_dataset$static_3 +
                m1$summary.fixed['annual_1', 'mean'] * proj_cov_dataset$annual_1 +
                m1$summary.fixed['annual_2', 'mean'] * proj_cov_dataset$annual_2 +
                m1$summary.fixed['annual_3', 'mean'] * proj_cov_dataset$annual_3 +
                m1$summary.fixed['annual_4', 'mean'] * proj_cov_dataset$annual_4 +
                m1$summary.fixed['annual_5', 'mean'] * proj_cov_dataset$annual_5 +
                m1$summary.fixed['annual_6', 'mean'] * proj_cov_dataset$annual_6 +
                m1$summary.fixed['annual_7', 'mean'] * proj_cov_dataset$annual_7 +
                m1$summary.fixed['annual_8', 'mean'] * proj_cov_dataset$annual_8 +
                m1$summary.fixed['annual_9', 'mean'] * proj_cov_dataset$annual_9 +
                m1$summary.fixed['annual_10', 'mean'] * proj_cov_dataset$annual_10 +
                m1$summary.fixed['annual_11', 'mean'] * proj_cov_dataset$annual_11 +
                m1$summary.fixed['annual_12', 'mean'] * proj_cov_dataset$annual_12 +
                m1$summary.fixed['annual_13', 'mean'] * proj_cov_dataset$annual_13 +
                m1$summary.fixed['annual_14', 'mean'] * proj_cov_dataset$annual_14 +
                m1$summary.fixed['annual_15', 'mean'] * proj_cov_dataset$annual_15 +
                field
  
  ###########################################
  # ####### Joint posterior draws
  ###########################################
  # Formula to evaluate when drawing posteriors
  inla_model_eval_fun <- function(){
    return(
      #Intercept +
        static_1 * proj_cov_dataset$static_1 +
        static_2 * proj_cov_dataset$static_2 +
        static_3 * proj_cov_dataset$static_3 +
        # static_4 * proj_cov_dataset$static_4 +
        annual_1 * proj_cov_dataset$annual_1 +
        annual_2 * proj_cov_dataset$annual_2 +
        annual_3 * proj_cov_dataset$annual_3 +
        annual_4 * proj_cov_dataset$annual_4 +
        annual_5 * proj_cov_dataset$annual_5 +
        annual_6 * proj_cov_dataset$annual_6 +
        annual_7 * proj_cov_dataset$annual_7 +
        annual_8 * proj_cov_dataset$annual_8 +
        annual_9 * proj_cov_dataset$annual_9 +
        annual_10 * proj_cov_dataset$annual_10 +
        annual_11 * proj_cov_dataset$annual_11 +
        annual_12 * proj_cov_dataset$annual_12 +
        annual_13 * proj_cov_dataset$annual_13 +
        annual_14 * proj_cov_dataset$annual_14 +
        annual_15 * proj_cov_dataset$annual_15 +
        Aprediction %*% field
    )
  }
  
  # Generate samples
  samples <- inla.posterior.sample(n = n_samples, m1, verbose = TRUE)
  eval_samples <- inla.posterior.sample.eval(inla_model_eval_fun, samples)
  
  # Calculate the posterior raster values and inverse transform as appropriate
  pred_samples <- matrix(0, nrow = dim(Aprediction)[1], ncol = n_samples)
  for (i in 1:n_samples){
    print(str_glue("Extracting sample {i} of {n_samples}..."))
    pred_samples[,i] <- eval_samples[i][[1]]@x
  }
  
  z_samples <- pred_samples 
  
  # Calculate the average and standard deviation raster
  z_mean <- pred_mean

  ###########################################
  # Save output raster
  ###########################################
  
  print("Saving...")
  # write results as a raster
  x <- as.matrix(reference_coordinates)
  
  pr.mdg.out_mean <- rasterFromXYZ(cbind(x, z_mean), crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs")

  # Save Raster
  save_filename_mean = str_glue("/mnt/efs/userdata/etan/mitn_outputs/outputs/INLA/rasters/inla_logmodel_npc/NPC_logmodel_{year}_mean.tif")
  writeRaster(pr.mdg.out_mean, save_filename_mean, NAflag = -9999, overwrite = TRUE)
  
  # Save sample draws for calculating quantiles later
  for (i in 1:n_samples_saved){
    pr.mdg.out_sample <- rasterFromXYZ(cbind(x, z_samples[,i]), crs = "+proj=longlat +datum=WGS84 +no_defs +type=crs")
    save_filename_sample = str_glue("/mnt/efs/userdata/etan/mitn_outputs/outputs/INLA/rasters/inla_logmodel_npc/NPC_logmodel_{year}_sample_{i}.tif")
    writeRaster(pr.mdg.out_sample, save_filename_sample, NAflag = -9999, overwrite = TRUE)
  }


# }

print("HURRAH! I FINISHED THE R SCRIPT THANK GOD")
