# Load all packages
library(INLA)
library(raster)
library(tidyverse)
library(sf)
library(terra)

# Load custom transformation functions
source("scripts/INLA/transforms.R")

# Import Models
# load("outputs/INLA/model1_npc_complete_logmodel.RData")
load("outputs/INLA/model1_npc_half_logmodel.RData")

# Check summary
summary(m1)

# Random fix
sf_use_s2(FALSE)

#########################
# Load training dataset
#########################

# load INLA regression data
inla_data <- read.csv('outputs/data_prep/INLA/inla_dataset_reduced.csv')
inla_data <- inla_data[seq(2,dim(inla_data)[1],2),]
# inla_data <- inla_data[1:1000,]
inla_data$yearidx <- (inla_data$monthidx %/% 12)+1#*12
inla_data$yearidx

# load Africa shapefile
global_shp <- read_sf("/mnt/s3/master_geometries/Admin_Units/Global/MAP/2023/MG_5K/admin2023_0_MG_5K.shp")

# Filter for required countries
ISO_list <- read.csv("datasets/ISO_list.csv")$ISO
exclusion_ISOs <- c("CPV","ZAF")
filt_ISOs <- setdiff(ISO_list, exclusion_ISOs)

test <- global_shp[global_shp$ISO %in% filt_ISOs,]

africa_geometry <- st_union(test$geometry)

# Load coordinates to regress on
coords <- cbind(inla_data$longitude, inla_data$latitude)

# Back calculate subnational estimates
inla_data$npc_subnat <- inla_data$npc - inla_data$npc_gap
inla_data$access_subnat <- inla_data$access - inla_data$access_gap

#########################
# Setup Covariates and Response Data
#########################
cov_data <- list(yearidx = inla_data$yearidx,
                 static_1 = inla_data$static_1,
                 static_2 = inla_data$static_2,
                 static_3 = inla_data$static_3,
                 annual_1 = inla_data$annual_1,
                 annual_2 = inla_data$annual_2,
                 annual_3 = inla_data$annual_3,
                 annual_4 = inla_data$annual_4,
                 annual_5 = inla_data$annual_5,
                 annual_6 = inla_data$annual_6,
                 annual_7 = inla_data$annual_7,
                 annual_8 = inla_data$annual_8,
                 annual_9 = inla_data$annual_9,
                 annual_10 = inla_data$annual_10,
                 annual_11 = inla_data$annual_11,
                 annual_12 = inla_data$annual_12,
                 annual_13 = inla_data$annual_13,
                 annual_14 = inla_data$annual_14,
                 annual_15 = inla_data$annual_15
)

###########################################
# Start Extracting required coefficients
###########################################
print("Prediction...")

# Get interpolation tensor matrix A
Aprediction <- inla.spde.make.A(mesh = africa_spde, loc = coords,
                                group = inla_data$yearidx,
                                group.mesh = temporal_mesh_annual)

# Calulate spatial structure
sfield_nodes <- m1$summary.random$field['mean']
field <- (Aprediction %*% as.data.frame(sfield_nodes)[, 1])
summary(m1)
# Calculate Predicted values using regression formula
pred <- m1$summary.fixed['Intercept', 'mean'] +
  # m1$summary.fixed['subnat_npc', 'mean'] * cov_data$subnat_npc +
  # m1$summary.fixed['subnat_access', 'mean'] * cov_data$subnat_access +
  m1$summary.fixed['static_1', 'mean'] * cov_data$static_1 +
  m1$summary.fixed['static_2', 'mean'] * cov_data$static_2 +
  m1$summary.fixed['static_3', 'mean'] * cov_data$static_3 +
  m1$summary.fixed['annual_1', 'mean'] * cov_data$annual_1 +
  m1$summary.fixed['annual_2', 'mean'] * cov_data$annual_2 +
  m1$summary.fixed['annual_3', 'mean'] * cov_data$annual_3 +
  m1$summary.fixed['annual_4', 'mean'] * cov_data$annual_4 +
  m1$summary.fixed['annual_5', 'mean'] * cov_data$annual_5 +
  m1$summary.fixed['annual_6', 'mean'] * cov_data$annual_6 +
  m1$summary.fixed['annual_7', 'mean'] * cov_data$annual_7 +
  m1$summary.fixed['annual_8', 'mean'] * cov_data$annual_8 +
  m1$summary.fixed['annual_9', 'mean'] * cov_data$annual_9 +
  m1$summary.fixed['annual_10', 'mean'] * cov_data$annual_10 +
  m1$summary.fixed['annual_11', 'mean'] * cov_data$annual_11 +
  m1$summary.fixed['annual_12', 'mean'] * cov_data$annual_12 +
  m1$summary.fixed['annual_13', 'mean'] * cov_data$annual_13 +
  m1$summary.fixed['annual_14', 'mean'] * cov_data$annual_14 +
  m1$summary.fixed['annual_15', 'mean'] * cov_data$annual_15 +
  # m1$summary.fixed['monthly_1', 'mean'] * cov_data$monthly_1 +
  # m1$summary.fixed['monthly_2', 'mean'] * cov_data$monthly_2 +
  field

# 
epsilon = 0.001
npc_subnat <- inla_data$npc - inla_data$npc_gap
npc_gap_ratio <- exp(pred)
npc_pred <- npc_gap_ratio*npc_subnat + (npc_gap_ratio - 1)*epsilon

# Calculate performance metrics
RMSE <- round(sqrt(mean((inla_data$npc - npc_pred)^2)), digits = 4)
MAE <- round(mean(abs(inla_data$npc - npc_pred)), digits = 4)
CORR <- round(cor(inla_data$npc, as.numeric(npc_pred)), digits = 4)

# Make Fit Plot
title <- str_glue("NPC Fit (Gap Model) Full Training Data\nMAE = {MAE}, RMSE = {RMSE}, CORR = {CORR}")
# title <- str_glue("NPC Fit (Gap Model) 50% Training, Out of Sample \nMAE = {MAE}, RMSE = {RMSE}, CORR = {CORR}")

plot(inla_data$npc, npc_pred, 
     col = rgb(red = 0, green = 0, blue = 1, alpha = 0.1), cex=0.1,
     xlim = c(0,1), ylim = c(0,1),
     xlab = "Survey NPC", ylab = "Fitted NPC", main = title)
abline(a=0, b=1)



