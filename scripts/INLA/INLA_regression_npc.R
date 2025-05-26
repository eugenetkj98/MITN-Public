# Set working directory
setwd("/mnt/efs/userdata/etan/map-itn")

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
library(lattice)     
library(gridExtra)
library(tomledit)

# Load custom transformation functions
source("scripts/INLA/transforms.R")

# Load TOML Config file
model_config = from_toml(read_toml("/mnt/efs/userdata/etan/map-itn/scripts/awsbatch/configs/model_config.toml"))

# Random fix
sf_use_s2(FALSE)

# load INLA regression data
inla_data <- read.csv('/mnt/efs/userdata/etan/map-itn/outputs/data_prep/INLA/inla_dataset_reduced.csv')
# inla_data <- inla_data[seq(1,dim(inla_data)[1],2),]
# inla_data <- inla_data[sample(dim(inla_data)[1], 25000, replace = FALSE),]
inla_data$yearidx <- (inla_data$monthidx %/% 12)+1#*12
inla_data$yearidx

# load Africa shapefile
global_shp <- read_sf("/mnt/s3/master_geometries/Admin_Units/Global/MAP/2023/MG_5K/admin2023_0_MG_5K.shp")

print("I got here :D")

# Filter for required countries
ISO_list <- model_config$ISO_LIST
exclusion_ISOs <- model_config$EXCLUSION_ISOS
filt_ISOs <- setdiff(ISO_list, exclusion_ISOs)

test <- global_shp[global_shp$ISO %in% filt_ISOs,]

africa_geometry <- st_union(test$geometry)

coords <- cbind(inla_data$longitude, inla_data$latitude)
africa_mesh <- inla.mesh.2d(loc = coords,
                            boundary = africa_geometry,
                            max.edge = c(1.25,3),
                            offset = c(1,5),
                            cutoff = 0.6)

africa_spde <- inla.spde2.matern(mesh = africa_mesh)

plot(africa_mesh)

# Construct temporal parts of model
start_year = model_config$YEAR_NAT_START
end_year = model_config$YEAR_NAT_END
n_years = (end_year-start_year+1)

# generate temporal mesh
temporal_mesh_annual <- inla.mesh.1d(seq(1,n_years+1,by=2),interval=c(1, n_years+1),degree=2)

# Make projection matrices
A_proj_annual <- inla.spde.make.A(mesh = africa_spde, loc = coords,
                                  group = inla_data$yearidx,
                                  group.mesh = temporal_mesh_annual)
S_index_annual <- inla.spde.make.index(name = "field",
                                       n.spde = africa_mesh$n,
                                       n.group = temporal_mesh_annual$m)

# Construct INLA stack for spatial
# Calculate NPC metric
epsilon <- 0.001
inla_data$npc_subnat <- inla_data$npc - inla_data$npc_gap
npc_gap_ratio <- (inla_data$npc + epsilon)/(inla_data$npc_subnat + epsilon)
res_npc_gap <- log(npc_gap_ratio)

# Setup response data
response_data <- list(res_npc_gap = res_npc_gap)

cov_data <- list(yearidx = inla_data$yearidx,
                 monthidx = inla_data$monthidx,
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

effects_data_annual <- list(c(S_index_annual, list(Intercept = 1)),cov_data)

africa_stack_annual <- inla.stack(data = response_data,
                                  A = list(A_proj_annual,1),
                                  effects = effects_data_annual,
                                  tag = "npc.data")

#############################
# Fit NPC GAP model
#############################
print("Fitting NPC gap spatio-temporal model...")
m1 <- inla(res_npc_gap ~ -1 + 
             static_1 +
             static_2 +
             static_3 +
             annual_1 +
             annual_2 +
             annual_3 +
             annual_4 +
             annual_5 +
             annual_6 +
             annual_7 +
             annual_8 +
             annual_9 +
             annual_10 +
             annual_11 +
             annual_12 +
             annual_13 +
             annual_14 +
             annual_15 +
             # annual_16 +
             # annual_17 +
             # annual_18 +
             f(field, model = spde, group = field.group, 
               control.group = list(model = 'ar1') ),
           data = inla.stack.data(africa_stack_annual, spde = africa_spde),
           family = "gaussian",
           control.predictor = list(A = inla.stack.A(africa_stack_annual), compute = TRUE),
           control.compute = list(cpo = TRUE, dic = TRUE, config = TRUE), 
           control.inla = list(strategy = "adaptive", int.strategy = "eb"),
           verbose = TRUE)

print("Saving NPC gap model outputs...")

save(africa_mesh, africa_spde, temporal_mesh_annual, m1, epsilon, file = "/mnt/efs/userdata/etan/map-itn/outputs/INLA/model1_npc_complete_logmodel.RData")

print("Saved NPC gap model")

summary(m1)

print("HURRAH! I FINISHED THE R SCRIPT THANK GOD")


