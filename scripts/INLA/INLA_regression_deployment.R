# Set working directory
setwd("/mnt/efs/userdata/etan/map-itn")

# Set timeout to allow enough time to install INLA
options(timeout=600)

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

# Load TOML Config file
model_config = from_toml(read_toml("/mnt/efs/userdata/etan/map-itn/scripts/awsbatch/configs/model_config.toml"))

# Load custom transformation functions
source("/mnt/efs/userdata/etan/map-itn/scripts/INLA/transforms.R")

# Random fix
sf_use_s2(FALSE)

# load INLA regression data
inla_data <- read.csv('/mnt/efs/userdata/etan/map-itn/outputs/data_prep/INLA/inla_dataset_reduced.csv')
# inla_data <- inla_data[seq(1,dim(inla_data)[1],2),]
# inla_data <- inla_data[which(inla_data$access > 0),]
inla_data <- inla_data[which(inla_data$npc > 0),]
inla_data$yearidx <- (inla_data$monthidx %/% 12)+1#*12
inla_data$yearidx

# load Africa shapefile
global_shp <- read_sf("/mnt/s3/master_geometries/Admin_Units/Global/MAP/2023/MG_5K/admin2023_0_MG_5K.shp")

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
n_years = (end_year-start_year + 1)
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
# Calculate optimum ihs theta and transformation for access response
inla_data$access_subnat <- inla_data$access - inla_data$access_gap
inla_data$npc_subnat <- inla_data$npc - inla_data$npc_gap
for (i in 1:length(inla_data$access)){
  inla_data$dep_subnat <- min(inla_data$access_subnat[i]/(2*inla_data$npc_subnat[i]),1)
  inla_data$dep[i] <- min(inla_data$access[i]/(2*inla_data$npc[i]),1)
  inla_data$p_dep[i] <- p_transform(inla_data$dep[i], inla_data$dep_subnat[i], n = 2)
}

dep_theta <- optimise(ihs_loglik, lower = 0.001, upper = 1000, x = gap_emplogit(inla_data$p_dep[which(inla_data$p_dep != 0)]), maximum = TRUE)$maximum
res_dep_gap <- ihs(gap_emplogit(inla_data$p_dep), dep_theta)

response_data <- list(res_dep_gap = res_dep_gap)

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
                 annual_15 = inla_data$annual_15,
                 monthly_1 = inla_data$monthly_1,
                 monthly_2 = inla_data$monthly_2
)

effects_data_annual <- list(c(S_index_annual, list(Intercept = 1)),cov_data)

africa_stack_annual <- inla.stack(data = response_data,
                                  A = list(A_proj_annual,1),
                                  effects = effects_data_annual,
                                  tag = "npc.data")

#############################
# Fit ACCESS GAP model
#############################
print("Fitting Deployment gap spatio-temporal model...")
m1 <- inla(res_dep_gap ~  -1 + #Intercept + 
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
             annual_11+
             annual_12 +
             annual_13 +
             annual_14 +
             annual_15 +
             f(field, model = spde, group = field.group, 
               control.group = list(model = 'ar1') ),
           data = inla.stack.data(africa_stack_annual, spde = africa_spde),
           family = "gaussian",
           control.predictor = list(A = inla.stack.A(africa_stack_annual), compute = TRUE),
           control.compute = list(cpo = TRUE, dic = TRUE, config = TRUE), 
           control.inla = list(strategy = "adaptive", int.strategy = "eb"),
           verbose = TRUE)

print("Saving Deployment gap model outputs...")

save(africa_mesh, africa_spde, temporal_mesh_annual, m1, dep_theta, file = "/mnt/efs/userdata/etan/map-itn/outputs/INLA/model1_deployment_complete_pmodel.RData")

print("Saved Deployment model")

print("HURRAH! I FINISHED THE R SCRIPT THANK GOD")