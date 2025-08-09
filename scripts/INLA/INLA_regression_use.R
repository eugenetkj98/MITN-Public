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
install.packages("remotes")
library(remotes)
remotes::install_version("INLA", version = "24.12.11",
repos = c(getOption("repos"), INLA = "https://inla.r-inla-download.org/R/stable"), dep = TRUE)

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

# Construct temporal parts of model
start_year = model_config$YEAR_NAT_START
end_year = model_config$YEAR_NAT_END
n_years = (end_year-start_year + 1)
n_months = n_years*12

# load INLA regression data
inla_data <- read.csv('/mnt/efs/userdata/etan/mitn_outputs/outputs/data_prep/INLA/inla_dataset_reduced.csv')
# inla_data <- inla_data[seq(1,dim(inla_data)[1],2),]
inla_data <- inla_data[which(inla_data$access > 0),]

# Need to replicate most recent year and append to allow INLA to extrapolate to final year
latest_data <- inla_data[which(inla_data$monthidx == max(inla_data$monthidx)),]
latest_data$monthidx <- n_months
inla_data <- rbind(inla_data, latest_data)
inla_data$yearidx <- ((inla_data$monthidx-1) %/% 12)+1#*12

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





# generate temporal mesh
temporal_mesh_monthly <- inla.mesh.1d(seq(1,n_months,by=12),interval=c(1, n_months),degree=2)

# Make projection matrices
A_proj_monthly <- inla.spde.make.A(mesh = africa_spde, loc = coords,
                                   group = inla_data$monthidx,
                                   group.mesh = temporal_mesh_monthly)

S_index_monthly <- inla.spde.make.index(name = "field",
                                        n.spde = africa_mesh$n,
                                        n.group = temporal_mesh_monthly$m)

# Construct INLA stack for spatial
# Calculate optimum ihs theta and transformation for use response
# Use gap
epsilon <- 1e-5
for (i in 1:length(inla_data$use)){
  inla_data$p_use[i] <- p_transform(inla_data$use[i], inla_data$access[i], n = 2)
}

use_theta <- optimise(ihs_loglik, lower = 0.001, upper = 200, x = gap_emplogit(inla_data$p_use[which(abs(inla_data$p_use) > 0.001)]), maximum = TRUE)$maximum
res_use_gap <- ihs(gap_emplogit(inla_data$p_use), use_theta)

response_data <- list(res_use_gap = res_use_gap)

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

effects_data_monthly <- list(c(S_index_monthly, list(Intercept = 1)),cov_data)

africa_stack_monthly <- inla.stack(data = response_data,
                                   A = list(A_proj_monthly,1),
                                   effects = effects_data_monthly,
                                   tag = "use.data")

#############################
# Fit USE GAP model
#############################
print("Fitting Use gap spatio-temporal model...")
m1 <- inla(res_use_gap ~ -1 +# Intercept + 
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
             monthly_1 +
             monthly_2 +
             f(field, model = spde, group = field.group, 
               control.group = list(model = 'ar1') ),
           data = inla.stack.data(africa_stack_monthly, spde = africa_spde),
           family = "gaussian",
           control.predictor = list(A = inla.stack.A(africa_stack_monthly), compute = TRUE),
           control.compute = list(cpo = TRUE, dic = TRUE, config = TRUE), 
           control.inla = list(strategy = "adaptive", int.strategy = "eb"),
           verbose = TRUE)

print("Saving Use gap model outputs...")

save(africa_mesh, africa_spde, temporal_mesh_monthly, m1, use_theta, file = "/mnt/efs/userdata/etan/mitn_outputs/outputs/INLA/model1_use_complete_logis.RData")

print("Saved Use model Part")

summary(m1)

print("HURRAH! I FINISHED THE R SCRIPT THANK GOD")

