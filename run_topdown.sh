#! /bin/bash
set -e

export DEFAULT_NUM_THREADS=32
export DISTRIBUTED_NUM_THREADS=5

julia --threads $DEFAULT_NUM_THREADS "scripts/RasterMaps/raster_timeseries_aggregation.jl" | tee -a outputs.txt

date | tee outputs.txt
# National MITN Regression

echo "Step 0: Data Extraction and Preparation" | tee -a outputs.txt
echo "Step 0.1: Survey Cleanup" | tee -a outputs.txt
julia --threads $DEFAULT_NUM_THREADS "scripts/DataProcessing/SurveyCleaning/DataEng_SurveyDataCleanup.jl" | tee -a outputs.txt

echo "Step 0.2: Extract Household Metrics for National" | tee -a outputs.txt
julia --threads $DEFAULT_NUM_THREADS "scripts/DataProcessing/National SNF/ExtractHouseholdMetrics.jl" | tee -a outputs.txt
echo "Completed data preparation from surveys"

echo "Step 1: Run National MITN regression" | tee -a outputs.txt
echo "Step 1.1: Net Crop extraction and regression" | tee -a outputs.txt
export JULIA_NUM_THREADS=$DISTRIBUTED_NUM_THREADS
julia scripts/NationalMITN/run_analysis_crop.jl | tee -a outputs.txt
export JULIA_NUM_THREADS=$DEFAULT_NUM_THREADS

echo "Step 1.2: Net Access extraction and regression"
julia --threads $DEFAULT_NUM_THREADS "scripts/NationalMITN/run_analysis_access.jl" | tee -a outputs.txt

echo "Step 1.3: Perform posterior samples of National MITN"
julia --threads $DEFAULT_NUM_THREADS "scripts/NationalMITN/generate_crop_access_draws.jl" | tee -a outputs.txt

echo "Step 1.4: Perform posterior samples of National MITN attrition parameters"
julia --threads $DEFAULT_NUM_THREADS "scripts/DataProcessing/National SNF/SummariseNationalMITN_attrition.jl" | tee -a outputs.txt

echo "Step 1.5: Perform posterior samples of National MITN net age demography"
julia --threads $DEFAULT_NUM_THREADS "scripts/DataProcessing/National SNF/SummariseNationalMITN_netage.jl" | tee -a outputs.txt
echo "Completed National MITN regression and sampling"

echo "Step 2: Run Subnational model"
echo "Step 2.1: Extract Household Metrics for Subnational"
julia --threads $DEFAULT_NUM_THREADS "scripts/DataProcessing/Subnational SNF/Subnational_ExtractHouseholdMetrics.jl" | tee -a outputs.txt

echo "Step 2.2: Subnational MITN regression"
export JULIA_NUM_THREADS=$DISTRIBUTED_NUM_THREADS
julia "scripts/SubnationalMITN/subnational_MITN_regression.jl" | tee -a outputs.txt
export JULIA_NUM_THREADS=$DEFAULT_NUM_THREADS

echo "Step 2.3: Perform posterior samples for Subnational MITN"
julia --threads $DEFAULT_NUM_THREADS "scripts/SubnationalMITN/generate_crop_access_draws.jl" | tee -a outputs.txt
echo "Completed Subnational MITN regression and sampling"

echo "Step 3: Prepare covariate data for spatial disaggregation"
echo "Step 3.1: Extract covariate data from rasters"
julia --threads $DEFAULT_NUM_THREADS "scripts/DataProcessing/Spatial Disaggregation/spatial_prep_hh_data.jl" | tee -a outputs.txt
echo "Step 3.2: Deal with collinearities"
julia --threads $DEFAULT_NUM_THREADS "scripts/DataProcessing/Spatial Disaggregation/spatial_prep_pca_collin.jl" | tee -a outputs.txt
echo "Completed covariate data prep."

echo "Step 4: Run INLA to regress INLA models and sample NPC and Access deviation rasters"
echo "Step 4.1: Train INLA NPC model"
Rscript "scripts/INLA/INLA_regression_npc.R"
echo "Step 4.2: Train INLA Access model"
Rscript "scripts/INLA/INLA_regression_access.R"
echo "Completed INLA model regression for NPC and access."
echo "Step 4.3: Beginning posterior sampling for mean NPC and Access deviation  rasters"
Rscript "vis_mondels_npc.R"
echo "Sampled NPC deviation Rasters."
Rscript "vis_models_access.R"
echo "Sampled Access deviation Rasters."

echo "Step 5: Construct adjust NPC and Access rasters and prep dataset to regress INLA Use model"
echo "Step 5.1: Prepare adjusted INLA training dataset (normalised w.r.t. to adjusted Access raster)"
julia --threads $DEFAULT_NUM_THREADS "scripts/DataProcessing/SpatialDisaggregation/spatial_prep_use_data.jl" | tee -a outputs.txt
echo "Constructed adjusted rasters."
echo "Step 5.2: Do INLA Use model regression"
Rscript "scripts/INLA/INLA_regression_use.R"
echo "Completed INLA model regression for use."
echo "Step 5.3: Beginning posterior sampling for use deviation rasters."
Rscript "vis_models_use.R"
echo "Complete Use deviation Rasters"
echo "Step 5.4: Construct Use rasters"
julia --threads $DEFAULT_NUM_THREADS "scripts/RasterMaps/generate_final_rasters.jl" | tee -a outputs.txt

export DEFAULT_NUM_THREADS=32
export DISTRIBUTED_NUM_THREADS=5
julia --threads $DEFAULT_NUM_THREADS "scripts/RasterMaps/raster_timeseries_aggregation.jl"
echo "Completed constructing all required rasters"

echo "Step 6: Extract time series for ITN coverage"
echo "Step 6.1: National population weighted estimates"
julia --threads $DEFAULT_NUM_THREADS "scripts/RasterMaps/raster_to_nat_timeseries.jl" | tee -a outputs.txt
echo "Step 6.2: Subnational population weighted estimates"
julia --threads $DEFAULT_NUM_THREADS "scripts/RasterMaps/raster_to_subnat_timeseries.jl" | tee -a outputs.txt
echo "Step 6.3: Compare model predictions against true survey values to calculate RMSE"
julia --threads $DEFAULT_NUM_THREADS "scripts/DataProcessing/RMSE Validation/compile_model_estimates.jl" | tee -a outputs.txt
