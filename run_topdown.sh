#! /bin/bash
set -e

export DEFAULT_NUM_THREADS=32
export DISTRIBUTED_NUM_THREADS=5

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
julia --threads $DEFAULT_NUM_THREADS "scripts/NationalMITN/run_analysis_access.jl"

echo "Step 1.3: Perform posterior samples of National MITN"
julia --threads $DEFAULT_NUM_THREADS "scripts/NationalMITN/generate_crop_access_draws.jl"

echo "Step 1.4: Perform posterior samples of National MITN attrition parameters"
julia --threads $DEFAULT_NUM_THREADS "scripts/DataProcesing/National SNF/SummariseNationalMITN_attrition.jl"

echo "Step 1.5: Perform posterior samples of National MITN net age demography"
julia --threads $DEFAULT_NUM_THREADS "scripts/DataProcesing/National SNF/SummariseNationalMITN_netage.jl"
echo "Completed National MITN regression and sampling"

echo "Step 2: Run Subnational model"
echo "Step 2.1: Extract Household Metrics for Subnational"
julia --threads $DEFAULT_NUM_THREADS "scripts/DataProcessing/Subnational SNF/Subnational_ExtractHouseholdMetrics.jl"
echo "Step 2.2: Subnational MITN regression"
export JULIA_NUM_THREADS=$DISTRIBUTED_NUM_THREADS
julia "scripts/SubnationalMITN/subnational_MITN_regression.jl"
export JULIA_NUM_THREADS=$DEFAULT_NUM_THREADS
echo "Step 2.3: Perform posterior samples for Subnational MITN"
julia --threads $DEFAULT_NUM_THREADS "scripts/SubnationalMITN/generate_crop_access_draws.jl"
echo "Completed Subnational MITN regression and sampling"

echo "Step 3: Prepare covariate data for spatial disaggregation"
echo "Step 3.1: Extract covariate data from rasters"
julia --threads $DEFAULT_NUM_THREADS "scripts/DataProcessing/Spatial Disaggregation/spatial_prep_hh_data.jl"
echo "Step 3.2: Deal with collinearities"
julia --threads $DEFAULT_NUM_THREADS "scripts/DataProcessing/Spatial Disaggregation/spatial_prep_pca_collin.jl"
echo "Completed covariate data prep."


4. Spatial Disaggregation (DIRECTORIES TO BE UPDATED IN FUTURE MERGE)
4.1. scripts/INLA/INLA_regression_npc.R
4.2. scripts/INLA/INLA_regression_access.R
4.3. scripts/INLA/INLA_regression_use.R
4.4. Optional Validation Checks
    - scripts/INLA/validation_models_annual_npc.R
    - scripts/INLA/validation_models_annual_access.R
    - scripts/INLA/validation_models_monthly_use.R
4.5. scripts/INLA/vis_models_npc.R
4.6. scripts/INLA/vis_models_access.R
4.7. scripts/INLA/vis_models_use.R
4.8. Export generated .tif rasters for final raster construction

5. ITN Coverage Raster Construction and Metrics
5.1. scripts/RasterMaps/generate_final_rasters.jl
5.2. scripts/RasterMaps/raster_to_nat_timeseries.jl
5.3. scripts/RasterMaps/raster_to_subnat_timeseries.jl