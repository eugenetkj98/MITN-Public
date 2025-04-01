Instructions to run current version of MITN model. v 1.0
Last Updated: 31/3/2025

1. National MITN
1.1. scripts/DataProcessing/SurveyCleaning/DataEng_SurveyDataCleanup.jl
1.2. scripts/DataProcessing/National SNF/ExtractHouseholdMetrics.jl
1.3. scripts/NationalMITN/run_analysis.jl
1.4. scripts/NationalMITN/generate_crop_access_draws.jl
1.5. scripts/DataProcesing/National SNF/SummariseNationalMITN_attrition.jl
1.6. scripts/DataProcesing/National SNF/SummariseNationalMITN_netage.jl

2. Subnational MITN
2.1. scripts/DataProcessing/Subnational SNF/Subnational_ExtractHouseholdMetrics.jl
2.2. scripts/SubnationalMITN/subnational_MITN_regression.jl
2.3. scripts/SubnationalMITN/generate_crop_access_draws.jl

3. Spatial Data Prep
3.1. scripts/DataProcessing/Spatial Disaggregation/spatial_prep_hh_data.jl
3.2. scripts/DataProcessing/Spatial Disaggregation/spatial_prep_pca_collin.jl
3.3. Transfer following files to INLA directory
    - covariate_normalisation_constants.csv
    - inla_dataset_reduced.csv
    - proj_matrix.csv
    - MAP_Regions_unclipped_5k.tif

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
