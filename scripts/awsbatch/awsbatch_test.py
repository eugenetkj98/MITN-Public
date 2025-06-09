from awsmap.workflow import Workflow, JuliaManifest, DockerImage

def create_workflow(wf: Workflow):
    # load config.toml file. Access variables using config["value"]
    config = wf.load_toml_config("configs/batch_config.toml")
    model_config = wf.load_toml_config("configs/model_config.toml")

    wf.init(
        name = "MITN_run",
        base = JuliaManifest(manifest_path="env_files/Manifest.toml")
    )

    # print("Commencing Step 1: Household Survey Cleanup\n")
    # step1 = wf.add_task(
    #     name = "Household_Survey_Cleanup",
    #     memory_mb = config["step1"]["memsize"],
    #     vcpus = config["step1"]["ncpus"],
    #     command = f"julia --threads={config["step1"]["ncpus"]} scripts/DataProcessing/SurveyCleaning/DataEng_SurveyDataCleanup.jl"
    # )
    # print("Completed Step 1 Succesfully!")

    # print("Commencing Step 2: Extract Household NPC Metrics\n")
    # step2 = wf.add_task(
    #     name = "Extract_Household_NPC_Metrics_National",
    #     memory_mb = config["step2"]["memsize"],
    #     vcpus = config["step2"]["ncpus"],
    #     command = f"julia --threads={config["step2"]["ncpus"]} \"scripts/DataProcessing/National SNF/ExtractHouseholdMetrics.jl\" ",
    #     after = [step1]
    # )
    # print("Completed Step 2 Succesfully!")

    # print("Commencing Step 3: Net Crop Model Regression\n")
    # step3 = wf.add_task(
    #     name = "National_NPC_SNF_Regression",
    #     memory_mb = config["step3"]["memsize"],
    #     vcpus = config["step3"]["ncpus"],
    #     command = f"julia --threads={config["step3"]["ncpus"]} scripts/NationalMITN/AWS_run_analysis_crop.jl "+"${ISO}",
    #     # after = [step2]
    #     array_parameters = {
    #         "ISO": model_config["ISO_LIST"]
    #     }
    # )
    # print("Completed Step 3 Succesfully!")

    # print("Commencing Step 4: Net Access Model Regression\n")
    # step4 = wf.add_task(
    #     name = "National_Access_SNF_Regression",
    #     memory_mb = config["step4"]["memsize"],
    #     vcpus = config["step4"]["ncpus"],
    #     command = f"julia --threads={config["step4"]["ncpus"]} scripts/NationalMITN/AWS_run_analysis_access.jl",
    #     # after = [step3]
    # )
    # print("Completed Step 4 Succesfully!")

    # print("Commencing Step 5: Net Crop SNF Posterior Draws\n")
    # step5 = wf.add_task(
    #     name = "National_NPC_SNF_Draws",
    #     memory_mb = config["step5"]["memsize"],
    #     vcpus = config["step5"]["ncpus"],
    #     command = f"julia --threads={config["step5"]["ncpus"]} scripts/NationalMITN/AWS_generate_crop_access_draws.jl "+"${ISO}",
    #     # after = [step4]
    #     array_parameters = {
    #         "ISO": model_config["ISO_LIST"]
    #     }
    # )
    # print("Completed Step 5 Succesfully!")

    # print("Commencing Step 6: Export Posterior attrition and net age draws for national\n")
    # step6a = wf.add_task(
    #     name = "National_NPC_SNF_Summarise_Attrition",
    #     memory_mb = config["step6"]["memsize"],
    #     vcpus = config["step6"]["ncpus"],
    #     command = f"julia --threads={config["step6"]["ncpus"]} \"scripts/DataProcessing/National SNF/SummariseNationalMITN_attrition.jl\" ",
    #     # after = [step5]
    # )
    # print("Completed Step 6 Succesfully!")

    # print("Commencing Step 6: Export Posterior attrition and net age draws for national\n")
    # step6b = wf.add_task(
    #     name = "National_NPC_SNF_Summarise_NetAge",
    #     memory_mb = config["step6"]["memsize"],
    #     vcpus = config["step6"]["ncpus"],
    #     command = f"julia --threads={config["step6"]["ncpus"]} \"scripts/DataProcessing/National SNF/SummariseNationalMITN_netage.jl\" ",
    #     # after = [step5]
    # )
    # print("Completed Step 6 Succesfully!")

    # print("Commencing Step 7: Extract and prep subnational household data\n")
    # step7 = wf.add_task(
    #     name = "Subnational_Extract_Household_Data",
    #     memory_mb = config["step7"]["memsize"],
    #     vcpus = config["step7"]["ncpus"],
    #     command = f"julia --threads={config["step7"]["ncpus"]} \"scripts/DataProcessing/Subnational SNF/Subnational_ExtractHouseholdMetrics.jl\" ",
    #     # after = [step5]
    # )
    # print("Completed Step 7 Succesfully!")

    # print("Commencing Step 8: Subnational MITN SNF Regression\n")
    # step8 = wf.add_task(
    #     name = "Subnational_NPC_SNF_regression",
    #     memory_mb = config["step8"]["memsize"],
    #     vcpus = config["step8"]["ncpus"],
    #     command = f"julia --threads={config["step8"]["ncpus"]} scripts/SubnationalMITN/AWS_subnational_MITN_regression.jl "+"${ISO}",
    #     # after = [step7]
    #     array_parameters = {
    #         "ISO": model_config["ISO_LIST"]
    #     }
    # )
    # print("Completed Step 8 Succesfully!")

    # print("Commencing Step 9: Subnational MITN SNF Draws\n")
    # step9 = wf.add_task(
    #     name = "Subnational_NPC_SNF_draws",
    #     memory_mb = config["step9"]["memsize"],
    #     vcpus = config["step9"]["ncpus"],
    #     command = f"julia --threads={config["step9"]["ncpus"]} scripts/SubnationalMITN/AWS_generate_crop_access_draws.jl "+"${ISO}",
    #     # after = [step8]
    #     array_parameters = {
    #         "ISO": model_config["ISO_LIST"]
    #     }
    # )
    # print("Completed Step 9 Succesfully!")


    # print("Commencing Step 9b: Extract Net Age summaries\n")
    # step9b = wf.add_task(
    #     name = "Subnational_Mean_Net_Age_draws",
    #     memory_mb = config["step9b"]["memsize"],
    #     vcpus = config["step9b"]["ncpus"],
    #     command = f"julia --threads={config["step9b"]["ncpus"]} \"scripts/DataProcessing/Subnational SNF/extract_netage_timeseries.jl\" ",
    #     # after = [step9]
    # )
    # print("Completed Step 9b Succesfully!")

    # print("Commencing Step 10: Prep covariates data for INLA regression\n")
    # step10 = wf.add_task(
    #     name = "ExtractCovariates",
    #     memory_mb = config["step10"]["memsize"],
    #     vcpus = config["step10"]["ncpus"],
    #     command = f"julia --threads={config["step10"]["ncpus"]} \"scripts/DataProcessing/Spatial Disaggregation/spatial_prep_hh_data.jl\" ",
    #     # after = [step9]
    # )
    # print("Completed Step 10 Succesfully!")

    # print("Commencing Step 11: Normalise covariates and deal with collinearities\n")
    # step11 = wf.add_task(
    #     name = "NormaliseCovariates_Collin",
    #     memory_mb = config["step11"]["memsize"],
    #     vcpus = config["step11"]["ncpus"],
    #     command = f"julia --threads={config["step11"]["ncpus"]} \"scripts/DataProcessing/Spatial Disaggregation/spatial_prep_pca_collin.jl\" ",
    #     after = [step10]
    # )
    # print("Completed Step 11 Succesfully!")

    ##########################################
    # NEED TO GET R CODE READY!!
    ##########################################

    month_year_argslist = []
    for year in range(model_config["YEAR_NAT_START"], model_config["YEAR_NAT_END"]+1):
        for month in range(1,13):
                month_year_argslist.append(f"{year} {month}")

    # step12 = wf.add_task(
    #     name = "INLA_NPC_Regression",
    #     command = "Rscript scripts/INLA/INLA_regression_npc.R",
    #     base = DockerImage("rocker/geospatial:4.4.3"),
    #     memory_mb = config["step12"]["memsize"],
    #     vcpus = config["step12"]["ncpus"],
    #     # after = [step11]
    # )

    # step12b = wf.add_task(
    #     name = "INLA_Access_Regression",
    #     command = "Rscript scripts/INLA/INLA_regression_access.R",
    #     base = DockerImage("rocker/geospatial:4.4.3"),
    #     memory_mb = config["step12"]["memsize"],
    #     vcpus = config["step12"]["ncpus"],
    #     # after = [step11]
    # )

    # step12b = wf.add_task(
    #     name = "INLA_Deployment_Regression",
    #     command = "Rscript scripts/INLA/INLA_regression_deployment.R",
    #     base = DockerImage("rocker/geospatial:4.4.3"),
    #     memory_mb = config["step12"]["memsize"],
    #     vcpus = config["step12"]["ncpus"],
    #     # after = [step11]
    # )

    # step12c = wf.add_task(
    #     name = "INLA_Use_Regression",
    #     command = "Rscript scripts/INLA/INLA_regression_use.R",
    #     base = DockerImage("rocker/geospatial:4.4.3"),
    #     memory_mb = config["step12"]["memsize"],
    #     vcpus = config["step12"]["ncpus"],
    #     # after = [step11]
    # )

    # step13 = wf.add_task(
    #     name = "INLA_NPC_Sampling",
    #     command = "Rscript scripts/INLA/vis_models_npc.R ${YEAR}",
    #     base = DockerImage("rocker/geospatial:4.4.3"),
    #     memory_mb = config["step13"]["memsize"],
    #     vcpus = config["step13"]["ncpus"],
    #     # after = [step12],
    #     array_parameters = {
    #         "YEAR": list(range(model_config["YEAR_NAT_START"], model_config["YEAR_NAT_END"] + 1)) 
    #     }
    # )

    # step13b = wf.add_task(
    #     name = "INLA_Access_Sampling",
    #     command = "Rscript scripts/INLA/vis_models_access.R ${YEAR}",
    #     base = DockerImage("rocker/geospatial:4.4.3"),
    #     memory_mb = config["step13"]["memsize"],
    #     vcpus = config["step13"]["ncpus"],
    #     after = [step12b],
    #     array_parameters = {
    #         "YEAR": list(range(model_config["YEAR_NAT_START"], model_config["YEAR_NAT_END"] + 1)) 
    #     }
    # )

    # step13b = wf.add_task(
    #     name = "INLA_Deployment_Sampling",
    #     command = "Rscript scripts/INLA/vis_models_deployment.R ${YEAR}",
    #     base = DockerImage("rocker/geospatial:4.4.3"),
    #     memory_mb = config["step13"]["memsize"],
    #     vcpus = config["step13"]["ncpus"],
    #     after = [step12b],
    #     array_parameters = {
    #         "YEAR": list(range(model_config["YEAR_NAT_START"], model_config["YEAR_NAT_END"] + 1)) 
    #     }
    # )

    month_year_argslist = []
    for year in range(model_config["YEAR_NAT_START"], model_config["YEAR_NAT_END"]+1):
        for month in range(1,13):
                month_year_argslist.append(f"{year} {month}")    

    step13c = wf.add_task(
        name = "INLA_Use_Sampling",
        command = "Rscript scripts/INLA/vis_models_use.R ${monthyear_arg}",
        base = DockerImage("rocker/geospatial:4.4.3"),
        memory_mb = config["step13"]["memsize"],
        vcpus = config["step13"]["ncpus"],
        # after = [step12c],
        array_parameters = {
            "monthyear_arg": month_year_argslist  
        }
    )

    # step13c = wf.add_task(
    #     name = "INLA_Use_Sampling",
    #     command = "Rscript scripts/INLA/vis_models_use_temp.R ${year}",
    #     base = DockerImage("rocker/geospatial:4.4.3"),
    #     memory_mb = config["step13"]["memsize"],
    #     vcpus = config["step13"]["ncpus"],
    #     # after = [step12c],
    #     array_parameters = {
    #         "YEAR": list(range(model_config["YEAR_NAT_START"], model_config["YEAR_NAT_END"] + 1)) 
    #     }
    # )


    ##########################################
    # NEED TO GET R CODE READY!!
    ##########################################
        
    # print("Commencing Step 18: Construct final rasters using INLA outputs\n")
    step18 = wf.add_task(
        name = "Construct_Final_ITN_Rasters",
        memory_mb = config["step18"]["memsize"],
        vcpus = config["step18"]["ncpus"],
        command = f"julia --threads={config["step18"]["ncpus"]} scripts/RasterMaps/samples_based_run/generate_final_rasters_samples.jl "+"${monthyear_arg}",
        after = [step13c],
        array_parameters = {
            "monthyear_arg": month_year_argslist  
        }
    )
    # print("Completed Step 18 Successfully!")

    # print("Commencing Step 19: Extract timeseries from generated rasters \n")
    step19 = wf.add_task(
        name = "Extract_Raster_Timeseries",
        memory_mb = config["step19"]["memsize"],
        vcpus = config["step19"]["ncpus"],
        command = f"julia --threads={config["step19"]["ncpus"]} scripts/RasterMaps/samples_based_run/raster_timeseries_aggregation_samples.jl "+"${monthyear_arg}",
        after = [step18],
        array_parameters = {
            "monthyear_arg": month_year_argslist  
        }
    )
    # print("Completed Step 19 Successfully!")

    # print("Commencing Step 19b: Aggregate raster time series .csv \n")
    step19b = wf.add_task(
        name = "Aggregate_Raster_Timeseries",
        memory_mb = config["step19b"]["memsize"],
        vcpus = config["step19b"]["ncpus"],
        command = f"julia --threads={config["step19b"]["ncpus"]} scripts/RasterMaps/samples_based_run/join_timeseries_aggregates_samples.jl",
        after = [step19]
    )
    # print("Completed Step 19 Successfully!")

    # print("Commencing Step 20: Construct mean monthly rasters \n")
    step20 = wf.add_task(
        name = "Construct_Mean_Monthly_Rasters",
        memory_mb = config["step20"]["memsize"],
        vcpus = config["step20"]["ncpus"],
        command = f"julia --threads={config["step20"]["ncpus"]} scripts/RasterMaps/samples_based_run/calculate_mean_rasters.jl "+"${monthyear_arg}",
        after = [step18],
        array_parameters = {
            "monthyear_arg": month_year_argslist  
        }
    )
    # print("Completed Step 20 Successfully!")

    # print("Commencing Step 20b: Construct mean annual rasters \n")
    step20b = wf.add_task(
        name = "Construct_Raster_Annual_Rasters",
        memory_mb = config["step20b"]["memsize"],
        vcpus = config["step20b"]["ncpus"],
        command = f"julia --threads={config["step20b"]["ncpus"]} scripts/RasterMaps/samples_based_run/calculate_annual_aggregate_rasters.jl "+"${YEAR}",
        after = [step20],
        array_parameters = {
            "YEAR": list(range(model_config["YEAR_NAT_START"], model_config["YEAR_NAT_END"] + 1))
        }
    )
    print("Completed Step 20b Successfully!")




    