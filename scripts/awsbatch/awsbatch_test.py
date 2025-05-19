from awsmap.workflow import Workflow, JuliaManifest

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
    step8 = wf.add_task(
        name = "Subnational_NPC_SNF_regression",
        memory_mb = config["step8"]["memsize"],
        vcpus = config["step8"]["ncpus"],
        command = f"julia --threads={config["step8"]["ncpus"]} scripts/SubnationalMITN/AWS_subnational_MITN_regression.jl "+"${ISO}",
        # after = [step5]
        array_parameters = {
            "ISO": model_config["ISO_LIST"]
        }
    )
    # print("Completed Step 8 Succesfully!")