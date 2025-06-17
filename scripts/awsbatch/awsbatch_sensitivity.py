from awsmap.workflow import Workflow, JuliaManifest, DockerImage

def create_workflow(wf: Workflow):
    # load config.toml file. Access variables using config["value"]
    config = wf.load_toml_config("configs/batch_config.toml")
    model_config = wf.load_toml_config("configs/model_config.toml")

    wf.init(
        name = "MITN_run",
        base = JuliaManifest(manifest_path="env_files/Manifest.toml")
    )

    # print("Ad-Hoc Test: Run Convergence Samples\n")
    step3 = wf.add_task(
        name = "AD_HOC_Sensitivity_test",
        memory_mb = config["step3"]["memsize"],
        vcpus = config["step3"]["ncpus"],
        command = f"julia --threads={config["step3"]["ncpus"]} scripts/AuxAnalysis/TechPaper_Convergence/sensitivity_analysis.jl "+"${data_size}",
        array_parameters = {
            "data_size": list(range(2, 39+2,2)) 
        }
    )
    # print("Completed Step 3 Succesfully!")



    