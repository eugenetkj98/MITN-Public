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
        name = "AD_HOC_Convergence_test",
        memory_mb = config["step3"]["memsize"],
        vcpus = config["step3"]["ncpus"],
        command = f"julia --threads={config["step3"]["ncpus"]} scripts/AuxAnalysis/TechPaper_Convergence/convergence_test.jl "+"${sample_idx}",
        array_parameters = {
            "sample_idx": list(range(1, 20 + 1)) 
            # "sample_idx": [6,11,13,14,15,16,17]
        }
    )
    # print("Completed Step 3 Succesfully!")



    