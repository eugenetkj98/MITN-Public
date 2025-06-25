from awsmap.workflow import Workflow, JuliaManifest, DockerImage

def create_workflow(wf: Workflow):
    # load config.toml file. Access variables using config["value"]
    config = wf.load_toml_config("configs/batch_config.toml")
    model_config = wf.load_toml_config("configs/model_config.toml")

    wf.init(
        name = "Admin2_Scenario_Convolution",
        base = JuliaManifest(manifest_path="env_files/Manifest.toml")
    )

    step3 = wf.add_task(
        name = "Admin2_Scenario_Convolution",
        memory_mb = 131072,
        vcpus = 32,
        command = f"julia --threads={config["step3"]["ncpus"]} scripts/AuxAnalysis/admin2_convolution.jl"
    )



    