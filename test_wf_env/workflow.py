from awsmap.workflow import Workflow, JuliaManifest

def create_workflow(wf: Workflow):
    # load config.toml file. Access variables using config["value"]
    config = wf.load_toml_config("config.toml")

    wf.init(
        name = "HELLO",
        base = JuliaManifest("Manifest.toml")
    )

    step1 = wf.add_task(
        name = "DataExtraction",
        command = "julia --threads 2 testrun.jl",
        memory_mb = 16384,
        vcpus = 4
    )
    
    step2 = wf.add_task(
        name = "$task2_name",
        command = "Rscript step2.R",
        after = [step1]
    )

    step3 = wf.add_task(
        name = "step3",
        command = "julia --threads 5 analysis.jl ${ISO}", # ${YEAR} refers to the year variable. Injected via array_parameters
        after = [step2], # this set of steps only starts running after step1 completes successfully
        vcpus = 5,
        memory_mb = 8192,
        array_parameters = {
            "ISO": list(readcsv(config["gen"]["inputs"]["iso_list_path"]).ISO)
        }
    )