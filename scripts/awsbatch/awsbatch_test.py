from awsmap.workflow import Workflow, JuliaManifest

def create_workflow(wf: Workflow):

    wf.init(
        name = "MITN_run",
        base = JuliaManifest(manifest_path="env_files/Manifest.toml")
    )

    step1 = wf.add_task(
        name = "DataExtract",
        command = "julia scripts/DataProcessing/SurveyCleaning/DataEng_SurveyDataCleanup.jl",
        memory_mb = 8192
    )