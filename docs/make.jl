using Pkg
include(pwd()*"/scripts/init_env.jl")
# Pkg.activate("docs/");
# Pkg.instantiate()
# Pkg.add("Missings")
Pkg.add("Documenter") # In case Documenter.jl is not installed
push!(LOAD_PATH,"src/");

# Function to recursively find all sub_directories and load all child directories
function recurse_loaddir(src_dir)
    push!(LOAD_PATH, src_dir)
    sub_dirs = src_dir.*readdir(src_dir)
    isdir_bool = isdir.(sub_dirs)
    isdir_idx = findall(isdir_bool)
    for idx in isdir_idx
        push!(LOAD_PATH,sub_dirs[idx]*"/")
        sum(isdir.(sub_dirs[idx]*"/".*readdir(sub_dirs[idx]*"/"))) > 0
        if sum(isdir.(sub_dirs[idx]*"/".*readdir(sub_dirs[idx]*"/"))) > 0
            subsub_dirs = sub_dirs[idx]*"/".*readdir(sub_dirs[idx]*"/")

            for subsub_dir in subsub_dirs
                recurse_loaddir(subsub_dir)
            end
        end
    end
end

# Load subdirectories from src folder
src_dir = pwd()*"/src/"
recurse_loaddir(src_dir)

# Load other packages and custom packages
using Documenter

using DateConversions,
        NetLoss,
        DataExtractions,
        NetCropModel,
        NetCropRegression,
        Logit,
        NetAccessModel,
        NetAccessRegression,
        PlottingFunctions

# Makefile

Description = "Model Description" => [
                "Net Crop SNF" => "lib/description_SNF.md",
                "Net Access" => "lib/description_access.md",
                "Spatial LGM" => "lib/description_spatial.md"
                ];
Pipeline = "Pipeline and Examples" => [
                "Data Extraction" => "lib/data_extractions.md",
                "Models" => "lib/models.md",
                "Regression" => "lib/regression.md",
                "Example Code" => "lib/example.md"
                ];

Model_IO = "Model IO" => [
                "Model Inputs" => "lib/inputs.md",
                "Model Outputs" => "lib/outputs.md"                        
                ];

Visualisation = "Visualisation" => [
                "Plotting Functions" => "lib/plotting_functions.md"                    
                ];

HelperFunctions = "Helper Functions" => [
                "Date Conversions" => "lib/date_conversions.md",
                "Net Attrition" => "lib/net_loss.md",
                "Logit Transformations" => "lib/logit.md"];

PAGES = [
    "SNF_Docs.md",
    Description,
    Pipeline,
    Model_IO,
    Visualisation,
    HelperFunctions
];

makedocs(sitename = "ITN Stock and Flow Documentation",
        remotes = nothing,
        format = Documenter.HTML(prettyurls = haskey(ENV, "CI")),
        pages = PAGES, pagesonly = true)