# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %%
using ProgressBars
using Plots
using CSV
using DataFrames

# %%
using DateConversions

# %%
survey_data = CSV.read("datasets/inla_dataset.csv", DataFrame)

YEAR_START = 2000
YEAR_END = 2023
n_months = (YEAR_END-YEAR_START+1)*12

# %%
@gif for i in ProgressBar(1:n_months)
    month, year_ref = monthidx_to_monthyear(i)
    year = year_ref + YEAR_START -1

    lats = survey_data[survey_data.monthidx .== i,"latitude"]
    longs = survey_data[survey_data.monthidx .== i,"longitude"]
    npc_gap = survey_data[survey_data.monthidx .== i,"npc_gap"]

    if length(lats) == 0
        scatter(markerstrokewidth = 0,
            xlims = (-20,60), ylims = (-35,35), clims = (-30,30),
            title = "NPC Gap\n$(year)-$(month)")
    else
        scatter(longs, lats, zcolor= npc_gap, markerstrokewidth = 0,
                xlims = (-20,60), ylims = (-35,35), clims = (-0.3,0.3),
                title = "NPC Gap\n$(year)-$(month)")
    end
end fps = 10
