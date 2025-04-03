"""
Author: Eugene Tan
Date Created: 4/2/2025
Last Updated: 4/2/2025
Plotter functions to visualise timeseries predictions of BV model and MITN model at raster level. 
With confidence intervals
"""
# %% Prep environment and subdirectories
include(pwd()*"/scripts/init_env.jl")

# %% Import filenames and directories from config file
include(pwd()*"/scripts/dir_configs.jl")

# %% Import packages
using CairoMakie
using JLD2
using CSV
using DataFrames
using Colors
using StatsBase

# %% Import Dataset for Plotting
model_nat_est_timeseries = load(OUTPUT_DIR*"coverage_timeseries/adj_nat_model_coverage.jld2")
survey_pred_data = CSV.read(OUTPUT_DIR*"coverage_timeseries/adj_model_prediction_comparisons.csv", DataFrame)


CairoMakie.plot(survey_pred_data.use, survey_pred_data.mitn_use, markersize = 3)
CairoMakie.plot(survey_pred_data.use, survey_pred_data.bv_use, markersize = 3)

idxs = findall(.!isnan.(survey_pred_data.use) .& 
                .!isnan.(survey_pred_data.bv_use) .&
                (survey_pred_data.admin_level .== 0))
sqrt(mean((survey_pred_data.use[idxs] .- survey_pred_data.mitn_use[idxs]).^2))

sqrt(mean((survey_pred_data.use[idxs] .- survey_pred_data.bv_use[idxs]).^2))


# %% Get ISO List
ISO_list = intersect(model_nat_est_timeseries["filt_ISOs"], survey_pred_data.ISO)
n_countries = length(ISO_list)
tiling_dim = (5,9)

# %% Set theme for plot
set_theme!(theme_light())
lw = 1.5
error_alpha = 0.1
linecol_bv = colorant"#2754B6"
linecol_mitn = colorant"#BA263A"
scat_size = 8
scat_col_1 = colorant"#5FBFA2"
scat_col_2 = colorant"#EC9F05"

# %% Calculate Time Labels
YEAR_START, YEAR_END = model_nat_est_timeseries["YEAR_LIST"][[1,end]]
YEAR_VALS = (YEAR_START:(1/12):YEAR_END+1)[1:end-1]

# %% Create Master Figure and subfigures
fig_npc = Figure(size = (2560,1300))
fig_access = Figure(size = (2560,1300))
fig_hh_use = Figure(size = (2560,1300))
fig_pop_use = Figure(size = (2560,1300))

npc_subfigs_list = []
access_subfigs_list = []
hh_use_subfigs_list = []
pop_use_subfigs_list = []
for i in 1:n_countries
    row_idx = ((i-1)÷tiling_dim[2])+1
    col_idx = mod(i-1,tiling_dim[2])+1
    npc_subfig = fig_npc[row_idx,col_idx] = GridLayout()
    access_subfig = fig_access[row_idx,col_idx] = GridLayout()
    hh_use_subfig = fig_hh_use[row_idx,col_idx] = GridLayout()
    pop_use_subfig = fig_pop_use[row_idx,col_idx] = GridLayout()
    push!(npc_subfigs_list, npc_subfig)
    push!(access_subfigs_list, access_subfig)
    push!(hh_use_subfigs_list, hh_use_subfig)
    push!(pop_use_subfigs_list, pop_use_subfig)
end

# %% Add a Master title
Label(fig_npc[0,1:tiling_dim[2]], "National NPC", fontsize = 30)
Label(fig_access[0,1:tiling_dim[2]], "National Household Access", fontsize = 30)
Label(fig_hh_use[0,1:tiling_dim[2]], "National Household Use Rate", fontsize = 30)
Label(fig_pop_use[0,1:tiling_dim[2]], "National Population Use Rate", fontsize = 30)

# %% For Plot NPC model Time series for each country
for i in 1:n_countries
    ISO = ISO_list[i]
    ax = Axis(npc_subfigs_list[i][1,1], title = "$ISO",
                xlabel = "Year", ylabel = "NPC",
                xticks = YEAR_START:2:YEAR_END, xticklabelrotation=π/2)
    CairoMakie.ylims!(ax, (-0.05,1.05))

    # Plot Survey observed values
    filt_data = survey_pred_data[(survey_pred_data.ISO .== ISO).&
                                        (survey_pred_data.admin_level .== 0),:]
                                        CairoMakie.plot!(ax, YEAR_VALS[filt_data.monthidx], filt_data.npc, color = scat_col_1, markersize = scat_size,
            label = "Survey")

    # Plot BV model time series prediction and errors

    CairoMakie.lines!(ax, YEAR_VALS,model_nat_est_timeseries["bv_nat_npc"][i,:,2], linewidth = lw, color = linecol_bv,
            label = "BV")
    CairoMakie.band!(ax, YEAR_VALS, model_nat_est_timeseries["bv_nat_npc"][i,:,1], model_nat_est_timeseries["bv_nat_npc"][i,:,3],
            color = linecol_bv, alpha = error_alpha)

    # Plot MITN model time series prediction and errors
    CairoMakie.lines!(ax, YEAR_VALS,model_nat_est_timeseries["mitn_nat_npc"][i,:,2], linewidth = lw, color = linecol_mitn,
            label = "MITN")
    CairoMakie.band!(ax, YEAR_VALS, model_nat_est_timeseries["mitn_nat_npc"][i,:,1], model_nat_est_timeseries["mitn_nat_npc"][i,:,3],
            color = linecol_mitn, alpha = error_alpha)

    # Add legend
    CairoMakie.axislegend(ax, position = :lt)
end

fig_npc

# %% For Plot Access model Time series for each country
for i in 1:n_countries
    ISO = ISO_list[i]
    ax = Axis(access_subfigs_list[i][1,1], title = "$ISO",
                xlabel = "Year", ylabel = "Household Access",
                xticks = YEAR_START:2:YEAR_END, xticklabelrotation=π/2)
    CairoMakie.ylims!(ax, (-0.05,1.05))

    # Plot Survey observed values
    filt_data = survey_pred_data[(survey_pred_data.ISO .== ISO).&
                                        (survey_pred_data.admin_level .== 0),:]
    if !(ISO ∈ ["BWA","DJI","GNQ","ETH","SOM","SSD"]) # Temporary fix. Exclude some countries from plotting survey access. There shouldn't be any available, (problem with data summary code)
        CairoMakie.plot!(ax, YEAR_VALS[filt_data.monthidx], filt_data.access, color = scat_col_1, markersize = scat_size,
                label = "Survey")
    end

    # Plot BV model time series prediction and errors

    CairoMakie.lines!(ax, YEAR_VALS,model_nat_est_timeseries["bv_nat_access"][i,:,2], linewidth = lw, color = linecol_bv,
            label = "BV")
    CairoMakie.band!(ax, YEAR_VALS, model_nat_est_timeseries["bv_nat_access"][i,:,1], model_nat_est_timeseries["bv_nat_access"][i,:,3],
            color = linecol_bv, alpha = error_alpha)

    # Plot MITN model time series prediction and errors
    CairoMakie.lines!(ax, YEAR_VALS,model_nat_est_timeseries["mitn_nat_access"][i,:,2], linewidth = lw, color = linecol_mitn,
            label = "MITN")
    CairoMakie.band!(ax, YEAR_VALS, model_nat_est_timeseries["mitn_nat_access"][i,:,1], model_nat_est_timeseries["mitn_nat_access"][i,:,3],
            color = linecol_mitn, alpha = error_alpha)

    # Add legend
    CairoMakie.axislegend(ax, position = :lt)
end

fig_access

# %% For Plot Use model Time series for each country

# LREG dirty conversion function from proper household use to population use
f_use(x) = 1.26813*x + 0.050331
inv_f_use(x) = (x - 0.050331)/1.26813

# For household use
for i in 1:n_countries
    ISO = ISO_list[i]
    ax = Axis(hh_use_subfigs_list[i][1,1], title = "$ISO",
                xlabel = "Year", ylabel = "Household Use",
                xticks = YEAR_START:2:YEAR_END, xticklabelrotation=π/2)
    CairoMakie.ylims!(ax, (-0.05,1.05))

    # Plot Survey observed values
    filt_data = survey_pred_data[(survey_pred_data.ISO .== ISO).&
                                        (survey_pred_data.admin_level .== 0),:]
    # CairoMakie.plot!(ax, YEAR_VALS[filt_data.monthidx], 2 .* filt_data.use, color = scat_col_1, markersize = scat_size,
    #         label = "Survey (HH)")
    CairoMakie.plot!(ax, YEAR_VALS[filt_data.monthidx], filt_data.use, color = scat_col_1, markersize = scat_size,
            label = "Survey (HH)")

    # Plot BV model time series prediction and errors

    # lines!(ax, YEAR_VALS, 2 .* inv_f_use.(model_nat_est_timeseries["bv_nat_use"][i,:,2]), linewidth = lw, color = linecol_bv,
    #         label = "BV")
    # band!(ax, YEAR_VALS, 2 .* inv_f_use.(model_nat_est_timeseries["bv_nat_use"][i,:,1]), 2 .* inv_f_use.(model_nat_est_timeseries["bv_nat_use"][i,:,3]),
    #         color = linecol_bv, alpha = error_alpha)

    CairoMakie.lines!(ax, YEAR_VALS, model_nat_est_timeseries["bv_nat_use"][i,:,2], linewidth = lw, color = linecol_bv,
            label = "BV")
    CairoMakie.band!(ax, YEAR_VALS, model_nat_est_timeseries["bv_nat_use"][i,:,1], model_nat_est_timeseries["bv_nat_use"][i,:,3],
            color = linecol_bv, alpha = error_alpha)

    # Plot MITN model time series prediction and errors
    # lines!(ax, YEAR_VALS, 2 .* model_nat_est_timeseries["mitn_nat_use"][i,:,2], linewidth = lw, color = linecol_mitn,
    #         label = "MITN")
    # band!(ax, YEAR_VALS, 2 .* model_nat_est_timeseries["mitn_nat_use"][i,:,1], 2 .* model_nat_est_timeseries["mitn_nat_use"][i,:,3],
    #         color = linecol_mitn, alpha = error_alpha)

    CairoMakie.lines!(ax, YEAR_VALS, model_nat_est_timeseries["mitn_nat_use"][i,:,2], linewidth = lw, color = linecol_mitn,
            label = "MITN")
    CairoMakie.band!(ax, YEAR_VALS, model_nat_est_timeseries["mitn_nat_use"][i,:,1], model_nat_est_timeseries["mitn_nat_use"][i,:,3],
            color = linecol_mitn, alpha = error_alpha)

    # Add legend
    CairoMakie.axislegend(ax, position = :lt, groupgap = 0)
end

fig_hh_use


# For population use
for i in 1:n_countries
    ISO = ISO_list[i]
    ax = Axis(pop_use_subfigs_list[i][1,1], title = "$ISO",
                xlabel = "Year", ylabel = "Household Use",
                xticks = YEAR_START:2:YEAR_END, xticklabelrotation=π/2)
    CairoMakie.ylims!(ax, (-0.05,1.05))

    # Plot Survey observed values
    filt_data = survey_pred_data[(survey_pred_data.ISO .== ISO).&
                                        (survey_pred_data.admin_level .== 0),:]

    # plot!(ax, YEAR_VALS[filt_data.monthidx], f_use.(filt_data.use), color = scat_col_2, markersize = scat_size,
    #         label = "Survey (Pop.)")
    CairoMakie.plot!(ax, YEAR_VALS[filt_data.monthidx], filt_data.use, color = scat_col_2, markersize = scat_size,
            label = "Survey (Pop.)")

    # Plot BV model time series prediction and errors

    CairoMakie.lines!(ax, YEAR_VALS, model_nat_est_timeseries["bv_nat_use"][i,:,2], linewidth = lw, color = linecol_bv,
            label = "BV")
    CairoMakie.band!(ax, YEAR_VALS, model_nat_est_timeseries["bv_nat_use"][i,:,1], model_nat_est_timeseries["bv_nat_use"][i,:,3],
            color = linecol_bv, alpha = error_alpha)

    # Plot MITN model time series prediction and errors
    # lines!(ax, YEAR_VALS, f_use.(model_nat_est_timeseries["mitn_nat_use"][i,:,2]), linewidth = lw, color = linecol_mitn,
    #         label = "MITN")
    # band!(ax, YEAR_VALS, f_use.(model_nat_est_timeseries["mitn_nat_use"][i,:,1]), f_use.(model_nat_est_timeseries["mitn_nat_use"][i,:,3]),
    #         color = linecol_mitn, alpha = error_alpha)

    CairoMakie.lines!(ax, YEAR_VALS, model_nat_est_timeseries["mitn_nat_use"][i,:,2], linewidth = lw, color = linecol_mitn,
            label = "MITN")
    CairoMakie.band!(ax, YEAR_VALS, model_nat_est_timeseries["mitn_nat_use"][i,:,1], model_nat_est_timeseries["mitn_nat_use"][i,:,3],
            color = linecol_mitn, alpha = error_alpha)

    # Add legend
    CairoMakie.axislegend(ax, position = :lt, groupgap = 0)
end

fig_pop_use

# %% Save plots
mkpath("output_plots/ITN Coverage summaries/timeseries/national/")
save("output_plots/ITN Coverage summaries/timeseries/national/nat_npc.pdf", fig_npc)
save("output_plots/ITN Coverage summaries/timeseries/national/nat_access.pdf", fig_access)
save("output_plots/ITN Coverage summaries/timeseries/national/nat_hh_use.pdf", fig_hh_use)
save("output_plots/ITN Coverage summaries/timeseries/national/nat_pop_use.pdf", fig_pop_use)

# %% Make national RMSE Plots
# Storage variable for RMSE
bv_npc_RMSE = (Vector{Float64}(undef, n_countries) .= NaN)
mitn_npc_RMSE = (Vector{Float64}(undef, n_countries) .= NaN)
bv_access_RMSE = (Vector{Float64}(undef, n_countries) .= NaN)
mitn_access_RMSE = (Vector{Float64}(undef, n_countries) .= NaN)
bv_hh_use_RMSE = (Vector{Float64}(undef, n_countries) .= NaN)
mitn_hh_use_RMSE = (Vector{Float64}(undef, n_countries) .= NaN)
bv_pop_use_RMSE = (Vector{Float64}(undef, n_countries) .= NaN)
mitn_pop_use_RMSE = (Vector{Float64}(undef, n_countries) .= NaN)


# Calculate RMSE values for each country
for i in 1:n_countries
    ISO = ISO_list[i]

    if (ISO ∈ ["BWA","DJI","GNQ","ETH","SOM","SSD"]) # Temporary fix. Exclude some countries from plotting survey access. There shouldn't be any available, (problem with data summary code)
        continue
    end

    # For NPC
    npc_filt_data = survey_pred_data[(survey_pred_data.ISO .== ISO).&
                                    (survey_pred_data.admin_level .== 0) .&
                                    (.!isnan.(survey_pred_data.npc)) .&
                                    (.!isnan.(survey_pred_data.bv_npc)) .&
                                    (.!isnan.(survey_pred_data.mitn_npc)),:]

    bv_npc_RMSE[i] = sqrt(mean((npc_filt_data.bv_npc .- npc_filt_data.npc).^2))
    mitn_npc_RMSE[i] = sqrt(mean((npc_filt_data.mitn_npc .- npc_filt_data.npc).^2))

    # For Access
    access_filt_data = survey_pred_data[(survey_pred_data.ISO .== ISO).&
                                    (survey_pred_data.admin_level .== 0) .&
                                    (.!isnan.(survey_pred_data.access)) .&
                                    (.!isnan.(survey_pred_data.bv_access)) .&
                                    (.!isnan.(survey_pred_data.mitn_access)),:]

    bv_access_RMSE[i] = sqrt(mean((access_filt_data.bv_access .- access_filt_data.access).^2))
    mitn_access_RMSE[i] = sqrt(mean((access_filt_data.mitn_access .- access_filt_data.access).^2))

    # For household use
    hh_use_filt_data = survey_pred_data[(survey_pred_data.ISO .== ISO).&
                                    (survey_pred_data.admin_level .== 0) .&
                                    (.!isnan.(survey_pred_data.use)) .&
                                    (.!isnan.(survey_pred_data.bv_use)) .&
                                    (.!isnan.(survey_pred_data.mitn_use)),:]

    # bv_hh_use_RMSE[i] = sqrt(mean(((2 .* inv_f_use.(hh_use_filt_data.bv_use)) .- (2 .* hh_use_filt_data.use)).^2))
    # mitn_hh_use_RMSE[i] = sqrt(mean(((2 .* hh_use_filt_data.mitn_use) .- (2 .* hh_use_filt_data.use)).^2))

    # bv_pop_use_RMSE[i] = sqrt(mean((hh_use_filt_data.bv_use .- f_use.(hh_use_filt_data.use)).^2))
    # mitn_pop_use_RMSE[i] = sqrt(mean((f_use.(hh_use_filt_data.mitn_use) .- f_use.(hh_use_filt_data.use)).^2))

    bv_hh_use_RMSE[i] = sqrt(mean(((hh_use_filt_data.bv_use) .- (hh_use_filt_data.use)).^2))
    mitn_hh_use_RMSE[i] = sqrt(mean(((hh_use_filt_data.mitn_use) .- (hh_use_filt_data.use)).^2))

    bv_pop_use_RMSE[i] = sqrt(mean((hh_use_filt_data.bv_use .- hh_use_filt_data.use).^2))
    mitn_pop_use_RMSE[i] = sqrt(mean((hh_use_filt_data.mitn_use .- hh_use_filt_data.use).^2))
end

# Calculate RMSE of full national level dataset
npc_filt_data = survey_pred_data[(survey_pred_data.admin_level .== 0) .&
                                (.!isnan.(survey_pred_data.npc)) .&
                                (.!isnan.(survey_pred_data.bv_npc)) .&
                                (.!isnan.(survey_pred_data.mitn_npc)),:]
bv_full_npc_RMSE = sqrt(mean((npc_filt_data.bv_npc .- npc_filt_data.npc).^2))
mitn_full_npc_RMSE = sqrt(mean((npc_filt_data.mitn_npc .- npc_filt_data.npc).^2))

access_filt_data = survey_pred_data[(survey_pred_data.admin_level .== 0) .&
                                (.!isnan.(survey_pred_data.access)) .&
                                (.!isnan.(survey_pred_data.bv_access)) .&
                                (.!isnan.(survey_pred_data.mitn_access)) .& 
                                (survey_pred_data.ISO .!= "BWA") .& 
                                (survey_pred_data.ISO .!= "DJI") .& 
                                (survey_pred_data.ISO .!= "GNQ") .& 
                                (survey_pred_data.ISO .!= "ETH") .& 
                                (survey_pred_data.ISO .!= "SOM") .& 
                                (survey_pred_data.ISO .!= "SSD"),:]
bv_full_access_RMSE = sqrt(mean((access_filt_data.bv_access .- access_filt_data.access).^2))
mitn_full_access_RMSE = sqrt(mean((access_filt_data.mitn_access .- access_filt_data.access).^2))

use_filt_data = survey_pred_data[(survey_pred_data.admin_level .== 0) .&
                                (.!isnan.(survey_pred_data.use)) .&
                                (.!isnan.(survey_pred_data.bv_use)) .&
                                (.!isnan.(survey_pred_data.mitn_use)),:]
# bv_full_hh_use_RMSE = sqrt(mean(( (2 .* inv_f_use.(use_filt_data.bv_use)) .- (2 .* use_filt_data.use)).^2))
# mitn_full_hh_use_RMSE = sqrt(mean(((2 .* use_filt_data.mitn_use) .- (2 .* use_filt_data.use)).^2))
bv_full_hh_use_RMSE = sqrt(mean(( use_filt_data.bv_use .- use_filt_data.use).^2))
mitn_full_hh_use_RMSE = sqrt(mean((use_filt_data.mitn_use .- use_filt_data.use).^2))


# bv_full_pop_use_RMSE = sqrt(mean((use_filt_data.bv_use .- f_use.(use_filt_data.use)).^2))
# mitn_full_pop_use_RMSE = sqrt(mean((f_use.(use_filt_data.mitn_use) .- f_use.(use_filt_data.use)).^2))
bv_full_pop_use_RMSE = sqrt(mean((use_filt_data.bv_use .- use_filt_data.use).^2))
mitn_full_pop_use_RMSE = sqrt(mean((use_filt_data.mitn_use .- use_filt_data.use).^2))

# %% Calculate indicator whether BV or MITN performed better
npc_rmse_indic = (Vector{Float64}(undef, n_countries) .= NaN)
access_rmse_indic = (Vector{Float64}(undef, n_countries) .= NaN)
hh_use_rmse_indic = (Vector{Float64}(undef, n_countries) .= NaN)
pop_use_rmse_indic = (Vector{Float64}(undef, n_countries) .= NaN)

for i in 1:n_countries
    if isnan(bv_npc_RMSE[i])
        if isnan(mitn_npc_RMSE[i])
            continue
        else
            npc_rmse_indic[i] = 2
        end
    else
        if bv_npc_RMSE[i] < mitn_npc_RMSE[i]
            npc_rmse_indic[i] = 1
        else 
            npc_rmse_indic[i] = 2
        end
    end

    if isnan(bv_access_RMSE[i])
        if isnan(mitn_access_RMSE[i])
            continue
        else
            access_rmse_indic[i] = 2
        end
    else
        if bv_access_RMSE[i] < mitn_access_RMSE[i]
            access_rmse_indic[i] = 1
        else 
            access_rmse_indic[i] = 2
        end
    end

    if isnan(bv_hh_use_RMSE[i])
        if isnan(mitn_hh_use_RMSE[i])
            continue
        else
            hh_use_rmse_indic[i] = 2
        end
    else
        if bv_hh_use_RMSE[i] < mitn_hh_use_RMSE[i]
            hh_use_rmse_indic[i] = 1
        else 
            hh_use_rmse_indic[i] = 2
        end
    end

    if isnan(bv_pop_use_RMSE[i])
        if isnan(mitn_pop_use_RMSE[i])
            continue
        else
            pop_use_rmse_indic[i] = 2
        end
    else
        if bv_pop_use_RMSE[i] < mitn_pop_use_RMSE[i]
            pop_use_rmse_indic[i] = 1
        else 
            pop_use_rmse_indic[i] = 2
        end
    end
end

# %% Make national RMSE plots
# Construct master plot
fig_rmse = Figure(size = (1920,1080))
labelsize = 17
titlesize = 22
scat_size = 13
sep = 0.15
highlight_width = 20
highlight_alpha = 0.5
alpha_mag = 5
highlight_cols = [linecol_bv, linecol_mitn]

subfig_npc_rmse = fig_rmse[1,1:8] = GridLayout()
subfig_access_rmse = fig_rmse[2,1:8] = GridLayout()
subfig_hh_use_rmse = fig_rmse[3,1:8] = GridLayout()
subfig_pop_use_rmse = fig_rmse[4,1:8] = GridLayout()
subfig_hh_pop_rmse = fig_rmse[5,1:8] = GridLayout()
subfig_full_rmse = fig_rmse[1:5,9] = GridLayout()

# Plot NPC rmse
ax_npc = Axis(subfig_npc_rmse[1,1], title = "NPC", titlesize = titlesize,
                xlabel = "Country", ylabel = "RMSE", 
                xlabelsize = labelsize, ylabelsize = labelsize, 
                xticks = (1:n_countries, ISO_list))
CairoMakie.xlims!(ax_npc, 0, n_countries + 1)
CairoMakie.ylims!(ax_npc, -0.01, 0.32)

### RMSE values
CairoMakie.plot!(ax_npc, (1:n_countries).-sep, bv_npc_RMSE, markersize = scat_size,
        label = "BV", color = linecol_bv)
CairoMakie.plot!(ax_npc, (1:n_countries).+sep, mitn_npc_RMSE, markersize = scat_size,
        label = "MITN", color = linecol_mitn)

### Highlights for better model
for i in 1:n_countries
    indic_val = npc_rmse_indic[i]
    if !isnan(indic_val)
        min_error_val = min(bv_npc_RMSE[i], mitn_npc_RMSE[i])
        max_error_val = max(bv_npc_RMSE[i], mitn_npc_RMSE[i])
        alpha_scale = min(alpha_mag*(max_error_val - min_error_val),1)#/max_error_val
        vlines!(ax_npc, i, color = highlight_cols[Int(indic_val)],
                linewidth = highlight_width, 
                alpha = highlight_alpha*alpha_scale)
    end
end

### Add legend
axislegend(ax_npc, position = :rt, groupgap = 0)

# Plot Access rmse
ax_access = Axis(subfig_access_rmse[1,1], title = "Household Access", titlesize = titlesize,
                xlabel = "Country", ylabel = "RMSE", 
                xlabelsize = labelsize, ylabelsize = labelsize, 
                xticks = (1:n_countries, ISO_list))
CairoMakie.xlims!(ax_access, 0, n_countries + 1)
CairoMakie.ylims!(ax_access, -0.01, 0.42)

### RMSE values
CairoMakie.plot!(ax_access, (1:n_countries).-sep, bv_access_RMSE, markersize = scat_size,
        label = "BV", color = linecol_bv)
CairoMakie.plot!(ax_access, (1:n_countries).+sep, mitn_access_RMSE, markersize = scat_size,
        label = "MITN", color = linecol_mitn)


### Highlights for better model
for i in 1:n_countries
    indic_val = access_rmse_indic[i]
    if !isnan(indic_val)
        min_error_val = min(bv_access_RMSE[i], mitn_access_RMSE[i])
        max_error_val = max(bv_access_RMSE[i], mitn_access_RMSE[i])
        alpha_scale = min(alpha_mag*(max_error_val - min_error_val),1)#/max_error_val
        vlines!(ax_access, i, color = highlight_cols[Int(indic_val)],
                linewidth = highlight_width, 
                alpha = highlight_alpha*alpha_scale)
    end
end

### Add legend
axislegend(ax_access, position = :rt, groupgap = 0)

# Plot Household Use rmse
ax_hh_use = Axis(subfig_hh_use_rmse[1,1], title = "Household Use", titlesize = titlesize,
                xlabel = "Country", ylabel = "RMSE", 
                xlabelsize = labelsize, ylabelsize = labelsize, 
                xticks = (1:n_countries, ISO_list))
CairoMakie.xlims!(ax_hh_use, 0, n_countries + 1)
CairoMakie.ylims!(ax_hh_use, -0.01, 0.42)

### RMSE values
CairoMakie.plot!(ax_hh_use, (1:n_countries).-sep, bv_hh_use_RMSE, markersize = scat_size,
        label = "BV", color = linecol_bv)
CairoMakie.plot!(ax_hh_use, (1:n_countries).+sep, mitn_hh_use_RMSE, markersize = scat_size,
        label = "MITN", color = linecol_mitn)


### Highlights for better model
for i in 1:n_countries
    indic_val = hh_use_rmse_indic[i]
    if !isnan(indic_val)
        min_error_val = min(bv_hh_use_RMSE[i], mitn_hh_use_RMSE[i])
        max_error_val = max(bv_hh_use_RMSE[i], mitn_hh_use_RMSE[i])
        alpha_scale = min(alpha_mag*(max_error_val - min_error_val),1)#/max_error_val
        vlines!(ax_hh_use, i, color = highlight_cols[Int(indic_val)],
                linewidth = highlight_width, 
                alpha = highlight_alpha*alpha_scale)
    end
end

### Add legend
axislegend(ax_hh_use, position = :rt, groupgap = 0)

# Plot Population Use rmse
ax_pop_use = Axis(subfig_pop_use_rmse[1,1], title = "Population Use", titlesize = titlesize,
                xlabel = "Country", ylabel = "RMSE", 
                xlabelsize = labelsize, ylabelsize = labelsize, 
                xticks = (1:n_countries, ISO_list))
CairoMakie.xlims!(ax_pop_use, 0, n_countries + 1)
CairoMakie.ylims!(ax_hh_use, -0.01, 0.42)

### RMSE values
CairoMakie.plot!(ax_pop_use, (1:n_countries).-sep, bv_pop_use_RMSE, markersize = scat_size,
        label = "BV", color = linecol_bv)
CairoMakie.plot!(ax_pop_use, (1:n_countries).+sep, mitn_pop_use_RMSE, markersize = scat_size,
        label = "MITN", color = linecol_mitn)


### Highlights for better model
for i in 1:n_countries
    indic_val = pop_use_rmse_indic[i]
    if !isnan(indic_val)
        min_error_val = min(bv_pop_use_RMSE[i], mitn_pop_use_RMSE[i])
        max_error_val = max(bv_pop_use_RMSE[i], mitn_pop_use_RMSE[i])
        alpha_scale = min(alpha_mag*(max_error_val - min_error_val),1)#/max_error_val
        vlines!(ax_pop_use, i, color = highlight_cols[Int(indic_val)],
                linewidth = highlight_width, 
                alpha = highlight_alpha*alpha_scale)
    end
end

### Add legend
axislegend(ax_pop_use, position = :rt, groupgap = 0)

# Plot comparative difference in gain between using BV and MITN in both hh_use and pop_use
ax_hh_pop = Axis(subfig_hh_pop_rmse[1,1], title = "Household vs Population Use Performance", titlesize = titlesize,
                xlabel = "Country", ylabel = "Normalised ΔRMSE", 
                xlabelsize = labelsize, ylabelsize = labelsize, 
                xticks = (1:n_countries, ISO_list))
CairoMakie.xlims!(ax_hh_pop, 0, n_countries + 1)
CairoMakie.ylims!(ax_hh_pop, -0.02, 1.02)

### RMSE values
bv_transfer_ratios = abs.(bv_pop_use_RMSE .- bv_hh_use_RMSE)./((bv_pop_use_RMSE .+ bv_hh_use_RMSE .+ mitn_pop_use_RMSE .+ mitn_hh_use_RMSE)./4)
mitn_transfer_ratios = abs.(mitn_pop_use_RMSE .- mitn_hh_use_RMSE)./((bv_pop_use_RMSE .+ bv_hh_use_RMSE .+ mitn_pop_use_RMSE .+ mitn_hh_use_RMSE)./4)
CairoMakie.plot!(ax_hh_pop, (1:n_countries).-sep, bv_transfer_ratios, markersize = scat_size,
        label = "BV", color = linecol_bv)
CairoMakie.plot!(ax_hh_pop, (1:n_countries).+sep, mitn_transfer_ratios, markersize = scat_size,
        label = "MITN", color = linecol_mitn)


### Highlights for better model
for i in 1:n_countries
    indic_val = NaN
    if isnan(bv_transfer_ratios[i])
        if isnan(mitn_transfer_ratios[i])
            continue
        else
            indic_val = 2
        end
    else
        if bv_transfer_ratios[i] < mitn_transfer_ratios[i]
            indic_val = 1
        else
            indic_val = 2
        end
    end
    if !isnan(indic_val)
        if indic_val == 1
            vlines!(ax_hh_pop, i, color = highlight_cols[Int(indic_val)],
                    linewidth = highlight_width, 
                    alpha = highlight_alpha*(mitn_transfer_ratios[i] - bv_transfer_ratios[i]))
        else
            vlines!(ax_hh_pop, i, color = highlight_cols[Int(indic_val)],
                    linewidth = highlight_width, 
                    alpha = highlight_alpha*(bv_transfer_ratios[i] - mitn_transfer_ratios[i]))
        end
    end
end

### Add legend
CairoMakie.axislegend(ax_hh_pop, position = :rt, groupgap = 0)

# Add side plot for metrics from full dataset
ax_full = Axis(subfig_full_rmse[1,1], title = "Full Survey Data", titlesize = titlesize,
                xlabel = "Metric", ylabel = "RMSE", 
                xlabelsize = labelsize, ylabelsize = labelsize, 
                xticks = (1:4, ["NPC", "Access", "HH\nUse", "Pop\nUse"]),
                yticks = (0:0.05:0.25))
CairoMakie.ylims!(ax_full, -0.002, 0.26)

CairoMakie.plot!(ax_full, (1:4) .- sep, [bv_full_npc_RMSE, bv_full_access_RMSE, bv_full_hh_use_RMSE, bv_full_pop_use_RMSE],
        markersize = scat_size, label = "BV", color = linecol_bv)
CairoMakie.plot!(ax_full, (1:4 .+ sep), [mitn_full_npc_RMSE, mitn_full_access_RMSE, mitn_full_hh_use_RMSE, mitn_full_pop_use_RMSE],
        markersize = scat_size, label = "MITN", color = linecol_mitn)

### Add legend
axislegend(ax_full, position = :rb, groupgap = 0)
fig_rmse
# %% Save RMSE Plot

mkpath("output_plots/ITN Coverage summaries/")
save("output_plots/ITN Coverage summaries/nat_RMSE.pdf", fig_rmse)