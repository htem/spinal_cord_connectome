module v

using JLD2
using Optim
using Random
using GLMakie
using StatsBase
using Parameters
using LinearAlgebra
using ProgressMeter
using DelimitedFiles
using OrderedCollections

@with_kw mutable struct Params

    T::Float64 = 85
    dt::Float64 = 0.01
    tstep::Int = floor(Int, T/dt)
    times::Array = 0:dt:T-dt

    stimonset::Float64 = 10
    stimdurat::Float64 = 75
    stimampli::Float64 = 1

    normalize_conn::Bool  = true

    #### from E
    β_EE::Float64 = 1
    β_IE::Float64 = 1

    #### from I
    β_EI::Float64 = 2
    β_II::Float64 = 2

    #####
    excit_cell_type::Array = [
                              "R",
                              "module-cltmr",
                              "module-adltmr",
                              "module-np",
                              "module-noci",
                              "V1",
                              "V2",
                              "P1",
                              "P2-poly",
                              "P2-noci",
                              "P2-cold",
                             ]

    #####
    inhib_cell_type::Array = [
                              "islet-cltmr",
                              "islet-adltmr",
                              "islet-mrgd",
                              "islet-adhtmr",
                             ]


    #####
    senso_cell_type::Array = [
                              "sst",
                              "cltmr",
                              "adltmr",
                              "mrgd",
                              "adhtmr",
                              "chtmr-p",
                              "ccold",
                             ]

    h0Dict::Dict{String,Bool} = Dict(
                       "sst"     => 0,
                       "cltmr"   => 0,
                       "adltmr"  => 0,
                       "mrgd"    => 0,
                       "adhtmr"  => 0,
                       "chtmr-p" => 0,
                       "ccold"   => 0,
                      )

    I0_pulseDict::Dict{String,Bool} = Dict(
                 "R"             => 0,
                 "module-cltmr"  => 0,
                 "module-adltmr" => 0,
                 "module-np"     => 0,
                 "module-noci"   => 0,
                 "V1"            => 0,
                 "V2"            => 0,

                 "P1"            => 0,
                 "P2-poly"       => 0,
                 "P2-noci"       => 0,
                 "P2-cold"       => 0,

                 "islet-cltmr"   => 0,
                 "islet-adltmr"  => 0,
                 "islet-mrgd"    => 0,
                 "islet-adhtmr"  => 0,
                )

    I0Dict = Dict(
                 "R"             => 0.0,
                 "module-cltmr"  => 0.0,
                 "module-adltmr" => 0.0,
                 "module-np"     => 0.0,
                 "module-noci"   => 0.0,
                 "V1"            => 0.0,
                 "V2"            => 0.0,

                 "P1"            => 0.0,
                 "P2-poly"       => 0.0,
                 "P2-noci"       => 0.0,
                 "P2-cold"       => 0.0,

                 "islet-cltmr"   => 0.0,
                 "islet-adltmr"  => 0.0,
                 "islet-mrgd"    => 0.0,
                 "islet-adhtmr"  => 0.0,
                 "islet-V"       => 0.0,
                )

    τDict = Dict(
                 "R"             => 1.0,
                 "module-cltmr"  => 1.0,
                 "module-adltmr" => 1.0,
                 "module-np"     => 1.0,
                 "module-noci"   => 1.0,
                 "V1"            => 1.0,
                 "V2"            => 1.0,

                 "P1"            => 1.0,
                 "P2-poly"       => 1.0,
                 "P2-noci"       => 1.0,
                 "P2-cold"       => 1.0,

                 "islet-cltmr"   => 1.0,
                 "islet-adltmr"  => 1.0,
                 "islet-mrgd"    => 1.0,
                 "islet-adhtmr"  => 1.0,
                 "islet-V"       => 1.0,
                )

    cell_type_name_mapping = OrderedDict(
                                         #{{{
                                         "R"    => "R",
                                         "module-cltmr"   => "C-LTMR mod",
                                         "module-np" => "C-HTMR/Heat mod",
                                         "module-adltmr"   => "Aδ-LTMR mod",
                                         "module-noci" => "Aδ-HTMR mod",
                                         "V1"   => "V1",
                                         "V2"   => "V2",
                                         "U"    => "unknown",
                                         "P1"   => "P1",
                                         "P2-poly"   => "P2-poly",
                                         "P2-noci"   => "P2-noci",
                                         "P2-cold"   => "P2-cold",
                                         "islet-cltmr"  => rich("islet",superscript("C-LTMR")),
                                         "islet-mrgd"   => rich("islet",superscript("C-HTMR/Heat")),
                                         "islet-adltmr" => rich("islet",superscript("Aδ-LTMR")),
                                         "islet-adhtmr" => rich("PSI-IhN",superscript("Aδ-HTMR")),
                                         "islet-V"      => "inhib vertical",
                                         "cltmr"   => "C-LTMR",
                                         "mrgd"    => rich("C-HTMR/Heat (MrgprD/A3/B4", superscript("+"), ")"),
                                         "sst"     => rich("C-HTMR/Heat (SST", superscript("+"), ")"),
                                         "chtmr-p" => "Peptidergic C-Nociceptor",
                                         "adltmr"  => "Aδ-LTMR",
                                         "adhtmr"  => "Aδ-HTMR",
                                         "ccold"   => "C-Cold",
                                         #}}}
                                        )

    colorDict = Dict(
                     #{{{
                     "R" => "#3F5482",
                     "module-cltmr" => "#7B2479",
                     "module-np" => "#00a168",
                     "module-adltmr" => "#D67300",
                     "module-noci" => "#008AA1",
                     "V1" => "#718600",
                     "V2" => "#730F1F",
                     "P1" => "#4D4182",
                     "P2-poly" => "#B4B753",
                     "P2-noci" => "#BB2D3B",
                     "P2-cold" => "#2E6CB2",
                     "islet-cltmr" => "#0864c3",
                     "islet-mrgd" => "#c36708",
                     "islet-adltmr" => "#c30864",
                     "islet-adhtmr" => "#08c367",
                     "islet-V" => "#c2c308",
                     #}}}
                    )
end

dp = Params()


##############################################################
function read_data(;p = dp)
    #{{{

    J = readdlm("data/J_recurrent.txt")
    W = readdlm("data/J_feedforward.txt")

    WpreDict = load("data/preDict_feedforward.jld2", "data")
    WpostDict = load("data/postDict_feedforward.jld2", "data")

    JpreDict = load("data/preDict_recurrent.jld2", "data")
    JpostDict = load("data/postDict_recurrent.jld2", "data")

    return J, W, JpreDict, JpostDict, WpreDict, WpostDict
    #}}}
end


##############################################################
function relu(x)
    return max(0, x)
end


##############################################################
function remap_dict_ranges(dict, kept_indices)
    #{{{
    index_map = Dict(old_idx => new_idx for (new_idx, old_idx) in enumerate(kept_indices))

    remapped_dict = OrderedDict()

    for (key, old_range) in dict
        old_indices = collect(old_range)
        new_indices = [index_map[idx] for idx in old_indices]
        
        if !isempty(new_indices)
            remapped_dict[key] = minimum(new_indices):maximum(new_indices)
        end
    end
    #}}}
    return remapped_dict
end


##############################################################
function set_connectivity(; p = dp)
    #{{{
    J0, W0, JpreDict, JpostDict, WpreDict, WpostDict = read_data(p = p)

    ############## filter to only include specified cell types #############
    valid_neuron_types = vcat(p.excit_cell_type, p.inhib_cell_type)
    valid_sensory_types = p.senso_cell_type

    JpreDict_filtered  = filter(q -> q.first in valid_neuron_types, JpreDict)
    JpostDict_filtered = filter(q -> q.first in valid_neuron_types, JpostDict)
    WpreDict_filtered  = filter(q -> q.first in valid_sensory_types, WpreDict)
    WpostDict_filtered = filter(q -> q.first in valid_neuron_types, WpostDict)

    neuron_indices_pre = sort(vcat([collect(inds) for inds in values(JpreDict_filtered)]...))
    neuron_indices_post = sort(vcat([collect(inds) for inds in values(JpostDict_filtered)]...))
    sensory_indices = sort(vcat([collect(inds) for inds in values(WpreDict_filtered)]...))

    J0 = J0[neuron_indices_post, neuron_indices_pre]
    W0 = W0[neuron_indices_post, sensory_indices]

    JpreDict  = remap_dict_ranges(JpreDict_filtered, neuron_indices_pre)
    JpostDict = remap_dict_ranges(JpostDict_filtered, neuron_indices_post)
    WpreDict  = remap_dict_ranges(WpreDict_filtered, sensory_indices)
    WpostDict = remap_dict_ranges(WpostDict_filtered, neuron_indices_post)

    ############# inhibitory (-) #############
    for pre_ct in p.inhib_cell_type
        J0[:,JpreDict[pre_ct]] .*= -1
    end

    ############# normalizing connectivity #############
    if p.normalize_conn
        indegree = dropdims(sum(abs.(J0),dims=2),dims=2) + dropdims(sum(abs.(W0),dims=2),dims=2)
        if any(iszero, indegree)
            indegree[indegree .== 0] .= 1
        end
        J = J0 ./ indegree
        W = W0 ./ indegree
    else
        J = J0
        W = W0
    end
    #}}}
    return J, W, JpreDict, JpostDict, WpreDict, WpostDict
end


##############################################################
function reduce_connectivity(J, preDict, postDict; figmode=false)
    #{{{
    J_bar = zeros(length(postDict), length(preDict))
    for (j, jtype) in enumerate(keys(preDict))
        for (i, itype) in enumerate(keys(postDict))
            J_bar[i,j] = sum(J[postDict[itype], preDict[jtype]]) / length(postDict[itype])
        end
    end

    preDict_bar = OrderedDict{String, Int}()
    for (i,k) in enumerate(keys(preDict))
        preDict_bar[k] = i
    end
    postDict_bar = OrderedDict{String, Int}()
    for (i,k) in enumerate(keys(postDict))
        postDict_bar[k] = i
    end

    if figmode
        fig = Figure()
        ax = Axis(fig[1,1])
        hm = heatmap!(ax, 1:size(J_bar,2), 1:size(J_bar,1), J_bar', colormap=:hot)
        Colorbar(fig[1,2], hm)
        ax.yticks = (1:length(postDict), [k for k in keys(postDict)])
        ax.xticks = (1:length(preDict), [k for k in keys(preDict)])
        ax.xticklabelrotation = π/4
    else
        return J_bar, preDict_bar, postDict_bar
    end
    #}}}
    return fig
end


##############################################################
function runsim(; p = dp)
    #{{{
    #####################################################
    J_full, W_full, JpreDict_full, JpostDict_full, WpreDict_full, WpostDict_full = set_connectivity(p = p)
    J, JpreDict, JpostDict = reduce_connectivity(J_full, JpreDict_full, JpostDict_full)
    W, WpreDict, WpostDict = reduce_connectivity(W_full, WpreDict_full, WpostDict_full)

    ###### E/I scaling #####
    for (pre_types, post_types, β) in [
        (p.excit_cell_type, p.excit_cell_type, p.β_EE),  # E→E
        (p.excit_cell_type, p.inhib_cell_type, p.β_IE),  # E→I
        (p.inhib_cell_type, p.excit_cell_type, p.β_EI),  # I→E
        (p.inhib_cell_type, p.inhib_cell_type, p.β_II)   # I→I
    ]
        for pre_ct in pre_types, post_ct in post_types
            J[JpostDict[post_ct], JpreDict[pre_ct]] *= β
        end
    end
    #####################################################

    ϕ = relu

    N = size(J, 1)  # number of neuron types
    M = size(W, 2)  # number of sensory types

    neuron_types = collect(keys(JpreDict))
    sensory_types = collect(keys(WpreDict))

    excit_indices = [JpreDict[ct] for ct in p.excit_cell_type if ct in neuron_types]
    inhib_indices = [JpreDict[ct] for ct in p.inhib_cell_type if ct in neuron_types]

    x = zeros(Float64, N)
    r = zeros(Float64, N)
    ηE = zeros(Float64, N)
    ηI = zeros(Float64, N)
    I0 = zeros(Float64, N)
    I0_pulse = zeros(Float64, N)

    h0 = zeros(Float64, M)
    hs = zeros(Float64, p.tstep, M)
    xs = zeros(Float64, p.tstep, N)
    rs = zeros(Float64, p.tstep, N)
    ηEs = zeros(Float64, p.tstep, N)
    ηIs = zeros(Float64, p.tstep, N)
    I0_pulses = zeros(Float64, p.tstep, N)

    ########## fill in inputs ##########
    for st in p.senso_cell_type
        active = p.h0Dict[st]
        idx = WpreDict[st]
        h0[idx] = active ? p.stimampli : 0.0
    end
    
    start_idx = round(Int, p.stimonset/p.dt) + 1
    end_idx   = round(Int, (p.stimonset + p.stimdurat)/p.dt)
    # write the same h0 pattern into hs between start and end
    for i = 1:M
        hs[start_idx:end_idx, i] .= h0[i]
    end

    ########## fill in bias current ##########
    for cell_type in neuron_types
        idx = JpreDict[cell_type]
        I0[idx] = p.I0Dict[cell_type]
        I0_pulse[idx] = p.I0_pulseDict[cell_type] ? p.stimampli : 0.0
    end
    for i = 1:N
        I0_pulses[start_idx:end_idx, i] .= I0_pulse[i]
    end

    ########### fill in dtoverτ ##########
    τ = zeros(Float64, N)
    for cell_type in neuron_types
        idx = JpreDict[cell_type]
        τ[idx] = p.τDict[cell_type]
    end
    dtoverτ = p.dt ./ τ

    ########## split J by columns (excitatory vs inhibitory) ##########
    JfromE = J[:, excit_indices]
    JfromI = J[:, inhib_indices]

    for i = 1:p.tstep
        ηE   = JfromE * r[excit_indices]
        ηI   = JfromI * r[inhib_indices]
        ηffw = W * hs[i,:]
        total_input = ηE + ηI + ηffw + I0 + I0_pulses[i,:]
        r = r + dtoverτ .* (-r + ϕ.(total_input))
        rs[i,:] = r
        ηEs[i,:] = ηE
        ηIs[i,:] = ηI
    end

    peak_rs = [maximum(rs[:,i]) for i = 1:size(rs,2)]

    ss_idx = min(end_idx, p.tstep)
    steady_rs = [rs[ss_idx, i] for i = 1:size(rs,2)]

    F = Dict()
    F[:J] = J
    F[:W] = W
    F[:xs] = xs
    F[:rs] = rs
    F[:hs] = hs
    F[:ηEs] = ηEs
    F[:ηIs] = ηIs
    F[:peak_rs] = peak_rs
    F[:steady_rs] = steady_rs
    F[:JpreDict] = JpreDict
    F[:JpostDict] = JpostDict
    F[:WpreDict] = WpreDict
    F[:WpostDict] = WpostDict

    last_window = round(Int, 1.0 / p.dt)
    drift = maximum(abs.(rs[end, :] - rs[end - last_window, :]))
    is_unstable = any(isnan.(rs)) || 
                  any(isinf.(rs)) || 
                  maximum(abs.(rs)) > 1e4 ||
                  drift > 1e-3
    F[:is_unstable] = maximum(abs.(rs)) > 1e4

    #}}}
    return F
end


########################
function plot_stimulations_sidebyside(; p = dp, 
                                      senso_type="cltmr", 
                                      show_types=["P1", "P2-poly", "P2-noci", "P2-cold"],
                                      fig = nothing,
                                      content_layout = nothing)
    #{{{
    standalone = isnothing(fig)
    if standalone
        set_my_theme!(fontsize=5, xlabelpadding=2, ylabelpadding=2, xticklabelpad=0, linewidth=0.5)
        fig = Figure(size=cm2px(4, 2))
        content_layout = GridLayout()
        fig[1, 1] = content_layout
    end

    axes = [Axis(content_layout[1, i]) for i = 1:length(show_types)]
    
    for key in keys(p.h0Dict)
        p.h0Dict[key] = 0
    end
    p.h0Dict[senso_type] = 1

    F = runsim(p = p)

    peak_sum = 0.0
    for cell_type in show_types
        cell_idx = F[:JpostDict][cell_type]
        peak_sum += F[:steady_rs][cell_idx] 
    end
    if peak_sum == 0.0
        peak_sum = 1.0
    end

    for (i, cell_type) in enumerate(show_types)
        cell_idx = F[:JpostDict][cell_type]
        response_normalized = F[:rs][:, cell_idx] / peak_sum

        lines!(axes[i], p.times, response_normalized, 
               label=p.cell_type_name_mapping[cell_type],
               linewidth=1.0,
               color=p.colorDict[cell_type])
        fill_between!(axes[i], p.times, zeros(length(p.times)), 
                      response_normalized, color=(p.colorDict[cell_type],0.05))
    end

    for ax in axes[2:end]
        hideydecorations!(ax, grid=false)
        ax.leftspinevisible = false
    end

    all_axes = [axes...]
    linkyaxes!(all_axes...)

    for ax in all_axes
        hidexdecorations!(ax, ticks=false)
    end
    axes[1].ylabel = "norm. response"
    axes[1].yticklabelsize = 5
    axes[1].ylabelsize = 5
    axes[1].yticks = LinearTicks(3)

    title = p.cell_type_name_mapping[senso_type]
    Label(content_layout[1, 1:length(show_types)], title, fontsize=6, 
          halign=:center, valign=:top, tellheight=false, tellwidth=false)

    colgap!(content_layout, 2)

    if standalone
        return fig
    end
    #}}}
end


########################
function plot_boxcar(; p = dp, 
                       fig = nothing, 
                       content_layout = nothing)

    #{{{
    standalone = isnothing(fig)
    if standalone
        set_my_theme!(fontsize=5, xlabelpadding=2, ylabelpadding=2, xticklabelpad=0, linewidth=0.5)
        fig = Figure(size=cm2px(2, 2))
        content_layout = GridLayout()
        fig[1, 1] = content_layout
    end

    Label(content_layout[1, 3], "sensory stimulation",
          fontsize=6, halign=:center, valign=:top,
          tellwidth = false, tellheight=false)

    axes = [Axis(content_layout[1, i]) for i = 1:4]

    boxcar = zeros(p.tstep)
    start_idx = round(Int, p.stimonset/p.dt) + 1
    end_idx   = round(Int, (p.stimonset + p.stimdurat)/p.dt)
    boxcar[start_idx:end_idx] .= p.stimampli

    lines!(axes[1], p.times, boxcar, color=:gray, linewidth=1.0)
    fill_between!(axes[1], p.times, zeros(length(p.times)), boxcar, color=(:gray,0.05))
    axes[1].ylabel = "sensory input"
    axes[1].xlabel = "time"

    for ax in axes[2:end]
        hideydecorations!(ax, grid=false)
        ax.leftspinevisible = false
        hidexdecorations!(ax, grid=false)
    end

    colgap!(content_layout, 2)

    if standalone
        return fig
    end
    #}}}
end


########################
function plot_heatmap_selectivity(; p = dp, 
                                  show_types=["P1", "P2-poly", "P2-noci", "P2-cold"],
                                  fig = nothing,
                                  content_layout = nothing)
    #{{{
    standalone = isnothing(fig)

    if standalone
        set_my_theme!(fontsize=5, xlabelpadding=2, ylabelpadding=2, xticklabelpad=0, linewidth=0.5)
        fig = Figure(size=cm2px(3.2, 6.5))
        content_layout = GridLayout()
        fig[1, 1] = content_layout
    end

    p_local = deepcopy(p)
    excit_cell_type = setdiff(p_local.excit_cell_type, show_types)

    n_projection  = length(show_types)
    n_sensory     = length(p_local.senso_cell_type)
    n_excitatory  = length(excit_cell_type)
    
    ##### sensory neurons #####
    peak_matrix_sensory = zeros(n_projection, n_sensory)
    
    for (j, senso_type) in enumerate(p_local.senso_cell_type)
        p_iter = deepcopy(p_local)
        for key in keys(p_iter.h0Dict)
            p_iter.h0Dict[key] = 0
        end
        p_iter.h0Dict[senso_type] = 1
        for key in keys(p_iter.I0Dict)
            p_iter.I0Dict[key] = 0.0
        end
        F = runsim(p = p_iter)
        for (i, proj_type) in enumerate(show_types)
            cell_idx = F[:JpostDict][proj_type]
            peak_matrix_sensory[i, j] = F[:steady_rs][cell_idx]
        end
    end
    
    peak_matrix_sensory_norm = peak_matrix_sensory ./ sum(peak_matrix_sensory, dims=1)
    peak_matrix_sensory_norm[isnan.(peak_matrix_sensory_norm)] .= 0

    selectivity_sensory = [calc_selectivity(peak_matrix_sensory_norm[:, j]) for j in 1:n_sensory]
    sorted_indices_sensory = sortperm(selectivity_sensory, rev=true)

    selectivity_sensory_sorted = selectivity_sensory[sorted_indices_sensory]
    peak_matrix_sensory_sorted = peak_matrix_sensory_norm[:, sorted_indices_sensory]
    sensory_labels_sorted = [p_local.cell_type_name_mapping[p_local.senso_cell_type[i]] for i in sorted_indices_sensory]
    
    ##### excitatory neurons #####
    peak_matrix_excitatory = zeros(n_projection, n_excitatory)
    
    for (j, excit_type) in enumerate(excit_cell_type)
        p_iter = deepcopy(p_local)
        for key in keys(p_iter.h0Dict)
            p_iter.h0Dict[key] = 0
        end
        for key in keys(p_iter.I0Dict)
            p_iter.I0Dict[key] = 0.0
        end
        p_iter.I0Dict[excit_type] = 1.0
        F = runsim(p = p_iter)
        for (i, proj_type) in enumerate(show_types)
            cell_idx = F[:JpostDict][proj_type]
            peak_matrix_excitatory[i, j] = F[:steady_rs][cell_idx]
        end
    end
    
    peak_matrix_excitatory_norm = peak_matrix_excitatory ./ sum(peak_matrix_excitatory, dims=1)
    peak_matrix_excitatory_norm[isnan.(peak_matrix_excitatory_norm)] .= 0

    selectivity_excitatory = [calc_selectivity(peak_matrix_excitatory_norm[:, j]) for j in 1:n_excitatory]
    sorted_indices_excitatory = sortperm(selectivity_excitatory, rev=true)

    selectivity_excitatory_sorted = selectivity_excitatory[sorted_indices_excitatory]
    peak_matrix_excitatory_sorted = peak_matrix_excitatory_norm[:, sorted_indices_excitatory]
    excitatory_labels_sorted = [p_local.cell_type_name_mapping[excit_cell_type[i]] for i in sorted_indices_excitatory]
    
    ##### figure #####
    ax1 = Axis(content_layout[1:2, 1])#, title="sensory")
    ax2 = Axis(content_layout[3:5, 1], ylabel="PN subtype\nselectivity")
    ax3 = Axis(content_layout[6:7, 1])#, title="excitatory")
    ax4 = Axis(content_layout[8:10, 1], ylabel="PN subtype\nselectivity")
    ax2.ylabelsize = 6
    ax4.ylabelsize = 6
    ax2.xlabelsize = 6
    ax4.xlabelsize = 6

  
    for a in [ax2, ax4]
        hlines!(a, [0.5, 1], color=(:gray,0.4), linewidth=0.5)
    end
    
    my_colormap = Reverse(:grayC)

    hm1 = heatmap!(ax1, 1:n_sensory, 1:n_projection, 
                   reverse(peak_matrix_sensory_sorted, dims=1)', 
                   colormap=my_colormap, colorrange=(0,1))
    ax1.yreversed = true
    ylims!(ax1, 0.5, n_projection + 0.5)
    ax1.yticks = (1:n_projection, reverse([p_local.cell_type_name_mapping[st] for st in show_types]))
    hidexdecorations!(ax1)
    
    scatter!(ax2, 1:n_sensory, selectivity_sensory_sorted, markersize=3.5, color="#696969")
    lines!(ax2, 1:n_sensory, selectivity_sensory_sorted, linewidth=1, color="#696969")
    ax2.xticks = (1:n_sensory, sensory_labels_sorted)
    ax2.yticks = (0:0.5:1, ["0", "0.5", "1"])
    ax2.xticklabelrotation = π/4
    ylims!(ax2, 0, 1)
    
    hm3 = heatmap!(ax3, 1:n_excitatory, 1:n_projection, 
                   reverse(peak_matrix_excitatory_sorted, dims=1)', 
                   colormap=my_colormap, colorrange=(0,1))
    ax3.yreversed = true
    ylims!(ax3, 0.5, n_projection + 0.5)
    ax3.yticks = (1:n_projection, reverse([p_local.cell_type_name_mapping[st] for st in show_types]))
    hidexdecorations!(ax3)
    
    scatter!(ax4, 1:n_excitatory, selectivity_excitatory_sorted, markersize=4, color="#696969")
    lines!(ax4, 1:n_excitatory, selectivity_excitatory_sorted, linewidth=1, color="#696969")
    ax4.xticks = (1:n_excitatory, excitatory_labels_sorted)
    ax4.yticks = (0:0.5:1, ["0", "0.5", "1"])
    ax4.xticklabelrotation = π/4
    ylims!(ax4, 0, 1)

    linkxaxes!(ax1, ax2)
    linkxaxes!(ax3, ax4)
    xlims!(ax2, 0.5, n_sensory + 0.5)
    xlims!(ax4, 0.5, n_excitatory + 0.5)
    
    ax1.title = "sensory"
    ax1.titlesize = 6
    ax3.title = "excitatory"
    ax3.titlesize = 6

    if standalone
        return fig
    end
    #}}}
end


#########################
function calc_selectivity(x)
    sorted_x = sort(x, rev=true)
    r_max    = sorted_x[1]
    r_second = sorted_x[2]
    if r_max + r_second == 0
        return 0.0
    end
    return (r_max - r_second) / (r_max + r_second)
end


##############################################################
function Plot_matrix(J, preDict, postDict;
                     pre_cell_type_name_mapping=nothing, 
                     post_cell_type_name_mapping=nothing,
                     maxcolor=nothing,
                     mincolor=nothing,
                     mycolormap=Reverse(:balance),
                     show_colorbar=true,
                     show_ylabel=true,
                     fig=nothing,
                     content_layout=nothing,
                     p=dp)

    #{{{
    standalone = isnothing(fig)
    if standalone
        set_my_theme!(fontsize=4, ticksize=1, xticklabelpad=0.5, yticklabelpad=1.0)
        fig = Figure(size=(cm2px(10,4.25)))
        content_layout = GridLayout()
        fig[1,1] = content_layout
    end

    if isnothing(pre_cell_type_name_mapping)
        pre_cell_type_name_mapping = p.cell_type_name_mapping
    end
    if isnothing(post_cell_type_name_mapping)
        post_cell_type_name_mapping = p.cell_type_name_mapping
    end

    N, M = size(J)
    ax = Axis(content_layout[1,1], 
              xticklabelsize=4.5, 
              yticklabelsize=4.5, 
              xlabelsize=6,
              ylabelsize=6,
             )
    ax.aspect = DataAspect() 

    J_hmARG = copy(J)
    
    if isnothing(maxcolor)
        hm = heatmap!(ax, 1:M, 1:N, J_hmARG', colormap=mycolormap)
        barlabel = "effective connection strength"
    else
        actual_mincolor = isnothing(mincolor) ? -maxcolor : mincolor
        hm = heatmap!(ax, 1:M, 1:N, J_hmARG', 
                     colorrange=(actual_mincolor, maxcolor), 
                     colormap=mycolormap)
        barlabel = "effective connection strength"
    end
    
    if show_colorbar
        cbar = Colorbar(content_layout[1,2], hm, label=barlabel)
        tickvals = cbar.ticks[]
        if tickvals isa Makie.Automatic || !(tickvals isa AbstractVector)
            actual_mincolor = isnothing(mincolor) ? -maxcolor : mincolor
            tickvals = range(actual_mincolor, maxcolor, length=5)
        end
        ticklabels = String[]
        for (i, val) in enumerate(tickvals)
            if i == 1 && !isnothing(mincolor)
                push!(ticklabels, "≤ " * string(round(val, digits=2)))
            elseif i == length(tickvals) && !isnothing(maxcolor)
                push!(ticklabels, "≥ " * string(round(val, digits=2)))
            else
                push!(ticklabels, string(round(val, digits=2)))
            end
        end
        cbar.ticks = (collect(tickvals), ticklabels)
    end
    
    ax.yreversed = true
    xlims!(ax, 0.5, M+0.5)
    ylims!(ax, 0.5, N+0.5)
    
    pre_indices  = [preDict[k]  for k in keys(preDict)]
    post_indices = [postDict[k] for k in keys(postDict)]
    
    ax.xticks = (pre_indices, [pre_cell_type_name_mapping[k] for k in keys(preDict)])
    ax.xticklabelrotation = π/4

    if show_ylabel
        ax.yticks = (post_indices, [post_cell_type_name_mapping[k] for k in keys(postDict)])
        ax.ylabel = "postsynaptic"
    else
        hideydecorations!(ax, grid=false)
    end
    Label(content_layout[0, 1], 
          "presynaptic excitatory\ninhibitory", 
          fontsize=6, halign=:center, valign=:top,
          tellwidth=false)
    #ax.title = "presynaptic excitatory\ninhibitory"
    #ax.titlesize = 6
    ax.rightspinevisible = true
    ax.topspinevisible = true
    
    if standalone
        return fig
    else
        return hm
    end
    #}}}
end


##############################################################
function plot_main_figure(; p = dp)
    #{{{
    set_my_theme!(fontsize=5, xlabelpadding=2, ylabelpadding=2, xticklabelpad=0, linewidth=0.5)

    fig = Figure(size=(cm2px(13.5,15)))

    left_frac  = 0.7
    right_frac = 1 - left_frac
    
    left_top_frac  = 0.4
    left_bot_frac  = 1 - left_top_frac
    right_top_frac = 0.5
    right_bot_frac = 1 - right_top_frac
    
    left_layout  = GridLayout()
    right_layout = GridLayout()
    fig[1, 1] = left_layout
    fig[1, 2] = right_layout
    
    #################################
    ############### A ###############
    #################################
    A_layout = GridLayout()
    left_layout[1, 1] = A_layout
    
    #################################
    J_full, W_full, JpreDict_full, JpostDict_full, WpreDict_full, WpostDict_full = set_connectivity(p=p)
    J, JpreDict, JpostDict = reduce_connectivity(J_full, JpreDict_full, JpostDict_full)
    W, WpreDict, WpostDict = reduce_connectivity(W_full, WpreDict_full, WpostDict_full)

    ###### E/I scaling #####
    for (pre_types, post_types, β) in [
        (p.excit_cell_type, p.excit_cell_type, p.β_EE),  # E→E
        (p.excit_cell_type, p.inhib_cell_type, p.β_IE),  # E→I
        (p.inhib_cell_type, p.excit_cell_type, p.β_EI),  # I→E
        (p.inhib_cell_type, p.inhib_cell_type, p.β_II)   # I→I
    ]
        for pre_ct in pre_types, post_ct in post_types
            J[JpostDict[post_ct], JpreDict[pre_ct]] *= β
        end
    end
    ################################

    proj_types = ["P1", "P2-poly", "P2-noci", "P2-cold"]
    JpreDict_filtered = filter(q -> !(q.first in proj_types), JpreDict)
    
    keep_cols = [JpreDict[k] for k in keys(JpreDict_filtered)]
    J_filtered = J[:, keep_cols]
    
    JpreDict_filtered = OrderedDict(k => i for (i,k) in enumerate(keys(JpreDict_filtered)))
    ###############################

    J_sublayout = GridLayout()
    A_layout[1, 1] = J_sublayout
    Plot_matrix(J_filtered, JpreDict_filtered, JpostDict,
                maxcolor=0.5, mincolor=-0.5,
                show_colorbar=false,
                show_ylabel=true,
                fig=fig, content_layout=J_sublayout, p=p)
    
    W_sublayout = GridLayout()
    A_layout[1, 2] = W_sublayout
    hm_ref = Plot_matrix(W, WpreDict, WpostDict,
                         maxcolor=0.5, mincolor=-0.5,
                         show_colorbar=false,
                         show_ylabel=false,
                         fig=fig, content_layout=W_sublayout, p=p)
    
    #~cbar = Colorbar(A_layout[1, 3], hm_ref, 
    #~                label="connection strength", 
    #~                width=8)
    #~tickvals  = collect(range(-0.5, 0.5, length=5))
    #~ticklabels = ["≤ -0.5", "-0.25", "0", "0.25", "≥ 0.5"]
    #~cbar.ticks = (tickvals, ticklabels)
    
    colgap!(A_layout, 0)
    rowgap!(A_layout, 0)

    #################################
    ############### B ###############
    #################################
    B_layout = GridLayout()
    left_layout[2, 1] = B_layout

    boxcar_layout = GridLayout()
    B_layout[1, 1] = boxcar_layout
    plot_boxcar(p=p, fig=fig, content_layout=boxcar_layout)

    senso_types_B = ["cltmr", "mrgd", "chtmr-p",
                     "ccold", "adltmr", "sst", "adhtmr"]

    for (k, senso_type) in enumerate(senso_types_B)
        if k <= 3
            row = k + 1
            col = 1
        else
            row = k - 3
            col = 2
        end
        sub_layout = GridLayout()
        B_layout[row, col] = sub_layout
        plot_stimulations_sidebyside(p = deepcopy(p),
                                     senso_type = senso_type,
                                     fig = fig,
                                     content_layout = sub_layout)
    end
    for col in 1:2
        colsize!(B_layout, col, Relative(0.5))
    end
    for row in 1:4
        rowsize!(B_layout, row, Relative(0.25))
    end


    #################################
    ############### C ###############
    #################################
    C_layout = GridLayout()
    right_layout[1, 1] = C_layout
    plot_heatmap_selectivity(p=p, fig=fig, content_layout=C_layout)


    #################################
    ############### D ###############
    #################################
    ax_D = Axis(right_layout[2, 1])
    
    colsize!(fig.layout, 1, Relative(left_frac))
    colsize!(fig.layout, 2, Relative(right_frac))
    
    rowsize!(left_layout,  1, Relative(left_top_frac))
    rowsize!(left_layout,  2, Relative(left_bot_frac))
    rowsize!(right_layout, 1, Relative(right_top_frac))
    rowsize!(right_layout, 2, Relative(right_bot_frac))
    

    for gl in (left_layout, right_layout, C_layout, A_layout)
        colgap!(gl, 0)
        rowgap!(gl, 0)
    end

    rowgap!(B_layout, 12)


    rowgap!(left_layout, 1, 8)
    rowgap!(right_layout, 1, 8)
    colgap!(fig.layout, 1, 40)
    rowgap!(fig.layout, 0)
   
    return fig
    #}}}
end


#########################
function plot_selectivity_combined(; p = dp,
                                    show_types = ["P1", "P2-poly", "P2-noci", "P2-cold"],
                                    #β_E_set = [0.5, 1.5],
                                    #β_I_set = [2.0, 2.0],
                                    β_E_set = [1.0, 1.0],
                                    β_I_set = [1.0, 4.0],
                                    colors  = ["#d99e73", "#96bfe6"],
                                    markers = [:diamond, :xcross])
    #{{{
    set_my_theme!(fontsize=5, xlabelpadding=2, ylabelpadding=2, xticklabelpad=0, linewidth=0.5)

    n_sensory    = length(p.senso_cell_type)
    n_projection = length(show_types)
    n_models     = 1 + length(β_E_set)

    fig = Figure(size=cm2px(4.5,4))

    ax1 = Axis(fig[1, 1], ylabel = "norm.\nresponse")
    ax2 = Axis(fig[2:3, 1], ylabel = "PN subtype\nselectivity")

    hlines!(ax2, [0.5, 1.0], color=(:gray, 0.4), linewidth=0.5)

    ########## shared computation ##########
    p_ref = deepcopy(p)
    peak_matrix_ref = zeros(n_projection, n_sensory)
    for (j, senso_type) in enumerate(p_ref.senso_cell_type)
        p_iter = deepcopy(p_ref)

        for key in keys(p_iter.h0Dict); p_iter.h0Dict[key] = 0; end
        p_iter.h0Dict[senso_type] = 1

        for key in keys(p_iter.I0Dict); p_iter.I0Dict[key] = 0.0; end
        F = runsim(p = p_iter)

        for (i, proj_type) in enumerate(show_types)
            cell_idx = F[:JpostDict][proj_type]
            peak_matrix_ref[i, j] = F[:steady_rs][cell_idx]
        end
    end
    peak_matrix_ref_norm = peak_matrix_ref ./ sum(peak_matrix_ref, dims=1)
    peak_matrix_ref_norm[isnan.(peak_matrix_ref_norm)] .= 0

    selectivity_ref        = [calc_selectivity(peak_matrix_ref_norm[:, j]) for j in 1:n_sensory]
    sorted_indices         = sortperm(selectivity_ref, rev=true)
    selectivity_ref_sorted = selectivity_ref[sorted_indices]
    sensory_labels_sorted  = [p_ref.cell_type_name_mapping[p_ref.senso_cell_type[i]]
                               for i in sorted_indices]

    test_norm = Vector{Matrix{Float64}}(undef, length(β_E_set))
    for (k, (βE, βI)) in enumerate(zip(β_E_set, β_I_set))
        p_test = deepcopy(p)
        p_test.β_EE = βE
        p_test.β_IE = βE
        p_test.β_EI = βI
        p_test.β_II = βI
        peak_matrix = zeros(n_projection, n_sensory)

        for (j, senso_type) in enumerate(p_test.senso_cell_type)
            p_iter = deepcopy(p_test)

            for key in keys(p_iter.h0Dict); p_iter.h0Dict[key] = 0; end
            p_iter.h0Dict[senso_type] = 1

            for key in keys(p_iter.I0Dict); p_iter.I0Dict[key] = 0.0; end
            F = runsim(p = p_iter)

            for (i, proj_type) in enumerate(show_types)
                cell_idx = F[:JpostDict][proj_type]
                peak_matrix[i, j] = F[:steady_rs][cell_idx]
            end
        end
        peak_matrix_norm = peak_matrix ./ sum(peak_matrix, dims=1)
        peak_matrix_norm[isnan.(peak_matrix_norm)] .= 0
        test_norm[k] = peak_matrix_norm
    end

    all_norm = Vector{Matrix{Float64}}(undef, n_models)
    all_norm[1] = test_norm[1]
    all_norm[2] = peak_matrix_ref_norm 
    all_norm[3] = test_norm[2]

    ref_label = rich("ref: ",
                     rich("β", font="Helvetica Neue Italic"), subscript("E"),
                     " = ", string(p.β_EE), ", ",
                     rich("β", font="Helvetica Neue Italic"), subscript("I"),
                     " = ", string(p.β_EI))
    model_labels = Any[ref_label]
    for (βE, βI) in zip(β_E_set, β_I_set)
        lbl = rich(rich("β", font="Helvetica Neue Italic"), subscript("E"),
                   " = ", string(βE), ", ",
                   rich("β", font="Helvetica Neue Italic"), subscript("I"),
                   " = ", string(βI))
        push!(model_labels, lbl)
    end

    ########## stacked bar ##########
    x_vals     = Int[]
    y_vals     = Float64[]
    dodge_vals = Int[]
    stack_vals = Int[]
    color_vals = Any[]

    for (col, j) in enumerate(sorted_indices)
        for model_k in 1:n_models
            for (pn_idx, pn_type) in enumerate(show_types)
                push!(x_vals,     col)
                push!(y_vals,     all_norm[model_k][pn_idx, j])
                push!(dodge_vals, model_k)
                push!(stack_vals, pn_idx)
                push!(color_vals, p.colorDict[pn_type])
            end
        end
    end

    barplot!(ax1, x_vals, y_vals,
             dodge       = dodge_vals,
             stack       = stack_vals,
             color       = color_vals,
             strokewidth = 0.0,
             gap         = 0.1,
             dodge_gap   = 0.04,
             width       = 0.7)

    ax1.yticks = (0:0.5:1, ["0", "0.5", "1"])
    ax1.xticks = 1:n_sensory
    ylims!(ax1, 0, 1)
    xlims!(ax1, 0.5, n_sensory + 0.5)
    hidexdecorations!(ax1, ticks=false)

    ########## selectivity overlay ##########
    scatter!(ax2, 1:n_sensory, selectivity_ref_sorted,
             markersize=4, color="#696969")
    lines!(ax2, 1:n_sensory, selectivity_ref_sorted,
           linewidth=1, color="#696969", label=ref_label)

    for (k, (βE, βI)) in enumerate(zip(β_E_set, β_I_set))
        norm_idx = k == 1 ? 1 : 3
        selectivity        = [calc_selectivity(all_norm[norm_idx][:, j]) for j in 1:n_sensory]
        selectivity_sorted = selectivity[sorted_indices]
        scatter!(ax2, 1:n_sensory, selectivity_sorted,
                 markersize=4, color=colors[k], marker=markers[k])
        lines!(ax2, 1:n_sensory, selectivity_sorted,
               linewidth=1, color=(colors[k], 0.8), label=model_labels[k+1])
    end

    ax2.xticks = (1:n_sensory, sensory_labels_sorted)
    ax2.xticklabelrotation = π/6
    ax2.yticks = (0:0.5:1, ["0", "0.5", "1"])
    ylims!(ax2, 0, 1)
    xlims!(ax2, 0.5, n_sensory + 0.5)

    linkxaxes!(ax1, ax2)

    rowgap!(fig.layout, 0)

    return fig
    #}}}
end


#########################
function compute_spearman_heatmap(; p = dp,
                                   show_types = ["P1", "P2-poly", "P2-noci", "P2-cold"],
                                   β_E_range = 0.1:0.1:2.0,
                                   β_I_range = 0.1:0.1:5.0)
    #{{{
    nβE = length(β_E_range)
    nβI = length(β_I_range)

    spearman_sensory    = zeros(nβE, nβI)
    instability_sensory = falses(nβE, nβI)

    ########## reference sensory stimulation ##########
    ref_peaks_sensory = zeros(length(show_types), length(p.senso_cell_type))
    for (j, senso_type) in enumerate(p.senso_cell_type)
        p_iter = deepcopy(p)
        for key in keys(p_iter.h0Dict); p_iter.h0Dict[key] = 0; end
        p_iter.h0Dict[senso_type] = 1
        for key in keys(p_iter.I0Dict); p_iter.I0Dict[key] = 0.0; end
        F = runsim(p = p_iter)
        for (i, proj_type) in enumerate(show_types)
            cell_idx = F[:JpostDict][proj_type]
            ref_peaks_sensory[i, j] = F[:steady_rs][cell_idx]
        end
    end

    ########## Scan over β_E and β_I ##########
    prog = Progress(nβE * nβI, 1, "computing spearman sensory (β_E vs β_I)... ")
    for (iE, βE) in enumerate(β_E_range)
        for (iI, βI) in enumerate(β_I_range)

            p_test = deepcopy(p)
            p_test.β_EE = βE
            p_test.β_IE = βE
            p_test.β_EI = βI
            p_test.β_II = βI

            test_peaks_sensory = zeros(length(show_types), length(p.senso_cell_type))
            for (j, senso_type) in enumerate(p.senso_cell_type)
                p_iter = deepcopy(p_test)

                for key in keys(p_iter.h0Dict); p_iter.h0Dict[key] = 0; end
                p_iter.h0Dict[senso_type] = 1

                for key in keys(p_iter.I0Dict); p_iter.I0Dict[key] = 0.0; end
                F = runsim(p = p_iter)
                if F[:is_unstable]
                    instability_sensory[iE, iI] = true
                end
                for (i, proj_type) in enumerate(show_types)
                    cell_idx = F[:JpostDict][proj_type]
                    test_peaks_sensory[i, j] = F[:steady_rs][cell_idx]
                end
            end

            spearman_sensory[iE, iI] = corspearman(vec(ref_peaks_sensory),
                                                    vec(test_peaks_sensory))
            next!(prog)
        end
    end

    spearman_sensory[instability_sensory] .= NaN

    return Dict(
        :spearman_sensory    => spearman_sensory,
        :instability_sensory => instability_sensory,
        :β_E_vals            => collect(β_E_range),
        :β_I_vals            => collect(β_I_range),
        :β_E_ref             => p.β_EE,
        :β_I_ref             => p.β_EI,
    )
    #}}}
end


#########################
function plot_spearman_heatmap(results_dict; p = dp, 
                               alongI_beta_E_set = [1.0, 1.0],
                               alongI_beta_I_set = [1.0, 4.0],
                               alongE_beta_E_set = [0.5, 1.5],
                               alongE_beta_I_set = [2.0, 2.0])
    #{{{
    set_my_theme!(fontsize=5, xlabelpadding=2, ylabelpadding=2, xticklabelpad=0, linewidth=0.5)
    fig = Figure(size=cm2px(4, 4))

    ax1 = Axis(fig[1, 1], 
               aspect = AxisAspect(1),
               title  = "sensory neurons",    
               xlabel = rich("excitatory strength (",
                             rich("β", font="Helvetica Neue Italic"),
                             subscript("E"),
                             ")"),
               ylabel = rich("inhibitory strength (",
                             rich("β", font="Helvetica Neue Italic"),
                             subscript("I"),
                             ")"))

    β_E_vals            = results_dict[:β_E_vals]
    β_I_vals            = results_dict[:β_I_vals]
    spearman_sensory    = results_dict[:spearman_sensory]
    instability_sensory = results_dict[:instability_sensory]
    β_E_ref             = results_dict[:β_E_ref]
    β_I_ref             = results_dict[:β_I_ref]

    heatmap!(ax1, β_E_vals, β_I_vals, Float64.(instability_sensory),
             colormap=[:transparent, :orange], colorrange=(0,1))

    hm1 = heatmap!(ax1, β_E_vals, β_I_vals, spearman_sensory,
                   colormap=Reverse(:tempo), colorrange=(0,1), nan_color=(:white,0))


    vlines!(ax1, [β_E_ref], color=(:gray, 0.6), linewidth=0.5, linestyle=:dash)
    hlines!(ax1, [β_I_ref], color=(:gray, 0.6), linewidth=0.5, linestyle=:dash)
    
    scatter!(ax1, [β_E_ref], [β_I_ref], color=:red, markersize=6, marker=:star5)

    #### along I
    scatter!(ax1, [alongI_beta_E_set[1]], [alongI_beta_I_set[1]], 
             color=:black, markersize=5, marker=:diamond)
    scatter!(ax1, [alongI_beta_E_set[2]], [alongI_beta_I_set[2]], 
             color=:black, markersize=5, marker=:xcross)
    
    #### along E
    scatter!(ax1, [alongE_beta_E_set[1]], [alongE_beta_I_set[1]], 
             color=:transparent, markersize=5, marker=:diamond,
             strokecolor=:black, strokewidth=0.5)
    scatter!(ax1, [alongE_beta_E_set[2]], [alongE_beta_I_set[2]], 
            color=:transparent, markersize=5, marker=:xcross,
            strokecolor=:black, strokewidth=0.5)

    cb_label = "Spearman correlation"
    Colorbar(fig[1, 2], hm1, label=cb_label, width=5)

    colgap!(fig.layout, 3)

    return fig
    #}}}
end


#########################
function set_my_theme!(;fontsize=5, 
                       linewidth=0.5,
                       tickwidth=0.3, 
                       ticksize=1, 
                       xlabelpadding=5, 
                       ylabelpadding=5, 
                       xticklabelpad=0.5, 
                       yticklabelpad=1.0, 
                       figure_padding=5,
                      )

    #{{{
    set_theme!(Theme(fontsize = fontsize, 
                     font = "Helvetica Neue", 
                     figure_padding = figure_padding, 
                     Lines = (linewidth = linewidth,), 
                     Axis = (xgridvisible = false,
                             ygridvisible = false,
                             rightspinevisible = false,
                             topspinevisible = false,
                             xtickwidth = tickwidth, 
                             ytickwidth = tickwidth, 
                             spinewidth = tickwidth, 
                             xticksize = ticksize, 
                             yticksize = ticksize, 
                             titlefont = "Helvetica Neue",
                             xlabelpadding = xlabelpadding, 
                             ylabelpadding = ylabelpadding,
                             xticklabelpad = xticklabelpad, 
                             yticklabelpad = yticklabelpad
                            ),
                     Colorbar = (tickwidth = tickwidth,
                             spinewidth = tickwidth,
                             ticksize = ticksize,
                             labelsize = fontsize,
                             ticklabelsize = fontsize,
                            )
                    )
              )
    #}}}
end


#########################
function cm2px(w, h; dpi = 96)
    w_px = floor(Int, w * dpi / 2.54) / 1.333327787
    h_px = floor(Int, h * dpi / 2.54) / 1.333327787
    return (w_px, h_px)
end


end

