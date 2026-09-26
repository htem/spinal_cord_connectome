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

    ##### optimization #####
    n_iters::Int    = 1000
    n_restarts::Int = 100
    ########################

    T::Float64    = 85
    dt::Float64   = 0.01
    tstep::Int    = floor(Int, T/dt)
    times::Array  = 0:dt:T-dt

    stimonset::Float64 = 1
    stimdurat::Float64 = 84
    stimampli::Float64 = 1

    normalize_conn::Bool = true

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

    # all active as default; h_input will scale each one
    h0Dict::Dict{String,Bool} = Dict(
                       "sst"     => 1,
                       "cltmr"   => 1,
                       "adltmr"  => 1,
                       "mrgd"    => 1,
                       "adhtmr"  => 1,
                       "chtmr-p" => 1,
                       "ccold"   => 1,
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
                                         "R"             => "R",
                                         "module-cltmr"  => "C-LTMR mod",
                                         "module-np"     => "C-HTMR/Heat mod",
                                         "module-adltmr" => "Aδ-LTMR mod",
                                         "module-noci"   => "Aδ-HTMR mod",
                                         "V1"            => "V1",
                                         "V2"            => "V2",
                                         "U"             => "unknown",
                                         "P1"            => "P1",
                                         "P2-poly"       => "P2-poly",
                                         "P2-noci"       => "P2-noci",
                                         "P2-cold"       => "P2-cold",
                                         "islet-cltmr"   => rich("islet", superscript("C-LTMR")),
                                         "islet-mrgd"    => rich("islet", superscript("C-HTMR/Heat")),
                                         "islet-adltmr"  => rich("islet", superscript("Aδ-LTMR")),
                                         "islet-adhtmr"  => rich("PSI-IhN", superscript("Aδ-HTMR")),
                                         "islet-V"       => "inhib vertical",
                                         "cltmr"         => "C-LTMR",
                                         "mrgd"          => rich("C-HTMR/Heat (MrgprD/A3/B4", superscript("+"), ")"),
                                         "sst"           => rich("C-HTMR/Heat (SST", superscript("+"), ")"),
                                         "chtmr-p"       => "Peptidergic C-Nociceptor",
                                         "adltmr"        => "Aδ-LTMR",
                                         "adhtmr"        => "Aδ-HTMR",
                                         "ccold"         => "C-Cold",
                                        )

    colorDict = Dict(
                     "R"            => "#3F5482",
                     "module-cltmr" => "#7B2479",
                     "module-np"    => "#00a168",
                     "module-adltmr"=> "#D67300",
                     "module-noci"  => "#008AA1",
                     "V1"           => "#718600",
                     "V2"           => "#730F1F",
                     "P1"           => "#4D4182",
                     "P2-poly"      => "#B4B753",
                     "P2-noci"      => "#BB2D3B",
                     "P2-cold"      => "#2E6CB2",
                     "islet-cltmr"  => "#0864c3",
                     "islet-mrgd"   => "#c36708",
                     "islet-adltmr" => "#c30864",
                     "islet-adhtmr" => "#08c367",
                     "islet-V"      => "#c2c308",
                    )
end

dp = Params()


##############################################################
function read_data(; p = dp)
    #{{{
    J  = readdlm("data/J_recurrent.txt")
    W  = readdlm("data/J_feedforward.txt")
    WpreDict  = load("data/preDict_feedforward.jld2", "data")
    WpostDict = load("data/postDict_feedforward.jld2", "data")
    JpreDict  = load("data/preDict_recurrent.jld2", "data")
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
    index_map     = Dict(old_idx => new_idx for (new_idx, old_idx) in enumerate(kept_indices))
    remapped_dict = OrderedDict()
    for (key, old_range) in dict
        old_indices = collect(old_range)
        new_indices = [index_map[idx] for idx in old_indices]
        if !isempty(new_indices)
            remapped_dict[key] = minimum(new_indices):maximum(new_indices)
        end
    end
    return remapped_dict
    #}}}
end


##############################################################
function set_connectivity(; p = dp)
    #{{{
    J0, W0, JpreDict, JpostDict, WpreDict, WpostDict = read_data(p = p)

    ############## filter to only include specified cell types #############
    valid_neuron_types  = vcat(p.excit_cell_type, p.inhib_cell_type)
    valid_sensory_types = p.senso_cell_type

    JpreDict_filtered  = filter(q -> q.first in valid_neuron_types,  JpreDict)
    JpostDict_filtered = filter(q -> q.first in valid_neuron_types,  JpostDict)
    WpreDict_filtered  = filter(q -> q.first in valid_sensory_types, WpreDict)
    WpostDict_filtered = filter(q -> q.first in valid_neuron_types,  WpostDict)

    neuron_indices_pre  = sort(vcat([collect(inds) for inds in values(JpreDict_filtered)]...))
    neuron_indices_post = sort(vcat([collect(inds) for inds in values(JpostDict_filtered)]...))
    sensory_indices     = sort(vcat([collect(inds) for inds in values(WpreDict_filtered)]...))

    J0 = J0[neuron_indices_post, neuron_indices_pre]
    W0 = W0[neuron_indices_post, sensory_indices]

    JpreDict  = remap_dict_ranges(JpreDict_filtered,  neuron_indices_pre)
    JpostDict = remap_dict_ranges(JpostDict_filtered, neuron_indices_post)
    WpreDict  = remap_dict_ranges(WpreDict_filtered,  sensory_indices)
    WpostDict = remap_dict_ranges(WpostDict_filtered, neuron_indices_post)

    ############# inhibitory (-) #############
    for pre_ct in p.inhib_cell_type
        J0[:, JpreDict[pre_ct]] .*= -1
    end

    ############# normalizing connectivity #############
    if p.normalize_conn
        indegree = dropdims(sum(abs.(J0), dims=2), dims=2) + dropdims(sum(abs.(W0), dims=2), dims=2)
        if any(iszero, indegree)
            indegree[indegree .== 0] .= 1
        end
        J = J0 ./ indegree
        W = W0 ./ indegree
    else
        J = J0
        W = W0
    end
    return J, W, JpreDict, JpostDict, WpreDict, WpostDict
    #}}}
end


##############################################################
function reduce_connectivity(J, preDict, postDict)
    #{{{
    J_bar = zeros(length(postDict), length(preDict))
    for (j, jtype) in enumerate(keys(preDict))
        for (i, itype) in enumerate(keys(postDict))
            J_bar[i,j] = sum(J[postDict[itype], preDict[jtype]]) / length(postDict[itype])
        end
    end
    preDict_bar  = OrderedDict{String,Int}(k => i for (i,k) in enumerate(keys(preDict)))
    postDict_bar = OrderedDict{String,Int}(k => i for (i,k) in enumerate(keys(postDict)))
    return J_bar, preDict_bar, postDict_bar
    #}}}
end


##############################################################
function project_h(x::AbstractVector)
    #{{{
    h   = max.(x, 0.0)       # non-negativity constraint
    nrm = norm(h)
    if nrm == 0
        h .= 1 / sqrt(length(h))
    else
        h ./= nrm
    end
    return h
    #}}}
end


##############################################################
function runsim_optim(h_input::Vector{Float64};
                      p         = dp,
                      target_cell = "P1",
                      J         = nothing,
                      W         = nothing,
                      JpreDict  = nothing,
                      JpostDict = nothing,
                      WpreDict  = nothing)
    #{{{

    ϕ = relu

    N = size(J, 1)
    M = size(W, 2)

    neuron_types = collect(keys(JpreDict))

    excit_indices = [JpreDict[ct] for ct in p.excit_cell_type if ct in neuron_types]
    inhib_indices = [JpreDict[ct] for ct in p.inhib_cell_type if ct in neuron_types]

    #### scale sensory input, h0[idx] = h_input[i] * stimampli
    h0 = zeros(Float64, M)
    for (i, st) in enumerate(p.senso_cell_type)
        idx = WpreDict[st]
        h0[idx] = h_input[i] * p.stimampli
    end

    start_idx = round(Int, p.stimonset / p.dt) + 1
    end_idx = round(Int, (p.stimonset + p.stimdurat) / p.dt)
    ss_idx = min(end_idx, p.tstep)

    Wh0 = W * h0
    zeros_N = zeros(Float64, N)

    I0 = zeros(Float64, N)
    for cell_type in neuron_types
        I0[JpreDict[cell_type]] = p.I0Dict[cell_type]
    end

    τ = zeros(Float64, N)
    for cell_type in neuron_types
        τ[JpreDict[cell_type]] = p.τDict[cell_type]
    end
    dtoverτ = p.dt ./ τ

    JfromE = J[:, excit_indices]
    JfromI = J[:, inhib_indices]

    x = zeros(Float64, N)
    r = zeros(Float64, N)

    steady_r_target = 0.0

    for i = 1:p.tstep
        ηE = JfromE * r[excit_indices]
        ηI = JfromI * r[inhib_indices]
        ηffw = (start_idx <= i <= end_idx) ? Wh0 : zeros_N
        total_input = ηE + ηI + ηffw + I0
        r = r + dtoverτ .* (-r + ϕ.(total_input))
        if i == ss_idx
            steady_r_target = r[JpostDict[target_cell]]
        end
    end

    return -steady_r_target   # negative for minimization
    #}}}
end


##############################################################
function optimize_h_neldermead(; target_cell = "P1",
                                 seed = 0,
                                 p = dp)
    #{{{
    Random.seed!(seed)

    M = length(p.senso_cell_type)

    ##### pre-compute connectivity once #####
    J_full, W_full, JpreDict_full, JpostDict_full, WpreDict_full, WpostDict_full = set_connectivity(p = p)
    J, JpreDict, JpostDict = reduce_connectivity(J_full, JpreDict_full, JpostDict_full)
    W, WpreDict, WpostDict = reduce_connectivity(W_full, WpreDict_full, WpostDict_full)

    for (pre_types, post_types, β) in [
        (p.excit_cell_type, p.excit_cell_type, p.β_EE),
        (p.excit_cell_type, p.inhib_cell_type, p.β_IE),
        (p.inhib_cell_type, p.excit_cell_type, p.β_EI),
        (p.inhib_cell_type, p.inhib_cell_type, p.β_II)
    ]
        for pre_ct in pre_types, post_ct in post_types
            J[JpostDict[post_ct], JpreDict[pre_ct]] *= β
        end
    end
    #########################################

    objective(x) = runsim_optim(project_h(x);
                                p           = p,
                                target_cell = target_cell,
                                J           = J,
                                W           = W,
                                JpreDict    = JpreDict,
                                JpostDict   = JpostDict,
                                WpreDict    = WpreDict)

    best_val = -Inf
    best_h   = nothing
    best_res = nothing

    for k in 1:p.n_restarts
        x0 = project_h(rand(M)) 
        res = optimize(objective, x0, NelderMead(),
                       Optim.Options(iterations = p.n_iters))

        h = project_h(Optim.minimizer(res))
        val = -Optim.minimum(res)

        println("restart $k/$(p.n_restarts): response = ", round(val, digits=4),
                " | iters = ", Optim.iterations(res))

        if val > best_val
            best_val = val
            best_h = h
            best_res = res
        end
    end

    input_keys = collect(keys(p.h0Dict))
    h_dict = OrderedDict(input_keys[i] => best_h[i] for i in 1:M)

    println("\n" * "="^50)
    println("Optimal input pattern for $target_cell:")
    println("="^50)
    for (key, val) in h_dict
        println("  ", rpad(key, 12), ": ", round(val, digits=4))
    end
    println("\nMaximum steady-state response: ", round(best_val, digits=4))

    return best_h, h_dict, best_val, best_res
    #}}}
end


##############################################################
function compute_optimal_h_multiple(; p = dp,
                                      seed = 0,
                                      measure_types = ["P1", "P2-poly", "P2-noci", "P2-cold"])
    #{{{
    M = length(p.senso_cell_type)
    n = length(measure_types)
    H = zeros(M, n)

    for (i, cell_type) in enumerate(measure_types)
        println("\n" * "="^50)
        println("optimizing for: $cell_type ($i/$n)")
        println("="^50)
        best_h, _, best_val, _ = optimize_h_neldermead(
            target_cell = cell_type,
            p = p,
            seed = seed,
        )
        H[:, i] = best_h
        println("best response: $(round(best_val, digits=4))")
    end

    return Dict(:H => H, :measure_types => measure_types)
    #}}}
end


# sensory_order = ["cltmr", "mrgd", "sst", "chtmr-p", "adltmr", "adhtmr", "ccold"]
# ##############################################################
function plot_optimal(results; p = dp,
                      sensory_order = nothing)
    #{{{
    set_my_theme!(fontsize=5)

    fig = Figure(size=cm2px(5, 5))
    ax  = Axis(fig[1,1])
    ax.aspect = DataAspect()

    H = results[:H]
    measure_types = results[:measure_types]

    M = size(H, 1)
    n = length(measure_types)

    sensory_types = p.senso_cell_type

    if !isnothing(sensory_order)
        reorder_indices = [findfirst(==(st), sensory_types) for st in sensory_order
                           if st in sensory_types]
        H_reordered = H[reorder_indices, :]
        sensory_labels = [p.cell_type_name_mapping[st] for st in sensory_order
                          if st in sensory_types]
        M_plot = length(reorder_indices)
    else
        preferred_pn = [argmax(H[i, :]) for i in 1:M]
        reorder_indices = sortperm(preferred_pn)
        H_reordered = H[reorder_indices, :]
        sensory_labels = [p.cell_type_name_mapping[sensory_types[i]]
                           for i in reorder_indices]
        M_plot = M
    end

    hm = heatmap!(ax, 1:n, 1:M_plot, H_reordered',
                  colormap   = :Oranges,
                  colorrange = (0, 1))

    ax.yreversed = true
    ax.ylabel = "sensory neurons"
    xlims!(ax, 0.5, n + 0.5)
    ylims!(ax, 0.5, M_plot + 0.5)

    ax.xticks = (1:n, [p.cell_type_name_mapping[k] for k in measure_types])
    ax.xticklabelrotation = π/4
    ax.yticks = (1:M_plot, sensory_labels)

    for i in 1:(n-1)
        vlines!(ax, i + 0.5, color=:black, linewidth=0.5)
    end
    ax.rightspinevisible = true
    ax.topspinevisible = true

    return fig
    #}}}
end


#########################
function set_my_theme!(; fontsize       = 5,
                         linewidth      = 0.5,
                         tickwidth      = 0.3,
                         ticksize       = 1,
                         xlabelpadding  = 5,
                         ylabelpadding  = 5,
                         xticklabelpad  = 0.5,
                         yticklabelpad  = 1.0,
                         figure_padding = 5)
    #{{{
    set_theme!(Theme(
        fontsize        = fontsize,
        font            = "Helvetica Neue",
        figure_padding  = figure_padding,
        Lines           = (linewidth = linewidth,),
        Axis            = (xgridvisible      = false,
                           ygridvisible      = false,
                           rightspinevisible = false,
                           topspinevisible   = false,
                           xtickwidth        = tickwidth,
                           ytickwidth        = tickwidth,
                           spinewidth        = tickwidth,
                           xticksize         = ticksize,
                           yticksize         = ticksize,
                           titlefont         = "Helvetica Neue",
                           xlabelpadding     = xlabelpadding,
                           ylabelpadding     = ylabelpadding,
                           xticklabelpad     = xticklabelpad,
                           yticklabelpad     = yticklabelpad),
        Colorbar        = (tickwidth      = tickwidth,
                           spinewidth     = tickwidth,
                           ticksize       = ticksize,
                           labelsize      = fontsize,
                           ticklabelsize  = fontsize),
    ))
    #}}}
end


#########################
function cm2px(w, h; dpi = 96)
    w_px = floor(Int, w * dpi / 2.54) / 1.333327787
    h_px = floor(Int, h * dpi / 2.54) / 1.333327787
    return (w_px, h_px)
end


end

