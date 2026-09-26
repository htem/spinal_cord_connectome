module v

using JLD2
using Optim
using Printf
using Random
using GLMakie
using StatsBase
using Parameters
using LinearAlgebra
using ProgressMeter
using DelimitedFiles
using OrderedCollections

@with_kw mutable struct Params

    min_linewidth::Float64 = 0.5
    max_linewidth::Float64 = 3.0
    weight_thresh::Float64 = 0.065

    normalize_conn::Bool = true

    β_EE::Float64 = 1
    β_IE::Float64 = 1

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
                              "P2-cold"
                             ]

    #####
    proje_cell_type::Array = [
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

    #####
    senso_pos = OrderedDict(
                            "sst"     => (1.0, 0.0),
                            "cltmr"   => (4.0, 0.0),
                            "adltmr"  => (7.0, 0.0),
                            "mrgd"    => (10.0, 0.0),
                            "adhtmr"  => (13.0, 0.0),
                            "chtmr-p" => (16.0, 0.0),
                            "ccold"   => (19.0, 0.0),
                           )

    #####
    excit_pos = OrderedDict(
                            "R"             => (2.0, 1.5),
                            "module-cltmr"  => (6.0, 1.8),
                            "module-adltmr" => (10.0, 1.6),
                            "module-np"     => (14.0, 1.7),
                            "module-noci"   => (18.0, 1.5),
                            "V1"            => (2.0, 2.0),
                            "V2"            => (6.0, 2.0),
                           )

    #####
    proje_pos = OrderedDict(
                            "P1"      => (2.5, 3.0),
                            "P2-poly" => (7.5, 3.0),
                            "P2-noci" => (12.5, 3.0),
                            "P2-cold" => (17.5, 3.0),
                           )
    #####
    inhib_pos = OrderedDict(
                            "islet-cltmr"  => (6.0, 1.0),
                            "islet-adltmr" => (10.0, 1.0),
                            "islet-mrgd"   => (14.0, 1.0),
                            "islet-adhtmr" => (18, 1.0),
                           )


    cell_type_name_mapping = OrderedDict(
                                         #{{{
                                         "R"    => "R",
                                         "module-cltmr" => "C-LTMR mod",
                                         "module-np" => "C-HTMR/Heat mod",
                                         "module-adltmr" => "Aδ-LTMR mod",
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
    #}}}
    return J, W, JpreDict, JpostDict, WpreDict, WpostDict
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
function plot_circuit(; p = dp, source_node = nothing)
    #{{{
    ################
    J_full, W_full, JpreDict_full, JpostDict_full, WpreDict_full, WpostDict_full = set_connectivity(p = p)
    J, JpreDict, JpostDict = reduce_connectivity(J_full, JpreDict_full, JpostDict_full)
    W, WpreDict, WpostDict = reduce_connectivity(W_full, WpreDict_full, WpostDict_full)

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
    J = clamp.(abs.(J), 0.0, 0.5)   # since signs are flipped
    W = clamp.(abs.(W), 0.0, 0.5)   # since signs are flipped
    ################

    fig = Figure(size=(cm2px(2* 3.2885, 2* 4.1482)))
    ax = Axis(fig[1, 1], xlabel = "", ylabel = "")
    
    node_size = 2
    
    positions = Point2f[]
    node_labels = String[]
    node_colors = []
    node_types = String[]
    
    #### sensory nodes
    for (st, pos) in p.senso_pos
        push!(positions, Point2f(pos[1], pos[2]))
        push!(node_labels, st)
        push!(node_colors, :orange)
        push!(node_types, "sensory")
    end
    
    #### excitatory nodes
    for (et, pos) in p.excit_pos
        push!(positions, Point2f(pos[1], pos[2]))
        push!(node_labels, et)
        push!(node_colors, haskey(p.colorDict, et) ? p.colorDict[et] : :steelblue)
        push!(node_types, "excitatory")
    end
    
    #### inhibitory nodes
    for (it, pos) in p.inhib_pos
        push!(positions, Point2f(pos[1], pos[2]))
        push!(node_labels, it)
        push!(node_colors, haskey(p.colorDict, it) ? p.colorDict[it] : :red)
        push!(node_types, "inhibitory")
    end
    
    #### projection nodes
    for (pt, pos) in p.proje_pos
        push!(positions, Point2f(pos[1], pos[2]))
        push!(node_labels, pt)
        push!(node_colors, haskey(p.colorDict, pt) ? p.colorDict[pt] : :purple)
        push!(node_types, "projection")
    end
    
    N = length(positions)
    
    name_to_idx = Dict(node_labels[i] => i for i in 1:N)
    
    #### build edge list
    edges = Tuple{Int, Int, Float64}[]
    
    #### add feedforward edges
    for (post_type, post_range) in WpostDict
        if post_type in node_labels
            post_idx = name_to_idx[post_type]

            for (pre_type, pre_range) in WpreDict
                if pre_type in node_labels
                    pre_idx = name_to_idx[pre_type]
                    weight = W[post_range, pre_range]
                    if abs(weight) > p.weight_thresh
                        push!(edges, (pre_idx, post_idx, abs(weight)))
                    end
                end
            end

        end
    end
    
    #### add recurrent edges
    for (post_type, post_range) in JpostDict
        if post_type in node_labels
        post_idx = name_to_idx[post_type]

        for (pre_type, pre_range) in JpreDict
            if pre_type in node_labels
            pre_idx = name_to_idx[pre_type]
            weight = J[post_range, pre_range]
            if abs(weight) > p.weight_thresh
                push!(edges, (pre_idx, post_idx, abs(weight)))
            end
        end

        end
    end
    end
    
    weights = [w for (_, _, w) in edges]

    for (src, dst, w) in edges
        if w < p.weight_thresh
            continue
        end
        if !isnothing(source_node)
            if node_labels[src] != source_node
                continue
            end
        end
        p_src = positions[src]
        p_dst = positions[dst]
        lw = weight_to_linewidth(w, weights)
        if src != dst
            v = p_dst - p_src
            edge_color = node_types[src] == "inhibitory" ? :red : (:black, 0.3)
            
            arrows!(ax,
                [p_src[1]], [p_src[2]],
                [v[1]], [v[2]],
                linewidth = lw,
                arrowsize = 10,
                linecolor = edge_color,
                arrowcolor = edge_color,
            )
        else
            center = p_src
            r = 0.15
            θs = range(0.3π, 1.7π, length=40)
            
            loopx = center[1] .+ r .* cos.(θs)
            loopy = center[2] .+ r .* sin.(θs)
            
            edge_color = node_types[src] == "inhibitory" ? :red : (:black, 0.3)
            lines!(ax, loopx, loopy; color = edge_color, linewidth = lw)
            
            p_before = Point2f(loopx[end-1], loopy[end-1])
            p_end = Point2f(loopx[end], loopy[end])
            v_arrow = p_end - p_before
            arrows!(ax,
                [p_before[1]], [p_before[2]],
                [v_arrow[1]], [v_arrow[2]],
                linewidth = lw,
                arrowsize = 10,
                linecolor = edge_color,
                arrowcolor = edge_color,
            )
        end
    end
    
    xs = [p[1] for p in positions]
    ys = [p[2] for p in positions]
    
    scatter!(ax, xs, ys;
        color = node_colors,
        strokecolor = :black,
        strokewidth = 0.0,
        markersize = node_size,
    )
    
    for i in 1:N
        label_text = get(p.cell_type_name_mapping, node_labels[i], node_labels[i])
        text!(ax, positions[i] .+ Point2f(0.1, 0);
              text = label_text,
              align = (:left, :center),
              fontsize = 3)
    end
    
    hidexdecorations!(ax)
    hideydecorations!(ax)
    hidespines!(ax)
    #}}}
    return fig
end


##############################################################
function weight_to_linewidth(w, ws; p = dp, use_log = false)
    if use_log
        log_ws = log10.(ws)
        log_w = log10(w)
        log_min = minimum(log_ws)
        log_max = maximum(log_ws)
        t = (log_w - log_min) / (log_max - log_min)
    else
        w_min = minimum(ws)
        w_max = maximum(ws)
        t = (w - w_min) / (w_max - w_min)
    end
    
    return p.min_linewidth + t * (p.max_linewidth - p.min_linewidth)
end


function cm2px(w, h; dpi = 96)
    w_px = floor(Int, w * dpi / 2.54)
    h_px = floor(Int, h * dpi / 2.54)
    w_px = w_px / 1.333327787
    h_px = h_px / 1.333327787
    return (w_px, h_px)
end


end

