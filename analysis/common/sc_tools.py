# -*- coding: utf-8 -*-
"""
Created on Sat Jan 25 17:36:00 2025

@author: Adrian Holtrup
"""

import itertools
from collections import defaultdict
import numpy as np
import scipy
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import SpectralClustering, AgglomerativeClustering
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def aff_bipartition(G, pre_1, pre_2, post):
    
    '''
    

    Parameters
    ----------
    G : networkx digraph
        Directed connectivity graph.
    pre_1 : list
        First list of pre neurons after which to partition post.
    pre_2 : list
        Second list of pre neurons after which to partition post.
    post : list
        List of post neurons that are partitioned based on pre_1, pre_2.

    Returns
    -------
    V1 : list
        Partition of post that is more affine to pre_1.
    V2 : list
        Partition of post that is more affine to pre_2.
    V1_affs : list
        Affinity indices of elements in V1.
    V2_affs : list
        Affinity indices of elements in V2.

    '''
    V1 = []
    V2 = []
    V1_affs = []
    V2_affs = []
    fac_r = len(pre_1)
    fac_c = len(pre_2)
    for v in post:
        con_r = sum([w['weight'] 
                     for n, w in G.pred[v].items() 
                     if n in pre_1])
        con_c = sum([w['weight'] 
                     for n, w in G.pred[v].items() 
                     if n in pre_2])
        try:
            aff = (fac_c * con_r - fac_r * con_c) / (fac_c * con_r + fac_r * con_c)
            if con_r > con_c:
                V1.append(v)
                V1_affs.append(aff)
            else:
                V2.append(v)
                V2_affs.append(aff)
        except:
            continue
    return V1, V2, V1_affs, V2_affs


def build_matrix(G, pre_list, post_list):
    
    '''
    

    Parameters
    ----------
    G : networkx digraph
        Directed connectivity graph.
    pre_list : list
        List of pre neuron lists.
    post_list : list
        List of post neuron lists.

    Returns
    -------
    con_mat : numpy array
        Connectivity matrix where rows are post, columns pre.
    pre_lens : list
        List of lengths of pre neuron lists.
    post_lens : list
        List of lengths of post neuron lists.

    '''
    
    pre_lens = list(map(len, pre_list))
    post_lens = list(map(len, post_list))
    
    n, m = sum(pre_lens), sum(post_lens)
    con_mat = np.zeros((m, n))
    i = 0
    for pre_type in pre_list:
        for pre in pre_type:
            j = 0
            for post_type in post_list:
                for post in post_type:
                    if G.has_edge(pre, post):
                        weight = G.edges[pre, post]['weight']
                        con_mat[j,i] = weight
                    j += 1
            i += 1        
    return con_mat, pre_lens, post_lens 


def cossim_matrix(G, pre_list, post_list):
    
    '''
    

    Parameters
    ----------
    G : networkx digraph
        Directed connectivity graph.
    pre_list : list
        List of pre neuron lists.
    post_list : list
        List of post neuron lists.

    Returns
    -------
    con_mat : numpy array
        Square matrix of output cosine similarities.
    pre_lens : list
        List of lengths of pre neuron lists.
    post_lens : list
        List of lengths of post neuron lists.

    '''
    
    pre_lens = list(map(len, pre_list))
    post_lens = list(map(len, post_list))
    
    pre_flat = list(itertools.chain.from_iterable(pre_list))
    post_flat = list(itertools.chain.from_iterable(post_list))
    
    n, m = sum(pre_lens), sum(post_lens)
    
    weight_matrix = np.zeros((n, m))
    for i, pre in enumerate(pre_flat):
        for j, post in enumerate(post_flat):
            if G.has_edge(pre, post):
                weight_matrix[i, j] = G.edges[pre, post]['weight']

    norms = np.linalg.norm(weight_matrix, axis=1, keepdims=True)
    sim_mat = (weight_matrix @ weight_matrix.T) / (norms @ norms.T + 1e-9)

    return sim_mat, pre_lens, post_lens, pre_flat, post_flat


def inp_cossim_matrix(G, post_list, pre_list):
    
    '''
    

    Parameters
    ----------
    G : networkx digraph
        Directed connectivity graph.
    post_list : list
        List of post neuron lists.    
    pre_list : list
        List of pre neuron lists.

    Returns
    -------
    con_mat : numpy array
        Square matrix of input cosine similarities.
    post_lens : list
        List of lengths of post neuron lists.    
    pre_lens : list
        List of lengths of pre neuron lists.

    '''
    
    post_lens = list(map(len, post_list))
    pre_lens = list(map(len, pre_list))

    post_flat = list(itertools.chain.from_iterable(post_list))
    pre_flat = list(itertools.chain.from_iterable(pre_list))
    
    m, n = sum(post_lens), sum(pre_lens)
    weight_matrix = np.zeros((n, m))
    for i, pre in enumerate(pre_flat):
        for j, post in enumerate(post_flat):
            if G.has_edge(pre, post):
                weight_matrix[i, j] = G.edges[pre, post]['weight']
                
    norms = np.linalg.norm(weight_matrix.T, axis=1, keepdims=True)
    sim_mat = (weight_matrix.T @ weight_matrix) / (norms @ norms.T + 1e-9)
    
    return sim_mat, post_lens, pre_lens, post_flat, pre_flat


def plot_dendrogram(model, **kwargs):
    
    '''
    
    
    Adapted from: https://scikit-learn.org/stable/auto_examples/cluster/plot_agglomerative_dendrogram.html

    Parameters
    ----------
    model : scikit-learn hierarchical clustering model
        Substrate for the hierachical tree.

    Returns
    -------
    order : list
        List of cluster and color tuples.

    '''
    
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)
    
    R = scipy.cluster.hierarchy.dendrogram(linkage_matrix, 
                                           get_leaves=True, 
                                           **kwargs)
    
    _, exms = np.unique(model.labels_, return_index=True)
    col_indices = [int(np.where(R['leaves']==exm)[0]) for exm in exms]
    col_dict = {f'{i+1}': R['leaves_color_list'][inx] 
                for i, inx in enumerate(col_indices)}

    return col_dict

def gen_clustering(
        G, 
        pre_list, 
        post_list, 
        MODEL_TYPE='spectral', 
        k=10, 
        AGG_THRESH=None, 
        N_SEEDS=10
        ):
    
    '''
    

    Parameters
    ----------
    G : networkx digraph
        Directed connectivity graph.
    post_list : list
        List of post neuron lists.    
    pre_list : list
        List of pre neuron lists.

    Returns
    -------
    clusts : list
        List of clusters.
        

    '''

    sim_mat, pre_lens, _, pre_flat, _ = cossim_matrix(
        G,
        pre_list, 
        post_list
        )
 
    
    
    if MODEL_TYPE == 'spectral':
        clust_model = SpectralClustering(n_clusters=k, affinity='rbf')
        clust_model.fit(sim_mat)

        

    elif MODEL_TYPE == 'spectral_averaged':
        clust_dict = defaultdict(list)
        for seed in range(N_SEEDS):
            clust_model = SpectralClustering(n_clusters=k, 
                                             affinity='rbf', 
                                             random_state=seed)
            clust_model.fit(sim_mat)
            
            clusts = [[] for _ in range(k)]
            for i, n in enumerate(pre_flat):
                clusts[clust_model.labels_[i]].append(n)   
                
            clust_dict[str(seed)] = clusts
        
        A = clust_dict['0']    
        for seed in range(1, N_SEEDS):
            B = clust_dict[str(seed)]
            cost_mat = np.zeros((k, k))
            
            for i in range(k):
                for j in range(i+1, k):
                    car_inter = len([n for n in A[i] if n in B[j]])
                    car_union = len(A[i]) + len(B[j]) - car_inter
                    J = car_inter / car_union
                    cost_mat[i, j] = cost_mat[j, i] = -J
                    
            _, assignment = linear_sum_assignment(cost_mat)
            B_sorted = [B[inx] for inx in assignment]
        
            clust_dict[str(seed)] = B_sorted
        
        
        clusts = [[] for _ in range(k)]
        for n in pre_flat:
            clust_col = np.empty((N_SEEDS))
            for seed in range(N_SEEDS):
                for clust in range(k):
                    if n in clust_dict[str(seed)][clust]:
                        clust_col[seed] = clust
                        
            clust = int(scipy.stats.mode(clust_col)[0])
            clusts[clust].append(n)            
    
    
    
    elif MODEL_TYPE == 'agglomerative':    
        clust_model = AgglomerativeClustering(n_clusters=None, 
                                              compute_full_tree=True,
                                              distance_threshold=AGG_THRESH)
        clust_model.fit(sim_mat)
        k = clust_model.n_clusters_
        fig, ax = plt.subplots()    
        col_dict = plot_dendrogram(clust_model, 
                                   color_threshold=AGG_THRESH, 
                                   above_threshold_color='black',
                                   ax=ax)
        patches = [mpatches.Patch(color=col_dict[f'{i+1}'], 
                                  label=f'{i+1}')
                   for i in range(k)]
        ax.legend(handles=patches)
        ax.set_xticks([])
        ax.set_ylabel('Euclidean distance') 
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_title('Output cosine similarity dendrogram')
        plt.show() 
    
    
    
    if MODEL_TYPE != 'spectral_averaged':    
        clusts = [[] for _ in range(k)]
        for i, n in enumerate(pre_flat):
            clusts[clust_model.labels_[i]].append(n) 
            
    return clusts        
