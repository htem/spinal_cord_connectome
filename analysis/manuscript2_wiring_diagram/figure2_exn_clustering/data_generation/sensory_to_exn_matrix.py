#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 13:04:16 2025

@author: wangchu
"""

# Libraries
import caveclient
import subprocess
import json
import numpy as np
import matplotlib.pyplot as plt
# Setup caveclient
client = caveclient.CAVEclient("wclee_mouse_spinalcord_cltmr", ) #server_address='https://global.daf-apis.com')

#%% load synapse table and cell label table

interneuron_df=client.materialize.query_table('excitatory_neuron_label')
interneuron_tags = ['c1', 'c2', 'c3', 'c4', 'c1, pep', 'r', 'u']
interneuron_ids = interneuron_df.loc[interneuron_df['tag'].isin(interneuron_tags), 'pt_root_id'].tolist()

sensory_df=client.materialize.query_table('sensory_neuron_label')
sensory_df.loc[sensory_df['pt_root_id'] == 0, 'pt_root_id'] = 720575940887704904
sensory_ids=sensory_df['pt_root_id'].tolist()

sensory_ids
#%% load synapse table
# then pull only synapses where both ends are in your allowed set
synapses_df = client.materialize.query_table(
    'synapses_v2',
    filter_in_dict={
        'pre_pt_root_id':  sensory_ids,
        'post_pt_root_id': interneuron_ids
    }
)


flip_df = client.materialize.query_table(
    'synapses_v2',
    filter_in_dict={
        'pre_pt_root_id':  interneuron_ids,
        'post_pt_root_id': sensory_ids
    }
)
# define the column‐pairs you want to swap
swap_pairs = [
    ('pre_pt_supervoxel_id', 'post_pt_supervoxel_id'),
    ('pre_pt_root_id',        'post_pt_root_id'),
    ('pre_pt_position',       'post_pt_position'),
]

for col_a, col_b in swap_pairs:
    # make copies so the swap doesn’t overwrite itself
    a = flip_df[col_a].copy()
    b = flip_df[col_b].copy()
    flip_df[col_a] = b
    flip_df[col_b] = a
    
import pandas as pd

records = synapses_df.to_dict('records') + flip_df.to_dict('records')
combined_df = pd.DataFrame.from_records(records, columns=synapses_df.columns)
combined_df

#%% plot the connectivity matrix
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# 1) Define precisely the tag-blocks you want, in order:
custom_order = [
    'ab-ltmr',
    'adltmr',
    'cltmr',
    'mrgd',
    'sst',
    'ad-htmr',
    'c-htmr-p',
    'c-cold'
]

# 2) Normalize your sensory_df tags to lowercase, stripped strings:
sensory_df['tag_norm'] = sensory_df['tag'].fillna('').astype(str).str.lower().str.strip()

# 3) Keep only rows with tags in your custom_order:
df_filt = sensory_df[sensory_df['tag_norm'].isin(custom_order)]

# 4) Now build two things in one pass: 
#    a) sorted_ids, by walking through each tag in custom_order  
#    b) tag_groups, which records (tag, start_idx, end_idx) for plotting
sorted_ids  = []
tag_groups  = []
start       = 0

for tag in custom_order:
    ids_for_tag = df_filt[df_filt['tag_norm'] == tag]['pt_root_id'].tolist()
    if not ids_for_tag:
        # skip blocks that happen to have zero members
        continue
    end = start + len(ids_for_tag)
    
    sorted_ids.extend(ids_for_tag)
    tag_groups.append((tag.upper(), start, end))
    
    start = end

# 5) Build your connectivity matrix just as before,
#    but now explicitly with this sorted_ids ordering:
matrix = (
    combined_df
      .groupby(['post_pt_root_id','pre_pt_root_id'])
      .size()
      .unstack(fill_value=0)
)
matrix = matrix.reindex(index=interneuron_ids,
                        columns=sorted_ids,
                        fill_value=0)

# 6) Plot, drawing a single separator & label per block:
plt.figure(figsize=(10,6))
plt.imshow(matrix.values, aspect='auto')
plt.colorbar(label='Number of connections')

plt.xticks([]); plt.yticks([])   # no raw ID ticks

for tag, start, end in tag_groups:
    # separator line between blocks (skip the very first)
    if start != 0:
        plt.axvline(start - 0.5, color='white', linewidth=2)
    # label once, centered
    mid = (start + end - 1) / 2
    plt.text(mid, -0.5, tag, ha='center', va='bottom', fontsize=12)

plt.xlabel('pre_pt_root_id (sensory → interneuron)')
plt.ylabel('post_pt_root_id')

plt.tight_layout()
plt.show()

#%% clustering based on the sensory inputs
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import SpectralClustering
import pandas as pd

# --- 1) Define your tag order and build sorted_ids/tag_groups as before ---
custom_order = [
    'ab-ltmr','adltmr','cltmr','mrgd',
    'sst','ad-htmr','c-htmr-p','c-cold'
]
sensory_df['tag_norm'] = sensory_df['tag'].fillna('').astype(str).str.lower().str.strip()

df_filt = sensory_df[sensory_df['tag_norm'].isin(custom_order)]
sorted_ids = []; tag_groups = []; start = 0
for tag in custom_order:
    ids_for_tag = df_filt[df_filt['tag_norm']==tag]['pt_root_id'].tolist()
    if not ids_for_tag: continue
    end = start + len(ids_for_tag)
    sorted_ids.extend(ids_for_tag)
    tag_groups.append((tag.upper(), start, end))
    start = end

# --- 2) Build & reorder your connectivity matrix so columns follow sorted_ids ---
matrix = (
    combined_df
      .groupby(['post_pt_root_id','pre_pt_root_id'])
      .size()
      .unstack(fill_value=0)
)
matrix = matrix.reindex(columns=sorted_ids, fill_value=0)

# --- 3) Spectral clustering on the rows (interneurons) ---
X = matrix.values
k = 4
sc = SpectralClustering(n_clusters=k, affinity='nearest_neighbors',
                        assign_labels='kmeans', random_state=0)
labels = sc.fit_predict(X)
order = np.argsort(labels)
clustered_X = X[order, :]
clustered_labels = labels[order]
boundaries = np.where(clustered_labels[1:] != clustered_labels[:-1])[0] + 1
clustered_X= np.clip(clustered_X, None, 30)
# --- 4) Plot with custom‐ordered x‐axis and clustered y‐axis ---
plt.figure(figsize=(12, 8))
plt.imshow(clustered_X, aspect='auto')
plt.colorbar(label='Number of connections')
plt.xticks([]); plt.yticks([])

# draw vertical separators & labels for each tag block
for tag, s, e in tag_groups:
    if s != 0:
        plt.axvline(s - 0.5, color='white', linewidth=2)
    mid = (s + e - 1) / 2
    plt.text(mid, -1.0, tag, ha='center', va='bottom', fontsize=10)

# draw horizontal separators & labels for each cluster
for b in boundaries:
    plt.axhline(b - 0.5, color='white', linewidth=2)
for cluster in np.unique(clustered_labels):
    idxs = np.where(clustered_labels == cluster)[0]
    start, end = idxs.min(), idxs.max()+1
    mid = (start + end - 1) / 2
    plt.text(-1.5, mid, f'Cluster {cluster}', ha='right', va='center', fontsize=10)

plt.xlabel('Sensory inputs (grouped by custom tags)')


plt.tight_layout()
plt.show()
#%% return the id in each clusters and tag
import numpy as np

# build a lookup from pt_root_id → tag
tag_lookup = interneuron_df.set_index('pt_root_id')['tag']

# recover the y-axis IDs in clustered order
y_axis_ids = matrix.index[order]

for cluster in np.unique(labels):
    # positions (row indices) in clustered_X that belong to this cluster
    pos = np.where(labels[order] == cluster)[0]
    # the actual pt_root_id values in that block
    cluster_ids  = y_axis_ids[pos]
    # look up their tags
    cluster_tags = tag_lookup.loc[cluster_ids]
    
    print(f"Cluster {cluster} ({len(cluster_ids)} cells):")
    for pid, tag in zip(cluster_ids, cluster_tags):
        print(f"  {pid}: {tag}")
    print()
    
    #%% clustering based on sensory input, ignore individual sensory neuron identity
import numpy as np
import pandas as pd
from sklearn.cluster import SpectralClustering

# --- PARAMETERS ---
custom_order = [
    'ab-ltmr', 'adltmr', 'cltmr', 'mrgd',
    'sst', 'ad-htmr', 'c-htmr-p', 'c-cold'
]
k = 5  # number of clusters

# --- 1) Normalize tags in sensory_df and build lookup dict ---
sensory_df['tag_norm'] = (
    sensory_df['tag']
    .fillna('')
    .astype(str)
    .str.lower()
    .str.strip()
)
tag_lookup = dict(zip(
    sensory_df['pt_root_id'],
    sensory_df['tag_norm']
))

# --- 2) Annotate combined_df with pre-synaptic tag ---
agg = combined_df.copy()
agg['pre_tag'] = agg['pre_pt_root_id'].map(tag_lookup).fillna('unknown')

# --- 3) Build aggregated matrix: rows=post_pt_root_id, cols=pre_tag counts ---
agg_matrix = (
    agg
    .groupby(['post_pt_root_id', 'pre_tag'])
    .size()
    .unstack(fill_value=0)
)

# --- 4) Reindex columns to custom_order (fill missing tags with zeros) ---
agg_matrix = agg_matrix.reindex(columns=custom_order, fill_value=0)

# --- 5) Spectral clustering on these 8-dimensional profiles ---
X = agg_matrix.values
sc = SpectralClustering(
    n_clusters=k,
    affinity='nearest_neighbors',
    assign_labels='kmeans',
    random_state=0
)
labels = sc.fit_predict(X)

# --- 6) Recover the y-axis ordering and cluster labels ---
order = np.argsort(labels)
y_axis_ids    = agg_matrix.index[order]
y_axis_labels = labels[order]

# --- 7) Print pt_root_id for each cluster in y-axis order ---
for cluster in np.unique(labels):
    positions   = np.where(y_axis_labels == cluster)[0]
    cluster_ids = y_axis_ids[positions].tolist()
    print(f"Cluster {cluster} ({len(cluster_ids)} cells):")
    print(cluster_ids)
    print()
#%% plot the matrix
import numpy as np
import matplotlib.pyplot as plt

# --- 1) Build the original connectivity matrix (post × pre) ---
conn_matrix = (
    combined_df
      .groupby(['post_pt_root_id', 'pre_pt_root_id'])
      .size()
      .unstack(fill_value=0)
)
# reorder columns by your custom‐tag order
conn_matrix = conn_matrix.reindex(columns=sorted_ids, fill_value=0)

# --- 2) Reorder rows by the clustering you already computed ---
# ‘labels’ is the array from SpectralClustering on the aggregated matrix
order = np.argsort(labels)
conn_clustered = conn_matrix.iloc[order, :]

# if you need the cluster boundaries for horizontal lines:
boundaries = np.where(labels[order][1:] != labels[order][:-1])[0] + 1
conn_clustered= np.clip(conn_clustered, None, 30)
# --- 3) Plot the heatmap with both tag‐grouped x‐axis and clustered y‐axis ---
plt.figure(figsize=(12, 8))
plt.imshow(conn_clustered.values, aspect='auto')
plt.colorbar(label='Number of connections')
plt.xticks([]); plt.yticks([])

# vertical separators & tag labels (x axis)
for tag, start, end in tag_groups:
    if start != 0:
        plt.axvline(start - 0.5, color='white', linewidth=2)
    mid = (start + end - 1) / 2
    plt.text(mid, -1.0, tag, ha='center', va='bottom', fontsize=10)

# horizontal separators & cluster labels (y axis)
for b in boundaries:
    plt.axhline(b - 0.5, color='white', linewidth=2)
for cluster in np.unique(labels):
    idxs = np.where(labels[order] == cluster)[0]
    start, end = idxs.min(), idxs.max() + 1
    mid = (start + end - 1) / 2
    plt.text(-1.5, mid, f'Cluster {cluster}', ha='right', va='center', fontsize=10)

plt.xlabel('Sensory inputs (grouped by tag)')
plt.ylabel('Post‐synaptic neurons (clustered)')

plt.tight_layout()
plt.show()
#%% print out the interneuron id
import numpy as np

# 1) Build a quick lookup from post‐pt_root_id to its tag
post_tag_lookup = interneuron_df.set_index('pt_root_id')['tag'].to_dict()

# 2) Recover the y‐axis IDs in clustered order
order        = np.argsort(labels)
y_axis_ids   = agg_matrix.index[order]
y_axis_lbls  = labels[order]

# 3) For each cluster, pull out the IDs and their tags
for cluster in np.unique(y_axis_lbls):
    # positions (row indices in clustered X) for this cluster
    pos         = np.where(y_axis_lbls == cluster)[0]
    cluster_ids = y_axis_ids[pos].tolist()
    
    print(f"\nCluster {cluster} ({len(cluster_ids)} cells):")
    for pt_id in cluster_ids:
        tag = post_tag_lookup.get(pt_id, "UNKNOWN")
        print(f"  {pt_id} → {tag}")


