# Dorsal horn connectome analysis

Analysis code and derived data tables for the two manuscripts described in the
[repository README](../README.md), studying the mouse lumbar spinal cord dorsal horn
reconstructed from the **dSC1** serial-section TEM volume:

- **Manuscript 1** — *Cell type-resolved connectomic reconstruction of the spinal dorsal horn reveals sensory neuron synaptic organization and modality-specific inhibition*
  (Xiang et al.) → [`manuscript1_sensory_synaptic_organization/`](manuscript1_sensory_synaptic_organization/)
- **Manuscript 2** — *A synaptic wiring diagram underlying somatosensory modality-specific
  ascending pathways of the superficial spinal dorsal horn* (Xiang et al.)
  → [`manuscript2_wiring_diagram/`](manuscript2_wiring_diagram/)

Each manuscript directory is organized by figure. Every figure directory contains the
notebooks/scripts that produce that figure's panels, numbered in the order they should be
run. See the per-manuscript README for a panel-by-panel map.

## What is *not* here

- **SPINE** (the sensory-neuron subtype classifier: SegCLR + GNN + MLP, Manuscript 1 Fig. 2)
  is maintained separately at **[htem/dorsalhorn-ml](https://github.com/htem/dorsalhorn-ml)**.
  The analyses here consume its per-neuron outputs
  (`cell_type_predictions_<root_id>.feather`).
- **Segmentation and synapse prediction** pipelines (CNN-based, run by collaborators).
- Two large data sets that some Manuscript 1 panels read directly — the per-neuron SPINE
  prediction tables and the meshwork caches — plus neuron meshes and full-resolution TIFF
  exports. These are at Harvard Dataverse, <https://doi.org/10.7910/DVN/NVXK7W>; see
  **Data availability** below.

## Data sources

### Browsing the volume — Neuroglancer

The EM data and reconstructed neurons can be browsed in a web browser, without installing
anything, through this Neuroglancer state:

<https://spelunker.cave-explorer.org/#!middleauth+https://global.daf-apis.com/nglstate/api/v1/5871146013032448>

Start here if you want to look up a particular neuron by root ID, inspect a synapse, or check
a reconstruction against the EM, rather than to re-run an analysis.

The volume EM data and segmentation — the dSC1 dataset and all dSC_APEX datasets — are also
deposited at [BossDB](https://bossdb.org).

### CAVE

Most analyses query the segmentation and annotation tables live through
[CAVE](https://github.com/CAVEconnectome/CAVEclient):

```python
import caveclient
client = caveclient.CAVEclient("wclee_mouse_spinalcord_cltmr")
```

Segmentation source:
`graphene://https://cave.fanc-fly.com/segmentation/table/wclee_mouse_spinalcord_cltmr`

#### Canonical tables

These are the authoritative tables. **Use these** when reproducing or extending any analysis
here, even where an individual notebook still references an older table name:

| Table | Contents | Use it for |
| --- | --- | --- |
| `synapses_version3` | Automated synapse predictions (cleft + pre/post partner assignment) | All synapse queries |
| `sensory_neuron_synapses` | Synapses of reconstructed sensory afferents | Sensory neuron synapse counts in **Manuscript 1, Figure 3** |
| `sensory_neuron_label` | Sensory afferent subtype annotations | Sensory neuron cell types |
| `spinalcord_neuron_clusters_2` | Final ExN and PN cluster labels | ExN and PN cell types (both manuscripts) |

Every neuron-level quantity reported in the papers can be regenerated from these four
tables — the cell-type lists, synapse counts, and per-neuron partner tables in `data/` are
convenience snapshots of CAVE queries, not independent primary data.

#### Superseded table names still present in the code

Individual notebooks were written over ~2 years and may query earlier versions. Where you
see these, substitute the canonical table above:

| Appears in code | Superseded by |
| --- | --- |
| `spinalcord_neuron_clusters` | `spinalcord_neuron_clusters_2` |
| `synapses_v2` | `synapses_version3` |
| `excitatory_neuron_label` | `spinalcord_neuron_clusters_2` (excitatory/inhibitory identity is carried in the cluster labels) |
| `all_neuron_label`, `manual_synapses` | still valid for specific sub-analyses; see the notebook that uses them |

Some analyses use additional tables beyond the four canonical ones. Rather than duplicating
a list that would drift, **the table name is always visible in the `client.materialize.query_table(...)`
call in the notebook itself** — grep for `query_table` in any figure directory to see exactly
what it reads.

#### Root IDs drift — pin the materialization

A CAVE root ID identifies a *version* of a neuron, not the neuron itself: it changes every
time the segment is edited. The neuron lists hard-coded in these notebooks were recorded
when the figures were made, so as proofreading continues some of them go stale. **A stale
root ID does not raise an error — the query simply returns nothing**, the neuron silently
contributes a zero, and the subtype mean drops.

Two defences are built into the CAVE-querying notebooks:

1. **The materialization is pinned**, e.g. `MATERIALIZATION = 669` (the permanent release
   version; the older versions the figures were originally made against have since expired).
   Set it to `None` to query the latest instead.
2. **A `resolve_current_ids()` helper** checks every published root ID with
   `chunkedgraph.is_latest_roots()` and substitutes the ChunkedGraph's suggested successor
   for any that are stale, printing each substitution. As of materialization 669, 3 of the
   79 sensory neurons need substituting.

If you add a new analysis, do the same — do not query hard-coded root IDs directly.

#### Known fix landing after materialization 669

Two `islet-cltmr` inhibitory neurons (Manuscript 1, Figure 5b) were still merged in version
669, leaving their axon and dendrite in a single segment. They have since been **split
manually**, so a materialization newer than 669 is preferred for that panel. Details and
the regeneration step are in
[the Figure 5 section](manuscript1_sensory_synaptic_organization/README.md#panel-b-needs-a-materialization-newer-than-669).
No other analysis in this repository is affected.

### CATMAID — the dSC_APEX volumes

The per-glomerulus quantifications in Manuscript 1 Fig. 3g use manually traced skeletons in
the **six dSC_APEX volumes**, hosted on a CATMAID server and accessed programmatically with
[`pymaid`](https://pymaid.readthedocs.io/).

<!-- LINK-TODO: replace with the public CATMAID project URL for the six dSC_APEX volumes -->

> **Link to be added.** A CATMAID instance hosting the six dSC_APEX volumes together with the
> glomerulus tracings behind Fig. 3g will be published here, so the traced skeletons can be
> browsed and checked against the EM directly. In the meantime the dSC_APEX image volumes
> themselves are available at [BossDB](https://bossdb.org).

Once the link is live, set the server and your API token as described under
[Credentials](#credentials); the notebooks in
[`figure3_glomerulus_quantification_catmaid/`](manuscript1_sensory_synaptic_organization/figure3_glomerulus_quantification_catmaid/)
connect through `pymaid` and expect both.

### Credentials

**No credentials are committed to this repository.** Every place a token was previously
hard-coded now contains a `REPLACE_WITH_YOUR_...` placeholder. Set them up once:

```bash
export CAVE_TOKEN="your-cave-token"
export CATMAID_API_TOKEN="your-catmaid-token"
```

and see [`common/credentials.py`](common/credentials.py) for helpers that read them. For
CAVE, the recommended route is the standard credentials file rather than an inline token:

```python
import caveclient
caveclient.auth.AuthClient().save_token(token="your-cave-token")
```

## Layout

```
common/                                  shared library + mesh export utilities
envs/                                    exact conda environment exports
dorsal_horn_model/                       connectome-constrained circuit model (Julia, M2 Fig. 5)
manuscript1_sensory_synaptic_organization/
  figureS1_synapse_prediction_validation/     precision / recall of synapse predictions
  figure3_sensory_synaptic_architecture/      synapse counts per sensory subtype
  figure3_glomerulus_quantification_catmaid/  per-glomerulus axoaxonic counts (dSC_APEX)
  figureS6_synapse_distribution/              synapse position along the axon arbor
  figure4_homotypic_psi/                      homotypic presynaptic inhibition circuits
  figureS9_lateral_inhibition/                lateral vs self inhibition by PSI-IhNs
  figure5_feedforward_inhibition/             heterotypic feedforward inhibition
  data/                                       tables and per-neuron files read by the above
manuscript2_wiring_diagram/
  figure1_reconstruction_overview/            reconstruction statistics, soma distributions
  figure2_exn_clustering/                     ExN sensory-input + morphology UMAP/HDBSCAN
  figure3_pn_clustering/                      PN identification (soma volume) and clustering
  figureS5_vertical_cell_subtypes/            Vert1 / Vert2 split
  figure4_wiring_diagram/                     connectivity, path overlap, network traversal
  data/
```

## Circuit model

[`dorsal_horn_model/`](dorsal_horn_model/) holds the connectome-constrained firing-rate model
of the dorsal horn microcircuit behind **Manuscript 2, Figure 5** — written in **Julia**, not
Python, and so independent of the conda environments below. Install its dependencies with
`julia install.jl` and see
[`dorsal_horn_model/README.md`](dorsal_horn_model/README.md) for how to produce each panel.
It is authored and maintained by Ramin Khajeh and carries its own MIT licence.

## Environment

Two conda environments, both exported directly from the environments the analyses were run
in (so versions are exact, not reconstructed):

```bash
conda env create -f envs/cave310.yml    # Python 3.10 — CAVE analyses, clustering, plotting
conda env create -f envs/pymaid.yml     # Python 3.7  — CATMAID, navis traversal, pymeshfix
```

Most notebooks run under `cave310`. See [`envs/README.md`](envs/README.md) for the
per-notebook breakdown and key package versions.

Scripts with `#%%` cell markers (the `*.py` files in the figure directories) are written to
be run cell-by-cell in Spyder or VS Code rather than end to end.

## Hard-coded paths

Notebooks were written against the original working directories
(`/Users/wangchu/connectivity_analysis`, `/Users/wangchu/mesh_analysis`). Absolute paths
remain in many cells; adjust them to your checkout, or set a shell variable and edit the
path constants at the top of each notebook. Paths to data files that ship with this repo
are listed in each manuscript's README.

## Data availability

The full EM volume, segmentation, and synapse predictions are served through CAVE (above).
To browse rather than query them, see
[Neuroglancer](#browsing-the-volume--neuroglancer) for dSC1 and
[CATMAID](#catmaid--the-dsc_apex-volumes) for the six dSC_APEX volumes and the glomerulus
tracings.

**What ships in this repository (~91 MB).** Every figure notebook for **Manuscript 2**, and
most of Manuscript 1, runs offline from a clone — no CAVE account and no extra download.
That covers the derived connectivity and overlap matrices, the per-neuron sensory-input
tables, cluster assignments, soma volumes, the curated SWC skeletons, and the per-PSI-IhN
SPINE tables, all under each manuscript's `data/`.

**What is too large for git** is deposited at Harvard Dataverse:

> **<https://doi.org/10.7910/DVN/NVXK7W>**

| Deposited | Size | Required by |
| --- | --- | --- |
| `psi_ihn_postsynaptic_feather/` | 132 MB | M1 Fig. 5f–j, S11, S12 |
| `meshwork_h5/` | 44 MB | M1 Fig. S6b, S6c, S6d–h, S6j |

These are **not** optional extras for those panels — the notebooks read them directly. To use
them, download the deposit and unpack the two directories into
`manuscript1_sensory_synaptic_organization/data/`, keeping those folder names. No other panel
needs anything beyond this repository.

The Dataverse deposit also contains a complete copy of this code, so it can be used standalone
without cloning. This repository is where corrections are made; prefer it for the code itself.

Neuron meshes (`.ply`) and full-resolution TIFF exports are not deposited: the mesh figures
(M1 Fig. 3h,i) fetch geometry from CloudVolume at runtime, so the meshes are outputs rather
than inputs.

The SPINE classifier that produced the prediction tables is at
[htem/dorsalhorn-ml](https://github.com/htem/dorsalhorn-ml).

## Citation

If you use this code, please cite the corresponding manuscript.
