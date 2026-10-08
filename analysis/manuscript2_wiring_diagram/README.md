# Manuscript 2 — A synaptic wiring diagram of the superficial dorsal horn

*A synaptic wiring diagram underlying somatosensory modality-specific ascending pathways of
the superficial spinal dorsal horn* (Xiang et al.)

Figure 5 is the connectome-constrained firing-rate model, written in Julia and kept in
[`../dorsal_horn_model/`](../dorsal_horn_model/) rather than here.

> **Which CAVE table to use.** ExN and PN cell-type labels come from
> **`spinalcord_neuron_clusters_2`** (this is where the final, expert-curated cluster
> assignments were deposited); sensory afferent subtypes from **`sensory_neuron_label`**;
> all synapse queries from **`synapses_version3`**. Notebooks written earlier in the project
> may still query `spinalcord_neuron_clusters`, `excitatory_neuron_label` or `synapses_v2`
> — substitute the tables above. See the
> [top-level README](../README.md#canonical-tables).
>
> The per-neuron tables in `data/` (sensory input counts, adjacency, overlap, cluster
> assignments) are snapshots of these CAVE queries, provided so the plotting notebooks run
> without a full re-query. They are not independent primary data.

---

## Figure 1 / S1 — Reconstruction overview

`figure1_reconstruction_overview/`

| File | Panel | Notes |
| --- | --- | --- |
| `01_reconstruction_summary_statistics.ipynb` | Fig. 1e,f; Methods | Neuron counts, total and per-neuron axon cable length (0.964 m total, 5.6 mm mean), ROI extent. |
| `02_exn_ihn_pn_soma_distribution.ipynb` | S1g, S1i | Soma positions of ExNs, IhNs and PNs. Reads `data/exn_ihn_soma_coords.json`. |
| `03_soma_depth_by_cluster.ipynb` | S2d | ExN soma depth from the white matter boundary, by cluster. |

---

## Figure 2 / S2 — ExN subtypes from sensory input + morphology

`figure2_exn_clustering/figure2bde_S2b_exn_sensory_and_clustering.ipynb` produces panels
**2b, 2d, 2e and S2b** in one pass. **Runs offline — no CAVE access.** Panels 2a and 2c are
a schematic and skeleton renderings.

All four panels rest on one table, `data/sensory_input_features/excitatory_neurons.csv`:
one row per laminae I–II excitatory interneuron (n = 177) with its per-subtype sensory
synapse counts, its non-sensory (`N/A`) count, and its subtype label.

| Panel | Content | Extra input |
| --- | --- | --- |
| **2b** | Per-ExN heatmap, fraction of sensory input from each afferent subtype | — |
| **2e** | Violin, sensory as a fraction of *total* input, by subtype | — |
| **2d** | UMAP of morphology + sensory-input features with HDBSCAN clusters | `data/swc_skeletons/` |
| **S2b** | Z-score heatmap of the 33 features used for the embedding | `data/swc_skeletons/` |

### Labels

The CSV's `cell_type` column is a snapshot of the CAVE table
**`spinalcord_neuron_clusters_2`** (the final expert-curated assignments). Verified against
CAVE at materialization 669: 156 of the 177 neurons are still present under their recorded
root ID and all 156 labels match; the remaining 21 have drifted root IDs and, once resolved
through the ChunkedGraph, **all 21 match as well**. Using the CSV is therefore equivalent to
querying CAVE and avoids silently dropping those 21. The notebook carries a
`verify_labels_against_cave()` helper that re-runs the comparison.

The seven stored tags collapse to six published subtypes — `v1` and `v2` are both *Vertical*
here, and are only separated in Figure S5:

| Tag(s) | Subtype | n |
| --- | --- | --- |
| `v1` + `v2` | Vertical | 77 |
| `r` | Radial | 24 |
| `np_module` | C-HTMR/Heat module | 27 |
| `cltmr_module` | C-LTMR module | 24 |
| `adltmr_module` | Aδ-LTMR module | 15 |
| `noci_module` | Aδ-HTMR module | 10 |

### Low-confidence synapses count as sensory

A synapse SPINE calls sensory but whose every subtype probability is ≤ 0.6 is its own
`low_confidence_sensory` category. It is a column in panel b **and is included in the
numerator of panel e**. This is not cosmetic: excluding it changes the Vertical-vs-Radial
comparison from *p* = 2.97 × 10⁻⁶ to 4.69 × 10⁻⁶ and no longer reproduces the paper.

### Reproduces

Rows in panel b keep their order in the source table within each group — the published
panel does no within-group sorting. Panel e's y axis is fixed to 0–0.5 with the
significance brackets stacked inside it.

Group sizes are exactly the published *n* (77 / 24 / 15 / 24 / 10 / 27), the feature matrix
is **177 × 33** (13 morphology + 20 connectivity) as stated in the Methods, and all five
panel-2e *p*-values match the legend exactly:

| Vertical vs | Computed | Published |
| --- | --- | --- |
| Radial | 2.97 × 10⁻⁶ | 2.97 × 10⁻⁶ |
| Aδ-LTMR module | 2.14 × 10⁻⁴ | 2.14 × 10⁻⁴ |
| C-LTMR module | 2.74 × 10⁻¹³ | 2.74 × 10⁻¹³ |
| Aδ-HTMR module | 7.06 × 10⁻⁷ | 7.06 × 10⁻⁷ |
| C-HTMR/Heat module | 1.81 × 10⁻²⁴ | 1.81 × 10⁻²⁴ |

### The CSV ↔ SWC root-ID map matters

The sensory CSV and the SWC exports were made at different times, so **19 of the 177
neurons appear under different root IDs in the two**. `data/sensory_input_features/exn_swc_id_map.csv`
records the correspondence; the notebook loads it and serves it to the ported cells in
place of a live ChunkedGraph lookup (`USE_CAVE = False`; set `True` for the real query).

This is not a convenience — it is load-bearing. Matching on the CSV IDs alone finds SWCs
for only 158 of 177 neurons; the other 19 get no morphology, are zero-filled, and HDBSCAN
then returns **4 clusters instead of 6**. With the map in place the notebook reproduces the
original exactly:

| | Original | This notebook |
| --- | --- | --- |
| neurons aligned | 177 | 177 |
| features retained | 33 | 33 |
| clusters (excl. noise) | 6 | 6 |
| cluster sizes | 15, 23, 80, 23, 9, 27 | 15, 23, 80, 23, 9, 27 |

Two details make the map correct. Only the **seven class folders** are searched
(`adltmr_module`, `cltmr_module`, `noci_module`, `np_module`, `r`, `v1`, `v2`) — they hold
exactly 177 SWCs; including `lamian_i_in`, a 413-file superset, pulls 28 `v1` neurons onto
the wrong file. And the **old ID is tried before its descendants**, because for two neurons
the SWC predates a later re-proofreading and only exists under the old ID. Every neuron's
SWC ends up in the folder matching its label.

Rebuild the map with `data_generation/build_exn_swc_id_map.ipynb`.

### `data_generation/`

| File | Purpose |
| --- | --- |
| `build_excitatory_neurons_table.ipynb` | Rebuilds `excitatory_neurons.csv` from CAVE + SPINE feathers |
| `skeletonize_exn.ipynb` | `pcg_skel` skeletons of the ExNs |
| `export_axon_dendrite_swc.ipynb` | Writes the SWCs in `data/swc_skeletons/` (type 2 = axon, 3 = dendrite) |
| `sensory_to_exn_matrix.py` | Standalone sensory → interneuron matrix with spectral clustering (`#%%` cells) |

---

## Figure 3 / S3 — Lamina I projection neurons

`figure3_pn_clustering/figure3bcd_S3fg_pn_sensory_and_clustering.ipynb` produces panels
**3b, 3c, 3d and S3f, S3g**. **Runs offline — no CAVE access.** Panel 3a is skeleton
renderings.

| Panel | Content |
| --- | --- |
| **S3g** | UMAP + HDBSCAN of the 40 P2 neurons on sensory input — the clustering that *defines* the subtypes |
| **3b** | Violin, sensory as a fraction of total input, by subtype |
| **3c** | Tri-axis scatter of absolute drive from C-Cold, Aδ-HTMR and peptidergic C-fibres |
| **3d** | Per-PN heatmap, fraction of sensory input by afferent subtype |
| **S3f** | Same rows as raw counts, log scale |

### The subtype labels are an output, not an input

`projection_neuron.csv` has a `cell_type` column, but it is the **pre-curation** grouping
(p1 6, p2_poly 18, p2_noci 14, p2_cold 8) and does not match the published *n*. The
clustering reassigns five neurons, giving **P1 6, P2-Cold 9, P2-Noci 16, P2-Poly 15** — the
legend's values, and what `spinalcord_neuron_clusters_2` holds. Checked against CAVE at
materialization 669: the HDBSCAN clusters map one-to-one onto the CAVE labels (0 → p2_cold,
1 → p2_noci, 2 → p2_poly, clean diagonal), while the CSV column disagrees for 5 of 46.
**Use the clusters, not the CSV column.**

### Row order is load-bearing

UMAP's result depends on the order of its input rows, and the original took that order from
an *unsorted* directory glob — filesystem order, which is not portable.
`data/sensory_input_features/pn_preclustering_order.csv` records the exact order that
reproduces the published clustering, and the notebook reads it rather than re-globbing.

The folders it refers to, `data/swc_skeletons_pn_preclustering/`, hold the **pre**-clustering
grouping — the input to this analysis — whereas `data/swc_skeletons/` holds the
post-curation grouping, its output. Both matter:

| Row order used | Clusters found |
| --- | --- |
| **recorded order (correct)** | **9, 16, 15 — ARI 1.0 vs published** |
| CSV order | 8, 32 |
| sorted filenames | 9, 8, 7, 8 + 8 noise |

The notebook asserts `ARI == 1.0` against the shipped assignment and that cluster numbering
is identical, so none of this can regress silently. Only the SWC *filenames* are used, for
ordering; the original's "drop skeletons with zero dendrite length" filter passes all 46,
so the file contents are not read.

### Clustering parameters

Sensory features only — PNs are frequently truncated at the volume boundary, so morphology
is unreliable for them (Methods). P1 is excluded, having been set aside on morphology
beforehand. Features are per-subtype fractions (CNS column dropped), log1p counts and the
CNS fraction — 18 in total; raw counts are deliberately not included. Aδ-HTMR and
peptidergic C-fibres are merged into one nociceptor category first. Then
UMAP(`n_neighbors=4`, `min_dist=0.2`, `random_state=42`) and
HDBSCAN(`min_cluster_size=7`, `min_samples=1`, `cluster_selection_method="leaf"`).

### Reproduces

| Subtype | n | median sensory % | vs P1 (computed) | vs P1 (published) |
| --- | --- | --- | --- | --- |
| P1 | 6 | 2.6 | — | — |
| P2-Poly | 15 | 11.7 | 0.072 | 0.072 |
| P2-Noci | 16 | 22.1 | 4.11 × 10⁻⁴ | 4.11 × 10⁻⁴ |
| P2-Cold | 9 | 42.5 | 6.62 × 10⁻⁷ | 6.62 × 10⁻⁷ |

---

## Figure S3a–e — soma volume and dendritic morphology

`figure3_pn_clustering/figureS3abcde_soma_volume_and_dendrite_morphology.ipynb` produces all
five panels. **Runs offline — no CAVE access.**

| Panel | Content | Input |
| --- | --- | --- |
| **a** | Soma-volume histogram, interneurons vs PN | `data/lamina1_soma_volumes_15um.csv` |
| **b** | Sholl profile (mean ± SEM) in the horizontal plane | `data/swc_skeletons/` |
| **c** | All dendritic arbors overlaid on the soma (rostrocaudal × medial-lateral) | `data/swc_skeletons/` |
| **d** | Dendritic territory — convex-hull area of the horizontal projection | `data/swc_skeletons/` |
| **e** | Dendritic branch density — branch points per µm | `data/swc_skeletons/` |

Soma volumes are measured on a 15 µm cube cropped around the soma centroid, with two manual
corrections: `720575940883198167` uses its 20 µm crop (the 15 µm cube clipped it), and
`720575940895281326` is excluded as a partial reconstruction. Neurons tagged `o` are dropped.
That leaves **45 PN and 411 interneurons** for panel a.

### Panels b–e use the 85 interneurons with a labelled axon

Panels b–e need the dendrite in isolation, which is only possible where the export labels the
axon (SWC `type == 2`). Of the 411 interneuron skeletons, **85** have one; the analysis is
restricted to those, and the notebook asserts the count. This is the set behind the published
panels — the violin scatter in d and e has exactly 85 interneuron and 46 PN points.

The restriction is load-bearing, not cosmetic. Including the other 326 would fold each one's
axon into its "dendrite" and inflate every measurement — median territory rises from
20,118 to 28,401 µm², and the Sholl curve peak from 9.7 to 11.2 intersections.

PNs are handled differently: their exports contain no axon, so every non-soma node is
dendrite, giving all 46.

### Reproduces

All statistics recorded in the original analysis are reproduced exactly, and asserted:

| Quantity | PN | Interneuron | Test |
| --- | --- | --- | --- |
| Dendritic territory (median, µm²) | 66,388 | 20,118 | U = 3658.0, *p* = 2.23 × 10⁻¹⁶ |
| Branch density (median, µm⁻¹) | 0.0106 | 0.0123 | U = 1401.0, *p* = 7.61 × 10⁻³ |
| Sholl AUC (median) | 1722.5 | 1240.0 | *p* = 2.30 × 10⁻⁵ |
| Sholl max radius (median, µm) | 415 | 210 | *p* = 4.63 × 10⁻¹³ |
| Sholl peak (median) | 9.5 | 13.0 | *p* = 1.47 × 10⁻³ |
| Sholl peak radius (median, µm) | 95 | 50 | *p* = 5.10 × 10⁻⁷ |

### `data_generation/`

| File | Purpose |
| --- | --- |
| `build_projection_neuron_table.ipynb` | Rebuilds `projection_neuron.csv` from CAVE + SPINE feathers |
| `quantify_soma_volume_15um.ipynb` | Measures soma volume from the neuron meshes on a 15 µm crop — the code behind `data/lamina1_soma_volumes_15um.csv` |
| `quantify_soma_volume_20um_first_pass.ipynb` | The earlier 20 µm-crop pass; its volumes are the `soma_volume_um3` column, still used for the one clipped soma |

---

## Figure S5 / S6 — Two connectivity-defined vertical cell subtypes

`figureS5_vertical_cell_subtypes/`

| File | Panel | Notes |
| --- | --- | --- |
| `01_vertical_cell_sensory_input.ipynb` | S5c, S5d, S5e | Sensory input fractions and raw Aδ-HTMR counts onto Vert1 vs Vert2. |
| `02_vert1_vert2_cosine_clustering.ipynb` | S5a, S5b | Cosine similarity of ExN input profiles across vertical cells; spectral clustering into Vert1 / Vert2. |
| `03_vertical_cell_spine_density.ipynb` | S6e | Dendritic spine density: NPFF⁺, GRPR⁺ (light microscopy) vs Vert1, Vert2 (EM). |

---

## Figure 4 / S7 — Wiring diagram

`figure4_wiring_diagram/` — four notebooks, one per group of panels. **All run offline; none
queries CAVE.** Panels 4a and 4d are schematics.

| File | Panels | Environment |
| --- | --- | --- |
| `figure4bcf_S7abc_connectivity_matrices.ipynb` | 4b, 4c, 4f, S7a, S7b, S7c | `cave310` |
| `figureS7def_sensory_to_exn_matrices.ipynb` | S7d, S7e, S7f | `cave310` |
| `figure4e_S7h_network_traversal.ipynb` | 4e, S7h | **`pymaid`** (needs `navis`) |
| `figure4g_S7g_wiring_diagram.ipynb` | S7g, 4g edge table | `cave310` |

Run the first two before the last one — it reads the cell-type tables they write into
`figure4_outputs/`.

### What each panel is

Everything rests on two neuron × neuron matrices in `data/adjacency_overlap/`
(`adjacency_matrix.csv`, synapse counts; `overlap_matrix.csv`, axon–dendrite path overlap in
µm; 177 presynaptic × 223 postsynaptic), plus `tag_of.csv` mapping each root ID to its cell
type. The sensory panels use the parallel pair in `data/sensory_exn_matrix/`.

| Panel | Content | Colour cap |
| --- | --- | --- |
| **4b** | Neuron-by-neuron synapse count | P95 of non-zero = **8** |
| **4c** | Normalized synapse density, synapses per 10 µm overlap | P95 = **0.806** |
| **4f** | Cell-type bubble grid: area = avg synapses per overlapping pair, colour = mean density | P95 |
| **S7a** | Axon–dendrite path overlap (µm) | P99 = **879** |
| **S7b** | Cell-type block: connection probability among overlapping pairs | — |
| **S7c** | Cell-type block: mean synapses per overlapping pair | — |
| **S7d** | Sensory → ExN/PN synapse count (all 161 afferents) | P80 of non-zero = **7** |
| **S7e** | Sensory → ExN/PN path overlap (136 afferents; Aβ-LTMR not run) | P95 = **591** |
| **S7f** | Sensory → ExN/PN bubble grid | fixed 0–0.14 syn/10 µm |
| **S7g** | Mean sensory input fraction per cell type (SPINE) | fixed 0–50 % |
| **S7h** | Cell-type mean traversal layer vs `THRESH` | — |
| **4e** | Violin of per-neuron mean traversal layer, one row per cell type | — |

Each colour cap is recomputed from the shipped data and `assert`ed against the value printed
on the published colorbar, so a silent change in the input cannot pass unnoticed.

### Block statistics are over *overlapping* pairs

S7b and S7c are computed over every ordered neuron pair in the block with `overlap > 0`:
S7b is the fraction of those that also have ≥ 1 synapse, S7c the mean synapse count over the
same set (so an overlapping-but-unconnected pair contributes a 0). **Self-pairs are not
excluded** from the within-type blocks — the original analysis includes the diagonal, and
the published numbers only reproduce with it in. Eight values quoted in the Results are
checked in the notebook:

| pre → post | p(connected) | paper | syn/pair | paper |
| --- | --- | --- | --- | --- |
| v1 → p1 | 0.723 | 0.727 | 2.01 | 2.01 |
| v1 → v1 | 0.426 | 0.428 | 0.91 | 0.91 |
| r → v1 | 0.599 | 0.600 | 1.51 | 1.51 |
| v2 → v1 | 0.420 | 0.420 | 1.21 | 1.21 |
| C-LTMR module → p2_poly | 0.164 | 0.160 | 0.36 | 0.36 |
| C-HTMR/Heat module → p2_poly | 0.240 | 0.240 | 0.43 | 0.42 |
| v2 → p2_poly | 0.167 | 0.170 | 0.31 | 0.31 |
| C-HTMR/Heat module → p2_noci | 0.267 | 0.270 | 0.44 | 0.44 |

### Network traversal reproduces exactly

`navis.models.TraversalModel` after Schlegel et al. (2021), run on the 136 sensory + 177 ExN
→ 177 ExN + 46 PN matrix, seeded from the sensory afferents, 10,000 iterations at a fixed
NumPy seed, at `THRESH` ∈ {0.30, 0.50}, plus a 7-point threshold sweep at 2,000 iterations.

Layer means are **bit-identical** to the shipped reference in `data/traversal_outputs/`
(max absolute difference 0.0 across all three output tables). Note the traversal layer is a
*mean first-encounter step*, not a shortest-path distance: a neuron reachable only through
weak edges lands deeper than a strongly driven one at the same hop count.

One portability fix: `model.summary` is indexed by `node` in `navis` 1.4.0 but returns
`node` as a column in other versions; the notebook normalizes both.

### Figure 4g — what is and is not reproduced

The diagram's **node layout was drawn by hand in Illustrator** and is not regenerated. What
the notebook does regenerate is the data-driven part — which arrows exist, their thickness
and their colour — written to `figure4_outputs/figure4g_wiring_edges.csv`.

Every arrow must clear the anatomical filter to be drawn solid: mean synapses per overlapping
pair > 0.25 **and** density > 0.0036 syn/µm (the Methods' "≥ 0.036 per 10 µm"). Both halves
apply to sensory and ExN arrows alike. Sensory arrows have an additional requirement — the
subtype must supply ≥ 3.5 % of the target's input (Fig. S7g); a sensory pair that clears that
but fails the anatomical filter is drawn **dashed**.

Getting the dashed set right depends on one label alignment: the peptidergic C-fibre
population is `Peptidergic C-Nociceptor` in the overlap tables and `Peptidergic C-fibers` in
the SPINE table. Without renaming, its four pairs drop out of the merge.

Thickness and colour bins are exactly those printed in the published legend. Note the
thickness ladder has **six** steps, not five — the sub-threshold `< 0.25` bin is its own
hairline width:

| Avg syn count | Thickness | | Syn. density (per 10 µm) | Colour |
| --- | --- | --- | --- | --- |
| < 0.25 | 0.3 pt (sub-threshold) | | < 0.036 | grey, dashed |
| 0.25–0.5 | 0.5 pt | | 0.036–0.07 | `#c9ddf0` |
| 0.5–1 | 1 pt | | 0.07–0.14 | `#77b5d9` |
| 1–2 | 2 pt | | > 0.14 | `#083e81` |
| 2–4 | 4 pt | | | |
| > 4 | 8 pt | | | |

### Verified against the published panel

The final panel was exported as vector SVG, so each arrow's stroke colour and width are
recoverable from it directly. Decoding that export gives **31 arrows — 28 solid + 3 dashed**,
and the generated edge table matches it in **all 15 (colour, thickness) buckets**:

| Colour | 0.3 | 0.5 | 1 | 2 | 4 | 8 |
| --- | --- | --- | --- | --- | --- | --- |
| `#c9ddf0` | — | 4 | 3 | 2 | 1 | — |
| `#77b5d9` | — | 3 | 2 | 3 | 1 | 1 |
| `#083e81` | — | — | 2 | 2 | 2 | 2 |
| grey, dashed | 2 | — | — | 1 | — | — |

The notebook asserts the per-colour solid counts (10 / 10 / 8) and the dashed count (3), so
neither threshold nor bin edge can drift silently.

### `data_generation/`

The notebooks that build the shipped matrices from CAVE. They still carry the original
absolute output paths — adjust them before re-running.

| File | Writes |
| --- | --- |
| `build_axon_dendrite_path_overlap.ipynb` | `overlap_matrix.csv`. `pcg_skel` skeletons, `pcg_meshwork` axon/dendrite assignment, resampled to ~100 nm internode spacing, `cKDTree` query at 5 µm, contiguous regions grouped by skeleton connectivity, geodesic diameter per region. |
| `build_neuron_by_neuron_adjacency.ipynb` | `adjacency_matrix.csv`, `tag_of.csv` |
| `build_sensory_to_exn_matrix.ipynb` | `data/sensory_exn_matrix/*` |
| `build_sensory_exn_to_exn_pn_matrix.ipynb` | `data/sensory_exn_to_exn_pn_matrix.csv`, `..._cell_tags.csv` |

---

## Data

`data/`

| Path | Contents |
| --- | --- |
| `sensory_input_features/` | `excitatory_neurons.csv` (per-ExN sensory synapse counts), `projection_neuron.csv` (per-PN), `pn_p2_hdbscan_clusters.csv`, `exn_swc_id_map.csv`, `pn_preclustering_order.csv`. |
| `swc_skeletons/` | Curated soma-corrected SWC skeletons, grouped by subtype (see above). 639 files, 25 MB. |
| `exn_combined_features_table.csv` | The combined morphology + sensory input feature table for the 177 ExNs. |
| `adjacency_overlap/` | `adjacency_matrix.csv` (synapse counts), `overlap_matrix.csv` (axon–dendrite path overlap, µm), cell-type summary tables, and the null-model outputs. |
| `adjacency_overlap_v1/` | Earlier version of the same two matrices, kept for provenance. |
| `interneuron_matrix/` | `neuron_pairs_with_axon_length.csv`. |
| `sensory_exn_matrix/` | Sensory → ExN cell-to-cell long-form tables, all possible pairs with overlap, node/edge tables for the wiring diagram. |
| `sensory_exn_to_exn_pn_matrix.csv`, `..._cell_tags.csv` | Input to the traversal model. |
| `traversal_matrix/`, `traversal_outputs/` | Traversal inputs and per-neuron / per-cell-type layer assignments, incl. the threshold sweep. |
| `projection_neuron/` | PN sensory input fractions by group. |
| `sensory_modules/` | Cosine-similarity orderings and the z-scored ExN clustering feature table. |
| `lamina1_soma_volumes_15um.csv` | Soma volumes for 456 lamina I neurons (15 µm crop, repaired mesh). |
| `pn_soma_coords.csv`, `ihn_soma_coords.csv`, `exn_ihn_soma_coords.json` | Soma coordinates. |
| `wiring_diagram_nodes.csv`, `wiring_diagram_edges.csv` | Node/edge tables behind Fig. 4g. |
| `projection_neuron_groups.csv`, `cluster_counts_summary.csv` | Final group assignments and counts. |

