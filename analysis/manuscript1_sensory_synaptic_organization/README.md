# Manuscript 1 — Sensory synaptic organization and modality-specific inhibition

*Cell type-resolved connectomic reconstruction of the spinal cord dorsal horn reveals
somatosensory neuron synaptic organization and modality-specific inhibition* (Xiang et al.)

Figures 1 and 2 (EM volume acquisition; the SPINE classifier) have no analysis code in this
repository — see the top-level [README](../README.md).

---

## Figure S1 — Validation of synapse predictions

`figureS1_synapse_prediction_validation/`

| File | Panel | Notes |
| --- | --- | --- |
| `synapse_prediction_accuracy.ipynb` | S1a | Precision / recall for sensory axodendritic, sensory axoaxonic, and spinal axodendritic predictions. Counts of manually reviewed predictions (2,216 predictions vs 1,932 confirmed synapses) are entered inline in the notebook. |

---

## Figure 3 — Synaptic architectures of somatosensory neuron subtypes

### In the dSC1 volume (CAVE)

`figure3_sensory_synaptic_architecture/`

One notebook per panel. Each is self-contained: style, neuron lists, query, statistics,
plot. Run top to bottom.

| Notebook | Panel | Data source | Statistics |
| --- | --- | --- | --- |
| `figure3d_axodendritic_synapse_counts.ipynb` | **3d** — axodendritic synapses per sensory neuron | CAVE `synapses_version3`, sensory neuron **presynaptic**, autapses dropped, cleft `size` ≥ 20 voxels. 79 neurons (10 per subtype, 9 C-Cold). | Kruskal–Wallis + Dunn (Holm) |
| `figure3e_downstream_neuron_number.ipynb` | **3e** — number of downstream neurons | `data/sensory_postsynaptic_partners.csv`; intact partners only (`num_supervoxels` ≥ 40,000), counted as distinct `post_pt_root_id` per sensory neuron. | Kruskal–Wallis + Dunn (Holm) |
| `figure3f_axoaxonic_synapse_counts.ipynb` | **3f** — axoaxonic synapses per sensory neuron | CAVE `sensory_neuron_synapses`, sensory neuron **postsynaptic**. The 16 manually validated neurons (2 per subtype). | Kruskal–Wallis + Dunn (Holm) |

Panels 3d and 3f are the same analysis with the query direction reversed; they differ only
in which table they read and which neurons they cover:

- **3d** uses the raw predictions (`synapses_version3`) because sensory→spinal axodendritic
  predictions are accurate (Figure S1a).
- **3f** uses `sensory_neuron_synapses`, the **manually curated** table, because axoaxonic
  predictions onto sensory axons have a high false-positive rate (Figure S1a–c). Only the
  two neurons per subtype that were fully hand-validated are quantified.

Panel 3e is the one panel here that reads a shipped table rather than querying CAVE: the
partner table carries a per-neuron `num_supervoxels` value used for the "intact" filter,
which is not obtainable from the synapse table alone.
`data_generation/build_sensory_postsynaptic_partner_table.ipynb` rebuilds that table from
CAVE (see below).

#### How Figure 3e counts "downstream neurons"

Each row of `data/sensory_postsynaptic_partners.csv` is one (sensory neuron, postsynaptic
object) pair. The panel:

1. keeps rows whose presynaptic neuron is one of the 79 reconstructed sensory neurons;
2. keeps rows with **`num_supervoxels` ≥ 40,000** — the "intact" filter;
3. counts **distinct `post_pt_root_id`** per sensory neuron;
4. assigns 0 to any sensory neuron left with no intact partner (it stays a data point).

There is deliberately **no minimum synapse threshold** — a single synapse makes a partner
count. Among retained pairs 42% are single-synapse contacts (mean 4.84 synapses per pair).

The intact filter does most of the work. Postsynaptic object size is sharply bimodal: 76%
of the 33,077 distinct objects have fewer than 100 supervoxels (unproofread fragments),
while a second population sits above 40,000. The threshold falls just past the trough and
retains 3,115 objects (9.4%), removing 72% of raw partner rows. Without it the median
sensory neuron would appear to contact 451 partners rather than 141, and the subtype
ordering would invert — Aδ-LTMR would outrank Aβ-LTMR purely by collecting more debris.

Absolute counts scale with the threshold but the **ordering is stable** from 5,000 to
100,000 supervoxels, so the conclusion (A-fibre LTMRs contact ~300 spinal neurons, far more
than C-HTMR/Heat, peptidergic and C-Cold subtypes) does not depend on the exact cut.

#### Rebuilding the partner table

`data_generation/build_sensory_postsynaptic_partner_table.ipynb`

The script that originally produced this table (February 2026) was not kept. The recipe was
reverse-engineered from the file and **verified against live CAVE**:

- `num_supervoxels` = `len(client.chunkedgraph.get_leaves(post_pt_root_id))` — reproduced
  the stored value exactly for 8/8 spot-checked partners;
- `num_synapses` = synapses from `synapses_version3` with the sensory neuron presynaptic,
  autapses dropped, cleft `size` ≥ 20 (identical filtering to panel 3d) — reproduced the
  stored counts for 100% of the 864 overlapping partners of a test sensory neuron;
- a 25-row end-to-end rebuild matched the shipped table exactly on both columns.

A full rebuild is ~33,000 ChunkedGraph calls at roughly 6/s (~1.5 h). Counts are cached to
disk and the notebook is resumable.

#### Root-ID drift

3d and 3f pin `MATERIALIZATION = 669` and pass every published root ID through
`resolve_current_ids()`, which substitutes a current successor for any ID made stale by
later proofreading (see the [top-level README](../README.md#root-ids-drift--pin-the-materialization)).
Three of the 79 sensory neurons currently need substituting:

| Published ID | Subtype | Current ID | Saved count | Count now |
| --- | --- | --- | --- | --- |
| 720575940922958792 | C-LTMR | 720575940859731351 | 2634 | 2659 |
| 720575940902363320 | Aδ-LTMR | 720575940955315288 | 2809 | 2821 |
| 720575940893392149 | Aδ-LTMR | 720575940897209813 | 4112 | 4112 |

With this in place, panel 3d reproduces the published subtype means (6 of 8 exactly; C-LTMR
2534.8 vs 2532.3 and Aδ-LTMR 2663.9 vs 2662.7 differ by <0.1% because of proofreading on
those three neurons) and **all three legend *p*-values exactly**
(4.07 × 10⁻⁴, 8.03 × 10⁻⁴, 4.52 × 10⁻²).

Panel 3e deliberately does **not** resolve IDs: it joins against a snapshot
(`sensory_postsynaptic_partners.csv`) that is keyed by the published IDs, so the two must
stay consistent.

| `figure3hi_mesh_synapse_density.ipynb` | **3h, 3i** | Axon mesh colour-coded by local synapse density: geodesic (Dijkstra-on-vertex-graph) neighbourhood of **R = 4 µm**, axodendritic in **h** and axoaxonic in **i**, on one shared `viridis` scale capped at 15. Exports colour-baked `.ply` meshes plus a colourbar. **CAVE required** — the mesh is fetched via CloudVolume and is far too large to ship. |
| `figure3jk_skeleton_synapse_maps.ipynb` | **3j, 3k** | Flat (x–y) skeleton trace per axon with synapses overlaid: axodendritic in **j**, axoaxonic in **k**. The 16 validated axons. Skeletons from `data/sensory_neuron_swc/`; synapse positions queried live, so **CAVE is required**. |

Panels 3j/3k are plotted in pixel coordinates aligned on the soma so traces are comparable
across neurons, and panel j runs a two-pass version that fixes common x/y limits across all
16 axons first. The synapse table and the SWC files use different voxel resolutions
(16 × 16 × 45 nm and 32 × 32 × 45 nm); the plotting functions convert between them, so
don't change one without the other.

#### `data_generation/export_sensory_neuron_swc.ipynb`

Skeletonizes each sensory axon with `pcg_skel.pcg_skeleton` and writes
`data/sensory_neuron_swc/<tag>_<root_id>.swc` (16 files, 968 KB, shipped). Only needed to
rebuild the skeletons; requires CAVE access.

`data_generation/remove_false_positive_axoaxonic.ipynb` records the manual removal of
false-positive axoaxonic predictions; the curated result lives in the CAVE
`sensory_neuron_synapses` table.

Every Figure 3 panel with analysis code is now covered (d, e, f, g, h, i, j, k).

### Per-glomerulus quantification in the dSC_APEX volumes (CATMAID)

`figure3_glomerulus_quantification_catmaid/`

| Notebook | Panel | Data source | Statistics |
| --- | --- | --- | --- |
| `figure3g_axoaxonic_synapses_per_glomerulus.ipynb` | **3g** — axoaxonic synapses per glomerulus, 10 subtypes | `data/catmaid_glomeruli/keyword_counts_*_*.json` (ships with the repo — no CATMAID access needed) | Kruskal–Wallis + Dunn (Holm) |

Violins are ranked by uncapped mean and capped at 15 for display (top tick reads `15+`);
the statistics use the uncapped values. The leading number in each JSON filename fixes the
published subtype order.

This notebook runs offline. Verified against the published legend: all ten glomerulus
counts (*n* = 42, 102, 79, 45, 118, 44, 51, 45, 34, 42) and all seven reported *p*-values
reproduce exactly.

#### `data_generation/` — how the JSON counts were produced

Only needed to rebuild the raw data from scratch; requires CATMAID access and the `pymaid`
environment. Skeletons were traced manually in CATMAID, one project per genetically
labeled dataset, then synapses grouped into glomeruli with an iterative 2 µm distance
threshold.

| Subtype | CATMAID project | Trace notebook | Glomerulus clustering |
| --- | --- | --- | --- |
| Aβ-RA-LTMR (Ret⁺) | 12 | `trace_abltmr.ipynb` | `cluster_abltmr.ipynb` |
| Aδ-LTMR (TrkB⁺) | 8, 9 | `trace_adltmr_trkb.ipynb` | `cluster_adltmr.ipynb` |
| C-LTMR | 6 | `trace_cltmr.ipynb` | `cluster_cltmr.ipynb` |
| C-HTMR/Heat (MrgprD⁺) | 13 | `trace_mrgd.ipynb` | `cluster_mrgd.ipynb` |
| C-HTMR/Heat (MrgprA3⁺) | 9 | `trace_mrgpra3.ipynb` | `cluster_mrgpra3.ipynb` |
| C-HTMR/Heat (Sst⁺, Cysltr2-labeled) | 6 | `trace_sst_cysltr2.ipynb` | `cluster_sst_cysltr2.ipynb` |
| Aδ-HTMR (Smr2⁺, dorsal) | 7 | `trace_smr2_dorsal.ipynb` | `cluster_smr2_dorsal.ipynb` |
| Aδ-HTMR (Smr2⁺, ventral) | 7 | `trace_smr2_ventral.ipynb` | `cluster_smr2_ventral.ipynb` |
| C-Heat (Sstr2⁺) | 15 | `trace_sstr2.ipynb` | `cluster_sstr2.ipynb` |
| C-Cold (Trpm8⁺) | 14 | `trace_trpm8.ipynb` | `cluster_trpm8.ipynb` |

Each `cluster_*` notebook writes one
`data/catmaid_glomeruli/keyword_counts_<n>_<subtype>.json`, whose first element is the
per-glomerulus axoaxonic synapse count array read by the panel 3g notebook.

---

## Figure S6 — Distribution of axodendritic and axoaxonic synapses

`figureS6_synapse_distribution/` — panels a–j. **All except panel a run offline**; the
meshwork `.h5` files and CATMAID `.npz` files ship with the repository.

| Notebook | Panels | Quantity | Input |
| --- | --- | --- | --- |
| `figureS6a_mesh_synapse_density.ipynb` | **a** | Example meshes coloured by synapse density | same geodesic R = 4 µm calculation as Fig. 3h,i, but with two single-hue ramps — axoaxonic PuRd → `#cb2f43`, axodendritic PuBu → `#1c4286` — instead of a shared viridis scale. **CAVE required.** |
| `figureS6b_axodendritic_distance_to_axon_tip.ipynb` | **b** | Path distance from each **axodendritic** synapse to the nearest axon tip, + grey cable-length reference | `data/meshwork_h5/` (80 axons, 10/subtype) |
| `figureS6c_axoaxonic_distance_to_axon_tip.ipynb` | **c** | Same for **axoaxonic** synapses, 5 subtypes | 16 validated axons (2/subtype) |
| `figureS6d_h_distance_to_nearest_axoaxonic.ipynb` | **d–h** | Path distance from each axodendritic synapse to its nearest **axoaxonic** synapse | 16 validated axons |
| `figureS6i_axoaxonic_within_4um_catmaid.ipynb` | **i** | Axoaxonic synapses within 4 µm of each axodendritic synapse, dSC_APEX volumes | `data/catmaid_psi_within_4um/*.npz` |
| `figureS6j_adhtmr_dorsal_vs_ventral.ipynb` | **j** | Paired dorsal vs ventral comparison, 8 Aδ-HTMRs | `data/meshwork_h5/adhtmr/` |

### Aggregation (panels b, c)

One CDF is computed **per neuron** and interpolated onto a shared grid; the plotted curve
is the **mean across neurons** of a subtype, not a pool of all synapses. This weights each
neuron equally regardless of synapse count. The grey reference curve weights each skeleton
edge by its physical cable length, giving the fraction of arbor within X µm of a tip — the
null expectation if synapses were spread evenly.

### Panel h / j split

Aδ-HTMR axons span laminae I–IV, so their axodendritic synapses are split into dorsal and
ventral by the mesh y coordinate at `Y_SPLIT = 140,800 nm`.

### Verified against the paper

- **b** — LTMRs 72.4 / 80.3 / 75.1% of axodendritic synapses within 25 µm of a tip
  (Results: "~75%"), against only **21.1%** of axon cable (Results: "~20%"). Aδ-HTMR,
  peptidergic and C-Cold reach ~75% only by 50 µm (70.8 / 69.2 / 78.3%), as stated.
- **d–h** — ≥5 µm co-localization: Aβ-LTMR 91.5%, Aδ-LTMR 97.7%, C-LTMR 86.5%
  (Results: "virtually all"); MrgprD⁺ C-HTMR/Heat 74.1%.
- **i** — all three legend *p*-values reproduce exactly (6.38 × 10⁻⁶⁰, 5.41 × 10⁻¹², 1).
- **j** — p = 0.01172, exactly the published 1.17 × 10⁻².

### Two statistical notes

- **Panel i uses Bonferroni correction**, whereas Figures 3d/3e/3g and the rest of this
  figure use Holm. Bonferroni is what reproduces the published values; Holm would give
  5.25 × 10⁻⁶⁰, 2.04 × 10⁻¹² and 0.714. The legend says only "Dunn's post hoc test".
- **Panel j's published *p* is one-sided** (two-sided p = 0.02344, one-sided = 0.01172).
  The notebook prints both.

### `data_generation/`

`build_meshwork_h5.ipynb` rebuilds the `.h5` exports from CAVE with
`pcg_skel.pcg_meshwork(..., synapse_table='sensory_neuron_synapses')`, and regenerates
`data/meshwork_h5_subtypes.csv` (the file → subtype map). Only needed to rebuild the
shipped data; requires CAVE access.

---

## Figure 4 — Homotypic presynaptic inhibition (PSI)

`figure4_homotypic_psi/`

Panels d–g. Both notebooks **run offline** — the SPINE feather tables ship with the
repository. Panels a–c are schematics and skeleton renderings, no analysis code.

| Notebook | Panels | Output |
| --- | --- | --- |
| `figure4de_psi_ihn_sensory_heatmaps.ipynb` | **d, e** | Per-neuron heatmaps: rows = the 34 PSI-IhNs grouped into 5 populations, columns = sensory subtypes, colour = that subtype's fraction of the neuron's sensory synapses |
| `figure4fg_psi_ihn_summed_fractions.ipynb` | **f, g** | Stacked bars: the same fractions pooled per population |

Panels **e/f** are *dendritic inputs* (sensory synapses onto the PSI-IhN); panels **d/g**
are *axoaxonic outputs* (PSI-IhN synapses onto sensory axons). The two directions use
identical code on two feather sets.

### Input data

`data/psi_ihn_feather/` — per-neuron SPINE prediction tables, one row per synapse with the
eight sensory subtype probabilities and a called `cell_type`. These are outputs of the SPINE
classifier, which is maintained separately at
[htem/dorsalhorn-ml](https://github.com/htem/dorsalhorn-ml):

| Path | Contents |
| --- | --- |
| `dendritic_inputs/` | 34 `.feather`, synapses onto each PSI-IhN dendrite (58,729 synapses) |
| `axoaxonic_outputs/` | 34 `.feather`, each PSI-IhN's axoaxonic synapses onto sensory axons (4,839 synapses) |
| `psi_ihn_subtypes.csv` | coordinate → root ID → population tag |

The 34 neurons split exactly as reported in the Methods: 9 PSI-IhN^Aδ-HTMR,
6 Islet-PSI-IhN^C-HTMR/Heat, 6 Islet-PSI-IhN^C-LTMR, 6 Islet-PSI-IhN^Aδ-LTMR,
7 PSI-IhN^Aβ-LTMR.

### Two details that matter

- **Low-confidence relabelling.** A synapse SPINE calls sensory but whose every subtype
  probability is ≤ 0.6 is reassigned to `low_confidence_sensory` (Methods, *Inhibitory
  neuron connectivity*) — it becomes its own column rather than being dropped or forced
  into a subtype.
- **Panels f/g pool, they do not average.** Counts are summed across the neurons of a
  population before normalising (`groupby(population).sum()`), so larger neurons contribute
  proportionally more. The per-neuron view is panels d/e.

Segments within each bar run in descending order, except the Aβ-LTMR row where the three
LTMR subtypes are pinned to the front (Aβ, Aδ, C) so the LTMR block reads together.

### Reproduces

The homotypic diagonal, mean per-neuron fraction from the neuron's own sensory subtype:

| Population | Input | Output |
| --- | --- | --- |
| PSI-IhN^Aδ-HTMR | 76.0% | 81.3% |
| Islet-PSI-IhN^C-HTMR/Heat | 91.6% | 92.2% |
| Islet-PSI-IhN^C-LTMR | 94.7% | 94.8% |
| Islet-PSI-IhN^Aδ-LTMR | 65.8% | 89.9% |
| PSI-IhN^Aβ-LTMR | 49.6% | 68.8% |

Pooled, Islet-PSI-IhN^C-LTMR receives **94.7%** of its sensory input from C-LTMRs, matching
the Results (">90%").

### `data_generation/`

`build_psi_ihn_subtype_map.ipynb` rebuilds `psi_ihn_subtypes.csv`: the feather files are
named by seed coordinate, so it reads the segmentation at each voxel (CloudVolume) to get a
root ID and joins to the CAVE `excitatory_neuron_label` table. Requires CAVE access; only
needed if the feather set changes.

---

## Figure S9 — Lateral inhibition mediated by PSI-IhNs

`figureS9_lateral_inhibition/figureS9cg_lateral_inhibition.ipynb` — panels **c–g**
(a and b are schematics). **Runs offline — no CAVE access.**

For each **PSI-IhN ↔ sensory neuron** pair the notebook computes an *incoming* fraction
(this afferent's share of the PSI-IhN's dendritic sensory input) and an *outgoing* fraction
(its share of the PSI-IhN's axoaxonic output), and joins them with a line. Pairs are
anchored on the axoaxonic connection — only pairs with a real PSI-IhN axon → sensory
synapse are included, and an afferent that makes no synapse back onto the dendrite gets an
incoming fraction of 0 rather than being dropped (Methods, *Lateral inhibition analysis*).

### Input data

`data/lateral_inhibition/` — four connectivity tables extracted from CAVE beforehand, which
is why the notebook needs no network access:

| File | Rows | Contents |
| --- | --- | --- |
| `nodes.csv` | 425 | reconstructed neurons → subtype tag and position |
| `edges_all.csv` | 10,225 | synapse counts between reconstructed neurons (gives the sensory → PSI-IhN dendrite direction) |
| `edges_axoaxonic.csv` | 2,964 | PSI-IhN axon → sensory axon counts |
| `psi_ihn_axon_dendrite.csv` | 28 | per PSI-IhN: axon and dendrite root IDs, and total sensory input onto the dendrite |

### Pair counts

| Panel | Pairing | Notebook | Published legend |
| --- | --- | --- | --- |
| c | Islet-PSI-IhN^C-LTMR / C-LTMR | 52 | 52 ✓ |
| d | Islet-PSI-IhN^C-HTMR/Heat / C-HTMR/Heat | **55** | 60 |
| e | Islet-PSI-IhN^Aδ-LTMR / Aδ-LTMR | **60** | 55 |
| f | PSI-IhN^Aδ-HTMR / Aδ-HTMR | 49 | 49 ✓ |
| g | PSI-IhN^Aβ-LTMR / Aβ-LTMR | 31 | 31 ✓ |

Panels d and e appear **transposed in the legend**. Counting the data lines in the
published SVGs (`lateral_in/homotypic_*_vs_*.svg`, where each pair draws exactly three
alpha-0.4 elements) gives 55 for the C-HTMR/Heat panel and 60 for the Aδ-LTMR panel —
matching the notebook, not the legend. The artwork is right; the two numbers in the caption
need swapping.

### An outlier in panel g

One pair reaches 70.4% incoming, which looks inconsistent with the Results. It is a
denominator effect: one `islet-abltmr` PSI-IhN has only **27** sensory input synapses on
its dendrite (vs 225–454 for the rest of that population, ~1,000–1,500 for the islet
C-LTMR and C-HTMR/Heat cells), so a single 19-synapse afferent is 70% of it. A second
neuron in that population has 66. Both are likely incompletely reconstructed dendrites. The
notebook prints the per-neuron denominators; nothing is filtered out.

---

## Figure 5 — Heterotypic feedforward inhibition

`figure5_feedforward_inhibition/` — panels b–e. Panel a is a schematic.

| Notebook | Panels | Needs CAVE? |
| --- | --- | --- |
| `figure5b_synapse_type_fractions.ipynb` | **b** — axoaxonic / axodendritic / dendritic output fractions per PSI-IhN | **yes** (live synapse counts) |
| `figure5c_homotypic_selectivity_ecdf.ipynb` | **c** — cumulative homotypic selectivity of postsynaptic neurons | no |
| `figure5de_postsynaptic_sensory_composition.ipynb` | **d, e** — per-neuron heatmap and pooled bars of postsynaptic sensory input | no |
| `figure5fj_S11_postsynaptic_sensory_modules_axon.ipynb` | **f–j** — sensory module composition of each population's postsynaptic neurons, **plus Fig. S11a–j** | no |
| `figureS12_postsynaptic_sensory_modules_dendrite.ipynb` | Fig. S12c–h — the same analysis for dendrodendritic (dendrite-presynaptic) targets | no |

### Panels f–j come from the supplement notebook

`figure5fj_S11_postsynaptic_sensory_modules_axon.ipynb` produces main-figure and
supplementary panels together, because they share one pass over the feather tables:

| Output | Panel |
| --- | --- |
| `violin_<tag>__axon_pre.svg` | Fig. S11a–e |
| `<tag>__selectivity__axon_pre.svg` — selective vs non-selective | Fig. S11f–j |
| `<tag>__composition__axon_pre.svg` — **module breakdown among selective neurons** | **Fig. 5f–j** |
| `ecdf_target_sensory_fraction__axon_pre.svg` | same quantity as Fig. 5c |
| `post_syn_sensory_composition__axon_pre_0910.svg` | same quantity as Fig. 5e |

A neuron belongs to a **sensory module** when one subtype supplies >50% of its sensory
input. Panels f–j show only those neurons, so their *n* is lower than S11a–e:

| Population | S11a–e *n* | 5f–j *n* |
| --- | --- | --- |
| PSI-IhN^Aδ-HTMR | 392 | 254 |
| Islet-PSI-IhN^C-HTMR/Heat | 159 | 121 |
| Islet-PSI-IhN^C-LTMR | 85 | 66 |
| Islet-PSI-IhN^Aδ-LTMR | 519 | 426 |
| PSI-IhN^Aβ-LTMR | 554 | 420 |
| **total** | **1709** | **1287** |

Both columns reproduce their legends exactly.

#### Where 1,553 and 77.1% come from

The Results quote "77.1% of the 1,553 dorsal horn neurons that are postsynaptic to
PSI-IhNs". That is a **per-neuron** statistic over **both** postsynaptic routes, so its
denominator is different from the per-connection, axon-only counts in the table above.
Both reconcile exactly:

| | count |
| --- | --- |
| axon-presynaptic neurons with a SPINE feather | 1435 |
| dendrite-only neurons with a feather (recovered partners, `_v2` list) | 118 |
| **distinct postsynaptic neurons** = feather files shipped | **1553** |

Of those 1,553, **1198 are selective** — one sensory subtype supplies ≥50% of their sensory
input — which is **77.1%**, matching the Results exactly. (The threshold is `>=` 50%; a
strict `>` gives 1194, or 76.9%.)

So there are three consistent numbers, each a different unit:

| Number | Unit | Where used |
| --- | --- | --- |
| **1553** neurons | distinct neurons, axodendritic ∪ dendrodendritic | Methods, Results |
| **1709** connections | PSI-IhN → neuron pairs, axodendritic only | Fig. 5c–e, S11a–e legends |
| **1287** connections | of those, selective | Fig. 5f–j legends |

A neuron contacted by several PSI-IhNs is one entry in the first row but several in the
second, which is why 1709 > 1553.

### The counting unit (panels c, d, e)

A spinal cord neuron contacted by more than one PSI-IhN population is counted under
**each** of them — per the Methods, "a single spinal cord neuron could belong to multiple
postsynaptic neuron groups". The unit is therefore the PSI-IhN → postsynaptic-neuron
**connection**, not the neuron. Deduplicating instead gives noticeably smaller groups and
does not reproduce the published *n* values.

From 1,786 axon-presynaptic connections (1,499 distinct neurons), 77 are dropped because
the neuron's feather file is missing or holds no classifiable sensory synapse, leaving
**1,709** — matching the legend exactly:

| Population | *n* connections |
| --- | --- |
| PSI-IhN^Aδ-HTMR | 392 |
| Islet-PSI-IhN^C-HTMR/Heat | 159 |
| Islet-PSI-IhN^C-LTMR | 85 |
| Islet-PSI-IhN^Aδ-LTMR | 519 |
| PSI-IhN^Aβ-LTMR | 554 |

### Reproduces

- **c** — median homotypic selectivity 0.26 / 0.72 / 0.50 / 0.34 / 0.46 across the five
  populations. Islet-PSI-IhN^C-HTMR/Heat is the outlier (75% of its partners are >50%
  homotypic), matching "most postsynaptic neurons show low homotypic selectivity, except
  those postsynaptic to Islet-PSI-IhN^C-HTMR/Heat".
- **e** — partners of Islet-PSI-IhN^C-LTMR draw **33.9%** of their sensory input from
  C-HTMR/Heat, the exact value quoted in the Results. Partners of PSI-IhN^Aδ-HTMR take more
  input from C-HTMR/Heat (34.1%) than from Aδ-HTMR (28.7%), as stated.

### Panel b needs a materialization **newer than 669**

Panel b splits each PSI-IhN's output by compartment, which only works where the axon and
dendrite are genuinely separate segments.

In materialization **669** — the newest available when this repository was assembled —
2 of the 34 neurons (both `islet-cltmr`, seeds `22476_4688_1997` and `24266_4298_1651`)
were still merged into one segment. Their "dendritic" count was therefore the whole
neuron's output counted a second time: 3427 + 149 = 3576 and 2775 + 155 = 2930, exactly.
That pins both neurons near 0.50/0.50 and drags the `islet-cltmr` mean with them.

**Both neurons have since been split manually**, so any materialization later than 669
should be correct. To use one:

1. run `data_generation/build_psi_ihn_compartment_map.ipynb` against the new version — it
   rewrites `data/psi_ihn_compartments.csv` with the new root IDs;
2. set `MATERIALIZATION` in the panel b notebook to that version;
3. check it reports **0 unsplit neurons**. The "split only" cell then becomes a no-op.

On version 669 the notebook plots the panel both ways:

| `islet-cltmr` | axoaxonic | axodendritic | dendritic |
| --- | --- | --- | --- |
| all 6 neurons (v669) | 0.072 | 0.246 | 0.682 |
| **4 split neurons only** | **0.096** | **0.131** | **0.773** |
| `islet-chtmr`, for comparison | 0.096 | 0.136 | 0.768 |

The split-only version puts `islet-cltmr` essentially on top of `islet-chtmr`, which is
what two islet populations should look like — so it is the better estimate until a newer
materialization is available.

### Data

| Path | Contents |
| --- | --- |
| `data/psi_ihn_postsynaptic_feather/` | 1,553 per-neuron SPINE tables, `cell_type_predictions_<root_id>.feather` (129 MB) — the postsynaptic neurons analysed in panels c–e |
| `data/psi_ihn_partners_filtered_gt400.csv` | PSI-IhN → postsynaptic partner edges after the ≥400-input and sensory-like filters; `kind` is `axon_pre` or `dendrite_pre`. **This is the version the figure notebooks read.** |
| `data/psi_ihn_partners_filtered_gt400_v2.csv` | Same, with additional recovered dendrite-presynaptic partners (630 vs 424 `dendrite_pre` rows; `axon_pre` identical). Used by the dendrite supplement. |
| `data/psi_ihn_compartments.csv` | Per PSI-IhN: axon and dendrite seed coordinates, root IDs **as of materialization 669**, and whether the two are separate segments. Regenerate with `data_generation/build_psi_ihn_compartment_map.ipynb` when moving to a newer version. |

---

## Data

`data/`

| Path | Contents |
| --- | --- |
| `con_graph17492931.pkl` | Pickled `networkx` connectivity graph snapshot (CAVE materialization 17492931). Nodes carry `type` (cell-type code) and `coords`; edges carry `weight` = synapse count. |
| `sensory_postsynaptic_partners.csv` | Per-sensory-neuron postsynaptic partner list (Figure 3e). |
| `psi_ihn_partners_filtered_gt400.csv`, `..._v2.csv` | PSI-IhN → postsynaptic partner edges (Figure 5c–e); see the Figure 5 section. |
| `psi_ihn_postsynaptic_feather/` | 1,553 per-neuron SPINE prediction tables for the postsynaptic neurons (129 MB). |
| `psi_ihn_compartments.csv` | PSI-IhN axon / dendrite compartment coordinates and root IDs (Figure 5b). |
| `psi_ihn_axon_partners.csv`, `psi_ihn_dendrite_partners.csv` | Raw partner lists before filtering. |
| `missing_dendrite_partners_*.csv` | Root-ID recovery intermediates. |
| `non_sensory_gt400_root_ids.csv` | Root IDs passing the ≥400-input filter. |
| `homotypic_lateral_inhibition_pairs.csv`, `psi_ihn_dendrite_input.csv` | Lateral inhibition pair fractions. |
| `overlap_ihn_axon_to_sensory.csv` | Axon–dendrite path overlap, IhN axon → sensory axon. |
| `triadic_motif_counts_by_type.csv`, `triadic_probabilities_SxIxE.csv` | Triadic motif counts and probabilities. |
| `ihn_connectivity/` | Per-compartment adjacency, overlap and density CSVs (`csv_axoaxonic/`, `csv_axodendritic/`, `csv_dendrodendritic/`), plus `triad_null_FINE_1000.csv` (1,000 null resamples). |
| `catmaid_glomeruli/` | `keyword_counts_*.json` per-glomerulus counts, `distance_matrix.csv`, `pairwise_dunn_axoaxonic.csv`. |
| `psi_ihn_feather/` | Per-neuron SPINE prediction tables for the 34 PSI-IhNs (Figure 4d–g): `dendritic_inputs/` and `axoaxonic_outputs/`, 34 `.feather` each, plus `psi_ihn_subtypes.csv`. 3.9 MB. |
| `lateral_inhibition/` | Four pre-extracted connectivity tables for Figure S9c–g (544 KB). |
| `catmaid_psi_within_4um/` | `*_pre-vaxon_final.npz`, each with a `psi_close` array of per-axodendritic-synapse axoaxonic counts (Figure S6i). `trpm8_deep` duplicates `trpm8` and is unused. |
| `meshwork_h5/` | 80 `pcg_skel` meshwork exports (`x_y_z.h5`, 42 MB) — level-2 skeleton, mesh and `sensory_neuron_synapses` per axon. Plus `adhtmr/` (8 files) for Figure S6j. |
| `meshwork_h5_subtypes.csv` | Maps each `.h5` to its sensory subtype, with per-file synapse counts. |
| `sensory_neuron_swc/` | SWC skeletons of the 16 validated sensory axons, `<tag>_<root_id>.swc` (Figure 3j,k). |

