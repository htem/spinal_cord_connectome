# Environments

Two conda environments were used. Both are exact `conda env export --no-builds` captures of
the environments the analyses were actually run in, so package versions are real rather than
reconstructed.

| File | Env name | Python | Use for |
| --- | --- | --- | --- |
| `cave310.yml` | `cave310` | 3.10.18 | Everything CAVE-based: connectivity matrices, path overlap, UMAP/HDBSCAN clustering, mesh synapse density, figure plotting. |
| `pymaid.yml` | `pymaid` | 3.7.16 | CATMAID tracing and glomerulus quantification (`pymaid`, `neuroboom`), the `navis` network traversal model, and `pymeshfix` soma volume repair. |

```bash
conda env create -f envs/cave310.yml
conda env create -f envs/pymaid.yml
```

## Which environment for which notebook

Most notebooks run under **`cave310`**. Use **`pymaid`** for:

- `manuscript1_sensory_synaptic_organization/figure3_glomerulus_quantification_catmaid/` —
  all of it (`pymaid`, `neuroboom` are only in this env).
- `manuscript2_wiring_diagram/figure4_wiring_diagram/figure4e_S7h_network_traversal.ipynb`
  — `navis==1.4.0` (and therefore `navis.models.TraversalModel`) is only in `pymaid`.
- `manuscript2_wiring_diagram/figure3_pn_clustering/data_generation/quantify_soma_volume_15um.ipynb`
  and `quantify_soma_volume_20um_first_pass.ipynb` — `pymeshfix==0.16.2` is only in `pymaid`.

Note that the two environments carry different CAVE stack versions
(`caveclient` 7.11.0 / 5.11.0, `cloud-volume` 12.5.0 / 8.27.0, `pcg-skel` 1.3.1 / 0.3.5).
`pymaid` is the older environment used earlier in the project; `cave310` is current.

## Key versions (`cave310`)

| Package | Version | | Package | Version |
| --- | --- | --- | --- | --- |
| python | 3.10.18 | | caveclient | 7.11.0 |
| numpy | 2.2.6 | | cloud-volume | 12.5.0 |
| scipy | 1.15.3 | | pcg-skel | 1.3.1 |
| pandas | 2.3.2 | | meshparty | 2.0.3 |
| matplotlib | 3.10.6 | | skeleton-plot | 0.1.0 |
| seaborn | 0.13.2 | | trimesh | 4.8.2 |
| networkx | 3.4.2 | | open3d | 0.19.0 |
| scikit-learn | 1.7.2 | | umap-learn | 0.5.9.post2 |
| statsmodels | 0.14.5 | | hdbscan | 0.8.40 |
| scikit-posthocs | 0.11.4 | | pyarrow | 19.0.0 |

## Key versions (`pymaid`)

| Package | Version |
| --- | --- |
| python | 3.7.16 |
| python-catmaid (`pymaid`) | 2.4.0 |
| neuroboom | 0.3.55 |
| navis | 1.4.0 |
| pymeshfix | 0.16.2 |
| caveclient | 5.11.0 |
| numpy | 1.21.5 |
