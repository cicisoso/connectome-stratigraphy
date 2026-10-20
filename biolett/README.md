# Neuronal birth order in the geometry of three *Drosophila* connectomes: analysis code

This is the code for Liu H, Zong X, Chen J, Xiong J, "Neuronal birth order leaves a reproducible trace in the geometry of three *Drosophila* connectomes" (submitted to *Biology Letters*).

It tests the site-of-addition rule: a neuron's age grows with its distance from the site where its lineage adds new cells. That site is the neuroblast for somata and the core of the lineage tract for fibres.

The tests use four connectomes:

- MaleCNS v1.0
- FlyWire v783
- BANC v888
- MANC v1.0 (connectivity only)

`run_all.sh` regenerates every result, figure and table of the paper and its electronic supplementary material from the public data.

## Requirements

Python 3.11 with:

| Package | Version |
|---|---|
| pandas | 2.3 |
| NumPy | 2.2 |
| SciPy | 1.16 |
| scikit-learn | 1.7 |
| PyArrow | 22 |
| Matplotlib | 3.10 |
| openpyxl | 3.1 |
| requests | 2.34 (only for `tools/`) |

`pip install -r requirements.txt` installs the exact versions used. Most of the run time is spent reading skeletons.

The figure scripts also run a panel-alignment check if the `audit_panel_alignment` module can be imported (from the path or from `$NATURE_FIGURE_SCRIPTS`). Without the module they skip the check; the figures are the same.

## Data

Set the environment variable `FLYBRAIN_DATA` to a folder with the layout below. The default is `data/` in this folder.

| File or folder under `FLYBRAIN_DATA` | Source |
|---|---|
| `body-annotations-male-cns-v1.0-minconf-0.5.feather`, `connectome-weights-male-cns-v1.0-minconf-0.5.feather` | MaleCNS v1.0 flat connectome, https://storage.googleapis.com/flyem-male-cns/v1.0/connectome-data/flat-connectome/ (https://male-cns.janelia.org) |
| `swc/<bodyId>.swc` | MaleCNS v1.0 skeletons, https://storage.googleapis.com/flyem-male-cns/v1.0/segmentation/skeletons-malecns/skeletons-swc/ (fetch with `tools/fetch_malecns_swc.py`) |
| `flywire/annotations/Supplemental_file1_neuron_annotations.tsv` | FlyWire annotations (Schlegel et al. 2024), https://github.com/flyconnectome/flywire_annotations |
| `flywire/connections.csv.gz`, `flywire/cell_stats.csv.gz` | FlyWire v783 downloads, https://codex.flywire.ai |
| `flywire/skeletons/sk_lod1_783_healed_ds2.parquet` | FlyWire v783 skeletons, https://zenodo.org/records/10877326 (5.4 GB; `tools/pdownload.py` downloads it in parallel) |
| `banc/banc_888_meta.feather`, `banc/banc_888_edgelist_simple_v3.feather` | BANC v888 compiled data (Bates et al. 2026), Harvard Dataverse https://doi.org/10.7910/DVN/7WTH1N |
| `banc/skel/swcs-from-pcg-skel/<root_id>.swc` | BANC v888 skeletons (`neuron_skeletons.zip` of the same release) |
| `manc/manc-v1.0-neuron-properties.feather`, `manc/traced-connections.csv` | MANC v1.0, https://storage.googleapis.com/flyem-manc-exports/v1.0 (https://neuprint.janelia.org) |

`external/` holds supplementary files 1 and 2 of Lee et al. 2020 (*eLife* 9, e53518; CC BY 4.0). These give the clone-induction windows used as birth order for the vPN and PAM lineages.

## Running

```
cd code
zsh ../run_all.sh
```

| Script | What it does |
|---|---|
| `common.py` | Loaders for the four connectomes, skeleton access, geometry (cell body fibre, tract cross-section, bundle measures, geometric age index) and statistics. |
| `a0_prepare.py` | Extracts the FlyWire skeletons that are needed into a compact cache. |
| `a1_kc.py` | Kenyon cells: peduncle and α/β-lobe fibre position, and soma distance from the clone axis, by class. |
| `a2_lineages.py` | Clonally dated lineages (adPN, lPN, vPN, PAM): soma depth, soma distance, fibre radial position and the geometric age index against birth rank. |
| `a3_hemilineages.py` | Hemilineage units: primary versus secondary, and early versus late secondary. |
| `a4_connectivity.py` | Synapse counts, size and partner classes of primary and secondary neurons. |
| `a5_summary.py` | Pooled estimates, concordance, sign tests and the summary tables in `results/`. |
| `a6_render_data.py` | Axes, centroids and sections for the renderings. |
| `fig1.py`, `fig2.py`, `fig3.py` | Main figures. |
| `figS_renderings.py`, `figS_more.py` | Supplementary figures S1–S6. |
| `make_tables.py` | Table 1, the supplementary workbook (tables S1–S11) and the per-neuron CSV files. |

Outputs go to `results/`, `figures/` and `tables/`. `a0_prepare.py` stores the FlyWire skeletons it needs in `results/flywire_skeleton_cache.parquet` (about 450 MB). The GitHub repository holds all these outputs except that cache and the TIFF and SVG copies of the figures.

## Archive

This code is release v2.0.0 of https://github.com/cicisoso/connectome-stratigraphy (folder `biolett/`), archived at Zenodo (all versions: https://doi.org/10.5281/zenodo.22760055).

## Licence

The code is released under the MIT licence. The connectome data are distributed by their creators under their own terms.
