# Connectome stratigraphy: neuronal birth order is written in the geometry of a static wiring diagram

Code, per-neuron tables and figures for the manuscript *"Connectome stratigraphy: neuronal birth order is written in the geometry of a static wiring diagram"* (Liu, Zong, Chen and Xiong; School of Information Technology, Zhejiang Financial College; submitted to *Development* as a Research Report).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22760056.svg)](https://doi.org/10.5281/zenodo.22760056)

We ask whether the order in which neurons were born can be read from the geometry of an adult connectome. In the MaleCNS v1.0 electron-microscopy connectome of a male *Drosophila* central nervous system we measure, for every neuron, the radial position and burial depth of its fibre within its lineage bundle and the position of its soma within its clone. These measures order the seven Kenyon cell classes by birth (youngest fibres in the core of the α/β lobes), recover the twin-spot MARCM birth orders of the anterodorsal and lateral projection neuron lineages, and separate embryonic-born from larval-born neurons in 155 hemilineage units of the nerve cord. Embryonic-born neurons make about four times more synapses than their larval-born lineage-mates and are 95% of the most connected neurons.

## Contents

| Folder | What it holds |
|---|---|
| `code/` | Analysis and figure scripts (Python 3.11) and `download_data.py`. |
| `results/` | Outputs of the analysis scripts used by the figures (parquet tables, per-unit AUCs, JSON statistics), plus the neuron selection (`selected_neurons.csv`, `swc_needed.txt`). |
| `data_tables/` | Table S1 (`Table_S1.xlsx`; sheets KC, ALPN, hemilineage, units, with CSV copies). |
| `figures/` | Figs 1-4 and Fig. S1 as submitted (PDF and PNG; 180 mm wide). |
| `qa/` | Automated figure audits (panel alignment, text collisions). |

## Data sources (public, not included)

MaleCNS v1.0 (Berg et al., 2026; CC-BY 4.0), https://male-cns.janelia.org. Files used, all in the public bucket `gs://flyem-male-cns/v1.0/`:

- `connectome-data/flat-connectome/body-annotations-male-cns-v1.0-minconf-0.5.feather` (cell class and type, hemilineage in the Ito-Lee and Truman nomenclatures, birth time, soma position, side, dimorphism);
- `connectome-data/flat-connectome/body-stats-male-cns-v1.0-minconf-0.5.feather` (presynapse, postsynapse and synaptic-weight counts; 778 MB; needed only for the connectivity analysis);
- `segmentation/skeletons-malecns/skeletons-swc/<bodyId>.swc` (skeletons; voxel units, 8 nm).

`code/download_data.py` fetches these files over HTTPS into `$FLYBRAIN_DATA` (default `data/`), downloading only the 17,225 skeletons listed in `results/swc_needed.txt` (Kenyon cells, antennal-lobe projection neurons and hemilineage neurons with soma positions; `results/selected_neurons.csv` holds their annotations).

## Setup

```bash
pip install -r requirements.txt
export FLYBRAIN_DATA=/path/to/data      # optional; default is ./data
python code/download_data.py            # add --no-stats to skip the 778 MB statistics table
```

## Pipeline (run from `code/`)

1. `kc_geometry.py` - Kenyon cells: peduncle segment per cell, per-hemisphere bundle axis, radial position and burial depth in 4 µm slabs, α/β-lobe radial position after k-means lobe assignment, soma distance from the clone axis. Writes `results/kc_geometry.parquet`, `kc_geometry_stats.json`, `kc_per_side.json`.
2. `pn_geometry.py` - anterodorsal and lateral projection neurons: literature birth ranks (Yu et al., 2010; Lin et al., 2012), soma measures within the lineage cluster and primary-neurite position within the bundle. Writes `results/pn_geometry.parquet`, `pn_geometry_stats.json`.
3. `hemilineage_geometry.py` - 155 hemilineage units (one hemilineage, neuromere and side; at least 5 primary and 10 secondary neurons): soma and neurite measures and per-unit AUCs. Writes `results/hemilineage_geometry.parquet`, `hemilineage_auc.csv`.
4. `birth_consequences.py` - synapse counts of primary versus secondary neurons per hemilineage, connectivity percentiles, dimorphism cross-tabulation. Writes `results/birth_consequences.json`, `dimorphism_birth.json`.
5. `dev_fig1.py` ... `dev_fig4.py`, `dev_figS1.py` - figures (`devfig.py` holds the shared styling). Figs 2, 3 and S1 need only `results/`; Fig. 1 re-reads the Kenyon cell skeletons for the cross-section and Fig. 4 needs the annotation and statistics tables. The optional render-time panel-alignment gate runs only if the `audit_panel_alignment` module is on the path; otherwise it is skipped with a message.

`strat.py` holds the shared geometry: skeleton parsing, geodesic paths, primary-neurite extraction, and `bundle_depth()` (bundle axis by PCA, slab cross-sections, convex-hull burial depth and radial distance to the geometric median, normalised by the 95th-percentile radius).

Steps 1-4 take about 20 minutes on a laptop once the skeletons are downloaded; the figure scripts take under two minutes.

## Citation

If you use these tables or code, please cite the paper (reference to be added on publication) and the Zenodo archive of this repository (https://doi.org/10.5281/zenodo.22760056), together with the MaleCNS connectome paper (Berg et al., 2026) and the lineage-tracing studies that provide the birth orders (Yu et al., 2010; Lin et al., 2012).

## License

Code: MIT License. Data tables and figures: CC BY 4.0.
