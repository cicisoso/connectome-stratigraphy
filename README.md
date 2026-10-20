# Connectome stratigraphy: neuronal birth order in the geometry of *Drosophila* connectomes

Code, per-neuron tables and figures for the manuscript *"Neuronal birth order leaves a reproducible trace in the geometry of three Drosophila connectomes"* (Liu, Zong, Chen and Xiong; School of Information Technology, Zhejiang Financial College; submitted to *Biology Letters*).

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22760055.svg)](https://doi.org/10.5281/zenodo.22760055)

We test a site-of-addition rule: a neuron's age grows with its distance from the site where its lineage adds new cells, which is the neuroblast for somata and the core of the lineage tract for fibres. In three electron-microscopy connectomes of adult *Drosophila* (MaleCNS v1.0, FlyWire v783 and BANC v888), Kenyon cells show the expected concentric order (Spearman ρ between peduncle radial position and class −0.28 to −0.45). In four lineages whose birth order was mapped by clonal analysis, an index combining soma depth, soma distance and fibre radial position correlates with birth rank (pooled ρ 0.20–0.47) and separates embryonic- from larval-born projection neurons (area under the curve 0.70–0.84). Primary (embryonic-born) somata look older than secondary somata in 230 of 261 hemilineage units. The rule is a tendency rather than a clock: of two neurons born at different times, the later-born has the younger-looking geometry in 55–73% of pairs. With the MANC v1.0 nerve cord added for connectivity, primary neurons are larger and carry 2.9–6.6 times more synapses than secondary neurons of the same unit, with little difference in partner classes.

## Versions

- **v2.0.0** (this version; *Biology Letters* submission): adds `biolett/`, which repeats every analysis in MaleCNS v1.0, FlyWire v783 and BANC v888 (MANC v1.0 for connectivity), uses the clone-induction windows of Lee et al. (2020) as birth order for the vPN and PAM lineages, extends the hemilineage test to the brain, reports concordance and areas under the curve, and compares synapse counts, size and partner classes of primary and secondary neurons. It also holds the figures, Table 1, tables S1–S11 and per-neuron data of the current manuscript.
- **v1.0.0** (https://doi.org/10.5281/zenodo.22760056): the MaleCNS-only analyses of an earlier version of the manuscript (`code/`, `results/`, `data_tables/`, `figures/`, `qa/`), kept unchanged. `biolett/` does not use them.

The connectivity numbers of v1.0.0 are superseded. v1.0.0 pooled each of 45 nerve-cord hemilineages and ranked neurons by synaptic weight; `biolett/` compares primary and secondary neurons within units of one hemilineage, side and neuromere and ranks neurons by synapse count. For the MaleCNS nerve cord the median ratio of synapse counts is then 6.6 rather than about 4, and primary neurons make up 71% of the 5% most connected neurons rather than 95%.

## Contents

| Path | What it holds |
|---|---|
| `biolett/` | v2.0.0: analysis and figure scripts (`code/`, run in order by `run_all.sh`), download helpers (`tools/`), supplementary files 1 and 2 of Lee et al. (2020) (`external/`), intermediate results (`results/`), Figs 1–3 and S1–S6 (`figures/`, PDF and PNG), and Table 1, tables S1–S11 (`ESM_tables.xlsx`) and the per-neuron CSV files (`tables/`). `biolett/README.md` lists the data sources, requirements and scripts. |
| `code/`, `results/`, `data_tables/`, `figures/`, `qa/`, `requirements.txt` | v1.0.0 (MaleCNS only; see below). |

## v2.0.0 quick start (`biolett/`)

```bash
pip install -r biolett/requirements.txt
export FLYBRAIN_DATA=/path/to/data      # file layout in biolett/README.md; default biolett/data
cd biolett/code
zsh ../run_all.sh                       # every result, figure and table
```

The figure scripts and `make_tables.py` read only `biolett/results/` and `biolett/external/`, so they run without the raw data, for example `python fig2.py`. The other steps need the raw connectome files listed in `biolett/README.md`.

## v1.0.0 (MaleCNS only)

### Data sources (public, not included)

MaleCNS v1.0 (Berg et al., 2026; CC-BY 4.0), https://male-cns.janelia.org. Files used, all in the public bucket `gs://flyem-male-cns/v1.0/`:

- `connectome-data/flat-connectome/body-annotations-male-cns-v1.0-minconf-0.5.feather` (cell class and type, hemilineage in the Ito-Lee and Truman nomenclatures, birth time, soma position, side, dimorphism);
- `connectome-data/flat-connectome/body-stats-male-cns-v1.0-minconf-0.5.feather` (presynapse, postsynapse and synaptic-weight counts; 778 MB; needed only for the connectivity analysis);
- `segmentation/skeletons-malecns/skeletons-swc/<bodyId>.swc` (skeletons; voxel units, 8 nm).

`code/download_data.py` fetches these files over HTTPS into `$FLYBRAIN_DATA` (default `data/`), downloading only the 17,225 skeletons listed in `results/swc_needed.txt` (Kenyon cells, antennal-lobe projection neurons and hemilineage neurons with soma positions; `results/selected_neurons.csv` holds their annotations).

### Setup

```bash
pip install -r requirements.txt
export FLYBRAIN_DATA=/path/to/data      # optional; default is ./data
python code/download_data.py            # add --no-stats to skip the 778 MB statistics table
```

### Pipeline (run from `code/`)

1. `kc_geometry.py` - Kenyon cells: peduncle segment per cell, per-hemisphere bundle axis, radial position and burial depth in 4 µm slabs, α/β-lobe radial position after k-means lobe assignment, soma distance from the clone axis. Writes `results/kc_geometry.parquet`, `kc_geometry_stats.json`, `kc_per_side.json`.
2. `pn_geometry.py` - anterodorsal and lateral projection neurons: literature birth ranks (Yu et al., 2010; Lin et al., 2012), soma measures within the lineage cluster and primary-neurite position within the bundle. Writes `results/pn_geometry.parquet`, `pn_geometry_stats.json`.
3. `hemilineage_geometry.py` - 155 hemilineage units (one hemilineage, neuromere and side; at least 5 primary and 10 secondary neurons): soma and neurite measures and per-unit AUCs. Writes `results/hemilineage_geometry.parquet`, `hemilineage_auc.csv`.
4. `birth_consequences.py` - synapse counts of primary versus secondary neurons per hemilineage, connectivity percentiles, dimorphism cross-tabulation. Writes `results/birth_consequences.json`, `dimorphism_birth.json`.
5. `dev_fig1.py` ... `dev_fig4.py`, `dev_figS1.py` - figures (`devfig.py` holds the shared styling). Figs 2, 3 and S1 need only `results/`; Fig. 1 re-reads the Kenyon cell skeletons for the cross-section and Fig. 4 needs the annotation and statistics tables. The optional render-time panel-alignment gate runs only if the `audit_panel_alignment` module is on the path; otherwise it is skipped with a message.

`strat.py` holds the shared geometry: skeleton parsing, geodesic paths, primary-neurite extraction, and `bundle_depth()` (bundle axis by PCA, slab cross-sections, convex-hull burial depth and radial distance to the geometric median, normalised by the 95th-percentile radius).

Steps 1-4 take about 20 minutes on a laptop once the skeletons are downloaded; the figure scripts take under two minutes.

## Citation

Please cite the paper (reference to be added on publication) and the Zenodo archive of the release you used (v1.0.0: https://doi.org/10.5281/zenodo.22760056; all versions: https://doi.org/10.5281/zenodo.22760055), together with the connectome papers and data deposits: MaleCNS (Berg et al., 2026), FlyWire (Dorkenwald et al., 2024; Schlegel et al., 2024), BANC (Bates et al., 2026, and the Dataverse deposit https://doi.org/10.7910/DVN/7WTH1N) and MANC (Takemura et al., 2024; Marin et al., 2024), and the lineage studies that provide the birth orders (Lee et al., 1999; Yu et al., 2010; Lin et al., 2012; Lee et al., 2020).

## License

Code: MIT License. Tables, results and figures: CC BY 4.0. `biolett/external/` holds supplementary files 1 and 2 of Lee et al. (2020, *eLife* 9, e53518), redistributed under their CC BY 4.0 licence.
