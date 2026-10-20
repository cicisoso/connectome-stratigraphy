"""Table 1 of the paper and the electronic supplementary material tables.
tables/Table1.csv, tables/Table1.md: summary of every test of the site-of-addition rule.
tables/ESM_tables.xlsx: all statistics, unit-level results, axes and cutting planes, and the birth-order ground truth.
tables/per_neuron/*.csv (zipped as tables/ESM_per_neuron_data.zip): the per-neuron measurements behind every figure."""
import os, json, zipfile, numpy as np, pandas as pd
import a2_lineages as A2

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
RES = f'{ROOT}/results'; TAB = f'{ROOT}/tables'; PN = f'{TAB}/per_neuron'
os.makedirs(TAB, exist_ok=True); os.makedirs(PN, exist_ok=True)
DS = ['MaleCNS', 'FlyWire', 'BANC']

KC = pd.read_csv(f'{RES}/summary_kc.csv'); L = pd.read_csv(f'{RES}/summary_lineages.csv'); P = pd.read_csv(f'{RES}/summary_lineages_pooled.csv')
E = pd.read_csv(f'{RES}/summary_adpn_epochs.csv'); H = pd.read_csv(f'{RES}/summary_hemilineages.csv'); C = pd.read_csv(f'{RES}/summary_connectivity.csv')
R = pd.read_csv(f'{RES}/summary_reproducibility.csv')
f2 = lambda x: '–' if x != x else f'{x:.2f}'.replace('-', '−')

# ------------------------------------------------------------------------------------------- Table 1
rows = []
def kc_row(measure, label):
    q = {ds: KC[(KC.ds == ds) & (KC.measure == measure)].iloc[0] for ds in DS}
    rows.append(['Kenyon cells: γ, α′β′, αβ', 'Lee et al. 1999', label] + [f'{f2(q[ds].rho_major)} ({f2(q[ds].C)})' for ds in DS] +
                [f'{sum(q[ds].rho_major < 0 for ds in DS)}/3 in the predicted direction'])   # both measures are predicted to fall with class
kc_row('ped_radial', 'peduncle fibre radial position, ρ (C)')
kc_row('soma_tang_clone', 'soma distance from clone axis, ρ (C)')
for lin, src, n in [('adPN', 'Yu et al. 2010', '33 types'), ('lPN', 'Lin et al. 2012', '12 types'), ('vPN', 'Lee et al. 2020', '8 types'), ('PAM', 'Lee et al. 2020', '15 types')]:
    q = {ds: L[(L.frame == 'hemilineage') & (L.lineage == lin) & (L.ds == ds) & (L.measure == 'geo_age')].iloc[0] for ds in DS}
    p = P[(P.frame == 'hemilineage') & (P.lineage == lin) & (P.measure == 'geo_age')].iloc[0]
    rows.append([f'{lin} lineage ({n})', src, 'geometric age index, ρ (C)'] + [f'{f2(q[ds].rho)} ({f2(q[ds].C)})' for ds in DS] +
                [f'pooled ρ = {f2(p.pooled_rho)} ({f2(p.ci_lo)} to {f2(p.ci_hi)})'])
rows.append(['adPN, embryonic- vs larval-born', 'Yu et al. 2010', 'geometric age index, AUC'] +
            [f2(E[E.ds == ds].auc_geo_age.iloc[0]) for ds in DS] + [f'{sum(E[E.ds == ds].auc_geo_age.iloc[0] > 0.5 for ds in DS)}/3 > 0.5'])
def hl_row(region, measure, label, lab2, ref):
    out = []
    for ds in DS:
        q = H[(H.ds == ds) & (H.region == region) & (H.measure == measure)]
        out.append('–' if not len(q) else f'{f2(q["median"].iloc[0])} ({int(q.n_consistent.iloc[0])}/{int(q.n_units.iloc[0])})')
    rows.append([label, ref, lab2] + out + [''])
hl_row('brain', 'auc_geo_age_soma', 'Brain hemilineages, primary vs secondary', 'soma index, median AUC (units > 0.5)', 'connectome annotations')
hl_row('vnc', 'auc_geo_age_soma', 'VNC hemilineages, primary vs secondary', 'soma index, median AUC (units > 0.5)', 'MANC birth classes')
hl_row('vnc', 'auc_early_vs_late_geo_age_soma', 'VNC hemilineages, early vs late secondary', 'soma index, median AUC (units > 0.5)', 'MANC birth classes')
T1 = pd.DataFrame(rows, columns=['Test', 'Birth-order evidence', 'Measure', 'MaleCNS', 'FlyWire', 'BANC', 'Across connectomes'])
for i in range(len(T1)):
    if T1.loc[i, 'Across connectomes'] == '' and T1.loc[i, 'Test'].startswith(('Brain', 'VNC')):
        vals = [float(v.split()[0]) for v in T1.loc[i, DS] if v != '–']; T1.loc[i, 'Across connectomes'] = f'{sum(v > 0.5 for v in vals)}/{len(vals)} medians > 0.5'
T1.to_csv(f'{TAB}/Table1.csv', index=False)
# the manuscript version cites the ground truths by citation key, numbered with the text when the manuscript is built
CITE = {'Lee et al. 1999': 'Lee et al. [@lee1999]', 'Yu et al. 2010': 'Yu et al. [@yu2010]', 'Lin et al. 2012': 'Lin et al. [@lin2012]',
        'Lee et al. 2020': 'Lee et al. [@lee2020]', 'MANC birth classes': 'MANC birth classes [@marin2024]'}
T1m = T1.assign(**{'Birth-order evidence': T1['Birth-order evidence'].map(lambda s: CITE.get(s, s))})
md = '| ' + ' | '.join(T1m.columns) + ' |\n|' + '---|' * len(T1m.columns) + '\n' + '\n'.join('| ' + ' | '.join(map(str, r)) + ' |' for r in T1m.values)
open(f'{TAB}/Table1.md', 'w').write(md + '\n'); print(md)

# ------------------------------------------------------------------------------------------- ESM workbook
def flat(d, prefix=''):
    out = {}
    for k, v in d.items():
        if isinstance(v, dict): out.update(flat(v, f'{prefix}{k}.'))
        elif isinstance(v, list): out[f'{prefix}{k}'] = ', '.join(f'{x:.4g}' if isinstance(x, float) else str(x) for x in v)
        else: out[f'{prefix}{k}'] = v
    return out
axes = []
for ds in DS:
    for kind, f in [('Kenyon cell bundle', 'axes_kc'), ('lineage frame', 'axes_lineages'), ('hemilineage unit', 'axes_hemilineages')]:
        for r in json.load(open(f'{RES}/{f}_{ds}.json')): axes.append(dict(kind=kind, **flat(r)))
AX = pd.DataFrame(axes)
lead = ['kind', 'ds', 'bundle', 'lineage', 'hl', 'unit', 'side', 'frame_type', 'frame']; AX = AX[[c for c in lead if c in AX] + [c for c in AX if c not in lead]]

pb, mixed = A2.pam_birth(); vb = A2.vpn_birth()
gt = [dict(lineage='adPN', source='Yu et al. 2010 (birth sequence of the anterodorsal PN lineage)', type=k, birth_rank=v,
           epoch='embryonic' if v <= A2.ADPN_EMBRYONIC_LAST else 'larval') for k, v in A2.ADPN.items()]
gt += [dict(lineage='lPN', source='Lin et al. 2012 (lateral PN lineage)', type=k, birth_rank=v) for k, v in A2.LPN.items()]
gt += [dict(lineage='vPN', source='Lee et al. 2020 (ALv1, mean clone-induction window)', type=k, birth_rank=v) for k, v in vb.items()]
for hl, g in pb.items():
    gt += [dict(lineage='PAM', hemilineage=hl, source='Lee et al. 2020 (clone-weighted mean induction window)', type=k, birth_rank=r.mean_window,
                n_clones=r.n_clones, n_lee_types=r.n_lee_types, mixed_windows=k in mixed) for k, r in g.iterrows()]
GT = pd.DataFrame(gt)

UH = pd.concat([pd.read_csv(f'{RES}/hemilineage_units_{ds}.csv') for ds in DS])
UC = pd.concat([pd.read_csv(f'{RES}/connectivity_units_{ds}.csv') for ds in DS + ['MANC']])
sheets = [('README', None), ('S1_Kenyon_cells', KC), ('S2_lineages', L), ('S3_lineages_pooled', P), ('S4_adPN_epochs', E),
          ('S5_type_reproducibility', R), ('S6_hemilineage_summary', H), ('S7_connectivity_summary', C), ('S8_hemilineage_units', UH),
          ('S9_connectivity_units', UC), ('S10_axes_and_sections', AX), ('S11_birth_order_ground_truth', GT)]
DESC = {
 'S1_Kenyon_cells': 'Kenyon cells: Spearman rho of each measure with major class (gamma 0, alpha-prime/beta-prime 1, alpha/beta 2), Fisher 95% CI, concordance C, adjacent-class AUCs, 7-class rho (MaleCNS, FlyWire), per-side rho and class medians.',
 'S2_lineages': 'Clonally dated lineages: per connectome, frame (whole hemilineage or projection neurons only), lineage and measure: n neurons, n types, rho with birth rank and Fisher CI, type-median rho, concordance C (predicted sign applied), and whether the sign was as predicted.',
 'S3_lineages_pooled': 'Fixed-effect pooled rho across the three connectomes (Fisher z, weights (n-3)/1.06^2), with 95% CI, number of connectomes in the predicted direction and mean C.',
 'S4_adPN_epochs': 'adPN neurons of embryonic-born (birth rank <= 17) versus larval-born types: AUC of the geometric age index and the primary/secondary labels of the connectome annotation in each group.',
 'S5_type_reproducibility': 'Spearman rho between type medians of the geometric age index in pairs of connectomes.',
 'S6_hemilineage_summary': 'Hemilineage units (hemilineage x side x neuromere or brain): median and quartiles of per-unit AUCs (primary neurons look older) and of the three-level rho, with the number of units in the predicted direction and the two-sided sign-test p.',
 'S7_connectivity_summary': 'Primary versus secondary neurons: median per-unit log2 ratio and AUC of input, output and total synapses, size (cable length; MANC: volume) and synapses per size; median difference in partner-class fractions; the fraction of primary neurons among the top 5% by synapse count. _k = units in the predicted direction, _p = sign-test p.',
 'S8_hemilineage_units': 'Per-unit results for the hemilineage analysis.',
 'S9_connectivity_units': 'Per-unit results for the connectivity analysis.',
 'S10_axes_and_sections': 'Every axis, centroid and cutting plane used: Kenyon-cell peduncle cylinders and slab ranges; for each lineage frame and hemilineage unit the soma centroid, fibre outflow centroid, soma depth axis, tract convergence point, tract axis and radius, and the cross-section basis. Coordinates in micrometres in each connectome\'s own space.',
 'S11_birth_order_ground_truth': 'Birth order used for each connectome cell type: adPN birth rank (Yu et al. 2010), lPN rank (Lin et al. 2012), vPN and PAM mean clone-induction windows computed from Lee et al. 2020 supplementary file 1.'}
with pd.ExcelWriter(f'{TAB}/ESM_tables.xlsx') as xw:
    pd.DataFrame({'sheet': [s for s, _ in sheets[1:]], 'contents': [DESC[s] for s, _ in sheets[1:]]}).to_excel(xw, 'README', index=False)
    for s, df in sheets[1:]: df.to_excel(xw, s, index=False)
print('ESM tables:', {s: (None if df is None else df.shape) for s, df in sheets})

# ------------------------------------------------------------------------------------------- per-neuron data
files = []
for ds in DS:
    for name, f in [('kenyon_cells', f'kc_{ds}'), ('lineages_hemilineage_frame', f'lineages_hemilineage_{ds}'),
                    ('lineages_pn_frame', f'lineages_pn_only_{ds}'), ('hemilineage_units', f'hemilineages_{ds}')]:
        out = f'{PN}/{ds}_{name}.csv'; pd.read_parquet(f'{RES}/{f}.parquet').to_csv(out, index=False, float_format='%.5g'); files.append(out)
for ds in DS + ['MANC']:
    out = f'{PN}/{ds}_connectivity.csv'; pd.read_parquet(f'{RES}/connectivity_{ds}.parquet').to_csv(out, index=False, float_format='%.5g'); files.append(out)
readme = f'{PN}/README.txt'
open(readme, 'w').write("""Per-neuron measurements behind every figure (one CSV per connectome and analysis).
Coordinates in micrometres in each connectome's own space; nid = the connectome's neuron identifier.
kenyon_cells: kctype, rank (0 gamma-d, 1 gamma-m, 2 alpha'/beta', 3-6 alpha/beta p, s, m, c), clone, ped_radial and ped_burial (fraction of
  peduncle radius), lobe_radial (alpha/beta lobe), soma_tang_clone (um from the clone axis).
lineages_*: lineage, birth rank (see ESM table S11), epoch (adPN), frame, in_section, soma_dist, soma_depth, soma_tang (um), radial and burial
  (fraction of tract radius), geo_age (mean of the z-scored, sign-aligned soma depth, soma distance and fibre radial position).
hemilineage_units: unit (hemilineage|side|neuromere), birth class (primary, early_secondary, secondary), the same measures, geo_age_soma
  (soma measures only).
connectivity: synapse counts, size (MaleCNS and BANC: cable length in um; FlyWire: skeleton length; MANC: volume), synapses per size and
  the fraction of input and output synapses from each partner class.
""")
with zipfile.ZipFile(f'{TAB}/ESM_per_neuron_data.zip', 'w', zipfile.ZIP_DEFLATED) as z:
    for f in files + [readme]: z.write(f, os.path.basename(f))
print('zip MB', round(os.path.getsize(f'{TAB}/ESM_per_neuron_data.zip') / 1e6, 1))
