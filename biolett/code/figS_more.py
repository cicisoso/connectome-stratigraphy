"""Electronic supplementary figures S3-S6.
S3: Kenyon cells by subclass, peduncle burial and alpha/beta lobe position. S4: lineage results by frame and measure, and
reproducibility of type medians across connectomes. S5: every measure in the hemilineage units, and example units.
S6: connectivity of primary versus secondary neurons, further measures."""
from figstyle import *
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
from scipy import stats

MAJOR = {0: 0, 1: 0, 2: 1, 3: 2, 4: 2, 5: 2, 6: 2}
SUBLAB = {0: 'γd', 1: 'γm', 2: "α′β′", 3: 'αβp', 4: 'αβs', 5: 'αβm', 6: 'αβc'}
SUBCOL = {0: '#D4D4D4', 1: '#9C9C9C', 2: '#56B4E9', 3: '#F3CDE2', 4: '#E3A3C8', 5: '#CC79A7', 6: '#8C4770'}
REG = {'brain': 'brain', 'vnc': 'VNC'}; MK = {'brain': 'o', 'vnc': 's'}
BCOL = {'primary': PRIM, 'early_secondary': '#CC79A7', 'secondary': '#BDBDBD'}
BLAB = {'primary': 'primary', 'early_secondary': 'early secondary', 'secondary': 'late secondary'}

def med_iqr(ax, x, v, w=0.32):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if not len(v): return
    q1, m, q3 = np.percentile(v, [25, 50, 75])
    ax.plot([x, x], [q1, q3], color='k', lw=1.4, solid_capstyle='butt', zorder=4); ax.plot([x - w, x + w], [m, m], color='k', lw=1.1, zorder=4)

def ds_labels(ax, centres, names, y=-0.07, second=None):
    ax.set_xticks(centres); ax.set_xticklabels([''] * len(centres)); ax.tick_params(axis='x', length=0)
    for i, (c, ds) in enumerate(zip(centres, names)):
        ax.text(c, y, ds, transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=6.5, color=DSC[ds], fontweight='bold')
        if second: ax.text(c, y - 0.09, second[i], transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=6.5)

# ======================================================================================== S3: Kenyon cells
KC = {ds: pd.read_parquet(f'{RES}/kc_{ds}.parquet').assign(major=lambda t: t['rank'].map(MAJOR)) for ds in DS}
fig = plt.figure(figsize=(174 * MM, 70 * MM))
gs = gridspec.GridSpec(1, 3, width_ratios=[1.5, 1, 1], left=0.07, right=0.99, top=0.8, bottom=0.17, wspace=0.38)
axA, axB, axC = [fig.add_subplot(gs[0, i]) for i in range(3)]
cent = []; pos = 0
for ds in DS:
    T = KC[ds].dropna(subset=['ped_radial']); ks = [0, 1, 2, 3, 4, 5, 6] if ds != 'BANC' else [0, 1, 2]
    key = 'rank' if ds != 'BANC' else 'major'; cols = [SUBCOL[k] for k in ks] if ds != 'BANC' else [KC_MAJOR_COL[k] for k in ks]
    xs = [pos + i for i in range(len(ks))]; boxes(axA, [T[T[key] == k].ped_radial for k in ks], xs, cols, width=0.7)
    r7 = stats.spearmanr(T[key], T.ped_radial)[0]
    axA.text(np.mean(xs), 1.0, f'ρ = {r7:.2f}'.replace('-', '−'), transform=axA.get_xaxis_transform(), ha='center', va='bottom', fontsize=6.5)
    cent.append(np.mean(xs)); pos += len(ks) + 1.2
ds_labels(axA, cent, DS, y=-0.04)
axA.set_ylabel('Peduncle radial position\n(fraction of radius)'); axA.set_ylim(0, 1.1)
axA.legend(handles=[Line2D([], [], marker='s', ls='', color=SUBCOL[k], ms=4, label=SUBLAB[k]) for k in range(7)], ncol=4, loc='lower center',
           bbox_to_anchor=(0.5, 1.07), fontsize=6.5, handletextpad=0.05, columnspacing=0.5, borderaxespad=0)
for g, ds in enumerate(DS):
    T = KC[ds].dropna(subset=['ped_burial']); xs = [g * 4 + k for k in [0, 1, 2]]
    boxes(axB, [T[T.major == k].ped_burial for k in [0, 1, 2]], xs, [KC_MAJOR_COL[k] for k in [0, 1, 2]], width=0.7)
    axB.text(g * 4 + 1, 1.0, f'ρ = {stats.spearmanr(T.major, T.ped_burial)[0]:.2f}'.replace('-', '−'), transform=axB.get_xaxis_transform(), ha='center', va='bottom', fontsize=6.5)
ds_labels(axB, [1, 5, 9], DS, y=-0.04); axB.set_ylabel('Peduncle burial depth\n(fraction of radius)')
axB.legend(handles=[Line2D([], [], marker='s', ls='', color=KC_MAJOR_COL[k], ms=4, label=KC_MAJOR_LAB[k]) for k in [0, 1, 2]], ncol=3, loc='lower center',
           bbox_to_anchor=(0.5, 1.07), fontsize=6.5, handletextpad=0.05, columnspacing=0.6, borderaxespad=0)
for g, ds in enumerate(['MaleCNS', 'FlyWire']):
    T = KC[ds].dropna(subset=['lobe_radial']); T = T[T['rank'] >= 3]; xs = [g * 5 + k for k in range(4)]
    boxes(axC, [T[T['rank'] == k].lobe_radial for k in [3, 4, 5, 6]], xs, [SUBCOL[k] for k in [3, 4, 5, 6]], width=0.7)
    axC.text(g * 5 + 1.5, 1.0, f'ρ = {stats.spearmanr(T["rank"], T.lobe_radial)[0]:.2f}'.replace('-', '−'), transform=axC.get_xaxis_transform(), ha='center', va='bottom', fontsize=6.5)
ds_labels(axC, [1.5, 6.5], ['MaleCNS', 'FlyWire'], y=-0.04); axC.set_ylabel('α/β lobe radial position\n(fraction of radius)')
axC.text(0.5, 1.08, 'definitional: αβ subclasses\nare named by lobe position', transform=axC.transAxes, ha='center', va='bottom', fontsize=6.5, color=GREY)
fig.canvas.draw(); panel_labels(fig, [axA, axB, axC], 'abc', dy_pt=3)
gate(fig, 'FigS3', axes=[axB, axC], panel_ids=['b', 'c']); export(fig, 'FigS3'); plt.close(fig)

# ======================================================================================== S4: lineages by frame and measure; reproducibility
L = pd.read_csv(f'{RES}/summary_lineages.csv')
SETS = ['adPN', 'lPN', 'vPN', 'PAM', 'PAM_single_window']; SLAB = {'PAM_single_window': 'PAM, single\nwindow'}
MS = [('geo_age', 1, 'geometric\nage'), ('soma_depth', 1, 'soma\ndepth'), ('soma_dist', -1, '−soma\ndistance'), ('radial', -1, '−fibre\nradial')]
rows = [(s, ds) for s in SETS for ds in DS]; colsM = [(fr, m) for fr in ['hemilineage', 'pn_only'] for m in MS]
H = np.full((len(rows), len(colsM)), np.nan)
for i, (s, ds) in enumerate(rows):
    for j, (fr, (m, sg, _)) in enumerate(colsM):
        q = L[(L.frame == fr) & (L.lineage == s) & (L.ds == ds) & (L.measure == m)]
        if len(q): H[i, j] = sg * q.rho.iloc[0]
fig = plt.figure(figsize=(174 * MM, 112 * MM))
gs = gridspec.GridSpec(2, 2, width_ratios=[1.35, 1], height_ratios=[1, 1], left=0.13, right=0.985, top=0.9, bottom=0.13, wspace=0.42, hspace=0.55)
axA = fig.add_subplot(gs[:, 0]); axB = fig.add_subplot(gs[0, 1]); axC = fig.add_subplot(gs[1, 1])
im = axA.imshow(H, cmap='PuOr_r', norm=TwoSlopeNorm(0, -0.8, 0.8), aspect='auto')
for i in range(H.shape[0]):
    for j in range(H.shape[1]):
        axA.text(j, i, '–' if np.isnan(H[i, j]) else f'{H[i, j]:.2f}'.replace('-', '−'), ha='center', va='center', fontsize=5.5,
                 color='white' if abs(np.nan_to_num(H[i, j])) > 0.55 else 'k')
axA.set_yticks(range(len(rows))); axA.set_yticklabels([ds for _, ds in rows], fontsize=6.5)
for t, (_, ds) in zip(axA.get_yticklabels(), rows): t.set_color(DSC[ds])
for k, s in enumerate(SETS):
    axA.text(-0.19, 1 - (k * 3 + 1.5) / len(rows), SLAB.get(s, s), transform=axA.transAxes, ha='center', va='center', fontsize=6.5, fontweight='bold', rotation=90)
    if k: axA.axhline(k * 3 - 0.5, color='white', lw=1.5)
axA.axvline(3.5, color='k', lw=0.8)
axA.set_xticks(range(len(colsM))); axA.set_xticklabels([lab for _, (_, _, lab) in colsM], fontsize=6)
axA.xaxis.tick_top(); axA.tick_params(length=0)
axA.text(1.5, -1.9, 'whole-hemilineage frame', ha='center', va='bottom', fontsize=6.5, fontweight='bold')
axA.text(5.5, -1.9, 'projection-neuron frame', ha='center', va='bottom', fontsize=6.5, fontweight='bold')
for s_ in axA.spines.values(): s_.set_visible(False)
cax = axA.inset_axes([0.15, -0.035, 0.7, 0.016]); cb = fig.colorbar(im, cax=cax, orientation='horizontal'); cb.outline.set_linewidth(0.4)
cb.ax.tick_params(labelsize=6.5, length=1.5); cb.set_label('Spearman ρ, signed so that the predicted direction is positive', fontsize=6.5)
TM = {}
for ds in DS:
    M = pd.read_parquet(f'{RES}/lineages_hemilineage_{ds}.parquet'); M = M[M['rank'].notna()]
    TM[ds] = M.groupby(['lineage', 'type', 'side']).geo_age.median().groupby(['lineage', 'type']).mean()
LC = {'adPN': '#0072B2', 'lPN': '#E69F00', 'vPN': '#009E73', 'PAM': '#CC79A7'}
REP = pd.read_csv(f'{RES}/summary_reproducibility.csv')
for ax, b in [(axB, 'FlyWire'), (axC, 'BANC')]:
    j = pd.concat([TM['MaleCNS'].rename('a'), TM[b].rename('b')], axis=1).dropna()
    for lin, g in j.groupby(level=0): ax.scatter(g.a, g.b, s=9, color=LC[lin], lw=0, alpha=0.9, label=lin)
    lim = [-2.6, 1.6]; ax.plot(lim, lim, color=GREY, lw=0.6, ls=(0, (2, 2))); ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect('equal')
    r = REP[(REP.pair == f'MaleCNS-{b}') & (REP.lineage == 'all')].iloc[0]
    ax.text(0.04, 0.96, f'ρ = {r.rho:.2f}, {int(r.n_types)} types'.replace('-', '−'), transform=ax.transAxes, ha='left', va='top', fontsize=6.5)
    ax.set_xlabel('MaleCNS type median'); ax.set_ylabel(f'{b} type median')
axB.legend(loc='lower right', fontsize=6.5, handletextpad=0.1, borderaxespad=0.1, markerscale=1.2)
fig.canvas.draw(); panel_labels(fig, [axA, axB, axC], 'abc', dy_pt=3)
export(fig, 'FigS4'); plt.close(fig)

# ======================================================================================== S5: hemilineage units, all measures, and example units
U = {ds: pd.read_csv(f'{RES}/hemilineage_units_{ds}.csv') for ds in DS}
HG = [('MaleCNS', 'brain'), ('FlyWire', 'brain'), ('BANC', 'brain'), ('MaleCNS', 'vnc'), ('BANC', 'vnc')]
HM = [('auc_soma_depth', 'soma depth'), ('auc_soma_dist', 'soma distance'), ('auc_soma_tang', 'soma tangential'), ('auc_geo_age_soma', 'soma index'),
      ('auc_radial', 'fibre radial*'), ('auc_burial', 'fibre burial*'), ('auc_geo_age', 'geometric age*')]
fig = plt.figure(figsize=(174 * MM, 132 * MM))
gs = gridspec.GridSpec(2, 4, height_ratios=[1, 1.05], left=0.075, right=0.99, top=0.92, bottom=0.08, wspace=0.12, hspace=0.42)
axA = fig.add_subplot(gs[0, :]); axs = [fig.add_subplot(gs[1, i]) for i in range(4)]
xt = []
for m, (col, lab) in enumerate(HM):
    for g, (ds, rg) in enumerate(HG):
        T = U[ds][U[ds].region == rg]; x = m * (len(HG) + 1.3) + g; v = T[col].dropna()
        strip(axA, v, x, DSC[ds], width=0.6, s=3 if rg == 'vnc' else 7, alpha=0.7, marker=MK[rg], seed=5 * m + g); med_iqr(axA, x, v, w=0.36)
    xt.append(m * (len(HG) + 1.3) + (len(HG) - 1) / 2)
axA.axhline(0.5, color=GREY, lw=0.6, ls=(0, (2, 2)), zorder=0); axA.set_ylim(0, 1); axA.set_xlim(-1, len(HM) * (len(HG) + 1.3) - 1.3)
axA.set_xticks(xt); axA.set_xticklabels([l for _, l in HM]); axA.tick_params(axis='x', length=0); axA.set_ylabel('AUC, primary neurons look older')
axA.legend(handles=[Line2D([], [], marker='s', ls='', color=DSC[ds], ms=3.5, label=ds) for ds in DS] +
                   [Line2D([], [], marker=MK[r], ls='', mfc='white', mec='k', ms=3.5, label=REG[r]) for r in ['brain', 'vnc']],
           loc='lower left', bbox_to_anchor=(0.0, 1.0), ncol=5, fontsize=6.5, handletextpad=0.1, columnspacing=0.8, borderaxespad=0.2)
axA.text(1.0, 1.02, '* uses fibre position, a criterion of the primary/secondary annotation', transform=axA.transAxes, ha='right', va='bottom', fontsize=6.5, color=GREY)

def unit_projection(ax, ds, hl, side):
    S = pd.read_parquet(f'{RES}/render_somata_{ds}.parquet'); S = S[(S.frame_hl == hl) & (S.frame_side == side)].copy()
    B = pd.read_parquet(f'{RES}/hemilineages_{ds}.parquet', columns=['nid', 'birth']).drop_duplicates('nid').set_index('nid').birth
    S['birth'] = S.nid.map(B)
    PA = pd.read_parquet(f'{RES}/render_paths_{ds}.parquet'); PA = PA[(PA.frame_hl == hl) & (PA.frame_side == side)]
    fr = [r for r in json.load(open(f'{RES}/render_frames_{ds}.json')) if r['hl'] == hl and r['side'] == side][0]
    cs = np.array(fr['soma_centroid']); axv = np.array(fr['soma_axis']); oc = np.array(fr['outflow_centroid'])
    X = S[['x', 'y', 'z']].values - cs; perp = X - np.outer(X @ axv, axv); e = np.linalg.svd(perp, full_matrices=False)[2][0]
    pr = lambda Q: np.c_[(Q - cs) @ e, (Q - cs) @ axv]
    for _, g in PA.groupby('nid'):
        q = pr(g[['x', 'y', 'z']].values); ax.plot(q[:, 0], q[:, 1], color='#D9D9D9', lw=0.3, zorder=1)
    for b in ['secondary', 'early_secondary', 'primary']:
        q = pr(S[S.birth == b][['x', 'y', 'z']].values); ax.scatter(q[:, 0], q[:, 1], s=7, color=BCOL[b], edgecolors='k', lw=0.15, zorder=3)
    o = pr(oc[None])[0]; ax.plot(*o, marker='s', color='k', ms=3, zorder=5); ax.plot(0, 0, marker='+', color='k', ms=5, mew=0.9, zorder=5)
    ax.annotate('', xy=(0, 15), xytext=(o[0], o[1]), arrowprops=dict(arrowstyle='->', lw=0.7), zorder=4)
    ax.set_aspect('equal', adjustable='datalim'); clean_3d_axes(ax)
    xl = ax.get_xlim(); yl = ax.get_ylim(); scalebar(ax, xl[0] + 1, yl[0] + 1, 20, '20 µm', dy=-1.2, fontsize=6.5)
    return S

def unit_section(ax, S):
    s = S[S.u.notna()]
    for b in ['secondary', 'early_secondary', 'primary']:
        q = s[s.birth == b]; ax.scatter(q.u, q.v, s=9, color=BCOL[b], edgecolors='k', lw=0.15, zorder=3)
    ax.plot(0, 0, marker='+', color='k', ms=5, mew=0.9, zorder=4)
    ax.set_aspect('equal', adjustable='datalim'); clean_3d_axes(ax)
    xl = ax.get_xlim(); yl = ax.get_ylim(); scalebar(ax, xl[0] + 0.2, yl[0] + 0.2, 2, '2 µm', dy=-0.15, fontsize=6.5)

for k, (ds, hl, side, title) in enumerate([('MaleCNS', '13A', 'L|T1', 'MaleCNS 13A, left T1'), ('BANC', '13A', 'L|T1', 'BANC 13A, left T1')]):
    S = unit_projection(axs[2 * k], ds, hl, side); unit_section(axs[2 * k + 1], S)
    axs[2 * k].set_title(f'{title}: somata', fontsize=7); axs[2 * k + 1].set_title(f'{title}: tract section', fontsize=7)
fig.legend(handles=[Line2D([], [], marker='o', ls='', color=BCOL[b], mec='k', mew=0.2, ms=3.5, label=BLAB[b]) for b in ['primary', 'early_secondary', 'secondary']],
           loc='lower center', bbox_to_anchor=(0.5, 0.005), ncol=3, fontsize=6.5, handletextpad=0.1, columnspacing=1.2)
fig.canvas.draw(); panel_labels(fig, [axA] + axs, 'abcde', dy_pt=3)
export(fig, 'FigS5'); plt.close(fig)

# ======================================================================================== S6: connectivity, further measures
CU = {ds: pd.read_csv(f'{RES}/connectivity_units_{ds}.csv') for ds in DS + ['MANC']}
SC = pd.read_csv(f'{RES}/summary_connectivity.csv')
CG = HG + [('MANC', 'vnc')]
CM = [('auc_n_in', 'input\nsynapses'), ('auc_n_out', 'output\nsynapses'), ('auc_n_syn', 'all\nsynapses'), ('auc_size', 'size'), ('auc_syn_per_size', 'synapses\nper size')]
fig = plt.figure(figsize=(174 * MM, 66 * MM))
gs = gridspec.GridSpec(1, 2, width_ratios=[2.2, 1], left=0.075, right=0.99, top=0.84, bottom=0.25, wspace=0.28)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
xt = []
for m, (col, lab) in enumerate(CM):
    for g, (ds, rg) in enumerate(CG):
        T = CU[ds][CU[ds].region == rg]; x = m * (len(CG) + 1.3) + g; v = T[col].dropna()
        strip(axA, v, x, DSC[ds], width=0.6, s=3 if rg == 'vnc' else 7, alpha=0.7, marker=MK[rg], seed=11 * m + g); med_iqr(axA, x, v, w=0.36)
    xt.append(m * (len(CG) + 1.3) + (len(CG) - 1) / 2)
axA.axhline(0.5, color=GREY, lw=0.6, ls=(0, (2, 2)), zorder=0); axA.set_ylim(0, 1.02)
axA.set_xticks(xt); axA.set_xticklabels([l for _, l in CM]); axA.tick_params(axis='x', length=0)
axA.set_ylabel('AUC, primary > secondary\n(per hemilineage unit)')
axA.legend(handles=[Line2D([], [], marker='s', ls='', color=DSC[ds], ms=3.5, label=ds) for ds in DS + ['MANC']] +
                   [Line2D([], [], marker=MK[r], ls='', mfc='white', mec='k', ms=3.5, label=REG[r]) for r in ['brain', 'vnc']],
           loc='lower center', bbox_to_anchor=(0.5, 1.01), ncol=6, fontsize=6.5, handletextpad=0.1, columnspacing=0.8, borderaxespad=0)
w = 0.38
for g, (ds, rg) in enumerate(CG):
    r = SC[(SC.ds == ds) & (SC.region == rg)].iloc[0]
    axB.bar(g - w / 2, r.top5_frac_primary_all, w, color='white', edgecolor=DSC[ds], lw=0.8)
    axB.bar(g + w / 2, r.top5_frac_primary_top, w, color=DSC[ds], edgecolor=DSC[ds], lw=0.8)
axB.set_xticks(range(len(CG))); axB.set_xticklabels([f'{ds} {REG[rg]}' for ds, rg in CG], rotation=40, ha='right', fontsize=6.5)
for t, (ds, rg) in zip(axB.get_xticklabels(), CG): t.set_color(DSC[ds])
axB.set_ylabel('Fraction primary'); axB.set_ylim(0, 1)
axB.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc='white', ec='k', lw=0.8, label='all neurons'),
                    plt.Rectangle((0, 0), 1, 1, fc='k', ec='k', lw=0.8, label='top 5% by synapses')], loc='upper right', fontsize=6.5, borderaxespad=0.1)
fig.canvas.draw(); panel_labels(fig, [axA, axB], 'ab')
gate(fig, 'FigS6', axes=[axA, axB], panel_ids=['a', 'b']); export(fig, 'FigS6'); plt.close(fig)
