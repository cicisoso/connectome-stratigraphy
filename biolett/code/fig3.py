"""Figure 3: hemilineages without clonal birth dating. Primary (embryonic-born) versus secondary (larval-born) neurons
in every hemilineage unit of three connectomes, and their connectivity in four."""
from figstyle import *
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm

U = {ds: pd.read_csv(f'{RES}/hemilineage_units_{ds}.csv') for ds in DS}
CU = {ds: pd.read_csv(f'{RES}/connectivity_units_{ds}.csv') for ds in DS + ['MANC']}
HG = [('MaleCNS', 'brain'), ('FlyWire', 'brain'), ('BANC', 'brain'), ('MaleCNS', 'vnc'), ('BANC', 'vnc')]
CG = HG + [('MANC', 'vnc')]
REG = {'brain': 'brain', 'vnc': 'VNC'}; MK = {'brain': 'o', 'vnc': 's'}

fig = plt.figure(figsize=(174 * MM, 128 * MM))
gs = gridspec.GridSpec(2, 2, width_ratios=[1.45, 1], height_ratios=[1, 1], left=0.075, right=0.885, top=0.925, bottom=0.12, wspace=0.42, hspace=0.66)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1]); axC = fig.add_subplot(gs[1, 0]); axD = fig.add_subplot(gs[1, 1])

def med_iqr(ax, x, v, w=0.32):
    v = np.asarray(v, float); v = v[np.isfinite(v)]
    if not len(v): return
    q1, m, q3 = np.percentile(v, [25, 50, 75])
    ax.plot([x, x], [q1, q3], color='k', lw=1.6, solid_capstyle='butt', zorder=4); ax.plot([x - w, x + w], [m, m], color='k', lw=1.2, zorder=4)

def group_labels(ax, groups, centres, y=-0.06):
    ax.set_xticks(centres); ax.set_xticklabels([''] * len(centres)); ax.tick_params(axis='x', length=0)
    for (ds, rg), c in zip(groups, centres):
        ax.text(c, y, ds, transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=6.5, color=DSC[ds], fontweight='bold')
        ax.text(c, y - 0.1, REG[rg], transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=6.5)

# (a) per-unit AUC, primary neurons look older: soma index (not an annotation criterion) and fibre radial position (one)
cent = []
for g, (ds, rg) in enumerate(HG):
    T = U[ds][U[ds].region == rg]; x0 = 3 * g
    for k, (col, hollow) in enumerate([('auc_geo_age_soma', False), ('auc_radial', True)]):
        v = T[col].dropna(); strip(ax=axA, values=v, x=x0 + k, color=DSC[ds], width=0.55, s=5 if rg == 'vnc' else 9, alpha=0.8,
                                     marker=MK[rg], facecolor='white' if hollow else None, seed=g * 10 + k)
        med_iqr(axA, x0 + k, v)
        axA.text(x0 + k, 1.02 + 0.075 * k, f'{int((v > 0.5).sum())}/{len(v)}', transform=axA.get_xaxis_transform(), ha='center', va='bottom', fontsize=6,
                 color='k' if k == 0 else GREY)
    cent.append(x0 + 0.5)
axA.axhline(0.5, color=GREY, lw=0.6, ls=(0, (2, 2)), zorder=0); axA.set_ylim(-0.32, 1.0); axA.set_xlim(-0.8, 3 * len(HG) - 1.2)
axA.set_yticks(np.arange(0, 1.01, 0.2)); axA.spines['left'].set_bounds(0, 1)
axA.set_ylabel('AUC, primary neurons\nlook older')
group_labels(axA, HG, cent)
axA.legend(handles=[Line2D([], [], marker='o', ls='', color='k', ms=3.5, label='soma index (soma depth and distance)'),
                    Line2D([], [], marker='o', ls='', mfc='white', mec='k', ms=3.5, label='fibre radial position (annotation criterion)')],
           loc='lower left', bbox_to_anchor=(0.0, 0.0), fontsize=6.5, handletextpad=0.1, borderaxespad=0.1, ncol=1)

# (b) the same rule within the secondary lineage: early versus late secondary neurons (VNC, MANC-derived birth classes)
XB = [0, 1.35, 3.75, 5.1]
for g, ds in enumerate(['MaleCNS', 'BANC']):
    T = U[ds][U[ds].region == 'vnc']
    for k, (col, lo) in enumerate([('auc_early_vs_late_geo_age_soma', 0.5), ('rho_level_geo_age_soma', 0.0)]):
        x = XB[2 * k + g]; v = T[col].dropna()
        strip(axB, v, x, DSC[ds], width=0.6, s=5, alpha=0.8, marker='s', seed=40 + g + k); med_iqr(axB, x, v)
        axB.text(x, 1.02, f'{int((v > lo).sum())}/{len(v)}', transform=axB.get_xaxis_transform(), ha='center', va='bottom', fontsize=6)
        axB.text(x, -0.03, ds, transform=axB.get_xaxis_transform(), ha='center', va='top', fontsize=6.5, color=DSC[ds], fontweight='bold')
axB.set_xlim(-0.8, 5.9); axB.set_ylim(-0.6, 1.0); axB.set_xticks(XB); axB.set_xticklabels([''] * 4); axB.tick_params(axis='x', length=0)
axB.plot([-0.8, 2.15], [0.5, 0.5], color=GREY, lw=0.6, ls=(0, (2, 2)), zorder=0); axB.plot([2.95, 5.9], [0, 0], color=GREY, lw=0.6, ls=(0, (2, 2)), zorder=0)
axB.text(0.675, -0.15, 'AUC, early secondary\nlook older than late', transform=axB.get_xaxis_transform(), ha='center', va='top', fontsize=6.5)
axB.text(4.425, -0.15, 'ρ, primary > early >\nlate secondary', transform=axB.get_xaxis_transform(), ha='center', va='top', fontsize=6.5)
axB.set_ylabel('Soma index statistic,\nper VNC unit')

# (c) connectivity: log2 ratio of primary to secondary medians, per hemilineage unit
MEAS = [('log2ratio_n_syn', 'synapses'), ('log2ratio_size', 'size'), ('log2ratio_syn_per_size', 'synapses per size')]
xt = []
for m, (col, lab) in enumerate(MEAS):
    for g, (ds, rg) in enumerate(CG):
        T = CU[ds][CU[ds].region == rg]; x = m * (len(CG) + 1.5) + g; v = T[col].dropna()
        strip(axC, v, x, DSC[ds], width=0.6, s=4 if rg == 'vnc' else 8, alpha=0.75, marker=MK[rg], seed=70 + 7 * m + g); med_iqr(axC, x, v, w=0.36)
    xt.append(m * (len(CG) + 1.5) + (len(CG) - 1) / 2)
axC.axhline(0, color=GREY, lw=0.6, ls=(0, (2, 2)), zorder=0)
axC.set_xticks(xt); axC.set_xticklabels([l for _, l in MEAS]); axC.tick_params(axis='x', length=0)
axC.set_ylabel('log$_2$ (primary / secondary),\nunit medians'); axC.set_ylim(-2.5, 6.5); axC.set_xlim(-1, 3 * (len(CG) + 1.5) - 1.5)
axC.legend(handles=[Line2D([], [], marker='s', ls='', color=DSC[ds], ms=3.5, label=ds) for ds in DS + ['MANC']] +
                   [Line2D([], [], marker=MK[r], ls='', mfc='white', mec='k', ms=3.5, label=REG[r]) for r in ['brain', 'vnc']],
           loc='upper right', ncol=3, fontsize=6.5, handletextpad=0.1, columnspacing=0.6, borderaxespad=0.1)

# (d) partner classes: difference in the fraction of input or output synapses (primary minus secondary), unit medians
CLS = ['intrinsic', 'sensory', 'ascending', 'descending', 'visual', 'motor/efferent']
rows = [(d, c) for d in ['in', 'out'] for c in CLS]
H = np.array([[CU[ds][CU[ds].region == rg][f'diff_{d}_{c}'].median() for (ds, rg) in CG] for d, c in rows])
keep = np.abs(H).max(1) >= 0.0015   # rows whose medians all lie below 0.15 percentage points are left out
rows = [r for r, k in zip(rows, keep) if k]; H = H[keep]; n_in = sum(d == 'in' for d, _ in rows)
im = axD.imshow(H * 100, cmap='PuOr_r', norm=TwoSlopeNorm(0, -8, 8), aspect='auto')
for i in range(H.shape[0]):
    for j in range(H.shape[1]):
        v = H[i, j] * 100
        axD.text(j, i, f'{v:+.1f}'.replace('-', '−') if abs(v) >= 0.05 else '0', ha='center', va='center', fontsize=5.5, color='white' if abs(v) > 5.5 else 'k')
axD.set_yticks(range(len(rows))); axD.set_yticklabels([f'{"input" if d == "in" else "output"}: {c}' for d, c in rows], fontsize=6.5)
axD.set_xticks(range(len(CG))); axD.set_xticklabels([f'{ds} {REG[rg]}' for ds, rg in CG], rotation=40, ha='right', fontsize=6.5)
for t, (ds, rg) in zip(axD.get_xticklabels(), CG): t.set_color(DSC[ds])
axD.axhline(n_in - 0.5, color='k', lw=0.6)
for s in axD.spines.values(): s.set_visible(False)
axD.tick_params(length=0)
cax = axD.inset_axes([1.03, 0.0, 0.04, 1.0]); cb = fig.colorbar(im, cax=cax); cb.outline.set_linewidth(0.4)
cb.ax.tick_params(labelsize=6.5, length=1.5); cb.set_label('primary − secondary (percentage points)', fontsize=6.5, labelpad=2)

fig.canvas.draw()
panel_labels(fig, [axA, axB, axC, axD], 'abcd')
gate(fig, 'Fig3', axes=[axA, axC], panel_ids=['a', 'c'])
export(fig, 'Fig3')
