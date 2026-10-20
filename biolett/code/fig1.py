"""Figure 1: the site-of-addition rule and its positive control, the Kenyon cells of three connectomes."""
from figstyle import *
import matplotlib.gridspec as gridspec
from matplotlib.patches import Circle, FancyBboxPatch, Wedge, FancyArrowPatch
from matplotlib.lines import Line2D
from scipy import stats

MAJOR = {0: 0, 1: 0, 2: 1, 3: 2, 4: 2, 5: 2, 6: 2}
KC = {ds: pd.read_parquet(f'{RES}/kc_{ds}.parquet').assign(major=lambda t: t['rank'].map(MAJOR)) for ds in DS}
SUM = pd.read_csv(f'{RES}/summary_kc.csv')

fig = plt.figure(figsize=(174 * MM, 118 * MM))
gs = gridspec.GridSpec(2, 2, width_ratios=[1.55, 1], height_ratios=[1, 1], left=0.075, right=0.985, top=0.955, bottom=0.115, wspace=0.22, hspace=0.5)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
gs2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, :], wspace=0.3)
axC = fig.add_subplot(gs2[0, 0]); axD = fig.add_subplot(gs2[0, 1])

# (a) schematic -------------------------------------------------------------------------------------------------
ax = axA; ax.set_xlim(0, 15.5); ax.set_ylim(0, 8); ax.set_aspect('equal'); ax.axis('off')
ax.add_patch(FancyBboxPatch((0.3, 0.25), 8.6, 1.5, boxstyle='round,pad=0.05,rounding_size=0.5', fc=LIGHT, ec='none'))
ax.text(0.6, 0.55, 'neuropil', fontsize=7, color=GREY, ha='left', va='center')
ax.plot([0.3, 8.9], [7.55, 7.55], color=GREY, lw=0.6, ls=(0, (2, 2))); ax.text(0.35, 7.7, 'brain surface', fontsize=7, color=GREY, va='bottom')
ax.add_patch(Circle((4.6, 7.0), 0.52, fc='white', ec='k', lw=0.9)); ax.text(5.3, 7.05, 'neuroblast', fontsize=7, va='center')
rng = np.random.default_rng(3); n = 22; somata = []
for k in range(n):   # k = birth order; older cells displaced towards the neuropil and to the sides
    depth = 6.15 - 3.4 * (1 - k / (n - 1)) + rng.normal(0, 0.12); spread = 0.35 + 1.55 * (1 - k / (n - 1))
    x = 4.6 + rng.uniform(-spread, spread); somata.append((x, depth, k))
entry = np.array([4.6, 1.9])
for x, y, k in sorted(somata, key=lambda s: s[2]):
    c = RANK_CMAP(k / (n - 1)); ax.plot([x, 4.6 + 0.25 * (x - 4.6) / 2, entry[0] + 0.08 * (x - 4.6)], [y, 2.6, entry[1]], color=c, lw=0.55, alpha=0.9, zorder=1)
for x, y, k in somata: ax.add_patch(Circle((x, y), 0.24, fc=RANK_CMAP(k / (n - 1)), ec='white', lw=0.4, zorder=3))
ax.annotate('', xy=(2.05, 6.1), xytext=(2.05, 2.95), arrowprops=dict(arrowstyle='->', lw=0.8, color='k'))
ax.text(1.85, 4.5, 'soma depth', rotation=90, fontsize=7, ha='right', va='center')
# tract cross-section: concentric layers, youngest in the core
cx, cy = 12.3, 4.4
ax.annotate('', xy=(cx - 2.2, cy - 0.2), xytext=(4.95, 2.45), arrowprops=dict(arrowstyle='->', lw=0.6, color=GREY, connectionstyle='arc3,rad=-0.25'))
for j, r in enumerate(np.linspace(2.0, 0.35, 6)): ax.add_patch(Circle((cx, cy), r, fc=RANK_CMAP(j / 5), ec='white', lw=0.5))
ax.text(cx, cy + 2.2, 'tract cross-section', fontsize=7, ha='center', va='bottom')
ax.annotate('', xy=(cx + 2.0, cy - 2.3), xytext=(cx, cy - 2.3), arrowprops=dict(arrowstyle='->', lw=0.8))
ax.text(cx + 1.0, cy - 2.45, 'fibre radial position', fontsize=7, ha='center', va='top')
# birth-order key
for j in range(40): ax.add_patch(plt.Rectangle((11.0 + j * 0.07, 0.35), 0.07, 0.3, fc=RANK_CMAP(j / 39), ec='none'))
ax.text(10.9, 0.5, 'early-born', fontsize=7, ha='right', va='center'); ax.text(13.9, 0.5, 'late-born', fontsize=7, ha='left', va='center')

# (b) Kenyon cell peduncle cross-section (MaleCNS, left) ----------------------------------------------------------
sec = pd.read_parquet(f'{RES}/kc_section_MaleCNS.parquet'); sec = sec[sec.side == 'L'].merge(KC['MaleCNS'][['nid', 'major']], on='nid')
for k in [0, 1, 2]:
    q = sec[sec.major == k]; axB.scatter(q.u, q.v, s=1.6, color=KC_MAJOR_COL[k], lw=0, rasterized=True, label=KC_MAJOR_LAB[k])
axB.set_aspect('equal', adjustable='datalim'); clean_3d_axes(axB)
x0 = sec.u.min(); y0 = sec.v.min() - 2.5; scalebar(axB, x0, y0, 10, '10 µm', dy=-0.8)
axB.set_title('Kenyon cell peduncle, MaleCNS (left)', fontsize=7)
axB.legend(handles=[Line2D([], [], marker='o', ls='', color=KC_MAJOR_COL[k], markersize=4, label=KC_MAJOR_LAB[k]) for k in [0, 1, 2]],
           loc='lower right', ncol=3, handletextpad=0.1, columnspacing=0.6, borderaxespad=0.0)

# (c,d) class medians in three connectomes ------------------------------------------------------------------------
def grouped(ax, col, ylabel, sign):
    xt, xl = [], []
    for g, ds in enumerate(DS):
        T = KC[ds].dropna(subset=[col]); base = g * 4
        boxes(ax, [T[T.major == k][col] for k in [0, 1, 2]], [base + k for k in [0, 1, 2]], [KC_MAJOR_COL[k] for k in [0, 1, 2]], width=0.7)
        s = SUM[(SUM.ds == ds) & (SUM.measure == col)].iloc[0]
        ax.text(base + 1, 1.02, f'ρ = {s.rho_major:.2f}\nC = {s.C:.2f}'.replace('-', '−'), transform=ax.get_xaxis_transform(), ha='center', va='bottom', fontsize=7)
        xt += [base + k for k in [0, 1, 2]]; xl += [KC_MAJOR_LAB[k] for k in [0, 1, 2]]
        ax.text(base + 1, -0.13, ds, transform=ax.get_xaxis_transform(), ha='center', va='top', fontsize=7, color=DSC[ds], fontweight='bold')
    ax.set_xticks(xt); ax.set_xticklabels(xl); ax.set_ylabel(ylabel); ax.set_xlim(-0.7, 10.7)
grouped(axC, 'ped_radial', 'Fibre radial position\n(fraction of peduncle radius)', -1)
grouped(axD, 'soma_tang_clone', 'Soma distance from\nclone axis (µm)', -1)
axD.text(9, 0.93, 'clones by\nk-means', transform=axD.get_xaxis_transform(), ha='center', va='top', fontsize=6.5, color=GREY)

fig.canvas.draw()
panel_labels(fig, [axA, axB, axC, axD], 'abcd')
gate(fig, 'Fig1', axes=[axC, axD], panel_ids=['c', 'd'])
export(fig, 'Fig1')
