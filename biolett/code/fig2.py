"""Figure 2: four lineages whose birth order was mapped by clonal analysis, in three connectomes."""
from figstyle import *
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy import stats

L = pd.read_csv(f'{RES}/summary_lineages.csv'); L = L[L.frame == 'hemilineage']
P = pd.read_csv(f'{RES}/summary_lineages_pooled.csv'); P = P[P.frame == 'hemilineage']

fig = plt.figure(figsize=(174 * MM, 132 * MM))
gs = gridspec.GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1.08], left=0.165, right=0.985, top=0.955, bottom=0.075, wspace=0.36, hspace=0.3)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1]); axC = fig.add_subplot(gs[1, 0]); axD = fig.add_subplot(gs[1, 1])

# (a,b) renderings: somata coloured by birth rank, cell body fibres (0-60 µm), soma depth axis and tract convergence point
def render(ax, ds, hl, side, norm, title, cbar_ticks, cbar_title, cbx=0.72):
    S = pd.read_parquet(f'{RES}/render_somata_{ds}.parquet'); S = S[(S.frame_hl == hl) & (S.frame_side == side)]
    PA = pd.read_parquet(f'{RES}/render_paths_{ds}.parquet'); PA = PA[(PA.frame_hl == hl) & (PA.frame_side == side)]
    fr = [r for r in json.load(open(f'{RES}/render_frames_{ds}.json')) if r['hl'] == hl and r['side'] == side][0]
    cs = np.array(fr['soma_centroid']); axv = np.array(fr['soma_axis']); oc = np.array(fr['outflow_centroid'])
    X = S[['x', 'y', 'z']].values - cs; perp = X - np.outer(X @ axv, axv); e = np.linalg.svd(perp, full_matrices=False)[2][0]
    proj = lambda Q: np.c_[(Q - cs) @ e, (Q - cs) @ axv]
    for nid, g in PA.groupby('nid'):
        q = proj(g[['x', 'y', 'z']].values); ax.plot(q[:, 0], q[:, 1], color='#CDCDCD', lw=0.3, zorder=1)
    rk = S['rank']; q = proj(S[['x', 'y', 'z']].values)
    ax.scatter(q[rk.isna(), 0], q[rk.isna(), 1], s=6, color='white', edgecolors=GREY, lw=0.4, zorder=2)
    sc = ax.scatter(q[rk.notna(), 0], q[rk.notna(), 1], s=9, c=norm(rk[rk.notna()]), cmap=RANK_CMAP, vmin=0, vmax=1, edgecolors='k', lw=0.2, zorder=3)
    o = proj(oc[None])[0]; cp = proj(np.array(fr['tract']['convergence_point'])[None])[0]
    ax.annotate('', xy=(0, 18), xytext=(o[0], o[1]), arrowprops=dict(arrowstyle='->', lw=0.9, color='k'), zorder=4)
    ax.plot(*o, marker='s', color='k', ms=3.5, zorder=5); ax.plot(0, 0, marker='+', color='k', ms=6, mew=1.0, zorder=5)
    ax.plot(*cp, marker='*', color=PRIM, ms=6, zorder=5)
    ax.set_aspect('equal', adjustable='datalim'); clean_3d_axes(ax)
    xl = ax.get_xlim(); yl = ax.get_ylim(); scalebar(ax, xl[0] + 2, yl[0] + 2, 20, '20 µm', dy=-1.5)
    ax.set_title(title, fontsize=7)
    cax = ax.inset_axes([cbx, 0.06, 0.26, 0.035]); cb = fig.colorbar(sc, cax=cax, orientation='horizontal'); cb.set_ticks([0, 1]); cb.set_ticklabels(cbar_ticks)
    cb.ax.tick_params(labelsize=6.5, length=1.5, pad=1); cb.outline.set_linewidth(0.4); cax.set_title(cbar_title, fontsize=6.5, pad=1.5)
render(axA, 'MaleCNS', 'ALad1', 'R', lambda r: r / 39, 'adPN lineage (ALad1), MaleCNS', ['first', 'last'], 'birth rank (Yu 2010)')
render(axB, 'FlyWire', 'CREa1_dorsal', 'L', lambda r: (r - 12) / 24, 'PAM lineage (CREa1 dorsal), FlyWire', ['early', 'late'], 'birth window (Lee 2020)', cbx=0.62)
axA.legend(handles=[Line2D([], [], marker='s', ls='', color='k', ms=3.5, label='fibre outflow centroid'),
                    Line2D([], [], marker='+', ls='', color='k', ms=6, mew=1.0, label='soma centroid'),
                    Line2D([], [], marker='*', ls='', color=PRIM, ms=6, label='tract convergence point'),
                    Line2D([], [], marker='o', ls='', mfc='white', mec=GREY, ms=3, label='unranked type')],
           loc='upper left', bbox_to_anchor=(-0.02, 1.0), fontsize=6.5, handletextpad=0.2, borderaxespad=0.0, labelspacing=0.3)

# (c) forest plot of the geometric age index, with its components
LIN = ['adPN', 'lPN', 'vPN', 'PAM']; LLAB = {'adPN': 'adPN (Yu 2010)', 'lPN': 'lPN (Lin 2012)', 'vPN': 'vPN (Lee 2020)', 'PAM': 'PAM (Lee 2020)'}
COMP = [('soma_depth', 1, '^'), ('soma_dist', -1, 's'), ('radial', -1, 'D')]
y = 0; yt, yl_ = [], []
for lin in LIN:
    p = P[(P.lineage == lin) & (P.measure == 'geo_age')].iloc[0]
    axC.plot([p.ci_lo, p.ci_hi], [y, y], color='k', lw=1.0); axC.plot(p.pooled_rho, y, marker='D', color='k', ms=5)
    yt.append(y); yl_.append(f'{LLAB[lin]}'); y += 1
    for ds in DS:
        r = L[(L.lineage == lin) & (L.ds == ds) & (L.measure == 'geo_age')].iloc[0]
        axC.plot([r.ci_lo, r.ci_hi], [y, y], color=DSC[ds], lw=0.9); axC.plot(r.rho, y, marker='o', color=DSC[ds], ms=4)
        for k, (m, sg, mk) in enumerate(COMP):
            rr = L[(L.lineage == lin) & (L.ds == ds) & (L.measure == m)]
            if len(rr) and rr.rho.iloc[0] == rr.rho.iloc[0]:
                axC.plot(sg * rr.rho.iloc[0], y + 0.28, marker=mk, ms=2.6, mfc='white', mec=DSC[ds], mew=0.6, ls='')
        axC.text(1.02, y, f'{r.C:.2f}', transform=axC.get_yaxis_transform(), fontsize=6.5, va='center', ha='left', color=DSC[ds])
        yt.append(y); yl_.append(f'{ds}  n = {int(r.n)}'); y += 1
    y += 0.6
axC.axvline(0, color=GREY, lw=0.6, ls=(0, (2, 2)))
axC.set_yticks(yt); axC.set_yticklabels(yl_, fontsize=6.5); axC.invert_yaxis(); axC.set_xlim(-0.5, 1.0)
axC.set_xlabel('Spearman ρ with birth rank (predicted direction > 0)')
axC.text(1.02, -1.2, 'C', transform=axC.get_yaxis_transform(), fontsize=6.5, ha='left', va='center', fontweight='bold')
axC.legend(handles=[Line2D([], [], marker='o', ls='', color='k', ms=4, label='geometric age'),
                    Line2D([], [], marker='^', ls='', mfc='white', mec='k', ms=3, label='soma depth'),
                    Line2D([], [], marker='s', ls='', mfc='white', mec='k', ms=3, label='−soma distance'),
                    Line2D([], [], marker='D', ls='', mfc='white', mec='k', ms=3, label='−fibre radial')],
           loc='lower right', fontsize=6.5, handletextpad=0.2, borderaxespad=0.2)

# (d) adPN type medians of the geometric age index against birth rank
for ds in DS:
    M = pd.read_parquet(f'{RES}/lineages_hemilineage_{ds}.parquet'); M = M[(M.lineage == 'adPN') & M['rank'].notna()]
    tm = M.groupby(['type', 'rank']).geo_age.median().reset_index()
    axD.scatter(tm['rank'], tm.geo_age, s=9, color=DSC[ds], lw=0, alpha=0.85, label=f'{ds}')
    b = np.polyfit(tm['rank'], tm.geo_age, 1); xx = np.array([0, 39]); axD.plot(xx, np.polyval(b, xx), color=DSC[ds], lw=0.9)
axD.axvspan(-1, 17.5, color=LIGHT, zorder=0, lw=0); axD.text(8.5, 0.98, 'embryonic-born', transform=axD.get_xaxis_transform(), ha='center', va='top', fontsize=7)
axD.text(28.5, 0.98, 'larval-born', transform=axD.get_xaxis_transform(), ha='center', va='top', fontsize=7)
axD.set_xlim(-1, 40); axD.set_xlabel('adPN birth rank (Yu et al. 2010)'); axD.set_ylabel('Geometric age index\n(type median, z units)')
r_t = [L[(L.lineage == 'adPN') & (L.ds == ds) & (L.measure == 'geo_age')].type_rho.iloc[0] for ds in DS]
axD.legend(handles=[Line2D([], [], marker='o', ls='', color=DSC[ds], ms=4, label=f'{ds}  ρ(types) = {r:.2f}') for ds, r in zip(DS, r_t)],
           loc='lower right', fontsize=6.5, handletextpad=0.2, borderaxespad=0.2)

fig.canvas.draw()
panel_labels(fig, [axA, axB, axC, axD], 'abcd')
gate(fig, 'Fig2', axes=[axC, axD], panel_ids=['c', 'd'])
export(fig, 'Fig2')
