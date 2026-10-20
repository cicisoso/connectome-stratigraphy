"""Electronic supplementary figures S1 and S2: renderings of every clonally dated lineage frame in the three connectomes.
S1: soma clusters and cell body fibres projected onto the plane of the soma depth axis; S2: tract cross-sections."""
from figstyle import *
import matplotlib.gridspec as gridspec

LIN = [('adPN', 'ALad1', 'L'), ('lPN', 'ALl1_ventral', 'L'), ('lPN, PNs only', 'ALl1_ventral', 'L|PN'), ('vPN', 'ALv1', 'L'),
       ('PAM', 'CREa1_dorsal', 'L'), ('PAM', 'CREa2_medial', 'L')]

def load(ds):
    S = pd.read_parquet(f'{RES}/render_somata_{ds}.parquet'); PA = pd.read_parquet(f'{RES}/render_paths_{ds}.parquet')
    FR = {(r['hl'], r['side']): r for r in json.load(open(f'{RES}/render_frames_{ds}.json'))}
    return S, PA, FR
DATA = {ds: load(ds) for ds in DS}

def norm_rank(r, hl):
    if hl == 'ALad1': return r / 39
    if hl == 'ALl1_ventral': return r / 12
    return (r - 12) / (36 - 12)   # mean clone-induction window (index) for ALv1 and PAM types

def projection(ax, ds, hl, side):
    S, PA, FR = DATA[ds]; fr = FR.get((hl, side))
    s = S[(S.frame_hl == hl) & (S.frame_side == side)]; pa = PA[(PA.frame_hl == hl) & (PA.frame_side == side)]
    if fr is None or not len(s): ax.axis('off'); return
    cs = np.array(fr['soma_centroid']); axv = np.array(fr['soma_axis']); oc = np.array(fr['outflow_centroid'])
    X = s[['x', 'y', 'z']].values - cs; perp = X - np.outer(X @ axv, axv); e = np.linalg.svd(perp, full_matrices=False)[2][0]
    pr = lambda Q: np.c_[(Q - cs) @ e, (Q - cs) @ axv]
    for _, g in pa.groupby('nid'):
        q = pr(g[['x', 'y', 'z']].values); ax.plot(q[:, 0], q[:, 1], color='#CFCFCF', lw=0.25, zorder=1)
    q = pr(s[['x', 'y', 'z']].values); rk = s['rank']
    ax.scatter(q[rk.isna(), 0], q[rk.isna(), 1], s=3, color='white', edgecolors=GREY, lw=0.3, zorder=2)
    ax.scatter(q[rk.notna(), 0], q[rk.notna(), 1], s=4, c=norm_rank(rk[rk.notna()], hl), cmap=RANK_CMAP, vmin=0, vmax=1, lw=0, zorder=3)
    o = pr(oc[None])[0]; ax.plot(*o, marker='s', color='k', ms=2.5, zorder=4); ax.plot(0, 0, marker='+', color='k', ms=5, mew=0.8, zorder=4)
    if 'tract' in fr:
        cp = pr(np.array(fr['tract']['convergence_point'])[None])[0]; ax.plot(*cp, marker='*', color=PRIM, ms=5, zorder=4)
    ax.annotate('', xy=(0, 15), xytext=(o[0], o[1]), arrowprops=dict(arrowstyle='->', lw=0.6))
    ax.set_aspect('equal', adjustable='datalim'); clean_3d_axes(ax)
    xl = ax.get_xlim(); yl = ax.get_ylim(); ax.plot([xl[0], xl[0] + 20], [yl[0], yl[0]], color='k', lw=0.9)

def section(ax, ds, hl, side):
    S, PA, FR = DATA[ds]; s = S[(S.frame_hl == hl) & (S.frame_side == side) & S.u.notna()]
    if not len(s): ax.axis('off'); ax.text(0.5, 0.5, 'no tract section', transform=ax.transAxes, ha='center', fontsize=6.5, color=GREY); return
    rk = s['rank']
    ax.scatter(s.u[rk.isna()], s.v[rk.isna()], s=4, color='white', edgecolors=GREY, lw=0.3)
    ax.scatter(s.u[rk.notna()], s.v[rk.notna()], s=6, c=norm_rank(rk[rk.notna()], hl), cmap=RANK_CMAP, vmin=0, vmax=1, edgecolors='k', lw=0.15)
    ax.plot(0, 0, marker='+', color='k', ms=5, mew=0.8)
    ax.set_aspect('equal', adjustable='datalim'); clean_3d_axes(ax)
    xl = ax.get_xlim(); yl = ax.get_ylim(); ax.plot([xl[0], xl[0] + 2], [yl[0], yl[0]], color='k', lw=0.9)

for name, fn, note in [('FigS1', projection, 'scale bar 20 µm'), ('FigS2', section, 'scale bar 2 µm')]:
    fig = plt.figure(figsize=(174 * MM, 228 * MM))
    gs = gridspec.GridSpec(len(LIN), 3, left=0.1, right=0.99, top=0.945, bottom=0.02, wspace=0.08, hspace=0.1)
    axes = []
    for i, (lin, hl, side) in enumerate(LIN):
        for j, ds in enumerate(DS):
            ax = fig.add_subplot(gs[i, j]); fn(ax, ds, hl, side); axes.append(ax)
            if i == 0: ax.set_title(ds, fontsize=8, color=DSC[ds], fontweight='bold')
            if j == 0: ax.text(-0.06, 0.5, f'{lin}\n{hl}', transform=ax.transAxes, rotation=90, ha='right', va='center', fontsize=7)
    fig.text(0.1, 0.988, f'Left-side frames; colour = birth rank or window (early dark, late yellow); open circles = types without a mapped birth rank; {note}.', ha='left', va='top', fontsize=6.5, color=GREY)
    if name == 'FigS1': fig.text(0.1, 0.976, 'Square = fibre outflow centroid; + = soma centroid; arrow = soma depth axis; star = tract convergence point.', ha='left', va='top', fontsize=6.5, color=GREY)
    export(fig, name)
