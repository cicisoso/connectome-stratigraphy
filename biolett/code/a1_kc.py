"""Kenyon cells: fibre position in the peduncle and the alpha/beta lobes, and soma position in the clone, against birth-order class.

Classes in birth order (Lee et al. 1999; Tanaka et al. 2008; Li et al. 2020): gamma-d, gamma-m, alpha'/beta', alpha/beta-p,
alpha/beta-s, alpha/beta-m, alpha/beta-c. Usage: python a1_kc.py MaleCNS|FlyWire|BANC
"""
import sys
from common import *
from sklearn.cluster import KMeans

RANK = {'KCg-d': 0, 'KCg-m': 1, "KCa'b'-ap1": 2, "KCa'b'-ap2": 2, "KCa'b'-m": 2, 'KCapbp-ap1': 2, 'KCapbp-ap2': 2, 'KCapbp-m': 2,
        "KCa'b'": 2, 'KCab-p': 3, 'KCab-s': 4, 'KCab-m': 5, 'KCab-c': 6}
LAB = {0: 'γd', 1: 'γm', 2: "α′β′", 3: 'αβp', 4: 'αβs', 5: 'αβm', 6: 'αβc'}

def kc_table(ds):
    d = LOADERS[ds]()
    if ds == 'MaleCNS':
        k = d[d.cls == 'Kenyon_Cell'].copy(); k['kctype'] = k.type
        a = pd.read_feather(f'{DATA}/body-annotations-male-cns-v1.0-minconf-0.5.feather', columns=['bodyId', 'itoleeHl']).set_index('bodyId').itoleeHl
        k['clone'] = k.nid.map(a)
    elif ds == 'FlyWire':
        k = d[d.cls == 'Kenyon_Cell'].copy(); k['kctype'] = k.hemibrain_type
        f = pd.read_csv(f'{DATA}/flywire/annotations/Supplemental_file1_neuron_annotations.tsv', sep='\t', usecols=['root_id', 'ito_lee_hemilineage']).set_index('root_id').ito_lee_hemilineage
        k['clone'] = k.nid.map(f)
    else:
        k = d[d.cls == 'kenyon_cell'].copy(); k['kctype'] = k.type
        k['clone'] = None
    k['rank'] = k.kctype.map(RANK)
    if ds == 'BANC':   # alpha/beta subclasses are mostly untyped in BANC: use the four major classes
        k.loc[k.kctype.isin(['KCab', 'KCab-p', 'KCab-s', 'KCab-m', 'KCab-c']), 'rank'] = 3
    return k[k['rank'].notna() & k.side.isin(['L', 'R'])].copy()

def side_branch_cable(sk, u, p_prev, p_next):
    seen = {u, p_prev, p_next}; tot = 0.0; q = deque([v for v in sk.adj[u] if v not in (p_prev, p_next)])
    for v in q: seen.add(v)
    stack = list(q); tot = sum(np.linalg.norm(sk.xyz[u] - sk.xyz[v]) for v in stack)
    while stack:
        w = stack.pop()
        for v in sk.adj[w]:
            if v not in seen: seen.add(v); tot += np.linalg.norm(sk.xyz[w] - sk.xyz[v]); stack.append(v)
    return tot

def peduncle(sk, soma):
    """Soma-to-tip path resampled every 1 µm and its longest branch-free stretch (side branches > 4 µm of cable count as
    branches). Returns (path points, stretch start, stretch end) or None."""
    r = soma_path(sk, soma)
    if r is None: return None
    path, gd = r; L = gd[-1]
    if L < 150: return None
    cuts = [0.0]
    for j in range(1, len(path) - 1):
        if len(sk.adj[path[j]]) > 2 and side_branch_cable(sk, path[j], path[j - 1], path[j + 1]) > 4: cuts.append(gd[j])
    cuts.append(L); gaps = [(cuts[i + 1] - cuts[i], cuts[i], cuts[i + 1]) for i in range(len(cuts) - 1)]
    g, s0, s1 = max(gaps)
    if g < 40: return None
    P = sk.xyz[path]; t = np.arange(0, L + 1e-9, 1.0)
    Q = np.c_[np.interp(t, gd, P[:, 0]), np.interp(t, gd, P[:, 1]), np.interp(t, gd, P[:, 2])]
    return Q, s0, s1

def peduncle_cylinder(paths, stretch, half=20.0):
    """Population peduncle: centre = geometric median of the points 25 µm into each fibre's longest branch-free stretch,
    axis = principal axis of the proximal stretches (5-45 µm) that start near that centre, radius = 99th percentile of their
    distance from the axis. Each fibre contributes its longest contiguous run of path points inside the cylinder."""
    ids = list(paths); mid = np.array([paths[i][min(int(stretch[i][0]) + 25, len(paths[i]) - 1)] for i in ids])
    c0 = geometric_median(mid); near = [i for i, m in zip(ids, mid) if np.linalg.norm(m - c0) < 20]
    X = np.vstack([paths[i][int(stretch[i][0]) + 5:int(stretch[i][0]) + 46] for i in near])
    a = np.linalg.svd(X - X.mean(0), full_matrices=False)[2][0]; c0 = X.mean(0)
    perp = np.linalg.norm((X - c0) - np.outer((X - c0) @ a, a), axis=1); R = np.quantile(perp, 0.99)
    pts = {}
    for i in ids:
        Q = paths[i]; t = (Q - c0) @ a; inside = (np.abs(t) <= half) & (np.linalg.norm((Q - c0) - np.outer(t, a), axis=1) <= R)
        if not inside.any(): continue
        k = np.flatnonzero(np.diff(np.r_[0, inside.astype(int), 0]))
        runs = k.reshape(-1, 2); best = runs[np.argmax(runs[:, 1] - runs[:, 0])]
        if best[1] - best[0] >= 1.5 * half: pts[i] = Q[best[0]:best[1]]
    return pts, dict(centre=c0.round(3).tolist(), axis=a.round(5).tolist(), radius_um=float(R), half_length_um=half,
                     n_near=len(near))

def run(ds):
    K = kc_table(ds); print(ds, 'KCs with class and side', len(K), flush=True)
    rows, paths, stretch, lobe_pts = [], {}, {}, {}
    for r in K.itertuples():
        sk = skeleton(ds, r.nid, getattr(r, 'skel', None))
        if sk is None: continue
        out = peduncle(sk, np.array([r.x, r.y, r.z]))
        if out is None: continue
        Q, s0, s1 = out; paths[r.nid] = Q; stretch[r.nid] = (s0, s1)
        if len(Q) > s1 + 20: lobe_pts[r.nid] = Q[int(s1) + 8:int(min(s1 + 50, len(Q) - 1)) + 1]
        rows.append(dict(nid=r.nid, side=r.side, kctype=r.kctype, rank=int(r.rank), clone=r.clone, x=r.x, y=r.y, z=r.z, ped_len=s1 - s0))
    T = pd.DataFrame(rows); print(ds, 'with soma path and branch-free stretch', len(T), flush=True)
    axes = []; T['ped_radial'] = np.nan; T['ped_burial'] = np.nan; T['lobe_radial'] = np.nan; T['lobe_burial'] = np.nan; T['lobe_id'] = np.nan
    sections = []
    for side in ['L', 'R']:
        m = T.side == side
        ped_pts, cyl = peduncle_cylinder({i: paths[i] for i in T.loc[m, 'nid']}, stretch)
        ids = np.array(list(ped_pts)); T.loc[m, 'in_cylinder'] = T.loc[m, 'nid'].isin(ped_pts)
        pts = np.vstack([ped_pts[i] for i in ids]); lab = np.concatenate([np.full(len(ped_pts[i]), i) for i in ids])
        bm, ax = bundle_measures(pts, lab, slab=4.0, min_cov=0.8, min_pts=50, axis=cyl['axis'])
        T.loc[m, 'ped_radial'] = T.loc[m, 'nid'].map(bm.radial); T.loc[m, 'ped_burial'] = T.loc[m, 'nid'].map(bm.burial)
        axes.append(dict(ds=ds, bundle=f'KC peduncle {side}', cylinder=cyl, **ax))
        # one cross-section near the middle of the peduncle for the rendering
        c = np.array(ax['centroid']); a = np.array(ax['axis']); e1 = np.array(ax['e1']); e2 = np.cross(a, e1)
        t = (pts - c) @ a; mid = np.median(t); sel = np.abs(t - mid) < 2.0
        sec = pd.DataFrame({'nid': lab[sel], 'u': (pts[sel] - c) @ e1, 'v': (pts[sel] - c) @ e2}).groupby('nid').mean().reset_index()
        sec['side'] = side; sec['ds'] = ds; sec['t_um'] = mid; sections.append(sec)
        # alpha/beta lobes: assign each alpha/beta fibre to the vertical or medial lobe by k-means of segment centroids
        ab = T[m & (T['rank'] >= 3) & T.nid.isin(lobe_pts.keys())]
        if ds != 'BANC' and len(ab) > 100:
            C = np.vstack([lobe_pts[i].mean(0) for i in ab.nid]); km = KMeans(2, n_init=10, random_state=0).fit(C)
            for L_ in [0, 1]:
                ii = ab.nid.values[km.labels_ == L_]
                p2 = np.vstack([lobe_pts[i] for i in ii]); l2 = np.concatenate([np.full(len(lobe_pts[i]), i) for i in ii])
                bm2, ax2 = bundle_measures(p2, l2, slab=4.0, min_cov=0.5, min_pts=50)
                T.loc[T.nid.isin(ii), 'lobe_radial'] = T.loc[T.nid.isin(ii), 'nid'].map(bm2.radial)
                T.loc[T.nid.isin(ii), 'lobe_burial'] = T.loc[T.nid.isin(ii), 'nid'].map(bm2.burial)
                T.loc[T.nid.isin(ii), 'lobe_id'] = L_
                if ax2: axes.append(dict(ds=ds, bundle=f'KC alpha/beta lobe {side}{L_}', **ax2))
        # soma: tangential distance from the axis running from the peduncle centroid through the clone's soma centroid
        if ds == 'BANC':   # BANC has no clone labels: split each mushroom-body soma cluster into four clones by k-means
            S = T.loc[m, ['x', 'y', 'z']].values; T.loc[m, 'clone'] = 'k' + pd.Series(KMeans(4, n_init=10, random_state=0).fit(S).labels_).astype(str).values
        pc = c
        for cl in T.loc[m, 'clone'].dropna().unique():
            mm = m & (T.clone == cl); S = T.loc[mm, ['x', 'y', 'z']].values; cs = S.mean(0); axv = (cs - pc) / np.linalg.norm(cs - pc)
            dep = (S - cs) @ axv; T.loc[mm, 'soma_tang_clone'] = np.linalg.norm((S - cs) - np.outer(dep, axv), axis=1)
    T['ds'] = ds
    T.to_parquet(f'{RES}/kc_{ds}.parquet'); pd.concat(sections).to_parquet(f'{RES}/kc_section_{ds}.parquet')
    st = {'ds': ds, 'n': int(len(T)), 'n_in_cylinder': int(T.in_cylinder.sum())}
    for col in ['ped_radial', 'ped_burial', 'soma_tang_clone']:
        v = T.dropna(subset=[col]); rho, p, n = spearman(v['rank'], v[col]); st[col] = {'rho': rho, 'p': p, 'n': n}
        st[col]['class_median'] = {LAB[k]: float(g[col].median()) for k, g in v.groupby('rank')}
        st[col]['adjacent_auc'] = {f'{LAB[k]}>{LAB[k + 1]}': auc(v[v['rank'] == k][col], v[v['rank'] == k + 1][col]) for k in range(6)
                                   if (v['rank'] == k).sum() > 5 and (v['rank'] == k + 1).sum() > 5}
        st[col]['per_side_rho'] = {s: spearman(g['rank'], g[col])[0] for s, g in v.groupby('side')}
    v = T[(T['rank'] >= 3)].dropna(subset=['lobe_radial'])
    if len(v) > 50:
        rho, p, n = spearman(v['rank'], v['lobe_radial']); st['lobe_radial'] = {'rho': rho, 'p': p, 'n': n,
            'class_median': {LAB[k]: float(g.lobe_radial.median()) for k, g in v.groupby('rank')},
            'auc_s_vs_c': auc(v[v['rank'] == 4].lobe_radial, v[v['rank'] == 6].lobe_radial),
            'auc_p_vs_c': auc(v[v['rank'] == 3].lobe_radial, v[v['rank'] == 6].lobe_radial),
            'auc_p_vs_s': auc(v[v['rank'] == 3].lobe_radial, v[v['rank'] == 4].lobe_radial)}
    save_json(st, f'kc_{ds}.json'); save_json(axes, f'axes_kc_{ds}.json'); print(json.dumps(st, indent=1, default=float)[:3000])

if __name__ == '__main__':
    run(sys.argv[1])
