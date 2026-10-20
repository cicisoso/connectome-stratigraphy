"""Shared loaders, skeleton access and geometry for the four-connectome stratigraphy analysis.

Datasets (all public; coordinates converted to micrometres):
  MaleCNS v1.0  - male brain and nerve cord (Berg et al. 2026), FIB-SEM, 8 nm voxels
  FlyWire v783  - female brain (Dorkenwald et al. 2024; Schlegel et al. 2024), 4 x 4 x 40 nm voxels
  BANC v888     - female brain and nerve cord (Bates et al. 2026), soma = nucleus position in nm
  MANC v1.0     - male nerve cord (Takemura et al. 2024; Marin et al. 2024), 8 nm voxels

Every loader returns one row per neuron with a soma position and the columns
  ds, nid, x, y, z, side, region, super, type, hl, birth, neuromere
where hl is the hemilineage with the '__prim' suffix stripped and birth is 'primary',
'early_secondary', 'secondary' or missing.
"""
import os, json, numpy as np, pandas as pd
from collections import defaultdict, deque
from scipy import stats
from scipy.spatial import ConvexHull

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.environ.get('FLYBRAIN_DATA', os.path.join(ROOT, 'data'))
RES = os.path.join(ROOT, 'results')
EXT = os.path.join(ROOT, 'external')
os.makedirs(RES, exist_ok=True)

VNC_NEUROMERES = ['T1', 'T2', 'T3'] + [f'A{i}' for i in range(1, 11)]
NOT_A_LINEAGE = {'putative_primary', 'primary', 'TBD', 'None', 'nan', 'NA', ''}

# ----------------------------------------------------------------------------- loaders
def _hl_birth(label):
    """Split an Ito-Lee hemilineage label into (hemilineage, birth)."""
    s = label.astype('string')
    prim = s.str.endswith('__prim').fillna(False)
    hl = s.str.replace('__prim$', '', regex=True)
    hl = hl.where(~hl.isin(NOT_A_LINEAGE))
    birth = pd.Series(np.where(prim, 'primary', np.where(hl.notna(), 'secondary', None)), index=label.index)
    birth[s.isin(['putative_primary', 'primary']).fillna(False)] = 'primary'
    return hl, birth

def _super_malecns(sc):
    m = {'cb_intrinsic': 'intrinsic', 'vnc_intrinsic': 'intrinsic', 'ol_intrinsic': 'optic',
         'descending_neuron': 'descending', 'ascending_neuron': 'ascending', 'vnc_motor': 'motor',
         'cb_motor': 'motor', 'visual_projection': 'visual_projection', 'visual_centrifugal': 'visual_centrifugal',
         'cb_endocrine': 'endocrine', 'vnc_endocrine': 'endocrine', 'vnc_efferent': 'efferent', 'cb_efferent': 'efferent',
         'efferent_ascending': 'efferent', 'efferent_descending': 'efferent'}
    out = sc.map(m)
    out[sc.astype(str).str.contains('sensory')] = 'sensory'
    return out.fillna('other')

def load_malecns():
    a = pd.read_feather(f'{DATA}/body-annotations-male-cns-v1.0-minconf-0.5.feather')
    a = a[a.somaLocation.notna()].copy()
    S = np.vstack(a.somaLocation.values).astype(float) * 0.008
    vnc = a.somaNeuromere.isin(VNC_NEUROMERES) | a.superclass.astype(str).str.startswith('vnc') | (a.superclass == 'ascending_neuron')
    hl_b, birth_b = _hl_birth(a.itoleeHl)
    # nerve cord: Truman hemilineage; birth time from the matched MANC neuron (primary / early secondary / secondary)
    m = pd.read_feather(f'{DATA}/manc/manc-v1.0-neuron-properties.feather', columns=['bodyId', 'birthtime']).set_index('bodyId').birthtime
    mb = a.mancBodyid.map(lambda v: m.get(int(v)) if pd.notna(v) else None)
    hl_v = a.trumanHl.astype('string').where(~a.trumanHl.astype('string').isin(NOT_A_LINEAGE))
    birth_v = mb.map({'primary': 'primary', 'early secondary': 'early_secondary', 'secondary': 'secondary'})
    df = pd.DataFrame({'ds': 'MaleCNS', 'nid': a.bodyId.astype('int64').values, 'x': S[:, 0], 'y': S[:, 1], 'z': S[:, 2],
                       'side': a.somaSide.where(a.somaSide.isin(['L', 'R'])).values,
                       'region': np.where(vnc, 'vnc', 'brain'), 'super': _super_malecns(a.superclass).values,
                       'cls': a['class'].values, 'type': a.type.astype('string').values,
                       'hl': np.where(vnc, hl_v, hl_b), 'birth': np.where(vnc, birth_v, birth_b),
                       'neuromere': a.somaNeuromere.where(a.somaNeuromere.isin(VNC_NEUROMERES)).values,
                       'malecns_early': (a.birthtime == 'early').values})
    df['hl'] = df.hl.astype(object).where(pd.notna(df.hl), None)
    df['birth'] = df.birth.astype(object).where(pd.notna(df.birth), None)
    return df.reset_index(drop=True)

def load_flywire():
    f = pd.read_csv(f'{DATA}/flywire/annotations/Supplemental_file1_neuron_annotations.tsv', sep='\t', low_memory=False)
    f = f[f.soma_x.notna()].copy()
    hl, birth = _hl_birth(f.ito_lee_hemilineage)
    smap = {'central': 'intrinsic', 'optic': 'optic', 'descending': 'descending', 'ascending': 'ascending',
            'motor': 'motor', 'visual_projection': 'visual_projection', 'visual_centrifugal': 'visual_centrifugal',
            'endocrine': 'endocrine', 'sensory': 'sensory'}
    df = pd.DataFrame({'ds': 'FlyWire', 'nid': f.root_id.astype('int64').values,
                       'x': f.soma_x.values * 0.004, 'y': f.soma_y.values * 0.004, 'z': f.soma_z.values * 0.040,
                       'side': f.side.map({'left': 'L', 'right': 'R'}).values, 'region': 'brain',
                       'super': f.super_class.map(smap).fillna('other').values, 'cls': f.cell_class.values,
                       'type': f.cell_type.astype('string').values, 'hemibrain_type': f.hemibrain_type.astype('string').values,
                       'hl': hl.values, 'birth': birth.values, 'neuromere': None})
    df['hl'] = df.hl.astype(object).where(pd.notna(df.hl), None)
    return df.reset_index(drop=True)

BANC_SWC = f'{DATA}/banc/skel/swcs-from-pcg-skel'
def load_banc():
    b = pd.read_feather(f'{DATA}/banc/banc_888_meta.feather')
    b = b[b.nucleus_id.notna() & b.root_position_nm.notna()].copy()
    P = np.array([[float(v) for v in str(s).split(',')] for s in b.root_position_nm]) / 1000.0
    files = set(int(f[:-4]) for f in os.listdir(BANC_SWC) if f.endswith('.swc'))
    sk = pd.Series(pd.NA, index=b.index, dtype='Int64')
    for c in ['root_888', 'root_890', 'root_850', 'root_626', 'root_id']:   # skeleton files are named by the root id of one release
        v = pd.to_numeric(b[c], errors='coerce').astype('Int64')
        hit = v.isin(files) & sk.isna()
        sk[hit] = v[hit]
    hl, birth = _hl_birth(b.hemilineage)
    vnc = b.region == 'ventral_nerve_cord'
    smap = {'central_brain_intrinsic': 'intrinsic', 'ventral_nerve_cord_intrinsic': 'intrinsic', 'optic_lobe_intrinsic': 'optic',
            'descending': 'descending', 'ascending': 'ascending', 'motor': 'motor', 'visual_projection': 'visual_projection',
            'visual_centrifugal': 'visual_centrifugal', 'sensory': 'sensory', 'sensory_ascending': 'sensory'}
    nm = b.neuromere.astype('string').str.upper().str.replace('PROTHORACIC', 'T1').str.replace('MESOTHORACIC', 'T2').str.replace('METATHORACIC', 'T3')
    df = pd.DataFrame({'ds': 'BANC', 'nid': pd.to_numeric(b.root_id, errors='coerce').astype('int64').values, 'skel': sk.values,
                       'x': P[:, 0], 'y': P[:, 1], 'z': P[:, 2], 'side': b.side.map({'left': 'L', 'right': 'R'}).values,
                       'region': np.where(vnc, 'vnc', 'brain'), 'super': b.super_class.map(smap).fillna('other').values,
                       'cls': b.cell_class.values, 'type': b.cell_type.astype('string').values,
                       'hl': hl.values, 'birth': np.where(vnc, None, birth.values), 'neuromere': nm.values,
                       'in_syn': b.input_connections.values, 'out_syn': b.output_connections.values,
                       'cable_um': b.l2_cable_length_um.values})
    df['hl'] = df.hl.astype(object).where(pd.notna(df.hl), None)
    return df.drop_duplicates('nid').reset_index(drop=True)   # a few root ids carry two nucleus records

def load_manc():
    m = pd.read_feather(f'{DATA}/manc/manc-v1.0-neuron-properties.feather')
    m = m[m.somaLocation.notna()].copy()
    S = np.vstack(m.somaLocation.values).astype(float) * 0.008
    smap = {'intrinsic neuron': 'intrinsic', 'ascending neuron': 'ascending', 'descending neuron': 'descending',
            'motor neuron': 'motor', 'sensory neuron': 'sensory', 'sensory ascending': 'sensory', 'efferent neuron': 'efferent',
            'efferent ascending': 'efferent'}
    hl = m.hemilineage.astype('string').where(~m.hemilineage.astype('string').isin(NOT_A_LINEAGE))
    df = pd.DataFrame({'ds': 'MANC', 'nid': m.bodyId.astype('int64').values, 'x': S[:, 0], 'y': S[:, 1], 'z': S[:, 2],
                       'side': m.somaSide.map({'LHS': 'L', 'RHS': 'R'}).values, 'region': 'vnc',
                       'super': m['class'].map(smap).fillna('other').values, 'cls': m['class'].values,
                       'type': m.type.astype('string').values, 'hl': hl.values,
                       'birth': m.birthtime.map({'primary': 'primary', 'early secondary': 'early_secondary', 'secondary': 'secondary'}).values,
                       'neuromere': m.somaNeuromere.where(m.somaNeuromere.isin(VNC_NEUROMERES)).values,
                       'pre': m.pre.values, 'post': m.post.values, 'size_vox': m['size'].values})
    df['hl'] = df.hl.astype(object).where(pd.notna(df.hl), None)
    return df.reset_index(drop=True)

LOADERS = {'MaleCNS': load_malecns, 'FlyWire': load_flywire, 'BANC': load_banc, 'MANC': load_manc}

# ----------------------------------------------------------------------------- skeletons
class Skeleton:
    """Nodes (N,3) in µm and an undirected adjacency list."""
    __slots__ = ('xyz', 'adj')
    def __init__(self, xyz, parent_idx):
        self.xyz = np.asarray(xyz, float)
        adj = [[] for _ in range(len(self.xyz))]
        for k, p in enumerate(parent_idx):
            if p >= 0: adj[k].append(p); adj[p].append(k)
        self.adj = adj

def _from_swc(path, scale):
    d = pd.read_csv(path, sep=r'\s+', comment='#', header=None, names=['id', 't', 'x', 'y', 'z', 'r', 'p'])
    idx = pd.Series(np.arange(len(d)), index=d.id.values)
    par = d.p.map(idx).fillna(-1).astype(int).values
    return Skeleton(d[['x', 'y', 'z']].values * scale, par)

_FW = {}
def flywire_cache_path(): return f'{RES}/flywire_skeleton_cache.parquet'
def build_flywire_cache(root_ids):
    """Extract the skeletons of the requested FlyWire neurons from the 5.4 GB Zenodo parquet once."""
    import pyarrow.parquet as pq, pyarrow as pa, pyarrow.compute as pc
    src = f'{DATA}/flywire/skeletons/sk_lod1_783_healed_ds2.parquet'
    pf = pq.ParquetFile(src); want = pa.array(sorted(set(int(r) for r in root_ids)), type=pa.int64())
    parts = []
    for g in range(pf.num_row_groups):
        t = pf.read_row_group(g)
        col = 'neuron'
        t = t.filter(pc.is_in(pc.cast(t[col], pa.int64()), value_set=want))
        if t.num_rows: parts.append(t)
    out = pa.concat_tables(parts); pq.write_table(out, flywire_cache_path()); return out.num_rows

def load_flywire_cache():
    if not _FW:
        t = pd.read_parquet(flywire_cache_path())
        col = 'neuron'
        for rid, g in t.groupby(col, sort=False): _FW[int(rid)] = g
    return _FW

def skeleton(ds, nid, skel_id=None):
    try:
        if ds == 'MaleCNS':
            p = f'{DATA}/swc/{int(nid)}.swc'
            return _from_swc(p, 0.008) if os.path.exists(p) else None
        if ds == 'BANC':
            if skel_id is None or pd.isna(skel_id): return None
            p = f'{BANC_SWC}/{int(skel_id)}.swc'
            return _from_swc(p, 1.0) if os.path.exists(p) else None
        if ds == 'FlyWire':
            g = load_flywire_cache().get(int(nid))
            if g is None: return None
            idx = pd.Series(np.arange(len(g)), index=g.node_id.values)
            par = g.parent_id.map(idx).fillna(-1).astype(int).values
            return Skeleton(g[['x', 'y', 'z']].values / 1000.0, par)
    except Exception:
        return None
    return None

# ----------------------------------------------------------------------------- geometry
def geodesic(sk, s):
    n = len(sk.xyz); dist = np.full(n, np.inf); prev = np.full(n, -1); dist[s] = 0; q = deque([s])
    while q:
        u = q.popleft()
        for v in sk.adj[u]:
            if dist[v] == np.inf:
                dist[v] = dist[u] + np.linalg.norm(sk.xyz[u] - sk.xyz[v]); prev[v] = u; q.append(v)
    return dist, prev

def soma_path(sk, soma, max_gap=12.0):
    """Path from the skeleton node nearest the soma to the most distant tip; None if the soma is not on the skeleton."""
    d0 = np.linalg.norm(sk.xyz - soma, axis=1); s = int(np.argmin(d0))
    if d0[s] > max_gap: return None
    dist, prev = geodesic(sk, s)
    tip = int(np.argmax(np.where(np.isfinite(dist), dist, -1))); path = []; u = tip
    while u != -1: path.append(u); u = prev[u]
    path = np.array(path[::-1]); return path, dist[path]

def neurite_segment(sk, soma, lo, hi, step=1.0):
    """Points of the cell body fibre between geodesic distances lo and hi from the soma, resampled every `step` µm."""
    r = soma_path(sk, soma)
    if r is None: return None
    path, gd = r
    if gd[-1] < hi: return None
    P = sk.xyz[path]; t = np.arange(lo, hi + 1e-9, step)
    return np.c_[np.interp(t, gd, P[:, 0]), np.interp(t, gd, P[:, 1]), np.interp(t, gd, P[:, 2])]

def geometric_median(P, it=100):
    gm = P.mean(0)
    for _ in range(it):
        w = 1.0 / np.maximum(np.linalg.norm(P - gm, axis=1), 1e-6); g2 = (P * w[:, None]).sum(0) / w.sum()
        if np.linalg.norm(g2 - gm) < 1e-5: break
        gm = g2
    return gm

def bundle_measures(points, ids, slab=3.0, min_cov=0.6, min_pts=20, axis=None):
    """Bundle of fibre points from many neurons.
    Fits the axis (first principal component), cuts slabs perpendicular to it and, within each slab, gives each fibre's
    distance from the geometric median of the cross-section (radial) and from the convex-hull boundary (burial depth),
    both divided by the slab's 95th-percentile radius. Per-neuron values are medians over slabs; only slabs that contain
    at least `min_cov` of the neurons are used. Returns (per-neuron table, axis record)."""
    ids = np.asarray(ids); X = np.asarray(points, float); c = X.mean(0)
    ax = np.linalg.svd(X - c, full_matrices=False)[2][0] if axis is None else np.asarray(axis, float)
    t = (X - c) @ ax; perp = (X - c) - np.outer(t, ax)
    e1 = np.linalg.svd(perp, full_matrices=False)[2][0]; e2 = np.cross(ax, e1); u = perp @ e1; v = perp @ e2
    bins = np.floor(t / slab).astype(int); n_total = len(np.unique(ids)); rec = defaultdict(list); used = []
    for b_ in np.unique(bins):
        m = bins == b_
        nid = np.unique(ids[m])
        if len(nid) < min_cov * n_total or m.sum() < min_pts: continue
        P = np.c_[u[m], v[m]]
        try: hull = ConvexHull(P)
        except Exception: continue
        dedge = np.min(-(P @ hull.equations[:, :2].T + hull.equations[:, 2]), axis=1)
        gm = geometric_median(P); rad = np.linalg.norm(P - gm, axis=1); r95 = np.quantile(rad, .95)
        if r95 <= 0: continue
        used.append((b_ * slab, (b_ + 1) * slab, int(len(nid)), float(r95)))
        for i, dd, rr in zip(ids[m], dedge, rad): rec[i].append((rr / r95, dd / r95, rr, dd))
    if not rec: return pd.DataFrame(columns=['radial', 'burial', 'radial_um', 'burial_um']), None
    out = pd.DataFrame({i: np.median(np.array(v), axis=0) for i, v in rec.items()}, index=['radial', 'burial', 'radial_um', 'burial_um']).T
    axrec = {'centroid': c.round(3).tolist(), 'axis': ax.round(5).tolist(), 'e1': e1.round(5).tolist(), 'slab_um': slab,
             'min_coverage': min_cov, 'n_neurons': int(n_total), 'slabs_used': len(used),
             'slab_range_um': [float(min(s[0] for s in used)), float(max(s[1] for s in used))] if used else None}
    return out, axrec

def side_branch_cable(sk, u, p_prev, p_next):
    """Cable (µm) of the side branch(es) leaving path node u other than towards p_prev and p_next."""
    seen = {u, p_prev, p_next}; stack = [v for v in sk.adj[u] if v not in (p_prev, p_next)]
    for v in stack: seen.add(v)
    tot = sum(np.linalg.norm(sk.xyz[u] - sk.xyz[v]) for v in stack)
    while stack:
        w = stack.pop()
        for v in sk.adj[w]:
            if v not in seen: seen.add(v); tot += np.linalg.norm(sk.xyz[w] - sk.xyz[v]); stack.append(v)
    return tot

def cell_body_fibre(sk, soma, skip=5.0, min_branch=4.0, max_len=150.0):
    """Cell body fibre: the soma-to-tip path from the soma to its first side branch carrying more than `min_branch` µm
    of cable (branches within `skip` µm of the soma are ignored), resampled every 1 µm and capped at `max_len` µm."""
    r = soma_path(sk, soma)
    if r is None: return None
    path, gd = r; end = min(gd[-1], max_len)
    for j in range(1, len(path) - 1):
        if gd[j] > end: break
        if gd[j] > skip and len(sk.adj[path[j]]) > 2 and side_branch_cable(sk, path[j], path[j - 1], path[j + 1]) > min_branch:
            end = gd[j]; break
    if end < 10: return None
    P = sk.xyz[path]; t = np.arange(0, end + 1e-9, 1.0)
    return np.c_[np.interp(t, gd, P[:, 0]), np.interp(t, gd, P[:, 1]), np.interp(t, gd, P[:, 2])]

_OFFS = np.array([(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)])
def tract_section(fibres, grid=2.0, half=6.0):
    """Cross-section of a lineage tract where its cell body fibres converge most.
    The convergence point is the centre of the 2-µm voxel (dilated by one voxel) crossed by the most fibres; the tract axis
    is the principal axis of the fibre points within 2*half µm (along each fibre) of their closest approach to that point;
    each fibre contributes its longest contiguous run of points inside the cylinder |t| <= half, distance <= R, where R is
    the 99th percentile of the distance of those near points from the axis. Returns (points per neuron, record)."""
    ids = [i for i, P in fibres.items() if P is not None and len(P) >= 5]
    if len(ids) < 10: return {}, None
    rows = []
    for k, i in enumerate(ids):
        V = np.unique(np.floor(fibres[i] / grid).astype(np.int64), axis=0)
        V = np.unique((V[:, None, :] + _OFFS[None]).reshape(-1, 3), axis=0); rows.append(np.c_[V, np.full(len(V), k)])
    A = pd.DataFrame(np.vstack(rows), columns=['x', 'y', 'z', 'k']); cnt = A.groupby(['x', 'y', 'z']).size()
    top = cnt[cnt == cnt.max()].index.to_frame().values.astype(float); c0 = (np.median(top, axis=0) + 0.5) * grid
    near = []
    for i in ids:
        P = fibres[i]; dd = np.linalg.norm(P - c0, axis=1); j = int(np.argmin(dd))
        if dd[j] <= 1.5 * grid: near.append(P[max(0, j - int(2 * half)):j + int(2 * half) + 1])
    X = np.vstack(near); a = np.linalg.svd(X - X.mean(0), full_matrices=False)[2][0]
    R = float(np.quantile(np.linalg.norm((X - c0) - np.outer((X - c0) @ a, a), axis=1), 0.99))
    pts = {}
    for i in ids:
        P = fibres[i]; t = (P - c0) @ a; inside = (np.abs(t) <= half) & (np.linalg.norm((P - c0) - np.outer(t, a), axis=1) <= R)
        if not inside.any(): continue
        runs = np.flatnonzero(np.diff(np.r_[0, inside.astype(int), 0])).reshape(-1, 2); best = runs[np.argmax(runs[:, 1] - runs[:, 0])]
        if best[1] - best[0] >= half: pts[i] = P[best[0]:best[1]]
    return pts, dict(convergence_point=c0.round(3).tolist(), tract_axis=a.round(5).tolist(), radius_um=round(R, 3), half_length_um=half,
                     max_fibres_per_voxel=int(cnt.max()), n_fibres=len(ids), n_near=len(near), n_in_section=len(pts))

def fibre_path(sk, soma, dmax=60.0):
    """Soma-to-tip path resampled every 1 µm from the soma out to `dmax` µm (None if the soma is not on the skeleton or the
    path is shorter than 31 µm)."""
    r = soma_path(sk, soma)
    if r is None: return None
    path, gd = r
    if gd[-1] < 31: return None
    P = sk.xyz[path]; t = np.arange(0, min(gd[-1], dmax) + 1e-9, 1.0)
    return np.c_[np.interp(t, gd, P[:, 0]), np.interp(t, gd, P[:, 1]), np.interp(t, gd, P[:, 2])]

def frame_geometry(F, paths, label='', slab=2.0, min_cov=0.6, grid=3.0):
    """Soma and fibre measures for the neurons of one frame (one hemilineage on one side), as defined in the original analysis:
    soma_dist   distance of the soma from the frame's soma centroid;
    soma_depth  position along the axis from the centroid of the cell body fibres (8-30 µm from each soma) to the soma
                centroid, i.e. the direction in which the fibres leave the cluster (positive = away from the neuropil);
    soma_tang   distance from that axis;
    radial / burial  position of the fibre in the cross-section of the tract where the fibres (soma paths 0-60 µm)
                converge most (tract_section, bundle_measures).
    Returns (per-neuron table, frame record)."""
    F = F.copy(); S = F[['x', 'y', 'z']].values.astype(float); cs = S.mean(0); F['soma_dist'] = np.linalg.norm(S - cs, axis=1)
    seg = [paths[i][8:31] for i in F.nid if paths.get(i) is not None and len(paths[i]) >= 31]
    if len(seg) < 5: return F, None
    oc = np.vstack(seg).mean(0); axv = cs - oc; axv /= np.linalg.norm(axv); dep = (S - cs) @ axv
    F['soma_depth'] = dep; F['soma_tang'] = np.linalg.norm((S - cs) - np.outer(dep, axv), axis=1)
    rec = dict(frame=label, n_somata=int(len(F)), n_paths=len(seg), soma_centroid=cs.round(3).tolist(), outflow_centroid=oc.round(3).tolist(),
               soma_axis=axv.round(5).tolist(), outflow_to_soma_um=round(float(np.linalg.norm(cs - oc)), 2))
    pts, sec = tract_section({i: paths.get(i) for i in F.nid}, grid=grid)
    F['in_section'] = F.nid.isin(pts)
    if sec is not None:
        rec['tract'] = sec
        if len(pts) >= 10:
            bm, ax = bundle_measures(np.vstack(list(pts.values())), np.concatenate([np.full(len(v), i) for i, v in pts.items()]),
                                     slab=slab, min_cov=min_cov, min_pts=10, axis=sec['tract_axis'])
            F['radial'] = F.nid.map(bm.radial); F['burial'] = F.nid.map(bm.burial)
            if ax: rec['section'] = ax
    return F, rec

def add_geo_age(F, soma_only=False):
    """Geometric age index: mean of the within-frame z-scores of soma_depth, -soma_dist and (unless soma_only) -radial."""
    z = []
    for c, sgn in [('soma_depth', 1), ('soma_dist', -1)] + ([] if soma_only else [('radial', -1)]):
        if c in F and F[c].notna().sum() > 2: z.append(sgn * (F[c] - F[c].mean()) / F[c].std())
    F['geo_age_soma' if soma_only else 'geo_age'] = pd.concat(z, axis=1).mean(1, skipna=True) if z else np.nan
    return F

def soma_measures(S, bundle_centroid=None):
    """Distance of each soma from the cluster centroid and, if a neurite-bundle centroid is given, position along the
    axis from the bundle centroid to the soma centroid (positive = away from the neuropil) and distance from that axis."""
    S = np.asarray(S, float); cs = S.mean(0); out = {'soma_dist': np.linalg.norm(S - cs, axis=1)}
    if bundle_centroid is not None:
        axv = cs - np.asarray(bundle_centroid, float); axv /= np.linalg.norm(axv)
        dep = (S - cs) @ axv; out['soma_depth'] = dep; out['soma_tang'] = np.linalg.norm((S - cs) - np.outer(dep, axv), axis=1)
        out['axis'] = axv
    out['centroid'] = cs
    return out

# ----------------------------------------------------------------------------- statistics
def auc(x1, x2):
    x1 = np.asarray(x1, float); x2 = np.asarray(x2, float); x1 = x1[np.isfinite(x1)]; x2 = x2[np.isfinite(x2)]
    if len(x1) < 2 or len(x2) < 2: return np.nan
    return stats.mannwhitneyu(x1, x2).statistic / (len(x1) * len(x2))

def spearman(x, y):
    x = np.asarray(x, float); y = np.asarray(y, float); m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 5: return np.nan, np.nan, int(m.sum())
    r = stats.spearmanr(x[m], y[m]); return float(r[0]), float(r[1]), int(m.sum())

def sign_test(v, null=0.5):
    v = np.asarray(v, float); v = v[np.isfinite(v)]; k = int((v > null).sum()); n = int((v != null).sum())
    return k, n, (stats.binomtest(k, n).pvalue if n else np.nan)

def zscore_within(df, cols, by):
    out = df.copy()
    for c in cols:
        g = out.groupby(by)[c]; out[c + '_z'] = (out[c] - g.transform('mean')) / g.transform('std')
    return out

def save_json(obj, name):
    def conv(o):
        if isinstance(o, (np.integer,)): return int(o)
        if isinstance(o, (np.floating,)): return None if not np.isfinite(o) else float(o)
        if isinstance(o, np.ndarray): return o.tolist()
        raise TypeError(type(o))
    json.dump(obj, open(f'{RES}/{name}', 'w'), indent=1, default=conv)
