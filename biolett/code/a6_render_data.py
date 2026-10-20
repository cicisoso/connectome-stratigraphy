"""Geometry for the renderings: soma clusters, cell body fibres, tract cross-sections, axes and centroids of every frame of
the clonally dated lineages, and of example hemilineage units.
Usage: python a6_render_data.py MaleCNS|FlyWire|BANC
"""
import sys
from common import *
import a2_lineages as A2

EXAMPLE_UNITS = {'MaleCNS': ['13A|L|T1', '19A|L|T2'], 'FlyWire': [], 'BANC': ['13A|L|T1']}

def section_coords(pts, rec):
    c0 = np.array(rec['convergence_point']); a = np.array(rec['tract_axis'])
    X = np.vstack(list(pts.values())); perp = (X - c0) - np.outer((X - c0) @ a, a)
    e1 = np.linalg.svd(perp, full_matrices=False)[2][0]; e2 = np.cross(a, e1); out = {}
    for i, P in pts.items():
        t = (P - c0) @ a; m = np.abs(t) <= 1.5
        if not m.any(): m = np.abs(t) == np.abs(t).min()
        q = (P[m] - c0); out[i] = (float((q @ e1).mean()), float((q @ e2).mean()))
    return out, e1, e2

def frame_render(ds, F, fib, label, rank=None):
    out, rec = frame_geometry(F, fib, label=label)
    if rec is None or 'tract' not in rec: return None, None, None
    pts, _ = tract_section({i: fib.get(i) for i in F.nid}, grid=3.0)
    uv, e1, e2 = section_coords(pts, rec['tract']) if pts else ({}, None, None)
    out['u'] = out.nid.map(lambda i: uv.get(i, (np.nan, np.nan))[0]); out['v'] = out.nid.map(lambda i: uv.get(i, (np.nan, np.nan))[1])
    if rank is not None: out['rank'] = out.nid.map(rank)
    rec['section_e1'] = None if e1 is None else e1.round(5).tolist(); rec['section_e2'] = None if e2 is None else e2.round(5).tolist()
    paths = [pd.DataFrame({'nid': i, 'k': np.arange(0, len(P), 2), 'x': P[::2, 0], 'y': P[::2, 1], 'z': P[::2, 2]}) for i, P in fib.items() if P is not None and i in set(F.nid)]
    return out, rec, pd.concat(paths) if paths else None

def run(ds):
    d = LOADERS[ds](); R = A2.lineage_frame(d, ds); rank = R.set_index('nid')['rank']; rank = rank[~rank.index.duplicated()]
    rows, recs, paths = [], [], []
    frames = [(hl, side, d[d.hl.eq(hl) & d.side.eq(side)]) for hls in A2.HL_OF.values() for hl in hls for side in ['L', 'R']]
    t = d.type.astype(str)   # lateral PNs leave ALl1 ventral in their own bundle, apart from the local neurons that dominate the hemilineage
    frames += [('ALl1_ventral', side + '|PN', d[d.hl.eq('ALl1_ventral') & d.side.eq(side) & t.str.endswith('_lPN')]) for side in ['L', 'R']]
    if ds in EXAMPLE_UNITS:
        for u in EXAMPLE_UNITS[ds]:
            hl, side, nm = u.split('|'); frames.append((hl, side + '|' + nm, d[d.hl.eq(hl) & d.side.eq(side) & d.neuromere.eq(nm)]))
    for hl, side, F in frames:
        if len(F) < 10: continue
        fib = {}
        for r in F.itertuples():
            sk = skeleton(ds, r.nid, getattr(r, 'skel', None)); fib[r.nid] = None if sk is None else fibre_path(sk, np.array([r.x, r.y, r.z]))
        out, rec, P = frame_render(ds, F.copy(), fib, f'{hl} {side}', rank)
        if out is None: continue
        out['ds'] = ds; out['frame_hl'] = hl; out['frame_side'] = side; rows.append(out)
        recs.append(dict(ds=ds, hl=hl, side=side, **rec))
        if P is not None: P['ds'] = ds; P['frame_hl'] = hl; P['frame_side'] = side; paths.append(P)
        print(ds, hl, side, len(F), flush=True)
    keep = ['ds', 'frame_hl', 'frame_side', 'nid', 'type', 'birth', 'rank', 'x', 'y', 'z', 'in_section', 'u', 'v', 'soma_dist', 'soma_depth', 'radial']
    pd.concat(rows)[keep].to_parquet(f'{RES}/render_somata_{ds}.parquet'); pd.concat(paths).to_parquet(f'{RES}/render_paths_{ds}.parquet')
    save_json(recs, f'render_frames_{ds}.json')

if __name__ == '__main__':
    run(sys.argv[1])
