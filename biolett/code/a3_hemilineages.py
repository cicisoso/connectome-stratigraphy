"""Primary (embryonic-born) versus secondary (larval-born) neurons of the same hemilineage.

Birth labels are connectome annotations: brain '__prim' labels (FlyWire, MaleCNS, BANC; Schlegel et al. 2024) and the
nerve-cord birthtimes of MANC (Marin et al. 2024), transferred to MaleCNS (mancBodyid) and BANC (manc_match).
Annotators used soma size, primary-neurite calibre and the course of the primary neurite relative to the hemilineage
tract (MANC: periphery or outside of the tract; FlyWire: small bundles next to the tract), so the fibre measures are not
independent of the labels. Soma position (depth, tangential distance) was not an annotation criterion.
A unit is one hemilineage on one side (and, in the nerve cord, in one neuromere); all of its neurons form the frame.
Usage: python a3_hemilineages.py MaleCNS|FlyWire|BANC
"""
import sys
from common import *

def banc_vnc_birth(d):
    b = pd.read_feather(f'{DATA}/banc/banc_888_meta.feather', columns=['root_id', 'manc_match'])
    m = pd.read_feather(f'{DATA}/manc/manc-v1.0-neuron-properties.feather', columns=['bodyId', 'birthtime']).set_index('bodyId').birthtime
    bt = pd.to_numeric(b.manc_match, errors='coerce').map(m).map({'primary': 'primary', 'early secondary': 'early_secondary', 'secondary': 'secondary'})
    s = pd.Series(bt.values, index=pd.to_numeric(b.root_id, errors='coerce').astype('int64').values); s = s[~s.index.duplicated()]
    vnc = d.region == 'vnc'; d.loc[vnc, 'birth'] = d.loc[vnc, 'nid'].map(s).values
    return d

LEVEL = {'primary': 0, 'early_secondary': 1, 'secondary': 2}
# P(a primary neuron looks older than a secondary neuron of the same unit): older = soma nearer the neuropil (smaller depth),
# farther from the cluster axis and centroid, fibre farther from the tract core, smaller geometric age
OLDER_IF_LARGER = {'soma_depth': False, 'soma_tang': True, 'soma_dist': True, 'radial': True, 'burial': False,
                   'geo_age_soma': False, 'geo_age': False}

def run(ds):
    d = LOADERS[ds]()
    if ds == 'BANC': d = banc_vnc_birth(d)
    d = d[d.hl.notna() & d.side.isin(['L', 'R'])].copy()
    d['unit'] = d.hl + '|' + d.side + '|' + np.where(d.region == 'vnc', d.neuromere.astype(str), 'brain')
    d = d[~((d.region == 'vnc') & d.neuromere.isna())]
    g = d.groupby('unit').birth.agg(prim=lambda s: (s == 'primary').sum(), sec=lambda s: (s == 'secondary').sum())
    keep = g[(g.prim >= 5) & (g.sec >= 10)].index; print(ds, 'units', len(keep), flush=True)
    rows, recs, per = [], [], []
    for k, u in enumerate(keep):
        F = d[d.unit == u].copy(); fib = {}
        for r in F.itertuples():
            sk = skeleton(ds, r.nid, getattr(r, 'skel', None))
            fib[r.nid] = None if sk is None else fibre_path(sk, np.array([r.x, r.y, r.z]))
        F['has_fibre'] = F.nid.map(lambda i: fib.get(i) is not None)
        out, rec = frame_geometry(F, fib, label=u)
        if rec is None: continue
        out = add_geo_age(add_geo_age(out, soma_only=True))
        rows.append(out); recs.append(dict(ds=ds, unit=u, **rec))
        P = out[out.birth == 'primary']; S = out[out.birth == 'secondary']
        res = dict(ds=ds, unit=u, hl=u.split('|')[0], side=u.split('|')[1], region='vnc' if u.split('|')[2] != 'brain' else 'brain',
                   n=len(out), n_prim=len(P), n_sec=len(S), n_early=int((out.birth == 'early_secondary').sum()),
                   frac_prim_in_section=float(P.in_section.mean()), frac_sec_in_section=float(S.in_section.mean()))
        for c, larger in OLDER_IF_LARGER.items():
            if c not in out: continue
            a = auc(P[c], S[c]); res[f'auc_{c}'] = a if larger else (1 - a if a == a else np.nan)
        E = out[out.birth == 'early_secondary']
        if len(E) >= 5:
            a = auc(E['geo_age_soma'], S['geo_age_soma']); res['auc_early_vs_late_geo_age_soma'] = 1 - a if a == a else np.nan
            lv = out.birth.map(LEVEL); res['rho_level_geo_age_soma'] = spearman(lv, out['geo_age_soma'])[0]
        per.append(res)
        if (k + 1) % 20 == 0: print(k + 1, 'units done', flush=True)
    A = pd.concat(rows); A['ds'] = ds
    keepc = ['ds', 'nid', 'unit', 'hl', 'side', 'region', 'neuromere', 'type', 'super', 'birth', 'x', 'y', 'z', 'has_fibre', 'in_section',
             'soma_depth', 'soma_tang', 'soma_dist', 'radial', 'burial', 'geo_age_soma', 'geo_age']
    A[[c for c in keepc if c in A.columns]].to_parquet(f'{RES}/hemilineages_{ds}.parquet')
    U = pd.DataFrame(per); U.to_csv(f'{RES}/hemilineage_units_{ds}.csv', index=False); save_json(recs, f'axes_hemilineages_{ds}.json')
    st = {}
    for reg, Q in U.groupby('region'):
        r = {'n_units': int(len(Q)), 'n_hemilineages': int(Q.hl.nunique()), 'n_neurons': int(Q.n.sum()), 'n_primary': int(Q.n_prim.sum())}
        for c in [x for x in Q.columns if x.startswith('auc_') or x.startswith('rho_')]:
            v = Q[c].dropna(); null = 0.0 if c.startswith('rho_') else 0.5; kk, nn, p = sign_test(v, null)
            r[c] = {'median': float(v.median()) if len(v) else None, 'n_units': int(len(v)), 'n_consistent': kk, 'sign_test_p': p}
        st[reg] = r
    save_json(st, f'hemilineages_{ds}.json'); print(json.dumps(st, indent=1, default=float))

if __name__ == '__main__':
    run(sys.argv[1])
