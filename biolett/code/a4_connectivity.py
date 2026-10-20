"""Connectivity of primary versus secondary neurons of the same hemilineage unit (reviewer 3): synapse numbers, synapses
per unit of cable (or volume), and the classes of their presynaptic and postsynaptic partners.
Usage: python a4_connectivity.py MaleCNS|FlyWire|BANC|MANC
"""
import sys
import pyarrow as pa, pyarrow.compute as pc, pyarrow.dataset as pds
from common import *
import common
from a3_hemilineages import banc_vnc_birth

PCLASS = {'sensory': 'sensory', 'ascending': 'ascending', 'descending': 'descending', 'motor': 'motor/efferent',
          'efferent': 'motor/efferent', 'endocrine': 'motor/efferent', 'intrinsic': 'intrinsic', 'optic': 'visual',
          'visual_projection': 'visual', 'visual_centrifugal': 'visual'}
CLASSES = ['intrinsic', 'sensory', 'ascending', 'descending', 'visual', 'motor/efferent', 'other']

def units(ds):
    d = LOADERS[ds]()
    if ds == 'BANC': d = banc_vnc_birth(d)
    d = d[d.hl.notna() & d.side.isin(['L', 'R']) & d.birth.isin(['primary', 'early_secondary', 'secondary'])].copy()
    d['unit'] = d.hl + '|' + d.side + '|' + np.where(d.region == 'vnc', d.neuromere.astype(str), 'brain')
    d = d[~((d.region == 'vnc') & d.neuromere.isna())]
    g = d.groupby('unit').birth.agg(prim=lambda s: (s == 'primary').sum(), sec=lambda s: (s == 'secondary').sum())
    return d[d.unit.isin(g[(g.prim >= 5) & (g.sec >= 10)].index)].copy()

def edges_and_classes(ds, ids):
    want = pa.array(sorted(set(int(i) for i in ids)), type=pa.int64())
    if ds == 'MaleCNS':
        dset = pds.dataset(f'{DATA}/connectome-weights-male-cns-v1.0-minconf-0.5.feather', format='feather')
        flt = pc.is_in(pc.field('body_pre'), value_set=want) | pc.is_in(pc.field('body_post'), value_set=want)
        E = dset.to_table(filter=flt).to_pandas().rename(columns={'body_pre': 'pre', 'body_post': 'post', 'weight': 'w'})
        a = pd.read_feather(f'{DATA}/body-annotations-male-cns-v1.0-minconf-0.5.feather', columns=['bodyId', 'superclass'])
        cls = pd.Series(common._super_malecns(a.superclass).values, index=a.bodyId.astype('int64').values)
    elif ds == 'FlyWire':
        E = pd.read_csv(f'{DATA}/flywire/connections.csv.gz', usecols=['pre_root_id', 'post_root_id', 'syn_count'])
        E = E.rename(columns={'pre_root_id': 'pre', 'post_root_id': 'post', 'syn_count': 'w'}).groupby(['pre', 'post'], as_index=False).w.sum()
        s = set(int(i) for i in ids); E = E[E.pre.isin(s) | E.post.isin(s)]
        c = pd.read_csv(f'{DATA}/flywire/classification.csv.gz', usecols=['root_id', 'super_class'])
        m = {'central': 'intrinsic', 'optic': 'optic', 'descending': 'descending', 'ascending': 'ascending', 'motor': 'motor',
             'visual_projection': 'visual_projection', 'visual_centrifugal': 'visual_centrifugal', 'endocrine': 'endocrine', 'sensory': 'sensory'}
        cls = pd.Series(c.super_class.map(m).fillna('other').values, index=c.root_id.astype('int64').values)
    elif ds == 'BANC':
        dset = pds.dataset(f'{DATA}/banc/banc_888_edgelist_simple_v3.feather', format='feather')
        t = dset.to_table(columns=['pre', 'post', 'count'])
        t = t.filter(pc.or_(pc.is_in(pc.cast(t['pre'], pa.int64()), value_set=want), pc.is_in(pc.cast(t['post'], pa.int64()), value_set=want)))
        E = t.to_pandas().rename(columns={'count': 'w'}); E['pre'] = E.pre.astype('int64'); E['post'] = E.post.astype('int64')
        meta = pd.read_feather(f'{DATA}/banc/banc_888_meta.feather', columns=['root_id', 'super_class'])
        smap = {'central_brain_intrinsic': 'intrinsic', 'ventral_nerve_cord_intrinsic': 'intrinsic', 'optic_lobe_intrinsic': 'optic',
                'descending': 'descending', 'ascending': 'ascending', 'motor': 'motor', 'visual_projection': 'visual_projection',
                'visual_centrifugal': 'visual_centrifugal', 'sensory': 'sensory', 'sensory_ascending': 'sensory', 'efferent': 'efferent',
                'endocrine': 'endocrine'}
        cls = pd.Series(meta.super_class.map(smap).fillna('other').values, index=pd.to_numeric(meta.root_id, errors='coerce').astype('int64').values)
    else:
        E = pd.read_csv(f'{DATA}/manc/traced-connections.csv').rename(columns={'bodyId_pre': 'pre', 'bodyId_post': 'post', 'weight': 'w'})
        s = set(int(i) for i in ids); E = E[E.pre.isin(s) | E.post.isin(s)]
        allm = pd.read_feather(f'{DATA}/manc/manc-v1.0-neuron-properties.feather', columns=['bodyId', 'class'])
        smap = {'intrinsic neuron': 'intrinsic', 'ascending neuron': 'ascending', 'descending neuron': 'descending', 'motor neuron': 'motor',
                'sensory neuron': 'sensory', 'sensory ascending': 'sensory', 'efferent neuron': 'efferent', 'efferent ascending': 'efferent'}
        cls = pd.Series(allm['class'].map(smap).fillna('other').values, index=allm.bodyId.astype('int64').values)
    cls = cls[~cls.index.duplicated()].map(PCLASS).fillna('other')
    return E, cls

def cable(ds, d):
    if ds == 'BANC': return d.set_index('nid').cable_um
    if ds == 'FlyWire':
        c = pd.read_csv(f'{DATA}/flywire/cell_stats.csv.gz', usecols=['root_id', 'length_nm']); return (c.set_index('root_id').length_nm / 1000.0)
    if ds == 'MANC': return d.set_index('nid').size_vox * 0.008 ** 3   # volume (µm^3), no skeletons
    out = {}
    for i in d.nid:
        p = f'{DATA}/swc/{int(i)}.swc'
        if not os.path.exists(p): continue
        x = pd.read_csv(p, sep=r'\s+', comment='#', header=None, names=['id', 't', 'x', 'y', 'z', 'r', 'p'])
        idx = pd.Series(np.arange(len(x)), index=x.id.values); par = x.p.map(idx)
        ok = par.notna().values; P = x[['x', 'y', 'z']].values * 0.008
        out[i] = float(np.linalg.norm(P[ok] - P[par[ok].astype(int).values], axis=1).sum())
    return pd.Series(out)

def run(ds):
    d = units(ds); ids = d.nid.values; print(ds, 'units', d.unit.nunique(), 'neurons', len(d), flush=True)
    E, cls = edges_and_classes(ds, ids); E = E[E.pre != E.post]
    s = set(int(i) for i in ids)
    out_e = E[E.pre.isin(s)].copy(); out_e['pc'] = out_e.post.map(cls).fillna('other')
    in_e = E[E.post.isin(s)].copy(); in_e['pc'] = in_e.pre.map(cls).fillna('other')
    fin = in_e.pivot_table(index='post', columns='pc', values='w', aggfunc='sum', fill_value=0)
    fout = out_e.pivot_table(index='pre', columns='pc', values='w', aggfunc='sum', fill_value=0)
    d['n_in'] = d.nid.map(fin.sum(1)).fillna(0); d['n_out'] = d.nid.map(fout.sum(1)).fillna(0); d['n_syn'] = d.n_in + d.n_out
    for c in CLASSES:
        d[f'in_{c}'] = d.nid.map(fin[c] if c in fin else pd.Series(dtype=float)).fillna(0) / d.n_in.replace(0, np.nan)
        d[f'out_{c}'] = d.nid.map(fout[c] if c in fout else pd.Series(dtype=float)).fillna(0) / d.n_out.replace(0, np.nan)
    cb = cable(ds, d); d['size'] = d.nid.map(cb); d['syn_per_size'] = d.n_syn / d['size'].replace(0, np.nan)
    d['ds'] = d.get('ds', ds)
    cols = ['ds', 'nid', 'unit', 'hl', 'side', 'region', 'neuromere', 'type', 'super', 'birth', 'n_in', 'n_out', 'n_syn', 'size', 'syn_per_size'] + \
           [f'{k}_{c}' for k in ['in', 'out'] for c in CLASSES]
    d[[c for c in cols if c in d.columns]].to_parquet(f'{RES}/connectivity_{ds}.parquet')
    per = []
    for u, Q in d.groupby('unit'):
        P = Q[Q.birth == 'primary']; S = Q[Q.birth == 'secondary']; r = dict(unit=u, region=Q.region.iloc[0], n_prim=len(P), n_sec=len(S))
        for c in ['n_in', 'n_out', 'n_syn', 'size', 'syn_per_size']:
            r[f'log2ratio_{c}'] = float(np.log2((P[c].median() + 1) / (S[c].median() + 1))) if P[c].notna().any() and S[c].notna().any() else np.nan
            r[f'auc_{c}'] = auc(P[c], S[c])
        for k in ['in', 'out']:
            for c in CLASSES: r[f'diff_{k}_{c}'] = float(P[f'{k}_{c}'].mean() - S[f'{k}_{c}'].mean())
        per.append(r)
    U = pd.DataFrame(per); U['ds'] = ds; U.to_csv(f'{RES}/connectivity_units_{ds}.csv', index=False)
    st = {}
    for reg, Q in U.groupby('region'):
        r = {'n_units': int(len(Q))}
        for c in [x for x in Q.columns if x.startswith(('log2ratio_', 'auc_', 'diff_'))]:
            v = Q[c].dropna(); null = 0.5 if c.startswith('auc_') else 0.0; k, n, p = sign_test(v, null)
            r[c] = {'median': float(v.median()) if len(v) else None, 'n_units': int(len(v)), 'n_positive': k, 'sign_test_p': p}
        D = d[d.region == reg]; thr = D.n_syn.quantile(0.95); top = D[D.n_syn >= thr]
        r['top5pct'] = {'frac_primary_top': float((top.birth == 'primary').mean()), 'frac_primary_all': float((D.birth == 'primary').mean()), 'n_top': int(len(top))}
        st[reg] = r
    save_json(st, f'connectivity_{ds}.json')
    for reg, r in st.items():
        print(ds, reg, r['n_units'], {c: (round(r[c]['median'], 3), f"{r[c]['n_positive']}/{r[c]['n_units']}") for c in ['log2ratio_n_syn', 'log2ratio_size', 'log2ratio_syn_per_size', 'auc_n_syn']}, r['top5pct'])
        print('   input-class differences', {c: round(r[f'diff_in_{c}']['median'], 3) for c in CLASSES})
        print('   output-class differences', {c: round(r[f'diff_out_{c}']['median'], 3) for c in CLASSES})

if __name__ == '__main__':
    run(sys.argv[1])
