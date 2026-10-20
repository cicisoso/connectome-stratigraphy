"""Lineages with birth order mapped by clonal analysis, independent of adult geometry.

adPN lineage (ALad1): twin-spot MARCM sequence of Yu et al. (2010)
lPN lineage (ALl1_ventral): twin-spot MARCM sequence of Lin et al. (2012)
vPN lineage (ALv1): mean clone-induction window of each type in Lee et al. (2020, supplementary file 1)
PAM dopaminergic neurons (CREa1_dorsal, CREa2_medial): mean clone-induction window in Lee et al. (2020)
For every lineage and side, the frame is the whole hemilineage (all of its neurons, ranked or not).
Usage: python a2_lineages.py MaleCNS|FlyWire|BANC
"""
import sys
from common import *

# Yu et al. 2010: positions in the full adPN sequence (polyglomerular types occupy positions 9 and 34-37; VM3 is made twice, at 14 and 16).
# Old glomerulus names: '1' = VM7v, VM7 = VM7d, VA1lm = VA1v, VM6+VP1 = VM6 (Grabe et al. 2015).
ADPN = {'DP1m': 0, 'VL2p': 1, 'DA4l': 2, 'VM4': 3, 'VM6': 4, 'VA6': 5, 'DL6': 6, 'VL2a': 7, 'DC1': 8, 'VA7l': 10, 'VA2': 11,
        'DM4': 12, 'DL5': 13, 'VM3': 15, 'DM3': 15, 'DL4': 17,
        'DL1': 18, 'DA3': 19, 'DC2': 20, 'D': 21, 'VA3': 22, 'DC3': 23, 'VA1d': 24, 'VM7v': 25, 'VM7d': 26, 'VM2': 27,
        'DM6': 28, 'VA1v': 29, 'VC4': 30, 'VM5d': 31, 'VM5v': 32, 'VC3': 33, 'DL2v': 38, 'DL2d': 39}
ADPN_EMBRYONIC_LAST = 17
# Lin et al. 2012: VA4-VC2-VC1-DM1-DM2-VA5-VA7m-DA1-DL3-DA1-VM1-DA2-DM5 (DA1 and DL3 are made recurrently)
LPN = {'VA4': 0, 'VC2': 1, 'VC1': 2, 'DM1': 3, 'DM2': 4, 'VA5': 5, 'VA7m': 6, 'DA1': 8, 'DL3': 8, 'VM1': 10, 'DA2': 11, 'DM5': 12}
# Lee et al. 2020 ALv1 types named by glomerulus -> connectome vPN type
VPN_MAP = {'P02(DA1)': 'DA1', 'P03(VL1)': 'VL1', 'P04(VA1lm)': 'VA1v', 'P05(VA1d)': 'VA1d', 'P06(DM1,4)': 'DM4',
           'P07(DL2d)': 'DL2d', 'P08(VL2p)': 'VL2p', 'P09(DP1l,VL2p)': 'DP1l'}
# Lee et al. 2020 PAM types named by mushroom-body zone -> connectome PAM type (instance names of MaleCNS, e.g. PAM11(a1))
PAM_ZONE = {'a1': 'PAM11', 'b1ped': 'PAM09', 'b1': 'PAM10', 'b2': 'PAM04', 'r5': 'PAM01', "b'2r5": 'PAM15', "b'2a": 'PAM02',
            "b'2m": 'PAM06', 'r4': 'PAM08', "b'2p": 'PAM05', "b'1m": 'PAM14', "b'1ap": 'PAM13', 'r3': 'PAM12',
            'r4<r2': 'PAM07', 'r4<r1r2': 'PAM07', 'b2b\'2': 'PAM03', 'r4/5': None}
LEE = f'{EXT}/elife-53518-supp1-v2.xlsx'

def lee_windows(sheet):
    d = pd.read_excel(LEE, sheet, header=0, index_col=0).apply(pd.to_numeric, errors='coerce').fillna(0)
    idx = np.arange(d.shape[1]); n = d.sum(1)
    return pd.DataFrame({'n_clones': n, 'mean_window': (d * idx).sum(1) / n.replace(0, np.nan), 'sum_w': (d * idx).sum(1)})

def pam_birth():
    """Clone-weighted mean induction window of each connectome PAM type, separately for CREa1A and CREa2A.
    Types made in more than one separate window (PAM01, PAM04, PAM06, PAM08) are flagged as mixed."""
    out = {}
    for lin, sheet in [('CREa1_dorsal', 'CREa1A_type'), ('CREa2_medial', 'CREa2A_type')]:
        w = lee_windows(sheet); w = w[w.index.str.startswith('H')]
        w['zone'] = w.index.str.extract(r'\((.*)\)')[0].values; w['pam'] = w.zone.map(PAM_ZONE)
        g = w.dropna(subset=['pam']).groupby('pam').agg(n_clones=('n_clones', 'sum'), sum_w=('sum_w', 'sum'), n_lee_types=('zone', 'size'))
        g['mean_window'] = g.sum_w / g.n_clones; out[lin] = g
    mixed = {'PAM01', 'PAM04', 'PAM06', 'PAM08'}
    return out, mixed

def vpn_birth():
    w = lee_windows('ALv1_type'); w = w[w.index.isin(VPN_MAP)]; w['glom'] = w.index.map(VPN_MAP)
    return w.set_index('glom').mean_window

def lineage_frame(d, ds):
    """Neurons of the five lineages with their birth rank."""
    t = d.type.astype(str)
    glom = t.str.replace(r'_(adPN|lPN|vPN)$', '', regex=True)
    ad = d[d.hl.eq('ALad1') & t.str.endswith('_adPN')].copy(); ad['rank'] = glom[ad.index].map(ADPN); ad['lineage'] = 'adPN'
    ad['epoch'] = np.where(ad['rank'].isna(), None, np.where(ad['rank'] <= ADPN_EMBRYONIC_LAST, 'embryonic', 'larval'))
    lp = d[d.hl.eq('ALl1_ventral') & t.str.endswith('_lPN')].copy(); lp['rank'] = glom[lp.index].map(LPN); lp['lineage'] = 'lPN'
    vb = vpn_birth(); vp = d[d.hl.eq('ALv1') & t.str.endswith('_vPN')].copy(); vp['rank'] = glom[vp.index].map(vb); vp['lineage'] = 'vPN'
    pb, mixed = pam_birth(); pm = []
    for lin in ['CREa1_dorsal', 'CREa2_medial']:
        q = d[d.hl.eq(lin) & t.str.match(r'^PAM\d\d')].copy(); q['pamtype'] = q.type.astype(str).str[:5]
        q['rank'] = q.pamtype.map(pb[lin].mean_window); q['mixed'] = q.pamtype.isin(mixed); q['lineage'] = 'PAM'; pm.append(q)
    pm = pd.concat(pm)
    R = pd.concat([ad, lp, vp, pm])
    return R

HL_OF = {'adPN': ['ALad1'], 'lPN': ['ALl1_ventral'], 'vPN': ['ALv1'], 'PAM': ['CREa1_dorsal', 'CREa2_medial']}
PN_SUFFIX = {'ALad1': '_adPN', 'ALl1_ventral': '_lPN', 'ALv1': '_vPN', 'CREa1_dorsal': 'PAM', 'CREa2_medial': 'PAM'}
MEASURES = ['soma_dist', 'soma_depth', 'radial', 'geo_age', 'soma_tang', 'burial']

def lineage_stats(M):
    st = {}
    for lineage, Q in M.groupby('lineage'):
        Qs = {'all': Q, 'single_window': Q[~Q['mixed'].astype(bool)]} if lineage == 'PAM' else {'all': Q}
        for sub, Q2 in Qs.items():
            key = lineage if sub == 'all' else f'{lineage}_{sub}'; Q2 = Q2[Q2['rank'].notna()]
            res = {'n_neurons': int(len(Q2)), 'n_types': int(Q2['rank'].nunique()), 'n_in_section': int(Q2.in_section.fillna(False).sum())}
            for c in MEASURES:
                rho, p, n = spearman(Q2['rank'], Q2[c]); res[c] = {'rho': rho, 'p': p, 'n': n}
                tm = Q2.groupby(['rank', 'type'])[c].median().reset_index()
                rt, pt, nt = spearman(tm['rank'], tm[c]); res[c].update({'type_rho': rt, 'type_p': pt, 'n_types': nt})
                res[c]['per_side_rho'] = {sd: spearman(g['rank'], g[c])[0] for sd, g in Q2.groupby('side')}
            if lineage == 'adPN':
                for c in MEASURES:
                    res[c]['auc_embryonic_vs_larval'] = auc(Q2[Q2.epoch == 'embryonic'][c], Q2[Q2.epoch == 'larval'][c])
                    lv = Q2[Q2.epoch == 'larval']; res[c]['rho_within_larval'] = spearman(lv['rank'], lv[c])[0]
                    em = Q2[Q2.epoch == 'embryonic']; res[c]['rho_within_embryonic'] = spearman(em['rank'], em[c])[0]
                res['annotation_check'] = {str(k): v for k, v in pd.crosstab(Q2.epoch, Q2.birth.fillna('NA')).to_dict().items()}
            st[key] = res
    return st

def run(ds):
    d = LOADERS[ds](); R = lineage_frame(d, ds)
    print(ds, 'ranked neurons', R[R['rank'].notna()].groupby('lineage').size().to_dict(), flush=True)
    frames = {'hemilineage': [], 'pn_only': []}; recs = []; fib = {}
    for lineage, hls in HL_OF.items():
        for hl in hls:
            for side in ['L', 'R']:
                F = d[d.hl.eq(hl) & d.side.eq(side)].copy()
                if len(F) < 10: continue
                for r in F.itertuples():
                    sk = skeleton(ds, r.nid, getattr(r, 'skel', None))
                    fib[r.nid] = None if sk is None else fibre_path(sk, np.array([r.x, r.y, r.z]))
                suf = PN_SUFFIX[hl]; t = F.type.astype(str); P = F[t.str.endswith(suf) | t.str.startswith(suf)]
                for fr, G in [('hemilineage', F), ('pn_only', P)]:
                    if len(G) < 10: continue
                    out, rec = frame_geometry(G, fib, label=f'{hl} {side} {fr}')
                    out = add_geo_age(out); frames[fr].append(out.assign(frame=fr))
                    if rec: recs.append(dict(ds=ds, lineage=lineage, hl=hl, side=side, frame_type=fr, **rec))
    keep = ['nid', 'side', 'hl', 'type', 'birth', 'x', 'y', 'z', 'n_frame', 'frame', 'in_section'] + MEASURES
    out = {}
    for fr, parts in frames.items():
        A = pd.concat(parts); A['n_frame'] = A.groupby(['hl', 'side']).nid.transform('size'); A = A[[c for c in keep if c in A.columns]]
        M = R[['nid', 'lineage', 'rank'] + [c for c in ['epoch', 'mixed', 'pamtype'] if c in R.columns]].merge(A, on='nid', how='inner'); M['ds'] = ds
        M.to_parquet(f'{RES}/lineages_{fr}_{ds}.parquet'); A.assign(ds=ds).to_parquet(f'{RES}/lineage_frames_{fr}_{ds}.parquet')
        out[fr] = lineage_stats(M)
    save_json(out, f'lineages_{ds}.json'); save_json(recs, f'axes_lineages_{ds}.json')
    for fr in out:
        print('==', fr)
        for k, v in out[fr].items():
            print(k, v['n_neurons'], 'in section', v['n_in_section'], {c: (None if v[c]['rho'] != v[c]['rho'] else round(v[c]['rho'], 2), None if v[c]['type_rho'] != v[c]['type_rho'] else round(v[c]['type_rho'], 2)) for c in MEASURES})

if __name__ == '__main__':
    run(sys.argv[1])
