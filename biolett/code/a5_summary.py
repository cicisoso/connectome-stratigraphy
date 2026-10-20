"""Summary statistics across connectomes: Kenyon-cell control, clonally dated lineages, hemilineage units and connectivity.
Predictive performance is the concordance C: the probability that, of two neurons of the same frame with different
birth ranks, the later-born one has the larger geometric age index (ties count one half).
Usage: python a5_summary.py
"""
from common import *

DS = ['MaleCNS', 'FlyWire', 'BANC']
TAB = os.path.join(ROOT, 'tables'); os.makedirs(TAB, exist_ok=True)

def concordance(df, rank, score, by):
    num = den = 0.0
    for _, g in df.dropna(subset=[rank, score]).groupby(by):
        r = g[rank].values.astype(float); s = g[score].values.astype(float)
        dr = np.sign(r[:, None] - r[None, :]); ds_ = np.sign(s[:, None] - s[None, :]); m = dr > 0
        num += ((ds_[m] > 0) + 0.5 * (ds_[m] == 0)).sum(); den += m.sum()
    return num / den if den else np.nan, int(den)

def fisher_ci(rho, n):
    if not (rho == rho) or n < 5: return np.nan, np.nan
    z = np.arctanh(rho); se = 1.06 / np.sqrt(n - 3); return float(np.tanh(z - 1.96 * se)), float(np.tanh(z + 1.96 * se))

def pooled(rhos, ns):
    z = np.array([np.arctanh(r) for r in rhos]); w = np.array([(n - 3) / 1.06 ** 2 for n in ns])
    zm = (w * z).sum() / w.sum(); se = 1 / np.sqrt(w.sum()); p = 2 * stats.norm.sf(abs(zm / se))
    return float(np.tanh(zm)), float(np.tanh(zm - 1.96 * se)), float(np.tanh(zm + 1.96 * se)), float(p)

# ----------------------------------------------------------------------------- Kenyon cells
MAJOR = {0: 0, 1: 0, 2: 1, 3: 2, 4: 2, 5: 2, 6: 2}; MLAB = {0: 'γ', 1: "α′β′", 2: 'αβ'}
kc_rows = []
for ds in DS:
    T = pd.read_parquet(f'{RES}/kc_{ds}.parquet'); T['major'] = T['rank'].map(MAJOR)
    for col, sign in [('ped_radial', -1), ('ped_burial', 1), ('soma_tang_clone', -1)]:
        v = T.dropna(subset=[col]); rho, p, n = spearman(v['major'], v[col]); lo, hi = fisher_ci(rho, n)
        C, npair = concordance(v.assign(score=sign * v[col]), 'major', 'score', 'side')
        aucs = {f'{MLAB[k]}>{MLAB[k + 1]}': auc(sign * v[v.major == k + 1][col], sign * v[v.major == k][col]) for k in range(2)}
        r7 = spearman(v['rank'], v[col])[0] if ds != 'BANC' else np.nan
        kc_rows.append(dict(ds=ds, measure=col, n=n, rho_major=rho, ci_lo=lo, ci_hi=hi, p=p, C=C, auc_g_apbp=aucs['γ>α′β′'],
                            auc_apbp_ab=aucs['α′β′>αβ'], rho_7class=r7, per_side={s: round(spearman(g['major'], g[col])[0], 3) for s, g in v.groupby('side')},
                            medians={MLAB[k]: round(float(g[col].median()), 3) for k, g in v.groupby('major')},
                            note='BANC: clones from k-means of soma positions (no clone labels)' if (ds == 'BANC' and col == 'soma_tang_clone') else ''))
KC = pd.DataFrame(kc_rows); KC.to_csv(f'{RES}/summary_kc.csv', index=False)
print(KC[['ds', 'measure', 'n', 'rho_major', 'ci_lo', 'ci_hi', 'C', 'auc_g_apbp', 'auc_apbp_ab', 'rho_7class']].round(3).to_string())

# ----------------------------------------------------------------------------- clonally dated lineages
SIGN = {'soma_dist': -1, 'soma_depth': 1, 'radial': -1, 'geo_age': 1, 'soma_tang': -1, 'burial': 1}
lin_rows = []
for fr in ['hemilineage', 'pn_only']:
    for ds in DS:
        M = pd.read_parquet(f'{RES}/lineages_{fr}_{ds}.parquet'); M = M[M['rank'].notna()]
        M['frame_id'] = M.hl + M.side
        subsets = [(l, M[M.lineage == l]) for l in ['adPN', 'lPN', 'vPN', 'PAM']] + [('PAM_single_window', M[(M.lineage == 'PAM') & ~M['mixed'].astype(bool)])]
        for lin, Q in subsets:
            for c, sg in SIGN.items():
                v = Q.dropna(subset=[c]); rho, p, n = spearman(v['rank'], v[c]); lo, hi = fisher_ci(rho, n)
                tm = v.groupby(['rank', 'type'])[c].median().reset_index(); rt, pt, nt = spearman(tm['rank'], tm[c])
                C, npair = concordance(v.assign(score=sg * v[c]), 'rank', 'score', 'frame_id')
                lin_rows.append(dict(frame=fr, ds=ds, lineage=lin, measure=c, predicted_sign=sg, n=n, n_types=nt, rho=rho, ci_lo=lo, ci_hi=hi, p=p,
                                     type_rho=rt, type_p=pt, C=C, n_pairs=npair, consistent=bool(rho == rho and np.sign(rho) == sg)))
L = pd.DataFrame(lin_rows); L.to_csv(f'{RES}/summary_lineages.csv', index=False)
pool = []
for (fr, lin, c), g in L[L.n >= 10].groupby(['frame', 'lineage', 'measure']):
    r, lo, hi, p = pooled(g.rho.values, g.n.values)
    pool.append(dict(frame=fr, lineage=lin, measure=c, k=len(g), k_consistent=int(g.consistent.sum()), pooled_rho=r, ci_lo=lo, ci_hi=hi, p=p,
                     mean_C=float(g.C.mean())))
P = pd.DataFrame(pool); P.to_csv(f'{RES}/summary_lineages_pooled.csv', index=False)
print(P[P.frame == 'hemilineage'].round(3).to_string())

# adPN embryonic versus larval
emb = []
for ds in DS:
    M = pd.read_parquet(f'{RES}/lineages_hemilineage_{ds}.parquet'); Q = M[(M.lineage == 'adPN') & M['rank'].notna()]
    emb.append(dict(ds=ds, n_emb=int((Q.epoch == 'embryonic').sum()), n_larval=int((Q.epoch == 'larval').sum()),
                    auc_geo_age=auc(Q[Q.epoch == 'larval'].geo_age, Q[Q.epoch == 'embryonic'].geo_age),
                    prim_label_in_embryonic=int(((Q.epoch == 'embryonic') & (Q.birth == 'primary')).sum()),
                    prim_label_in_larval=int(((Q.epoch == 'larval') & (Q.birth == 'primary')).sum())))
EMB = pd.DataFrame(emb); EMB.to_csv(f'{RES}/summary_adpn_epochs.csv', index=False); print(EMB.round(3).to_string())

# cross-dataset reproducibility of type medians (geometric age index, hemilineage frame)
rep = []
TM = {}
for ds in DS:
    M = pd.read_parquet(f'{RES}/lineages_hemilineage_{ds}.parquet'); M = M[M['rank'].notna()]
    TM[ds] = M.groupby(['lineage', 'type', 'side']).geo_age.median().groupby(['lineage', 'type']).mean()
for i, a in enumerate(DS):
    for b in DS[i + 1:]:
        j = pd.concat([TM[a].rename('a'), TM[b].rename('b')], axis=1).dropna()
        for lin, g in j.groupby(level=0):
            r, p, n = spearman(g.a, g.b); rep.append(dict(pair=f'{a}-{b}', lineage=lin, n_types=n, rho=r, p=p))
        r, p, n = spearman(j.a, j.b); rep.append(dict(pair=f'{a}-{b}', lineage='all', n_types=n, rho=r, p=p))
REP = pd.DataFrame(rep); REP.to_csv(f'{RES}/summary_reproducibility.csv', index=False); print(REP.round(3).to_string())

# ----------------------------------------------------------------------------- hemilineage units
hl_rows = []
for ds in DS:
    U = pd.read_csv(f'{RES}/hemilineage_units_{ds}.csv')
    for reg, Q in U.groupby('region'):
        for c in [x for x in Q.columns if x.startswith(('auc_', 'rho_'))]:
            v = Q[c].dropna()
            if not len(v): continue
            k, n, p = sign_test(v, 0.0 if c.startswith('rho_') else 0.5)
            hl_rows.append(dict(ds=ds, region=reg, measure=c, n_units=n, n_hemilineages=int(Q.loc[v.index, 'hl'].nunique()), median=float(v.median()),
                                q25=float(v.quantile(.25)), q75=float(v.quantile(.75)), n_consistent=k, sign_test_p=p,
                                n_neurons=int(Q.loc[v.index, 'n'].sum()), n_primary=int(Q.loc[v.index, 'n_prim'].sum())))
H = pd.DataFrame(hl_rows); H.to_csv(f'{RES}/summary_hemilineages.csv', index=False)
print(H[H.measure.isin(['auc_soma_dist', 'auc_soma_depth', 'auc_geo_age_soma', 'auc_radial', 'auc_early_vs_late_geo_age_soma'])].round(3).to_string())

# ----------------------------------------------------------------------------- connectivity
con = []
for ds in DS + ['MANC']:
    f = f'{RES}/connectivity_{ds}.json'
    if not os.path.exists(f): continue
    s = json.load(open(f))
    for reg, r in s.items():
        row = dict(ds=ds, region=reg, n_units=r['n_units'], **{f'top5_{k}': v for k, v in r['top5pct'].items()})
        for c, v in r.items():
            if isinstance(v, dict) and 'median' in v: row[c] = v['median']; row[c + '_k'] = f"{v['n_positive']}/{v['n_units']}"; row[c + '_p'] = v['sign_test_p']
        con.append(row)
CON = pd.DataFrame(con); CON.to_csv(f'{RES}/summary_connectivity.csv', index=False)
print(CON[['ds', 'region', 'n_units', 'log2ratio_n_syn', 'log2ratio_size', 'log2ratio_syn_per_size', 'log2ratio_syn_per_size_k', 'top5_frac_primary_top', 'top5_frac_primary_all']].round(3).to_string())
