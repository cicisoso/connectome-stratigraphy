"""Primary (embryonic-born) versus secondary (larval-born) neurons across hemilineages: soma position and primary-neurite bundle depth."""
import os, pandas as pd, numpy as np, json
from strat import *
R=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
a=pd.read_feather(f'{D}/body-annotations-male-cns-v1.0-minconf-0.5.feather')
sel=pd.read_csv(f'{R}/selected_neurons.csv')
a=a[a.bodyId.isin(sel.bodyId)&a.somaLocation.notna()].copy()
# define hemilineage id and class
a['hl']=np.where(a.trumanHl.notna()&(a.trumanHl!='TBD'),'VNC:'+a.trumanHl.astype(str)+':'+a.somaNeuromere.astype(str),'CB:'+a.itoleeHl.astype(str).str.replace('__prim$','',regex=True))
a['primary']=np.where(a.hl.str.startswith('VNC'),a.birthtime=='early',a.itoleeHl.astype(str).str.endswith('__prim'))
a=a[~a['class'].isin(['Kenyon_Cell','ALPN'])]
g=a.groupby(['hl','somaSide']).agg(n=('bodyId','size'),prim=('primary','sum')); g['sec']=g.n-g.prim; keep=g[(g.prim>=5)&(g.sec>=10)].index
print('hemilineage-side units with >=5 primary and >=10 secondary:',len(keep),' neurons',g.loc[keep].n.sum(),flush=True)
rows=[]; per=[]
for (hl,side) in keep:
    Q=a[(a.hl==hl)&(a.somaSide==side)]; pts=[]; ids=[]; S=[]; okid=[]
    for r in Q.itertuples():
        t=tree(r.bodyId)
        if t is None: continue
        xyz,adj=t; sm=np.array(r.somaLocation,dtype=float)*VOX
        P,path,dist=primary_neurite(xyz,adj,sm,lo=8,hi=30)
        if len(P)<3: continue
        pts.append(P); ids.append(np.full(len(P),r.bodyId)); S.append(sm); okid.append(r.bodyId)
    if len(okid)<15: continue
    S=np.vstack(S); prim=Q.set_index('bodyId').loc[okid,'primary'].values
    # soma metrics: distance from cluster centroid; depth along axis from neurite-bundle centroid to soma centroid (positive = away from neuropil)
    cs=S.mean(0); cb=np.vstack(pts).mean(0); axv=cs-cb; axv/=np.linalg.norm(axv); dep=(S-cs)@axv; dc=np.linalg.norm(S-cs,axis=1)
    try: bd=bundle_depth(np.vstack(pts),np.concatenate(ids),slab=3.0,min_cov=0.6)
    except Exception as e: bd=pd.DataFrame()
    df=pd.DataFrame({'bodyId':okid,'hl':hl,'side':side,'primary':prim,'soma_dist':dc,'soma_depth':dep}).set_index('bodyId').join(bd[['depth_norm','radial_norm']] if len(bd) else pd.DataFrame(columns=['depth_norm','radial_norm']))
    rows.append(df.reset_index())
    rec={'hl':hl,'side':side,'n_prim':int(prim.sum()),'n_sec':int((~prim).sum())}
    for col in ['soma_dist','soma_depth','depth_norm','radial_norm']:
        v=df.dropna(subset=[col])
        if v.primary.sum()>=5 and (~v.primary).sum()>=10: rec['AUC_'+col]=round(auc(v[v.primary][col],v[~v.primary][col]),3)
    per.append(rec); print(rec,flush=True)
A=pd.concat(rows); A.to_parquet(f'{R}/hemilineage_geometry.parquet'); P=pd.DataFrame(per); P.to_csv(f'{R}/hemilineage_auc.csv',index=False)
print(P.describe().round(3).to_string())
for col in ['AUC_soma_dist','AUC_soma_depth','AUC_depth_norm','AUC_radial_norm']:
    if col in P: v=P[col].dropna(); print(col,'median',round(v.median(),3),' frac >0.5',round((v>0.5).mean(),3),' frac <0.5',round((v<0.5).mean(),3),' n units',len(v),' sign test p',float(f'{stats.binomtest(int((v>0.5).sum()),len(v)).pvalue:.2g}'))
