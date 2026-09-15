"""Antennal-lobe projection neurons: literature birth order (twin-spot MARCM) versus soma position and neurite-bundle position."""
import pandas as pd, numpy as np, os, json
from scipy import stats
D=os.environ.get('FLYBRAIN_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')); R=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results'); VOX=0.008
a=pd.read_feather(f'{D}/body-annotations-male-cns-v1.0-minconf-0.5.feather')
# adPN birth order (Yu et al. 2010 PLoS Biol): embryonic uniglomerular (15, in order), then larval block A (in order), then larval block B (in order)
emb=['DP1m','VL2p','DA4l','VM4','VA6','VL2a','DC1','VA7l','VA2','DM4','DL5','VM3','DM3','DL4']
larA=['DL1','DA3','DC2','D','VA3','DC3','VA1d','VM7d','VM2','DM6','VA1v']   # VM7 -> VM7d assumed; VA1lm = VA1v
larB=['VC4','VM5d','VM5v','VC3','DL2v','DL2d']
ad_order={g:i for i,g in enumerate(emb+larA+larB)}; ad_epoch={g:('embryonic' if g in emb else 'larval') for g in ad_order}
# lPN birth order (Lin et al. 2012 PLoS Biol): VA4-VC2-VC1-DM1-DM2-VA5-VA7m-DA1-DL3-VM1-DA2-DM5
l_order={g:i for i,g in enumerate(['VA4','VC2','VC1','DM1','DM2','VA5','VA7m','DA1','DL3','VM1','DA2','DM5'])}
pn=a[(a['class']=='ALPN')&a.somaLocation.notna()].copy()
pn['glom']=pn.type.str.replace(r'_(adPN|lPN|lvPN|vPN|ilPN|l2PN|vPN.*)$','',regex=True)
pn['lineage']=np.where(pn.type.str.endswith('_adPN'),'adPN',np.where(pn.type.str.endswith('_lPN'),'lPN','other'))
pn['rank']=np.where(pn.lineage=='adPN',pn.glom.map(ad_order),np.where(pn.lineage=='lPN',pn.glom.map(l_order),np.nan))
pn['epoch']=np.where(pn.lineage=='adPN',pn.glom.map(ad_epoch),None)
pn=pn[pn['rank'].notna()].copy(); print('PNs with literature birth rank:',pn.groupby('lineage').size().to_dict())
print('annotation check: adPN hemilineage by epoch'); print(pd.crosstab(pn[pn.lineage=='adPN'].epoch,pn[pn.lineage=='adPN'].itoleeHl).to_string())
def read_swc(b):
    p=f'{D}/swc/{b}.swc'
    return pd.read_csv(p,sep=r'\s+',comment='#',header=None,names=['id','t','x','y','z','r','p']) if os.path.exists(p) else None
rows=[]
for i,r in pn.iterrows():
    d=read_swc(r.bodyId)
    if d is None: continue
    xyz=d[['x','y','z']].values*VOX; soma=np.array(r.somaLocation,dtype=float)*VOX
    # primary neurite: nodes within 10-30 um geodesic of soma approximated by euclidean shell 10-25 um from soma along the tree (use euclidean for simplicity)
    dd=np.linalg.norm(xyz-soma,axis=1); shell=xyz[(dd>10)&(dd<25)]
    rows.append(dict(bodyId=r.bodyId,type=r.type,glom=r.glom,lineage=r.lineage,rank=float(r['rank']),epoch=r.epoch,side=r.somaSide,hl=r.itoleeHl,soma=soma,shell=shell))
P=pd.DataFrame(rows); print('with skeletons',len(P))
out={}
for lin in ['adPN','lPN']:
    for side in ['L','R']:
        m=(P.lineage==lin)&(P.side==side); Q=P[m]
        if len(Q)<10: continue
        S=np.vstack(Q.soma.values); cs=S.mean(0)
        # bundle: pooled shell nodes -> axis by PCA; radial distance of each neuron's shell nodes to local centroid
        allsh=np.vstack([np.c_[np.full(len(s),b),s] for b,s in zip(Q.bodyId,Q.shell) if len(s)>0]); X=allsh[:,1:]; c=X.mean(0); ax=np.linalg.svd(X-c,full_matrices=False)[2][0]
        t=(X-c)@ax; perp=(X-c)-np.outer(t,ax); df=pd.DataFrame({'b':allsh[:,0],'bin':np.floor(t/3).astype(int),'x':perp[:,0],'y':perp[:,1],'z':perp[:,2]})
        cent=df.groupby('bin')[['x','y','z']].transform('mean'); df['rloc']=np.linalg.norm(df[['x','y','z']].values-cent.values,axis=1)
        per=df.groupby('b').rloc.median(); P.loc[m,'bundle_radial']=P.loc[m,'bodyId'].map(per)
        # soma: distance from cluster centroid, and from the cluster apex (point farthest from AL entry along cluster axis)
        entry=X.mean(0); axv=cs-entry; axv/=np.linalg.norm(axv); dep=(S-cs)@axv
        P.loc[m,'soma_dist_centroid']=np.linalg.norm(S-cs,axis=1); P.loc[m,'soma_depth']=dep; P.loc[m,'soma_tang']=np.linalg.norm((S-cs)-np.outer(dep,axv),axis=1)
for lin in ['adPN','lPN']:
    Q=P[P.lineage==lin]; res={}
    for col in ['bundle_radial','soma_dist_centroid','soma_depth','soma_tang']:
        v=Q.dropna(subset=[col]); rho=stats.spearmanr(v['rank'],v[col]); res[col]={'spearman':round(rho[0],3),'p':float(f'{rho[1]:.1e}'),'n':int(len(v))}
        if lin=='adPN':
            e=v[v.epoch=='embryonic'][col]; l=v[v.epoch=='larval'][col]; res[col]['AUC_emb_vs_larval']=round(stats.mannwhitneyu(e,l).statistic/(len(e)*len(l)),3)
            lv=v[v.epoch=='larval']; res[col]['spearman_within_larval']=round(stats.spearmanr(lv['rank'],lv[col])[0],3)
    out[lin]=res
print(json.dumps(out,indent=1)); json.dump(out,open(f'{R}/pn_geometry_stats.json','w'),indent=1)
P.drop(columns=['soma','shell']).to_parquet(f'{R}/pn_geometry.parquet')
