"""Kenyon cell stratigraphy: radial position in the peduncle and concentric soma position versus birth-order class."""
import pandas as pd, numpy as np, os, json
from collections import defaultdict, deque
from scipy import stats
D=os.environ.get('FLYBRAIN_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')); R=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
a=pd.read_feather(f'{D}/body-annotations-male-cns-v1.0-minconf-0.5.feather')
kc=a[(a['class']=='Kenyon_Cell')&a.somaLocation.notna()].copy()
order={'KCg-d':0,'KCg-m':1,"KCa'b'-ap1":2,"KCa'b'-ap2":2,"KCa'b'-m":2,'KCab-p':3,'KCab-s':4,'KCab-m':5,'KCab-c':6}
kc['rank']=kc.type.map(order); kc=kc[kc['rank'].notna()]
VOX=0.008
def read_swc(b):
    p=f'{D}/swc/{b}.swc'
    if not os.path.exists(p): return None
    d=pd.read_csv(p,sep=r'\s+',comment='#',header=None,names=['id','t','x','y','z','r','p'])
    return d
def soma_to_tip_path(d,soma_xyz):
    idx={i:k for k,i in enumerate(d.id.values)}; xyz=d[['x','y','z']].values*VOX; par=d.p.values
    adj=defaultdict(list)
    for k in range(len(d)):
        if par[k]!=-1 and par[k] in idx: j=idx[par[k]]; adj[k].append(j); adj[j].append(k)
    s=int(np.argmin(np.linalg.norm(xyz-soma_xyz,axis=1)))
    # BFS distances (geodesic) from soma
    dist=np.full(len(d),np.inf); dist[s]=0; prev=np.full(len(d),-1); q=deque([s])
    while q:
        u=q.popleft()
        for v in adj[u]:
            if dist[v]==np.inf: dist[v]=dist[u]+np.linalg.norm(xyz[u]-xyz[v]); prev[v]=u; q.append(v)
    tip=int(np.argmax(np.where(np.isfinite(dist),dist,-1)))
    path=[]; u=tip
    while u!=-1: path.append(u); u=prev[u]
    path=path[::-1]
    deg=np.array([len(adj[k]) for k in range(len(d))])
    return xyz,path,dist,deg,adj
rows=[]; ped_nodes=[]
for i,r in kc.iterrows():
    d=read_swc(r.bodyId)
    if d is None: continue
    soma=np.array(r.somaLocation,dtype=float)*VOX
    xyz,path,dist,deg,adj=soma_to_tip_path(d,soma)
    pxyz=xyz[path]; pdist=dist[path]; L=pdist[-1]
    if L<150: continue
    # branch points along path with substantial side branches (>4 um cable)
    def side_len(u,p_prev,p_next):
        seen={u,p_prev,p_next}; tot=0.0; q=deque([v for v in adj[u] if v not in (p_prev,p_next)])
        for v in q: seen.add(v)
        while q:
            w=q.popleft()
            for v in adj[w]:
                if v not in seen: seen.add(v); tot+=np.linalg.norm(xyz[w]-xyz[v]); q.append(v)
                    
        return tot
    bp=[]
    for j in range(1,len(path)-1):
        if deg[path[j]]>2:
            sl=side_len(path[j],path[j-1],path[j+1])
            if sl>4: bp.append((pdist[j],sl))
    # peduncle = longest branch-free stretch between consecutive branch points, located after the calyx (first branch cluster)
    cuts=[0.0]+[b[0] for b in bp]+[L]
    gaps=[(cuts[k+1]-cuts[k],cuts[k],cuts[k+1]) for k in range(len(cuts)-1)]
    g=max(gaps)
    if g[0]<40: continue
    s0,s1=g[1],g[2]
    sel=(pdist>=s0+0.15*(s1-s0))&(pdist<=s0+0.85*(s1-s0))
    P=pxyz[sel]
    rows.append(dict(bodyId=r.bodyId,type=r.type,rank=int(r['rank']),side=r.somaSide,nb=r.itoleeHl,soma=soma,ped_len=s1-s0,ped_start=s0,path_len=L,n_ped=int(sel.sum())))
    ped_nodes.append(np.c_[np.full(len(P),r.bodyId),P])
K=pd.DataFrame(rows); PN=np.vstack(ped_nodes); print('KCs analysed',len(K),' of',len(kc),'; peduncle length median',round(K.ped_len.median(),1),'um')
# per side: bundle axis by PCA on all peduncle nodes; slabs; radial position normalized
out={}
K['radial']=np.nan; K['radial_norm']=np.nan
for side in ['L','R']:
    ids=set(K[K.side==side].bodyId); m=np.isin(PN[:,0],list(ids)); X=PN[m,1:]; b=PN[m,0]
    c=X.mean(0); u,s,vt=np.linalg.svd(X-c,full_matrices=False); ax=vt[0]
    t=(X-c)@ax; perp=(X-c)-np.outer(t,ax); rad=np.linalg.norm(perp,axis=1)
    # slabs of 4 um along axis; local centroid + radius normalisation
    bins=np.floor(t/4).astype(int); df=pd.DataFrame({'b':b,'bin':bins,'x':perp[:,0],'y':perp[:,1],'z':perp[:,2]})
    cent=df.groupby('bin')[['x','y','z']].transform('mean'); rloc=np.linalg.norm(df[['x','y','z']].values-cent.values,axis=1)
    df['rloc']=rloc; r95=df.groupby('bin').rloc.transform(lambda v: v.quantile(.95)); df['rn']=df.rloc/r95
    cnt=df.groupby('bin').b.transform('nunique'); df=df[cnt>=200]   # keep slabs where most KCs present
    per=df.groupby('b').agg(radial=('rloc','median'),radial_norm=('rn','median'))
    K.loc[K.bodyId.isin(per.index),'radial']=K.bodyId.map(per.radial); K.loc[K.bodyId.isin(per.index),'radial_norm']=K.bodyId.map(per.radial_norm)
# soma geometry per side: axis from calyx-side (peduncle start centroid) to soma cluster centroid; tangential radius and depth
K['soma_x']=[s[0] for s in K.soma]; K['soma_y']=[s[1] for s in K.soma]; K['soma_z']=[s[2] for s in K.soma]
K['soma_tang']=np.nan; K['soma_depth']=np.nan
for side in ['L','R']:
    m=K.side==side; S=K.loc[m,['soma_x','soma_y','soma_z']].values; cs=S.mean(0)
    ids=set(K[m].bodyId); mm=np.isin(PN[:,0],list(ids)); cp=PN[mm,1:].mean(0)
    axv=(cs-cp); axv/=np.linalg.norm(axv)
    dep=(S-cs)@axv; perp=(S-cs)-np.outer(dep,axv); K.loc[m,'soma_depth']=dep; K.loc[m,'soma_tang']=np.linalg.norm(perp,axis=1)
    # per neuroblast clone: tangential radius relative to the clone's own centroid
    for nb in K[m].nb.unique():
        mn=m&(K.nb==nb); Sn=K.loc[mn,['soma_x','soma_y','soma_z']].values; cn=Sn.mean(0); pn=(Sn-cn)-np.outer((Sn-cn)@axv,axv); K.loc[mn,'soma_tang_clone']=np.linalg.norm(pn,axis=1)
K.drop(columns=['soma']).to_parquet(f'{R}/kc_geometry.parquet')
lab={0:'γd',1:'γm',2:"α'β'",3:'αβp',4:'αβs',5:'αβm',6:'αβc'}
summ=K.groupby('rank').agg(n=('bodyId','size'),radial_norm=('radial_norm','median'),radial_um=('radial','median'),soma_tang=('soma_tang','median'),soma_tang_clone=('soma_tang_clone','median'),soma_depth=('soma_depth','median')).round(3); summ.index=[lab[i] for i in summ.index]; print(summ.to_string())
for col in ['radial_norm','soma_tang','soma_tang_clone','soma_depth']:
    v=K.dropna(subset=[col]); rho=stats.spearmanr(v['rank'],v[col]); out[col]={'spearman':round(rho[0],3),'p':float(f'{rho[1]:.1e}'),'n':int(len(v))}
    # AUC for adjacent classes
    aucs={}
    for k in range(6):
        x1=v[v['rank']==k][col]; x2=v[v['rank']==k+1][col]
        if len(x1)>5 and len(x2)>5:
            u=stats.mannwhitneyu(x1,x2).statistic; aucs[f'{lab[k]}>{lab[k+1]}']=round(u/(len(x1)*len(x2)),3)
    out[col]['adjacent_AUC']=aucs
    # coarse: gamma vs a'b' vs ab
    g=v[v['rank']<=1][col]; ab=v[v['rank']>=3][col]; out[col]['AUC_gamma_vs_ab']=round(stats.mannwhitneyu(g,ab).statistic/(len(g)*len(ab)),3)
print(json.dumps(out,indent=1)); json.dump(out,open(f'{R}/kc_geometry_stats.json','w'),indent=1)
