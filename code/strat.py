"""Shared geometry helpers for connectome stratigraphy."""
import pandas as pd, numpy as np, os
from collections import defaultdict, deque
from scipy.spatial import ConvexHull
from scipy import stats
D=os.environ.get('FLYBRAIN_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')); VOX=0.008
def tree(b):
    p=f'{D}/swc/{b}.swc'
    if not os.path.exists(p): return None
    d=pd.read_csv(p,sep=r'\s+',comment='#',header=None,names=['id','t','x','y','z','r','p'])
    idx={i:k for k,i in enumerate(d.id.values)}; xyz=d[['x','y','z']].values*VOX; par=d.p.values; adj=defaultdict(list)
    for k in range(len(d)):
        if par[k]!=-1 and par[k] in idx: j=idx[par[k]]; adj[k].append(j); adj[j].append(k)
    return xyz,adj
def geodesic_from(xyz,adj,s):
    dist=np.full(len(xyz),np.inf); dist[s]=0; prev=np.full(len(xyz),-1); q=deque([s])
    while q:
        u=q.popleft()
        for v in adj[u]:
            if dist[v]==np.inf: dist[v]=dist[u]+np.linalg.norm(xyz[u]-xyz[v]); prev[v]=u; q.append(v)
    return dist,prev
def primary_neurite(xyz,adj,soma_xyz,lo=8.0,hi=30.0):
    """Nodes on the path from the soma outward between geodesic distances lo and hi (the cell-body fibre / primary neurite)."""
    s=int(np.argmin(np.linalg.norm(xyz-soma_xyz,axis=1))); dist,prev=geodesic_from(xyz,adj,s)
    tip=int(np.argmax(np.where(np.isfinite(dist),dist,-1))); path=[]; u=tip
    while u!=-1: path.append(u); u=prev[u]
    path=path[::-1]; pd_=dist[path]; sel=(pd_>=lo)&(pd_<=hi); return xyz[path][sel], path, dist
def bundle_depth(points, ids, slab=4.0, min_cov=0.8, n_total=None):
    """points: (N,3) fibre nodes from many neurons; ids: (N,) neuron id. Fits the bundle axis, slices it, and for each neuron
    returns burial depth (distance to the convex-hull boundary of the slab cross-section) and radial distance to the geometric median,
    both normalised by the slab's 95th-percentile radius. Only slabs containing >= min_cov of the neurons are used."""
    ids=np.asarray(ids); X=np.asarray(points,float); c=X.mean(0); ax=np.linalg.svd(X-c,full_matrices=False)[2][0]
    t=(X-c)@ax; perp=(X-c)-np.outer(t,ax)
    e1=np.linalg.svd(perp,full_matrices=False)[2][0]; e2=np.cross(ax,e1); u=perp@e1; v=perp@e2
    bins=np.floor(t/slab).astype(int); n_total=n_total or len(np.unique(ids))
    rec=defaultdict(list)
    for b_ in np.unique(bins):
        m=bins==b_; nid=np.unique(ids[m])
        if len(nid)<min_cov*n_total or m.sum()<30: continue
        P=np.c_[u[m],v[m]]
        try: hull=ConvexHull(P)
        except Exception: continue
        # distance of each point to hull boundary (min distance to hull edges)
        eq=hull.equations  # a x + b y + c <= 0 inside
        dist_edge=np.min(-(P@eq[:,:2].T+eq[:,2]),axis=1)  # positive inside
        # geometric median (Weiszfeld)
        gm=P.mean(0)
        for _ in range(50):
            w=1/np.maximum(np.linalg.norm(P-gm,axis=1),1e-6); gm=(P*w[:,None]).sum(0)/w.sum()
        rad=np.linalg.norm(P-gm,axis=1); r95=np.quantile(rad,.95)
        for i,dd,rr in zip(ids[m],dist_edge,rad): rec[i].append((dd/r95,rr/r95,dd,rr))
    out=pd.DataFrame({i:np.median(np.array(v),axis=0) for i,v in rec.items()},index=['depth_norm','radial_norm','depth_um','radial_um']).T
    return out
def auc(x1,x2):
    x1=np.asarray(x1); x2=np.asarray(x2)
    return stats.mannwhitneyu(x1,x2).statistic/(len(x1)*len(x2)) if len(x1)>2 and len(x2)>2 else np.nan
