"""Consequences of birth epoch: connectivity (synapse counts) and sexual dimorphism, primary vs secondary within hemilineages."""
import os, pandas as pd, numpy as np, json, pyarrow.feather as pf, pyarrow.compute as pc, pyarrow as pa
from scipy import stats
D=os.environ.get('FLYBRAIN_DATA', os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')); R=os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
a=pd.read_feather(f'{D}/body-annotations-male-cns-v1.0-minconf-0.5.feather')
a['hl']=np.where(a.trumanHl.notna()&(a.trumanHl!='TBD'),'VNC:'+a.trumanHl.astype(str),np.where(a.itoleeHl.notna(),'CB:'+a.itoleeHl.astype(str).str.replace('__prim$','',regex=True),None))
a['primary']=np.where(a.hl.astype(str).str.startswith('VNC'),a.birthtime=='early',a.itoleeHl.astype(str).str.endswith('__prim'))
a=a[a.hl.notna()&~a.itoleeHl.isin(['putative_primary','primary'])&a.superclass.isin(['cb_intrinsic','vnc_intrinsic','vnc_motor','descending_neuron','ascending_neuron','visual_projection','visual_centrifugal'])]
t=pf.read_table(f'{D}/body-stats-male-cns-v1.0-minconf-0.5.feather',columns=['body','pre','post','downstream','synweight'])
bs=t.filter(pc.is_in(t['body'],value_set=pa.array(a.bodyId.astype('int64').tolist()))).to_pandas()
a=a.merge(bs,left_on='bodyId',right_on='body',how='inner'); a['syn']=a.pre+a.post
g=a.groupby('hl').agg(n=('bodyId','size'),prim=('primary','sum')); units=g[(g.prim>=5)&(g.n-g.prim>=10)].index; a=a[a.hl.isin(units)]
print('hemilineages',len(units),' neurons',len(a),' primary',int(a.primary.sum()))
out={}
for col in ['pre','post','syn','downstream']:
    ratios=[]; aucs=[]
    for hl,Q in a.groupby('hl'):
        p=Q[Q.primary][col]; s=Q[~Q.primary][col]
        if len(p)>=5 and len(s)>=10: ratios.append(np.log2((p.median()+1)/(s.median()+1))); aucs.append(stats.mannwhitneyu(p,s).statistic/(len(p)*len(s)))
    ratios=np.array(ratios); aucs=np.array(aucs)
    out[col]={'n_hl':int(len(ratios)),'median_log2_ratio_prim_over_sec':round(float(np.median(ratios)),3),'frac_hl_prim_larger':round(float((ratios>0).mean()),3),'median_AUC':round(float(np.median(aucs)),3),'wilcoxon_p':float(f'{stats.wilcoxon(ratios).pvalue:.2g}')}
    print(col,out[col])
# rich club proxy: top 5% synweight overall -> fraction primary vs expected
thr=a.synweight.quantile(.95); top=a[a.synweight>=thr]; out['top5pct_synweight']={'frac_primary_in_top5':round(float(top.primary.mean()),3),'frac_primary_overall':round(float(a.primary.mean()),3),'n_top':int(len(top))}
print(out['top5pct_synweight'])
# dimorphism vs birth epoch (brain hemilineages only, where primary is annotated by suffix)
b=a[a.hl.str.startswith('CB')].copy(); b['dim']=b.dimorphism.fillna('isomorphic').replace({'potentially male-specific':'male-specific','potentially sexually dimorphic':'sexually dimorphic'})
ct=pd.crosstab(b.dim,b.primary); print(ct.to_string()); 
fr=ct.div(ct.sum(1),axis=0); out['dimorphism_frac_primary']=fr[True].round(3).to_dict()
from scipy.stats import fisher_exact
if 'male-specific' in ct.index: out['fisher_male_specific_vs_iso']=float(f"{fisher_exact(ct.loc[['isomorphic','male-specific']].values)[1]:.2g}")
json.dump(out,open(f'{R}/birth_consequences.json','w'),indent=1); print(json.dumps(out,indent=1))
