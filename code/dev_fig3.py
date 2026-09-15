from devfig import *
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','Helvetica','DejaVu Sans','Liberation Sans']
matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})
fig_width_mm = 180
import matplotlib.gridspec as gridspec
from strat import D
R=os.path.join(ROOT,'results'); A=pd.read_parquet(f'{R}/hemilineage_geometry.parquet'); U=pd.read_csv(f'{R}/hemilineage_auc.csv')
fig=plt.figure(figsize=(180*MM,135*MM)); gs=gridspec.GridSpec(2,2,wspace=0.42,hspace=0.6,left=0.09,right=0.985,top=0.93,bottom=0.1)
axA,axB,axC,axD=[fig.add_subplot(gs[i//2,i%2]) for i in range(4)]
u=U.sort_values('n_sec',ascending=False).iloc[0]; ex=A[(A.hl==u.hl)&(A.side==u.side)]
ann=pd.read_feather(f'{D}/body-annotations-male-cns-v1.0-minconf-0.5.feather'); sl=dict(zip(ann.bodyId,ann.somaLocation))
S=np.vstack([np.array(sl[b],dtype=float)*0.008 for b in ex.bodyId]); S-=S.mean(0); pc=np.linalg.svd(S,full_matrices=False)[2][:2]; P2=S@pc.T
axA.scatter(P2[~ex.primary.values,0],P2[~ex.primary.values,1],s=12,color=C['light'],edgecolor=C['grey'],lw=0.4,label='Secondary (larval-born)'); axA.scatter(P2[ex.primary.values,0],P2[ex.primary.values,1],s=16,color=C['DN'],edgecolor='k',lw=0.4,label='Primary (embryonic-born)')
axA.set_aspect('equal',adjustable='datalim'); axA.set_xlabel('Soma position, PC1 (µm)'); axA.set_ylabel('Soma position, PC2 (µm)')
axA.legend(fontsize=8,loc='lower center',bbox_to_anchor=(0.5,1.0),ncol=2,handletextpad=0.3,borderaxespad=0.1,columnspacing=1.0,title=f"Hemilineage {u.hl.replace('VNC:','')}, {'left' if u.side=='L' else 'right'} side",title_fontsize=8)
cols=[('AUC_soma_dist','Centroid\ndistance'),('AUC_soma_depth','Axis\nposition'),('AUC_radial_norm','Neurite\nradial'),('AUC_depth_norm','Neurite\nburial')]
data=[U[c].dropna() for c,_ in cols]; print('units before',len(U),'after dropping missing per metric',[len(d) for d in data])
# jitter below is for display only; it does not alter data values
boxes(axB,data,range(4),[C['DN_light']]*4)
rng=np.random.default_rng(0)
for i,d in enumerate(data): axB.scatter(i+rng.uniform(-0.2,0.2,len(d)),d,s=5,color=C['grey'],alpha=0.5,lw=0,zorder=3)
axB.axhline(0.5,color='k',lw=0.6,ls='--'); axB.set_xticks(range(4)); axB.set_xticklabels([f'{l}\n(n = {len(d)})' for (_,l),d in zip(cols,data)]); axB.set_ylabel('AUC (primary versus secondary)'); axB.set_ylim(0,1.02)
axB.set_title('155 hemilineage units',fontsize=8)
frac=[(U[c].dropna()>0.5).mean() if c not in ('AUC_soma_depth','AUC_depth_norm') else (U[c].dropna()<0.5).mean() for c,_ in cols]
axC.bar(range(4),frac,color=C['DN'],width=0.6); axC.axhline(0.5,color='k',lw=0.6,ls='--'); axC.set_xticks(range(4)); axC.set_xticklabels([l for _,l in cols]); axC.set_ylabel('Fraction of units with\nprimary neurons outside'); axC.set_ylim(0,1.08)
for i,f in enumerate(frac): axC.text(i,f+0.02,f'{100*f:.0f}%',ha='center',fontsize=8)
A['z']=A.groupby(['hl','side']).soma_dist.transform(lambda v:(v-v.mean())/v.std())
bins=np.linspace(-3,3,31); axD.hist(A[~A.primary].z.dropna(),bins=bins,histtype='step',color=C['grey'],lw=1.2,label='Secondary'); axD.hist(A[A.primary].z.dropna(),bins=bins,histtype='step',color=C['DN'],lw=1.2,label='Primary')
axD.set_xlabel('Soma distance from centroid\n(z-score within unit)'); axD.set_ylabel('Neurons'); axD.legend(fontsize=8,loc='upper right',handletextpad=0.4,borderaxespad=0.2)
panel_labels_tight(fig,[axA,axB,axC,axD],'ABCD')
require_matplotlib_panel_alignment(fig,json_out=f'{QA}/Fig3.alignment.json',overlay_svg=f'{QA}/Fig3.alignment.svg',tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=False,strict=True,axes=[axA,axB,axC,axD],panel_ids=['A','B','C','D']); base=f'{OUT}/Fig3'; fig.savefig(base+'.pdf'); fig.savefig(base+'.svg'); fig.savefig(base+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'}); fig.savefig(base+'.png',dpi=300); print('saved',base)
