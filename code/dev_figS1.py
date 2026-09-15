from devfig import *
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','Helvetica','DejaVu Sans','Liberation Sans']
matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})
fig_width_mm = 180
import matplotlib.gridspec as gridspec, json
from matplotlib.lines import Line2D
R=os.path.join(ROOT,'results'); K=pd.read_parquet(f'{R}/kc_geometry.parquet'); P=pd.read_parquet(f'{R}/pn_geometry.parquet'); J=json.load(open(f'{R}/kc_geometry_stats.json'))
lab=KC_LAB; cols=KC_COLS
fig=plt.figure(figsize=(180*MM,130*MM)); gs=gridspec.GridSpec(2,2,wspace=0.42,hspace=0.55,left=0.09,right=0.985,top=0.95,bottom=0.1)
axA,axB,axC,axD=[fig.add_subplot(gs[i//2,i%2]) for i in range(4)]
g=K.groupby(['rank','side']).radial3_norm.median().unstack(); axA.scatter(g.L,g.R,s=30,c=[cols[i] for i in g.index],edgecolor='k',lw=0.5); axA.plot([0.2,1],[0.2,1],'k--',lw=0.6)
axA.legend(handles=[Line2D([],[],marker='o',ls='',color=cols[i],markeredgecolor='k',markeredgewidth=0.5,markersize=5,label=lab[i]) for i in g.index],fontsize=8,loc='lower right',ncol=2,handletextpad=0.3,columnspacing=0.7)
axA.set_xlim(0.2,1); axA.set_ylim(0.2,1); axA.set_xlabel('Radial position, left peduncle\n(class median)'); axA.set_ylabel('Radial position, right peduncle\n(class median)')
data=[K[K['rank']==k].depth_norm.dropna() for k in range(7)]; print('KCs before',len(K),'after dropping missing depth',sum(len(d) for d in data))
boxes(axB,data,range(7),[cols[k] for k in range(7)])
axB.set_xticks(range(7)); axB.set_xticklabels([lab[k] for k in range(7)]); axB.set_ylabel('Burial depth in peduncle\n(fraction of bundle radius)'); axB.set_xlabel('Birth order →')
axB.set_title(f"ρ = {J['depth_norm']['spearman']:.2f}, n = {J['depth_norm']['n']:,}",fontsize=8)
ad=P[P.lineage=='adPN'].copy(); ad['annot']=np.where(ad.hl.astype(str).str.endswith('__prim'),'annotated primary','annotated secondary'); ct=pd.crosstab(ad.epoch,ad.annot)
x=np.arange(2); e=ct.loc['embryonic']; l=ct.loc['larval']
axC.bar(x-0.18,[e.get('annotated primary',0),l.get('annotated primary',0)],width=0.34,color=C['DN'],label='Annotated primary'); axC.bar(x+0.18,[e.get('annotated secondary',0),l.get('annotated secondary',0)],width=0.34,color=C['light'],edgecolor=C['grey'],label='Annotated secondary')
for xx,val in zip([x[0]-0.18,x[0]+0.18,x[1]-0.18,x[1]+0.18],[e.get('annotated primary',0),e.get('annotated secondary',0),l.get('annotated primary',0),l.get('annotated secondary',0)]): axC.text(xx,val+1.5,str(val),ha='center',fontsize=8)
axC.set_xticks(x); axC.set_xticklabels(['Embryonic-born\n(MARCM)','Larval-born\n(MARCM)']); axC.set_ylabel('adPN neurons'); axC.set_ylim(0,115); axC.legend(fontsize=8,loc='upper left',borderaxespad=0.2)
axD.hist(K.ped_len,bins=40,color=C['DN_light'],lw=0); axD.set_xlabel('Branch-free peduncle stretch (µm)'); axD.set_ylabel('Kenyon cells')
panel_labels_tight(fig,[axA,axB,axC,axD],'ABCD')
require_matplotlib_panel_alignment(fig,json_out=f'{QA}/FigS1.alignment.json',overlay_svg=f'{QA}/FigS1.alignment.svg',tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=False,strict=True); base=f'{OUT}/FigS1'; fig.savefig(base+'.pdf'); fig.savefig(base+'.svg'); fig.savefig(base+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'}); fig.savefig(base+'.png',dpi=300); print('saved',base); print(ct)
