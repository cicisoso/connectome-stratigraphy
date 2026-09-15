from devfig import *
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','Helvetica','DejaVu Sans','Liberation Sans']
matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})
fig_width_mm = 180
import matplotlib.gridspec as gridspec, json
from matplotlib.lines import Line2D
from strat import tree, primary_neurite, VOX, D
R=os.path.join(ROOT,'results'); K=pd.read_parquet(f'{R}/kc_geometry.parquet'); J=json.load(open(f'{R}/kc_geometry_stats.json'))
a=pd.read_feather(f'{D}/body-annotations-male-cns-v1.0-minconf-0.5.feather'); soma=dict(zip(a.bodyId,a.somaLocation))
CL={k:KC_COLS[k] for k in range(7)}; lab=KC_LAB
side='R'; Ks=K[K.side==side]; pts=[]; ids=[]
for r in Ks.itertuples():
    xyz,adj=tree(r.bodyId); sm=np.array(soma[r.bodyId],dtype=float)*VOX; P,path,dist=primary_neurite(xyz,adj,sm,lo=r.ped_start,hi=r.ped_start+r.ped_len); pts.append(P); ids.append(np.full(len(P),r.bodyId))
X=np.vstack(pts); ids=np.concatenate(ids); c=X.mean(0); ax_=np.linalg.svd(X-c,full_matrices=False)[2][0]; t=(X-c)@ax_; perp=(X-c)-np.outer(t,ax_)
e1=np.linalg.svd(perp,full_matrices=False)[2][0]; e2=np.cross(ax_,e1); u=perp@e1; v=perp@e2
S=np.vstack([np.array(soma[i],dtype=float)*VOX for i in Ks.bodyId])
if ((S-c)@ax_).mean()>0: t=-t
m=(t>-2)&(t<2)
rk=dict(zip(K.bodyId,K['rank'])); rr=np.array([rk[i] for i in ids[m]])
fig=plt.figure(figsize=(180*MM,140*MM)); gs=gridspec.GridSpec(2,2,wspace=0.42,hspace=0.62,left=0.09,right=0.985,top=0.94,bottom=0.1)
axA,axB,axC,axD=[fig.add_subplot(gs[i//2,i%2]) for i in range(4)]
order=[0,1,2,3,4,5,6]
for k in [1,0,2,4,3,5,6]:
    mm=rr==k; axA.scatter(u[m][mm],v[m][mm],s=2.5,color=CL[k],lw=0,rasterized=True)
axA.set_aspect('equal',adjustable='datalim'); axA.set_xticks([]); axA.set_yticks([]); [s.set_visible(False) for s in axA.spines.values()]
x0=u[m].min()+1; y0=v[m].min()-4
axA.plot([x0,x0+10],[y0,y0],color='k',lw=1.5); axA.text(x0+12,y0,'10 µm',ha='left',va='center',fontsize=8)
axA.set_title('Peduncle cross-section, right hemisphere',fontsize=8)
axA.legend(handles=[Line2D([],[],marker='o',ls='',color=CL[k],markersize=5,label=lab[k]) for k in order],loc='upper center',bbox_to_anchor=(0.5,-0.02),ncol=4,fontsize=8,handletextpad=0.3,columnspacing=0.9,frameon=False,title='Kenyon cell class (oldest to youngest)',title_fontsize=8)
def box(ax,col,ranks,ylabel,title):
    data=[K[K['rank']==k][col].dropna() for k in ranks]; print(col,'cells before',int(K['rank'].isin(ranks).sum()),'after dropping missing',sum(len(d) for d in data))
    boxes(ax,data,range(len(ranks)),[CL[k] for k in ranks])
    ax.set_xticks(range(len(ranks))); ax.set_xticklabels([lab[k] for k in ranks]); ax.set_ylabel(ylabel); ax.set_title(title,fontsize=8); ax.set_xlabel('Birth order →')
box(axB,'radial3_norm',order,'Radial position in peduncle\n(fraction of bundle radius)',f"Peduncle: ρ = {J['radial3_norm']['spearman']:.2f}, n = {J['radial3_norm']['n']:,}")
box(axC,'lobe_radial_norm',[3,4,5,6],'Radial position in α/β lobe\n(fraction of bundle radius)',f"α/β lobe: ρ = {J['lobe_radial_norm']['spearman']:.2f}, n = {J['lobe_radial_norm']['n']:,}")
box(axD,'soma_tang_clone',order,'Soma distance from clone axis (µm)',f"Somata: ρ = {J['soma_tang_clone']['spearman']:.2f}, n = {J['soma_tang_clone']['n']:,}")
panel_labels_tight(fig,[axA,axB,axC,axD],'ABCD')
require_matplotlib_panel_alignment(fig,json_out=f'{QA}/Fig1.alignment.json',overlay_svg=f'{QA}/Fig1.alignment.svg',tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=False,strict=True,axes=[axA,axB,axC,axD],panel_ids=['A','B','C','D']); base=f'{OUT}/Fig1'; fig.savefig(base+'.pdf'); fig.savefig(base+'.svg'); fig.savefig(base+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'}); fig.savefig(base+'.png',dpi=300); print('saved',base)
