from devfig import *
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','Helvetica','DejaVu Sans','Liberation Sans']
matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})
fig_width_mm = 180
import matplotlib.gridspec as gridspec, json
R=os.path.join(ROOT,'results'); P=pd.read_parquet(f'{R}/pn_geometry.parquet'); J=json.load(open(f'{R}/pn_geometry_stats.json'))
fig=plt.figure(figsize=(180*MM,135*MM)); gs=gridspec.GridSpec(2,2,wspace=0.4,hspace=0.55,left=0.09,right=0.985,top=0.94,bottom=0.15)
axA,axB,axC,axD=[fig.add_subplot(gs[i//2,i%2]) for i in range(4)]
ad=P[P.lineage=='adPN']; lp=P[P.lineage=='lPN']
def rankplot(ax,df,col,ylabel,title,color):
    g=df.groupby(['glom','rank'])[col].median().reset_index().sort_values('rank'); ax.scatter(df['rank'],df[col],s=9,color=C['light'],lw=0); ax.plot(g['rank'],g[col],'o-',ms=4,lw=1.0,color=color)
    ax.set_xlabel('Birth rank (twin-spot MARCM)'); ax.set_ylabel(ylabel); ax.set_title(title,fontsize=8)
rankplot(axA,ad,'soma_dist_centroid','Soma distance from\ncluster centroid (µm)',f"adPN lineage: ρ = {J['adPN']['soma_dist_centroid']['spearman']:.2f}, n = {J['adPN']['soma_dist_centroid']['n']}",C['DN']); yl=axA.get_ylim(); axA.axvspan(-0.5,13.5,color=C['light'],alpha=0.5,lw=0,zorder=0); axA.set_ylim(yl[0],yl[1]*1.12)
axA.text(6.5,axA.get_ylim()[1]*0.985,'Embryonic-born',ha='center',va='top',fontsize=8,color=C['grey']); axA.text(23,axA.get_ylim()[1]*0.985,'Larval-born',ha='center',va='top',fontsize=8,color=C['grey'])
rankplot(axB,ad,'bundle_radial2','Neurite radial position in bundle\n(fraction of bundle radius)',f"adPN lineage: ρ = {J['adPN']['bundle_radial2']['spearman']:.2f}, n = {J['adPN']['bundle_radial2']['n']}",C['DN'])
rankplot(axC,lp,'soma_depth','Soma position along\ncluster axis (µm)',f"lPN lineage: ρ = {J['lPN']['soma_depth']['spearman']:.2f}, n = {J['lPN']['soma_depth']['n']}",C['AN'])
rankplot(axD,lp,'bundle_radial2','Neurite radial position in bundle\n(fraction of bundle radius)',f"lPN lineage: ρ = {J['lPN']['bundle_radial2']['spearman']:.2f}, n = {J['lPN']['bundle_radial2']['n']}",C['AN'])
from matplotlib.lines import Line2D
fig.legend(handles=[Line2D([],[],marker='o',ls='',color=C['light'],markersize=4,label='Individual neuron'),Line2D([],[],marker='o',ls='-',color=C['DN'],markersize=4,label='Type median, adPN'),Line2D([],[],marker='o',ls='-',color=C['AN'],markersize=4,label='Type median, lPN')],loc='lower center',bbox_to_anchor=(0.5,0.005),ncol=3,fontsize=8,frameon=False,handletextpad=0.4,columnspacing=1.5)
panel_labels_tight(fig,[axA,axB,axC,axD],'ABCD')
require_matplotlib_panel_alignment(fig,json_out=f'{QA}/Fig2.alignment.json',overlay_svg=f'{QA}/Fig2.alignment.svg',tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=False,strict=True); base=f'{OUT}/Fig2'; fig.savefig(base+'.pdf'); fig.savefig(base+'.svg'); fig.savefig(base+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'}); fig.savefig(base+'.png',dpi=300); print('saved',base)
