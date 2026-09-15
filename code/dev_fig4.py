from devfig import *
from strat import D
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','Helvetica','DejaVu Sans','Liberation Sans']
matplotlib.rcParams.update({"svg.fonttype": "none", "pdf.fonttype": 42})
fig_width_mm = 180
import matplotlib.gridspec as gridspec, json, pyarrow.feather as pf, pyarrow.compute as pc, pyarrow as pa
R=os.path.join(ROOT,'results'); J=json.load(open(f'{R}/birth_consequences.json')); DM=json.load(open(f'{R}/dimorphism_birth.json'))
a=pd.read_feather(f'{D}/body-annotations-male-cns-v1.0-minconf-0.5.feather')
a['hl']=np.where(a.trumanHl.notna()&(a.trumanHl!='TBD'),'VNC:'+a.trumanHl.astype(str),np.where(a.itoleeHl.notna(),'CB:'+a.itoleeHl.astype(str).str.replace('__prim$','',regex=True),None))
a['primary']=np.where(a.hl.astype(str).str.startswith('VNC'),a.birthtime=='early',a.itoleeHl.astype(str).str.endswith('__prim'))
a=a[a.hl.notna()&~a.itoleeHl.isin(['putative_primary','primary'])&a.superclass.isin(['cb_intrinsic','vnc_intrinsic','vnc_motor','descending_neuron','ascending_neuron','visual_projection','visual_centrifugal'])]
t=pf.read_table(f'{D}/body-stats-male-cns-v1.0-minconf-0.5.feather',columns=['body','pre','post','synweight'])
bs=t.filter(pc.is_in(t['body'],value_set=pa.array(a.bodyId.astype('int64').tolist()))).to_pandas(); a=a.merge(bs,left_on='bodyId',right_on='body'); a['syn']=a.pre+a.post
g=a.groupby('hl').agg(n=('bodyId','size'),prim=('primary','sum')); units=g[(g.prim>=5)&(g.n-g.prim>=10)].index; a=a[a.hl.isin(units)]
print('neurons in analysed hemilineages',len(a),'hemilineages',len(units))
fig=plt.figure(figsize=(180*MM,130*MM)); gs=gridspec.GridSpec(2,2,wspace=0.42,hspace=0.55,left=0.1,right=0.985,top=0.95,bottom=0.1)
axA,axB,axC,axD=[fig.add_subplot(gs[i//2,i%2]) for i in range(4)]
med=a.groupby(['hl','primary']).syn.median().unstack(); axA.scatter(med[False],med[True],s=14,color=C['DN'],lw=0); lim=[50,30000]; axA.plot(lim,lim,'k--',lw=0.6); axA.set_xscale('log'); axA.set_yscale('log'); plain_log(axA); axA.set_xlim(lim); axA.set_ylim(lim)
axA.set_xlabel('Larval-born neurons: median synapses'); axA.set_ylabel('Embryonic-born neurons:\nmedian synapses'); axA.text(0.04,0.96,f"{len(med)} hemilineages\nmedian ratio {2**J['syn']['median_log2_ratio_prim_over_sec']:.1f}×",transform=axA.transAxes,va='top',fontsize=8)
a['pct']=a.synweight.rank(pct=True); bins=np.linspace(0,1,21); a['bin']=pd.cut(a.pct,bins,labels=False); fr=a.groupby('bin').primary.mean(); axB.plot((bins[:-1]+bins[1:])/2*100,fr.values,'o-',ms=4,lw=1.0,color=C['DN']); axB.axhline(a.primary.mean(),color=C['grey'],lw=0.6,ls='--'); axB.set_xlabel('Connectivity percentile (synaptic weight)'); axB.set_ylabel('Fraction embryonic-born'); axB.set_ylim(0,1); axB.text(2,a.primary.mean()+0.03,'All neurons',fontsize=8,color=C['grey'])
ratios={}  # log2 of (median+1)/(median+1): pseudocount of 1 guards against zero
# jitter below is for display only
for col in ['pre','post']:
    rr=[]
    for hl,Q in a.groupby('hl'):
        p=Q[Q.primary][col]; s=Q[~Q.primary][col]; rr.append(np.log2((p.median()+1)/(s.median()+1)))
    ratios[col]=np.array(rr)
boxes(axC,[ratios['pre'],ratios['post']],[0,1],[C['DN_light']]*2,width=0.5)
rng=np.random.default_rng(1)
for i,k in enumerate(['pre','post']): axC.scatter(i+rng.uniform(-0.15,0.15,len(ratios[k])),ratios[k],s=5,color=C['grey'],alpha=0.6,lw=0,zorder=3)
axC.axhline(0,color='k',lw=0.6,ls='--'); axC.set_xticks([0,1]); axC.set_xticklabels(['Presynapses','Postsynapses']); axC.set_ylabel('log$_2$ ratio per hemilineage\n(embryonic-born / larval-born)')
cb=DM['frac_primary']; vn=DM['vnc_frac_primary']; cats=['isomorphic','male-specific']
x=np.arange(2); axD.bar(x-0.18,[cb[c] for c in cats],width=0.34,color=C['DN'],label='Central brain'); axD.bar(x+0.18,[vn[c] for c in cats],width=0.34,color=C['AN'],label='Nerve cord')
axD.set_xticks(x); axD.set_xticklabels(['Isomorphic','Male-specific']); axD.set_ylabel('Fraction embryonic-born'); axD.legend(fontsize=8,loc='upper left',borderaxespad=0.2); axD.set_ylim(0,0.8)
panel_labels_tight(fig,[axA,axB,axC,axD],'ABCD')
require_matplotlib_panel_alignment(fig,json_out=f'{QA}/Fig4.alignment.json',overlay_svg=f'{QA}/Fig4.alignment.svg',tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=False,strict=True); base=f'{OUT}/Fig4'; fig.savefig(base+'.pdf'); fig.savefig(base+'.svg'); fig.savefig(base+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'}); fig.savefig(base+'.png',dpi=300); print('saved',base)
