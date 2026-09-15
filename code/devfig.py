"""Shared styling for Development (Company of Biologists) figures.
Specs applied: max 180 x 210 mm, panel letters 12 pt bold uppercase, other labels 8 pt Arial,
RGB colour, editable vector text (svg.fonttype none, pdf.fonttype 42)."""
import sys, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
ROOT=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
try:  # optional render-time panel-alignment gate (QA tool); skipped when the module is absent
    from audit_panel_alignment import require_matplotlib_panel_alignment
except ImportError:
    def require_matplotlib_panel_alignment(fig,**kw): print('audit_panel_alignment not installed; alignment gate skipped')
plt.rcParams['font.family']='sans-serif'; plt.rcParams['font.sans-serif']=['Arial','Helvetica','DejaVu Sans','Liberation Sans']
plt.rcParams['svg.fonttype']='none'; plt.rcParams['pdf.fonttype']=42
plt.rcParams.update({'mathtext.fontset':'custom','mathtext.rm':'Arial','mathtext.it':'Arial:italic','mathtext.bf':'Arial:bold','mathtext.default':'regular'})
plt.rcParams.update({'font.size':8,'axes.labelsize':8,'axes.titlesize':8,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,
    'axes.spines.right':False,'axes.spines.top':False,'axes.linewidth':0.7,'xtick.major.width':0.7,'ytick.major.width':0.7,
    'xtick.major.size':3,'ytick.major.size':3,'legend.frameon':False,'lines.linewidth':1.0})
MM=1/25.4
OUT=os.path.join(ROOT,'figures'); QA=os.path.join(ROOT,'qa'); os.makedirs(OUT,exist_ok=True); os.makedirs(QA,exist_ok=True)
C={'DN':'#0F4D92','DN_light':'#7FA6D6','AN':'#D9782D','AN_light':'#F0B98A','grey':'#767676','light':'#E6E6E6'}
# birth-order palette, oldest -> youngest (grey -> blue -> orange -> dark red); no red/green pairing
KC_COLS=['#D9D9D9','#B8B8B8','#7FB3D5','#F4B183','#E58E5C','#C9472B','#7A1F14']
KC_LAB={0:'γd',1:'γm',2:"α′β′",3:'αβp',4:'αβs',5:'αβm',6:'αβc'}
def panel_labels_tight(fig,axes,letters,dx_pt=0,dy_pt=3,fontsize=12):
    fig.canvas.draw(); r=fig.canvas.get_renderer()
    for ax,l in zip(axes,letters):
        bb=ax.get_tightbbox(r).transformed(fig.transFigure.inverted())
        fig.text(max(0.004,bb.x0+dx_pt/72/fig.get_figwidth()),min(1-14/72/fig.get_figheight(),bb.y1+dy_pt/72/fig.get_figheight()),l,fontsize=fontsize,fontweight='bold',ha='left',va='bottom')
def gate(fig,name,**kw):
    base=f'{OUT}/{name}'
    require_matplotlib_panel_alignment(fig,json_out=f'{QA}/{name}.alignment.json',overlay_svg=f'{QA}/{name}.alignment.svg',tolerance_pt=1.5,gutter_tolerance_pt=1.5,require_panel_labels=False,strict=True,**kw)
    return base
def export(fig,base):
    fig.savefig(base+'.pdf'); fig.savefig(base+'.svg'); fig.savefig(base+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'}); fig.savefig(base+'.png',dpi=300); print('saved',base)
def plain_log(ax,axis='both'):
    fmt=matplotlib.ticker.FuncFormatter(lambda v,p: f'{v:g}')
    if axis in ('x','both'): ax.xaxis.set_major_formatter(fmt)
    if axis in ('y','both'): ax.yaxis.set_major_formatter(fmt)
def boxes(ax,data,positions,colors,width=0.6):
    bp=ax.boxplot(data,positions=positions,widths=width,whis=(5,95),showfliers=False,patch_artist=True,medianprops=dict(color='k',lw=0.9),whiskerprops=dict(lw=0.6),capprops=dict(lw=0.6))
    for p,c in zip(bp['boxes'],colors): p.set_facecolor(c); p.set_edgecolor('k'); p.set_linewidth(0.6)
    return bp
