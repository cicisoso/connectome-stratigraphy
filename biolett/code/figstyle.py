"""Shared figure style for the Biology Letters package.
Royal Society artwork: figures at final size (double column 174 mm wide, at most 230 mm tall), sans-serif text of 7-8 pt,
panel labels (a), (b), ... in bold, RGB colour, vector PDF/SVG with editable text and 600 dpi TIFF; colour-blind-safe palette.
"""
import os, sys, json, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
try:   # optional panel-alignment gate; skipped when the module is not installed
    sys.path.insert(0, os.environ.get('NATURE_FIGURE_SCRIPTS', os.path.expanduser('~/.claude/skills/nature-figure/scripts')))
    from audit_panel_alignment import require_matplotlib_panel_alignment
except ImportError:
    def require_matplotlib_panel_alignment(fig, **kw): print('audit_panel_alignment not installed; alignment gate skipped')
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'], 'svg.fonttype': 'none',
                     'pdf.fonttype': 42, 'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
                     'mathtext.bf': 'Arial:bold', 'mathtext.default': 'regular', 'font.size': 7, 'axes.labelsize': 7, 'axes.titlesize': 7,
                     'xtick.labelsize': 7, 'ytick.labelsize': 7, 'legend.fontsize': 7, 'axes.spines.right': False, 'axes.spines.top': False,
                     'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'xtick.major.size': 2.5,
                     'ytick.major.size': 2.5, 'legend.frameon': False, 'lines.linewidth': 0.9})
MM = 1 / 25.4
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE)
RES = f'{ROOT}/results'; OUT = f'{ROOT}/figures'; QA = f'{ROOT}/figures/qa'
os.makedirs(OUT, exist_ok=True); os.makedirs(QA, exist_ok=True)
DS = ['MaleCNS', 'FlyWire', 'BANC']
DSC = {'MaleCNS': '#0072B2', 'FlyWire': '#E69F00', 'BANC': '#009E73', 'MANC': '#56B4E9'}
DSL = {'MaleCNS': 'MaleCNS (♂ CNS)', 'FlyWire': 'FlyWire (♀ brain)', 'BANC': 'BANC (♀ CNS)', 'MANC': 'MANC (♂ VNC)'}
PRIM = '#D55E00'; SEC = '#9A9A9A'; GREY = '#767676'; LIGHT = '#E6E6E6'
KC_MAJOR_COL = {0: '#B9B9B9', 1: '#56B4E9', 2: '#CC79A7'}
KC_MAJOR_LAB = {0: 'γ', 1: "α′β′", 2: 'αβ'}
RANK_CMAP = plt.cm.viridis

def panel_labels(fig, axes, letters, dx_pt=0, dy_pt=3, fontsize=8, align_rows=True):
    """Bold (a), (b), ... above the top-left corner of each panel's tight box; panels in one grid row share a baseline."""
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    bbs = [ax.get_tightbbox(r).transformed(fig.transFigure.inverted()) for ax in axes]; tops = [b.y1 for b in bbs]
    if align_rows:
        y1 = [ax.get_position().y1 for ax in axes]
        tops = [max(bbs[j].y1 for j in range(len(axes)) if abs(y1[j] - y1[i]) < 0.03) for i in range(len(axes))]
    for bb, top, l in zip(bbs, tops, letters):
        fig.text(max(0.004, bb.x0 + dx_pt / 72 / fig.get_figwidth()), min(1 - 10 / 72 / fig.get_figheight(), top + dy_pt / 72 / fig.get_figheight()),
                 f'({l})', fontsize=fontsize, fontweight='bold', ha='left', va='bottom')

def gate(fig, name, **kw):
    require_matplotlib_panel_alignment(fig, json_out=f'{QA}/{name}.alignment.json', overlay_svg=f'{QA}/{name}.alignment.svg',
                                       tolerance_pt=1.5, gutter_tolerance_pt=1.5, require_panel_labels=False, strict=True, **kw)

def export(fig, name):
    base = f'{OUT}/{name}'
    fig.savefig(base + '.pdf'); fig.savefig(base + '.svg'); fig.savefig(base + '.png', dpi=300)
    fig.savefig(base + '.tiff', dpi=600, pil_kwargs={'compression': 'tiff_lzw'}); print('saved', base)

def boxes(ax, data, positions, colors, width=0.6, lw=0.5):
    bp = ax.boxplot(data, positions=positions, widths=width, whis=(5, 95), showfliers=False, patch_artist=True,
                    medianprops=dict(color='k', lw=0.8), whiskerprops=dict(lw=lw), capprops=dict(lw=lw))
    for p, c in zip(bp['boxes'], colors): p.set_facecolor(c); p.set_edgecolor('k'); p.set_linewidth(lw)
    return bp

def strip(ax, values, x, color, width=0.28, s=4, alpha=0.75, seed=0, marker='o', facecolor=None):
    v = np.asarray(values, float); v = v[np.isfinite(v)]; rng = np.random.default_rng(seed)
    ax.scatter(x + rng.uniform(-width / 2, width / 2, len(v)), v, s=s, color=facecolor or color, edgecolors=color if facecolor else 'none',
               lw=0.4 if facecolor else 0, alpha=alpha, marker=marker, zorder=3)

def scalebar(ax, x0, y0, length, label, lw=1.2, fontsize=7, dy=0.0):
    ax.plot([x0, x0 + length], [y0, y0], color='k', lw=lw, solid_capstyle='butt')
    ax.text(x0 + length / 2, y0 + dy, label, ha='center', va='top', fontsize=fontsize)

def clean_3d_axes(ax):
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values(): s.set_visible(False)
