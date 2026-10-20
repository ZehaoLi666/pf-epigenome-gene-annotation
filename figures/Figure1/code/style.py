"""
Shared style for the Figure 1 panel scripts.

rcParams, colours and helper functions are copied verbatim from the approved generating code
(v4/common.py = annotation/common.py).  Only additions:
  * svg.fonttype='none' / pdf.fonttype=42 / ps.fonttype=42 so text stays editable in Illustrator
  * PANELS: per-panel canvas size and the exact axes rectangles measured from the approved
    composite (7.2 x 10.9 in, 300 dpi), so every panel is drawn at 100 % of its size in the composite.
"""
import os
import matplotlib as mpl
mpl.use('Agg')
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyBboxPatch, FancyArrowPatch

mpl.rcParams.update({'font.family': 'Arial', 'font.size': 6.5, 'axes.titlesize': 7, 'axes.labelsize': 6.5, 'xtick.labelsize': 5.5, 'ytick.labelsize': 5.5,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5, 'xtick.major.size': 2, 'ytick.major.size': 2,
    'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})

INK = '#222222'
COL = dict(agree='#4a78c2', ol_unit='#a6400a', ol_lnc='#e08214', ol_anti='#f5b041', ol_none='#fbd9a6', short='#7b4fa0', noutr='#2a9d8f',
           between='#c4a46b', na='#e3e3e3', ref='#1f3b73', cds='#1f3b73', utr='#8fb3e0', lnc='#e08214', cage='#c23b6e', daft='#4a78c2',
           tss='#7b4fa0', tts='#2a9d8f', body='#4a78c2', inter='#b5b5b5', trunc='#9a9a9a', art='#d9822b')

def title(ax, t, y=1.03): ax.text(0, y, t, transform=ax.transAxes, fontsize=7, va='bottom', ha='left')
def box(ax, x, y, w, h, fc='white', ec=INK, lw=0.8, ls='-'): ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.0", fc=fc, ec=ec, lw=lw, ls=ls))
def arrow(ax, a, b, c=INK): ax.add_patch(FancyArrowPatch(a, b, arrowstyle='-|>', mutation_scale=7, color=c, lw=0.8))

# ------------------------------------------------------------------ geometry
# size_px: panel canvas in pixels at 300 dpi (size in inches = px / 300; mm = px / 300 * 25.4).
# origin_px_in_composite: top-left corner of this canvas in the approved composite Figure1.png
# (7.2 x 10.9 in = 2160 x 3270 px; pixels from the top-left).  axes_px: [left, bottom, width, height]
# of each axes, in pixels from the canvas's bottom-left corner, measured from the composite's GridSpec.
# Panels in one figure row share the same canvas top and height, so aligning their top edges in
# Illustrator reproduces the composite.
DPI = 300
PANELS = {
    'A': dict(size_px=(1417, 771), origin_px_in_composite=(26, 31),
        axes_px=[[200.79999999999998, 591.7650352469277, 1203.84, 112.13496475307228], [506.2519402985075, 97.16974428248295, 898.3880597014925, 427.1808181069425]]),
    'B': dict(size_px=(836, 771), origin_px_in_composite=(1304, 31),
        axes_px=[[316.72, 498.22872687541803, 506.88000000000005, 205.67127312458194], [316.72, 97.16974428248295, 506.88000000000005, 205.67127312458194]]),
    'C': dict(size_px=(717, 433), origin_px_in_composite=(29, 873),
        axes_px=[[197.79999999999998, 83.45895168384004, 506.88000000000005, 303.3651278587587]]),
    'D': dict(size_px=(669, 433), origin_px_in_composite=(774, 873),
        axes_px=[[149.76, 83.45895168384004, 506.88000000000005, 303.3651278587587]]),
    'E': dict(size_px=(787, 433), origin_px_in_composite=(1353, 873),
        axes_px=[[267.72, 83.45895168384004, 506.88000000000005, 303.3651278587587]]),
    'F': dict(size_px=(1926, 166), origin_px_in_composite=(214, 1393),
        axes_px=[[12.799999999999983, 12.901840914802733, 1900.8, 109.21144602915268]]),
    'G': dict(size_px=(680, 652), origin_px_in_composite=(66, 1671),
        axes_px=[[160.79999999999998, 122.67674077971833, 506.88000000000005, 439.8794353951999]]),
    'H': dict(size_px=(621, 652), origin_px_in_composite=(822, 1671),
        axes_px=[[101.75999999999999, 122.67674077971833, 506.88000000000005, 439.8794353951999]]),
    'I': dict(size_px=(540, 652), origin_px_in_composite=(1600, 1671),
        axes_px=[[20.720000000000027, 122.67674077971833, 506.88000000000005, 439.8794353951999]]),
    'J': dict(size_px=(590, 925), origin_px_in_composite=(332, 2327),
        axes_px=[[12.443688238353502, 728.2069293280505, 564.9156698120217, 109.12414671178364], [12.443688238353502, 504.1204940554029, 564.9156698120217, 204.60777508459427], [12.443688238353502, 382.3379463250524, 564.9156698120217, 102.30388754229713], [12.443688238353502, 199.17306606932362, 564.9156698120217, 163.6862200676754], [12.443688238353502, 63.75000000000007, 564.9156698120217, 115.94440588127014]]),
    'K': dict(size_px=(590, 925), origin_px_in_composite=(941, 2327),
        axes_px=[[12.564009213165946, 728.2069293280505, 564.9156698120217, 109.12414671178364], [12.564009213165946, 504.1204940554029, 564.9156698120217, 204.60777508459427], [12.564009213165946, 382.3379463250524, 564.9156698120217, 102.30388754229713], [12.564009213165946, 199.17306606932362, 564.9156698120217, 163.6862200676754], [12.564009213165946, 63.75000000000007, 564.9156698120217, 115.94440588127014]]),
    'L': dict(size_px=(590, 925), origin_px_in_composite=(1550, 2327),
        axes_px=[[12.684330187978503, 728.2069293280505, 564.9156698120217, 109.12414671178364], [12.684330187978503, 504.1204940554029, 564.9156698120217, 204.60777508459427], [12.684330187978503, 382.3379463250524, 564.9156698120217, 102.30388754229713], [12.684330187978503, 199.17306606932362, 564.9156698120217, 163.6862200676754], [12.684330187978503, 63.75000000000007, 564.9156698120217, 115.94440588127014]]),
    'JKL_track_labels': dict(size_px=(170, 925), origin_px_in_composite=(164, 2327),
        axes_px=[[180.4436882383535, 728.2069293280505, 564.9156698120217, 109.12414671178364], [180.4436882383535, 504.1204940554029, 564.9156698120217, 204.60777508459427], [180.4436882383535, 382.3379463250524, 564.9156698120217, 102.30388754229713], [180.4436882383535, 199.17306606932362, 564.9156698120217, 163.6862200676754], [180.4436882383535, 63.75000000000007, 564.9156698120217, 115.94440588127014]]),
}
for _p in PANELS.values():
    _p['size_in'] = (_p['size_px'][0] / DPI, _p['size_px'][1] / DPI)
    _p['size_mm'] = (_p['size_px'][0] / DPI * 25.4, _p['size_px'][1] / DPI * 25.4)

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')
OUT = os.path.join(HERE, '..', 'panels')

def data(name): return os.path.join(DATA, name)
def read(name, **kw): return pd.read_csv(data(name), sep='\t', **kw)

def panel_figure(key):
    """Figure with the panel's exact canvas size and axes placed at their composite positions."""
    p = PANELS[key]; W, H = p['size_px']
    fig = plt.figure(figsize=(W / DPI, H / DPI), dpi=DPI)
    axs = [fig.add_axes([l / W, b / H, w / W, h / H]) for l, b, w, h in p['axes_px']]
    return fig, axs

def save_panel(fig, key):
    """Figure1_<key>.pdf / .svg (vector, transparent background, editable text) and .png (300 dpi, white)."""
    os.makedirs(OUT, exist_ok=True); stem = os.path.join(OUT, f"Figure1_{key}")
    fig.savefig(stem + '.pdf', dpi=DPI, transparent=True)
    fig.savefig(stem + '.svg', dpi=DPI, transparent=True)
    fig.savefig(stem + '.png', dpi=DPI, facecolor='white')
    plt.close(fig)
    W, H = PANELS[key]['size_mm']; print(f"Figure1_{key}: {W:.1f} x {H:.1f} mm")
