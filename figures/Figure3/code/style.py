"""Shared style for Figure 3 panel scripts.

rcParams, colours and drawing helpers are copied verbatim from the code that
produced the approved composite Figure 3 (v4/common.py), plus the two settings
that keep text editable in Adobe Illustrator (pdf.fonttype 42, svg.fonttype none).

Geometry
--------
Every panel is drawn on its own canvas whose size equals the panel's footprint
in the composite (7.2 x 9.0 in, 300 dpi).  Panel boxes are given in composite
pixels at 300 dpi: (x0, y0, width, height), origin at the bottom-left corner of
the composite.  Axes rectangles are given in composite inches (x0, y0, w, h).
Because each panel canvas is a window onto the composite coordinate system,
fonts, line widths, marker sizes and axes sizes are identical to the composite,
so the panel files can be placed in Illustrator at 100 % scale.
"""
import os
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.normpath(os.path.join(HERE, "..", "data"))
PANELS = os.path.normpath(os.path.join(HERE, "..", "panels"))

# ---- rcParams (identical to the composite; savefig.bbox stays 'standard' so canvas size = panel size)
mpl.rcParams.update({
    'font.family': 'Arial', 'font.size': 6.5, 'axes.titlesize': 7, 'axes.labelsize': 6.5,
    'xtick.labelsize': 5.5, 'ytick.labelsize': 5.5,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'xtick.major.size': 2, 'ytick.major.size': 2,
    'axes.spines.top': False, 'axes.spines.right': False,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'savefig.bbox': 'standard', 'svg.hashsalt': 'figure3',
})

# ---- colours (as used in the composite)
INK = '#222222'
COL = dict(ref='#1f3b73', cds='#1f3b73', utr='#8fb3e0', cage='#c23b6e', daft='#4a78c2',
           tss='#7b4fa0', tts='#2a9d8f', body='#4a78c2',
           validation='#f5b041', training='#d9d9d9', chromatin_track='#6a51a3',
           alt='#7f96bd', longer='#e08214', shorter='#f5c27a', nocall='#d9d9d9',
           nochrom='#7fa7d9', grey_call='#aaaaaa', model_full='#555555', model_trunc='#b0b0b0')

FIG_W, FIG_H, DPI = 7.2, 9.0, 300          # composite size (in) and raster resolution

# ---- panel footprints in composite pixels @300 dpi (x0, y0, w, h); origin bottom-left
BOX = {
    'A': (159, 2277, 1209, 340),
    'B': (1490, 2277, 621, 368),
    'C': (104, 1228, 2025, 916),
    'legend_key': (213, 1156, 1941, 44),
    'D': (0, 578, 1043, 540),
    'E': (1238, 575, 828, 555),
    'F': (83, 31, 936, 497),
    'G': (1138, 56, 957, 459),
}


class Panel:
    """A figure whose canvas is the panel's window onto the composite."""

    def __init__(self, key):
        x0, y0, w, h = BOX[key]
        self.key = key
        self.x0, self.y0, self.w, self.h = x0 / DPI, y0 / DPI, w / DPI, h / DPI   # inches
        self.fig = plt.figure(figsize=(self.w, self.h))

    def rect(self, r):
        """composite-inch rectangle (x0, y0, w, h) -> panel figure fraction."""
        x, y, w, h = r
        return [(x - self.x0) / self.w, (y - self.y0) / self.h, w / self.w, h / self.h]

    def axes(self, r):
        return self.fig.add_axes(self.rect(r))

    def text(self, x, y, s, **kw):
        """Figure-level text at composite-inch position (x, y)."""
        return self.fig.text((x - self.x0) / self.w, (y - self.y0) / self.h, s, **kw)

    def save(self, name=None):
        name = name or f"Figure3_{self.key}"
        os.makedirs(PANELS, exist_ok=True)
        base = os.path.join(PANELS, name)
        self.fig.savefig(base + ".pdf", metadata={'CreationDate': None, 'Producer': None})
        self.fig.savefig(base + ".svg", metadata={'Date': None})
        self.fig.savefig(base + ".png", dpi=DPI)
        plt.close(self.fig)
        return base


def data(name):
    return os.path.join(DATA, name)


# ---- drawing helpers (verbatim from the composite code)
def box(ax, x, y, w, h, fc='white', ec=INK, lw=0.8, ls='-'):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.0", fc=fc, ec=ec, lw=lw, ls=ls))


def arrow(ax, a, b, c=INK):
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle='-|>', mutation_scale=7, color=c, lw=0.8))


def title(ax, t, y=1.03):
    ax.text(0, y, t, transform=ax.transAxes, fontsize=7, va='bottom', ha='left')
