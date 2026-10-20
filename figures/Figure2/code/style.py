"""Shared style, colours and geometry for Figure 2 panels.

Every panel script imports this module. The rcParams below are the ones in effect
when the approved composite (Figure2.png, 7.2 x 6.6 in) was rendered: matplotlib
defaults + the project's common.py settings. Re-running the original make_fig2.py
with exactly these settings reproduces the approved PNG pixel-for-pixel.
Only 'svg.fonttype' (editable SVG text) and two determinism settings are added;
none of them changes the PDF/PNG appearance.

Panels are drawn on fixed-size canvases cut from the composite layout, with every
axes at exactly the size/position (inches) it has in the composite, so the files can
be placed in Illustrator at 100 % scale. Panel letters are NOT drawn.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
PANELS = os.path.join(HERE, "..", "panels")

mpl.rcdefaults()
RC = {  # from v4/common.py (unchanged)
    "font.family": "Arial", "font.size": 6.5, "axes.titlesize": 7, "axes.labelsize": 6.5,
    "xtick.labelsize": 5.5, "ytick.labelsize": 5.5,
    "axes.linewidth": 0.6, "xtick.major.width": 0.5, "ytick.major.width": 0.5,
    "xtick.major.size": 2, "ytick.major.size": 2,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42,
    # added for Illustrator hand-off / reproducibility (no visual effect)
    "ps.fonttype": 42, "svg.fonttype": "none", "svg.hashsalt": "Figure2",
}
mpl.rcParams.update(RC)

INK = "#222222"
MASK = "#dcdcdc"
COL = dict(ref="#1f3b73", cds="#1f3b73", utr="#8fb3e0",
           tss="#7b4fa0", tts="#2a9d8f", body="#4a78c2", inter="#b5b5b5", trunc="#9a9a9a")

DPI = 300


def data(name, **kw):
    """Read a tab-separated source-data file from ../data."""
    return pd.read_csv(os.path.join(DATA, name), sep="\t", **kw)


def new_panel(key):
    """Create the panel canvas and its axes.

    Returns (fig, axes_dict, to_fig) where axes_dict maps element name -> Axes placed at
    the exact composite size/position, and to_fig(x_in, y_in) converts composite
    coordinates (inches from the composite's lower-left corner) to canvas fractions.
    """
    lay = data("Figure2_panel_layout.tsv")
    lay = lay[lay.panel == key].sort_values("order")
    cv = lay[lay.element == "canvas"].iloc[0]
    W, H = cv.w_in, cv.h_in
    fig = plt.figure(figsize=(W, H))
    to_fig = lambda x, y: ((x - cv.x0_in) / W, (y - cv.y0_in) / H)
    axes = {}
    for r in lay[lay.kind == "axes"].itertuples():
        fx, fy = to_fig(r.x0_in, r.y0_in)
        axes[r.element] = fig.add_axes([fx, fy, r.w_in / W, r.h_in / H])
    return fig, axes, to_fig, lay


def anchor(lay, element):
    r = lay[lay.element == element].iloc[0]
    return r.x0_in, r.y0_in


def save(fig, name):
    """Write <name>.pdf/.svg/.png (300 dpi) into ../panels at the fixed canvas size."""
    os.makedirs(PANELS, exist_ok=True)
    fig.savefig(os.path.join(PANELS, f"{name}.pdf"), dpi=DPI, metadata={"CreationDate": None})
    fig.savefig(os.path.join(PANELS, f"{name}.svg"), dpi=DPI, metadata={"Date": None})
    fig.savefig(os.path.join(PANELS, f"{name}.png"), dpi=DPI)
    plt.close(fig)
