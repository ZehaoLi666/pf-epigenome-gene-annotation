"""Shared style and panel geometry for Figure 5.

Every panel script imports this module. Panels are drawn with the same physical
axes size (inches) as in the assembled draft figure (7.2 x 7.6 in), so they can be
placed in Illustrator at 100 % scale. Panel letters are NOT drawn.
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

RC = {
    "axes.labelsize": 8.0, "axes.linewidth": 0.6, "axes.spines.right": False,
    "axes.spines.top": False, "axes.titlelocation": "left", "axes.titlesize": 8.0,
    "figure.dpi": 200.0, "font.size": 8.0, "legend.fontsize": 7.0,
    "legend.frameon": False, "lines.linewidth": 1.2, "patch.linewidth": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "savefig.bbox": "tight", "savefig.dpi": 300.0,
    "xtick.labelsize": 6.0, "xtick.major.size": 3.0, "xtick.major.width": 0.6,
    "ytick.labelsize": 6.0, "ytick.major.size": 3.0, "ytick.major.width": 0.6,
}
mpl.rcParams.update(RC)

NAVY = "#1f3a68"   # refined annotation
ORG = "#e08214"    # PlasmoDB v71
GRY = "#4d4d4d"    # experimental 5' ends
LT = "#7fa7d9"     # light blue (secondary series)
LG = "#bdbdbd"     # light grey (controls / null)
ACOL = {"refined": NAVY, "v71": ORG, "experimental": GRY}
ALAB = {"refined": "Refined", "v71": "PlasmoDB v71", "experimental": "Experimental 5\u2032 end"}


def data(name, **kw):
    """Read a tab-separated source-data file from ../data."""
    return pd.read_csv(os.path.join(DATA, name), sep="\t", **kw)


def smooth(v, w):
    """Centred moving average (same length), as used for all metagene curves."""
    return np.convolve(v, np.ones(w) / w, mode="same")


def new_panel(key, margin=0.9):
    """Create a figure whose axes have exactly the size/arrangement they have in the
    assembled figure. Returns (fig, [axes...]) in layout order."""
    lay = data("Figure5_panel_layout.tsv")
    lay = lay[lay.panel == key].sort_values("order")
    x0, y0 = lay.x0_in.min(), lay.y0_in.min()
    x1, y1 = (lay.x0_in + lay.w_in).max(), (lay.y0_in + lay.h_in).max()
    W, H = (x1 - x0) + 2 * margin, (y1 - y0) + 2 * margin
    fig = plt.figure(figsize=(W, H))
    axes = [fig.add_axes([(r.x0_in - x0 + margin) / W, (r.y0_in - y0 + margin) / H,
                          r.w_in / W, r.h_in / H]) for r in lay.itertuples()]
    return fig, axes


def save(fig, name):
    os.makedirs(PANELS, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(PANELS, f"{name}.{ext}"), dpi=300, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
