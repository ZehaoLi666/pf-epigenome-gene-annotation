import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

"""Shared line key for panels A-D (place below row A-D in Illustrator)."""
from matplotlib.lines import Line2D
n = data("Figure5_A_atac_n.tsv", index_col=0)["n_genes"]
fig = plt.figure(figsize=(5.0, 0.3))
h = [Line2D([], [], color=ACOL[a], lw=1.1 if a != "experimental" else 0.9, ls="-" if a != "experimental" else "--",
            label=f"{ALAB[a]} (n={n[a]})") for a in ["refined", "v71", "experimental"]]
fig.legend(handles=h, loc="center", ncol=3, frameon=False, fontsize=6)
save(fig, "Figure5_legend_AD")
