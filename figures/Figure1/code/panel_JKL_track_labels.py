"""Figure 1 J-L - shared track row labels (sit left of panel J in the composite).
Standalone: python panel_JKL_track_labels.py  (no data; writes ../panels/Figure1_JKL_track_labels.pdf/.svg/.png)
The canvas has the same top edge and height as Figure1_J/K/L, so top-aligning it with them restores the composite."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from track_panel import *
fig, axs = panel_figure('JKL_track_labels')   # the five (invisible) track axes of panel J
for ax, l in zip(axs, ROWLAB):
    ax.set_axis_off()
    ax.text(-0.04, 0.5, l, transform=ax.transAxes, fontsize=5.4, ha='right', va='center', linespacing=0.95)
save_panel(fig, 'JKL_track_labels')
