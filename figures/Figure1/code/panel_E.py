"""Figure 1E - Internal-priming signature at 3′ ends
Standalone: python panel_E.py   (reads ../data/Figure1_E_internal_priming.tsv; writes ../panels/Figure1_E.pdf/.svg/.png)
Drawing code copied from the approved make_fig1.py; panel letter is NOT drawn (added in Illustrator)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

IP = read("Figure1_E_internal_priming.tsv")
fig, (axE,) = panel_figure('E')
ye = np.arange(len(IP))[::-1]
axE.barh(ye, IP.frac * 100, color=[COL['art'], COL['art'], '#e8b27a', COL['utr'], '#cfcfcf'], height=0.62)
for yi, v in zip(ye, IP.frac): axE.text(v * 100 + 0.6, yi, f"{v:.0%}", va='center', fontsize=5.3)
axE.set_yticks(ye); axE.set_yticklabels(list(IP.label), fontsize=5.0); axE.tick_params(axis='y', length=0)
axE.set_xlim(0, 34); axE.set_xticks([0, 10, 20, 30]); axE.set_xlabel("3′ ends with ≥80% A in the next 20 nt (%)", fontsize=5.6, labelpad=1)
title(axE, "False 3′ ends from internal priming", y=1.02)
save_panel(fig, 'E')
