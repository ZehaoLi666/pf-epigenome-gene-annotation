"""Figure 1D - Non-templated 5′ G (cap signature) at long-read TSS clusters
Standalone: python panel_D.py   (reads ../data/Figure1_D_tss_capG.tsv; writes ../panels/Figure1_D.pdf/.svg/.png)
Drawing code copied from the approved make_fig1.py; panel letter is NOT drawn (added in Illustrator)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

TG = read("Figure1_D_tss_capG.tsv")
fig, (axD,) = panel_figure('D')
xx = np.arange(4); axD.bar(xx, TG.ONT_extraG_frac, color=[COL['tss'], COL['tss'], COL['art'], COL['art']], width=0.62)
for xi, v in enumerate(TG.ONT_extraG_frac): axD.text(xi, v + 0.015, f"{v:.0%}", ha='center', fontsize=5.3)
axD.set_xticks(xx); axD.set_xticklabels(['Upstream,\nCAGE high', 'Upstream,\nCAGE low', 'Inside\nCDS', 'Inside CDS,\nweak'], fontsize=4.8)   # = TG.label, line-broken
axD.set_ylim(0, 0.68); axD.set_yticks([0, 0.2, 0.4, 0.6]); axD.set_yticklabels(['0', '20%', '40%', '60%']); axD.set_ylabel("ONT reads with a\nnon-templated 5′ G", fontsize=5.6)
title(axD, "Degraded 5′ ends lack the cap", y=1.02)
save_panel(fig, 'D')
