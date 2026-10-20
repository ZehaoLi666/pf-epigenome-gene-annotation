"""Figure 1C - Completeness of stage-confirmed nuclear mRNA models
Standalone: python panel_C.py   (reads ../data/Figure1_C_model_completeness.tsv; writes ../panels/Figure1_C.pdf/.svg/.png)
Drawing code copied from the approved make_fig1.py; panel letter is NOT drawn (added in Illustrator)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

CM = read("Figure1_C_model_completeness.tsv").set_index('cls'); tot = int(CM.n_models.sum())
fig, (axC,) = panel_figure('C')
rowsC = [('full CDS', 'Complete CDS', COL['ref']), ('5′ truncated', '5′ end inside CDS', COL['trunc']), ('3′ truncated', '3′ end inside CDS', '#c8c8c8')]
yc_ = np.arange(3)[::-1]
for yi, (k, l, cc) in zip(yc_, rowsC):
    n = int(CM.loc[k, 'n_models']); axC.barh(yi, n / tot * 100, color=cc, height=0.6)
    axC.text(n / tot * 100 + 1, yi, f"{n/tot:.0%}", va='center', fontsize=5.3)
axC.set_yticks(yc_); axC.set_yticklabels([r[1] for r in rowsC], fontsize=5.1); axC.tick_params(axis='y', length=0)
axC.set_xlim(0, 55); axC.set_xticks([0, 25, 50]); axC.set_xlabel(f"Stage-confirmed mRNA models (%, n = {tot:,})", fontsize=5.6, labelpad=1)
title(axC, "Most models are truncated", y=1.02)
save_panel(fig, 'C')
