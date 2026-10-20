"""Figure 4A - TSS accuracy against the consensus of three 5′-end datasets (1,649 genes).
Dot = fraction of genes whose end lies within 200 bp of the truth position; bar = bootstrap 95% CI.
Data: ../data/Figure4_A_tss_consensus_acc200.tsv   Output: ../panels/Figure4_A.pdf|.svg|.png"""
import pandas as pd, numpy as np
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))  # find style.py
from style import Panel, DATA, COLM, SM

D = pd.read_csv(DATA / 'Figure4_A_tss_consensus_acc200.tsv', sep='\t')
n = int(D.n.iloc[0])
pn = Panel('A')
ax = pn.add_axes([0.2, 0.80, 0.26, 0.16])          # composite figure-fraction rect (7.2 x 8.9 in canvas)
y = np.arange(len(D))[::-1]
for yy in y: ax.axhline(yy, color='#f0f0f0', lw=0.5, zorder=0)
for yy, r in zip(y, D.itertuples()):
    c = COLM[r.key]; ax.plot([r.ci_lo, r.ci_hi], [yy, yy], color=c, lw=0.9, solid_capstyle='butt'); ax.plot(r.acc200, yy, 'o', ms=4, color=c, zorder=3)
    ax.text(r.ci_hi + 0.025, yy, f'{r.acc200:.2f}', va='center', fontsize=SM, color=c, fontweight='bold' if r.key == 'Ours_chromatin' else 'normal')
ax.set_yticks(y); ax.set_yticklabels(D.method); ax.set_xlim(0, 1); ax.set_xticks([0, 0.5, 1]); ax.set_ylim(-0.6, len(D) - 0.4)
for tl, k in zip(ax.get_yticklabels(), D.key): tl.set_color(COLM[k])
ax.set_xlabel(f'Genes with TSS ≤200 bp from the consensus\nof three 5′-end datasets (n = {n:,})', linespacing=1.05)
pn.save()
