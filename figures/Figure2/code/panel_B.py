"""Figure 2B - number of training labels per class.
Reads ../data/Figure2_B_label_counts.tsv (derived by prep_data.py from the per-label
tables, also supplied as Figure2_B_source_*.tsv). Output: ../panels/Figure2_B.{pdf,svg,png}
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # find style.py from any cwd
from style import *

fig, AX, to_fig, lay = new_panel("B")
axB = AX["counts"]
T = data("Figure2_B_label_counts.tsv").sort_values(["bar_order", "segment_order"])
SEGCOL = {("TSS", "terminal, cap-G"): COL['tss'], ("TSS", "terminal, no cap-G"): '#b49bcc',
          ("TSS", "truncated-model end"): '#d6c6e6', ("TSS", "inside CDS (w 0.5)"): '#ece4f3',
          ("TTS", "terminal, poly(A) reads"): COL['tts'], ("TTS", "terminal, no poly(A) reads"): '#7fc4bb',
          ("TTS", "truncated-model end"): '#bfe3de', ("TTS", "inside CDS"): '#e3f3f1',
          ("Flank", "tier A"): '#8a8a8a', ("Flank", "tier B"): '#b5b5b5', ("Flank", "tier C"): '#d9d9d9',
          ("Masked gaps", "spanned by reads"): MASK}
bars = []  # [(bar name, [(segment, n, colour), ...]), ...] top to bottom
for _, g in T.groupby("bar_order", sort=True):
    nm = g.bar.iloc[0]
    bars.append((nm, [(r.segment, int(r.n), SEGCOL[(nm, r.segment)]) for r in g.itertuples()]))
yb = np.arange(len(bars))[::-1]
for yi, (nm, parts) in zip(yb, bars):
    left = 0
    for l, n, c in parts:
        axB.barh(yi, n, left=left, color=c, height=0.6, ec='white', lw=0.4, hatch='/////' if c == MASK else None)
        if n >= 700:
            axB.text(left + n / 2, yi, f"{n:,}", ha='center', va='center', fontsize=4.8,
                     color='white' if c in (COL['tss'], COL['tts'], '#8a8a8a') else INK)
        left += n
    axB.text(left + 80, yi, f"{left:,}", va='center', fontsize=5.3, fontweight='bold')
axB.set_yticks(yb); axB.set_yticklabels([b[0] for b in bars], fontsize=5.6); axB.tick_params(axis='y', length=0)
axB.set_xlim(0, 8400); axB.set_xticks([0, 2500, 5000]); axB.set_xticklabels(['0', '2,500', '5,000']); axB.spines['bottom'].set_bounds(0, 5000)
axB.set_xlabel("Labelled ends or gaps", fontsize=5.8, labelpad=1)
axB.text(0.98, 0.03, "Gene bodies:\n2,389 reference genes", transform=axB.transAxes, fontsize=5.2, color=COL['body'], ha='right', va='bottom')
save(fig, "Figure2_B")
