"""Figure 3G - where predictions fall at +/-200 bp: main long-read end, alternative long-read
end of the same gene, longer / shorter UTR (>200 bp from every long-read end), or no call."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # find style.py next to this script
import pandas as pd
from matplotlib.patches import Rectangle
from style import Panel, COL, INK, title, data

G = pd.read_csv(data("Figure3_G_outcome_by_gene.tsv"), sep="\t")
cats = [('main end', 'main', COL['ref']), ('alternative end', 'alternative', COL['alt']), ('longer', 'longer', COL['longer']),
        ('shorter', 'shorter', COL['shorter']), ('no call', 'no call', COL['nocall'])]

P = Panel('G')
axE = P.axes((4.031265, 0.495, 2.844735, 1.096217))
for yi, end in enumerate(['TSS', 'TTS']):
    left = 0; v = G['cat_' + end].value_counts()
    for k, l, c in cats:
        n = int(v.get(k, 0)); f = n / len(G); axE.barh(yi, f, left=left, color=c, height=0.6, ec='white', lw=0.4)
        if f >= 0.06:
            axE.text(left + f / 2, yi, f"{f:.0%}", ha='center', va='center', fontsize=5.2, color='white' if c in (COL['ref'], COL['alt'], COL['longer']) else INK)
        left += f
axE.set_yticks([0, 1]); axE.set_yticklabels(['TSS', 'TTS']); axE.invert_yaxis(); axE.set_xlim(0, 1); axE.set_xticks([0, 0.5, 1]); axE.set_xticklabels(['0', '50%', '100%'])
axE.set_xlabel("Reference genes")
for i, (k, l, c) in enumerate(cats):
    xx_ = [0.0, 0.17, 0.40, 0.58, 0.80][i]; yy_ = 0.5
    axE.add_patch(Rectangle((xx_, yy_ - 0.03), 0.05, 0.06, transform=axE.transAxes, fc=c, ec='none', clip_on=False))
    axE.text(xx_ + 0.07, yy_, l, transform=axE.transAxes, fontsize=4.8, va='center')
title(axE, "Where predictions fall (±200 bp)", y=1.02)
P.save()
