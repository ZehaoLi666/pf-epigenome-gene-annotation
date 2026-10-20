"""Figure 3D - cumulative distance from each predicted end to the nearest long-read end
cluster of the same gene (main or alternative). Genes without a call count as misses."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # find style.py next to this script
import numpy as np, pandas as pd
from style import Panel, COL, title, data

NE = pd.read_csv(data("Figure3_D_nearest_end.tsv"), sep="\t")

P = Panel('D')
axD = P.axes((0.54, 2.237985, 2.844735, 1.357221))
xg = np.arange(0, 1001, 5)
for end, c in [('TSS', COL['tss']), ('TTS', COL['tts'])]:
    d = NE[f'{end}_nearest_lr_end_bp'].values; d = np.where(np.isnan(d), np.inf, d)   # no call = miss
    cum = [(d <= v).sum() / len(d) for v in xg]; axD.plot(xg, cum, color=c, lw=1.2)
    f200 = (d <= 200).sum() / len(d); med = np.median(d[np.isfinite(d)])
    axD.text(560, 0.42 if end == 'TSS' else 0.58, f"{end}: {f200:.1%} ≤200 bp\nmedian {med:.0f} bp", fontsize=5.3, color=c, va='center')
axD.axvline(200, color='#999999', lw=0.5, ls=(0, (2, 2)))
axD.set_xlim(0, 1000); axD.set_ylim(0, 1); axD.set_yticks([0, 0.5, 1]); axD.set_yticklabels(['0', '50%', '100%'])
axD.set_xlabel("Distance to the nearest long-read end of the gene (bp)"); axD.set_ylabel(f"Reference genes\n(n = {len(NE):,})")
title(axD, "Predicted ends fall on long-read ends", y=1.02)
P.save()
