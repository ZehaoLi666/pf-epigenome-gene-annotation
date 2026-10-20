"""Figure 3A - boundary model schematic (no source data; drawn from constants)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # find style.py next to this script
import numpy as np
from matplotlib.patches import Rectangle
from style import Panel, COL, box, arrow

P = Panel('A')
axA = P.axes((0.54, 7.603185, 4.008490, 0.991815))      # composite inches
axA.set_xlim(0, 100); axA.set_ylim(0, 40); axA.axis('off')
axA.text(2, 41, "Boundary model (applied to each strand)", fontsize=7, va='bottom')
box(axA, 1, 10, 21, 24, ec='#888888'); axA.text(11.5, 30.5, "Input", fontsize=6.2, fontweight='bold', ha='center')
axA.text(11.5, 21.5, "90 channels\nchromatin · RNA\nDNA · gLM\n(Figure 2)", fontsize=5.2, ha='center', va='center', linespacing=1.2, color='#444444')
arrow(axA, (22.5, 22), (27.5, 22))
box(axA, 28, 10, 18, 24, fc='#eef2f9', ec=COL['ref']); axA.text(37, 30.5, "Dilated CNN", fontsize=6.2, fontweight='bold', ha='center', color=COL['ref'])
for k in range(4):
    axA.add_patch(Rectangle((30.5 + k * 3.4, 14 + k * 1.6), 2.6, 12 - k * 3.2, fc=COL['ref'], alpha=0.25 + 0.15 * k, ec='none'))
axA.text(37, 11.8, "50-bp output bins", fontsize=4.9, ha='center', color='#444444')
arrow(axA, (46.5, 22), (51.5, 22))
x0 = 52; xs = np.linspace(0, 1, 200)
curves = [("TSS", COL['tss'], lambda x: np.exp(-((x - 0.15) / 0.03) ** 2) + 0.35 * np.exp(-((x - 0.3) / 0.03) ** 2)),
          ("gene body", COL['body'], lambda x: 1 / (1 + np.exp(-(x - 0.18) * 60)) * (1 - 1 / (1 + np.exp(-(x - 0.85) * 60)))),
          ("TTS", COL['tts'], lambda x: np.exp(-((x - 0.86) / 0.03) ** 2))]
for j, (nm, c, f) in enumerate(curves):
    yb = 27.5 - j * 8.0
    axA.plot(x0 + xs * 20, yb + f(xs) * 5.2, color=c, lw=0.9)
    axA.plot([x0, x0 + 20], [yb, yb], color='#cccccc', lw=0.4)
    axA.text(x0 + 20.8, yb + 1.5, nm, fontsize=5.2, color=c, va='center')
axA.text(x0 + 10, 6.2, "sigmoid outputs; flank = all 0", fontsize=4.9, ha='center', color='#444444')
arrow(axA, (80.5, 22), (84.0, 22))
box(axA, 84.5, 7, 15, 30, ec='#888888'); axA.text(92, 33, "Gene calls", fontsize=6.2, fontweight='bold', ha='center')
axA.text(92, 20.5, "peaks: one per\n300 bp\n\nTSS ≤3 kb\nupstream of ATG\nTTS ≤2 kb\ndownstream of stop", fontsize=4.8, ha='center', va='center', linespacing=1.12, color='#444444')
P.save()
