"""Figure 3 shared symbol key (sits below panel C in the composite). No source data."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # find style.py next to this script
from style import Panel, COL

P = Panel('legend_key')
lgx = P.axes((0.72, 3.863974, 6.336, 0.126)); lgx.axis('off'); lgx.set_xlim(0, 100); lgx.set_ylim(0, 1)
xx = 0.5
lgx.plot([xx, xx + 2], [0.5, 0.5], color=COL['model_full'], lw=2); lgx.plot([xx + 2.6, xx + 4.6], [0.5, 0.5], color=COL['model_trunc'], lw=2)
lgx.text(xx + 5.2, 0.5, "complete-CDS / truncated model", fontsize=4.9, va='center'); xx += 23
lgx.plot(xx, 0.5, 'o', ms=2.6, color=COL['tss']); lgx.plot(xx + 1.2, 0.5, 'o', ms=2.6, color=COL['tts'])
lgx.text(xx + 2.2, 0.5, "model end within 200 bp of a call", fontsize=4.9, va='center'); xx += 21
lgx.plot([xx, xx], [0.2, 0.8], color=COL['tss'], lw=1.2); lgx.plot([xx + 0.8, xx + 0.8], [0.2, 0.8], color=COL['tts'], lw=1.2)
lgx.text(xx + 1.8, 0.5, "long-read end clusters", fontsize=4.9, va='center'); xx += 15
lgx.plot(xx, 0.5, marker='v', ms=4, color=COL['tss'], ls='none'); lgx.text(xx + 1.2, 0.5, "main call", fontsize=4.9, va='center'); xx += 8
lgx.plot(xx, 0.5, marker='v', ms=3.6, mfc='white', mec=COL['tss'], ls='none'); lgx.plot(xx + 1.2, 0.5, marker='v', ms=3.6, mfc='white', mec=COL['grey_call'], ls='none')
lgx.text(xx + 2.4, 0.5, "other calls (grey: no long-read end)", fontsize=4.9, va='center'); xx += 22
lgx.plot([xx, xx + 1.6], [0.5, 0.5], color=COL['tss'], lw=0.6, ls=(0, (2, 2))); lgx.text(xx + 2.2, 0.5, "main long-read TSS", fontsize=4.9, va='center')
P.save("Figure3_legend_key")
