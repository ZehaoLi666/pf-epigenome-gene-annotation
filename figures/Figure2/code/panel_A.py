"""Figure 2A - label-definition schematic (no source data; geometry is illustrative).
Drawing code copied verbatim from the approved make_fig2.py; only the axes placement
(and removal of the panel letter) differs. Output: ../panels/Figure2_A.{pdf,svg,png}
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # find style.py from any cwd
from style import *
from matplotlib.patches import Rectangle

fig, AX, to_fig, lay = new_panel("A")
axA = AX["schematic"]
axA.set_xlim(0, 100); axA.set_ylim(0, 30); axA.axis('off')
axA.text(3, 30.5, "Labels come only from the long-read reference (per strand)", fontsize=7, va='bottom')
X = lambda v: 3 + v * 0.95
yg = 21.5; axA.plot([X(0), X(100)], [yg, yg], color='#bbbbbb', lw=0.5)
genes = [(6, 36, 11, 32), (58, 92, 64, 88)]
for a, b, c0, c1 in genes:
    axA.add_patch(Rectangle((X(a), yg - 0.8), X(c0) - X(a), 1.6, fc=COL['utr'], ec='none'))
    axA.add_patch(Rectangle((X(c0), yg - 1.3), X(c1) - X(c0), 2.6, fc=COL['ref'], ec='none'))
    axA.add_patch(Rectangle((X(c1), yg - 0.8), X(b) - X(c1), 1.6, fc=COL['utr'], ec='none'))
# reads: full-CDS models and a truncated model
for k, (a, b) in enumerate([(6.3, 36), (7.6, 35.6), (58.3, 92), (59.5, 91.6)]):
    axA.plot([X(a), X(b)], [yg + 2.7 + (k % 2) * 1.0] * 2, color='#555555', lw=0.9)
axA.plot([X(74), X(91.8)], [yg + 4.7] * 2, color='#9a9a9a', lw=0.9)
axA.text(X(47), yg + 3.2, "complete-CDS models", fontsize=5.0, ha='center', color='#555555')
axA.text(X(83), yg + 5.6, "truncated model", fontsize=4.8, ha='center', color='#888888')
axA.annotate("", xy=(X(36.4), yg - 2.4), xytext=(X(57.6), yg - 2.4), arrowprops=dict(arrowstyle='<->', lw=0.6, color=INK))
axA.text(X(47), yg - 3.0, "no read spans the gap", fontsize=4.9, ha='center', va='top', style='italic')
# label track
ys = 10.5; CL = {'TSS': COL['tss'], 'body': COL['body'], 'TTS': COL['tts'], 'inter': COL['inter'], 'mask': MASK}
segs = [(0, 2, 'inter'), (2, 4.5, 'mask'), (4.5, 7.5, 'TSS'), (7.5, 34.5, 'body'), (34.5, 37.5, 'TTS'), (37.5, 40, 'mask'),
        (40, 54, 'inter'), (54, 56.5, 'mask'), (56.5, 59.5, 'TSS'), (59.5, 72.5, 'body'), (72.5, 75.5, 'TSS'),
        (75.5, 90.5, 'body'), (90.5, 93.5, 'mask'), (93.5, 100, 'inter')]
for a, b, k in segs:
    axA.add_patch(Rectangle((X(a), ys - 1.1), X(b) - X(a), 2.2, fc=CL[k], ec='white', lw=0.3, hatch='/////' if k == 'mask' else None))
for x, l, k in [(6, "TSS", 'TSS'), (21, "gene body", 'body'), (36, "TTS", 'TTS'), (47, "flank (intergenic)", 'inter'),
                (74, "internal TSS\n(cap-G only, w 0.5)", 'TSS'), (92, "masked", 'mask')]:
    axA.text(X(x), ys - 2.0, l, fontsize=5.0, ha='center', va='top', color=CL[k] if k not in ('inter', 'mask') else '#666666', linespacing=0.95)
axA.text(0, ys, "labels", fontsize=5.4, ha='left', va='center', color='#444444')
notes = ["TSS: 5\u2032 end clusters \u226450 bp; weight = read share \u00d7 (1 with cap-G, 0.5 without)",
         "TTS: masked when \u226560% A follows the end (internal priming)",
         "Flank: unspanned same-strand gap, tier A/B/C by end evidence (w 1 / 0.6 / 0.3); spanned gaps masked"]
for i, t in enumerate(notes):
    axA.text(X(0), 3.6 - i * 2.2, "\u2022 " + t, fontsize=5.0, va='center', color='#333333')
save(fig, "Figure2_A")
