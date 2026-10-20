"""Figure 1H - v71 5′ UTR length vs farthest reference TSS
Standalone: python panel_H.py   (reads ../data/Figure1_H_5UTR_lengths.tsv; writes ../panels/Figure1_H.pdf/.svg/.png)
Drawing code copied from the approved make_fig1.py; panel letter is NOT drawn (added in Illustrator)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

T = read("Figure1_H_5UTR_lengths.tsv")
fig, (ax,) = panel_figure('H')
k = 0; lbl = "5′ UTR"
Lv = T.v71_UTR_length_bp; Dm = T.farthest_reference_end_bp; cc = T['class']
gcol = {'agree': COL['agree'], 'over-long': '#d9822b', 'too short': COL['short'], 'no annotated UTR': COL['noutr'], 'between ends': COL['between']}
tr = lambda x: np.log10(1 + np.clip(np.asarray(x, float), 0, None)); ticks = [0, 100, 300, 1000, 3000]
for gname in ['agree', 'between ends', 'no annotated UTR', 'too short', 'over-long']:
    m = (cc == gname) & Dm.notna()
    jit = np.where(Lv[m] == 0, np.random.default_rng(1).uniform(-0.06, 0.06, m.sum()), 0)   # x-jitter for genes without an annotated UTR
    ax.scatter(tr(Lv[m]) + jit, tr(Dm[m]), s=1.3, lw=0, color=gcol[gname], alpha=0.7, rasterized=True)
lim = (-0.15, np.log10(6001)); ax.plot(lim, lim, color='#777777', lw=0.5, ls=(0, (2, 2)), zorder=0)
ax.set_xlim(lim); ax.set_ylim(lim); tl = ['0', '100', '300', '1k', '3k']
ax.set_xticks(tr(ticks)); ax.set_xticklabels(tl); ax.set_yticks(tr(ticks)); ax.set_yticklabels(tl if k == 0 else [])
ax.set_xlabel(f"PlasmoDB v71 {lbl} length (bp)", fontsize=5.8)
if k == 0: ax.set_ylabel("Farthest reference end (bp)", fontsize=5.8, labelpad=1)
title(ax, lbl, y=1.02)
if k == 1:
    for i, g_ in enumerate(['agree', 'over-long', 'too short', 'no annotated UTR', 'between ends']):
        ax.text(0.98, 0.40 - i * 0.08, g_, transform=ax.transAxes, ha='right', fontsize=5.0, color=gcol[g_])
save_panel(fig, 'H')
