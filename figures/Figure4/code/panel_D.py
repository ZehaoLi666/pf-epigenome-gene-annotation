"""Figure 4D - change of each refined UTR end relative to PlasmoDB v71, by source (long reads / model).
agrees: |difference| <=100 bp; trimmed: refined UTR >100 bp shorter; extended: >100 bp longer; UTR added: no UTR in v71.
Numbers in bars are percentages (segments <6 % unlabelled). The colour key above the bars is part of this panel.
Data: ../data/Figure4_D_change_vs_v71.tsv   Output: ../panels/Figure4_D.pdf|.svg|.png"""
import pandas as pd, numpy as np, matplotlib as mpl
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))  # find style.py
from style import Panel, DATA, COL, SM

TD = pd.read_csv(DATA / 'Figure4_D_change_vs_v71.tsv', sep='\t')
_d = TD.set_index(['end', 'source'])
CHG = ['agrees', 'trimmed', 'extended', 'UTR added']
CHC = {'agrees': COL['agree'], 'trimmed': COL['ol_lnc'], 'extended': COL['short'], 'UTR added': COL['noutr']}
CHL = {'agrees': 'agrees (≤100 bp)', 'trimmed': 'trimmed (v71 too long)', 'extended': 'extended (v71 too short)', 'UTR added': 'UTR added (none in v71)'}

pn = Panel('D')
axD = pn.add_axes([0.15, 0.325, 0.33, 0.125])
rows = [('5′ end', 'long-read'), ('5′ end', 'model'), ('3′ end', 'long-read'), ('3′ end', 'model')]; yD = np.array([3.3, 2.3, 1.0, 0.0])
for y, (e, s) in zip(yD, rows):
    r = _d.loc[(e, s)]; left = 0
    for c in CHG:
        axD.barh(y, r[c], left=left, color=CHC[c], height=0.72, lw=0)
        if r[c] >= 0.06: axD.text(left + r[c] / 2, y, f'{r[c] * 100:.0f}', ha='center', va='center', fontsize=SM - 0.3, color='white')
        left += r[c]
axD.set_yticks(yD); axD.set_yticklabels([f'{e[:2]} {"long reads" if s == "long-read" else "model"} ({int(_d.loc[(e, s), "n"]):,})' for e, s in rows], fontsize=SM + 0.2)
axD.set_xlim(0, 1); axD.set_xticks([0, 0.5, 1]); axD.set_xticklabels(['0', '50', '100'])
axD.set_xlabel(f'Refined ends (%); no refined end: 5′ {int(_d.loc[("5′ end", "none"), "n"]):,}, 3′ {int(_d.loc[("3′ end", "none"), "n"]):,} genes', fontsize=6)
axD.set_ylim(-0.5, 3.8); axD.spines['left'].set_visible(False); axD.tick_params(axis='y', length=0)
hD = [mpl.patches.Patch(color=CHC[c], label=CHL[c]) for c in CHG]
axD.legend(handles=hD, frameon=False, fontsize=SM - 0.2, loc='lower left', bbox_to_anchor=(-0.32, 1.0), ncol=2, handlelength=0.9, columnspacing=0.8, borderaxespad=0.1)
pn.save()
