"""Figure 3F - accuracy of the full model (navy) and the model without chromatin (light blue)
by DAFT-seq expression quintile (same run as E). Bars: bootstrap 95% CI; * p < 0.05 (paired sign test)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # find style.py next to this script
import numpy as np, pandas as pd, matplotlib as mpl
from style import Panel, COL, data

TF3 = pd.read_csv(data("Figure3_F_chromatin_by_expression.tsv"), sep="\t")
QS = ['Q1 (lowest)', 'Q2', 'Q3', 'Q4', 'Q5 (highest)']
C_CH, C_NC = COL['ref'], COL['nochrom']

P = Panel('F')
axF2 = [P.axes((0.756, 0.495, 1.206367, 0.931784)), P.axes((2.178367, 0.495, 1.206367, 0.931784))]
for i, (ax, tl) in enumerate(zip(axF2, ['Long-read TSS', 'Long-read TTS'])):
    d = TF3[TF3.truth == tl].set_index('stratum').loc[QS]; x = np.arange(5)
    for m, cc, dx in (('A_full', C_CH, -0.08), ('L_nochrom', C_NC, 0.08)):
        ax.errorbar(x + dx, d[f'acc_{m}'], yerr=[d[f'acc_{m}'] - d[f'lo_{m}'], d[f'hi_{m}'] - d[f'acc_{m}']],
                    fmt='-o', ms=2.3, lw=0.9, elinewidth=0.6, color=cc, capsize=0)
    for xi, r in zip(x, d.itertuples()):
        if r.p_sign < 0.05:
            ax.text(xi, 0.99, '*', ha='center', va='top', fontsize=7, transform=mpl.transforms.blended_transform_factory(ax.transData, ax.transAxes))
    ax.set_ylim(0.5, 0.95); ax.set_xlim(-0.5, 4.5); ax.set_xticks(x)
    ax.set_xticklabels([f'Q{k + 1}\n{v:.0f}' if v >= 10 else f'Q{k + 1}\n{v:.1f}' for k, v in enumerate(d.daft_median)], fontsize=4.8, linespacing=1.0)
    ax.text(0.03, 0.03, tl, transform=ax.transAxes, fontsize=5.4, color=COL['tss'] if 'TSS' in tl else COL['tts'], va='bottom')
    if i: ax.tick_params(labelleft=False)
axF2[0].set_ylabel('Accuracy within 200 bp')
axF2[1].text(4.4, 0.62, 'full model', color=C_CH, fontsize=5.0, ha='right'); axF2[1].text(4.4, 0.585, '\u2212 chromatin', color=C_NC, fontsize=5.0, ha='right')
P.text(2.070367, 0.135, 'DAFT-seq expression quintile (median CPM)', ha='center', fontsize=5.6)
P.text(0.288, 1.658717, 'Chromatin gain by expression', fontsize=7, ha='left', va='bottom')
P.save()
