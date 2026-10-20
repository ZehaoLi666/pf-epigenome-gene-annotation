"""Figure 4F - accuracy within 200 bp by DAFT-seq expression quintile for reference genes (long-read) and model-predicted genes,
TSS (5'-end consensus truth) and TTS (major Siegel 2014 poly(A) site truth). Open circles: quintiles with <30 genes; error bars
(model only): bootstrap 95% CI; small grey numbers: genes per quintile; median CPM under each quintile.
The line key above the four plots and the explanatory note below them are part of this panel.
Data: ../data/Figure4_F_accuracy_by_expression.tsv   Output: ../panels/Figure4_F.pdf|.svg|.png"""
import numpy as np, pandas as pd, matplotlib as mpl
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))  # find style.py
from style import Panel, DATA, REF, V71, EXPT

DSF = pd.read_csv(DATA / 'Figure4_F_accuracy_by_expression.tsv', sep='\t')
QS = ['Q1 (lowest)', 'Q2', 'Q3', 'Q4', 'Q5 (highest)']
FM = {'Model (out-of-fold)': ('Our model', REF, '-', 1.3), 'Final annotation (reference / model)': ('Long-read reference', EXPT, '-', 1.0),
      'PlasmoDB v71': ('PlasmoDB v71', V71, '-', 1.0), 'Fixed UTR length': ('Fixed UTR length', '#9a9a9a', (0, (2.5, 1.5)), 0.9)}
FP = [('Reference genes (long-read)', 'TSS consensus', 'Reference genes: TSS'), ('Reference genes (long-read)', 'Siegel 2014 poly(A)', 'Reference genes: TTS'),
      ('Model-predicted genes', 'TSS consensus', 'Model-predicted genes: TSS'), ('Model-predicted genes', 'Siegel 2014 poly(A)', 'Model-predicted genes: TTS')]

pn = Panel('F'); fig = pn.fig
axF = []; FX = [0.085, 0.32, 0.555, 0.79]
for i, (grp, tru, ttl) in enumerate(FP):
    ax = pn.add_axes([FX[i], 0.062, 0.185, 0.115]); axF.append(ax); x = np.arange(5)
    for m, (lab, cc, ls, lw) in FM.items():
        if grp.startswith('Model') and m.startswith('Final'): continue
        d = DSF[(DSF.group == grp) & (DSF.truth == tru) & (DSF.method == m)].set_index('stratum').loc[QS]
        ax.plot(x, d.acc200, color=cc, ls=ls, lw=lw, zorder=4 if m.startswith('Model') else 3)
        if m.startswith('Model'): ax.errorbar(x, d.acc200, yerr=[d.acc200 - d.lo, d.hi - d.acc200], fmt='none', ecolor=cc, elinewidth=0.5, capsize=0, zorder=2)
        for xi, r in zip(x, d.itertuples()): ax.plot(xi, r.acc200, 'o', ms=2.6, mfc=cc if r.n >= 30 else 'white', mec=cc, mew=0.7, zorder=5 if m.startswith('Model') else 3)
    nn = DSF[(DSF.group == grp) & (DSF.truth == tru) & (DSF.method == 'Model (out-of-fold)')].set_index('stratum').loc[QS]
    for xi, r in zip(x, nn.itertuples()): ax.text(xi, 0.215, f'{r.n:,}', ha='center', va='bottom', fontsize=4.3, color='#777777')
    ax.set_ylim(0.2, 1.0); ax.set_xlim(-0.4, 4.4); ax.set_xticks(x)
    ax.set_xticklabels([f'Q{k + 1}\n{v:.0f}' if v >= 10 else f'Q{k + 1}\n{v:.1f}' for k, v in enumerate(nn.DAFT_median)], fontsize=4.9, linespacing=1.0)
    ax.text(0.0, 1.03, ttl, transform=ax.transAxes, fontsize=5.6, va='bottom', ha='left')
    if i: ax.tick_params(labelleft=False)
    if i == 0: ax.set_ylabel('Accuracy within 200 bp')
fig.text(*pn.xy(0.5, 0.008), 'Expression quintile (DAFT-seq; median CPM below each quintile). TSS truth: 5′-end consensus; TTS truth: major Siegel 2014 poly(A) site. Small numbers, genes per quintile',
         ha='center', fontsize=5.3, color='#444444')
hF = [mpl.lines.Line2D([], [], color=cc, ls=ls, lw=lw + 0.2, label=lab) for m, (lab, cc, ls, lw) in FM.items()] + \
     [mpl.lines.Line2D([], [], marker='o', ls='', ms=2.6, mfc='white', mec='#555555', label='n < 30')]
fig.legend(handles=hF, frameon=False, fontsize=5.3, loc='lower left', ncol=5, handlelength=1.6, columnspacing=1.2, bbox_to_anchor=pn.xy(FX[0] - 0.005, 0.2), borderaxespad=0)
pn.save()
