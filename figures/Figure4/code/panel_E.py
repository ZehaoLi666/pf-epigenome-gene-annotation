"""Figure 4E - 5'UTR and 3'UTR length distributions, refined annotation (navy, solid) vs PlasmoDB v71 (orange, dashed).
Gaussian KDE of log10 length (bandwidth factor 0.18); dotted lines = medians.
Data: ../data/Figure4_E_utr_lengths_per_gene.tsv (per-gene lengths), Figure4_E_utr_length.tsv (medians / IQR / Mann-Whitney p)
Output: ../panels/Figure4_E.pdf|.svg|.png"""
import numpy as np, pandas as pd
from scipy import stats
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))  # find style.py
from style import Panel, DATA, REF, V71, SM

LG = pd.read_csv(DATA / 'Figure4_E_utr_lengths_per_gene.tsv', sep='\t')
TE = pd.read_csv(DATA / 'Figure4_E_utr_length.tsv', sep='\t')
ALAB = {'Refined': 'Refined annotation', 'v71': 'PlasmoDB v71'}
UL = {(u, k): LG[(LG.utr == u) & (LG.set == ALAB[k])].length_bp for u in ('5′', '3′') for k in ALAB}

pn = Panel('E'); axE = []
for i, u in enumerate(('5′', '3′')):
    ax = pn.add_axes([0.62 + i * 0.195, 0.325, 0.165, 0.125]); axE.append(ax); xs = np.linspace(1, 4, 300)
    for k, c, ls in (('v71', V71, '--'), ('Refined', REF, '-')):
        v = np.log10(UL[(u, k)][UL[(u, k)] > 0]); kd = stats.gaussian_kde(v, bw_method=0.18); ax.plot(xs, kd(xs), color=c, ls=ls, lw=1.2 if k == 'Refined' else 1.0)
        md = UL[(u, k)].median(); ax.axvline(np.log10(md), color=c, lw=0.5, ls=':')
    ax.set_xticks([1, 2, 3, 4]); ax.set_xticklabels(['10', '100', '1k', '10k']); ax.set_xlim(1.3, 3.9); ax.set_ylim(0, None); ax.set_yticks([])
    ax.spines['left'].set_visible(False); ax.set_xlabel(f'{u}UTR length (bp)')
    t = TE[TE.utr == u].set_index('set')
    ax.text(0.98, 0.98, f'median\n{t.loc["Refined annotation", "median"]:.0f} bp', transform=ax.transAxes, ha='right', va='top', fontsize=SM, color=REF, linespacing=1.0)
    ax.text(0.98, 0.70, f'{t.loc["PlasmoDB v71", "median"]:.0f} bp', transform=ax.transAxes, ha='right', va='top', fontsize=SM, color=V71)
    ax.text(0.02, 0.98, u + 'UTR', transform=ax.transAxes, ha='left', va='top', fontsize=SM + 0.4)
axE[0].set_ylabel('Density')
pn.save()
