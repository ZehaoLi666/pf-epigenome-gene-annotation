"""Figure 4C - independent 5'-end (Adjalley 2016, Chappell 2020, Shaw 2021) and poly(A) (Siegel 2014) data centred on the
refined (navy) and PlasmoDB v71 (orange) TSS/TTS. Mean per-gene-normalised signal summed in 10-bp bins, shown within ±250 bp;
numbers at upper right = % of signal (or genes) within ±100 bp. The line key above the four plots is part of this panel.
Data: ../data/Figure4_C_profiles.tsv, Figure4_C_refined_vs_v71_signal.tsv, Figure4_C_anchor_gene_counts.tsv
Output: ../panels/Figure4_C.pdf|.svg|.png"""
import numpy as np, pandas as pd, matplotlib as mpl
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))  # find style.py
from style import Panel, DATA, REF, V71, SM

PR = pd.read_csv(DATA / 'Figure4_C_profiles.tsv', sep='\t')
TC = pd.read_csv(DATA / 'Figure4_C_refined_vs_v71_signal.tsv', sep='\t')
NS = pd.read_csv(DATA / 'Figure4_C_anchor_gene_counts.tsv', sep='\t').set_index(['anchor', 'end']).n_genes
OFF = PR.offset_bp.values
DSETS = ['Adjalley 2016', 'Chappell 2020', 'Shaw 2021', 'Siegel 2014']
ASSAY = {'Adjalley 2016': 'cap-seq', 'Chappell 2020': '5′UTR-seq', 'Shaw 2021': '5CAPture-seq', 'Siegel 2014': 'major poly(A) site'}
ALAB = {'Refined': 'Refined annotation', 'v71': 'PlasmoDB v71'}
BIN, XL = 10, 250
def binned(p):
    k = np.abs(OFF) <= XL; e = np.arange(-XL, XL + BIN, BIN); idx = np.clip(np.digitize(OFF[k], e) - 1, 0, len(e) - 2)
    return (e[:-1] + e[1:]) / 2, np.bincount(idx, p[k], len(e) - 1) * 100

pn = Panel('C'); fig = pn.fig
for i, ds in enumerate(DSETS):
    ax = pn.add_axes([0.085 + i * 0.235, 0.565, 0.18, 0.135])
    for k, c, lw in (('v71', V71, 1.0), ('Refined', REF, 1.3)):
        xc, y = binned(PR[f'{ds}|{ALAB[k]}'].values); ax.plot(xc, y, color=c, lw=lw)
    ax.axvline(0, color='#bbbbbb', lw=0.5, ls=':', zorder=0); ax.set_xlim(-XL, XL); ax.set_ylim(0, None)
    ax.set_xticks([-250, -100, 0, 100, 250]); ax.set_xticklabels(['−250', '', '0', '', '+250'])
    ax.text(0.03, 0.98, f'{ds}\n{ASSAY[ds]}', transform=ax.transAxes, fontsize=SM + 0.2, va='top', ha='left', linespacing=1.05)
    s = TC[TC.dataset == ds].set_index('anchor')
    for j, (k, c) in enumerate((('Refined', REF), ('v71', V71))):
        ax.text(0.98, 0.98 - j * 0.09, f'{s.loc[ALAB[k], "frac_pm100"] * 100:.0f}%', transform=ax.transAxes, fontsize=SM, color=c, va='top', ha='right',
                fontweight='bold' if k == 'Refined' else 'normal')
    ax.set_xlabel('Distance from TSS (bp)' if i < 3 else 'Distance from TTS (bp)', fontsize=6)
    if i == 0: ax.set_ylabel('5′-end signal\n(% per 10 bp)', linespacing=1.05)
    if i == 3: ax.set_ylabel('Genes with site\n(% per 10 bp)', linespacing=1.05)
hC = [mpl.lines.Line2D([], [], color=REF, lw=1.5, label=f'Refined annotation (TSS n = {NS[("Refined annotation", "TSS")]:,}; TTS n = {NS[("Refined annotation", "TTS")]:,})'),
      mpl.lines.Line2D([], [], color=V71, lw=1.2, label=f'PlasmoDB v71, all genes (TSS n = {NS[("PlasmoDB v71", "TSS")]:,}; TTS n = {NS[("PlasmoDB v71", "TTS")]:,})')]
fig.legend(handles=hC, frameon=False, fontsize=SM, loc='lower left', ncol=2, handlelength=1.4, columnspacing=1.5, bbox_to_anchor=pn.xy(0.08, 0.706), borderaxespad=0)
pn.save()
