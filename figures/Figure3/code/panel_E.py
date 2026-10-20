"""Figure 3E - input ablation scored against long-read ends on held-out chromosomes.
Points: change in accuracy within 200 bp vs the full model retrained in the same run;
bars: paired bootstrap 95% CI; filled: p < 0.05 (two-sided sign test)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # find style.py next to this script
import numpy as np, pandas as pd, matplotlib as mpl
from style import Panel, COL, INK, data

AB3 = pd.read_csv(data("Figure3_E_input_ablation.tsv"), sep="\t")
ABLM = [('L_nochrom', 'Chromatin'), ('L_nodaft', 'DAFT-seq'), ('L_nostrand', 'Inferred stranded RNA'), ('N_noshortread', 'DAFT-seq + inferred stranded'),
        ('L_nogro', 'GRO-seq'), ('L_nopoly', 'Polysome'), ('L_noonehot', 'DNA one-hot'), ('L_noglm', 'gLM embeddings'),
        ('S_chrom', 'Chromatin only'), ('S_chromdna', 'Chromatin + DNA only'), ('S_rna', 'RNA only'), ('S_dna', 'DNA + gLM only')]
C_CH = COL['ref']

P = Panel('E')
axE2 = [P.axes((5.219265, 2.237985, 0.738367, 1.167210)), P.axes((6.137633, 2.237985, 0.738367, 1.167210))]
yA = np.array(list(range(13, 5, -1)) + [3.9, 2.9, 1.9, 0.9], float)
tr = mpl.transforms.blended_transform_factory(axE2[0].transAxes, axE2[0].transData)
for ax, tl, xl in ((axE2[0], 'Long-read TSS', (-0.10, 0.03)), (axE2[1], 'Long-read TTS', (-0.24, 0.03))):
    d = AB3[AB3.truth == tl].set_index('model').loc[[m for m, _ in ABLM]]
    ax.axvline(0, color='#999999', lw=0.5, zorder=0); ax.axhline(5.45, color='#e5e5e5', lw=0.5)
    for y, (m, lab) in zip(yA, ABLM):
        r = d.loc[m]; cc = C_CH if 'hrom' in lab else INK
        ax.plot([r.lo, r.hi], [y, y], color=cc, lw=0.8, solid_capstyle='butt')
        ax.plot(r.delta, y, 'o', ms=2.8, mfc=cc if r.p_sign < 0.05 else 'white', mec=cc, mew=0.8, zorder=3)
    ax.set_xlim(*xl); ax.set_ylim(0.3, 14.8); ax.set_yticks([])
    ax.xaxis.set_major_formatter(mpl.ticker.FuncFormatter(lambda v, _: f'{v * 100:+.0f}'.replace('+0', '0').replace('-', '\u2212')))
    ax.text(0.5, 1.0, f'{tl}\nfull {d.acc_full.iloc[0]:.3f} (n = {d.n.iloc[0]:,})', transform=ax.transAxes, ha='center', va='bottom', fontsize=5.0, linespacing=1.05,
            color=COL['tss'] if 'TSS' in tl else COL['tts'])
for y, (m, lab) in zip(yA, ABLM):
    axE2[0].text(-0.06, y, lab, transform=tr, ha='right', va='center', fontsize=5.0, color=C_CH if 'hrom' in lab else INK)
axE2[0].text(-0.06, 14.3, 'Remove one class', transform=tr, ha='right', va='center', fontsize=5.0, fontweight='bold')
axE2[0].text(-0.06, 4.85, 'Use only', transform=tr, ha='right', va='center', fontsize=5.0, fontweight='bold')
P.text(6.047633, 1.949985, 'Δ accuracy within 200 bp vs full model (points)', ha='center', fontsize=5.6)
axE2[1].plot([], [], 'o', ms=2.8, mfc=INK, mec=INK, ls='', label='p < 0.05'); axE2[1].plot([], [], 'o', ms=2.8, mfc='white', mec=INK, ls='', label='p \u2265 0.05')
axE2[1].legend(loc='lower left', frameon=False, fontsize=4.8, handletextpad=0.2, borderaxespad=0.1, labelspacing=0.25)
P.text(4.139265, 3.662706, 'Input ablation: chromatin selects the TSS', fontsize=7, ha='left', va='bottom')
P.save()
