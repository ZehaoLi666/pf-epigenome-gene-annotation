"""Figure 1B - Read-level long-read coverage (top: genome; bottom: reads per gene)
Standalone: python panel_B.py   (reads ../data/Figure1_B_genome_coverage.tsv, Figure1_B_gene_longread_counts.tsv; writes ../panels/Figure1_B.pdf/.svg/.png)
Drawing code copied from the approved make_fig1.py; panel letter is NOT drawn (added in Illustrator)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

GC = read("Figure1_B_genome_coverage.tsv"); PG = read("Figure1_B_gene_longread_counts.tsv", index_col=0)
fig, (axB, axB2) = panel_figure('B')
lab = dict(zip(GC.source, GC.label))
order = ['ALL', 'Shaw2022_PacBio', 'Shaw2022_ONT', 'Lee2021', 'Hoshizaki2022', 'Yang2021']; yb = np.arange(len(order))[::-1]
g1 = GC[GC.min_depth == 1].set_index('source').genome; g5 = GC[GC.min_depth == 5].set_index('source').genome
axB.barh(yb, [g1[o] * 100 for o in order], color='#c6d3e8', height=0.66)
axB.barh(yb, [g5[o] * 100 for o in order], color=COL['ref'], height=0.66)
for yi, o in zip(yb, order): axB.text(g1[o] * 100 + 1.5, yi, f"{g1[o]:.0%}", va='center', fontsize=4.8)
axB.set_yticks(yb); axB.set_yticklabels([lab[o] for o in order], fontsize=5.0); axB.tick_params(axis='y', length=0)
axB.get_yticklabels()[0].set_fontweight('bold')
axB.set_xlim(0, 112); axB.set_xticks([0, 50, 100]); axB.spines['bottom'].set_bounds(0, 100); axB.set_xlabel("Nuclear genome covered (%)", fontsize=5.6, labelpad=1)
axB.text(0.99, 1.02, "≥1 read", transform=axB.transAxes, color='#7f96bd', fontsize=4.8, ha='right', va='bottom'); axB.text(0.70, 1.02, "≥5 reads", transform=axB.transAxes, color=COL['ref'], fontsize=4.8, ha='right', va='bottom')
title(axB, "Long reads cover most of the genome", y=1.13)
grp = [('confirmed', 'Reference', COL['ref']), ('truncated-only', 'Only truncated models', '#8c8c8c'), ('HC model, unconfirmed', 'Not stage-confirmed', '#a9a9a9'), ('no HC model', 'No supported model', '#cfcfcf')]
yg = np.arange(len(grp))[::-1]
for yi, (k, l, cc) in zip(yg, grp):
    v = np.log10(PG.loc[PG.lr_group == k, 'pooled_reads_on_CDS'] + 1)
    axB2.boxplot([v], positions=[yi], orientation='horizontal', widths=0.58, showfliers=False, patch_artist=True, medianprops=dict(color='white' if cc == COL['ref'] else INK, lw=0.9),
                 boxprops=dict(fc=cc, ec=INK, lw=0.5), whiskerprops=dict(lw=0.5), capprops=dict(lw=0.5))
    axB2.text(4.62, yi, f"{int(10**np.median(v)-0.5):,}", va='center', fontsize=4.9, ha='right')
axB2.set_yticks(yg); axB2.set_yticklabels([f"{l} ({(PG.lr_group==k).sum():,})" for k, l, _ in grp], fontsize=5.0); axB2.tick_params(axis='y', length=0)
axB2.set_xlim(0, 4.7); axB2.set_xticks([0, 1, 2, 3, 4]); axB2.set_xticklabels(['1', '10', '100', '1k', '10k']); axB2.spines['bottom'].set_bounds(0, 4)
axB2.set_xlabel("Long reads on the CDS per gene (pooled)", fontsize=5.6, labelpad=1)
title(axB2, "Non-reference genes also have reads", y=1.02)
save_panel(fig, 'B')
