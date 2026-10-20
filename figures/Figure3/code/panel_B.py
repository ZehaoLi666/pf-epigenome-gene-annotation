"""Figure 3B - three-fold chromosome cross-validation (chromosome lengths to scale)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # find style.py next to this script
import pandas as pd
from matplotlib.patches import Rectangle
from style import Panel, COL, data

CH = pd.read_csv(data("Figure3_B_chromosome_folds.tsv"), sep="\t")
NG = pd.read_csv(data("Figure3_B_test_gene_counts.tsv"), sep="\t")
ROLE_COL = {'test': COL['ref'], 'validation': COL['validation'], 'training': COL['training']}

P = Panel('B')
axB = P.axes((5.195020, 7.603185, 1.680980, 0.991815))
for i, fr in enumerate(NG.itertuples()):
    y = 2 - i; x = 0; role = CH[fr.fold.replace(' ', '').lower() + '_role']
    for L, ro in zip(CH.length_bp, role):
        wd = L / 1e6
        axB.add_patch(Rectangle((x, y - 0.32), wd, 0.64, fc=ROLE_COL[ro], ec='white', lw=0.6)); x += wd
    axB.text(x + 0.25, y, f"{fr.n_reference_genes:,}", fontsize=5.2, va='center')
    axB.text(-0.25, y, fr.fold, fontsize=5.4, ha='right', va='center')
x = 0
for k, L in zip(CH.chromosome, CH.length_bp):
    if k in (1, 5, 6, 10, 11, 14): axB.text(x + L / 2e6, -0.75, str(k), fontsize=4.8, ha='center', va='top')
    x += L / 1e6
axB.set_xlim(0, x + 2.4); axB.set_ylim(-1.3, 2.7); axB.axis('off')
for j, (l, c) in enumerate([("test", COL['ref']), ("validation", COL['validation']), ("training", COL['training'])]):
    # key swatches sit at the top edge of the axes and are clipped to a thin bar, as in the composite
    axB.add_patch(Rectangle((j * 6.4, 2.62), 0.9, 0.3, fc=c, ec='none')); axB.text(j * 6.4 + 1.2, 2.77, l, fontsize=5.0, va='center')
axB.text(x + 0.25, 2.77, "genes", fontsize=4.8, va='center', color='#555555')
axB.text(x / 2, -1.25, "Chromosome (length to scale)", fontsize=5.2, ha='center')
axB.text(0, 1.12, "Three-fold chromosome cross-validation", transform=axB.transAxes, fontsize=7, va="bottom")
P.save()
