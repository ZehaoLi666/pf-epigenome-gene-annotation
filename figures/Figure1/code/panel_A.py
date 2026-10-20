"""Figure 1A - Long-read datasets and reference construction (funnel)
Standalone: python panel_A.py   (reads ../data/Figure1_A_datasets.tsv, Figure1_A_reference_funnel.tsv; writes ../panels/Figure1_A.pdf/.svg/.png)
Drawing code copied from the approved make_fig1.py; panel letter is NOT drawn (added in Illustrator)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

DS = read("Figure1_A_datasets.tsv"); FUN = read("Figure1_A_reference_funnel.tsv")
fig, (axA, axF0) = panel_figure('A')
axA.set_xlim(0, 100); axA.set_ylim(0, 100); axA.axis('off')
axA.text(0, 118, "Long-read data and reference construction", fontsize=7, va='bottom')
for i, (t1, l1, l2) in enumerate(zip(DS.dataset, DS.line1, DS.line2)):
    x = 1 + i * 25; box(axA, x, 22, 23, 76, ec='#888888')
    axA.text(x + 11.5, 78, t1, fontsize=6.0, fontweight='bold', ha='center', va='center'); axA.text(x + 11.5, 43, l1 + "\n" + l2, fontsize=4.9, ha='center', va='center', color='#555555', linespacing=1.05)
    arrow(axA, (x + 11.5, 20), (x + 11.5, -6))
steps = [(r.row_label, int(r.n_loci_or_genes), int(r.n_models), None, int(r.step) - 1, r.unit) for r in FUN.itertuples()]
GRP = list(dict.fromkeys(FUN.step_label))
yy = []; yc = 0.0; prev = None
for st in steps:
    if st[4] != prev: yc -= 0.55 if prev is not None else 0; yy.append(None); prev = st[4]
    yc -= 1; yy.append(yc)
hdr = []; ys = []
for k, v in enumerate(yy):
    if v is None: hdr.append(yy[k + 1] + 0.72)
    else: ys.append(v)
y = np.array(ys)
cols = ['#cfcfcf', '#9e9e9e', '#9e9e9e'] + ['#8fa3c7'] * 3 + [COL['ref']]
axF0.barh(y, [s_[1] for s_ in steps], color=cols, height=0.62)
for yi, (nm, l, m, note, g, unit) in zip(y, steps):
    axF0.text(l + 150, yi, f"{l:,} {unit} · {m:,} models", va='center', fontsize=5.1)
for hy, t in zip(hdr, GRP): axF0.text(-0.55, hy, t, transform=axF0.get_yaxis_transform(), fontsize=5.3, fontweight='bold', va='center', ha='left', color='#333333')
axF0.set_yticks(y); axF0.set_yticklabels([s_[0] for s_ in steps], fontsize=5.3); axF0.set_xlim(0, 17000); axF0.set_xticks([0, 5000, 10000]); axF0.spines['bottom'].set_bounds(0, 10000)
axF0.set_ylim(y.min() - 0.6, hdr[0] + 0.5)
axF0.set_xlabel("Loci (first two rows) or genes", fontsize=6); axF0.tick_params(axis='y', length=0)
save_panel(fig, 'A')
