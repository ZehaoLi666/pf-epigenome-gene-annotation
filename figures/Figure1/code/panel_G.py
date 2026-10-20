"""Figure 1G - 5′ x 3′ end classes of the 2,389 reference genes
Standalone: python panel_G.py   (reads ../data/Figure1_G_gene_classes.tsv (Figure1_G_class_matrix.tsv = the plotted counts); writes ../panels/Figure1_G.pdf/.svg/.png)
Drawing code copied from the approved make_fig1.py; panel letter is NOT drawn (added in Illustrator)."""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

S = read("Figure1_G_gene_classes.tsv").set_index('gene_id'); NREF = len(S); assert NREF == 2389
c5, c3 = S.c5_100, S.c3_100; FIX = ['over-long', 'too short', 'no annotated UTR', 'between ends']
ndiff = int((c5.isin(FIX) | c3.isin(FIX)).sum())
fig, (axG,) = panel_figure('G')
R5 = ['agree', 'over-long', 'too short', 'no annotated UTR', 'between ends']; C3 = R5 + ['not assessable']
v = pd.crosstab(pd.Categorical(c5, R5), pd.Categorical(c3, C3)).values
assert (v == read("Figure1_G_class_matrix.tsv").iloc[:, 1:].values).all()
axG.imshow(np.log10(v + 1), cmap='Greys', vmin=0, vmax=np.log10(v.max() + 1) * 1.25, aspect='auto')
for i in range(v.shape[0]):
    for j in range(v.shape[1]): axG.text(j, i, f"{v[i,j]:,}", ha='center', va='center', fontsize=5.2, color='white' if v[i, j] > 250 else INK)
lab = ['Agrees', 'Over-\nlong', 'Too\nshort', 'No\nUTR', 'Between\nends', 'Not\nassessed']
axG.set_xticks(range(6)); axG.set_xticklabels(lab, fontsize=4.8); axG.set_yticks(range(5)); axG.set_yticklabels(lab[:5], fontsize=5)
axG.set_xlabel("3′ end vs reference"); axG.set_ylabel("5′ end vs reference")
for k in ['top', 'right']: axG.spines[k].set_visible(True)
axG.add_patch(Rectangle((-0.5, -0.5), 1, 1, fill=False, ec=COL['agree'], lw=1.2))
title(axG, f"{ndiff:,} of {NREF:,} reference genes ({ndiff/NREF:.0%})\ndiffer from v71 at ≥1 end")
save_panel(fig, 'G')
