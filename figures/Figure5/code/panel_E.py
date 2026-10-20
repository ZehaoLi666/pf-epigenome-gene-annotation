import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

"""Panel E: DAFT-seq co-expression (Pearson r, 7 time points) of neighbouring genes by orientation,
with random inter-chromosomal pairs as null. Horizontal bars and numbers = medians."""
d = data("Figure5_E_neighbour_coexpression.tsv")
groups = ["divergent", "tandem", "convergent", "random_interchromosomal"]
vals = [d.pearson_r[d.group == g].values for g in groups]
fig, (ax,) = new_panel("E")
vp = ax.violinplot(vals, positions=range(4), showextrema=False, widths=0.8)
for b_, c in zip(vp["bodies"], [NAVY, LT, LT, LG]):
    b_.set_facecolor(c); b_.set_alpha(0.8); b_.set_edgecolor("none")
for i, v in enumerate(vals):
    ax.plot([i - 0.25, i + 0.25], [np.median(v)] * 2, color="k", lw=1)
    ax.text(i, 1.08, f"{np.median(v):.2f}", ha="center", fontsize=5.5)
ax.set_xticks(range(4)); ax.set_xticklabels(["Head-\nto-head", "Tan-\ndem", "Tail-\nto-tail", "Ran-\ndom"], fontsize=5.5)
ax.set_ylabel("Pearson r (DAFT-seq)"); ax.set_ylim(-1.05, 1.2); ax.set_title("Neighbour co-expression", fontsize=6.5)
save(fig, "Figure5_E")
