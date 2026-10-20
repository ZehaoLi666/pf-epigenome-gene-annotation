import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

"""Panel I: distal-TSS proportion (CAGE solid, 5'UTR-seq dashed) and DAFT-seq expression (bars,
divided by the maximum) in ring (R), trophozoite (T) and schizont (S) for PfA-M1, LDH and eIF3i.
Titles: distal / proximal 5'UTR length (nt)."""
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
d = data("Figure5_I_focal_genes.tsv"); st = ["ring", "troph", "schizont"]
fig, axes = new_panel("I")
for k, (ax, r) in enumerate(zip(axes, d.itertuples())):
    r = r._asdict(); cpm = np.array([r["DAFT_CPM_" + s] for s in st])
    ax.bar(range(3), cpm / cpm.max(), color="0.85", width=0.6)
    ax.plot(range(3), [r["CAGE_prop_" + s] for s in st], "o-", color=NAVY, ms=3, lw=1)
    ax.plot(range(3), [r["U5_prop_" + s] for s in st], "s--", color=LT, ms=3, lw=1)
    ax.set_xticks(range(3)); ax.set_xticklabels(["R", "T", "S"], fontsize=6); ax.set_ylim(0, 1.05)
    ax.set_title(f"{r['gene_name']}\n{int(r['distal_utr5_len'])}/{int(r['proximal_utr5_len'])} nt", fontsize=5.8)
    if k == 0: ax.set_ylabel("Distal-TSS proportion /\nrelative expression")
    else: ax.set_yticklabels([])
    if k == 1:
        ax.legend(handles=[Line2D([], [], color=NAVY, marker="o", ms=3, label="CAGE"), Line2D([], [], color=LT, marker="s", ls="--", ms=3, label="5′UTR-seq"),
                           Patch(color="0.85", label="DAFT (÷ max)")], frameon=False, fontsize=5, loc="upper center",
                  bbox_to_anchor=(0.5, -0.14), ncol=3, handlelength=1.4, columnspacing=0.8)
save(fig, "Figure5_I")
