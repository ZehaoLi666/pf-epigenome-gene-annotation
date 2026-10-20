import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

"""Panel H: change in distal-TSS proportion (CAGE) vs change in gene expression (DAFT-seq log2FC) for
65 stage-specific TSS pairs (MSP5 excluded: MSP4 read-through) and 888 single-TSS controls."""
d = data("Figure5_H_gradeAB_pairs.tsv"); c = data("Figure5_H_singleTSS_controls.tsv")
fig, (ax,) = new_panel("H")
ax.scatter(c.DAFT_log2FC.abs(), c.CAGE_dprop.abs(), s=2, c=LG, lw=0, label=f"Single-TSS controls (n={len(c)})", rasterized=True)
CC = {"2_TSS_switch_only": (NAVY, "Switch only"), "3_both": (LT, "Switch + expression"), "1_expression_only": (GRY, "Expression only"),
      "4_neither": ("white", "Neither"), "5_unclassified": ("white", "Unclassified")}
for cs, (col, lb) in CC.items():
    x = d[d.case == cs]
    ax.scatter(x.DAFT_log2FC.abs(), x.CAGE_dprop.abs(), s=10 if cs[0] in "23" else 6, c=col,
               edgecolors="k" if col == "white" else "none", lw=0.4, label=f"{lb} ({len(x)})", zorder=3 if cs[0] in "23" else 2)
ax.axhline(0.2, color="k", lw=0.5, ls=":"); ax.axvline(1, color="k", lw=0.5, ls=":")
for g, nm in [("PF3D7_1311800", "PfA-M1"), ("PF3D7_1324900", "LDH"), ("PF3D7_0716800", "eIF3i")]:
    r = d[d.gene_id == g].iloc[0]
    ax.annotate(nm, (abs(r.DAFT_log2FC), abs(r.CAGE_dprop)), xytext=(8, -2), textcoords="offset points", fontsize=5.5, arrowprops=dict(arrowstyle="-", lw=0.4))
ax.set_xlabel("|Expression change| (DAFT log2FC)"); ax.set_ylabel("|Δ distal proportion| (CAGE)"); ax.set_xlim(0, 4.6); ax.set_ylim(0, 1.12)
ax.legend(frameon=False, fontsize=5, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, handletextpad=0.2, columnspacing=0.8)
ax.set_title("Stage-specific TSS changes (65 pairs)", fontsize=6.5)  # 65 Grade A/B pairs; MSP5 (excluded) not drawn
save(fig, "Figure5_H")
