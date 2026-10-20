import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

"""Panel G: mean DAFT-seq co-expression (+/- 95% CI = 1.96 x SEM) of head-to-head pairs by TSS
distance, split by the number of MNase-defined NDRs between the two TSSs (bins with n < 8 hidden)."""
d = data("Figure5_G_pairs_by_distance_and_NDR.tsv"); st = data("Figure5_statistics.tsv").set_index("panel")
lab = ["≤200", "200–400", "400–600", "600–800", "0.8–1k", "1–1.5k", "1.5–2k", "2–3k", ">3k"]
d = d[d.pearson_r.notna() & d.distance_bin.notna()]
fig, (ax,) = new_panel("G")
for cl, col, mk, off in [("1 NDR", NAVY, "o", -0.1), ("≥2 NDRs", LT, "s", 0.1)]:
    q = d[d.ndr_class == cl].groupby("distance_bin").pearson_r.agg(["mean", "sem", "size"]); q = q[q["size"] >= 8]
    q = q.loc[[b for b in lab if b in q.index]]
    xi = np.array([lab.index(i) for i in q.index]) + off
    ax.errorbar(xi, q["mean"], yerr=1.96 * q["sem"], fmt=mk + "-", color=col, ms=3, lw=1, elinewidth=0.6, capsize=1.5, label=cl + " between TSSs")
ax.set_xticks(np.arange(len(lab))); ax.set_xticklabels(lab, rotation=50, ha="right", fontsize=5.5)
ax.set_xlabel("Distance between TSSs (bp)"); ax.set_ylabel("Mean r ± 95% CI"); ax.legend(frameon=False, fontsize=5.5, loc="lower left")
ax.set_title(f"Co-expression follows distance,\nnot NDR number (300–800 bp: p = {st.loc['G','p']:.2f})", fontsize=6.5)
save(fig, "Figure5_G")
