import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

"""Panel D: GTGCA(C) motif around refined, PlasmoDB v71 and experimental TSSs (1,236 trimmed genes)."""
d = data("Figure5_D_motif_GTGCA_per_kb.tsv"); n = data("Figure5_D_motif_n.tsv", index_col=0)["n_genes"]
fig, (ax,) = new_panel("D")
for an in ["experimental", "v71", "refined"]:
    y = smooth(d[an].values, 101)
    ax.plot(d.bp_from_TSS, y, color=ACOL[an], lw=1.1 if an != "experimental" else 0.9,
            ls="-" if an != "experimental" else "--", label=f"{ALAB[an]} (n={n[an]})")
ax.axvline(0, color="k", lw=0.5, ls=":"); ax.set_xlim(-1000, 500)
ax.set_title("GTGCA(C) motif", fontsize=6.5); ax.set_ylabel("Motifs per kb"); ax.set_xlabel("bp from TSS")

save(fig, "Figure5_D")
