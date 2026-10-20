import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

"""Panel A: ATAC accessibility around refined, PlasmoDB v71 and experimental TSSs (1,236 trimmed genes)."""
d = data("Figure5_A_atac_metagene.tsv"); n = data("Figure5_A_atac_n.tsv", index_col=0)["n_genes"]
fig, (ax,) = new_panel("A")
for an in ["experimental", "v71", "refined"]:
    y = smooth(d[an].values, 51)
    ax.plot(d.bp_from_TSS, y, color=ACOL[an], lw=1.1 if an != "experimental" else 0.9,
            ls="-" if an != "experimental" else "--", label=f"{ALAB[an]} (n={n[an]})")
ax.axvline(0, color="k", lw=0.5, ls=":"); ax.set_xlim(-1000, 500)
ax.set_title("ATAC accessibility", fontsize=6.5); ax.set_ylabel("Signal / genome mean"); ax.set_xlabel("bp from TSS")

save(fig, "Figure5_A")
