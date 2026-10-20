import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

"""Panel C: MNase occupancy around refined, PlasmoDB v71 and experimental TSSs (1,236 trimmed genes)."""
d = data("Figure5_C_mnase_metagene.tsv"); n = data("Figure5_C_mnase_n.tsv", index_col=0)["n_genes"]
fig, (ax,) = new_panel("C")
for an in ["experimental", "v71", "refined"]:
    y = smooth(d[an].values, 11)
    ax.plot(d.bp_from_TSS, y, color=ACOL[an], lw=1.1 if an != "experimental" else 0.9,
            ls="-" if an != "experimental" else "--", label=f"{ALAB[an]} (n={n[an]})")
ax.axvline(0, color="k", lw=0.5, ls=":"); ax.set_xlim(-300, 300)
ax.set_title("MNase occupancy", fontsize=6.5); ax.set_ylabel("Occupancy / median"); ax.set_xlabel("bp from TSS")
ax.axhline(0.5, color="0.7", lw=0.5, ls="--")
save(fig, "Figure5_C")
