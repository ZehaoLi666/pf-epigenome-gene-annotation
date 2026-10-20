import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

"""Panel J: upstream AUGs in distal vs proximal 5' leaders for the 24 TSS-proportion-switch pairs
(dark: switch only; light: switch + expression). Paired Wilcoxon signed-rank test."""
d = data("Figure5_J_uAUG_pairs.tsv"); st = data("Figure5_statistics.tsv").set_index("panel")
fig, (ax,) = new_panel("J")
for r in d.itertuples():
    ax.plot([0, 1], [r.distal_n_uAUG, r.proximal_n_uAUG], color=NAVY if r.case.startswith("2") else LT, lw=0.6, marker="o", ms=2)
ax.set_xticks([0, 1]); ax.set_xticklabels(["Distal", "Proximal"], fontsize=6); ax.set_xlim(-0.3, 1.3); ax.set_ylabel("uAUGs per 5′ leader")
ax.set_title(f"uAUGs (n = {len(d)})\np = {st.loc['J','p']:.0e}".replace("e-0", "×10⁻"), fontsize=6)
save(fig, "Figure5_J")
