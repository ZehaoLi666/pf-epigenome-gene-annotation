import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from style import *

"""Panel F: MNase occupancy (log2 occupancy / genome median, Kensche 2016 via nucDetective, mean of
8 time points) across head-to-head pairs centred on the midpoint between the two refined TSSs,
sorted by TSS distance. Black dots = refined TSSs."""
z = np.load(os.path.join(DATA, "Figure5_F_mnase_heatmap_matrix.npz"))
M = z["log2_occ_over_median"].astype(np.float32); W = int(z["bp_from_midpoint"].max()) + 1
o = data("Figure5_F_pair_order.tsv")
fig, (ax,) = new_panel("F")
im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2, extent=[-W, W, len(M), 0], interpolation="nearest")
hd = (o.refined_TSS_right - o.refined_TSS_left).values / 2; yy = np.arange(len(o)) + 0.5
ax.scatter(-hd, yy, s=0.3, c="k", lw=0); ax.scatter(hd, yy, s=0.3, c="k", lw=0)
ax.set_xlim(-W, W); ax.set_xticks([-1500, 0, 1500]); ax.set_xticklabels(["−1.5", "mid", "+1.5"])
ax.set_xlabel("kb from midpoint"); ax.set_ylabel(f"Head-to-head pairs (n = {len(M)}),\nsorted by TSS distance")
ax.set_title("MNase occupancy between refined TSSs", fontsize=6.5)
cax = ax.inset_axes([1.04, 0.55, 0.04, 0.45]); cb = fig.colorbar(im, cax=cax, orientation="vertical")
cb.set_ticks([-2, 0, 2]); cb.ax.tick_params(labelsize=5, pad=1); cb.set_label("log2 occ.\n/ median", fontsize=5, labelpad=1)
save(fig, "Figure5_F")
