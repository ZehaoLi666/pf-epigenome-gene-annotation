"""Figure 2C - input-track coverage along the reference genes (metagene).
Reads ../data/Figure2_C_metagene_trimmed_means.tsv (5 % trimmed means per bin, computed by
prep_data.py) and ../data/Figure2_C_metadata.tsv. The colour key (bottom-right cell) is part
of this panel, as in the composite. Output: ../panels/Figure2_C.{pdf,svg,png}
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # find style.py from any cwd
from style import *
from matplotlib.patches import Rectangle

fig, AX, to_fig, lay = new_panel("C")
M = data("Figure2_C_metagene_trimmed_means.tsv")
META = dict(data("Figure2_C_metadata.tsv").values)
NG = int(META["n_genes"])
PR = {c: M[c].to_numpy() for c in M.columns if c not in ("bin", "x", "region")}
panels = [('H2A.Z', '#6a51a3'), ('H2B.Z', '#6a51a3'), ('H3K9ac', '#6a51a3'), ('H3K27ac', '#6a51a3'), ('H3K4me3', '#6a51a3'),
          ('ATAC-seq', '#6a51a3'), ('MNase-seq', '#6a51a3'), ('Pol II Ser2P', '#6a51a3'),
          ('GRO-seq', '#2a9d8f'), ('Polysome profiling', '#2a9d8f'), ('DAFT-seq', '#4a78c2')]
xb = M["x"].to_numpy()
seg = [(0, 18, COL['inter']), (18, 22, COL['tss']), (22, 58, COL['body']), (58, 62, COL['tts']), (62, 80, COL['inter'])]
for i, (nm, c) in enumerate(panels):
    ax = AX[f"cell{i:02d}"]
    if nm == 'DAFT-seq':
        g = PR['DAFT-seq, gene strand']; o = PR['DAFT-seq, opposite strand']
        ax.plot(xb, g, color=c, lw=1.0); ax.plot(xb, o, color='#9db3d9', lw=1.0, ls=(0, (2, 1)))
        ax.text(0.98, 0.93, "gene strand", transform=ax.transAxes, ha='right', va='top', fontsize=4.8, color=c)
        ax.text(0.98, 0.78, "opposite", transform=ax.transAxes, ha='right', va='top', fontsize=4.8, color='#7f96bd')
        ymax = max(g.max(), o.max())
    else:
        m = PR[nm]; ax.plot(xb, m, color=c, lw=1.0); ymax = m.max()
    ax.axhline(1, color='#999999', lw=0.5, ls=(0, (2, 2)))
    YM = 4.2; ax.set_ylim(0, YM); ax.set_xlim(0, 80); ax.set_yticks([0, 1, 2, 3, 4])
    if ymax > YM:  # clipped peak (GRO-seq): label its maximum
        k = int(np.argmax(m)); ax.text(xb[k] + 1.5, YM * 0.93, f"peak {ymax:.1f}\u00d7", fontsize=4.8, color=c, va='top')
    for a, b, cc in seg:
        ax.add_patch(Rectangle((a, -0.10 * YM), b - a, 0.06 * YM, fc=cc, ec='none', clip_on=False))
    ax.set_xticks([0, 20, 60, 80]); ax.set_xticklabels(['\u22121 kb', 'TSS', 'TTS', '+1 kb'] if i >= 7 else []); ax.tick_params(axis='x', pad=4.5)
    for xv in (20, 60):
        ax.axvline(xv, color='#cccccc', lw=0.4, zorder=0)
    ax.set_title(nm, fontsize=6.2, color=c, pad=2, loc='left')
    if i % 4 == 0:
        ax.set_ylabel("Fold over\ngenome mean", fontsize=5.4, labelpad=1)
axk = AX["cell11"]; axk.axis('off'); axk.set_xlim(0, 10); axk.set_ylim(0, 10)
for j, (l, cc) in enumerate([("TSS label (\u00b1100 bp)", COL['tss']), ("gene-body label (scaled)", COL['body']),
                             ("TTS label (\u00b1100 bp)", COL['tts']), ("1 kb outside the gene", COL['inter'])]):
    axk.add_patch(Rectangle((0.3, 8.6 - j * 1.9), 1.0, 0.8, fc=cc, ec='none')); axk.text(1.7, 9.0 - j * 1.9, l, fontsize=5.1, va='center')
axk.text(0.3, 0.6, f"{NG:,} reference genes\n5% trimmed mean; 3-stage mean\n(chromatin, GRO-seq, polysome)", fontsize=4.6, va='center', color='#555555', linespacing=1.1)
fx, fy = to_fig(*anchor(lay, "title"))
fig.text(fx, fy, "Input coverage along reference genes", fontsize=7, va='bottom')
save(fig, "Figure2_C")
