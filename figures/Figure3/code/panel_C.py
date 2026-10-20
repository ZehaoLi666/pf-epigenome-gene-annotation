"""Figure 3C - three example loci on test chromosomes, drawn 5'->3'.

Rows: PlasmoDB v71 models | stage-confirmed long-read models | model calls |
DAFT-seq | H2A.Z | H3K9ac | H3K4me3 | ATAC-seq | CAGE.
The shared symbol key below the panel is a separate file (legend_key.py).
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))   # find style.py next to this script
import numpy as np, pandas as pd
from matplotlib.patches import Rectangle
from style import Panel, COL, INK, data

EXM = pd.read_csv(data("Figure3_C_examples.tsv"), sep="\t")
FEAT = pd.read_csv(data("Figure3_C_v71_features.tsv"), sep="\t")
MOD = pd.read_csv(data("Figure3_C_longread_models.tsv"), sep="\t")
CLU = pd.read_csv(data("Figure3_C_longread_end_clusters.tsv"), sep="\t")
CAL = pd.read_csv(data("Figure3_C_model_calls.tsv"), sep="\t")
TRK = pd.read_csv(data("Figure3_C_signal_tracks.tsv"), sep="\t")

# axes geometry in composite inches: one column per example, nine rows
COL_X = [0.991187, 3.010920, 5.030653]; COL_W = 1.845347
ROWS = [(6.744967, 0.211449), (6.231697, 0.480566), (5.929876, 0.269117), (5.628055, 0.269117), (5.356990, 0.238361),
        (5.085925, 0.238361), (4.814860, 0.238361), (4.543795, 0.238361), (4.241974, 0.269117)]
rowlab = ["PlasmoDB v71", "Stage-confirmed\nlong-read models", "Model calls", "Stranded RNA\n(DAFT-seq)", "H2A.Z", "H3K9ac", "H3K4me3", "ATAC-seq", "CAGE\n(not an input)"]
SIG = [(3, 'DAFTseq_log10_1p_CPM', COL['daft'], 0), (4, 'H2AZ_minusInput_mean_RTS_clip0', COL['chromatin_track'], 0),
       (5, 'H3K9ac_minusInput_mean_RTS_clip0', COL['chromatin_track'], 0), (6, 'H3K4me3_minusInput_mean_RTS_clip0', COL['chromatin_track'], 0),
       (7, 'ATAC_minusInput_mean_RTS_clip0', COL['chromatin_track'], 0), (8, 'CAGE_CPM', COL['cage'], 0.5)]
RANGE = {3: (COL['daft'], 1), 4: (COL['chromatin_track'], 0), 5: (COL['chromatin_track'], 0), 6: (COL['chromatin_track'], 0),
         7: (COL['chromatin_track'], 0), 8: (COL['cage'], 0)}   # colour, decimals of the "0-max" label

P = Panel('C')
ROWAX = {}
for j, e in enumerate(EXM.itertuples()):
    a, b, s, gname = e.window_start, e.window_end, e.strand, e.gene_id
    axs = [P.axes((COL_X[j], y, COL_W, h)) for y, h in ROWS]
    for ax in axs:
        ax.set_xlim((b, a) if s == '-' else (a, b)); ax.set_xticks([]); ax.set_yticks([])
        for sp in ['left', 'bottom', 'top', 'right']: ax.spines[sp].set_visible(False)
        ax.axvline(e.main_longread_TSS, color=COL['tss'], lw=0.5, ls=(0, (2, 2)), alpha=0.8, zorder=0)
    # --- PlasmoDB v71 gene models
    ax = axs[0]; ax.set_ylim(-1, 1.3)
    gg = FEAT[FEAT.example == e.example]
    for gn_, sb in gg.groupby('gene'):
        foc = gn_ == gname; alp = 1 if foc else 0.35; ex_ = sb[sb.type == 'exon']
        if len(ex_) > 1: ax.plot([ex_.start.min(), ex_.end.max()], [0, 0], color=COL['cds'], lw=0.5, alpha=alp)
        for f in sb[sb.type != 'exon'].itertuples():
            h = 0.9 if f.type == 'CDS' else 0.5
            ax.add_patch(Rectangle((f.start, -h / 2), f.end - f.start, h, fc=COL['cds'] if f.type == 'CDS' else COL['utr'], ec='none', alpha=alp))
        if foc: ax.text((max(sb.start.min(), a) + min(sb.end.max(), b)) / 2, 0.75, gn_, ha='center', fontsize=5.3, fontweight='bold')
    # --- long-read models (dark: complete CDS; light: truncated); dots = end within 200 bp of a call
    ax = axs[1]; lr = MOD[MOD.example == e.example].sort_values('row')
    for t_ in lr.itertuples():
        k = t_.row; cc = COL['model_full'] if t_.full_cds else COL['model_trunc']
        ax.plot([t_.start, t_.end], [-k, -k], color=cc, lw=0.4)
        for bb in t_.blocks.split(';'):
            ax.plot([int(v) for v in bb.split('-')], [-k, -k], color=cc, lw=2.0, solid_capstyle='butt')
        if t_.end5_within_200bp_of_TSS_call: ax.plot(t_.end5, -k, 'o', ms=2.6, color=COL['tss'], zorder=5)
        if t_.end3_within_200bp_of_TTS_call: ax.plot(t_.end3, -k, 'o', ms=2.6, color=COL['tts'], zorder=5)
    yb = -max(len(lr), 1) - 0.2
    cl = CLU[CLU.example == e.example]
    for p in cl[cl.end == 'TSS'].position: ax.plot([p, p], [yb - 0.35, yb + 0.35], color=COL['tss'], lw=1.2)
    for p in cl[cl.end == 'TTS'].position: ax.plot([p, p], [yb - 0.35, yb + 0.35], color=COL['tts'], lw=1.2)
    ax.set_ylim(yb - 0.8, 0.8)
    # --- model calls (filled: main; open: other, coloured if a long-read end lies within 200 bp)
    ax = axs[2]; ax.set_ylim(-1, 1.6); ax.axhline(0, color='#bbbbbb', lw=0.4)
    ca = CAL[CAL.example == e.example]
    for end in ('TSS', 'TTS'):
        for q in ca[(ca.end == end) & (ca.call == 'other')].itertuples():
            c_ = COL[end.lower()] if str(q.longread_end_within_200bp) == 'True' else COL['grey_call']
            ax.plot(q.position, 0.45, marker='v', ms=3.4, mfc='white', mec=c_, mew=0.8, ls='none')
    ax.plot(e.main_TSS_call, 0.45, marker='v', ms=4.6, color=COL['tss'], ls='none')
    ax.plot(e.main_TTS_call, 0.45, marker='v', ms=4.6, color=COL['tts'], ls='none')
    ax.text(e.main_TSS_call, 0.95, f"TSS p={e.main_TSS_prob:.2f}", fontsize=4.7, ha='center', color=COL['tss'])
    # --- signal tracks
    tk = TRK[TRK.example == e.example].sort_values('bin'); x = tk.position.values
    for row, colname, cc, lw in SIG:
        v = tk[colname].values; ax = axs[row]
        ax.fill_between(x, 0, v, color=cc, lw=lw); ROWAX.setdefault(row, []).append((ax, v.max()))
    # --- scale bar and locus label under the CAGE row
    ax = axs[8]
    sb0 = b - 120 if s == '-' else a + 120; sb1 = sb0 - 500 if s == '-' else sb0 + 500
    ax.plot([sb0, sb1], [-0.15] * 2, color=INK, lw=0.8, clip_on=False, transform=ax.get_xaxis_transform())
    ax.text((sb0 + sb1) / 2, -0.27, "500 bp", transform=ax.get_xaxis_transform(), fontsize=4.8, ha='center', va='top', clip_on=False)
    ax.text(0.99, -0.2, f"{e.chrom.replace('Pf3D7_', 'chr').replace('_v3', '')} ({s}), test chromosome", transform=ax.transAxes, fontsize=4.6, ha='right', va='top', color='#666666')
    axs[0].text(0, 1.25, e.title, transform=axs[0].transAxes, fontsize=6.0, va='bottom')
    if j == 0:
        for i, l in enumerate(rowlab):
            axs[i].text(-0.04, 0.5, l, transform=axs[i].transAxes, fontsize=5.2, ha='right', va='center', linespacing=0.95)
# same y scale per row across the three examples; range printed at the right of each row
for row, lst in ROWAX.items():
    ym = max(m for _, m in lst) * 1.1
    for ax_, _ in lst: ax_.set_ylim(0, ym)
    cc, dec = RANGE[row]
    lst[-1][0].text(1.005, 0.95, f"0–{ym / 1.1:.{dec}f}", transform=lst[-1][0].transAxes, fontsize=4.6, ha='left', va='top', color=cc)
P.save()
