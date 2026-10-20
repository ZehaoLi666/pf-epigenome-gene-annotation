"""Shared drawing routine for the genome-track example panels J, K and L (copied from make_fig1.py).
Track row labels (PlasmoDB v71, Long-read models, ...) are drawn separately by panel_JKL_track_labels.py."""
from style import *

ROWLAB = ["PlasmoDB v71", "Long-read\nmodels", "lncRNA", "Stranded RNA\n(DAFT-seq)", "CAGE"]

def draw_track_panel(P):
    M = read(f"Figure1_{P}_meta.tsv").iloc[0]
    gff = read(f"Figure1_{P}_v71_gene_models.tsv").rename(columns={'gene_id': 'gene'})
    LRm = read(f"Figure1_{P}_longread_models.tsv"); UN = read(f"Figure1_{P}_lncRNA.tsv"); CV = read(f"Figure1_{P}_coverage_10bp_bins.tsv")
    gname = M.gene_id; a, b = int(M.window_start), int(M.window_end); s = M.strand; chrom_label = M.chrom_label
    shade = (M.shade_start, M.shade_end); shc = M.shade_color; al = float(M.shade_alpha)
    fig, axs = panel_figure(P)
    for ax in axs:
        ax.axvspan(*shade, color=shc, alpha=al, lw=0, zorder=0); ax.set_xlim((b, a) if s == '-' else (a, b)); ax.set_xticks([]); ax.set_yticks([])
        for sp in ['left', 'bottom', 'top', 'right']: ax.spines[sp].set_visible(False)
    # PlasmoDB v71 gene models (focal gene opaque and labelled; neighbours at 40 %)
    ax = axs[0]; ax.set_ylim(-1.0, 1.5)
    for gn_, sb in gff.groupby('gene'):
        foc = gn_ == gname; alp = 1 if foc else 0.4; ex_ = sb[sb.type == 'exon']
        if len(ex_) > 1: ax.plot([ex_.start.min(), ex_.end.max()], [0, 0], color=COL['cds'], lw=0.5, alpha=alp)
        for _, f in sb[sb.type != 'exon'].iterrows():
            h = 0.9 if f.type == 'CDS' else 0.5; ax.add_patch(Rectangle((f.start, -h / 2), f.end - f.start, h, fc=COL['cds'] if f.type == 'CDS' else COL['utr'], ec='none', alpha=alp))
        if foc: ax.text((max(sb.start.min(), a) + min(sb.end.max(), b)) / 2, 0.85, gn_, ha='center', fontsize=5.4, fontweight='bold')
    # stage-confirmed full-CDS long-read models
    ax = axs[1]
    for k, t_ in enumerate(LRm.itertuples()):
        ax.plot([t_.start, t_.end], [-k, -k], color='#555555', lw=0.4)
        for bb in t_.blocks.split(','): ax.plot([int(v) for v in bb.split('-')], [-k, -k], color='#555555', lw=2.0, solid_capstyle='butt')
    ax.set_ylim(-max(len(LRm), 5) + 0.3, 0.8)
    # lncRNA units (filled = same strand as the gene, outlined = opposite strand)
    ax = axs[2]
    for k, (u0, u1, us) in enumerate(zip(UN.start, UN.end, UN.strand)): ax.add_patch(Rectangle((u0, -k * 0.9 - 0.3), u1 - u0, 0.6, fc=COL['lnc'] if us == s else 'none', ec=COL['lnc'], lw=0.5))
    ax.set_ylim(-max(len(UN), 2) * 0.9 - 0.1, 0.6)
    # stranded RNA-seq, log10(1 + CPM)
    ax = axs[3]; x = CV.position.values
    ds = np.log10(1 + CV.daft_gene_strand_cpm.values); do = np.log10(1 + CV.daft_opposite_strand_cpm.values)
    ax.fill_between(x, 0, ds, color=COL['daft'], lw=0); ax.fill_between(x, 0, -do, color='#b8c7e0', lw=0); ax.axhline(0, color='#999999', lw=0.3)
    ym = max(ds.max(), do.max()) * 1.05; ax.set_ylim(-ym, ym)
    ax.text(0.01, 0.93, "gene strand", transform=ax.transAxes, fontsize=4.8, va='top', color=COL['daft']); ax.text(0.01, 0.07, "opposite", transform=ax.transAxes, fontsize=4.8, va='bottom', color='#8aa0c8')
    # CAGE (CTSS CPM, gene strand; sum of the ring, trophozoite and schizont tracks as in the approved script)
    ax = axs[4]; cg = CV.cage_ctss_cpm_sum_3_stages.values
    ax.fill_between(x, 0, cg, color=COL['cage'], lw=0.6); ax.plot(x, cg, color=COL['cage'], lw=0.5); ax.set_ylim(0, cg.max() * 1.1); ax.axhline(0, color='#999999', lw=0.3)
    ax.text(0.99, 0.9, f"max {cg.max():.0f} CPM", transform=ax.transAxes, fontsize=4.8, ha='right', va='top', color=COL['cage'])
    sb0 = b - 150 if s == '-' else a + 150; sb1 = sb0 - 500 if s == '-' else sb0 + 500
    ax.plot([sb0, sb1], [-0.18 * cg.max()] * 2, color=INK, lw=0.8, clip_on=False); ax.text((sb0 + sb1) / 2, -0.3 * cg.max(), "500 bp", fontsize=5, ha='center', va='top', clip_on=False)
    ax.text(0.99, -0.2, f"{chrom_label}:{a:,}–{b:,} ({s})", transform=ax.transAxes, fontsize=4.8, ha='right', va='top', color='#666666')
    axs[0].text(0, 1.45, M.title, transform=axs[0].transAxes, fontsize=6.5, va='bottom')
    save_panel(fig, P)
