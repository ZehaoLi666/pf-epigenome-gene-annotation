"""Gather / derive every input of Figure 3 and write per-panel source tables to ../data/.

This is the only script that touches upstream files.  It reproduces, line for line,
the data handling of the code that rendered the approved composite (make_fig3.py,
abl_panels.py and tracks.py).

Upstream inputs
---------------
* Intermediate result tables (names below), read from RESULTS_DIR (config/paths.yaml),
  or from a folder given with:  python prep_data.py --upstream /path/to/folder
* Local genome tracks (paths below; needs pyBigWig), only for panel C.

Panel scripts never read upstream files; they read only ../data/Figure3_*.tsv.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import os, sys, argparse
import numpy as np, pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.normpath(os.path.join(HERE, "..", "data"))

ARTIFACTS = {   # file names in RESULTS_DIR
    "gene_level_eval.tsv": "gene_level_eval.tsv",            # per-gene CV evaluation (2,389 reference genes)
    "example_confirmed_models.tsv": "example_confirmed_models.tsv",   # stage-confirmed long-read models, 3 example loci
    "fig3E_ablation_longread.tsv": "fig3E_ablation_longread.tsv",    # input ablation vs long-read ends
    "fig3F_chromatin_by_expression_longread.tsv": "fig3F_chromatin_by_expression_longread.tsv",  # full vs no-chromatin by DAFT quintile
}
LAB = LAB_DIR
IGV = f"{LAB}/IGV_longread_transcript_models"
GFF = f"{LAB}/PlasmoDB/PlasmoDB-68/PlasmoDB-68_Pfalciparum3D7.gff"
TRK = f"{LAB}/lncRNA project tracks"
CAGE = f"{TRK}/cage_tss"


def resolve(name, upstream):
    if upstream:
        return os.path.join(upstream, name)
    return result(ARTIFACTS[name])


# ---------------- track readers (verbatim from fig1/tracks.py) ----------------
def load_gff():
    g = pd.read_csv(GFF, sep="\t", comment="#", header=None, names="chrom src type start end score strand phase attr".split(), dtype={'start': int, 'end': int})
    g = g[g.type.isin(['protein_coding_gene', 'ncRNA_gene', 'mRNA', 'CDS', 'five_prime_UTR', 'three_prime_UTR', 'exon', 'gene'])].copy()
    g['ID'] = g.attr.str.extract(r'ID=([^;]+)')
    g['Parent'] = g.attr.str.extract(r'Parent=([^;]+)')
    g['gene'] = np.where(g.type.str.endswith('gene'), g.ID, None)
    g.loc[g.type == 'mRNA', 'gene'] = g.loc[g.type == 'mRNA', 'Parent']
    tx2g = dict(zip(g[g.type == 'mRNA'].ID, g[g.type == 'mRNA'].Parent))
    m = ~g.type.str.endswith('gene') & (g.type != 'mRNA')
    g.loc[m, 'gene'] = g.loc[m, 'Parent'].map(lambda p: tx2g.get(str(p).split(',')[0], p))
    return g


def bw_vals(path, chrom, a, b, nb):
    import pyBigWig
    bw = pyBigWig.open(path); v = np.array(bw.stats(chrom, a, b, type='mean', nBins=nb), dtype=float); bw.close()
    return np.nan_to_num(v)


def daft(chrom, a, b, strand, nb):
    return bw_vals(f"{IGV}/DAFTseq_pool_{'plus' if strand == '+' else 'minus'}_CPM_10bp.bw", chrom, a, b, nb)


def cage(chrom, a, b, strand, nb):
    s = 'pos' if strand == '+' else 'neg'
    return sum(bw_vals(f"{CAGE}/CTSS_{st}_{s}_CPM.bw", chrom, a, b, nb) for st in ['ring', 'trophozoite', 'schizont'])


def pars(x):
    return [float(v) for v in str(x).split(';')] if isinstance(x, str) and x else []


def w(df, name):
    df.to_csv(os.path.join(OUT, name), sep="\t", index=False)
    print(f"  {name}: {len(df)} rows")


def main(upstream=None):
    os.makedirs(OUT, exist_ok=True)
    C = pd.read_csv(resolve("gene_level_eval.tsv", upstream), sep="\t", index_col=0)
    C.index.name = 'gene_id'

    # ---------------- B: chromosome folds ----------------
    CS = [640851, 947102, 1067971, 1200490, 1343557, 1418242, 1445207, 1472805, 1541735, 1687656, 2038340, 2271494, 2925236, 3291936]
    folds = [("Fold 1", range(1, 6), [6, 7]), ("Fold 2", range(6, 11), [11]), ("Fold 3", range(11, 15), [10])]
    rows = []
    for k, L in enumerate(CS, 1):
        r = dict(chromosome=k, length_bp=L)
        for nm, te, va in folds:
            r[nm.replace(' ', '').lower() + '_role'] = 'test' if k in te else ('validation' if k in va else 'training')
        rows.append(r)
    w(pd.DataFrame(rows), "Figure3_B_chromosome_folds.tsv")
    chrn = C.chrom.str[6:8].astype(int)
    w(pd.DataFrame([dict(fold=nm, test_chromosomes=f"{min(te)}-{max(te)}", validation_chromosomes=",".join(map(str, va)),
                         n_reference_genes=int(chrn.isin(list(te)).sum())) for nm, te, va in folds]),
      "Figure3_B_test_gene_counts.tsv")

    # ---------------- C: example loci ----------------
    gff = load_gff()
    MOD = pd.read_csv(resolve("example_confirmed_models.tsv", upstream), sep="\t")
    MOD['blk'] = MOD.blocks.apply(lambda x: [tuple(map(int, b.split('-'))) for b in x.split(';')])
    EX = [("PF3D7_1302800", "Main TSS call matches the main long-read TSS"),
          ("PF3D7_0716900", "Main TSS call matches an alternative TSS"),
          ("PF3D7_1323700", "Main TSS call off; another call matches")]
    ex_rows, feat, mods, clus, calls, trk = [], [], [], [], [], []
    for j, (gname, ttl) in enumerate(EX, 1):
        r = C.loc[gname]; s = r.strand; sg = 1 if s == '+' else -1; chrom = r.chrom
        refT = [r.atg - sg * d for d in pars(r.ref_TSS_dists)]; refE = [r.stop + sg * d for d in pars(r.ref_TTS_dists)]
        alt = [float(v) for v in str(r.model_alt_TSS).split(',')] if isinstance(r.model_alt_TSS, str) else []
        altE = [float(v) for v in str(r.model_alt_TTS).split(',')] if isinstance(r.model_alt_TTS, str) else []
        mm = MOD[MOD.locus == gname]
        pts = refT + refE + [r.model_TSS, r.model_TTS, r.lr_TSS, r.lr_TTS] + alt + altE + list(mm.start) + list(mm.end)
        a = int(np.nanmin(pts)) - 350; b = int(np.nanmax(pts)) + 350; nb = int((b - a) / 10)
        ex_rows.append(dict(example=j, gene_id=gname, title=ttl, chrom=chrom, strand=s, window_start=a, window_end=b, n_bins=nb,
                            main_longread_TSS=r.lr_TSS, main_TSS_call=r.model_TSS, main_TSS_prob=r.model_TSS_prob, main_TTS_call=r.model_TTS))
        # v71 track: GFF features overlapping the window (file order preserved)
        gg = gff[(gff.chrom == chrom) & (gff.end > a) & (gff.start < b) & gff.type.isin(['CDS', 'five_prime_UTR', 'three_prime_UTR', 'exon'])]
        for _, f in gg.iterrows():
            feat.append(dict(example=j, gene=f.gene, type=f.type, start=f.start, end=f.end, focal_gene=f.gene == gname))
        # long-read models, in plotted order (row 0 at top)
        lr = mm.copy(); lr['p5'] = np.where(lr.strand == '+', lr.start, lr.end); lr['p3'] = np.where(lr.strand == '+', lr.end, lr.start)
        lr = lr.sort_values('p5', ascending=(s == '+'))
        tss_calls = [r.model_TSS] + alt; tts_calls = [r.model_TTS] + altE
        for k, (_, t_) in enumerate(lr.iterrows()):
            mods.append(dict(example=j, row=k, transcript_id=t_.tid, chrom=t_.chrom, start=t_.start, end=t_.end, strand=t_.strand,
                             blocks=t_.blocks, full_cds=bool(t_.full_cds), end5=t_.p5, end3=t_.p3,
                             end5_within_200bp_of_TSS_call=bool(min(abs(t_.p5 - q) for q in tss_calls) <= 200),
                             end3_within_200bp_of_TTS_call=bool(min(abs(t_.p3 - q) for q in tts_calls) <= 200)))
        for p in refT: clus.append(dict(example=j, end='TSS', position=p))
        for p in refE: clus.append(dict(example=j, end='TTS', position=p))
        sup = lambda q, ends, refs: min([abs(e - q) for e in ends] + [abs(e - q) for e in refs] + [1e9]) <= 200
        calls.append(dict(example=j, end='TSS', call='main', position=r.model_TSS, probability=r.model_TSS_prob, longread_end_within_200bp=np.nan))
        calls.append(dict(example=j, end='TTS', call='main', position=r.model_TTS, probability=r.model_TTS_prob, longread_end_within_200bp=np.nan))
        for q in alt: calls.append(dict(example=j, end='TSS', call='other', position=q, probability=np.nan, longread_end_within_200bp=sup(q, list(lr.p5), refT)))
        for q in altE: calls.append(dict(example=j, end='TTS', call='other', position=q, probability=np.nan, longread_end_within_200bp=sup(q, list(lr.p3), refE)))
        # signal tracks, nb bins of 10 bp across the window
        x = np.linspace(a, b, nb)
        T = dict(example=j, bin=np.arange(nb), position=x, DAFTseq_log10_1p_CPM=np.log10(1 + daft(chrom, a, b, s, nb)))
        for nm, pat in [("H2AZ", "histone_marks/H2AZ_{s}_minusInput.bw"), ("H3K9ac", "histone_marks/H3K9ac_{s}_minusInput.bw"),
                        ("H3K4me3", "histone_marks/H3K4me3_{s}_minusInput.bw"), ("ATAC", "accessibility_nucleosome/ATAC_{s}_minusInput.bw")]:
            v = np.mean([bw_vals(f"{TRK}/{pat.format(s=st)}", chrom, a, b, nb) for st in ('Ring', 'Trophozoite', 'Schizont')], axis=0)
            T[nm + "_minusInput_mean_RTS_clip0"] = np.clip(v, 0, None)
        T['CAGE_CPM'] = np.abs(cage(chrom, a, b, s, nb))
        trk.append(pd.DataFrame(T))
    w(pd.DataFrame(ex_rows), "Figure3_C_examples.tsv")
    w(pd.DataFrame(feat), "Figure3_C_v71_features.tsv")
    w(pd.DataFrame(mods), "Figure3_C_longread_models.tsv")
    w(pd.DataFrame(clus), "Figure3_C_longread_end_clusters.tsv")
    w(pd.DataFrame(calls), "Figure3_C_model_calls.tsv")
    w(pd.concat(trk, ignore_index=True), "Figure3_C_signal_tracks.tsv")

    # ---------------- D: distance to nearest long-read end ----------------
    def nearest(row, end):
        m = row['u5_m'] if end == 'TSS' else row['u3_m']; lr_ = row['u5_lr'] if end == 'TSS' else row['u3_lr']
        refs = pars(row['ref_TSS_dists'] if end == 'TSS' else row['ref_TTS_dists']) + [lr_]
        return np.nan if pd.isna(m) else min(abs(m - q) for q in refs)
    NE = {end: C.apply(lambda r_: nearest(r_, end), axis=1).values for end in ('TSS', 'TTS')}
    w(pd.DataFrame({'gene_id': C.index, 'TSS_nearest_lr_end_bp': NE['TSS'], 'TTS_nearest_lr_end_bp': NE['TTS']}), "Figure3_D_nearest_end.tsv")
    xg = np.arange(0, 1001, 5); cum = {'distance_bp': xg}
    for end in ('TSS', 'TTS'):
        d = np.where(np.isnan(NE[end]), np.inf, NE[end]); cum[f'{end}_fraction_within'] = [(d <= v).sum() / len(d) for v in xg]
    w(pd.DataFrame(cum), "Figure3_D_cumulative_curve.tsv")

    # ---------------- E: input ablation ----------------
    AB3 = pd.read_csv(resolve("fig3E_ablation_longread.tsv", upstream), sep="\t")
    w(AB3, "Figure3_E_input_ablation.tsv")

    # ---------------- F: chromatin effect by expression quintile ----------------
    TF3 = pd.read_csv(resolve("fig3F_chromatin_by_expression_longread.tsv", upstream), sep="\t")
    w(TF3, "Figure3_F_chromatin_by_expression.tsv")

    # ---------------- G: outcome categories ----------------
    w(C[['cat_TSS', 'cat_TTS']].reset_index(), "Figure3_G_outcome_by_gene.tsv")
    cats = ['main end', 'alternative end', 'longer', 'shorter', 'no call']
    w(pd.DataFrame([dict(end=end, category=k, n_genes=int((C['cat_' + end] == k).sum()), n_total=len(C),
                         fraction=(C['cat_' + end] == k).sum() / len(C)) for end in ('TSS', 'TTS') for k in cats]),
      "Figure3_G_outcome_summary.tsv")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--upstream", default=None)
    args, _ = ap.parse_known_args()
    main(args.upstream)
