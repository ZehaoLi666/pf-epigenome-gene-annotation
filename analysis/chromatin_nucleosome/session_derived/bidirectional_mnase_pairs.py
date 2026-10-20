"""Session-derived computation of bidirectional_mnase_pairs.tsv.

Recovered from the record of the interactive analysis session in which bidirectional_mnase_pairs.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_boundaries_v3_all_genes.tsv, gene_daft_expression_truth.tsv.gz, bidirectional_pairs.tsv, divergent_pairs_coexpression_classes.tsv, MNase-seq_T5_qNorm.bw, MNase-seq_T10_qNorm.bw, MNase-seq_T15_qNorm.bw, MNase-seq_T20_qNorm.bw, MNase-seq_T25_qNorm.bw, MNase-seq_T30_qNorm.bw, MNase-seq_T35_qNorm.bw, MNase-seq_T40_qNorm.bw, confirmed_genes_stage_TSS_switches.tsv, chappell2020_5UTRseq_tagclusters.tsv, fig1_neighbour_pairs.png, bioreport_altTSS_families.png, fig_isoform2_alt_tss_prevalence.png, refined_utr_annotation.gff3, daft_gene_expression_v68.tsv.gz, stage_alt_tss_all_candidates.tsv, hc_reference.bed.gz, gene_transcript_end_catalogue.tsv, fig_isoform1_stage_switch_funnel.png, ablation_cv_input_classes.tsv, alt_tss_stage_usage_fullCDS.tsv, alt_tss_stage_usage_allmodels.tsv, fig2_te_vs_utr.png, Figure4_legend.md, divergent_pairs_mnase_ndr.tsv, report_stage_specific_alt_tss.md, fig_tss_predictor_accuracy.png, fig1_tss_metagene_apiap2_atac.png, fig_ndr_edge.png, fig3_model_comparison.png, fig2_termination_metagene.png, confirmed_alternative_TSSs.tsv, tss_cluster_positions.tsv, fig2_motif_position_tss_context.png, tss_cluster_cage.tsv, fig3_bidirectional_promoters.png, benchmark_method_predictions.tsv, gene2go_propagated.tsv, tss_clusters.tsv.gz, gene_vs_readend_v48.tsv, report_promoter_architecture.md, gro_flank_per_gene.tsv, transcript_stage_confirmation.tsv.gz, biology_analysis_report.md, fig1_uorf_census.png, fig_nucdetective_plus1.png, shaw2021_5CAPture_5CN_clusters.tsv, report_isoforms_families.md, go_terms.tsv, fig_hmm_states_tss_contrast.png, fig3_pas_motifs.png, refined_utr_annotation.tsv, fig_isoform3_family_forest.png
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import numpy as np
import pandas as pd
import pyBigWig
from scipy import stats
import statsmodels.formula.api as smf

# Load data
gb = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep='\t')
dft = pd.read_csv(result('gene_daft_expression_truth.tsv.gz'), sep='\t')
dv = pd.read_csv(result('divergent_pairs_coexpression_classes.tsv'), sep='\t')
T = pd.read_csv(result('divergent_pairs_mnase_ndr.tsv'), sep='\t')

AID = {
    "T5": result('MNase-seq_T5_qNorm.bw'),
    "T10": result('MNase-seq_T10_qNorm.bw'),
    "T15": result('MNase-seq_T15_qNorm.bw'),
    "T20": result('MNase-seq_T20_qNorm.bw'),
    "T25": result('MNase-seq_T25_qNorm.bw'),
    "T30": result('MNase-seq_T30_qNorm.bw'),
    "T35": result('MNase-seq_T35_qNorm.bw'),
    "T40": result('MNase-seq_T40_qNorm.bw'),
}
BW = {k: v for k, v in AID.items()}

b = pyBigWig.open(BW["T5"])
CH = {c: l for c, l in b.chroms().items()}
b.close()

raw = {c: np.zeros(L_, np.float32) for c, L_ in CH.items()}
for tp, p in BW.items():
    b = pyBigWig.open(p)
    for c in CH:
        raw[c] += np.nan_to_num(np.array(b.values(c, 0, CH[c]), dtype=np.float32)) / 8
    b.close()


def smooth(v, w=21):
    cs = np.concatenate([[0], np.cumsum(v, dtype=np.float64)])
    h = w // 2
    n = len(v)
    lo = np.clip(np.arange(n) - h, 0, n)
    hi = np.clip(np.arange(n) + h + 1, 0, n)
    return ((cs[hi] - cs[lo]) / (hi - lo)).astype(np.float32)


MN = {c: smooth(raw[c]) for c in CH}
del raw
MED = float(np.median(np.concatenate(list(MN.values()))))


def call_ndr(sig, thr, minw=50, gap=20):
    out = {}
    for c in sig:
        b_ = np.concatenate([[False], sig[c] < thr, [False]]).astype(np.int8)
        d = np.diff(b_)
        s = np.where(d == 1)[0]
        e = np.where(d == -1)[0]
        if len(s):
            ks = [s[0]]
            ke = [e[0]]
            for a_, z in zip(s[1:], e[1:]):
                if a_ - ke[-1] <= gap:
                    ke[-1] = z
                else:
                    ks.append(a_)
                    ke.append(z)
            s, e = np.array(ks), np.array(ke)
            k = (e - s) >= minw
            s, e = s[k], e[k]
        out[c] = (s, e)
    return out


NDR = call_ndr(MN, 0.5 * MED)

# Merge T with dv to get TSS columns from dv
k = ['gene_left', 'gene_right']
ex = T[k + ['TSS_left_final', 'TSS_right_final']].rename(columns={'TSS_left_final': 'L_ex', 'TSS_right_final': 'R_ex'})
d2 = dv.merge(ex, on=k, how='left')
for s, e in [('TSS_left_final', 'L_ex'), ('TSS_right_final', 'R_ex')]:
    ok = d2[e].notna()
    d2.loc[ok, s] = d2.loc[ok, e]
d2 = d2.drop(columns=['L_ex', 'R_ex'], errors='ignore')
P = d2.copy()

# Add TSS_left_v71 and TSS_right_v71 from T
T_merge = T[['gene_left', 'gene_right', 'TSS_left_v71', 'TSS_right_v71', 'atg_left', 'atg_right',
             'n_NDR_interTSS', 'n_nuc_interTSS_allNucs', 'mnase_P',
             'ndr_ATAC_refinedRule_anyStage', 'ATAC_meanStages_interTSS_refinedRule']].copy()
# Also get columns needed from T
T_extra = T[['gene_left', 'gene_right']].copy()
for col in T.columns:
    if col not in ['gene_left', 'gene_right'] and col not in P.columns:
        T_extra[col] = T[col]
P = P.merge(T[['gene_left', 'gene_right', 'TSS_left_v71', 'TSS_right_v71', 'atg_left', 'atg_right',
               'n_NDR_interTSS', 'n_nuc_interTSS_allNucs', 'mnase_P',
               'ndr_ATAC_refinedRule_anyStage']], on=k, how='left')

P['dist_v71'] = P.TSS_right_v71 - P.TSS_left_v71

# Load ATAC
AP = f"{LAB_DIR}/lncRNA project tracks/accessibility_nucleosome/"
ATAC = {}
for st in ["Ring", "Trophozoite", "Schizont"]:
    b = pyBigWig.open(AP + f"ATAC_{st}_minusInput.bw")
    ch = b.chroms()
    for c in CH:
        v = np.nan_to_num(np.array(b.values(c, 0, CH[c]), dtype=np.float32)) if c in ch else np.zeros(CH[c], np.float32)
        ATAC[c] = ATAC.get(c, 0) + v / 3
    b.close()

mo = []
at = []
for r in P.itertuples():
    if not (pd.notna(r.TSS_left_final) and pd.notna(r.TSS_right_final)):
        mo.append(np.nan)
        at.append(np.nan)
        continue
    a, b_ = int(r.TSS_left_final), int(r.TSS_right_final)
    if b_ - a < 300:
        m = (a + b_) // 2
        a, b_ = m - 150, m + 150
    if a < 0 or b_ > CH.get(r.chrom, 0):
        mo.append(np.nan)
        at.append(np.nan)
        continue
    mo.append(MN[r.chrom][a:b_].mean() / MED)
    at.append(ATAC[r.chrom][a:b_].mean())
P['mn_occ_rel'] = mo
P['atac_mean'] = at

P = P[P.TSS_left_final.notna() & P.TSS_right_final.notna()].copy()


def exit_metric(c, x, atg, strand):
    s, e = NDR[c]
    if strand == '+':
        k = (e >= atg - 3000) & (e <= atg)
        ss, ee = s[k], e[k]
        edge = ee
        sg = 1
    else:
        k = (s >= atg) & (s <= atg + 3000)
        ss, ee = s[k], e[k]
        edge = ss - 1
        sg = -1
    if len(ss) == 0 or not np.isfinite(x):
        return np.nan, False, False
    dd = (x - edge) * sg
    j = np.argmin(np.abs(dd))
    inN = np.any((ss <= x) & (x < ee))
    return dd[j], (0 <= dd[j] <= 150), inN


def utr_spans_ndr(c, x, atg, strand):
    s, e = NDR[c]
    if not np.isfinite(x):
        return np.nan
    if strand == '+':
        return bool(np.any((s > x) & (e <= atg)))
    return bool(np.any((s > atg) & (e <= x)))


Pi = P.set_index(P.gene_left + '|' + P.gene_right)

rows = []
for r in P.itertuples():
    for side, strand in [('left', '-'), ('right', '+')]:
        atg = getattr(r, f'atg_{side}')
        for an, col in [('refined', f'TSS_{side}_final'), ('v71', f'TSS_{side}_v71')]:
            tss_val = getattr(r, col)
            d_, h_, i_ = exit_metric(r.chrom, tss_val, atg, strand)
            sp = utr_spans_ndr(r.chrom, tss_val, atg, strand)
            rows.append(dict(pair=f"{r.gene_left}|{r.gene_right}", side=side, anchor=an,
                             d=d_, hit=h_, inNDR=i_, spans=sp, tss=tss_val, atg=atg))
EX = pd.DataFrame(rows)

# Coexpression models
H = ["h0", "h8", "h16", "h24", "h32", "h40", "h48"]
dft_idx = dft.set_index('gene_id')
L_daft = np.log2(dft_idx[H] + 1)

bins = [-1e9, 200, 400, 600, 800, 1000, 1500, 2000, 3000, 1e9]
lab = ['≤200', '200–400', '400–600', '600–800', '0.8–1k', '1–1.5k', '1.5–2k', '2–3k', '>3k']
P['dbin'] = pd.cut(P.TSS_dist_final, bins, labels=lab)
P['ncls'] = np.where(P.n_NDR_interTSS <= 1, '1 NDR', np.where(P.n_NDR_interTSS >= 2, '≥2 NDRs', '?'))

SB = P.groupby('dbin', observed=False).agg(
    n=('r', 'size'),
    mn_shared=('mnase_P', lambda s: (s == 'shared').mean()),
    atac_shared=('ndr_ATAC_refinedRule_anyStage', lambda s: (s == 'shared').mean()),
    nNDR=('n_NDR_interTSS', 'median')
)

D = P[(P.TSS_dist_final > 0) & P.r.notna() & (P.min_DAFT_max > 0)].copy()
D['z'] = np.arctanh(D.r.clip(-0.999, 0.999))
D['ld'] = np.log10(D.TSS_dist_final)
D['le'] = np.log10(D.min_DAFT_max)
D['atac_z'] = (D.atac_mean - D.atac_mean.mean()) / D.atac_mean.std()
D['mn_z'] = (D.mn_occ_rel - D.mn_occ_rel.mean()) / D.mn_occ_rel.std()
D['two'] = (P.loc[D.index, 'n_NDR_interTSS'] >= 2).astype(int)

m2 = smf.ols('z~ld+le+two', D).fit(cov_type='HC3')
ci2 = m2.conf_int().loc['two'].values
m3 = smf.ols('z~ld+le+atac_z+mn_z', D).fit(cov_type='HC3')
ci3 = m3.conf_int()

# Save output
out = P[['gene_left', 'gene_right', 'chrom', 'TSS_left_final', 'TSS_right_final', 'TSS_dist_final',
         'TSS_left_v71', 'TSS_right_v71', 'dist_v71', 'r', 'min_DAFT_max', 'cat',
         'mnase_P', 'ndr_ATAC_refinedRule_anyStage', 'n_NDR_interTSS', 'n_nuc_interTSS_allNucs',
         'mn_occ_rel', 'atac_mean']].copy()
ex2 = EX.pivot_table(index='pair', columns=['side', 'anchor'], values=['hit', 'spans', 'd'], aggfunc='first')
ex2.columns = ['_'.join(map(str, c)) for c in ex2.columns]
out = out.merge(ex2, left_on=out.gene_left + '|' + out.gene_right, right_index=True, how='left').drop(columns='key_0', errors='ignore')
out.to_csv('bidirectional_mnase_pairs.tsv', sep='\t', index=False)
