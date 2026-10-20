"""Session-derived computation of divergent_pairs_coexpression_classes.tsv.

Recovered from the record of the interactive analysis session in which divergent_pairs_coexpression_classes.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_boundaries_v3_all_genes.tsv, gene_daft_expression_truth.tsv.gz, bidirectional_pairs.tsv, divergent_pairs_mnase_ndr.tsv, confirmed_genes_stage_TSS_switches.tsv, chappell2020_5UTRseq_tagclusters.tsv, fig1_neighbour_pairs.png, bioreport_altTSS_families.png, fig_isoform2_alt_tss_prevalence.png, refined_utr_annotation.gff3, daft_gene_expression_v68.tsv.gz, stage_alt_tss_all_candidates.tsv, hc_reference.bed.gz, gene_transcript_end_catalogue.tsv, fig_isoform1_stage_switch_funnel.png, ablation_cv_input_classes.tsv, alt_tss_stage_usage_fullCDS.tsv, alt_tss_stage_usage_allmodels.tsv, fig2_te_vs_utr.png, Figure4_legend.md, report_stage_specific_alt_tss.md, fig_tss_predictor_accuracy.png, fig1_tss_metagene_apiap2_atac.png, fig_ndr_edge.png, fig3_model_comparison.png, fig2_termination_metagene.png, confirmed_alternative_TSSs.tsv, tss_cluster_positions.tsv, fig2_motif_position_tss_context.png, tss_cluster_cage.tsv, fig3_bidirectional_promoters.png, benchmark_method_predictions.tsv, gene2go_propagated.tsv, tss_clusters.tsv.gz, gene_vs_readend_v48.tsv, report_promoter_architecture.md, gro_flank_per_gene.tsv, transcript_stage_confirmation.tsv.gz, biology_analysis_report.md, fig1_uorf_census.png, fig_nucdetective_plus1.png, shaw2021_5CAPture_5CN_clusters.tsv, report_isoforms_families.md, go_terms.tsv, fig_hmm_states_tss_contrast.png, fig3_pas_motifs.png, refined_utr_annotation.tsv, fig_isoform3_family_forest.png
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import pandas as pd
import numpy as np

# Load data
dft = pd.read_csv(result('gene_daft_expression_truth.tsv.gz'), sep='\t').set_index('gene_id')
H = ['h0', 'h8', 'h16', 'h24', 'h32', 'h40', 'h48']
L = np.log2(dft[H] + 1)

bp = pd.read_csv(result('bidirectional_pairs.tsv'), sep='\t')

# Build divergent pairs subset
dv = bp[bp.orientation == 'divergent'].dropna(subset=['pearson_DAFT_log2']).copy()

# Compute null distribution for category thresholds
rng = np.random.default_rng(0)
ids = L.index[L.var(axis=1) > 0].values
gb = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep='\t').set_index('gene_id')
a = rng.choice(ids, 20000)
b = rng.choice(ids, 20000)
keep = gb.loc[a, 'chrom'].values != gb.loc[b, 'chrom'].values

def rows_r(x, y):
    X = L.loc[x].values
    Y = L.loc[y].values
    X = X - X.mean(1, keepdims=True)
    Y = Y - Y.mean(1, keepdims=True)
    return (X * Y).sum(1) / np.sqrt((X**2).sum(1) * (Y**2).sum(1))

rnull = rows_r(a[keep], b[keep])

bins = [-1.01, -0.3, 0.3, 0.7, 1.01]
labs = ['anti (r≤−0.3)', 'uncorrelated (−0.3–0.3)', 'moderate (0.3–0.7)', 'strong (≥0.7)']

ph_t = np.array([0, 8, 16, 24, 32, 40])

def phase(g):
    v = L.loc[g, H[:-1]].values
    return ph_t[np.argmax(v)]

def amp(g):
    v = L.loc[g, H].values
    return v.max() - v.min()

dv['phase_L'] = dv.gene_left.map(phase)
dv['phase_R'] = dv.gene_right.map(phase)
d_ = np.abs(dv.phase_L - dv.phase_R)
dv['phase_diff_h'] = np.minimum(d_, 48 - d_)
dv['min_amp'] = [min(amp(x), amp(y)) for x, y in zip(dv.gene_left, dv.gene_right)]
dv['min_expr'] = dv.min_DAFT_max

het = r'pfemp1|erythrocyte membrane protein 1|\brifin\b|stevor|surfin|pfmc-2tm'
desc = gb.desc.fillna('').str.replace('%2C', ',').str.lower()
hetset = set(desc[desc.str.contains(het)].index)
dv['any_het'] = dv.gene_left.isin(hetset) | dv.gene_right.isin(hetset)
dv['ndr'] = dv.NDR_refinedRule_anyStage
dv['cat'] = pd.cut(dv.pearson_DAFT_log2, bins, labels=labs)

# Merge updated TSS coordinates from MNase NDR file
m = pd.read_csv(result('divergent_pairs_mnase_ndr.tsv'), sep='\t')
k = ['gene_left', 'gene_right']
ex = m[k + ['TSS_left_final', 'TSS_right_final']].rename(columns={'TSS_left_final': 'L_ex', 'TSS_right_final': 'R_ex'})
d2 = dv.merge(ex, on=k, how='left')
for s, e in [('TSS_left_final', 'L_ex'), ('TSS_right_final', 'R_ex')]:
    ok = d2[e].notna()
    d2.loc[ok, s] = d2.loc[ok, e]
d2 = d2.drop(columns=['L_ex', 'R_ex', 'TSS_left_refinedRule', 'TSS_right_refinedRule'])

d2.to_csv('divergent_pairs_coexpression_classes.tsv', sep='\t', index=False)
