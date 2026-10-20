from pfconfig import FIGURES_DIR, result  # paths: config/paths.yaml
import pandas as pd, numpy as np
import os
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'common.py')).read())
FIG = FIGURES_DIR
rd = lambda v: pd.read_csv(result(v), sep='\t', compression='gzip' if result(v).endswith('.gz') else 'infer')
rf = lambda f: pd.read_csv(f"{FIG}/{f}", sep='\t')

# common metric descriptions
MET = {'n': 'Genes evaluated (genes with a truth position).', 'n_truth': 'Genes with a truth position.', 'n_genes': 'Genes evaluated.', 'acc200': 'Fraction of genes whose predicted end is within 200 bp of the truth (genes without a prediction count as misses).',
       'ci_lo': 'Bootstrap 95% CI lower bound (2,000 gene resamples, seed 0).', 'ci_hi': 'Bootstrap 95% CI upper bound.', 'lo': 'Bootstrap 95% CI lower bound.', 'hi': 'Bootstrap 95% CI upper bound.',
       'median_abs_err': 'Median |predicted − truth| (bp) among genes with a prediction.', 'median_err': 'Median |predicted − truth| (bp) among genes with a prediction.', 'coverage': 'Fraction of evaluated genes with a prediction.',
       'acc@50': 'Fraction within 50 bp.', 'acc@100': 'Fraction within 100 bp.', 'acc@200': 'Fraction within 200 bp.', 'acc@500': 'Fraction within 500 bp.', 'key': 'Method key.', 'method': 'Method.',
       'truth': 'Truth set.', 'end': 'TSS or TTS.', 'tolerance_bp': 'Distance tolerance (bp).', 'acc': 'Fraction of genes within the tolerance.', 'group': 'Gene group (long-read status).'}
MKEY = ('Methods: Ours_chromatin = final model (V3_top8); Ours_nochromatin = same model without chromatin inputs; FixedOffset = fixed UTR length from the ATG/stop; '
        'v71all = PlasmoDB v71 as annotated, genes without a UTR take the ATG/stop (used in Figure 4A/B and the manuscript); PlasmoDB_v71 = v71 UTR ends only (genes without a v71 UTR count as no prediction); '
        'ChromHMM = in-house 10-state ChromHMM promoter state on the model chromatin tracks; ChromatinPeak = strongest chromatin peak upstream of the ATG; StringTie_short / StringTie_long = StringTie transcript ends from short-read / long-read assemblies; '
        'LongRead_v3 = long-read reference ends; Exp_* = the other experimental 5′-end sets used as predictors (truth-vs-truth agreement).')

# ======================= S5 =======================
s5 = []
A = rf('Figure4/data/Figure4_A_tss_consensus_acc200.tsv'); Bt = rf('Figure4/data/Figure4_B_tts_siegel_acc200.tsv')
assert abs(A.set_index('key').loc['Ours_chromatin', 'acc200'] - 0.7544) < 1e-3 and abs(Bt.set_index('key').loc['Ours_chromatin', 'acc200'] - 0.7589) < 1e-3
s5.append(dict(name='TSS_consensus_Fig4A', df=A, desc=MET, source='Figure4/data/Figure4_A_tss_consensus_acc200.tsv (approved table Figure4_A_tss_consensus_acc200.tsv)', note='1,649 genes; v71 = v71all'))
s5.append(dict(name='TTS_Siegel_Fig4B', df=Bt, desc=MET, source='Figure4/data/Figure4_B_tts_siegel_acc200.tsv (approved table Figure4_B_tts_siegel_acc200.tsv)', note='3,439 genes; v71 = v71all'))
T = rd('Figure4_C_consensus_tolerance.tsv'); s5.append(dict(name='TSS_consensus_tolerance', df=T, desc=MET, source='Figure4_C_consensus_tolerance.tsv (Figure4_C_consensus_tolerance.tsv)', note='25–500 bp; v71 here = PlasmoDB_v71 (UTR-only) variant'))
T = rd('Figure4_H_tts_siegel_tolerance.tsv'); s5.append(dict(name='TTS_Siegel_tolerance', df=T, desc=MET, source='Figure4_H_tts_siegel_tolerance.tsv (Figure4_H_tts_siegel_tolerance.tsv)', note='25–500 bp; v71 here = PlasmoDB_v71 (UTR-only) variant'))
T = rd('Figure4_D_per_truth_set_acc200.tsv'); s5.append(dict(name='TSS_per_truth_set', df=T, desc=MET, source='Figure4_D_per_truth_set_acc200.tsv (Figure4_D_per_truth_set_acc200.tsv)', note='Adjalley / Chappell / Shaw separately; v71 = UTR-only variant'))
T = rd('benchmark_results_all.tsv'); s5.append(dict(name='All_truths_groups_tols', df=T, desc=MET, source='benchmark_results_all.tsv (benchmark_results_all.tsv)', note='every method × truth set (TSS: Adjalley, Chappell, Shaw, consensus, long-read; TTS: long-read) × gene group × 50/100/200/500 bp'))
T = rd('benchmark_acc200_bootstrap.tsv'); s5.append(dict(name='Bootstrap_CI_acc200', df=T, desc=MET, source='benchmark_acc200_bootstrap.tsv (benchmark_acc200_bootstrap.tsv)', note='bootstrap CIs for benchmark_results_all'))
T = rd('siegel_pas_tts_benchmark.tsv')
s5.append(dict(name='TTS_Siegel_by_subset', df=T, desc={**MET, 'subset': 'Gene subset.', 'median_err_major': 'Median |predicted − major PAS| (bp).', 'acc200_lo': 'Bootstrap 95% CI lower bound of acc@200.', 'acc200_hi': 'Upper bound.',
    'acc200_any_PAS': 'Fraction within 200 bp of any Siegel PAS cluster of the gene.', 'signed_median': 'Median signed error (bp; negative = upstream of the PAS).'}, source='siegel_pas_tts_benchmark.tsv (siegel_pas_tts_benchmark.tsv)', note='v71 = PlasmoDB v71 as in that analysis'))
T = rf('Figure4/data/Figure4_F_accuracy_by_expression.tsv')
s5.append(dict(name='Accuracy_by_expression_Fig4F', df=T, desc={**MET, 'stratum': 'DAFT-seq expression quintile (Q1 lowest).', 'DAFT_median': 'Median DAFT-seq max CPM of the stratum.', 'acc50': 'Fraction within 50 bp.',
    'has_call': 'Fraction of genes with a prediction.', 'group': 'Gene group: reference genes (long-read ends) / model-predicted genes / no prediction.'}, source='Figure4/data/Figure4_F_accuracy_by_expression.tsv (approved Figure4_F_accuracy_by_expression.tsv)', note='Figure 4F'))
T = rf('Figure3/data/Figure3_D_cumulative_curve.tsv')
s5.append(dict(name='LongRead_truth_cumulative', df=T, desc={'distance_bp': 'Distance to the nearest long-read end cluster of the same gene (bp).', 'TSS_fraction_within': 'Fraction of the 2,389 reference genes whose model TSS is within the distance.', 'TTS_fraction_within': 'Same for TTS.'},
    source='Figure3/data/Figure3_D_cumulative_curve.tsv', note='Figure 3D (out-of-fold, long-read truth)'))
T = rf('Figure3/data/Figure3_G_outcome_summary.tsv')
s5.append(dict(name='LongRead_truth_outcomes', df=T, desc={'end': 'TSS or TTS.', 'category': 'Outcome of the model call relative to the gene\'s long-read ends (main end / alternative end / longer / shorter / no call ...).', 'n_genes': 'Genes.', 'n_total': 'Reference genes.', 'fraction': 'Fraction.'},
    source='Figure3/data/Figure3_G_outcome_summary.tsv', note='Figure 3G'))
T = rd('hmm_state_tss_predictor_benchmark.tsv')
s5.append(dict(name='HMM_chromatin_states', df=T, desc={**MET, 'gene_group': 'Gene group.', 'segmentation': 'Segmentation (Brown 2026 11-state; in-house 10-state ChromHMM) or comparator.', 'states_or_description': 'States treated as TSS-associated.',
    'stage': 'Segmentation stage.', 'rule': 'Rule converting a state block to a TSS (proximal edge, midpoint, ...).', 'acc@200_CI95_lo': 'Bootstrap 95% CI lower bound.', 'acc@200_CI95_hi': 'Upper bound.', 'median_abs_err_bp': 'Median absolute error (bp).',
    'vs_Ours_chromatin_diff_acc@200': 'Difference in acc@200 vs our model.', 'vs_Ours_chromatin_mcnemar_p': 'McNemar p vs our model.', 'vs_inhouse_ChromHMM_baseline_diff_acc@200': 'Difference vs original ChromHMM baseline.',
    'vs_inhouse_ChromHMM_baseline_mcnemar_p': 'McNemar p vs ChromHMM baseline.', 'vs_Ours_chromatin_mcnemar_BH_q': 'BH-adjusted q.', 'vs_inhouse_ChromHMM_baseline_mcnemar_BH_q': 'BH-adjusted q.', 'note': 'Note.'},
    source='hmm_state_tss_predictor_benchmark.tsv (hmm_state_tss_predictor_benchmark.tsv)', note='consensus TSS truth; 30 Brown rule combinations are post hoc except the primary rule'))
T = rd('patseq_signal_at_tts_summary.tsv')
s5.append(dict(name='PATseq_TTS_signal', df=T, desc={'subset': 'Gene subset.', 'method': 'TTS source.', 'n': 'Genes.', 'median_pct': 'Median percentile of PAT-seq signal at the predicted TTS within the window.', 'frac_PAT_ge1': 'Fraction of predicted TTS with ≥1 PAT-seq CPM.',
    'median_PAT_cpm': 'Median PAT-seq CPM at the predicted TTS.', 'within100_of_peak': 'Fraction within 100 bp of the window PAT-seq peak.', 'within200_of_peak': 'Fraction within 200 bp.', 'median_dist_peak': 'Median distance to the PAT-seq peak (bp).'},
    source='patseq_signal_at_tts_summary.tsv (patseq_signal_at_tts_summary.tsv)', note='supplementary only; PAT-seq has strong internal priming'))
P = rd('benchmark_method_predictions.tsv'); TS = rd('benchmark_truth_sets.tsv')
SP = rd('siegel2014_pas_per_gene.tsv')[['gene_id', 'n_clusters', 'PAS_major']].rename(columns={'n_clusters': 'Siegel2014_n_PAS_clusters', 'PAS_major': 'Siegel2014_major_PAS__TTS'})
PG = TS.merge(SP, on='gene_id', how='left').merge(P, on='gene_id', how='left')
dpg = {'gene_id': 'PlasmoDB gene ID.', 'Siegel2014_n_PAS_clusters': 'Number of Siegel 2014 PAS clusters for the gene.', 'Siegel2014_major_PAS__TTS': 'Major Siegel 2014 PAS (0-based; source table is 1-based, converted).'}
for c in PG.columns:
    if c in dpg: continue
    k, e = c.split('__'); dpg[c] = (f'Truth {e} (0-based) from {k}.' if c in TS.columns else f'{e} predicted by {k} (0-based).')
s5.append(dict(name='Per_gene_predictions', df=PG, desc=dpg, source='benchmark_truth_sets.tsv (benchmark_truth_sets.tsv); benchmark_method_predictions.tsv (benchmark_method_predictions.tsv); siegel2014_pas_per_gene.tsv (siegel2014_pas_per_gene.tsv)', note='5,285 genes; inputs to all S5 sheets'))
intro5 = ['Benchmark of TSS/TTS prediction methods against independent 5′-end and poly(A) data and against the long-read reference.',
          'Primary truth sets: TSS = consensus of ≥2 of Adjalley 2016 cap-seq, Chappell 2020 5′UTR-seq and Shaw 2021 5CAPture-seq within 100 bp (1,649 genes); TTS = Siegel 2014 major poly(A) site (3,439 genes). Per-dataset TSS sets: Adjalley 4,096, Chappell 3,029, Shaw 1,806 genes.',
          MKEY,
          'IMPORTANT: two PlasmoDB v71 variants exist in the saved results. Figure 4A/B and the manuscript (TSS 0.55, TTS 0.54) use v71all (ATG/stop fallback; TSS 0.5476, TTS 0.5406). The tolerance curves, per-truth-set and all-truth tables use PlasmoDB_v71 (UTR-only; TSS 0.5397, TTS 0.5379 at 200 bp). Harmonise before submission.',
          'Model predictions are out-of-fold (each chromosome predicted by the fold that did not see it in training).']
p5, m5 = write_table('TableS5_benchmarks.xlsx', 'Table S5. Benchmarks of boundary prediction methods', intro5, s5)
print(p5, m5)

# ======================= S6 =======================
s6 = []
tasks = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'evaluation', 'ablation', 'tasks.tsv'), sep='\t', header=None, names=['model', 'input_groups', 'chromatin_tracks', 'fold', 'n_runs'])
GRP = {'chrom': 'chromatin (24 ch)', 'onehot': 'DNA one-hot (4)', 'glm': 'Caduceus gLM (32)', 'gro': 'GRO-seq (3)', 'poly': 'polysome (3)', 'strand': 'inferred stranded RNA (18)', 'daft': 'DAFT-seq (6)'}
tasks['input_groups_readable'] = tasks.input_groups.map(lambda s: '; '.join(GRP.get(x, x) for x in s.split(',')))
tasks['fold_test_chromosomes'] = tasks.fold.map({'A': '1–5', 'B': '6–10', 'C': '11–14'})
cfg = tasks.groupby(['model', 'input_groups', 'input_groups_readable', 'chromatin_tracks']).agg(folds=('fold', lambda x: ','.join(sorted(set(x)))), n_runs_per_fold=('n_runs', 'max')).reset_index()
cfg['note'] = np.where(cfg.model.isin(['A_full']), 'fold A uses the original V3_top8 runs (not re-trained)', np.where(cfg.model.str.startswith('TS_'), 'trained on fold A only', ''))
ABL = ('Settings: A_full = full model (all 7 input groups, top-8 chromatin tracks "T8" = H2A.Z, H2B.Z, H3K9ac, H3K27ac, H3K4me3, ATAC, MNase, Pol II Ser2P); A_full2 = full model re-trained in the single-track batch (reference for TL/TA comparisons); '
       'L_no* = leave one input group out; S_* = single input class only; N_noshortread = without DAFT-seq and inferred stranded RNA; TL_* = full model leaving one of the 8 chromatin tracks out; TA_* = full model adding one extra track; TS_* = chromatin-only model with a single track (fold A only).')
dcfg = {'model': 'Setting name.', 'input_groups': 'Input groups (code names).', 'input_groups_readable': 'Input groups with channel counts.', 'chromatin_tracks': 'Chromatin tracks used (T8 = the 8 tracks of the final model; "-" = none).',
        'folds': 'Chromosome folds trained (A: test chr 1–5; B: 6–10; C: 11–14).', 'n_runs_per_fold': 'Training runs per fold combined into the evaluated ensemble (3 = base seed 0 + seed 1 + dropout 0.3/wd 0.01 variant, per abl_train.sbatch; 2 and 1 as listed, to verify which runs).', 'note': 'Note.'}
s6.append(dict(name='Settings', df=cfg, desc=dcfg, source='evaluation/ablation/tasks.tsv and evaluation/ablation/abl_train.sbatch', note=ABL))
s6.append(dict(name='Fold_definition', df=rf('Figure3/data/Figure3_B_chromosome_folds.tsv'), desc={'chromosome': 'Chromosome.', 'length_bp': 'Length (bp).', 'fold1_role': 'Role in fold 1 (A).', 'fold2_role': 'Role in fold 2 (B).', 'fold3_role': 'Role in fold 3 (C).'}, source='Figure3/data/Figure3_B_chromosome_folds.tsv', note='Figure 3B'))
s6.append(dict(name='Fold_gene_counts', df=rf('Figure3/data/Figure3_B_test_gene_counts.tsv'), desc={'fold': 'Fold.', 'test_chromosomes': 'Test chromosomes.', 'validation_chromosomes': 'Validation chromosomes.', 'n_reference_genes': 'Reference genes on the test chromosomes.'}, source='Figure3/data/Figure3_B_test_gene_counts.tsv', note=''))
DEL = {'truth': 'Truth set.', 'model': 'Setting (see Settings).', 'label': 'Input removed / kept.', 'n': 'Genes.', 'acc_full': 'acc@200 of the full model (A_full).', 'acc': 'acc@200 of the setting.', 'delta': 'acc − acc_full.',
       'lo': 'Bootstrap 95% CI lower bound of delta.', 'hi': 'Upper bound.', 'better': 'Genes correct only with the setting.', 'worse': 'Genes correct only with the full model.', 'p_sign': 'Two-sided exact sign (McNemar) test p.'}
s6.append(dict(name='Fig3E_longread_truth', df=rf('Figure3/data/Figure3_E_input_ablation.tsv'), desc=DEL, source='Figure3/data/Figure3_E_input_ablation.tsv (= fig3E_ablation_longread.tsv fig3E_ablation_longread.tsv)', note='pooled out-of-fold (folds A+B+C)'))
DIC = {'model': 'Setting.', 'folds': 'Folds pooled (ABC = all three, out-of-fold over chr 1–14; A = fold A only, chr 1–5).', 'truth': 'Truth set.', 'end': 'TSS or TTS.', 'subset': 'Gene subset: all / low_expression_Q1 / heterochromatic.',
       'n': 'Genes.', 'acc200': 'Fraction within 200 bp.', 'lo': 'Bootstrap 95% CI lower bound.', 'hi': 'Upper bound.', 'coverage': 'Fraction with a prediction.', 'delta_vs_ref': 'acc200 − reference setting (A_full; A_full2 for single-track settings).',
       'better': 'Genes correct only with the setting.', 'worse': 'Genes correct only with the reference.', 'p_vs_ref': 'Exact sign (McNemar) test p vs reference.', 'q_BH': 'BH-adjusted q.'}
s6.append(dict(name='Input_classes_all_truths', df=rd('ablation_cv_input_classes.tsv'), desc=DIC, source='ablation_cv_input_classes.tsv (ablation_cv_input_classes.tsv)', note='independent (consensus TSS, Siegel PAS) and long-read truths; e.g. L_nochrom TSS consensus delta −0.062'))
s6.append(dict(name='Single_tracks', df=rd('ablation_single_track.tsv'), desc=DIC, source='ablation_single_track.tsv (ablation_single_track.tsv)', note='TL/TA (folds ABC) and TS (fold A) settings'))
PMd = {'model': 'Setting.', 'folds': 'Folds pooled (ABC or A).', 'end_type': 'TSS or TTS.', 'scheme': 'Matching scheme: one_to_one (each truth end matched to ≤1 call) or merged_region.', 'n_chrom': 'Chromosomes evaluated.', 'n_calls': 'Number of calls.'}
for t in (50, 100, 200, 500): PMd.update({f'P@{t}': f'Precision at {t} bp (calls vs all long-read end labels).', f'R@{t}': f'Recall at {t} bp.', f'F1@{t}': f'F1 at {t} bp.'})
s6.append(dict(name='Point_metrics_PRF', df=rd('ablation_cv_point_metrics.tsv'), desc=PMd, source='ablation_cv_point_metrics.tsv (ablation_cv_point_metrics.tsv; identical to cv_point_metrics.tsv cv_point_metrics.tsv)', note='genome-wide point calls, long-read labels'))
CE = rf('Figure3/data/Figure3_F_chromatin_by_expression.tsv')
s6.append(dict(name='Chromatin_by_expression', df=CE, desc={'truth': 'Truth set.', 'stratum': 'DAFT-seq expression quintile.', 'n': 'Genes.', 'daft_median': 'Median DAFT max CPM.', 'acc_A_full': 'acc@200, full model.', 'lo_A_full': 'CI lower.', 'hi_A_full': 'CI upper.',
    'acc_L_nochrom': 'acc@200, no-chromatin model.', 'lo_L_nochrom': 'CI lower.', 'hi_L_nochrom': 'CI upper.', 'delta': 'Difference.', 'better': 'Genes correct only without chromatin.', 'worse': 'Genes correct only with chromatin.', 'p_sign': 'Sign test p.'},
    source='Figure3/data/Figure3_F_chromatin_by_expression.tsv (= fig3F_chromatin_by_expression_longread.tsv fig3F_chromatin_by_expression_longread.tsv)', note='Figure 3F'))
s6.append(dict(name='Thresholds_lncRNA_by_fold', df=rd('ablation_cv_lncRNA_strand.tsv'), desc={'model': 'Setting.', 'folds': 'Folds pooled.', 'lncRNA_TSS_called': 'Fraction of long-read lncRNA TSSs receiving a TSS call.', 'n_lncRNA': 'lncRNA TSSs evaluated.',
    'mirror_TSS_opposite': 'Fraction of reference TSSs with a spurious call on the opposite strand at the mirrored position.', 'n_lr': 'Reference genes evaluated.',
    'thr_A_TSS': 'TSS probability threshold selected on the validation chromosomes of fold A.', 'thr_A_TTS': 'TTS threshold, fold A.', 'thr_B_TSS': 'TSS threshold, fold B.', 'thr_B_TTS': 'TTS threshold, fold B.', 'thr_C_TSS': 'TSS threshold, fold C.', 'thr_C_TTS': 'TTS threshold, fold C.'},
    source='ablation_cv_lncRNA_strand.tsv (ablation_cv_lncRNA_strand.tsv; identical to cv_extra.tsv cv_extra.tsv)', note='per-fold calling thresholds; column meanings of lncRNA/mirror metrics to verify'))
intro6 = ['Ablation of input classes and individual chromatin tracks. Every setting was re-trained from scratch with the same chromosome folds (A: test chr 1–5, B: 6–10, C: 11–14); evaluation is on pooled out-of-fold predictions unless the folds column says A.',
          ABL,
          'MISSING: per-seed metrics and per-fold metrics for the pooled (ABC) settings were not saved as tables; each fold was evaluated as an ensemble of its runs (Settings, n_runs_per_fold). Per-seed checkpoints exist (v3_model_code_weights.tar.gz v3_model_code_weights.tar.gz for the final model; ablation runs on HPCC) but per-seed scores would require re-evaluation.',
          'Excluded: abl_point_metrics.tsv / ablation_point_metrics_longread.tsv / ablation_independent_truth.tsv / abl_independent_truth.tsv (earlier 360-gene / fold-A evaluation round whose fold provenance is not recorded).']
p6, m6 = write_table('TableS6_ablation.xlsx', 'Table S6. Ablation of model inputs', intro6, s6, coord=False)
print(p6, m6)
