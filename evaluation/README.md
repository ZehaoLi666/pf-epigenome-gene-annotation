# Evaluation

## benchmarks/
Independent truth sets (5′-end consensus of Adjalley 2016 cap-seq, Chappell 2020 5′UTR-seq and Shaw 2021
5CAPture; Siegel 2014 poly(A) sites; long-read ends; PAT-seq) and alternative methods.

- `run_stringtie.sh`, `stringtie_to_genes.py` — StringTie baselines (short reads, long reads), transcripts assigned to genes.
- `chromhmm_10state/` — ChromHMM 10-state model (emissions / transitions) used for the ChromHMM promoter baseline.
- `bench_eval.py` — gene-level accuracy at fixed tolerances for a methods × truths table.
- `pat_anchor.sh` — PAT-seq cleavage sites from poly(T)-anchored reads (HPCC).
- `session_derived/` — derivations of `benchmark_truth_sets.tsv`, `benchmark_method_predictions.tsv`,
  `benchmark_results_all.tsv`, `benchmark_acc200_bootstrap.tsv`, the Siegel / PAT-seq / 5′-end truth tables and the
  tolerance and per-truth-set tables used in Table S5.

## ablation/
Three-fold chromosome cross-validation of input-class and single-track ablations (Figure 3E/F, Table S6), run as
SLURM arrays on the HPCC.

- `tasks.tsv` — every setting (model name, input groups, chromatin tracks, fold, number of runs); `configs.tsv` — input groups.
- `env.sh`, `abl_train.sbatch` (training array), `pred_prob.py`, `abl_eval.sbatch`, `abl_eval.py` (fold-wise
  calling and scoring), `repro.sbatch`, `repro2.sbatch` (reproduction of the reference model).
- `score_local.py`, `cv_score.py` — scoring of per-gene CV calls against independent truths and long-read ends,
  overall and in subsets (lowest expression quintile, heterochromatic genes), paired sign tests.
- `session_derived/` — derivations of the ablation summary tables.
