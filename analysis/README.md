# Downstream analyses

| Folder | Scripts | Paper |
|---|---|---|
| `chromatin_nucleosome/` | `common.py` (shared inputs), `part1.py`, `part1c.py` (HMM chromatin states at TSS/TTS; state-based TSS predictor benchmark), `part2.py`, `part2b.py`, `part2c.py`, `part2d.py` (nucleosome positions, meta-profiles and dyad density around refined vs v71 TSSs, nucDetective calls), `ndr_edge.py` (NDR edges vs TSS), `divergent_mnase.py` (MNase NDR in head-to-head pairs), `input_coverage.py` (input signal percentiles) | Figure 5, Figures S7–S8, Tables S5, S8 |
| `alternative_tss/` | `model_alt_tss.py` (validation of held-out model TSS calls with 3 independent 5′-end datasets; gene-level alternative TSSs, family / GO enrichment) | Figure 5, Figure S9, Table S7 |
| `translation/` | `te_analysis.py` (uORF / uAUG census and polysome translation efficiency models), `riboseq_map.sbatch`, `riboseq_trim.py`, `ribo_orf.py` (Caro 2014 Ribo-seq), `daft_gene_expr.py` (per-gene DAFT-seq expression) | Figure 5, Figure S10, Table S7 |
| `novel_units/` | `novel_scan.py` (genome-wide out-of-fold transcription units) | optional supplement |
| `*/session_derived/` | derivations of tables computed interactively (pairs, co-expression, TSS decomposition, ApiAP2 track classes, ...) | |

`te_analysis.py` expects a per-gene covariate table `cov` to be defined before it is executed (see its first line).
Scripts in `chromatin_nucleosome/` are run from a working folder that contains `common.py` and the nucDetective
output under `ext/nuc/` (Holzinger et al.; Zenodo 16779899).
