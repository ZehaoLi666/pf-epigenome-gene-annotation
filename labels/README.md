# Training labels (v3, CAGE-free)

1. `gene_set_v3.py <transcript_stage_confirmation.tsv.gz>` — gene set: loci with ≥ 1 stage-confirmed full-CDS model
   (confirmation per stage = detected & structure OK & splice junctions supported & 3′ end QC); writes
   `models_v3.tsv.gz`, `gaps_v3.tsv`, `gene_side_gaps_v3.tsv`.
2. `flank_gaps/` — same-strand gaps next to confirmed genes: `daft_map.sbatch` + `daft_cov.py` (stranded DAFT-seq
   coverage), `strand_check.py` (library orientation), `gap_stats.py` (per-gap evidence), `gap_classify.py`
   (tiers A/B/C or spanned).
3. `daft_stage20.py` — stage-binned stranded DAFT-seq (20-bp bins) used as model input and in label QC.
4. `labels_v3.py` — 20-bp TSS / TTS / body / flank labels with weights (cap-G and poly(A)-tail support from
   `read_end_qc.tsv`; internal-priming flags from `end_flags_min.tsv`); writes `Lml3_pool.npz`,
   `ml3_end_positives.tsv` and `flank_gaps_v3_longread_tiers.tsv` (Table S3, Figure 2).
5. `labels_soft.py` — Gaussian-softened end targets: `SOFT_BASE=data/Lml3_pool.npz SOFT_POS=data/ml3_end_positives.tsv python labels_soft.py 40 data/Lml3soft40_pool.npz`
   → `Lml3soft40_pool.npz` (training targets used by `model/run_v3.sh`).

Scripts use the project working layout (`v3/`, `work/data/`); run them from a folder with that layout.
