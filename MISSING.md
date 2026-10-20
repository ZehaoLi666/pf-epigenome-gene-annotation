# Steps not yet in this repository

The scripts below were run on the UCR HPCC cluster (or produced files outside the project's saved code) and could
not be recovered for this release. Their outputs are used by the code here and will be deposited with the data;
the scripts themselves are **to add**. Nothing in this list was reconstructed or guessed.

| # | Step | Output used here | Used by |
|---|---|---|---|
| 1 | Alignment and TAMA collapse of Lee 2021 and Hoshizaki 2022 (ONT direct RNA) and Yang 2021 (PacBio Iso-Seq) long reads | per-dataset BED / `*_trans_read.bed` in `HPCC_LAB/SRA_toolkit/Long_read_RNA_datasets/` | `longread_reference/05_transcript_catalogue/tama_models.py` |
| 2 | SQANTI3 classification of the long-read models | SQANTI3 classification table | `longread_reference/03_short_read_support/end_contexts.py` (optional input), Figure S1 |
| 3 | Short-read RNA-seq (López-Barragán 2011, 7 stages) alignment (hisat2) and featureCounts | BAMs / counts in `HPCC_LAB/SRA_toolkit/RNAseq_7stages/processed_data/` | `03_short_read_support/sr_*.py`, `model/features/stranded_evidence.py` |
| 4 | Gene-vs-read-end table of the GSE186109 reads | `gene_vs_readend_v48.tsv` | `01_gse186109_processing/offset_stats.py`, several `session_derived/` scripts |
| 5 | Internal-priming / no-tail flags of long-read 3′ ends | `end_flags_min.tsv` | `labels/labels_v3.py` |
| 6 | Job wrapper (ends file + BAM list) for read-level end QC | `read_end_qc.tsv` | `labels/labels_v3.py` (the script itself is `longread_reference/04_end_qc/read_end_qc.py`) |
| 7 | Invocation of `gap_stats.py` / `gap_classify.py` on the v3 gap set | `gap_stats_v3.tsv` | `labels/labels_v3.py` |
| 8 | Processing of the ChIP-seq, ATAC-seq, MNase-seq, GRO-seq and polysome data into bigWig tracks (`*_minusInput.bw`, `*_CPM.bw`, `*_mean_normalized.bw`) | `LAB_DIR/lncRNA project tracks/` | `model/features/build_features.py`, figures, analyses |
| 9 | Processing of Adjalley 2016 cap-seq into CTSS bigWigs (`CTSS_<stage>_{pos,neg}_CPM.bw`) | `LAB_DIR/lncRNA project tracks/cage_tss/` | validation only (`model/cage_validate.py`, Figures 1, 3, 4) |
| 10 | Pooled DAFT-seq bigWigs for genome browsing (`DAFTseq_pool_{plus,minus}_CPM_10bp.bw`) | `LAB_DIR/IGV_longread_transcript_models/` | Figures 1 and 3 (track panels) |
| 11 | Read evidence (non-templated G) per TSS group | `tss_read_evidence_by_group.tsv` | Figure 1D |
| 12 | Conda environment used for training / ablation on the cluster (`enformer_pure`) | — | `evaluation/ablation/env.sh` |
| 13 | Download of the Caduceus DNA language-model checkpoint used by `glm_embed.py` | local checkpoint | `model/features/glm_embed.py` |
| 14 | Exact command lines of `gene_boundaries.py` for the final fold-A/B/C ensembles (`oof_top8`) | `oof_top8_genes.tsv.gz`, `oof_top8_oof_prob.npz` | `model/session_derived/gene_boundaries_v3_all_genes.py` |
| 15 | Derivation of the per-gene long-read coverage table and of the distal-TSS-share vs TE table (computed interactively; the recovered session code contained placeholders and was not included) | `gene_longread_coverage.tsv.gz`, `altTSS_distal_share_vs_TE.tsv` | Figure 1B, Table S7 |

External processing not done in this project: nucleosome calls and qNorm MNase tracks (nucDetective; Holzinger et
al., Zenodo 16779899); HMM chromatin states (Brown, Llinás & Mahony; GitHub seqcode/pfal-chrom-states-trackhub).
