# Long-read reference

Builds the long-read transcript reference (2,389 genes / 9,142 stage-confirmed models in the paper) from four
long-read studies (Shaw 2022 GSE186109 PacBio Iso-Seq + ONT cDNA; Lee 2021 and Hoshizaki 2022 ONT direct RNA;
Yang 2021 PacBio Iso-Seq). Most steps ran on the UCR HPCC cluster (SLURM); paths come from `config/hpcc.env`.

| Step | Folder / script | What it does |
|---|---|---|
| 1 | `01_gse186109_processing/` `download2.sh`, `download3.sh`, `organize.sh`, `qc.sh`, `launch_all.sh` | SRA download, per-sample organisation, read QC (seqkit) |
| 1 | `align_array.sbatch`, `make_bed12.sh` | minimap2 splice alignment to 3D7 (PlasmoDB v48 genome), BAM → BED12 |
| 1 | `cluster_ends.py`, `decay_filter.py`, `offset_stats.py`, `make_summary.py` | read-end clustering, filter of 5′-truncated (decay) reads, end-offset statistics vs v48 genes |
| 2 | `02_tama_collapse_hc_models/` `ont_orient.py`, `orient_check.py`, `ont_orient_align.sbatch`, `ont_ub_array.sbatch`, `ont_c14_cap.sbatch`, `ont_tama_chr.sbatch` | ONT read orientation and alignment, per-chromosome TAMA collapse of ONT reads |
| 2 | `make_tama_env.sh`, `tama_collapse_up.sbatch`, `tama_collapse_relaxed.sbatch`, `fl_counts.sh` | TAMA collapse of the PacBio / capped sets, full-length read counts per isoform |
| 2 | `merge_evidence.sbatch` (+ `merge_filelist.txt`), `evidence_tiers.py` | TAMA merge into the master reference; evidence layers, chimera guard, high-confidence (HC) models. Its CAGE-based tiers / stage tests are **superseded** by the CAGE-free confirmation in `labels/gene_set_v3.py`; the HC model set is still the starting point. |
| 3 | `03_short_read_support/` `sr_run.sh`, `sr_junctions.py`, `sr_pool_tpm.py`, `end_contexts.py` | short-read (7-stage RNA-seq) junction support and TPM; transcript classes and end contexts |
| 4 | `04_end_qc/` `read_end_qc.py`, `lr_reprocess.sbatch`, `lrcov2.sbatch` | read-level 5′ cap-G and 3′ poly(A)-tail evidence per end (`read_end_qc.tsv`, used by `labels/labels_v3.py`); re-processing / coverage of the four datasets |
| 5 | `05_transcript_catalogue/` `tama_models.py`, `utr_stats.py`, `utr_analysis.py`, `utr_classify.py`, `transcript_catalog.py`, `run_catalog.py`, `longread_gene_assign_hpcc.py` | export of collapsed models, read connectivity across annotated UTRs, per-gene catalogue of transcript-end situations vs PlasmoDB v68/v71 |
| – | `session_derived/` | derivations of the reference summary tables (models_v3, v71 comparison, funnel, coverage, internal priming, cap-G groups) |

`run_catalog.py` is a fragment executed inside `transcript_catalog.py`'s namespace (it expects the objects listed in
its first line).
