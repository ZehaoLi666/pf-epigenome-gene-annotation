# Chromatin-informed transcript boundaries in *Plasmodium falciparum*

Code by **Zehao Li** — a long-read–anchored, chromatin-informed model of transcription start
(TSS) and end (TTS) sites across the *P. falciparum* 3D7 genome, and the refined UTR annotation and biology derived
from it.



## Pipeline overview

```
 raw long reads (PacBio Iso-Seq, ONT cDNA, ONT direct RNA; 4 studies)        short-read RNA-seq (7 stages)
        │ longread_reference/01_gse186109_processing  (download, minimap2, BED12)     │ 03_short_read_support
        │ longread_reference/02_tama_collapse_hc_models (ONT orientation, TAMA,       │ (junctions, TPM)
        │                                    merge, high-confidence (HC) models)      │
        ▼                                                                             ▼
 HC transcript models ── 04_end_qc (5′ cap-G / poly(A)-tail read evidence) ── 05_transcript_catalogue
        │                     (per-gene catalogue of transcript ends vs PlasmoDB v68/v71)
        ▼
 long-read reference (2,389 genes / 9,142 stage-confirmed models)       ──► Figure 1, Table S2
        │ labels/  gene_set_v3.py → labels_v3.py (CAGE-free TSS/TTS/body/flank labels, 20-bp bins)
        │          flank_gaps/ (same-strand gap classification from long reads + DAFT-seq)
        │          labels_soft.py (Gaussian-softened targets)                   ──► Figure 2, Table S3
        ▼
 inputs (model/features/): 8 chromatin tracks (ChIP/ATAC/MNase), DNA one-hot, gLM (Caduceus)
        embeddings, GRO-seq, polysome RNA, inferred stranded RNA, DAFT-seq  (bigWig → 20-bp matrices)
        │ model/  exp.py / exp_ml.py (multi-task CNN, segmodel.py + boundary_model/bm)
        │         run_v3.sh (final models), run_v3cv.sh (chromosome-fold CV)
        │         predict_logit.py → peak_eval.py (logit-peak calling, thresholds on validation chromosomes)
        │         gene_boundaries.py (assignment of calls to genes → per-gene TSS/TTS)
        ▼
 per-gene refined boundaries (gene_boundaries_v3_all_genes.tsv)          ──► Figure 3, Table S4
        │ evaluation/benchmarks  (independent 5′-end and 3′-end truth sets; StringTie, ChromHMM, fixed offset, v71)
        │ evaluation/ablation     (three-fold chromosome CV, input-class and single-track ablation; HPCC)
        │                                                                         ──► Figure 4, Tables S5–S6
        │ annotation/            (refined GFF3 / TSV resource)                   ──► Data S1
        ▼
 analysis/  chromatin_nucleosome (promoter NDRs, nucleosome positions, chromatin states, head-to-head pairs)
            alternative_tss      (model-predicted alternative TSSs validated with 3 independent 5′-end sets)
            translation          (uORF / uAUG census, polysome TE, Ribo-seq)
            novel_units          (genome-wide out-of-fold transcription units; optional supplement)
                                                                                 ──► Figure 5, Tables S7–S9
```

Folder by folder (each folder has a short README):

| Folder | Contents |
|---|---|
| `longread_reference/` | Long-read processing on the HPCC cluster (SLURM), HC transcript models, read-end QC, transcript-end catalogue |
| `labels/` | Gene set and training labels (v3, CAGE-free), flank-gap tiers, soft labels |
| `model/` | Feature construction, network, training / CV launch scripts, prediction, peak calling, gene assignment |
| `evaluation/` | Benchmarks against independent truth sets and alternative methods; ablation (HPCC job scripts + scoring) |
| `annotation/` | Builds and validates the refined UTR annotation (GFF3 + TSV) |
| `analysis/` | Downstream biology: promoter/nucleosome, head-to-head pairs, alternative TSS, translation, novel units |
| `figures/Figure1..Figure5/` | Per-panel figure code: `prep_data.py` (inputs → source tables), `panel_*.py`, `run_all.py` |
| `supplementary/tables/` | Builders of Supplementary Tables S1–S9 (Excel + TSV) |
| `*/session_derived/` | Derivations of intermediate tables that were computed interactively (see below) |
| `config/` | Path configuration (`paths.example.yaml`, `hpcc.env.example`) |

## Installation

```bash
git clone https://github.com/ZehaoLi666/pf-epigenome-gene-annotation.git
cd pf-epigenome-gene-annotation
conda env create -f environment.yml          # analysis, figures, tables (Python 3.12)
conda env create -f environment-model.yml    # model training / prediction (PyTorch)
conda activate pf-boundaries
pip install -e .                             # makes `pfconfig` importable from every folder
cp config/paths.example.yaml config/paths.yaml   # then edit the paths
```

All scripts read external locations from `config/paths.yaml` through `pfconfig.py` (`LAB_DIR`, `RESULTS_DIR`,
`HPCC_BIGDATA`, `HPCC_LAB`, `HPCC_HOME`, `FIGURES_DIR`, `TABLES_OUT_DIR`; each can be overridden with an
environment variable `PF_<KEY>`). Cluster job scripts read the same locations from `config/hpcc.env`
(copy `config/hpcc.env.example`) and need `PF_REPO` to point at the repository checkout.
Cluster software paths in the job scripts are those of the UCR HPCC (`/opt/linux/rocky/8.x/...`); change them for
another site.

## Reproducing the figures

Each `figures/FigureN/` folder is self-contained:

```bash
cd figures/Figure1
python code/prep_data.py   # reads RESULTS_DIR + LAB_DIR, writes data/Figure1_*.tsv (per-panel source data)
python code/run_all.py     # draws every panel into panels/ (PDF, SVG, PNG) from data/ only
```

`run_all.py` needs only numpy, pandas, matplotlib (and scipy for some panels); `prep_data.py` also needs pyBigWig.
The per-panel source tables (`data/`) are deposited with the paper (see *Data availability*), so `run_all.py` can be
run without any upstream input. Panels are drawn at final size; panel letters and the final layout are added in
Adobe Illustrator.

Checks run for this release (October 2026): with the deposited source tables, `run_all.py` for Figures 1–5
regenerated all 41 PNG panels pixel-identical to the submitted panels; with the deposited intermediate tables,
`prep_data.py` for Figures 1–5 regenerated all 74 per-panel source tables byte-identical.

Supplementary tables: `python supplementary/tables/build_S1.py` (likewise `build_S2_S4.py`, `build_S5_S6.py`,
`build_S7_S9.py`) — they read RESULTS_DIR and the figure source tables in FIGURES_DIR and write to TABLES_OUT_DIR
(requires openpyxl).

## Reproducing the model

1. Long-read reference (HPCC): `longread_reference/01_*` → `02_*` → `03_*` → `04_*`, then
   `05_transcript_catalogue/` (see its README for the order).
2. Labels: `labels/gene_set_v3.py`, `labels/flank_gaps/*`, `labels/daft_stage20.py`, `labels/labels_v3.py`,
   then `labels/labels_soft.py` (SOFT_BASE=data/Lml3_pool.npz → `Lml3soft40_pool.npz`).
3. Features: `model/features/build_features.py` (bigWig → `data/raw_*.npz`, `data/meta.json`),
   `model/features/rebuild_inputs.py`, `model/features/glm_embed.sbatch` (GPU), `model/features/stranded_evidence.py`.
4. Training: `model/run_v3.sh` (final configuration `V3_top8`: 8 chromatin tracks + DNA + gLM + GRO + polysome +
   stranded RNA + DAFT-seq; three runs per setting: seed 0, seed 1, dropout) and `model/run_v3cv.sh`
   (folds B and C); fold A = chromosomes 1–5 test / 6–7 validation / 8–14 training.
5. Calling: `model/predict_logit.py` → `model/peak_eval.py` → `model/gene_boundaries.py`.
6. Evaluation: `evaluation/benchmarks/`, `evaluation/ablation/` (HPCC array jobs, `tasks.tsv` lists every setting).
7. Annotation: `annotation/build_resource.py` then `annotation/validate_resource.py`.

The run scripts expect the working-directory layout used during the project (`data/`, `exp_ml/`, `boundary_model/`
next to the scripts); run them from `model/` with `PYTHONPATH=.:boundary_model`.

## Session-derived code

Some intermediate tables (e.g. the benchmark truth sets, the method-prediction table, the per-gene boundary table,
the Figure 1 summary tables) were computed interactively rather than by a saved script. For these, the
computation was recovered from the recorded analysis session and saved as `session_derived/<table>.py` next to the
code it depends on. Each file states the table it produces and its inputs. These files document the exact computation
but may rely on objects defined earlier in the session; check them before re-running. Three of them were re-run for
this release (the Table S5 tolerance and per-truth-set tables) and reproduced the deposited tables byte-identically;
the others have not been re-run.

## Input data layout

`LAB_DIR` (local data root) is expected to contain:

```
PlasmoDB/PlasmoDB-48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta
PlasmoDB/PlasmoDB-68/PlasmoDB-68_Pfalciparum3D7.gff          (v68 gene structures = v71)
PlasmoDB/PlasmoDB-71/PlasmoDB-71_Pfalciparum3D7.gff
lncRNA project tracks/histone_marks/<mark>_<Stage>_minusInput.bw        (H2AZ, H2BZ, H3, H3K27ac, H3K4me1, H3K4me3, H3K9ac, H3K9me3)
lncRNA project tracks/accessibility_nucleosome/{ATAC_<Stage>_minusInput.bw, MNase_<Stage>_CPM.bw}
lncRNA project tracks/machinery/{HP1_<Stage>_minusInput.bw, PolII_<stage>_minusInput.bw}
lncRNA project tracks/apiap2_tf/<AP2>_<Stage>_minusInput.bw
lncRNA project tracks/nascent_transcription/GRO_<Stage>_mean_normalized.bw
lncRNA project tracks/translatome/Polysome_<stage>_CPM_final.bw
lncRNA project tracks/cage_tss/CTSS_<stage>_{pos,neg}_CPM.bw            (validation only)
IGV_longread_transcript_models/  (long-read model BEDs, DAFT-seq pooled bigWigs, genome FASTA + .fai)
```

`RESULTS_DIR` holds the intermediate result tables named in the scripts (`result("<name>")`), e.g.
`gene_boundaries_v3_all_genes.tsv`, `models_v3.tsv.gz`, `benchmark_truth_sets.tsv`,
`benchmark_method_predictions.tsv`, `ml3_end_positives.tsv`.

## Data availability

Large files are not stored in this repository (no file > 1 MB, no BAM/FASTQ/bigWig, no model weights).

| Item | Where |
|---|---|
| Raw sequencing data (public) | accessions listed in Table S1 (e.g. GEO GSE186109, GSE85478, GSE104075, GSE184658/GSE184659, GSE68982, GSE103036, GSE58402, GSE66185; SRA SRP021890, SRP026365, SRP252482; ArrayExpress E-MTAB-11766; ENA ERP001570) |
| Long-read reference (GTF/BED), refined annotation (GFF3/TSV), model prediction tracks (bigWig) | Data S1–S3 — to be deposited (Zenodo, DOI to be added) |
| Intermediate result tables (`RESULTS_DIR`) and per-panel figure source tables (`figures/*/data`) | to be deposited (Zenodo, DOI to be added) |
| Trained model weights (final and CV models) | to be deposited (Zenodo, DOI to be added) |
| Processed bigWig input tracks | to be deposited (Zenodo / GEO, accession to be added) |

## Software

Analysis and figures: Python 3.12 with numpy 2.5, pandas 3.0, scipy 1.18, matplotlib 3.11, statsmodels 0.15,
pyBigWig 0.3.26, openpyxl 3.1 (`environment.yml`). Model: PyTorch 2.13, scikit-learn 1.9, lightgbm 4.7
(`environment-model.yml`); the HPCC training runs used a separate conda environment whose exact versions are not
recorded here. Command-line tools on the cluster: minimap2, samtools, bedtools, TAMA (Python 2.7), seqkit, SRA toolkit,
hisat2/STAR, featureCounts, bowtie2, deepTools, StringTie, ChromHMM, SQANTI3, SLURM.

## Not yet included

See [`MISSING.md`](MISSING.md) for pipeline steps whose scripts exist only on the HPCC cluster and still need to be
added.

## AI assistance

Code was developed with AI assistance; all code was reviewed by the authors.

## License

MIT (see `LICENSE`).
