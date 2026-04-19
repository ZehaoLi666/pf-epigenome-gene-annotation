# Epigenome-based Gene Boundary Prediction in Plasmodium falciparum

This repository contains the published epigenomic boundary-prediction workflow used to identify candidate transcript starts and ends in *Plasmodium falciparum*. The supported public workflow is the `window1000_bin100` pipeline, which builds binned chromatin features around annotated gene boundaries, trains classifiers, scores intergenic sites, and exports candidate loci for downstream validation.


## What Is Supported

The actively supported workflow is:

1. build positive training features from known gene boundaries,
2. build or regenerate negative examples,
3. train the TSS/TTS classifiers,
4. build intergenic feature matrices for genome-wide scoring,
5. score and refine candidate boundaries,
6. optionally compare predicted versus annotated UTR lengths.

Everything under `archive/` is retained for provenance only and is not part of the supported public method.

## Repository Layout

- `window1000_bin100/`: main published workflow.
- `validation_plot/`: lightweight downstream validation plots.
- `archive/`: legacy figure scripts, one-off helpers, and placeholders.
- `window500_bin50/`, `Figure*`, `model_comparison/`, `other_model/`: exploratory or alternative analyses kept for reference.

## Software Requirements

### Python

Create an environment and install:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements-python.txt
```

Core Python packages used by the supported workflow include:

- `numpy`
- `pandas`
- `scikit-learn`
- `xgboost`
- `joblib`
- `matplotlib`
- `statsmodels`
- `tqdm`


### Command-Line Tools

Some optional regeneration and filtering steps expect:

- `bedtools`
- `samtools`
- Slurm (`sbatch`) for the provided cluster wrapper scripts

## External Inputs

The repository does not track large per-base chromatin coverage files. The active workflow expects marker files in plain text format with three whitespace-delimited columns:

```text
chromosome  position  value
```

Positions are assumed to be 1-based. One file is needed per marker, for example:

- `ATAC-seq.txt`
- `H3K4me3.txt`
- `H3K27ac.txt`
- `MNseq_sub.txt`

### Marker Configuration

The tracked marker JSON files are templates:

- `window1000_bin100/training_features_values/marker_files.json`
- `window1000_bin100/training_features_values/marker_files_ATAC.json`

They use `${LNC_CHIP_MARKER_DIR}` placeholders. You can configure inputs in either of two ways:

1. Set one environment variable:

```bash
export LNC_CHIP_MARKER_DIR=/absolute/path/to/per_base_marker_tracks
```

2. Create ignored local overrides by copying the tracked JSON templates to:

- `window1000_bin100/training_features_values/marker_files.local.json`
- `window1000_bin100/training_features_values/marker_files_ATAC.local.json`

The shell wrappers automatically prefer `*.local.json` when present.

### Lightweight Reference Inputs Already Tracked

These files are already in the repository and are sufficient to run the documented workflow once marker tracks are configured:

- `window1000_bin100/training_features_values/reference_gene_boundaries.csv`
- `window1000_bin100/training_features_values/negatives_2k.tsv`
- `window1000_bin100/prediction/intergenic_midpoints_100bp.csv`
- `window1000_bin100/prediction/results_intergenic_TSS/PlasmoDB-48_Pfalciparum3D7_Genes.bed`

Optional regeneration helpers additionally use:

- `window1000_bin100/training_features_values/Pf3D7.chrom.sizes`
- `window1000_bin100/training_features_values/telomere_fixed.bed`
- `window1000_bin100/training_features_values/Merged.bed`
- a genome `.fai` file for chromosome sizes during final TSS filtering

## Main Workflow

### 1. Build Positive Training Features

Use the tracked reference gene boundaries and the configured marker panel:

```bash
cd window1000_bin100/training_features_values
sbatch assign_value.sh
python3 merge_feature_chunks.py --input-dir positive_values
```

Expected outputs:

- `positive_values/features_chunk_*_TSS_features_±1000bp_20bins.csv`
- `positive_values/features_chunk_*_TTS_features_±1000bp_20bins.csv`
- `positive_values/merged_TSS_500bp_20bins.csv`
- `positive_values/merged_TTS_500bp_20bins.csv`

The merged filenames are kept for compatibility with the original training scripts even though the active window is `±1000 bp` with `20` bins.

### 2. Build Negative Training Features

If you want to reuse the tracked negatives:

```bash
cd window1000_bin100/training_features_values
sbatch assign_value_Neg.sh
python3 negative_values/split.py \
  --in-csv ATAC_negatives_2k_features.csv \
  --out-prefix negative_values/ATAC_negatives_2k_features
```

Expected outputs:

- `ATAC_negatives_2k_features.csv`
- `negative_values/ATAC_negatives_2k_features_TSS_only.csv`
- `negative_values/ATAC_negatives_2k_features_TTS_only.csv`

If you want to regenerate negatives from BED references instead of reusing `negatives_2k.tsv`, use:

- `window1000_bin100/training_features_values/make_quiet_intergenic_from_beds.py`
- `window1000_bin100/training_features_values/generate_banded_negatives_v2.py`

### 3. Train the Classifier

The default training wrapper trains from the TSS feature tables:

```bash
cd window1000_bin100/model_training
sbatch model.training.sh
```

Key outputs:

- trained model directory under `model_training/`
- per-run metrics CSVs and summary JSON
- ROC/PR/confusion-matrix figures

The script defaults to:

- positives: `training_features_values/positive_values/merged_TSS_500bp_20bins.csv`
- negatives: `training_features_values/negative_values/ATAC_negatives_2k_features_TSS_only.csv`

You can override `POS`, `NEG`, and `OUT` via environment variables.

### 4. Build Intergenic Feature Matrices for Prediction

For the full marker panel:

```bash
cd window1000_bin100/prediction
sbatch run_bins_per_marker.sh
python3 merge_intergenic_bins.py --input-dir intergenic_bins
```

For ATAC only:

```bash
cd window1000_bin100/prediction
sbatch assign_value.sh
```

Expected outputs:

- `intergenic_bins/*_intergenic_bins_pm1000bp_20bins.csv`
- `merged_intergenic_features.csv`

### 5. Score and Refine Candidate Boundaries

The repository retains the original scoring notebooks for the final interactive candidate-calling step:

- `window1000_bin100/prediction/prediction_TSS.ipynb`
- `window1000_bin100/prediction/prediction_TTS.ipynb`
- `window1000_bin100/prediction/candidates_cleaning.ipynb`

These consume:

- `merged_intergenic_features.csv`
- the trained model output from `model_training/`
- the TSS/TTS result directories in `prediction/`

The notebooks are preserved as interactive analysis artifacts. Their first configuration cells may still need local path updates for model files, GFF annotation, or genome index files that are not distributed in this repository.

The final strand-aware TSS refinement step is also available as a script:

```bash
cd window1000_bin100/prediction/results_intergenic_TSS
CHROM_SIZES=/absolute/path/to/genome.fai bash filter_TSS.sh
```

Expected outputs include:

- `prediction/results_intergenic_TSS/predictions_all.csv`
- `prediction/results_intergenic_TSS/refined_TSS_intergenic.tsv`
- `prediction/results_intergenic_TTS/predictions_all_TTS_thr0.7.csv`
- `prediction/novel_units.tsv`
- `prediction/novel_units.bed`



## Notes For Reuse

- The active scripts no longer require the original `/rhome` or `/bigdata` lab paths.
- Marker file paths may be absolute, relative, `~`-expanded, or `${ENV_VAR}`-expanded.
- Relative marker paths inside the JSON configs are resolved relative to the JSON file location.
- Large binaries, joblib models, compressed matrices, and logs remain ignored by `.gitignore`.

## Citation

If you publish work derived from this repository, cite the associated paper and note any modifications you made to the workflow or marker set.
