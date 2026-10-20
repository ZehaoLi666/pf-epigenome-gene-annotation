# Model

Multi-task 1-D convolutional network that predicts, per strand and 20-bp bin, TSS, TTS and gene-body probabilities
from chromatin, DNA and RNA inputs.

| File | Role |
|---|---|
| `features/build_features.py` | bigWig tracks → per-stage 20-bp matrices `data/raw_<kind>_<stage>.npz`, `data/meta.json` |
| `features/rebuild_inputs.py` | normalised input arrays `data/X_chrom_<stage>.npz`, `data/X_rna_<stage>.npz` |
| `features/glm_embed.py`, `features/glm_embed.sbatch` | Caduceus DNA language-model embeddings (PCA-16, 20 bp; GPU, HPCC) |
| `features/stranded_evidence.py` | strand-specific evidence from unstranded RNA-seq BAMs (HPCC) |
| `pfprep.py`, `host_paths.py` | normalisation, chromosome splits, annotation helpers, result-file locations |
| `seg_features.py` | assembles input channels for an input specification (`chrom,onehot,glm,gro,poly,strand,daft`) |
| `segmodel.py`, `boundary_model/bm/` | network blocks (dilated residual CNN / U-Net), decoding and evaluation helpers |
| `exp.py`, `exp_ml.py` | training (multi-label sigmoid heads, strand-flip augmentation) |
| `run_v3.sh` | final settings (`V3_top8` = 8 chromatin tracks; plus `V3_all`, `V3_noAP2`, `V3_nochrom`, `V3_old_inputs`), 3 runs each |
| `run_v3cv.sh` | chromosome-fold cross-validation, folds B (test 6–10, val 11) and C (test 11–14, val 10) |
| `predict_logit.py` | averaged forward + strand-flipped logits (`exp_ml/logit*_<name>.npz`) |
| `peak_eval.py` | logit-peak calling, thresholds / NMS tuned on validation chromosomes, one-to-one matching |
| `gene_boundaries.py` | assigns out-of-fold calls to genes → per-gene TSS / TTS / UTRs |
| `gene_truth.py`, `strata_eval.py`, `cage_validate.py`, `bio_profiles.py` | long-read truth per gene, recall by expression stratum, independent CAGE check, chromatin read-outs at anchors |
| `session_derived/` | derivation of the final per-gene table `gene_boundaries_v3_all_genes.tsv` |

Run from this folder with `export PYTHONPATH=.:boundary_model POOL=1 ML_LABELS=data/Lml3soft40_pool.npz`.
`pfprep.py` is the latest saved version (the version in the archived v3 bundle was not stored separately);
`boundary_model/bm` and `segmodel.py` are the versions saved with the segmentation code. Weights are not in the
repository (see top-level README, *Data availability*).
