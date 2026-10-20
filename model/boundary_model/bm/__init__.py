"""bm — chromatin-based transcript-boundary model scaffold.

Modules
-------
schema    reference-contract schema + validator (what the reference must deliver)
tracks    bigWig -> binned matrix, per-stage normalisation, modality-partition guard
labels    soft (Gaussian) stranded label tracks + ignore masks from a boundary table
splits    chromosome hold-out and window tiling
models    dilated residual 1D CNN, multi-task, stage-conditioned (FiLM), strand-flip aware
losses    masked soft-label BCE / focal
baselines window-feature RF / LightGBM / logistic (paper-comparable rows)
decode    peak calling + gene-grammar decoding (weighted interval scheduling)
evaluate  AUPRC, positional error, recall@tolerance, track-shift null FDR, stratification
"""
__version__ = "0.1.0"
