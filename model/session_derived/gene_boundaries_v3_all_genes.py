"""Session-derived computation of gene_boundaries_v3_all_genes.tsv.

Recovered from the record of the interactive analysis session in which gene_boundaries_v3_all_genes.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs: oof_top8_cage.tsv, oof_nochrom_cage.tsv (outputs of model/gene_boundaries.py + model/cage_validate.py)
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
import numpy as np
import pandas as pd

# Out-of-fold per-gene calls of the final model (V3_top8) and of the chromatin-free model, each annotated with CAGE
# support by model/cage_validate.py (oof_<name>_genes.tsv.gz -> oof_<name>_cage.tsv); files in the working directory.
A_ = pd.read_csv('oof_top8_cage.tsv', sep='\t'); B_ = pd.read_csv('oof_nochrom_cage.tsv', sep='\t')

G = A_.merge(B_[['gene_id', 'pred_TSS', 'pred_TTS']].rename(columns={'pred_TSS': 'nc_TSS', 'pred_TTS': 'nc_TTS'}), on='gene_id')

conf = G.lr_group == 'confirmed'
G['final_TSS'] = np.where(conf, G.lr_TSS, G.pred_TSS)
G['final_TTS'] = np.where(conf, G.lr_TTS, G.pred_TTS)
G['final_src'] = np.where(conf, 'long-read', np.where(G.pred_TSS.notna(), 'model', 'none'))
G['anchor_TSS'] = np.where(conf, G.lr_TSS, G.nc_TSS)
G['anchor_TTS'] = np.where(conf, G.lr_TTS, G.nc_TTS)
G['utr5_final'] = (G.atg - G.final_TSS).abs()
G['utr3_final'] = (G.final_TTS - G.stop).abs()
G['utr5_anchor'] = (G.atg - G.anchor_TSS).abs()

cols = ['gene_id', 'chrom', 'strand', 'desc', 'atg', 'stop', 'lr_group', 'final_src', 'final_TSS', 'final_TTS', 'utr5_final', 'utr3_final',
        'pred_TSS', 'pred_TSS_prob', 'alt_TSS', 'n_TSS_calls', 'pred_TTS', 'pred_TTS_prob', 'alt_TTS', 'n_TTS_calls', 'nc_TSS', 'nc_TTS',
        'lr_TSS', 'lr_TTS', 'annot_tss', 'annot_tts', 'utr5_annot', 'utr3_annot', 'cage_pred_TSS', 'pct_pred_TSS', 'pct_annot_tss', 'pct_dist_matched', 'd_pred_to_cage_peak']

OUT = G[cols].rename(columns={
    'pred_TSS': 'model_TSS',
    'pred_TSS_prob': 'model_TSS_prob',
    'alt_TSS': 'model_alt_TSS',
    'pred_TTS': 'model_TTS',
    'pred_TTS_prob': 'model_TTS_prob',
    'alt_TTS': 'model_alt_TTS',
    'nc_TSS': 'chromfree_model_TSS',
    'nc_TTS': 'chromfree_model_TTS',
    'cage_pred_TSS': 'CAGE_cpm_at_model_TSS',
    'pct_pred_TSS': 'CAGE_pct_model_TSS',
    'pct_annot_tss': 'CAGE_pct_v68_TSS',
    'pct_dist_matched': 'CAGE_pct_distance_matched_ctrl',
    'd_pred_to_cage_peak': 'dist_model_TSS_to_CAGE_peak'
})

OUT.to_csv('gene_boundaries_v3_all_genes.tsv', sep='\t', index=False)
