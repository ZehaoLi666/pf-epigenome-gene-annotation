"""Classify same-strand gaps next to long-read confirmed genes (UTR_BOUNDARY_MODEL_DESIGN.md v2, section 2).
Inputs: gap_stats.tsv (per gap x long-read dataset), gap_stats_daft.tsv (stranded DAFT-seq), gap_local_evidence.tsv (CAGE, A content, AT).
Criteria (same strand):
  1 no spanning long read       : span_total == 0, or span_total <= 2 and span_total < 2% of the weaker side's reads ('near zero')
  2 end evidence                : upstream-unit TTS = direct-RNA 3' end within 50 bp, or A-tailed read end, or >=2 non-A-rich cDNA 3' ends
                                  with genomic downstream A < 60%;  downstream-unit TSS = CAGE >= 1 CPM (+-100 bp, any stage) or cap G read
  3 stranded coverage drop      : min 100-bp DAFT-seq window in the gap < 25% of the weaker 300-bp unit edge (pooled CPM; edge >= 0.5 CPM)
  4 enough reads on both sides  : >= 5 reads touching each unit edge (summed over datasets)
Tiers: Not flank (reads span) | Uncertain (too few reads) | Flank strong (1+2 both+3+4) | Flank weak (1+4, one end or no drop) | Uncertain (no end evidence)."""
import numpy as np, pandas as pd
DRNA = {'LeeVV2021_dRNA', 'Hoshizaki2022_dRNA'}
def classify(stats, daft, local):
    S = stats.copy(); S['is_drna'] = S.dataset.isin(DRNA)
    agg = S.groupby('gap_id').agg(left=('left', 'sum'), right=('right', 'sum'), span=('span', 'sum'), n_ds_span=('span', lambda x: int((x > 0).sum())),
                                  n3=('n3', 'sum'), n3_tail=('n3_tail', 'sum'), n3_arich=('n3_arich', 'sum'), n5=('n5', 'sum'), n5_G1=('n5_G1', 'sum'))
    agg['n3_drna'] = S[S.is_drna].groupby('gap_id').n3.sum().reindex(agg.index).fillna(0)
    pool = daft[daft.track == 'pool'].set_index('gap_id')[['gap_mean', 'gap_min100', 'left_edge', 'right_edge', 'antisense_gap_mean']]
    A = local.set_index('gap_id').join(agg).join(pool)
    A['weak_side'] = A[['left', 'right']].min(axis=1)
    A['c1_no_span'] = (A.span == 0) | ((A.span <= 2) & (A.span < 0.02 * A.weak_side))
    A['c1_strict'] = A.span == 0
    A['tts_ok'] = (A.n3_drna >= 1) | (A.n3_tail >= 1) | (((A.n3 - A.n3_arich) >= 2) & (A.downA_tts < 0.6))
    A['tss_ok'] = (A.cage_tss >= 1) | (A.n5_G1 >= 1)
    A['n_end'] = A.tts_ok.astype(int) + A.tss_ok.astype(int)
    edge = A[['left_edge', 'right_edge']].min(axis=1); A['edge_min'] = edge
    A['cov_ratio'] = A.gap_min100 / edge.where(edge > 0)
    A['cov_informative'] = edge >= 0.5
    A['c3_drop'] = A.cov_informative & (A.cov_ratio < 0.25)
    A['c4_enough'] = (A.left >= 5) & (A.right >= 5)
    A['AT_flag'] = A.AT95_frac > 0.5
    def tier(r):
        if not r.c1_no_span: return 'Not flank: reads span the gap'
        if not r.c4_enough: return 'Uncertain: too few reads'
        if r.n_end == 2 and r.c3_drop: return 'Flank (strong)'
        if r.n_end >= 1: return 'Flank (weak)'
        return 'Uncertain: no end evidence'
    A['tier'] = A.apply(tier, axis=1)
    return A
TIERS = ['Flank (strong)', 'Flank (weak)', 'Uncertain: no end evidence', 'Uncertain: too few reads', 'Not flank: reads span the gap']
