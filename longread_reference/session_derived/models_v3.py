"""Session-derived computation of models_v3.tsv.gz.

Recovered from the record of the interactive analysis session in which models_v3.tsv.gz was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: transcript_stage_confirmation.tsv.gz
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import glob, numpy as np, pandas as pd, sys

LAB = LAB_DIR
S = ['ring', 'troph', 'schizont']

tsc = pd.read_csv(result('transcript_stage_confirmation.tsv.gz'), sep='\t', low_memory=False)

for s in S:
    tsc[f'conf3_{s}'] = (tsc[f'detect_{s}'] == 'detected') & tsc.struct_ok & tsc[f'sj_{s}'] & tsc.tes_ok
tsc['conf3_any'] = tsc[[f'conf3_{s}' for s in S]].any(axis=1)

g = pd.read_csv(glob.glob(f'{LAB}/PlasmoDB/PlasmoDB-68/*.gff')[0], sep='\t', comment='#', header=None, low_memory=False, names=['chrom','src','type','start','end','score','strand','phase','attr'])
cds = g[g.type == 'CDS'].copy(); cds['gene'] = cds.attr.str.extract(r'Parent=([^;.]+)')[0]
env = cds.groupby('gene').agg(cds_lo=('start', 'min'), cds_hi=('end', 'max')); env['cds_lo'] -= 1

T = tsc.join(env, on='locus'); mr = T.unit_class == 'mRNA'
T['full_cds'] = mr & T.cds_lo.notna() & (T.start <= T.cds_lo) & (T.end >= T.cds_hi)
T['trunc5'] = mr & T.cds_lo.notna() & ~T.full_cds & np.where(T.strand == '+', T.start > T.cds_lo, T.end < T.cds_hi)
T['trunc3'] = mr & T.cds_lo.notna() & ~T.full_cds & ~T.trunc5

T.to_csv('models_v3.tsv.gz', sep='\t', index=False)
