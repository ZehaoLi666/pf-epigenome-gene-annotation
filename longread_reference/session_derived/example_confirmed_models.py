"""Session-derived computation of example_confirmed_models.tsv.

Recovered from the record of the interactive analysis session in which example_confirmed_models.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_boundaries_v3_all_genes.tsv, hc_reference.bed.gz, gene_transcript_end_catalogue.tsv, gene_vs_readend_v48.tsv, transcript_stage_confirmation.tsv.gz
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import gzip
import numpy as np
import pandas as pd

T = pd.read_csv(result('transcript_stage_confirmation.tsv.gz'), sep="\t", low_memory=False)
S3 = ['ring', 'troph', 'schizont']
T['conf'] = np.any([(T[f'detect_{s}'] == 'detected') & T.struct_ok & T[f'sj_{s}'] & T.tes_ok for s in S3], axis=0)

p = result('hc_reference.bed.gz')
rows = []
with gzip.open(p, 'rt') as fh:
    for ln in fh:
        if ln.startswith(('track', '#')): continue
        f = ln.rstrip('\n').split('\t'); s0 = int(f[1]); sz = [int(x) for x in f[10].strip(',').split(',')]; st = [int(x) for x in f[11].strip(',').split(',')]
        rows.append((f[0], s0, int(f[2]), f[3], f[5], [(s0+a, s0+a+b) for a, b in zip(st, sz)]))
HB = pd.DataFrame(rows, columns=['chrom', 'start', 'end', 'name', 'strand', 'blocks'])
HB['tid'] = HB.name.str.split(';').str[1]
HB['locus'] = HB.name.str.split(';').str[0]

G = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep="\t").set_index('gene_id')

conf = set(T.loc[T.conf, 'transcript_id'])
EXg = ["PF3D7_1302800", "PF3D7_0716900", "PF3D7_1323700"]
M = HB[HB.tid.isin(conf) & HB.locus.isin(EXg)].copy()
lo = np.minimum(G.atg, G.stop)
hi = np.maximum(G.atg, G.stop)
M['cds_lo'] = M.locus.map(lo) - 1
M['cds_hi'] = M.locus.map(hi)
M['full_cds'] = (M.start <= M.cds_lo) & (M.end >= M.cds_hi)
M['blocks'] = M.blocks.apply(lambda b: ';'.join(f"{a}-{c}" for a, c in b))
M.drop(columns=['name']).to_csv("example_confirmed_models.tsv", sep="\t", index=False)
