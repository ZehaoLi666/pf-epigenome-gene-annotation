"""Session-derived computation of shaw2021_5CAPture_5CN_clusters.tsv.

Recovered from the record of the interactive analysis session in which shaw2021_5CAPture_5CN_clusters.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: GSE103036_combine_enrichment_score.txt.gz
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import gzip
import numpy as np
import pandas as pd

f = result('GSE103036_combine_enrichment_score.txt.gz')

rows = []
with gzip.open(f, 'rt') as fh:
    next(fh)
    for line in fh:
        x = line.rstrip('\n').split('\t')
        try: p = float(x[-1])
        except ValueError: continue
        if p >= 0.01: continue
        body = '\t'.join(x[3:-1]).split('\t*'); a = b = 0.0; ns = 0
        for g in body:
            t = [y for y in g.split('\t') if y]
            if t and t[0] != 'n/a': a += float(t[0]); b += float(t[1]); ns += 1
        rows.append((x[2], x[1], int(x[0]), p, a, b, ns))
C5 = pd.DataFrame(rows, columns=['chrom', 'strand', 'pos', 'p', 'd1', 'd2', 'n_samples'])

def cluster(df, gap=20, wcol='d1'):
    out = []
    for (c, s), q in df.sort_values('pos').groupby(['chrom', 'strand']):
        p = q.pos.values; w = q[wcol].values; br = np.r_[0, np.where(np.diff(p) > gap)[0] + 1, len(p)]
        for i, j in zip(br[:-1], br[1:]):
            k = i + int(np.argmax(w[i:j])); out.append((c, s, int(p[i]), int(p[j - 1]), int(p[k]), float(w[i:j].sum()), j - i))
    return pd.DataFrame(out, columns=['chrom', 'strand', 'start', 'end', 'dominant', 'score', 'n_pos'])

CAP = cluster(C5)
CAP['dominant'] -= 1   # 1-based -> 0-based
CAP = CAP[CAP.chrom.str.match(r'Pf3D7_\d\d_v3')]

CAP.to_csv('shaw2021_5CAPture_5CN_clusters.tsv', sep='\t', index=False)
