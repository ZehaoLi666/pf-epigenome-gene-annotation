"""Session-derived computation of gene_heterochromatin_groups.tsv.

Recovered from the record of the interactive analysis session in which gene_heterochromatin_groups.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: meta.json
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import sys, numpy as np, pandas as pd
sys.path[:0] = ['.', 'boundary_model']
import pfprep as P, seg_discover as SD
import json

RES = 20; S = P.STAGES

T, g48, g68, lnc = SD.intervals()
g = g68[g68.biotype.isin(['protein_coding_gene', 'pseudogene'])].copy()
d = g.desc.fillna('').str.lower()
g['family'] = np.select([d.str.contains('erythrocyte membrane protein 1'), d.str.contains('rifin'), d.str.contains('stevor'), d.str.contains('pfmc-2tm')],
                        ['var', 'rifin', 'stevor', 'Pfmc-2TM'], 'other')

META = json.load(open(result('meta.json')))
TR = META['chrom_tracks']; iH, iK = TR.index('HP1'), TR.index('H3K9me3')
hp1 = np.zeros(len(g)); k9 = np.zeros(len(g))
for s in S:
    z = np.load(f'data/raw_chrom_{s}.npz')
    for i, r in enumerate(g.itertuples(index=False)):
        a = z[r.chrom]; lo, hi = r.start // RES, max(r.start // RES + 1, r.end // RES)
        hp1[i] = max(hp1[i], a[iH, lo:hi].mean()); k9[i] = max(k9[i], a[iK, lo:hi].mean())
g['HP1_body'] = hp1; g['H3K9me3_body'] = k9
thrH, thrK = np.quantile(hp1, 0.9), np.quantile(k9, 0.9)
g['HP1_domain'] = (hp1 >= thrH) & (k9 >= thrK)
g['group'] = np.where(g.family != 'other', 'var/rifin/stevor/Pfmc-2TM', np.where(g.HP1_domain, 'other HP1+H3K9me3-high', 'euchromatic'))
g.to_csv('gene_heterochromatin_groups.tsv', sep='\t', index=False)
print(g.group.value_counts().to_dict(), 'thr HP1', round(thrH, 3), 'K9', round(thrK, 3))
