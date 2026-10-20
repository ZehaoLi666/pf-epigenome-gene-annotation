"""Session-derived computation of siegel2014_pas_per_gene.tsv.

Recovered from the record of the interactive analysis session in which siegel2014_pas_per_gene.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: 12864_2013_7013_MOESM10_ESM.xlsx
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import numpy as np
import pandas as pd
import openpyxl

wb = openpyxl.load_workbook(result('12864_2013_7013_MOESM10_ESM.xlsx'), read_only=True)
ws = wb[wb.sheetnames[0]]
rows = list(ws.iter_rows(values_only=True))
h = [i for i, r in enumerate(rows) if r and r[0] == 'gene ID'][0]
SG = pd.DataFrame([r[:11] for r in rows[h + 1:] if r[0]], columns=['gene_id', 'strand', 'orf3', 'orf5', 'n_clusters', 'd_closest', 'pos_closest', 'd_furthest', 'pos_furthest', 'd_major', 'pos_major'])
for c_ in ['orf3', 'orf5', 'n_clusters', 'pos_closest', 'pos_furthest', 'pos_major', 'd_major']:
    SG[c_] = pd.to_numeric(SG[c_], errors='coerce')
SG['PAS_major'] = SG.pos_major - 1
SG['PAS_closest'] = SG.pos_closest - 1
SG['PAS_furthest'] = SG.pos_furthest - 1
SG.to_csv('siegel2014_pas_per_gene.tsv', sep='\t', index=False)
