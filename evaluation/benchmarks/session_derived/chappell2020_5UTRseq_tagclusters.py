"""Session-derived computation of chappell2020_5UTRseq_tagclusters.tsv.

Recovered from the record of the interactive analysis session in which chappell2020_5UTRseq_tagclusters.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: 12864_2020_6787_MOESM10_ESM.xlsx
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import openpyxl
import numpy as np
import pandas as pd

f = result('12864_2020_6787_MOESM10_ESM.xlsx')
wb = openpyxl.load_workbook(f, read_only=True)
parts = []
for sh in wb.sheetnames:
    ws = wb[sh]
    rows = list(ws.iter_rows(values_only=True))
    h = [i for i, r in enumerate(rows) if r and r[0] == 'Chromosome'][0]
    d = pd.DataFrame(rows[h + 1:], columns=rows[h])
    d['tp'] = sh
    parts.append(d)
S9 = pd.concat(parts).dropna(subset=['Chromosome'])
S9 = S9[S9.Chromosome.astype(str).str.match(r'Pf3D7_\d\d_v3')]
for c_ in ('Quantile_10', 'Quantile_90', 'Dominant_tag'):
    S9[c_] = S9[c_].astype(int)
S9['Quantile_TPM'] = S9.Quantile_TPM.astype(float)
S9.to_csv('chappell2020_5UTRseq_tagclusters.tsv', sep='\t', index=False)
