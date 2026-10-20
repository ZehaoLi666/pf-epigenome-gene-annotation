from pfconfig import TABLES_OUT_DIR  # paths: config/paths.yaml
import os, re, pandas as pd, numpy as np
from openpyxl import load_workbook
from openpyxl.styles import Font, Alignment
from openpyxl.utils import get_column_letter

OUT = TABLES_OUT_DIR
COORD = ("Coordinates: all genomic positions in this workbook are 0-based single-nucleotide positions on the P. falciparum 3D7 "
         "nuclear genome (PlasmoDB v48 = v68 = v71 assembly; chromosome names Pf3D7_01_v3 ... Pf3D7_14_v3). Add 1 to obtain the "
         "1-based coordinate shown in PlasmoDB/GFF3 browsers (checked: PF3D7_0100400 ATG = 50362 here, CDS start = 50363 in the GFF3). "
         "Interval columns explicitly marked 'BED-style' are half-open [start, end) with a 0-based start. TSS = first transcribed base; "
         "TTS = last transcribed base; ATG = first base of the start codon; stop = last base of the stop codon (strand-aware). "
         "Model-called ends are placed on a 20-bp grid (bin centre) by the calling script (gene_boundaries.py, RES = 20). "
         "Gene IDs are PlasmoDB IDs (PF3D7_xxxxxxx).")
B = Font(bold=True)

def sname(s):
    return re.sub(r'[\[\]\*\?/\\:]', '_', s)[:31]

def write_table(fname, title, intro, sheets, coord=True):
    """sheets: list of dict(name, df, desc{col:str}, source:str, note:str)."""
    os.makedirs(OUT, exist_ok=True)
    stem = fname.replace('.xlsx', '')
    tsvdir = os.path.join(OUT, stem + '_tsv'); os.makedirs(tsvdir, exist_ok=True)
    rows = [[title], [""]] + [[l] for l in intro] + [[""]]
    if coord: rows += [[COORD], [""]]
    rows.append(["Sheet", "Rows", "Source (result file name and/or script)", "Notes"])
    for s in sheets: rows.append([s['name'], len(s['df']), s['source'], s.get('note', '')])
    rows += [[""], ["Column dictionary"], ["Sheet", "Column", "Description"]]
    missing = []
    for s in sheets:
        for c in s['df'].columns:
            d = s['desc'].get(c)
            if d is None:
                d = 'NO DESCRIPTION (to verify)'; missing.append((s['name'], c))
            rows.append([s['name'], c, d])
    path = os.path.join(OUT, fname)
    with pd.ExcelWriter(path, engine='openpyxl') as xw:
        pd.DataFrame(rows).to_excel(xw, sheet_name='README', index=False, header=False)
        for s in sheets: s['df'].to_excel(xw, sheet_name=sname(s['name']), index=False)
    wb = load_workbook(path); ws = wb['README']
    for col, w in zip('ABCD', (34, 40, 110, 80)): ws.column_dimensions[col].width = w
    ws['A1'].font = Font(bold=True, size=13)
    for r in ws.iter_rows():
        if r[0].value in ('Sheet', 'Column dictionary'):
            for c in r: c.font = B
    for s in sheets:
        w = wb[sname(s['name'])]; w.freeze_panes = 'A2'
        for j, col in enumerate(s['df'].columns, 1):
            w.cell(1, j).font = B
            w.column_dimensions[get_column_letter(j)].width = min(max(10, len(str(col)) + 2), 45)
    wb.save(path)
    for s in sheets:
        s['df'].to_csv(os.path.join(tsvdir, f"{stem}_{sname(s['name'])}.tsv"), sep='\t', index=False)
    return path, missing
