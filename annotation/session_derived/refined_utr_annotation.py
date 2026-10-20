"""Session-derived computation of refined_utr_annotation.gff3.

Recovered from the record of the interactive analysis session in which refined_utr_annotation.gff3 was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_boundaries_v3_all_genes.tsv, genes_without_utr_prediction.tsv, common.py, tracks.py
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import pandas as pd
import numpy as np
import re
import os

# helpers: annotation/common.py and annotation/tracks.py
# Load the tracks helper
exec(open(result('common.py')).read())
exec(open(result('tracks.py')).read())

GENE_TABLE = result('gene_boundaries_v3_all_genes.tsv')
NOPRED_TABLE = result('genes_without_utr_prediction.tsv')

G = pd.read_csv(GENE_TABLE, sep='\t')
NW = pd.read_csv(NOPRED_TABLE, sep='\t')[['gene_id', 'no_TSS_reason', 'no_TTS_reason', 'het_group', 'IDC_TPM']]
gff = load_gff()

GROUP = {'confirmed': 'reference_stage_confirmed', 'truncated-only': 'truncated_long_reads_only',
         'HC model, unconfirmed': 'supported_model_not_stage_confirmed', 'no HC model': 'no_supported_model'}

def introns(blocks):
    b = sorted(blocks); return [(b[i][1] + 1, b[i + 1][0] - 1) for i in range(len(b) - 1) if b[i + 1][0] - b[i][1] > 1]

cds = gff[gff.type == 'CDS'].copy(); cds['tx'] = cds.Parent.str.split(',').str[0]
cdsr = cds.groupby('tx').agg(s=('start', 'min'), e=('end', 'max'), gene=('gene', 'first'))
utr = gff[gff.type.isin(['five_prime_UTR', 'three_prime_UTR'])].copy(); utr['tx'] = utr.Parent.str.split(',').str[0]
CDSB = {tx: d.sort_values('start')[['start', 'end', 'phase']].values.tolist() for tx, d in cds.groupby('tx')}
UTRB = {(tx, ty): sorted(d[['start', 'end']].values.tolist()) for (tx, ty), d in utr.groupby(['tx', 'type'])}
EXB = {tx: d.sort_values('start')[['start', 'end']].values.tolist() for tx, d in gff[gff.type == 'exon'].assign(tx=lambda x: x.Parent.str.split(',').str[0]).groupby('tx')}
TXI = {tx: introns(b) for tx, b in EXB.items()}

def pick_tx(r):
    a1 = r.atg + 1; s1 = r.stop + 1
    lo, hi = (a1, s1) if r.strand == '+' else (s1, a1)
    c = cdsr[(cdsr.gene == r.gene_id) & (cdsr.s == lo) & (cdsr.e == hi)]
    assert len(c) >= 1, r.gene_id
    return sorted(c.index)[0], lo, hi

def utr_blocks(lo, hi, v71_introns, end_is_lo):
    """UTR span [lo,hi] (1-based closed) minus PlasmoDB v71 introns inside the span. An intron may abut the CDS side
    (intron directly before the ATG / after the stop) but must not touch the transcript end (TSS/TTS)."""
    keep = [(a, b) for a, b in v71_introns if ((a > lo and b <= hi) if end_is_lo else (a >= lo and b < hi))]
    out, cur = [], lo
    for a, b in keep: out.append([cur, a - 1]); cur = b + 1
    out.append([cur, hi]); return [b for b in out if b[1] >= b[0]], len(keep)

def fmt_alt(x):
    if pd.isna(x) or str(x).strip() == '': return ''
    return ','.join(str(int(float(v)) + 1) for v in str(x).split(','))

def esc(v): return re.sub(r'([;=&,])', lambda m: '%{:02X}'.format(ord(m.group(1))), str(v))

rows, lines, qc = [], ['##gff-version 3'], dict(utr5_intron_retained=0, utr3_intron_retained=0)
for r in G.itertuples():
    tx, clo, chi = pick_tx(r); st = r.strand; plus = st == '+'
    rec = dict(gene_id=r.gene_id, transcript_id=tx, chrom=r.chrom, strand=st, gene_group=GROUP[r.lr_group],
               cds_start=clo, cds_end=chi)
    feats = []
    for end, ty in (('5', 'five_prime_UTR'), ('3', 'three_prime_UTR')):
        fin = r.final_TSS if end == '5' else r.final_TTS
        ann_len = r.utr5_annot if end == '5' else r.utr3_annot
        v71b = UTRB.get((tx, ty), [])
        v71i = [iv for iv in TXI[tx] if (iv[1] < clo or iv[0] > chi)]   # v71 introns outside the CDS (UTR introns)
        if pd.notna(fin):
            ev = 'long-read' if r.lr_group == 'confirmed' else 'model'
            p1 = int(fin) + 1
            # UTR span on the genome
            if end == '5': lo, hi = (p1, clo - 1) if plus else (chi + 1, p1)
            else: lo, hi = (chi + 1, p1) if plus else (p1, clo - 1)
            end_is_lo = (end == '5') == plus
            if hi >= lo: blocks, nk = utr_blocks(lo, hi, v71i, end_is_lo)
            else: blocks, nk = [], 0
            qc[f'utr{end}_intron_retained'] += nk > 0
        elif ann_len > 0:
            ev = 'PlasmoDB_v71'; blocks = v71b
            assert len(blocks) > 0, (r.gene_id, ty)
            p1 = (min(b[0] for b in blocks) if plus else max(b[1] for b in blocks)) if end == '5' else \
                 (max(b[1] for b in blocks) if plus else min(b[0] for b in blocks))
        else:
            ev = 'none'; blocks = []; p1 = (clo if plus else chi) if end == '5' else (chi if plus else clo)
        L = sum(b[1] - b[0] + 1 for b in blocks); span = abs(p1 - ((clo if plus else chi) if end == '5' else (chi if plus else clo)))
        k = 'TSS' if end == '5' else 'TTS'
        rec[k] = p1; rec[f'{k}_evidence'] = ev; rec[f'utr{end}_genomic_len'] = span if ev != 'none' else 0
        rec[f'utr{end}_exonic_len'] = L; rec[f'utr{end}_n_exons'] = len(blocks)
        prob = getattr(r, f'model_{k}_prob'); rec[f'{k}_model_prob'] = round(prob, 3) if ev == 'model' and pd.notna(prob) else np.nan
        rec[f'{k}_alternatives'] = fmt_alt(getattr(r, f'model_alt_{k}')) if ev == 'model' else ''
        feats += [(ty, b, ev) for b in blocks]
    rec['TSS_CAGE_percentile'] = round(r.CAGE_pct_model_TSS, 3) if rec['TSS_evidence'] == 'model' and pd.notna(r.CAGE_pct_model_TSS) else np.nan
    rec['v71_TSS'] = int(r.annot_tss) + 1; rec['v71_TTS'] = int(r.annot_tts) + 1
    rec['v71_utr5_len'] = int(r.utr5_annot); rec['v71_utr3_len'] = int(r.utr3_annot)
    rows.append(rec)
    # ---- GFF3 lines ----
    cb = CDSB[tx]
    allb = [b for _, b, _ in feats] + [[a, b] for a, b, _ in cb]
    ms, me = min(b[0] for b in allb), max(b[1] for b in allb)
    gattr = f"ID={r.gene_id};gene_group={GROUP[r.lr_group]};description={esc(str(r.desc).replace('%2C', ','))}"
    mattr = [f"ID={tx}", f"Parent={r.gene_id}", f"TSS_evidence={rec['TSS_evidence']}", f"TTS_evidence={rec['TTS_evidence']}"]
    for k in ('TSS', 'TTS'):
        if pd.notna(rec[f'{k}_model_prob']): mattr.append(f"{k}_model_prob={rec[f'{k}_model_prob']}")
        if rec[f'{k}_alternatives']: mattr.append(f"{k}_alternatives={rec[f'{k}_alternatives']}")
    if pd.notna(rec['TSS_CAGE_percentile']): mattr.append(f"TSS_CAGE_percentile={rec['TSS_CAGE_percentile']}")
    mattr += [f"PlasmoDB_v71_TSS={rec['v71_TSS']}", f"PlasmoDB_v71_TTS={rec['v71_TTS']}"]
    lines.append('\t'.join(map(str, [r.chrom, 'refinedUTR', 'protein_coding_gene', ms, me, '.', st, '.', gattr])))
    lines.append('\t'.join(map(str, [r.chrom, 'refinedUTR', 'mRNA', ms, me, '.', st, '.', ';'.join(mattr)])))
    # exons = union of UTR and CDS blocks, merging abutting blocks
    ex = sorted(allb); mer = [list(ex[0])]
    for a, b in ex[1:]:
        if a <= mer[-1][1] + 1: mer[-1][1] = max(mer[-1][1], b)
        else: mer.append([a, b])
    if not plus: mer = mer[::-1]
    out = []
    for i, (a, b) in enumerate(mer): out.append((a, 'exon', a, b, '.', f"ID={tx}-E{i+1};Parent={tx}"))
    n5 = n3 = 0
    for ty, (a, b), ev in feats:
        if ty == 'five_prime_UTR': n5 += 1; i = n5
        else: n3 += 1; i = n3
        out.append((a, ty, a, b, '.', f"ID={tx}-{'5UTR' if ty[0]=='f' else '3UTR'}{i};Parent={tx};evidence={ev}"))
    for i, (a, b, ph) in enumerate(cb if plus else cb[::-1]):
        out.append((a, 'CDS', a, b, ph, f"ID={tx}-CDS{i+1};Parent={tx};evidence=PlasmoDB_v71"))
    order = {'exon': 0, 'five_prime_UTR': 1, 'CDS': 2, 'three_prime_UTR': 3}
    for _, ty, a, b, ph, at in sorted(out, key=lambda x: (x[0], order[x[1]])):
        lines.append('\t'.join(map(str, [r.chrom, 'refinedUTR', ty, a, b, '.', st, ph, at])))

T = pd.DataFrame(rows).merge(NW, on='gene_id', how='left')

open('refined_utr_annotation.gff3', 'w').write('\n'.join(lines) + '\n')
