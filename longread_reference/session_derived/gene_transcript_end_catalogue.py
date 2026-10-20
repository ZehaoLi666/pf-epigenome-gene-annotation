"""Session-derived computation of gene_transcript_end_catalogue.tsv.

Recovered from the record of the interactive analysis session in which gene_transcript_end_catalogue.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: tss_clusters.tsv.gz, utr_analysis.py, utr_classify.py, confirmed_genes_stage_TSS_switches.tsv, cluster_eval_original_model.tsv, hc_reference.bed.gz, gene_boundaries_by_stage.tsv.gz, tts_read_evidence_by_position.tsv, tss_error_detail.tsv.gz, gene_heterochromatin_groups.tsv, confusion_test_bins.tsv, cluster_eval_new_ens.tsv, confirmed_alternative_TSSs.tsv, gene_stage_confirmation.tsv, utr5_connectivity_per_dataset.tsv.gz, tss_read_evidence_by_group.tsv, transcript_stage_confirmation.tsv.gz, tts_error_detail.tsv.gz, stage_TPM_9stages.tsv, pr_curve_old_vs_new.tsv, utr5_connectivity_summary.tsv, loci_transcript_classes.tsv, end_transcription_contexts.tsv.gz, hc_model_evidence.tsv.gz
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import os, sys, numpy as np, pandas as pd, glob, importlib, shutil

# Load utr_analysis and utr_classify modules
shutil.copy(result('utr_analysis.py'), 'utr_analysis.py')
shutil.copy(result('utr_classify.py'), 'utr_classify.py')

import utr_analysis as UA

# Load data
D = UA.load('work/utr_stats_v68_fixed.tsv.gz')
M = pd.read_csv('work/tama_models.tsv.gz', sep='\t')
M = M[M.chrom.str.match(r'Pf3D7_\d\d_v3')]

TCL = pd.read_csv(result('tss_clusters.tsv.gz'), sep='\t')
cc = [x for x in TCL.columns if x.startswith('cage_')]
TCL['cage'] = TCL[cc].mean(1)

LAB = LAB_DIR

def read_fa(f):
    seqs, name, buf = {}, None, []
    for line in open(f):
        if line.startswith('>'):
            if name: seqs[name] = ''.join(buf)
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip().upper())
    seqs[name] = ''.join(buf); return seqs

G = read_fa(glob.glob(f'{LAB}/PlasmoDB/PlasmoDB-68/*Genome.fasta')[0])
gff = pd.read_csv(glob.glob(f'{LAB}/PlasmoDB/PlasmoDB-68/*.gff')[0], sep='\t', comment='#', header=None, low_memory=False,
                  names=['chrom','src','type','start','end','score','strand','phase','attr'])
gff = gff[gff.chrom.str.match(r'Pf3D7_\d\d_v3')]; gff['start'] -= 1
cds = gff[gff.type == 'CDS'].copy(); cds['gene'] = cds.attr.str.extract(r'Parent=([^;.]+)')[0]
genes_all = gff[gff.type.isin(['protein_coding_gene','ncRNA_gene','pseudogene'])].copy()
genes_all['gene'] = genes_all.attr.str.extract(r'ID=([^;]+)')[0]

lnc = pd.concat([pd.read_csv(f, sep='\t', header=None, usecols=[0,1,2,5], names=['chrom','start','end','strand'])
                 for f in glob.glob(f'{LAB}/3 lncRNA datasets/*.bed')])

import transcript_catalog as TC
importlib.reload(TC)

ovl_cds = TC.cds_mask_factory(cds)

tci = {k: (v.tss.values, v.cage.values) for k, v in TCL.groupby(['chrom','strand'])}
def cage_at(c, s, p, w=100):
    t = tci.get((c, s)); return 0.0 if t is None else float(t[1][np.abs(t[0] - p) <= w].sum())

Mi = TC.Index(M)

GN = D[(D.dataset == D.dataset.iloc[0]) & D.transcript.str.endswith('.1')].drop(columns=['dataset']).reset_index(drop=True)

flip = {'+': '-', '-': '+'}
gidx = {k: v for k, v in genes_all.groupby(['chrom', 'strand'])}

def other_gene_overlap(c, s, a, b, me):
    q = gidx.get((c, s))
    if q is None: return ''
    q = q[(q.start < b) & (q.end > a) & (q.gene != me)]; return ','.join(q.gene)

lncI = {k: v for k, v in lnc.groupby(['chrom','strand'])}
def lnc_ov(c, s, a, b):
    q = lncI.get((c, s)); return False if q is None else bool(((q.start < b) & (q.end > a)).any())

def indep_units(c, s, a, b, gene):
    q = Mi.query(c, s, a, b)
    if len(q) == 0: return q
    ov = (np.minimum(q.end, b) - np.maximum(q.start, a)) / (q.end - q.start)
    q = q[ov >= 0.5]
    keep = [not ovl_cds(c, s, *TC.exons(r)) for r in q.itertuples(index=False)]
    return q[keep]

def antisense_units(c, s, a, b, g=None):
    o = flip[s]; q = Mi.query(c, o, a, b)
    if len(q) == 0 or g is None: return q
    lo, hi = min(g.cds5, g.cds3), max(g.cds5, g.cds3) + 1
    keep = [(not ovl_cds(c, o, *TC.exons(r))) and not (r.start < hi and r.end > lo) for r in q.itertuples(index=False)]
    return q[keep]

def chimeric(g, a, b):
    lo, hi = min(g.cds5, g.cds3), max(g.cds5, g.cds3) + 1
    pa, pb = [], []
    for x, y in zip(a, b):
        if x < lo: pa.append(x); pb.append(min(y, lo))
        if y > hi: pa.append(max(x, hi)); pb.append(y)
    return ovl_cds(g.chrom, g.strand, pa, pb) if pa else False

rows = []
for g in GN.itertuples(index=False):
    r = dict(transcript=g.transcript, gene=g.gene, chrom=g.chrom, strand=g.strand, cds5=g.cds5, cds3=g.cds3,
             L5=g.L5 if g.has5 else 0, L3=g.L3 if g.has3 else 0)
    for side, L in ((5, r['L5']), (3, r['L3'])):
        cl, nds, rd, nch, ext = TC.ends_side(Mi, g, side, chimeric=chimeric); P = f'p{side}_'
        sup = [x for x in cl if TC.supported(x)]
        r[P + 'n_ds'] = nds; r[P + 'reads'] = rd; r[P + 'chimeric_reads'] = nch; r[P + 'n_clusters_supported'] = len(sup)
        a, b = TC.interval(g, side, 0, max(L, 1))
        r[P + 'other_gene_same_strand'] = other_gene_overlap(g.chrom, g.strand, a, b, g.gene) if L > 0 else ''
        r[P + 'other_gene_opposite'] = other_gene_overlap(g.chrom, flip[g.strand], a, b, g.gene) if L > 0 else ''
        r[P + 'lnc_sense'] = lnc_ov(g.chrom, g.strand, a, b) if L > 0 else False
        r[P + 'lnc_antisense'] = lnc_ov(g.chrom, flip[g.strand], a, b) if L > 0 else False
        if nds < 2 or not sup: r[P + 'cat'] = 'Insufficient long-read models'; continue
        major = max(sup, key=lambda x: x['reads']); distal = max(sup, key=lambda x: x['peak'])
        ud = ext; r[P + 'major_end'] = major['peak']; r[P + 'distal_end'] = ud; r[P + 'major_share'] = major['reads'] / rd
        alt = TC.alt_ends(cl, rd); r[P + 'n_alt'] = len(alt)
        if side == 5:
            cg = [cage_at(g.chrom, g.strand, (g.cds5 - x['peak']) if g.strand == '+' else (g.cds5 + x['peak'])) for x in alt]
            r['p5_n_alt_cage'] = int(sum(v >= 1 for v in cg)); r['p5_alt_positions'] = ','.join(str(x['peak']) for x in alt)
        if L == 0:
            r[P + 'cat'] = 'No annotated UTR; reads end at CDS (≤50 bp)' if major['peak'] <= 50 else 'No annotated UTR; long reads define a UTR'
            continue
        tol = max(200, 0.25 * L)
        r[P + 'cat'] = 'Long reads extend beyond annotated end' if ud > L + 200 else ('Annotated end supported' if ud >= L - tol else 'Annotated UTR longer than any supported end')
        r[P + 'excess_bp'] = max(0, L - ud)
        U = indep_units(g.chrom, g.strand, a, b, g.gene); okU, ndU, rdU = TC.unit_support(U, G)
        AS = antisense_units(g.chrom, g.strand, a, b, g); okA = bool(len(AS) and AS.dataset.nunique() >= 2)
        r[P + 'sense_unit'] = okU; r[P + 'sense_unit_ds'] = ndU; r[P + 'sense_unit_reads'] = rdU; r[P + 'antisense_unit'] = okA
        if okU:
            far = [((g.cds5 - x.start) if g.strand == '+' else (x.end - 1 - g.cds5)) if side == 5 else ((x.end - 1 - g.cds3) if g.strand == '+' else (g.cds3 - x.start)) for x in U.itertuples()]
            near = [((g.cds5 - x.end + 1) if g.strand == '+' else (x.start - g.cds5)) if side == 5 else ((x.start - g.cds3) if g.strand == '+' else (g.cds3 - x.end + 1)) for x in U.itertuples()]
            r[P + 'unit_far'] = int(np.median(far)); r[P + 'unit_near'] = int(np.median(near))
    rows.append(r)

C = pd.DataFrame(rows)

def classify_side(r, P):
    cat = r[P + 'cat']; L = r['L5' if P == 'p5_' else 'L3']
    if cat == 'Insufficient long-read models': return '0 No annotated UTR; insufficient long-read models' if L == 0 else '1 Insufficient long-read models'
    if cat.startswith('No annotated UTR; long'): return '8 No annotated UTR; long reads define a UTR'
    if cat.startswith('No annotated UTR; reads'): return '9 No annotated UTR; reads start/end at CDS'
    ud = r[P + 'distal_end']
    if cat == 'Annotated UTR longer than any supported end':
        if r[P + 'sense_unit'] and r[P + 'unit_far'] > ud + 100: return '2a Over-long: separate same-strand unit in the unsupported part'
        if r[P + 'other_gene_same_strand']: return '2b Over-long: overlaps a neighbouring same-strand gene'
        if r[P + 'antisense_unit']: return '2c Over-long: antisense unit only'
        return '2d Over-long: no other unit detected'
    if cat == 'Long reads extend beyond annotated end':
        return '6 Main end beyond annotated end (annotation too short)' if r[P + 'major_end'] > L + 200 else '7 Annotated end is main; minor ends further out'
    if r[P + 'sense_unit']: return '3 End supported; separate same-strand unit nested inside UTR'
    if r[P + 'n_alt'] >= 2: return '4 End supported; alternative ends'
    return '5 End supported; single end'

def split2d(r, P):
    c = r[P + 'class']
    if c.startswith('2d'): return '2d Over-long: catalogued lncRNA inside, no separate long-read unit' if r[P + 'lnc_sense'] else '2e Over-long: no other unit detected'
    return c

for P in ('p5_', 'p3_'): C[P + 'class'] = C.apply(lambda r: classify_side(r, P), axis=1)
for P in ('p5_', 'p3_'): C[P + 'class'] = C.apply(lambda r: split2d(r, P), axis=1)

C['p5_alt'] = C.p5_n_alt.fillna(0) >= 2
C['p5_alt_cage2'] = C.p5_n_alt_cage.fillna(0) >= 2
C['p5_readthrough_upstream_gene'] = C.p5_chimeric_reads >= 3
C['p3_readthrough_downstream_gene'] = C.p3_chimeric_reads >= 3
C['p5_opposite_gene_in_UTR'] = C.p5_other_gene_opposite.fillna('') != ''

C.to_csv('gene_transcript_end_catalogue.tsv', sep='\t', index=False)
