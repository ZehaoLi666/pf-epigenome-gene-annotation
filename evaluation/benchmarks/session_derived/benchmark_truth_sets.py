"""Session-derived computation of benchmark_truth_sets.tsv.

Recovered from the record of the interactive analysis session in which benchmark_truth_sets.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: PlasmoDB-71_Pfalciparum3D7.gff, gene_boundaries_v3_all_genes.tsv, oof_top8_genes.tsv.gz, oof_top8_cage.tsv, oof_nochrom_cage.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import os, sys, glob, gzip, numpy as np, pandas as pd
import pyBigWig
from scipy import stats

os.chdir('.')

# Load GFF files
cols = ['chrom','src','type','start','end','score','strand','phase','attr']
def rd(p):
    g = pd.read_csv(p, sep='\t', comment='#', header=None, names=cols, low_memory=False)
    return g

g68 = rd(glob.glob(f'{LAB_DIR}/PlasmoDB/PlasmoDB-68/*.gff')[0])
g71 = rd(result('PlasmoDB-71_Pfalciparum3D7.gff'))
hdr71 = [l for l in open(result('PlasmoDB-71_Pfalciparum3D7.gff')) if l.startswith('##')][:30]
hdr68 = [l for l in open(glob.glob(f'{LAB_DIR}/PlasmoDB/PlasmoDB-68/*.gff')[0]) if l.startswith('##')][:30]

# Get description tables
from urllib.parse import unquote
def gdesc(g):
    q = g[g.type == 'protein_coding_gene']
    return pd.DataFrame({'gene_id': q.attr.str.extract(r'ID=([^;]+)')[0].values, 'desc': q.attr.str.extract(r'description=([^;]+)')[0].fillna('').map(unquote).values,
                         'name': q.attr.str.extract(r'Name=([^;]+)')[0].fillna('').values})
D = gdesc(g68).merge(gdesc(g71), on='gene_id', suffixes=('_68', '_71'))

# Load main gene boundaries file
B = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep='\t').merge(D[['gene_id', 'desc_71']], on='gene_id', how='left')

# Load top8 genes file
M8 = pd.read_csv(result('oof_top8_genes.tsv.gz'), sep='\t')[['gene_id', 'tss_win_lo', 'tss_win_hi', 'tts_win_lo', 'tts_win_hi', 'max_TSS_prob_window', 'max_TTS_prob_window']]
B = B.merge(M8, on='gene_id')

# Build gene truth table
T = pd.read_csv('work/models_v3.tsv.gz', sep='\t', low_memory=False) if os.path.exists('work/models_v3.tsv.gz') else pd.read_csv('v3/models_v3.tsv.gz', sep='\t', low_memory=False)
T = T[T.chrom.str.match(r'Pf3D7_\d\d_v3')]
S = ['ring', 'troph', 'schizont']; T['FL'] = T[[f'FL_{s}' for s in S]].sum(axis=1) + 0.01
T['tss'] = np.where(T.strand == '+', T.start, T.end - 1); T['tes'] = np.where(T.strand == '+', T.end - 1, T.start)
def dom(p, w):
    o = np.argsort(p); p, w = p[o], w[o]; cl = [[p[0], w[0], [p[0]], [w[0]]]]
    for x, y in zip(p[1:], w[1:]):
        if x - cl[-1][2][-1] <= 50: cl[-1][2].append(x); cl[-1][3].append(y); cl[-1][1] += y
        else: cl.append([x, y, [x], [y]])
    b = max(cl, key=lambda c: c[1]); return int(np.average(b[2], weights=b[3])), len(cl)
rows = []
F = T[T.full_cds & T.conf3_any]
for loc, d in F.groupby('locus'):
    t, n5 = dom(d.tss.values, d.FL.values); e, n3 = dom(d.tes.values, d.FL.values)
    rows.append(dict(gene_id=loc, lr_TSS=t, lr_TTS=e, lr_n_TSS_clusters=n5, lr_n_TTS_clusters=n3))
R = pd.DataFrame(rows)
grp = pd.Series('no long-read model', index=pd.Index(T.locus.unique()))
grp[T[~T.conf3_any].locus.unique()] = 'HC model, unconfirmed'
grp[T[(T.trunc5 | T.trunc3) & T.conf3_any].locus.unique()] = 'truncated-only'
grp[F.locus.unique()] = 'confirmed'
GT_truth = grp.rename('lr_group').rename_axis('gene_id').reset_index().merge(R, on='gene_id', how='left')

# Build GT with fold and window info
gff = g71.copy()
cds = gff[gff.type == 'CDS'].copy(); cds['gene_id'] = cds.attr.str.extract(r'Parent=([^;.]+)')[0]
env = cds.groupby('gene_id').agg(cds_lo=('start', 'min'), cds_hi=('end', 'max')); env['cds_lo'] -= 1
gene = gff[gff.type.isin(['protein_coding_gene'])].copy()
gene['gene_id'] = gene.attr.str.extract(r'ID=([^;]+)')[0]
G_gff = gene[['gene_id', 'chrom', 'strand', 'start', 'end']].merge(env, left_on='gene_id', right_index=True)
G_gff['start'] -= 1; G_gff['atg'] = np.where(G_gff.strand == '+', G_gff.cds_lo, G_gff.cds_hi - 1)
G_gff['stop'] = np.where(G_gff.strand == '+', G_gff.cds_hi - 1, G_gff.cds_lo)
G_gff['annot_tss'] = np.where(G_gff.strand == '+', G_gff.start, G_gff.end - 1)
G_gff['annot_tts'] = np.where(G_gff.strand == '+', G_gff.end - 1, G_gff.start)

rows2 = []
for (c, s), q in G_gff.groupby(['chrom', 'strand']):
    q = q.sort_values('cds_lo'); lo_, hi_ = q.cds_lo.values, q.cds_hi.values
    for i, r in enumerate(q.itertuples()):
        prev_hi = hi_[i - 1] if i > 0 else -10**9; next_lo = lo_[i + 1] if i + 1 < len(q) else 10**12
        if s == '+': w5 = (max(r.cds_lo - 3000, prev_hi + 50), r.cds_lo + 1); w3 = (r.cds_hi - 1, min(r.cds_hi + 2000, next_lo - 50))
        else: w5 = (r.cds_hi - 1, min(r.cds_hi + 3000, next_lo - 50)); w3 = (max(r.cds_lo - 2000, prev_hi + 50), r.cds_lo + 1)
        rows2.append(dict(gene_id=r.gene_id, chrom=c, strand=s, atg=int(r.atg), stop=int(r.stop),
                          annot_tss=int(r.annot_tss), annot_tts=int(r.annot_tts),
                          tss_win_lo=int(w5[0]), tss_win_hi=int(w5[1]), tts_win_lo=int(w3[0]), tts_win_hi=int(w3[1])))
GT = pd.DataFrame(rows2)
GT = GT.merge(GT_truth[['gene_id', 'lr_group', 'lr_TSS', 'lr_TTS']], on='gene_id', how='left')
fold = lambda c: 0 if int(c[6:8]) <= 5 else (1 if int(c[6:8]) <= 10 else 2)
GT['fold'] = GT.chrom.map(fold)

# Load CAGE bigwigs
CB = {(s, st): pyBigWig.open(f"{LAB_DIR}/lncRNA project tracks/cage_tss/CTSS_{s}_{st}_CPM.bw") for s in ('ring', 'trophozoite', 'schizont') for st in ('pos', 'neg')}
def track(c, strand, a, b):
    k = 'pos' if strand == '+' else 'neg'
    return np.max([np.nan_to_num(np.array(CB[(s, k)].values(c, max(0, a), b), float)) for s in ('ring', 'trophozoite', 'schizont')], axis=0)

rng = np.random.default_rng(0)
def adj(r, minw=5):
    lo, hi = int(r.tss_win_lo), int(r.tss_win_hi)
    if hi - lo < 50: return np.nan
    v = track(r.chrom, r.strand, lo - 100, hi + 100); sm = np.convolve(v, np.ones(201), 'valid')
    if sm.max() < minw: return np.nan
    c0 = int(np.argmax(sm)); seg = v[c0:c0 + 201]; return lo - 100 + c0 + int(np.argmax(seg))
T_ad = pd.Series({r.gene_id: adj(r) for r in GT.itertuples()}).dropna()

# Load Shaw 2021 5CAPture data
f_shaw = 'bench/indep/GSE103036_combine_enrichment_score.txt.gz'
rows_shaw = []
with gzip.open(f_shaw, 'rt') as fh:
    next(fh)
    for line in fh:
        x = line.rstrip('\n').split('\t')
        try: p = float(x[-1])
        except ValueError: continue
        if p >= 0.01: continue
        body = '\t'.join(x[3:-1]).split('\t*'); a = b = 0.0; ns = 0
        for g in body:
            t_v = [y for y in g.split('\t') if y]
            if t_v and t_v[0] != 'n/a': a += float(t_v[0]); b += float(t_v[1]); ns += 1
        rows_shaw.append((x[2], x[1], int(x[0]), p, a, b, ns))
C5 = pd.DataFrame(rows_shaw, columns=['chrom', 'strand', 'pos', 'p', 'd1', 'd2', 'n_samples'])

def cluster(df, gap=20, wcol='d1'):
    out = []
    for (c, s), q in df.sort_values('pos').groupby(['chrom', 'strand']):
        p = q.pos.values; w = q[wcol].values; br = np.r_[0, np.where(np.diff(p) > gap)[0] + 1, len(p)]
        for i, j in zip(br[:-1], br[1:]):
            k = i + int(np.argmax(w[i:j])); out.append((c, s, int(p[i]), int(p[j - 1]), int(p[k]), float(w[i:j].sum()), j - i))
    return pd.DataFrame(out, columns=['chrom', 'strand', 'start', 'end', 'dominant', 'score', 'n_pos'])
CAP = cluster(C5); CAP['dominant'] -= 1
CAP = CAP[CAP.chrom.str.match(r'Pf3D7_\d\d_v3')]

def truth_from_points(P, poscol, wcol, minw, gap=50):
    idx = {k: v.sort_values(poscol) for k, v in P.groupby(['chrom', 'strand'])}; out = {}
    for r in GT.itertuples():
        q = idx.get((r.chrom, r.strand))
        if q is None: continue
        p = q[poscol].values; i, j = np.searchsorted(p, r.tss_win_lo), np.searchsorted(p, r.tss_win_hi, side='right')
        if j <= i: continue
        pp, ww = p[i:j], q[wcol].values[i:j]; br = np.r_[0, np.where(np.diff(pp) > gap)[0] + 1, len(pp)]
        best = max(((ww[a_:b_].sum(), pp[a_ + int(np.argmax(ww[a_:b_]))]) for a_, b_ in zip(br[:-1], br[1:])))
        if best[0] >= minw: out[r.gene_id] = best[1]
    return pd.Series(out)

T_sh = truth_from_points(CAP.rename(columns={'dominant': 'pos'}), 'pos', 'score', minw=100)

# Load Chappell 2020 5UTR-seq data
import openpyxl
f_chap = glob.glob('bench/indep/chappell/*MOESM10_ESM.xlsx')[0]
wb = openpyxl.load_workbook(f_chap, read_only=True); parts = []
for sh in wb.sheetnames:
    ws = wb[sh]; rows_c = list(ws.iter_rows(values_only=True)); h = [i for i, r in enumerate(rows_c) if r and r[0] == 'Chromosome'][0]
    d = pd.DataFrame(rows_c[h + 1:], columns=rows_c[h]); d['tp'] = sh; parts.append(d)
S9 = pd.concat(parts).dropna(subset=['Chromosome']); S9 = S9[S9.Chromosome.astype(str).str.match(r'Pf3D7_\d\d_v3')]
for c_ in ('Quantile_10', 'Quantile_90', 'Dominant_tag'): S9[c_] = S9[c_].astype(int)
S9['Quantile_TPM'] = S9.Quantile_TPM.astype(float)
S9p = S9.assign(chrom=S9.Chromosome, strand=S9.Strand, pos=S9.Dominant_tag - 1, w=S9.Quantile_TPM)
T_ch = truth_from_points(S9p, 'pos', 'w', minw=5)

# Load cov 3UTR and 5UTR from Chappell MOESM files
def gff_tab(m):
    f = glob.glob(f'bench/indep/chappell/*MOESM{m}_ESM.xlsx')[0]; wb = openpyxl.load_workbook(f, read_only=True); ws = wb[wb.sheetnames[0]]
    rows_g = [r for r in ws.iter_rows(values_only=True)]; hdr = [r[0] for r in rows_g[:3]]
    d = pd.DataFrame([r[:9] for r in rows_g if r[0] and str(r[0]).startswith('Pf3D7_') and r[3] is not None], columns=['chrom', 'src', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attr'])
    return hdr, d

h2, U5 = gff_tab(3); U5['gene_id'] = U5.attr.str.extract(r'Parent=([^;.]+)')[0]; U5['start'] = U5.start.astype(int) - 1; U5['end'] = U5.end.astype(int)
U5['tss'] = np.where(U5.strand == '+', U5.start, U5.end - 1)
h4, U3 = gff_tab(5); U3['gene_id'] = U3.attr.str.extract(r'Parent=([^;.]+)')[0]; U3['start'] = U3.start.astype(int) - 1; U3['end'] = U3.end.astype(int)
U3['tts'] = np.where(U3.strand == '+', U3.end - 1, U3.start)

# Build truth table
TR_ = GT[['gene_id']].copy()
TR_['LongRead__TSS'] = GT.gene_id.map(GT_truth.set_index('gene_id').lr_TSS)
TR_['LongRead__TTS'] = GT.gene_id.map(GT_truth.set_index('gene_id').lr_TTS)
TR_['Adjalley2016_capseq__TSS'] = TR_.gene_id.map(T_ad)
TR_['Shaw2021_5CAPture__TSS'] = TR_.gene_id.map(T_sh)
TR_['Chappell2020_5UTRseq__TSS'] = TR_.gene_id.map(T_ch)

# Compute consensus
def consensus(r, tol=100):
    v = [r[f'{k}__TSS'] for k in ('Adjalley2016_capseq', 'Shaw2021_5CAPture', 'Chappell2020_5UTRseq') if pd.notna(r[f'{k}__TSS'])]
    best = None
    for a_ in v:
        s_v = [b_ for b_ in v if abs(b_ - a_) <= tol]
        if len(s_v) >= 2 and (best is None or len(s_v) > len(best)): best = s_v
    return float(np.median(best)) if best else np.nan
TR_['Consensus_5end__TSS'] = TR_.apply(consensus, axis=1)

TRb = TR_.copy()

TRb.to_csv('benchmark_truth_sets.tsv', sep='\t', index=False)
