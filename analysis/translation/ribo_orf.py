"""Ribo-seq evidence for ORFs (Caro 2014 footprints). Inputs: sites/<run>.tsv.gz (chrom, read strand, 5' pos, length, count),
scripts/cds_exons_v68.tsv, scripts/query_orfs_genomic.tsv. Steps:
1. Strand check (footprints on CDS sense vs antisense).
2. Per read length 25-35: P-site offset = mode distance from 5' end to annotated start codon (reads whose 5' end lies 6-20 nt upstream);
   keep lengths whose CDS P-sites show >= 45% in frame 0.
3. Per query ORF: sense P-site count, frame distribution relative to the ORF, RPKM (per million CDS P-sites), binomial test frame 0 > 1/3;
   plus antisense count. Pooled over all 8 runs and per stage."""
from pfconfig import HPCC_BIGDATA  # paths: config/paths.yaml
import glob, gzip, collections, numpy as np, pandas as pd
from scipy.stats import binomtest
W = f'{HPCC_BIGDATA}/riboseq_caro2014'
RUNS = {l.split('\t')[0]: l.split('\t')[1] for l in open(f'{W}/scripts/runs.tsv').read().strip().split('\n')}
EX = pd.read_csv(f'{W}/scripts/cds_exons_v68.tsv', sep='\t')
# CDS position -> (frame) arrays per chrom/strand: build per gene transcript-coordinate map
cdsmap = {}   # (chrom, strand) -> dict genomic pos -> frame
starts = []
for g, q in EX.groupby('gene_id'):
    c, s = q.chrom.iloc[0], q.strand.iloc[0]; iv = sorted(zip(q.start, q.end)); pos = np.concatenate([np.arange(a, b) for a, b in iv])
    if s == '-': pos = pos[::-1]
    d = cdsmap.setdefault((c, s), {})
    for k, p in enumerate(pos): d[int(p)] = k % 3
    starts.append((c, s, int(pos[0])))
SITES = {}
for f in sorted(glob.glob(f'{W}/sites/*.tsv.gz')):
    run = f.split('/')[-1].split('.')[0]; SITES[run] = pd.read_csv(f, sep='\t', header=None, names=['chrom', 'strand', 'pos', 'len', 'n'])
A = pd.concat(SITES.values()).groupby(['chrom', 'strand', 'pos', 'len']).n.sum().reset_index()
# 1. strand check
def in_cds(df, flipstrand=False):
    m = np.zeros(len(df), bool)
    for (c, s), q in df.groupby(['chrom', 'strand']):
        s2 = {'+': '-', '-': '+'}[s] if flipstrand else s; d = cdsmap.get((c, s2), {})
        m[q.index.values] = [p in d for p in q.pos.values]
    return m
A = A.reset_index(drop=True); sense = in_cds(A); anti = in_cds(A, True)
print('reads in CDS sense / antisense:', int(A.n[sense].sum()), int(A.n[anti].sum()), flush=True)
# 2. offsets
OFF = {}
for L in range(25, 36):
    q = A[A.len == L]; D = dict(zip(zip(q.chrom, q.strand, q.pos), q.n)); h = collections.Counter()
    for c, s, p in starts:
        for d in range(6, 21): h[d] += D.get((c, s, p - d if s == '+' else p + d), 0)
    if sum(h.values()): OFF[L] = h.most_common(1)[0][0]
print('offsets', OFF, flush=True)
FR5 = collections.defaultdict(lambda: np.zeros(3))
for (c, s), q in A.groupby(['chrom', 'strand']):
    d = cdsmap.get((c, s), {})
    for p, L, n in zip(q.pos.values, q.len.values, q.n.values):
        f = d.get(int(p))
        if f is not None: FR5[L][f] += n
print('5prime frame dist', {int(L): np.round(FR5[L] / max(FR5[L].sum(), 1), 3).tolist() for L in sorted(FR5)}, flush=True)
def psites(df):
    o = df.len.map(OFF); df = df[o.notna()].copy(); o = o[o.notna()].astype(int)
    df['p'] = np.where(df.strand == '+', df.pos + o, df.pos - o); return df
P = psites(A)
fr = collections.defaultdict(lambda: np.zeros(3))
for (c, s), q in P.groupby(['chrom', 'strand']):
    d = cdsmap.get((c, s), {})
    for p, L, n in zip(q.p.values, q.len.values, q.n.values):
        f = d.get(int(p))
        if f is not None: fr[L][f] += n
good = [L for L in OFF if fr[L].sum() > 0 and fr[L][0] / fr[L].sum() >= 0.45] or [L for L in OFF if fr[L].sum() > 0]
print('frame0 by length', {L: round(fr[L][0] / max(fr[L].sum(), 1), 3) for L in sorted(fr)}, 'kept', good, flush=True)
P = P[P.len.isin(good)]; tot_cds = sum(fr[L].sum() for L in good) / 1e6
IDX = {k: (q.p.values, q.n.values) for k, q in P.groupby(['chrom', 'strand'])}
Q = pd.read_csv(f'{W}/scripts/query_orfs_genomic.tsv', sep='\t'); rows = []
flip = {'+': '-', '-': '+'}
for r in Q.itertuples():
    d = dict(qid=r.qid)
    if pd.isna(r.orf_g0) or r.orf_aa < 10: rows.append(d); continue
    a, b = int(r.orf_g0), int(r.orf_g1); pp, nn = IDX.get((r.chrom, r.strand), (np.array([]), np.array([])))
    m = (pp >= a) & (pp < b); rel = (pp[m] - a) if r.strand == '+' else (b - 1 - pp[m]); f3 = np.bincount((rel % 3).astype(int), weights=nn[m], minlength=3)
    k = int(f3[0]); n = int(f3.sum()); pa, na = IDX.get((r.chrom, flip[r.strand]), (np.array([]), np.array([]))); ma = (pa >= a) & (pa < b)
    d.update(psites=n, frame0=k, frame0_frac=k / n if n else np.nan, p_frame=binomtest(k, n, 1 / 3, alternative='greater').pvalue if n else np.nan,
             rpkm=n / ((b - a) / 1000) / tot_cds, antisense_psites=int(na[ma].sum()))
    mu = (pp >= r.start) & (pp < r.end); n_unit = float(nn[mu].sum()); L_out = (r.end - r.start) - (b - a)
    d.update(unit_psites=n_unit, outside_psites=n_unit - n, outside_len=L_out, rpkm_outside=(n_unit - n) / max(L_out / 1000, 0.05) / tot_cds)
    rows.append(d)
R = pd.DataFrame(rows); R.to_csv(f'{W}/ribo_orf_metrics.tsv', sep='\t', index=False)
pd.DataFrame([dict(len=L, offset=OFF.get(L), frame0=fr[L][0] / max(fr[L].sum(), 1), n=fr[L].sum(), kept=L in good) for L in range(25, 36)]).to_csv(f'{W}/ribo_offsets.tsv', sep='\t', index=False)
print('DONE', len(R))
