"""Labels v3 (CAGE-free). 20-bp bins, pooled stages. Outputs work/data/Lml3_pool.npz (Y, W) and work/data/ml3_end_positives.tsv.
Gene set: loci with >=1 v3-confirmed full-CDS model (gene_set_v3.py); truncated-only v3-confirmed loci contribute body labels and supported ends.
TSS clusters (<=50 bp) of full-CDS models: positive, weight max(0.3, read share) x (1.0 if cap-G supported else 0.5); cap-G = >=2 cDNA reads
  with a non-templated 5' G and >=20% of reads (read_end_qc.tsv, +-20 bp).
TTS clusters of full-CDS models: positive if downstream 20 nt < 60% A and not flagged internal-no-tail; else masked.
Ends of truncated models that fall outside full-CDS end clusters: positive only if supported (TSS cap-G rule; TTS >=2 tailed reads, >=5%,
  downstream A < 60%); weight 0.5 if internal to the CDS envelope, 1.0 otherwise; unsupported -> masked.
Flank: same-strand gaps (gaps_v3.tsv) classified from long reads only (gap_stats_v3.tsv): spanned -> masked; A/B/C -> 1.0/0.6/0.3."""
from pfconfig import LAB_DIR  # paths: config/paths.yaml
import sys, glob, numpy as np, pandas as pd
sys.path[:0] = ['work']
RES, WIN, TRANS, EDGE = 20, 100, 300, 300
LAB = LAB_DIR; S = ['ring', 'troph', 'schizont']
T = pd.read_csv('v3/models_v3.tsv.gz', sep='\t', low_memory=False); T = T[T.chrom.str.match(r'Pf3D7_\d\d_v3')]
G = {}; name = None; buf = []
for line in open(glob.glob(f'{LAB}/PlasmoDB/PlasmoDB-48/*Genome.fasta')[0]):
    if line[0] == '>':
        if name: G[name] = ''.join(buf)
        name = line[1:].split()[0].split('|')[0].strip(); buf = []
    else: buf.append(line.strip().upper())
G[name] = ''.join(buf); comp = str.maketrans('ACGT', 'TGCA')
CH = sorted(c for c in G if c.startswith('Pf3D7_') and c[6:8].isdigit()); nb = {c: len(G[c]) // RES for c in CH}
def downA(c, s, p, n=20):
    w = G[c][p + 1:p + 1 + n] if s == '+' else G[c][max(0, p - n):p].translate(comp)[::-1]
    return w.count('A') / max(len(w), 1)
Q = pd.read_csv('work/data/read_end_qc.tsv', sep='\t'); Q['G1'] = Q.pacbio_G1 + Q.ont_G1; Q['n'] = Q.pacbio_n + Q.ont_n; Q['tail'] = Q.pacbio_tail + Q.ont_tail
Qi = {k: (v['pos'].values, v['G1'].values, v['n'].values, v['tail'].values) for k, v in Q.sort_values('pos').groupby(['chrom', 'strand', 'end_type'])}
def qc(c, s, et, p, tol=20):
    t = Qi.get((c, s, et))
    if t is None: return 0, 0, 0
    j = np.searchsorted(t[0], p); cand = [k for k in (j - 1, j) if 0 <= k < len(t[0]) and abs(t[0][k] - p) <= tol]
    if not cand: return 0, 0, 0
    k = min(cand, key=lambda k: abs(t[0][k] - p)); return t[1][k], t[2][k], t[3][k]
cap_ok = lambda c, s, p: (lambda g, n, _: g >= 2 and g / max(n, 1) >= 0.2)(*qc(c, s, 'TSS', p))
tail_ok = lambda c, s, p: (lambda _, n, t: t >= 2 and t / max(n, 1) >= 0.05)(*qc(c, s, 'TTS', p)) and downA(c, s, p) < 0.6
FLG = pd.read_csv('work/data/end_flags_min.tsv', sep='\t'); FLG = FLG[(FLG.end_type == 'TTS') & FLG.suspect]   # TTS flag = inside gene, no tail (no CAGE)
fl = {k: np.sort(v.pos.values) for k, v in FLG.groupby(['chrom', 'strand'])}
def tts_suspect(c, s, p):
    a = fl.get((c, s)); return a is not None and len(a) and np.abs(a[np.clip(np.searchsorted(a, p), 0, len(a) - 1)] - p) <= 50
conf = T.conf3_any; F = T[T.full_cds & conf].copy(); TR = T[(T.trunc5 | T.trunc3) & conf].copy()
geneset = set(F.locus); tr_only = set(TR.locus) - geneset
Y = {c: np.zeros((2, 3, n), np.uint8) for c, n in nb.items()}; W = {c: np.zeros((2, 3, n), np.float32) for c, n in nb.items()}
si = lambda s: 0 if s == '+' else 1; Bn = lambda p: int(p) // RES
def setr(c, s, ch, lo, hi, y=None, w=None, maxw=False):
    lo, hi = max(0, Bn(lo)), min(nb[c], Bn(hi - 1) + 1)
    if hi <= lo: return
    if y is not None: Y[c][si(s), ch, lo:hi] = y
    if w is not None: W[c][si(s), ch, lo:hi] = np.maximum(W[c][si(s), ch, lo:hi], w) if maxw else w
for d in (F, TR): d['tss'] = np.where(d.strand == '+', d.start, d.end - 1); d['tes'] = np.where(d.strand == '+', d.end - 1, d.start); d['FL'] = d[[f'FL_{s}' for s in S]].sum(axis=1)
pos, masked = [], []; st = dict(tss_cap=0, tss_nocap=0, tts_pos=0, tts_masked=0, trunc_added_tss=0, trunc_added_tts=0, trunc_masked=0)
for loc, d in F.groupby('locus'):
    c, s = d.chrom.iloc[0], d.strand.iloc[0]
    for r in d.itertuples(): setr(c, s, 2, r.start, r.end, y=1, w=1.0)
    for col, et in (('tss', 'TSS'), ('tes', 'TTS')):
        cl = []
        for r in d.sort_values(col).itertuples():
            p = getattr(r, col)
            if cl and p - cl[-1]['hi'] <= 50: cl[-1]['hi'] = p; cl[-1]['ps'].append(p); cl[-1]['w'].append(r.FL + 0.01)
            else: cl.append(dict(hi=p, ps=[p], w=[r.FL + 0.01]))
        mx = max(sum(k['w']) for k in cl)
        for k in cl:
            p = int(np.average(k['ps'], weights=k['w'])); share = max(0.3, sum(k['w']) / mx)
            if et == 'TSS':
                ok = cap_ok(c, s, p); pos.append((c, s, et, p, share * (1.0 if ok else 0.5), loc, 'terminal', bool(ok))); st['tss_cap' if ok else 'tss_nocap'] += 1
            elif downA(c, s, p) < 0.6 and not tts_suspect(c, s, p): pos.append((c, s, et, p, share, loc, 'terminal', bool(tail_ok(c, s, p)))); st['tts_pos'] += 1
            else: masked.append((c, s, et, p)); st['tts_masked'] += 1
POS = pd.DataFrame(pos, columns=['chrom', 'strand', 'end_type', 'pos', 'weight', 'locus', 'kind', 'end_evidence'])
known = {k: np.sort(v.pos.values) for k, v in POS.groupby(['chrom', 'strand', 'end_type'])}
def near_known(c, s, et, p, d=50):
    a = known.get((c, s, et)); return a is not None and len(a) and np.abs(a[np.clip(np.searchsorted(a, p), 0, len(a) - 1)] - p) <= d
add = []; seen = set()
for r in TR.itertuples():
    if r.locus in tr_only: setr(r.chrom, r.strand, 2, r.start, r.end, y=1, w=1.0)
    for et, p in (('TSS', r.tss), ('TTS', r.tes)):
        key = (r.chrom, r.strand, et, int(p) // 50)
        if key in seen or near_known(r.chrom, r.strand, et, p): continue
        seen.add(key); internal = (r.cds_lo <= p < r.cds_hi) if not np.isnan(r.cds_lo) else False
        ok = cap_ok(r.chrom, r.strand, p) if et == 'TSS' else tail_ok(r.chrom, r.strand, p)
        if ok: add.append((r.chrom, r.strand, et, int(p), 0.5 if internal else 1.0, r.locus, 'internal' if internal else 'truncated-model end', True)); st['trunc_added_tss' if et == 'TSS' else 'trunc_added_tts'] += 1
        else: masked.append((r.chrom, r.strand, et, int(p))); st['trunc_masked'] += 1
POS = pd.concat([POS, pd.DataFrame(add, columns=POS.columns)], ignore_index=True)
for ch in (0, 1):
    for c in Y: W[c][:, ch] = np.where(Y[c][:, 2] == 1, 1.0, 0.0)
other = T[T.locus.isin(geneset | tr_only)]; opl = other.strand == '+'
for c_, s_, p_ in zip(other.chrom, other.strand, np.where(opl, other.start, other.end - 1)): setr(c_, s_, 0, p_ - TRANS, p_ + TRANS + 1, w=0.0)
for c_, s_, p_ in zip(other.chrom, other.strand, np.where(opl, other.end - 1, other.start)): setr(c_, s_, 1, p_ - TRANS, p_ + TRANS + 1, w=0.0)
for c_, s_, et, p_ in masked: setr(c_, s_, 0 if et == 'TSS' else 1, p_ - TRANS, p_ + TRANS + 1, w=0.0)
for r in POS.itertuples(): setr(r.chrom, r.strand, 0 if r.end_type == 'TSS' else 1, r.pos - TRANS, r.pos + TRANS + 1, w=0.0)
for r in POS.itertuples(): setr(r.chrom, r.strand, 0 if r.end_type == 'TSS' else 1, r.pos - WIN, r.pos + WIN + 1, y=1, w=r.weight, maxw=True)
# flank (long reads only)
gs = pd.read_csv('v3/gap_stats_v3.tsv', sep='\t'); gp = pd.read_csv('v3/gaps_v3.tsv', sep='\t')
DRNA = {'LeeVV2021_dRNA', 'Hoshizaki2022_dRNA'}; gs['is_drna'] = gs.dataset.isin(DRNA)
agg = gs.groupby('gap_id').agg(left=('left', 'sum'), right=('right', 'sum'), span=('span', 'sum'), n3=('n3', 'sum'), n3_tail=('n3_tail', 'sum'), n3_arich=('n3_arich', 'sum'), n5_G1=('n5_G1', 'sum'))
agg['n3_drna'] = gs[gs.is_drna].groupby('gap_id').n3.sum().reindex(agg.index).fillna(0)
A = gp.set_index('gap_id').join(agg); A['downA_tts'] = [downA(r.chrom, r.strand, int(r.tts_pos)) for r in A.itertuples()]
weak = A[['left', 'right']].min(axis=1); A['no_span'] = (A.span == 0) | ((A.span <= 2) & (A.span < 0.02 * weak))
A['tts_ok'] = (A.n3_drna >= 1) | (A.n3_tail >= 1) | (((A.n3 - A.n3_arich) >= 2) & (A.downA_tts < 0.6)); A['tss_ok'] = A.n5_G1 >= 1
A['enough'] = (A.left >= 5) & (A.right >= 5); n_end = A.tts_ok.astype(int) + A.tss_ok.astype(int)
A['tier'] = np.select([~A.no_span, ~A.enough, n_end == 2, n_end == 1], ['spanned', 'C', 'A', 'B'], 'C')
TW = {'A': 1.0, 'B': 0.6, 'C': 0.3}; nfl = 0
for g in A.reset_index().itertuples():
    if g.tier == 'spanned':
        for ch in range(3): setr(g.chrom, g.strand, ch, g.a, g.b, y=0, w=0.0)
        continue
    lo, hi = g.a + EDGE, g.b - EDGE
    if hi - lo < RES: continue
    a_, b_ = Bn(lo), Bn(hi - 1) + 1; k = si(g.strand); free = Y[g.chrom][k, 2, a_:b_] == 0
    for ch in range(3):
        W[g.chrom][k, ch, a_:b_] = np.where(free, TW[g.tier], W[g.chrom][k, ch, a_:b_]); Y[g.chrom][k, ch, a_:b_] = np.where(free, 0, Y[g.chrom][k, ch, a_:b_])
    nfl += int(free.sum())
np.savez_compressed('work/data/Lml3_pool.npz', **{f'{c}__Y': Y[c] for c in Y}, **{f'{c}__W': W[c].astype(np.float16) for c in W})
POS.to_csv('work/data/ml3_end_positives.tsv', sep='\t', index=False); A.reset_index().to_csv('v3/flank_gaps_v3_longread_tiers.tsv', sep='\t', index=False)
st.update(genes=len(geneset), truncated_only_genes=len(tr_only), flank_bins=nfl, positives=POS.groupby('end_type').size().to_dict(), tiers=A.tier.value_counts().to_dict())
print(st)
