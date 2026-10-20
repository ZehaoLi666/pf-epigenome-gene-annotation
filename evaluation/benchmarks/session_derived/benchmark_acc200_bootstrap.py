"""Session-derived computation of benchmark_acc200_bootstrap.tsv.

Recovered from the record of the interactive analysis session in which benchmark_acc200_bootstrap.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: (none recorded)
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR  # paths: config/paths.yaml
import os, sys, glob, numpy as np, pandas as pd, pyBigWig, openpyxl, gzip, json
from scipy import stats

os.chdir('.')
sys.path[:0] = ['work', 'work/boundary_model']

# ── bench_eval helper ────────────────────────────────────────────────────────
TOL = (50, 100, 200, 500)

def evaluate(M, T, groups):
    meth = sorted({c.split('__')[0] for c in M if '__' in c})
    tru  = sorted({c.split('__')[0] for c in T if '__' in c})
    X = M.merge(T, on='gene_id', how='inner').merge(groups, on='gene_id', how='left')
    rows = []
    for end in ('TSS', 'TTS'):
        for t in tru:
            tc = f'{t}__{end}'
            if tc not in X: continue
            for g, q0 in [('all', X)] + list(X.groupby('group')):
                q = q0[q0[tc].notna()]
                if len(q) < 10: continue
                for m in meth:
                    mc = f'{m}__{end}'
                    if mc not in q: continue
                    d = (q[mc] - q[tc]).abs()
                    r = dict(end=end, truth=t, group=g, method=m, n_truth=len(q),
                             coverage=q[mc].notna().mean(), median_err=d.median())
                    for tol in TOL: r[f'acc@{tol}'] = (d <= tol).mean()
                    rows.append(r)
    return pd.DataFrame(rows)

# ── stringtie_to_genes helper ────────────────────────────────────────────────
def read_gtf(f):
    g = pd.read_csv(f, sep='\t', comment='#', header=None,
                    names=['chrom','src','type','start','end','score','strand','frame','attr'])
    t = g[g.type == 'transcript'].copy(); t['start'] -= 1
    t['cov'] = t.attr.str.extract(r'cov "([0-9.eE+-]+)"')[0].astype(float).fillna(1.0)
    t['tss'] = np.where(t.strand == '+', t.start, t.end - 1)
    t['tes'] = np.where(t.strand == '+', t.end - 1, t.start)
    return t[t.strand.isin(['+', '-'])]

def assign(t, G):
    out = []; idx = {k: v for k, v in t.groupby(['chrom', 'strand'])}
    for r in G.itertuples():
        q = idx.get((r.chrom, r.strand))
        if q is None: out.append((r.gene_id, np.nan, np.nan)); continue
        lo, hi = min(r.atg, r.stop), max(r.atg, r.stop) + 1
        ov = (np.minimum(q.end.values, hi) - np.maximum(q.start.values, lo)).clip(0) / (hi - lo)
        q = q[ov >= 0.5]
        a = q[(q.tss >= r.tss_win_lo) & (q.tss <= r.tss_win_hi)]
        b = q[(q.tes >= r.tts_win_lo) & (q.tes <= r.tts_win_hi)]
        out.append((r.gene_id,
                    a.loc[a['cov'].idxmax(), 'tss'] if len(a) else np.nan,
                    b.loc[b['cov'].idxmax(), 'tes'] if len(b) else np.nan))
    return pd.DataFrame(out, columns=['gene_id', 'TSS', 'TTS'])

# ── load gene boundaries and truth ───────────────────────────────────────────
B = pd.read_csv('gene_boundaries_v3_all_genes.tsv', sep='\t')
GT_full = pd.read_csv('work/gene_truth_v3.tsv', sep='\t')
B = B.merge(GT_full[['gene_id', 'lr_TSS', 'lr_TTS', 'lr_group']], on='gene_id', how='left')

# load gene windows from oof_top8
M8 = pd.read_csv('work/oof_top8_genes.tsv.gz', sep='\t')[
    ['gene_id', 'tss_win_lo', 'tss_win_hi', 'tts_win_lo', 'tts_win_hi',
     'atg', 'stop', 'chrom', 'strand', 'annot_tss', 'annot_tts',
     'utr5_annot', 'utr3_annot']]
B = B.merge(M8, on='gene_id', how='left', suffixes=('', '_m8'))
# fill missing atg/stop/chrom/strand from M8
for col in ('atg', 'stop', 'chrom', 'strand', 'tss_win_lo', 'tss_win_hi', 'tts_win_lo', 'tts_win_hi', 'annot_tss', 'annot_tts', 'utr5_annot', 'utr3_annot'):
    if col not in B.columns:
        B[col] = B[f'{col}_m8']
    else:
        B[col] = B[col].fillna(B.get(f'{col}_m8', np.nan))

GT = B[['gene_id', 'chrom', 'strand', 'atg', 'stop',
        'tss_win_lo', 'tss_win_hi', 'tts_win_lo', 'tts_win_hi',
        'lr_TSS', 'lr_TTS', 'lr_group']].copy()
GT['fold'] = GT.chrom.map(lambda c: 0 if int(c[6:8]) <= 5 else (1 if int(c[6:8]) <= 10 else 2))

# ── CAGE bigwig tracks ────────────────────────────────────────────────────────
CB = {(s, st): pyBigWig.open(
        f"{LAB_DIR}/lncRNA project tracks/cage_tss/CTSS_{s}_{st}_CPM.bw")
      for s in ('ring', 'trophozoite', 'schizont') for st in ('pos', 'neg')}

def track(c, strand, a, b):
    k = 'pos' if strand == '+' else 'neg'
    return np.max([np.nan_to_num(np.array(CB[(s, k)].values(c, max(0, a), b), float))
                   for s in ('ring', 'trophozoite', 'schizont')], axis=0)

# ── Shaw 2021 5CAPture ────────────────────────────────────────────────────────
f_shaw = 'bench/indep/GSE103036_combine_enrichment_score.txt.gz'
hdr = gzip.open(f_shaw, 'rt').readline().rstrip('\n').split('\t')
d1 = [h for h in hdr if h.endswith('(depth_1)')]
d2 = [h for h in hdr if h.endswith('(depth_2)')]
use = ['position', 'strand', 'chromosome', 'combine p-value'] + d1 + d2
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
            t = [y for y in g.split('\t') if y]
            if t and t[0] != 'n/a': a += float(t[0]); b += float(t[1]); ns += 1
        rows_shaw.append((x[2], x[1], int(x[0]), p, a, b, ns))
C5 = pd.DataFrame(rows_shaw, columns=['chrom', 'strand', 'pos', 'p', 'd1', 'd2', 'n_samples'])

def cluster(df, gap=20, wcol='d1'):
    out = []
    for (c, s), q in df.sort_values('pos').groupby(['chrom', 'strand']):
        p = q.pos.values; w = q[wcol].values
        br = np.r_[0, np.where(np.diff(p) > gap)[0] + 1, len(p)]
        for i, j in zip(br[:-1], br[1:]):
            k = i + int(np.argmax(w[i:j]))
            out.append((c, s, int(p[i]), int(p[j-1]), int(p[k]), float(w[i:j].sum()), j-i))
    return pd.DataFrame(out, columns=['chrom','strand','start','end','dominant','score','n_pos'])

CAP = cluster(C5); CAP['dominant'] -= 1
CAP = CAP[CAP.chrom.str.match(r'Pf3D7_\d\d_v3')]

# ── Chappell 2020 5UTR-seq tag clusters ──────────────────────────────────────
f_chap = glob.glob('bench/indep/chappell/*MOESM10_ESM.xlsx')[0]
wb = openpyxl.load_workbook(f_chap, read_only=True); parts_ch = []
for sh in wb.sheetnames:
    ws = wb[sh]; rows_c = list(ws.iter_rows(values_only=True))
    h = [i for i, r in enumerate(rows_c) if r and r[0] == 'Chromosome'][0]
    d = pd.DataFrame(rows_c[h+1:], columns=rows_c[h]); d['tp'] = sh; parts_ch.append(d)
S9 = pd.concat(parts_ch).dropna(subset=['Chromosome'])
S9 = S9[S9.Chromosome.astype(str).str.match(r'Pf3D7_\d\d_v3')]
for c_ in ('Quantile_10', 'Quantile_90', 'Dominant_tag'): S9[c_] = S9[c_].astype(int)
S9['Quantile_TPM'] = S9.Quantile_TPM.astype(float)

# ── truth helper: cluster points into gene TSS windows ───────────────────────
def truth_from_points(P, poscol, wcol, minw, gap=50):
    idx = {k: v.sort_values(poscol) for k, v in P.groupby(['chrom', 'strand'])}
    out = {}
    for r in GT.itertuples():
        q = idx.get((r.chrom, r.strand))
        if q is None: continue
        p = q[poscol].values
        i, j = np.searchsorted(p, r.tss_win_lo), np.searchsorted(p, r.tss_win_hi, side='right')
        if j <= i: continue
        pp, ww = p[i:j], q[wcol].values[i:j]
        br = np.r_[0, np.where(np.diff(pp) > gap)[0]+1, len(pp)]
        best = max(((ww[a:b].sum(), pp[a + int(np.argmax(ww[a:b]))]) for a, b in zip(br[:-1], br[1:])))
        if best[0] >= minw: out[r.gene_id] = best[1]
    return pd.Series(out)

S9p = S9.assign(chrom=S9.Chromosome, strand=S9.Strand, pos=S9.Dominant_tag - 1, w=S9.Quantile_TPM)
T_ch = truth_from_points(S9p, 'pos', 'w', minw=5)
T_sh = truth_from_points(CAP.rename(columns={'dominant': 'pos'}), 'pos', 'score', minw=100)

rng_adj = np.random.default_rng(0)
def adj(r, minw=5):
    lo, hi = int(r.tss_win_lo), int(r.tss_win_hi)
    if hi - lo < 50: return np.nan
    v = track(r.chrom, r.strand, lo-100, hi+100)
    sm = np.convolve(v, np.ones(201), 'valid')
    if sm.max() < minw: return np.nan
    c0 = int(np.argmax(sm)); seg = v[c0:c0+201]
    return lo - 100 + c0 + int(np.argmax(seg))

T_ad = pd.Series({r.gene_id: adj(r) for r in GT.itertuples()}).dropna()

TR_ = GT[['gene_id']].copy()
TR_['LongRead__TSS'] = GT.lr_TSS.values; TR_['LongRead__TTS'] = GT.lr_TTS.values
TR_['Adjalley2016_capseq__TSS'] = TR_.gene_id.map(T_ad)
TR_['Shaw2021_5CAPture__TSS']   = TR_.gene_id.map(T_sh)
TR_['Chappell2020_5UTRseq__TSS']= TR_.gene_id.map(T_ch)

def consensus(r, tol=100):
    v = [r[f'{k}__TSS'] for k in ('Adjalley2016_capseq', 'Shaw2021_5CAPture', 'Chappell2020_5UTRseq')
         if pd.notna(r[f'{k}__TSS'])]
    best = None
    for a_ in v:
        s = [b_ for b_ in v if abs(b_ - a_) <= tol]
        if len(s) >= 2 and (best is None or len(s) > len(best)): best = s
    return float(np.median(best)) if best else np.nan

TR_['Consensus_5end__TSS'] = TR_.apply(consensus, axis=1)

# ── chromatin feature matrix ──────────────────────────────────────────────────
TRK = json.load(open('work/data/meta.json'))['chrom_tracks']
XC = {s: {c: v.astype(np.float32) for c, v in np.load(f'work/data/X_chrom_{s}.npz').items()}
      for s in ('ring', 'troph', 'schizont')}
for s in XC:
    allv = np.concatenate([XC[s][c] for c in XC[s]], 1)
    mu, sd = allv.mean(1, keepdims=True), allv.std(1, keepdims=True) + 1e-6
    for c in XC[s]: XC[s][c] = (XC[s][c] - mu) / sd

ii = [TRK.index(m) for m in ('H2AZ', 'H3K9ac', 'ATAC')]; im = TRK.index('MNase')
PSC = {c: np.convolve(np.mean([XC[s][c][ii].sum(0) - XC[s][c][im] for s in XC], 0),
                      np.ones(5)/5, 'same')
       for c in XC['troph']}

def argmax_win(r):
    lo, hi = int(r.tss_win_lo), int(r.tss_win_hi)
    if hi - lo < 100: return np.nan
    v = PSC[r.chrom][lo//20:hi//20+1]
    return lo//20*20 + int(np.argmax(v))*20 + 10

GT['peak_raw'] = [argmax_win(r) for r in GT.itertuples()]
sgn = np.where(GT.strand == '+', 1, -1)
GT['peak_off'] = (GT.lr_TSS - GT.peak_raw) * sgn
GT['utr5_lr']  = (GT.atg - GT.lr_TSS).abs()
GT['utr3_lr']  = (GT.lr_TTS - GT.stop).abs()

M_ = GT[['gene_id']].copy()
for k in (0, 1, 2):
    tr = GT[(GT.fold != k) & GT.lr_TSS.notna()]
    te = GT.fold == k
    off = tr.peak_off.median()
    u5, u3 = tr.utr5_lr.median(), tr.utr3_lr.median()
    s_ = np.where(GT.strand[te] == '+', 1, -1)
    M_.loc[te, 'ChromatinPeak__TSS'] = GT.peak_raw[te] + s_ * off
    M_.loc[te, 'FixedOffset__TSS']   = GT.atg[te]  - s_ * u5
    M_.loc[te, 'FixedOffset__TTS']   = GT.stop[te] + s_ * u3

Bm = B.set_index('gene_id')
M_['PlasmoDB_v71__TSS']      = M_.gene_id.map(Bm['annot_tss'].where(Bm.get('utr5_annot', pd.Series(dtype=float)) > 0))
M_['PlasmoDB_v71__TTS']      = M_.gene_id.map(Bm['annot_tts'].where(Bm.get('utr3_annot', pd.Series(dtype=float)) > 0))
M_['Ours_chromatin__TSS']    = M_.gene_id.map(Bm['model_TSS'])
M_['Ours_chromatin__TTS']    = M_.gene_id.map(Bm['model_TTS'])
M_['Ours_nochromatin__TSS']  = M_.gene_id.map(Bm['chromfree_model_TSS'])
M_['Ours_nochromatin__TTS']  = M_.gene_id.map(Bm['chromfree_model_TTS'])

# ── ChromHMM states ───────────────────────────────────────────────────────────
SEG = {}
for s in XC:
    d = pd.read_csv(f'bench/chromhmm/out10/{s}_10_segments.bed', sep='\t', header=None,
                    names=['chrom','start','end','state'])
    d['st'] = d.state.str[1:].astype(int)
    for c, q in d.groupby('chrom'):
        a = np.zeros(int(q.end.max())//200, np.int8)
        for r in q.itertuples(): a[r.start//200:r.end//200] = r.st
        SEG.setdefault(c, []).append(a)
SEG = {c: np.stack([v[:min(map(len, Ls))] for v in Ls]) for c, Ls in SEG.items()}

def enrich(genes):
    cnt = np.zeros(11); bg = np.zeros(11)
    for c, a in SEG.items(): bg += np.bincount(a.ravel(), minlength=11)
    for r in genes.itertuples():
        b = int(r.lr_TSS)//200
        if b < SEG[r.chrom].shape[1]: cnt += np.bincount(SEG[r.chrom][:, b], minlength=11)
    return (cnt/cnt.sum()) / (bg/bg.sum())

def hmm_anchor(r, P):
    lo, hi = int(r.tss_win_lo)//200, int(r.tss_win_hi)//200+1
    a = SEG[r.chrom][:, lo:hi]
    if a.shape[1] == 0: return np.nan
    m = np.isin(a, P).any(0)
    if not m.any(): return np.nan
    idx = np.where(m)[0]; br = np.r_[0, np.where(np.diff(idx)>1)[0]+1, len(idx)]
    segs = [(idx[i], idx[j-1]) for i, j in zip(br[:-1], br[1:])]
    atgb = int(r.atg)//200 - lo
    sg = min(segs, key=lambda x: min(abs(x[0]-atgb), abs(x[1]-atgb)))
    return (lo + (sg[0]+sg[1])/2)*200 + 100

for k in (0, 1, 2):
    tr = GT[(GT.fold != k) & GT.lr_TSS.notna()]
    e = enrich(tr); P = [i for i in range(1, 11) if e[i] >= 2]
    anc_tr = np.array([hmm_anchor(r, P) for r in tr.itertuples()])
    off = np.nanmedian((tr.lr_TSS.values - anc_tr) * np.where(tr.strand == '+', 1, -1))
    te = GT.fold == k; anc = np.array([hmm_anchor(r, P) for r in GT[te].itertuples()])
    M_.loc[te, 'ChromHMM__TSS'] = anc + np.where(GT.strand[te] == '+', 1, -1) * off

# ── StringTie baselines ───────────────────────────────────────────────────────
GTw = GT.copy()
for nm, lab in (('ST_short_merged', 'StringTie_short'), ('ST_long_merged', 'StringTie_long')):
    t = read_gtf(f'bench/{nm}.gtf.gz')
    a = assign(t, GTw)
    M_[f'{lab}__TSS'] = M_.gene_id.map(a.set_index('gene_id').TSS)
    M_[f'{lab}__TTS'] = M_.gene_id.map(a.set_index('gene_id').TTS)

Mx = M_.copy()
Mx['LongRead_v3__TSS'] = Mx.gene_id.map(GT.set_index('gene_id').lr_TSS)
Mx['LongRead_v3__TTS'] = Mx.gene_id.map(GT.set_index('gene_id').lr_TTS)
for k in ('Adjalley2016_capseq', 'Shaw2021_5CAPture', 'Chappell2020_5UTRseq'):
    Mx[f'Exp_{k}__TSS'] = Mx.gene_id.map(TR_.set_index('gene_id')[f'{k}__TSS'])

GRP = GT[['gene_id', 'lr_group']].rename(columns={'lr_group': 'group'})
TRb = TR_.copy()

# ── bootstrap confidence intervals ───────────────────────────────────────────
METH = [('Ours_chromatin',   'Ours: CNN, chromatin + DNA + RNA',   '#1f5fa8'),
        ('Ours_nochromatin', 'Ours: CNN, no chromatin',             '#7fb2e5'),
        ('PlasmoDB_v71',     'PlasmoDB v71 UTRs',                   '#d08c2e'),
        ('FixedOffset',      'Fixed UTR length (median)',            '#bdbdbd'),
        ('StringTie_short',  'StringTie, stranded short reads',      '#8c6bb1'),
        ('StringTie_long',   'StringTie, long reads',                '#c2a5cf'),
        ('ChromatinPeak',    'Chromatin peak (H2A.Z+H3K9ac+ATAC)',  '#74c476'),
        ('ChromHMM',         'ChromHMM promoter states',             '#238b45')]
TRUTHS = [('Consensus_5end',        'Consensus of 3 cap/5′-end datasets'),
          ('Adjalley2016_capseq',   'Adjalley 2016 cap-seq'),
          ('Chappell2020_5UTRseq',  'Chappell 2020 5′UTR-seq'),
          ('Shaw2021_5CAPture',     'Shaw 2021 5CAPture'),
          ('LongRead',              'Long-read reference (training truth)')]

rng = np.random.default_rng(1)
X = Mx.merge(TRb, on='gene_id').merge(GRP, on='gene_id')

def boot(end, truth, m, tol=200, nb=500, sub=None):
    q = X[X[f'{truth}__{end}'].notna()] if sub is None else X[X[f'{truth}__{end}'].notna() & sub]
    h = ((q[f'{m}__{end}'] - q[f'{truth}__{end}']).abs() <= tol).values.astype(float)
    bs = [h[rng.integers(0, len(h), len(h))].mean() for _ in range(nb)]
    return h.mean(), np.percentile(bs, 2.5), np.percentile(bs, 97.5), len(h)

CI = pd.DataFrame([
    dict(end=e, truth=t, method=m,
         **dict(zip(['acc200','lo','hi','n'], boot(e, t, m))))
    for e, ts in (('TSS', [t for t, _ in TRUTHS]), ('TTS', ['LongRead']))
    for t in ts
    for m, _, _ in METH
    if f'{m}__{e}' in X
])

CI.to_csv('benchmark_acc200_bootstrap.tsv', sep='\t', index=False)
