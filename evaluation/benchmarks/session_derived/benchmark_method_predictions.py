"""Session-derived computation of benchmark_method_predictions.tsv.

Recovered from the record of the interactive analysis session in which benchmark_method_predictions.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: meta.json, oof_top8_genes.tsv.gz, oof_top8_cage.tsv, oof_nochrom_cage.tsv, ST_short_merged.gtf.gz, ST_long_merged.gtf.gz, X_chrom_ring.npz, X_chrom_troph.npz, X_chrom_schizont.npz
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import os, sys, glob, json, numpy as np, pandas as pd, pyBigWig
from scipy import stats

# Load gene boundary predictions
A_ = pd.read_csv(result('oof_top8_cage.tsv'), sep='\t')
B_ = pd.read_csv(result('oof_nochrom_cage.tsv'), sep='\t')

G = A_.merge(B_[['gene_id', 'pred_TSS', 'pred_TTS']].rename(columns={'pred_TSS': 'nc_TSS', 'pred_TTS': 'nc_TTS'}), on='gene_id')
conf = G.lr_group == 'confirmed'
G['final_TSS'] = np.where(conf, G.lr_TSS, G.pred_TSS)
G['final_TTS'] = np.where(conf, G.lr_TTS, G.pred_TTS)
G['final_src'] = np.where(conf, 'long-read', np.where(G.pred_TSS.notna(), 'model', 'none'))
G['anchor_TSS'] = np.where(conf, G.lr_TSS, G.nc_TSS)
G['anchor_TTS'] = np.where(conf, G.lr_TTS, G.nc_TTS)
G['utr5_final'] = (G.atg - G.final_TSS).abs()
G['utr3_final'] = (G.final_TTS - G.stop).abs()
G['utr5_anchor'] = (G.atg - G.anchor_TSS).abs()

# Load M8 for window info
M8 = pd.read_csv(result('oof_top8_genes.tsv.gz'), sep='\t')[['gene_id', 'tss_win_lo', 'tss_win_hi', 'tts_win_lo', 'tts_win_hi', 'max_TSS_prob_window', 'max_TTS_prob_window']]
B = G.merge(M8, on='gene_id')
Bm = B.set_index('gene_id')

# Load chromatin data
TRK = json.load(open(result('meta.json')))['chrom_tracks']
XC = {'ring': {c: v.astype(np.float32) for c, v in np.load(result('X_chrom_ring.npz')).items()},
      'troph': {c: v.astype(np.float32) for c, v in np.load(result('X_chrom_troph.npz')).items()},
      'schizont': {c: v.astype(np.float32) for c, v in np.load(result('X_chrom_schizont.npz')).items()}}
for s in XC:
    allv = np.concatenate([XC[s][c] for c in XC[s]], 1); mu, sd = allv.mean(1, keepdims=True), allv.std(1, keepdims=True) + 1e-6
    for c in XC[s]: XC[s][c] = (XC[s][c] - mu) / sd

ii = [TRK.index(m) for m in ('H2AZ', 'H3K9ac', 'ATAC')]; im = TRK.index('MNase')
PSC = {c: np.convolve(np.mean([XC[s][c][ii].sum(0) - XC[s][c][im] for s in XC], 0), np.ones(5) / 5, 'same') for c in XC['troph']}

# Build GT table from B
GT = B[['gene_id', 'chrom', 'strand', 'atg', 'stop', 'tss_win_lo', 'tss_win_hi', 'tts_win_lo', 'tts_win_hi', 'lr_group', 'lr_TSS', 'lr_TTS']].copy()
fold = lambda c: 0 if int(c[6:8]) <= 5 else (1 if int(c[6:8]) <= 10 else 2)
GT['fold'] = GT.chrom.map(fold)

def argmax_win(r):
    lo, hi = int(r.tss_win_lo), int(r.tss_win_hi)
    if hi - lo < 100: return np.nan
    v = PSC[r.chrom][lo // 20:hi // 20 + 1]; return lo // 20 * 20 + int(np.argmax(v)) * 20 + 10
GT['peak_raw'] = [argmax_win(r) for r in GT.itertuples()]
sgn = np.where(GT.strand == '+', 1, -1)
GT['peak_off'] = (GT.lr_TSS - GT.peak_raw) * sgn
GT['utr5_lr'] = (GT.atg - GT.lr_TSS).abs(); GT['utr3_lr'] = (GT.lr_TTS - GT.stop).abs()

M_ = GT[['gene_id']].copy()
for k in (0, 1, 2):
    tr = GT[(GT.fold != k) & GT.lr_TSS.notna()]; te = GT.fold == k
    off = tr.peak_off.median(); u5, u3 = tr.utr5_lr.median(), tr.utr3_lr.median()
    s_ = np.where(GT.strand[te] == '+', 1, -1)
    M_.loc[te, 'ChromatinPeak__TSS'] = GT.peak_raw[te] + s_ * off
    M_.loc[te, 'FixedOffset__TSS'] = GT.atg[te] - s_ * u5; M_.loc[te, 'FixedOffset__TTS'] = GT.stop[te] + s_ * u3

M_['PlasmoDB_v71__TSS'] = M_.gene_id.map(Bm.annot_tss.where(Bm.utr5_annot > 0))
M_['PlasmoDB_v71__TTS'] = M_.gene_id.map(Bm.annot_tts.where(Bm.utr3_annot > 0))
M_['Ours_chromatin__TSS'] = M_.gene_id.map(Bm.pred_TSS)
M_['Ours_chromatin__TTS'] = M_.gene_id.map(Bm.pred_TTS)
M_['Ours_nochromatin__TSS'] = M_.gene_id.map(Bm.nc_TSS)
M_['Ours_nochromatin__TTS'] = M_.gene_id.map(Bm.nc_TTS)

# ChromHMM
MK12 = ['H2AZ', 'H2BZ', 'H3', 'H3K27ac', 'H3K4me1', 'H3K4me3', 'H3K9ac', 'H3K9me3', 'HP1', 'PolII_S2P', 'ATAC', 'MNase']
os.makedirs('bench/chromhmm/bin', exist_ok=True)
raw = {'ring': np.load(result('X_chrom_ring.npz')),
       'troph': np.load(result('X_chrom_troph.npz')),
       'schizont': np.load(result('X_chrom_schizont.npz'))}
def b200(a): n = a.shape[1] // 10; return a[:, :n * 10].reshape(a.shape[0], n, 10).mean(2)
Bn = {s: {c: b200(raw[s][c][[TRK.index(m) for m in MK12]].astype(np.float32)) for c in raw[s].files} for s in XC}
thr = {s: np.percentile(np.concatenate([Bn[s][c] for c in Bn[s]], 1), 90, axis=1) for s in XC}
for s in XC:
    for c in Bn[s]:
        bb = (Bn[s][c] > thr[s][:, None]).astype(int)
        with open(f'bench/chromhmm/bin/{s}_{c}_binary.txt', 'w') as fh:
            fh.write(f'{s}\t{c}\n' + '\t'.join(MK12) + '\n'); np.savetxt(fh, bb.T, fmt='%d', delimiter='\t')
with open('bench/chromhmm/chromsizes.txt', 'w') as fh:
    for c in Bn['troph']: fh.write(f"{c}\t{Bn['troph'][c].shape[1] * 200}\n")

os.makedirs('bench/chromhmm/out10', exist_ok=True)
JV = os.path.join(os.environ.get('CONDA_PREFIX', ''), 'lib/jvm/bin/java')
CH = os.path.join(os.environ.get('CONDA_PREFIX', ''), 'share/chromhmm-1.28-0/ChromHMM.jar')
os.system(f"{JV} -Djava.awt.headless=true -Xmx6g -jar {CH} LearnModel -p 4 -noautoopen -nobrowser -b 200 -l bench/chromhmm/chromsizes.txt bench/chromhmm/bin bench/chromhmm/out10 10 Pf3D7 > bench/chromhmm/learn.log 2>&1")

SEG = {}
for s in XC:
    d = pd.read_csv(f'bench/chromhmm/out10/{s}_10_segments.bed', sep='\t', header=None, names=['chrom', 'start', 'end', 'state'])
    d['st'] = d.state.str[1:].astype(int)
    for c, q in d.groupby('chrom'):
        a = np.zeros(int(q.end.max()) // 200, np.int8)
        for r in q.itertuples(): a[r.start // 200:r.end // 200] = r.st
        SEG.setdefault(c, []).append(a)
SEG = {c: np.stack([v[:min(map(len, Ls))] for v in Ls]) for c, Ls in SEG.items()}

def enrich(genes):
    cnt = np.zeros(11); bg = np.zeros(11)
    for c, a in SEG.items(): bg += np.bincount(a.ravel(), minlength=11)
    for r in genes.itertuples():
        b = int(r.lr_TSS) // 200
        if b < SEG[r.chrom].shape[1]: cnt += np.bincount(SEG[r.chrom][:, b], minlength=11)
    return (cnt / cnt.sum()) / (bg / bg.sum())

def hmm_anchor(r, P):
    lo, hi = int(r.tss_win_lo) // 200, int(r.tss_win_hi) // 200 + 1; a = SEG[r.chrom][:, lo:hi]
    if a.shape[1] == 0: return np.nan
    m = np.isin(a, P).any(0)
    if not m.any(): return np.nan
    idx = np.where(m)[0]; br = np.r_[0, np.where(np.diff(idx) > 1)[0] + 1, len(idx)]; segs = [(idx[i], idx[j - 1]) for i, j in zip(br[:-1], br[1:])]
    atgb = int(r.atg) // 200 - lo; sg = min(segs, key=lambda x: min(abs(x[0] - atgb), abs(x[1] - atgb)))
    return (lo + (sg[0] + sg[1]) / 2) * 200 + 100

for k in (0, 1, 2):
    tr = GT[(GT.fold != k) & GT.lr_TSS.notna()]; e = enrich(tr); P = [i for i in range(1, 11) if e[i] >= 2]
    anc_tr = np.array([hmm_anchor(r, P) for r in tr.itertuples()]); off = np.nanmedian((tr.lr_TSS.values - anc_tr) * np.where(tr.strand == '+', 1, -1))
    te = GT[GT.fold == k]; anc = np.array([hmm_anchor(r, P) for r in te.itertuples()])
    M_.loc[te.index, 'ChromHMM__TSS'] = anc + np.where(te.strand == '+', 1, -1) * off

# StringTie
exec(open('bench/stringtie_to_genes.py').read()) if os.path.exists('bench/stringtie_to_genes.py') else None

def read_gtf(f):
    g = pd.read_csv(f, sep='\t', comment='#', header=None, names=['chrom','src','type','start','end','score','strand','frame','attr'])
    t = g[g.type == 'transcript'].copy(); t['start'] -= 1
    t['cov'] = t.attr.str.extract(r'cov "([0-9.eE+-]+)"')[0].astype(float).fillna(1.0)
    t['tss'] = np.where(t.strand == '+', t.start, t.end - 1); t['tes'] = np.where(t.strand == '+', t.end - 1, t.start)
    return t[t.strand.isin(['+', '-'])]

def assign(t, G):
    out = []; idx = {k: v for k, v in t.groupby(['chrom', 'strand'])}
    for r in G.itertuples():
        q = idx.get((r.chrom, r.strand))
        if q is None: out.append((r.gene_id, np.nan, np.nan)); continue
        lo, hi = min(r.atg, r.stop), max(r.atg, r.stop) + 1
        ov = (np.minimum(q.end.values, hi) - np.maximum(q.start.values, lo)).clip(0) / (hi - lo); q = q[ov >= 0.5]
        a = q[(q.tss >= r.tss_win_lo) & (q.tss <= r.tss_win_hi)]; b = q[(q.tes >= r.tts_win_lo) & (q.tes <= r.tts_win_hi)]
        out.append((r.gene_id, a.loc[a['cov'].idxmax(), 'tss'] if len(a) else np.nan, b.loc[b['cov'].idxmax(), 'tes'] if len(b) else np.nan))
    return pd.DataFrame(out, columns=['gene_id', 'TSS', 'TTS'])

GTw = GT.copy()
for nm, lab in (('ST_short_merged', 'StringTie_short'), ('ST_long_merged', 'StringTie_long')):
    f = result('ST_short_merged.gtf.gz') if nm == 'ST_short_merged' else result('ST_long_merged.gtf.gz')
    t = read_gtf(f); a = assign(t, GTw)
    M_[f'{lab}__TSS'] = M_.gene_id.map(a.set_index('gene_id').TSS)
    M_[f'{lab}__TTS'] = M_.gene_id.map(a.set_index('gene_id').TTS)

Mx = M_.copy()
Mx['LongRead_v3__TSS'] = Mx.gene_id.map(GT.set_index('gene_id').lr_TSS)
Mx['LongRead_v3__TTS'] = Mx.gene_id.map(GT.set_index('gene_id').lr_TTS)

# CAGE-based truths
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
Mx['Exp_Adjalley2016_capseq__TSS'] = Mx.gene_id.map(T_ad)

# Shaw 2021 5CAPture
import gzip
f_cap = 'bench/indep/GSE103036_combine_enrichment_score.txt.gz'
rows = []
with gzip.open(f_cap, 'rt') as fh:
    next(fh)
    for line in fh:
        x = line.rstrip('\n').split('\t')
        try: p = float(x[-1])
        except ValueError: continue
        if p >= 0.01: continue
        body = '\t'.join(x[3:-1]).split('\t*'); a = b = 0.0; ns = 0
        for g in body:
            t_ = [y for y in g.split('\t') if y]
            if t_ and t_[0] != 'n/a': a += float(t_[0]); b += float(t_[1]); ns += 1
        rows.append((x[2], x[1], int(x[0]), p, a, b, ns))
C5 = pd.DataFrame(rows, columns=['chrom', 'strand', 'pos', 'p', 'd1', 'd2', 'n_samples'])
C5 = C5[C5.chrom.str.match(r'Pf3D7_\d\d_v3')]

def cluster(df, gap=20, wcol='d1'):
    out = []
    for (c, s), q in df.sort_values('pos').groupby(['chrom', 'strand']):
        p = q.pos.values; w = q[wcol].values; br = np.r_[0, np.where(np.diff(p) > gap)[0] + 1, len(p)]
        for i, j in zip(br[:-1], br[1:]):
            k = i + int(np.argmax(w[i:j])); out.append((c, s, int(p[i]), int(p[j - 1]), int(p[k]), float(w[i:j].sum()), j - i))
    return pd.DataFrame(out, columns=['chrom', 'strand', 'start', 'end', 'dominant', 'score', 'n_pos'])
CAP = cluster(C5); CAP['dominant'] -= 1

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
Mx['Exp_Shaw2021_5CAPture__TSS'] = Mx.gene_id.map(T_sh)

# Chappell 2020 5UTR-seq
import openpyxl
f_chap = 'bench/indep/chappell/12864_2020_6787_MOESM10_ESM.xlsx'
wb = openpyxl.load_workbook(f_chap, read_only=True); parts = []
for sh in wb.sheetnames:
    ws = wb[sh]; rows_wb = list(ws.iter_rows(values_only=True)); h = [i for i, r in enumerate(rows_wb) if r and r[0] == 'Chromosome'][0]
    d = pd.DataFrame(rows_wb[h + 1:], columns=rows_wb[h]); d['tp'] = sh; parts.append(d)
S9 = pd.concat(parts).dropna(subset=['Chromosome']); S9 = S9[S9.Chromosome.astype(str).str.match(r'Pf3D7_\d\d_v3')]
for c_ in ('Quantile_10', 'Quantile_90', 'Dominant_tag'): S9[c_] = S9[c_].astype(int)
S9['Quantile_TPM'] = S9.Quantile_TPM.astype(float)
S9p = S9.assign(chrom=S9.Chromosome, strand=S9.Strand, pos=S9.Dominant_tag - 1, w=S9.Quantile_TPM)
T_ch = truth_from_points(S9p, 'pos', 'w', minw=5)
Mx['Exp_Chappell2020_5UTRseq__TSS'] = Mx.gene_id.map(T_ch)

# Consensus
TR_ = GT[['gene_id']].copy()
TR_['Adjalley2016_capseq__TSS'] = TR_.gene_id.map(T_ad)
TR_['Shaw2021_5CAPture__TSS'] = TR_.gene_id.map(T_sh)
TR_['Chappell2020_5UTRseq__TSS'] = TR_.gene_id.map(T_ch)
def consensus(r, tol=100):
    v = [r[f'{k}__TSS'] for k in ('Adjalley2016_capseq', 'Shaw2021_5CAPture', 'Chappell2020_5UTRseq') if pd.notna(r[f'{k}__TSS'])]
    best = None
    for a_ in v:
        s_ = [b_ for b_ in v if abs(b_ - a_) <= tol]
        if len(s_) >= 2 and (best is None or len(s_) > len(best)): best = s_
    return float(np.median(best)) if best else np.nan
TR_['Consensus_5end__TSS'] = TR_.apply(consensus, axis=1)
Mx['Exp_Consensus_5end__TSS'] = Mx.gene_id.map(TR_.set_index('gene_id').Consensus_5end__TSS)

Mx.to_csv('benchmark_method_predictions.tsv', sep='\t', index=False)
