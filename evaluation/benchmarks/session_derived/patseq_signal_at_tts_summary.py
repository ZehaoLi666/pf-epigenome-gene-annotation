"""Session-derived computation of patseq_signal_at_tts_summary.tsv.

Recovered from the record of the interactive analysis session in which patseq_signal_at_tts_summary.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: gene_boundaries_v3_all_genes.tsv, pat_anchored_sites.tsv.gz, PAT_Ring_pos.bw, PAT_Ring_neg.bw, PAT_Troph_pos.bw, PAT_Troph_neg.bw, PAT_Schizont_pos.bw, PAT_Schizont_neg.bw
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import numpy as np
import pandas as pd
import pyBigWig
import glob
from scipy.stats import wilcoxon, binomtest

# Load gene boundaries
G = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep='\t')

# Load genome FASTA for downstream A-richness calculation
Gfa = {}; name = None; buf = []
for line in open(glob.glob(f'{LAB_DIR}/PlasmoDB/PlasmoDB-68/*Genome.fasta')[0]):
    if line[0] == '>':
        if name: Gfa[name] = ''.join(buf)
        name = line[1:].split()[0]; buf = []
    else: buf.append(line.strip().upper())
Gfa[name] = ''.join(buf)
COMP = str.maketrans('ACGTN', 'TGCAN')

def downA(c, s, p, n=20):
    w = Gfa[c][p + 1:p + 1 + n] if s == '+' else Gfa[c][max(0, p - n):p].translate(COMP)[::-1]
    return w.count('A') / max(len(w), 1)

# Load PAT-seq bigWigs
BW = {
    ('Ring', 'pos'): pyBigWig.open(result('PAT_Ring_pos.bw')),
    ('Ring', 'neg'): pyBigWig.open(result('PAT_Ring_neg.bw')),
    ('Troph', 'pos'): pyBigWig.open(result('PAT_Troph_pos.bw')),
    ('Troph', 'neg'): pyBigWig.open(result('PAT_Troph_neg.bw')),
    ('Schizont', 'pos'): pyBigWig.open(result('PAT_Schizont_pos.bw')),
    ('Schizont', 'neg'): pyBigWig.open(result('PAT_Schizont_neg.bw')),
}
CL = BW[('Troph', 'pos')].chroms()

def track(c, strand, a, b):
    sd = 'pos' if strand == '+' else 'neg'; a = max(0, a); b = min(CL[c], b)
    return np.vstack([np.nan_to_num(np.array(BW[(s, sd)].values(c, a, b), float)) for s in ('Ring', 'Troph', 'Schizont')])

# Compute TTS windows
G['lo_cds'] = G[['atg', 'stop']].min(axis=1); G['hi_cds'] = G[['atg', 'stop']].max(axis=1)
win = {}
for (c, s), q in G.groupby(['chrom', 'strand']):
    q = q.sort_values('lo_cds'); lo, hi, ids, stp = q.lo_cds.values, q.hi_cds.values, q.gene_id.values, q.stop.values
    for i in range(len(q)):
        if s == '+':
            nxt = lo[i + 1:][lo[i + 1:] > stp[i]]; e = min(stp[i] + 2000, (nxt.min() - 50) if len(nxt) else 10**9); win[ids[i]] = (int(stp[i]), int(e))
        else:
            prv = hi[:i][hi[:i] < stp[i]]; b = max(stp[i] - 2000, (prv.max() + 50) if len(prv) else 0); win[ids[i]] = (int(b), int(stp[i]))
G['tts_lo'] = G.gene_id.map(lambda g: win[g][0]); G['tts_hi'] = G.gene_id.map(lambda g: win[g][1])

sg = np.where(G.strand == '+', 1, -1)
med3 = (G.lr_TTS - G.stop).abs().median(); G['fixed_TTS'] = G.stop + sg * med3
G['v71_TTS'] = G.annot_tts.where(G.utr3_annot > 0)

METH = ['model_TTS', 'chromfree_model_TTS', 'v71_TTS', 'fixed_TTS', 'lr_TTS']

# Load anchored PAT-seq cleavage sites
A = pd.read_csv(result('pat_anchored_sites.tsv.gz'), sep='\t', header=None, names=['srr', 'chrom', 'strand', 'pos', 'n'])
STG = {'SRR39981131': 'Ring', 'SRR39981130': 'Ring', 'SRR39981126': 'Ring', 'SRR39981125': 'Troph', 'SRR39981124': 'Troph', 'SRR39981123': 'Troph', 'SRR39981122': 'Schizont', 'SRR39981121': 'Schizont', 'SRR39981120': 'Schizont'}
A['stage'] = A.srr.map(STG)
A = A[A.chrom.isin(Gfa.keys())].copy()
A['downA'] = [downA(c_, s_, p_) for c_, s_, p_ in zip(A.chrom, A.strand, A.pos)]

Aall = A[A.chrom.str.match(r'Pf3D7_\d\d_v3')].groupby(['chrom', 'strand', 'pos']).agg(n=('n', 'sum'), n_lib=('srr', 'nunique'), downA=('downA', 'first')).reset_index()
IDXa = {k: (v.pos.values, v.n.values) for k, v in Aall.groupby(['chrom', 'strand'])}

# Per-gene PAT-seq signal computation
W = 50; rng = np.random.default_rng(0)
rows = []
for r in G.itertuples():
    lo, hi = r.tts_lo, r.tts_hi
    if hi - lo < 200: continue
    v = track(r.chrom, r.strand, lo - W, hi + W + 1); sm = np.vstack([np.convolve(x, np.ones(2 * W + 1), 'valid') for x in v])
    best = sm.max(0); tot = float(best.max()); d = dict(gene_id=r.gene_id, win_max=tot)
    if tot > 0:
        pk = lo + int(np.argmax(best)); d['PAT_peak'] = pk; d['peak_downA'] = downA(r.chrom, r.strand, pk)
    sc = lambda p: np.nan if pd.isna(p) or not (lo <= p <= hi) else float(best[int(p) - lo])
    pc = lambda p: np.nan if pd.isna(p) or not (lo <= p <= hi) or tot <= 0 else float((best < best[int(p) - lo]).mean() + 0.5 * (best == best[int(p) - lo]).mean())
    for m in METH: p = getattr(r, m); d[f'PAT_{m}'] = sc(p); d[f'pct_{m}'] = pc(p)
    rp = int(rng.integers(lo, hi + 1)); d['PAT_random'] = sc(rp); d['pct_random'] = pc(rp)
    if pd.notna(r.model_TTS) and lo <= r.model_TTS <= hi:
        alt = [p for p in (r.model_TTS - 300, r.model_TTS + 300) if lo <= p <= hi]
        if alt: p = alt[int(rng.integers(len(alt)))]; d['PAT_shift300'] = sc(p); d['pct_shift300'] = pc(p)
    rows.append(d)
P = pd.DataFrame(rows); P = G[['gene_id', 'chrom', 'strand', 'stop', 'lr_group', 'final_src'] + METH].merge(P, on='gene_id')

LAB = {'model_TTS': 'Ours (chromatin)', 'chromfree_model_TTS': 'Ours (no chromatin)', 'v71_TTS': 'PlasmoDB v71', 'fixed_TTS': 'Fixed 3′UTR (384 bp)', 'shift300': 'Ours shifted ±300 bp', 'random': 'Random in window', 'lr_TTS': 'Long-read TTS'}

def summarise(Q, tag):
    out = []
    for m in ['model_TTS', 'chromfree_model_TTS', 'v71_TTS', 'fixed_TTS', 'shift300', 'random', 'lr_TTS']:
        k = Q[f'pct_{m}'].notna(); q = Q[k]; pos = Q[m] if m in Q else None
        d = (q[m] - q.PAT_peak).abs() if m in METH else pd.Series(dtype=float)
        out.append(dict(subset=tag, method=LAB[m], n=int(k.sum()), median_pct=q[f'pct_{m}'].median(), frac_PAT_ge1=(q[f'PAT_{m}'] >= 1).mean(),
                        median_PAT_cpm=q[f'PAT_{m}'].median(), within100_of_peak=(d <= 100).mean() if len(d) else np.nan, within200_of_peak=(d <= 200).mean() if len(d) else np.nan,
                        median_dist_peak=d.median() if len(d) else np.nan))
    return out

lrA = [downA(r.chrom, r.strand, int(r.lr_TTS)) for r in G[G.lr_TTS.notna()].itertuples()]
inf = P.win_max >= 1
S = pd.DataFrame(summarise(P[inf], 'PAT ≥1 CPM in window') + summarise(P[inf & (P.peak_downA < 0.6)], 'PAT peak not A-rich') + summarise(P[P.win_max >= 5], 'PAT ≥5 CPM')
                 + summarise(P[inf & P.lr_TTS.isna()], 'No long-read TTS') + summarise(P[inf & P.lr_TTS.notna()], 'Long-read TTS available'))

S.to_csv('patseq_signal_at_tts_summary.tsv', sep='\t', index=False)
