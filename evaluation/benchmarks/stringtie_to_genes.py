"""Assign StringTie transcripts to genes: same strand, overlapping >=50% of the CDS envelope; the transcript with the highest coverage
whose 5' end lies in the gene's TSS search window gives the TSS (and, separately, the highest-coverage one whose 3' end lies in the TTS window gives the TTS)."""
import sys, re, numpy as np, pandas as pd
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
