"""Per-gene catalogue of transcript-end situations from TAMA collapsed models (7 long-read model sets) vs PlasmoDB v68.

Coordinates: cds5 / cds3 = 0-based first / last CDS base in transcription direction (from utr_stats.py).
u = distance upstream of the CDS start (5' side), v = distance downstream of the CDS end (3' side).

Connected model (5'): same strand, an exon block contains cds5, and the model extends >= 100 bp into the CDS.
Connected model (3'): same strand, an exon block contains cds3, and the model starts >= 100 bp inside the CDS.
5'/3' end clusters: pooled over datasets, single linkage (gap <= 100 bp in u / v).
A cluster is SUPPORTED if seen in >= 2 datasets or carried by >= 3 model reads; alternative TSS/TES require
>= 2 supported clusters >= 200 bp apart, each with >= 10% of the gene's connected reads.
Independent same-strand unit: model on the gene strand with NO exon overlapping any annotated CDS (any gene), supported
(>= 2 datasets or >= 3 reads), lying (>= 50% of its span) inside the annotated UTR or the window up to 2 kb upstream.
Antisense unit: opposite-strand model overlapping the annotated UTR, supported, not overlapping an opposite-strand annotated CDS.
"""
import numpy as np, pandas as pd
GAP, MIN_DS, MIN_READS, MIN_SHARE, ALT_SEP = 100, 2, 3, 0.10, 200
def exons(r):
    s = np.fromstring(r.block_starts.rstrip(','), sep=',', dtype=np.int64) + r.start
    return s, s + np.fromstring(r.block_sizes.rstrip(','), sep=',', dtype=np.int64)
def cluster(pos, w, ds):
    """pos, weights, dataset labels -> list of dict(lo, hi, peak, reads, n_ds)."""
    if len(pos) == 0: return []
    o = np.argsort(pos); pos, w, ds = pos[o], w[o], ds[o]
    cut = np.where(np.diff(pos) > GAP)[0] + 1; out = []
    for a, b in zip(np.r_[0, cut], np.r_[cut, len(pos)]):
        p, ww = pos[a:b], w[a:b]
        out.append(dict(lo=int(p.min()), hi=int(p.max()), peak=int(p[np.argmax(ww)]), reads=int(ww.sum()), n_ds=len(set(ds[a:b]))))
    return out
def supported(c): return c['n_ds'] >= MIN_DS or c['reads'] >= MIN_READS
class Index:
    def __init__(self, M):
        self.g = {k: (v.sort_values('start').reset_index(drop=True)) for k, v in M.groupby(['chrom', 'strand'])}
        self.st = {k: v.start.values for k, v in self.g.items()}
        self.mx = {k: np.maximum.accumulate(v.end.values - v.start.values).max() for k, v in self.g.items()}
    def query(self, c, s, a, b):
        k = (c, s)
        if k not in self.g: return self.g[next(iter(self.g))].iloc[0:0]
        st = self.st[k]; i0 = np.searchsorted(st, a - self.mx[k]); i1 = np.searchsorted(st, b)
        q = self.g[k].iloc[i0:i1]; return q[(q.end > a) & (q.start < b)]
def cds_mask_factory(cds):
    idx = {k: (np.sort(v.start.values), v.sort_values('start').end.values) for k, v in cds.groupby(['chrom', 'strand'])}
    mxl = {k: (e - s).max() for k, (s, e) in idx.items()}
    def overlaps_cds(c, s, starts, ends):
        t = idx.get((c, s))
        if t is None: return False
        cs, ce = t
        for a, b in zip(starts, ends):
            i0 = np.searchsorted(cs, a - mxl[(c, s)]); i1 = np.searchsorted(cs, b)
            if ((ce[i0:i1] > a) & (cs[i0:i1] < b)).any(): return True
        return False
    return overlaps_cds

DRNA = {'LeeVV2021_dRNA', 'Hoshizaki2022_dRNA'}
def read_fa(f):
    seqs, name, buf = {}, None, []
    for line in open(f):
        if line.startswith('>'):
            if name: seqs[name] = ''.join(buf)
            name = line[1:].split()[0]; buf = []
        else: buf.append(line.strip().upper())
    seqs[name] = ''.join(buf); return seqs
COMP = str.maketrans('ACGT', 'TGCA')
def downA(G, c, s, p3, n=20):
    w = G[c][p3 + 1:p3 + 1 + n] if s == '+' else G[c][max(0, p3 - n):p3].translate(COMP)[::-1]
    return w.count('A') / max(len(w), 1)
def unit_support(q, G):
    """q: candidate independent models -> (supported?, n_ds, reads) after removing A-rich-ended cDNA models (internal priming)."""
    if len(q) == 0: return False, 0, 0
    p3 = np.where(q.strand == '+', q.end - 1, q.start)
    ar = np.array([downA(G, c, s, p) for c, s, p in zip(q.chrom, q.strand, p3)])
    keep = q.dataset.isin(DRNA).values | (ar < 0.6)
    k = q[keep]; nd, rd = k.dataset.nunique(), int(k.n_reads.sum())
    return (nd >= MIN_DS or rd >= MIN_READS), nd, rd
def ends_side(Mi, g, side, chimeric=None, max_intron=5000):
    """connected-model ends on one side -> clusters, n datasets, reads, connected models."""
    c, s = g.chrom, g.strand; sg = 1 if s == '+' else -1
    anchor = g.cds5 if side == 5 else g.cds3
    q = Mi.query(c, s, anchor, anchor + 1)
    pos, w, ds, rows, nchim = [], [], [], [], 0
    for r in q.itertuples(index=False):
        a, b = exons(r)
        if not ((a <= anchor) & (b > anchor)).any(): continue
        if len(a) > 1 and (a[1:] - b[:-1]).max() > max_intron: continue
        if chimeric is not None and chimeric(g, a, b): nchim += r.n_reads; continue
        if side == 5:
            into = (r.end - 1 - anchor) if sg == 1 else (anchor - r.start)
            if into < 100: continue
            pos.append((anchor - r.start) if sg == 1 else (r.end - 1 - anchor))
        else:
            into = (anchor - r.start) if sg == 1 else (r.end - 1 - anchor)
            if into < 100: continue
            pos.append((r.end - 1 - anchor) if sg == 1 else (anchor - r.start))
        w.append(r.n_reads); ds.append(r.dataset); rows.append(r)
    cl = cluster(np.array(pos), np.array(w), np.array(ds, dtype=object))
    return cl, len(set(ds)), int(np.sum(w)) if w else 0, nchim, extent(np.array(pos), np.array(w), np.array(ds, dtype=object))
def extent(pos, w, ds):
    """farthest distance p such that models ending at >= p come from >= 2 datasets or carry >= 3 reads."""
    if len(pos) == 0: return np.nan
    o = np.argsort(-pos); pos, w, ds = pos[o], w[o], ds[o]
    seen, cw = set(), 0
    for p, ww, d in zip(pos, w, ds):
        seen.add(d); cw += ww
        if len(seen) >= MIN_DS or cw >= MIN_READS: return int(p)
    return np.nan
def alt_ends(cl, total):
    sup = [x for x in cl if supported(x) and x['reads'] >= MIN_SHARE * max(total, 1)]
    sup = sorted(sup, key=lambda x: x['peak']); keep = []
    for x in sup:
        if not keep or x['peak'] - keep[-1]['peak'] >= ALT_SEP: keep.append(x)
        elif x['reads'] > keep[-1]['reads']: keep[-1] = x
    return keep
def interval(g, side, u_lo, u_hi):
    """genomic [a, b) for distances u_lo..u_hi from the CDS anchor on the given side."""
    if side == 5:
        return (g.cds5 - u_hi, g.cds5 - u_lo + 1) if g.strand == '+' else (g.cds5 + u_lo, g.cds5 + u_hi + 1)
    return (g.cds3 + u_lo, g.cds3 + u_hi + 1) if g.strand == '+' else (g.cds3 - u_hi, g.cds3 - u_lo + 1)
