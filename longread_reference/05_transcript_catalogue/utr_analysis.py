"""Descriptive statistics + truncation null model for annotated 5' UTRs (v68), per long-read dataset.
Null model (per dataset): S(d) = P(read 5' end reaches >= d bp upstream of the CDS start | read anchored >= 100 bp in the CDS),
estimated by Kaplan-Meier over ALL annotated 5' UTRs: a read whose 5' end falls inside the annotated UTR is a 'stop' at that
distance; reads reaching the annotated TSS are censored there (the molecule ended). Genes whose UTR is really shorter than
annotated add extra stops, so the pooled hazard is, if anything, too high -> expected counts too low -> the test is conservative.
Test per UTR: k = anchored reads reaching f*L5; expected = N_anchor_reachcds * S(f*L5)/S(0+); binomial lower tail p.
"""
import numpy as np, pandas as pd
from scipy import stats
NB = 63; EDGES = np.arange(-300, 6001, 100)       # bin i covers [EDGES[i], EDGES[i]+100); bin 63 = >= 6000
def hist(s):
    h = np.zeros(NB + 1)
    if isinstance(s, str) and s:
        for kv in s.split(','):
            k, n = kv.split(':'); h[int(k)] += int(n)
    return h
def load(path):
    D = pd.read_csv(path, sep='\t')
    for c in ('has5', 'has3'): D[c] = D[c].astype(str) == 'True'
    return D
def survival_from_reference(D, ds):
    """S(d) for d = bin lower edges >= 0: fraction of anchored reads with u5 >= d, using genes with full-UTR support,
    conditioned on reads reaching the CDS start (u5 > 0). Uses each gene's reads only up to its own TSS (censoring at L5)."""
    d = D[(D.dataset == ds) & D.has5 & (D.L5 > 0)]
    at_risk = np.zeros(NB + 1); events = np.zeros(NB + 1)
    for h, L in zip(d.hist_anchored_5end.map(hist), d.L5):
        lb = int((L + 300) // 100)                  # bin containing the TSS (censoring point)
        h = h.copy(); h[:3] = 0                       # bins < 0 (5' end inside CDS) excluded: condition u5 > 0
        n = h[3:].sum()
        if n == 0: continue
        # reads ending (5' end) in bin b "stop" there; genes contribute at-risk only up to their TSS bin
        stop = h.copy(); stop[lb:] = 0                 # stops at/after TSS are the true molecule end, not truncation
        cum = n - np.concatenate([[0], np.cumsum(h)[:-1]])
        cum[:3] = n
        at_risk[3:lb] += cum[3:lb]; events[3:lb] += stop[3:lb]
    haz = np.divide(events, at_risk, out=np.zeros_like(events), where=at_risk > 0)
    S = np.cumprod(1 - haz); S[:3] = 1.0
    return pd.Series(S, index=EDGES[:NB + 1]), at_risk
def test_utrs(D, ds, S):
    d = D[(D.dataset == ds) & D.has5 & (D.L5 > 0)].copy()
    Sv = S.values
    def s_at(x): return float(Sv[min(int((x + 300) // 100), NB)])
    out = {}
    for f in (50, 75, 100):
        e = np.array([n * s_at(f / 100 * L) for n, L in zip(d.reach_cds, d.L5)])
        k = d[f'reach_f{f}'].values
        d[f'exp_f{f}'] = e
        d[f'p_f{f}'] = [stats.binom.cdf(kk, n, min(ee / n, 1)) if n > 0 else np.nan for kk, n, ee in zip(k, d.reach_cds, e)]
    return d
