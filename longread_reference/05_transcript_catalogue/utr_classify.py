"""Per-UTR classification from utr_stats.py output (same rules as utr5_connectivity_summary.tsv).
classify(path) -> per-transcript DataFrame with the flags full_length_supported, partial_only, utr_only_reads, candidate_overlong."""
import numpy as np, pandas as pd
from scipy import stats
import utr_analysis as UA
DROP = {'GSE186109_PacBio_schiztroph_cyto'}
CDNA = {'Yang2021_IsoSeq', 'GSE186109_PacBio_ring', 'GSE186109_PacBio_troph', 'GSE186109_PacBio_schizont', 'GSE186109_ONT_cDNA'}
def fisher(p):
    p = np.clip(np.asarray(p, float), 1e-300, 1); return stats.chi2.sf(-2 * np.log(p).sum(), 2 * len(p)) if len(p) else np.nan
def extent(hs):
    h = np.sum([UA.hist(x) for x in hs], axis=0); tail = h[::-1].cumsum()[::-1]
    ok = np.where(tail >= 2)[0]; return UA.EDGES[ok.max()] if len(ok) else np.nan
def reach_tol(hs, L, tol=100):
    h = UA.hist(hs); b = int(np.floor((max(L - tol, 0) + 300) / 100)); return h[min(b, UA.NB):].sum()
def classify(path):
    D = UA.load(path); D = D[D.transcript.str.endswith('.1') & ~D.dataset.isin(DROP)]
    d5 = D[D.has5 & (D.L5 > 0)].copy(); DS = list(dict.fromkeys(d5.dataset))
    hcds, SURV = {}, {}
    for ds in DS:
        x = d5[d5.dataset == ds]; h = np.sum([UA.hist(v) for v in x.hist_anchored_5end], axis=0)
        risk = h[::-1].cumsum()[::-1]; haz = h / np.maximum(risk, 1); hcds[ds] = (haz[0] + haz[1]) / 2
        SURV[ds], _ = UA.survival_from_reference(D, ds)
    T = d5
    T['reach_tol'] = [reach_tol(h, L) for h, L in zip(T.hist_anchored_5end, T.L5)]
    T['exp_tol_cds'] = T.reach_cds * (1 - T.dataset.map(hcds)) ** (np.maximum(T.L5 - 100, 0) / 100)
    T['exp_tol_km'] = [n * float(SURV[ds].values[min(int((max(L - 100, 0) + 300) // 100), UA.NB)]) for n, ds, L in zip(T.reach_cds, T.dataset, T.L5)]
    for k in ('cds', 'km'):
        T[f'p_tol_{k}'] = [stats.binom.cdf(a, n, min(e / n, 1)) if n > 0 else np.nan for a, n, e in zip(T.reach_tol, T.reach_cds, T[f'exp_tol_{k}'])]
    fr = np.where(T.utr_only > 0, T.utr_only_Arich / T.utr_only.replace(0, np.nan), 0)
    T['utr_only_clean'] = np.where(T.dataset.isin(CDNA), T.utr_only_in_utr * (1 - np.nan_to_num(fr)), T.utr_only_in_utr)
    T['ds_utr_only'] = (T.utr_only_clean >= 2).astype(int); T['ds_full'] = (T.reach_tol >= 2).astype(int)
    sub = T[T.reach_cds >= 5]
    A = T.groupby('transcript').agg(gene=('gene', 'first'), L5=('L5', 'first'), reach_cds=('reach_cds', 'sum'), reach_tol=('reach_tol', 'sum'),
            n_ds_cov=('reach_cds', lambda x: int((x >= 5).sum())), n_ds_full=('ds_full', 'sum'), n_ds_utr_only=('ds_utr_only', 'sum'),
            extent_bp=('hist_anchored_5end', extent))
    A = A.join(sub.groupby('transcript').agg(p_cdsnull=('p_tol_cds', fisher), p_kmnull=('p_tol_km', fisher)))
    cov = A.n_ds_cov >= 2
    A['full_length_supported'] = cov & (A.n_ds_full >= 2)
    A['partial_only'] = cov & (A.extent_bp < A.L5 - np.maximum(200, 0.25 * A.L5))
    A['utr_only_reads'] = cov & (A.n_ds_utr_only >= 2)
    A['candidate_overlong'] = A.partial_only | A.utr_only_reads
    A['group'] = np.select([~cov, A.full_length_supported & ~A.candidate_overlong, A.candidate_overlong & A.full_length_supported, A.candidate_overlong],
                           ['Insufficient coverage', 'Full length supported', 'Candidate, also full-length reads', 'Candidate, no full-length reads'], 'Covered, ambiguous')
    return A, pd.Series(hcds), d5.groupby('dataset').reach_cds.sum()
