"""Compare end-calling and matching schemes on the same predictions (multi-label models, labels v2 truth, test chr 1-5).
Callers  : prob_nms300 (current: local maxima of P, >300 bp apart)
           plateau     (one call per contiguous run of P >= thr, at the maximum of the smoothed logit)
           logit_nmsD  (local maxima of the 3-bin-smoothed logit, P >= thr, NMS D bp; D tuned on validation)
Matching : one_to_one  (each call <-> one truth cluster; current)
           one_to_many (each truth -> its nearest call within tol; a call is correct if >=1 truth is assigned; redundant calls are FP)
           many        (a truth is recalled if any call is within tol; a call is correct if any truth is within tol)
           merged      (truth clusters merged when <=300 bp apart; one-to-one to the merged site, distance to the site span)
Unmatched calls in masked bins of the channel are neutral. Threshold (and D) tuned on chr 6-7 for F1 at 200 bp with the same scheme."""
import os, sys, numpy as np, pandas as pd
sys.path[:0] = ['.', 'boundary_model']
from bm.evaluate import match
RES = 20; num = lambda c: int(c[6:8])
Z = np.load(os.environ.get('PE_LABELS', 'data/Lml2_pool.npz')); CH = sorted({k.split('__')[0] for k in Z.files}); LW = {c: Z[f'{c}__W'].astype(np.float32) for c in CH}
POS = pd.concat([pd.read_csv(f, sep='\t')[['chrom', 'strand', 'pos', 'end_type']] for f in os.environ.get('PE_POS', 'data/ml_end_positives.tsv,data/ml2_added_end_positives.tsv').split(',')]).drop_duplicates()
TR = {k: np.sort(g.pos.values) for k, g in POS.groupby(['chrom', 'strand', 'end_type'])}
def merged_sites(p, gap=300):
    if len(p) == 0: return np.zeros((0, 2))
    out = [[p[0], p[0]]]
    for x in p[1:]:
        if x - out[-1][1] <= gap: out[-1][1] = x
        else: out.append([x, x])
    return np.array(out, float)
def load(names, kind):
    ps = [{k.split('__')[-1]: v.astype(np.float32) for k, v in np.load(f'exp_ml/{kind}_{n}.npz').items()} for n in names.split('+')]
    keys = set.intersection(*[set(p) for p in ps]); return {k: sum(p[k] for p in ps) / len(ps) for k in keys}
def smooth(v, w=3): return np.convolve(v, np.ones(w) / w, 'same')
def nms(idx, sc, d):
    keep = []; order = np.argsort(-sc); taken = np.zeros(len(idx), bool)
    for i in order:
        if taken[i]: continue
        keep.append(i); taken |= np.abs(idx - idx[i]) <= d
    keep = np.sort(np.array(keep, int)); return idx[keep], sc[keep]
def call(prob, logit, caller, thr, d=300):
    if caller == 'prob_nms300':
        from bm.decode import call_peaks          # the caller used in all earlier evaluations
        return call_peaks(prob, thr, 15)[0]
    if caller == 'plateau':
        on = prob >= thr; e = np.diff(np.r_[0, on.astype(int), 0]); st, en = np.where(e == 1)[0], np.where(e == -1)[0]; s = smooth(logit)
        return np.array([a + int(np.argmax(s[a:b])) for a, b in zip(st, en)], int)
    s = smooth(logit); m = (s[1:-1] >= s[:-2]) & (s[1:-1] > s[2:]) & (prob[1:-1] >= thr); idx = np.where(m)[0] + 1
    return nms(idx, s[idx], max(1, d // RES))[0]
def score(calls, chroms, ch, scheme, tol):
    et = 'TSS' if ch == 0 else 'TTS'; tp_c = fp = tp_t = nt = 0
    for c in chroms:
        for k, st in ((0, '+'), (1, '-')):
            cc = calls.get((c, k), np.zeros(0, int)) * RES + RES // 2; t = TR.get((c, st, et), np.zeros(0))
            if scheme == 'merged':
                S = merged_sites(t); nt += len(S)
                if len(cc) and len(S):
                    d = np.maximum(0, np.maximum(S[None, :, 0] - cc[:, None], cc[:, None] - S[None, :, 1]))
                    ii, jj = np.where(d <= tol); o = np.argsort(d[ii, jj], kind='stable'); uc, us = set(), set()
                    for q in o:
                        a, b = ii[q], jj[q]
                        if a in uc or b in us: continue
                        uc.add(a); us.add(b)
                    matched_c = np.array(sorted(uc), int); tp_t += len(us)
                else: matched_c = np.zeros(0, int)
            elif scheme == 'merged_region':
                S = merged_sites(t); nt += len(S)
                if len(cc) and len(S):
                    d = np.maximum(0, np.maximum(S[None, :, 0] - cc[:, None], cc[:, None] - S[None, :, 1]))
                    tp_t += int((d.min(0) <= tol).sum()); matched_c = np.where(d.min(1) <= tol)[0]
                else: matched_c = np.zeros(0, int)
            elif scheme == 'one_to_many':
                nt += len(t)
                if len(cc) and len(t):
                    dd = np.abs(cc[:, None] - t[None]); near = dd.argmin(0); ok = dd.min(0) <= tol
                    tp_t += int(ok.sum()); matched_c = np.unique(near[ok])
                else: matched_c = np.zeros(0, int)
            elif scheme == 'many':
                nt += len(t)
                if len(cc) and len(t):
                    dd = np.abs(cc[:, None] - t[None]); tp_t += int((dd.min(0) <= tol).sum()); matched_c = np.where(dd.min(1) <= tol)[0]
                else: matched_c = np.zeros(0, int)
            else:
                nt += len(t); pi, ti, _ = match(cc, t, tol); tp_t += len(pi); matched_c = pi
            tp_c += len(matched_c); un = np.setdiff1d(np.arange(len(cc)), matched_c)
            if len(un): fp += int((LW[c][k, ch, np.minimum(cc[un] // RES, LW[c].shape[-1] - 1)] > 0).sum())
    P = tp_c / max(tp_c + fp, 1); R = tp_t / max(nt, 1); return P, R, 2 * P * R / max(P + R, 1e-9), tp_c + fp
def run(label, names, callers=('prob_nms300', 'plateau', 'logit_nms'), schemes=('one_to_one', 'one_to_many', 'merged', 'many')):
    PR = load(names, 'pred'); LG = load(names, 'logit'); rows = []
    va = [c for c in PR if num(c) in (6, 7)]; te = [c for c in PR if num(c) <= 5]
    for ch, et in ((0, 'TSS'), (1, 'TTS')):
        for caller in callers:
            Ds = (60, 100, 160, 300) if caller == 'logit_nms' else (300,)
            cache = {}
            def calls_for(chroms, thr, d):
                key = (tuple(chroms), thr, d)
                if key not in cache: cache[key] = {(c, k): call(PR[c][k, ch], LG[c][k, ch][:PR[c].shape[-1]], caller, thr, d) for c in chroms for k in (0, 1)}
                return cache[key]
            for scheme in schemes:
                best = max(((score(calls_for(va, t, d), va, ch, scheme, 200)[2], t, d) for t in np.round(np.arange(0.1, 0.96, 0.05), 2) for d in Ds))
                _, thr, d = best; Ct = calls_for(te, thr, d)
                row = dict(model=label, end_type=et, caller=caller if caller != 'logit_nms' else f'logit_nms{d}', scheme=scheme, thr=thr)
                for tol in (50, 100, 200, 500):
                    P, R, F, n = score(Ct, te, ch, scheme, tol); row.update({f'P@{tol}': P, f'R@{tol}': R, f'F1@{tol}': F}); row['n_calls'] = n
                rows.append(row); print(label, et, row['caller'], scheme, round(row['F1@200'], 3), flush=True)
    return rows
if __name__ == '__main__':
    out = []
    for spec in sys.argv[2:]:
        lab, nm = spec.split('::'); out += run(lab, nm)
    pd.DataFrame(out).to_csv(sys.argv[1], sep='\t', index=False)

def report_fixed(label, names, caller='logit_nms', Ds=(60, 100, 160, 300), schemes=('one_to_one', 'one_to_many', 'merged', 'merged_region', 'many')):
    """Calls fixed once per end type (threshold, NMS distance tuned on chr 6-7 for one-to-one F1@200); every matching rule scores the SAME calls."""
    PR = load(names, 'pred'); LG = load(names, 'logit'); rows = []
    va = [c for c in PR if num(c) in (6, 7)]; te = [c for c in PR if num(c) <= 5]
    for ch, et in ((0, 'TSS'), (1, 'TTS')):
        mk = lambda chroms, t, d: {(c, k): call(PR[c][k, ch], LG[c][k, ch][:PR[c].shape[-1]], caller, t, d) for c in chroms for k in (0, 1)}
        _, thr, d = max((score(mk(va, t, d), va, ch, 'one_to_one', 200)[2], t, d) for t in np.round(np.arange(0.1, 0.96, 0.05), 2) for d in (Ds if caller == 'logit_nms' else (300,)))
        Ct = mk(te, thr, d); errs = []
        for c in te:
            for k, st in ((0, '+'), (1, '-')):
                _, _, e = match(Ct[(c, k)] * RES + RES // 2, TR.get((c, st, et), np.zeros(0)), 500); errs += list(np.abs(e))
        for scheme in schemes:
            row = dict(model=label, end_type=et, caller=caller, nms_bp=d, thr=thr, scheme=scheme, median_err_bp=float(np.median(errs)))
            for tol in (50, 100, 200, 500):
                P, R, F, n = score(Ct, te, ch, scheme, tol); row.update({f'P@{tol}': P, f'R@{tol}': R, f'F1@{tol}': F}); row['n_calls'] = n
            rows.append(row)
    return rows
