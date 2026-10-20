"""Session-derived computation of cv_point_metrics.tsv.

Recovered from the record of the interactive analysis session in which cv_point_metrics.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: stage_TPM_9stages.tsv, gene_heterochromatin_groups.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
import sys, os, numpy as np, pandas as pd
sys.path[:0] = ['.', 'boundary_model']
import peak_eval as PE

RES = 20
FOLDS = {'A': ('', [6, 7], range(1, 6), 'logit'), 'B': ('_cvB', [11], range(6, 11), 'logitall'), 'C': ('_cvC', [10], range(11, 15), 'logitall')}
T = pd.read_csv('hpcc/tasks.tsv', sep='\t', header=None, names=['name', 'spec', 'tracks', 'fold', 'ns'])
P1 = list(dict.fromkeys(T[T.ns == 3].name))
CFG = [(n, ['', '_s1', '_drop']) for n in P1] + [('A_full2', None)] + [(n, ['', '_s1']) for n in dict.fromkeys(T[T.ns == 2].name)] + [(n, ['']) for n in T[T.ns == 1].name]
W = pd.read_csv('hpcc/gene_windows_truth.tsv', sep='\t')
LN = pd.read_csv('hpcc/genefree_lncRNA_tss.tsv', sep='\t')
num = lambda c: int(c[6:8])
ex = lambda f: os.path.exists(f'exp_ml/{f}.npz')
pts, genes, extra = [], [], []
for name, seeds in CFG:
    base = 'A_full' if name == 'A_full2' else name; seeds = ['', '_s1'] if name == 'A_full2' else seeds
    C, used = {}, []; par = {}
    for f, (sfx, va, te, lk) in FOLDS.items():
        nm = [f'{base}{s}{sfx}' for s in seeds]
        if not all(ex(f'pred_{x}') and ex(f'{lk}_{x}') for x in nm): continue
        PR = PE.load('+'.join(nm), 'pred'); LG = PE.load('+'.join(nm), lk); used.append(f)
        vac = [c for c in PR if num(c) in va]; tec = [c for c in PR if num(c) in te]
        for ch in (0, 1):
            mk = lambda chroms, t, d: {(c, k): PE.call(PR[c][k, ch], LG[c][k, ch][:PR[c].shape[-1]], 'logit_nms', t, d) for c in chroms for k in (0, 1)}
            _, thr, d = max((PE.score(mk(vac, t, d), vac, ch, 'one_to_one', 200)[2], t, d) for t in np.round(np.arange(0.1, 0.96, 0.05), 2) for d in (60, 100, 160, 300))
            par[(f, ch)] = (thr, d)
            for (c, k), idx in mk(tec, thr, d).items(): C[(c, k, ch)] = (idx, PR[c][k, ch, idx])
    if not used: print('missing', name, flush=True); continue
    tag = ''.join(used); chroms = sorted({c for c, _, _ in C})
    for ch, et in ((0, 'TSS'), (1, 'TTS')):
        calls = {(c, k): C[(c, k, ch)][0] for c in chroms for k in (0, 1)}
        for scheme in ('one_to_one', 'merged_region'):
            row = dict(model=name, folds=tag, end_type=et, scheme=scheme, n_chrom=len(chroms))
            for tol in (50, 100, 200, 500):
                P, R, F, n = PE.score(calls, chroms, ch, scheme, tol); row.update({f'P@{tol}': P, f'R@{tol}': R, f'F1@{tol}': F}); row['n_calls'] = n
            pts.append(row)
    for r in W[W.chrom.isin(chroms)].itertuples():
        k = 0 if r.strand == '+' else 1; d = dict(model=name, folds=tag, gene_id=r.gene_id)
        for ch, et, lo, hi in ((0, 'TSS', r.tss_win_lo, r.tss_win_hi), (1, 'TTS', r.tts_win_lo, r.tts_win_hi)):
            p, s = C[(r.chrom, k, ch)]; p = p * RES + RES // 2; m = (p >= lo) & (p <= hi)
            d[et] = int(p[m][np.argmax(s[m])]) if m.any() else np.nan
        genes.append(d)
    def near(c, k, x, tol=200):
        p = C[(c, k, 0)][0] * RES + RES // 2; return bool(len(p) and np.abs(p - x).min() <= tol)
    L_ = LN[LN.chrom.isin(chroms)]; q = W[W.chrom.isin(chroms) & W.lr_TSS.notna()]
    extra.append(dict(model=name, folds=tag, lncRNA_TSS_called=np.mean([near(r.chrom, 0 if r.strand == '+' else 1, r.tss) for r in L_.itertuples()]), n_lncRNA=len(L_),
                      mirror_TSS_opposite=np.mean([near(r.chrom, 1 if r.strand == '+' else 0, r.lr_TSS) for r in q.itertuples()]), n_lr=len(q),
                      **{f'thr_{f}_{"TSS" if ch == 0 else "TTS"}': v[0] for (f, ch), v in par.items()}))
    print('evaluated', name, tag, flush=True)
    pd.DataFrame(pts).to_csv('cv_point_metrics.tsv', sep='\t', index=False)
    pd.DataFrame(genes).to_csv('cv_gene_calls.tsv.gz', sep='\t', index=False)
    pd.DataFrame(extra).to_csv('cv_extra.tsv', sep='\t', index=False)
print('EVAL_DONE')
