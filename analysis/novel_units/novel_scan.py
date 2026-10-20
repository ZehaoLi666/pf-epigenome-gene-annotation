"""Genome-wide out-of-fold transcription units from the final model configuration (A_full: 3-seed ensemble per chromosome fold).
Fold A predicts chr 1-5 (val 6-7), B chr 6-10 (val 11), C chr 11-14 (val 10). TSS/TTS calls: logit-peak caller with threshold / NMS tuned on
each fold's validation chromosomes (one-to-one F1@200, as in the evaluation). Units: runs of smoothed body probability >= 0.5 (>= 300 bp);
each unit gets the strongest TSS call within 300 bp of its 5' edge and TTS call within 300 bp of its 3' edge.
Outputs: novel_scan_units.tsv, novel_scan_calls.tsv (all TSS/TTS calls genome-wide)."""
import sys, numpy as np, pandas as pd
sys.path[:0] = ['.', 'boundary_model']
import peak_eval as PE
RES = 20; num = lambda c: int(c[6:8])
FOLDS = {'A': ('', [6, 7], range(1, 6), 'logit'), 'B': ('_cvB', [11], range(6, 11), 'logitall'), 'C': ('_cvC', [10], range(11, 15), 'logitall')}
units, calls = [], []
for f, (sfx, va, te, lk) in FOLDS.items():
    nm = '+'.join(f'A_full{s}{sfx}' for s in ('', '_s1', '_drop'))
    PR = PE.load(nm, 'pred'); LG = PE.load(nm, lk)
    vac = [c for c in PR if num(c) in va]; tec = sorted(c for c in PR if num(c) in te)
    CL = {}
    for ch in (0, 1):
        mk = lambda chroms, t, d: {(c, k): PE.call(PR[c][k, ch], LG[c][k, ch][:PR[c].shape[-1]], 'logit_nms', t, d) for c in chroms for k in (0, 1)}
        _, thr, d = max((PE.score(mk(vac, t, d), vac, ch, 'one_to_one', 200)[2], t, d) for t in np.round(np.arange(0.1, 0.96, 0.05), 2) for d in (60, 100, 160, 300))
        print(f, ch, thr, d, flush=True)
        for (c, k), idx in mk(tec, thr, d).items():
            CL[(c, k, ch)] = idx
            for i in idx: calls.append((c, '+' if k == 0 else '-', 'TSS' if ch == 0 else 'TTS', int(i) * RES + RES // 2, float(PR[c][k, ch, i]), f))
    for c in tec:
        for k in (0, 1):
            st = '+' if k == 0 else '-'; b = np.convolve(PR[c][k, 2], np.ones(5) / 5, 'same') >= 0.5
            e = np.diff(np.r_[0, b.astype(int), 0]); s_, t_ = np.where(e == 1)[0], np.where(e == -1)[0]
            for a, z in zip(s_, t_):
                if (z - a) * RES < 300: continue
                lo, hi = a * RES, z * RES; five, three = (lo, hi) if st == '+' else (hi, lo)
                d = dict(chrom=c, strand=st, start=lo, end=hi, len_bp=hi - lo, body_mean=float(PR[c][k, 2, a:z].mean()), fold=f)
                for ch, et, edge in ((0, 'TSS', five), (1, 'TTS', three)):
                    idx = CL[(c, k, ch)]; m = np.abs(idx * RES + RES // 2 - edge) <= 300
                    if m.any():
                        j = idx[m][np.argmax(PR[c][k, ch, idx[m]])]; d[f'{et}_pos'] = int(j) * RES + RES // 2; d[f'{et}_prob'] = float(PR[c][k, ch, j])
                    else: d[f'{et}_pos'] = np.nan; d[f'{et}_prob'] = float(PR[c][k, ch, max(0, edge // RES - 15):edge // RES + 16].max())
                units.append(d)
pd.DataFrame(units).to_csv('novel_scan_units.tsv', sep='\t', index=False)
pd.DataFrame(calls, columns=['chrom', 'strand', 'end_type', 'pos', 'prob', 'fold']).to_csv('novel_scan_calls.tsv', sep='\t', index=False)
print('DONE', len(units), len(calls))
