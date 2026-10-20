import sys, os; sys.path.append(os.getcwd())
import numpy as np, pandas as pd, pyBigWig, re, pickle
from scipy.stats import binomtest, wilcoxon
from statsmodels.stats.multitest import multipletests
from common import *

NP = 'ext/nuc/Pf_nucleosome_annotation_of_nucDetective/nucleosome_profiles/'
ND = 'ext/nuc/output_nucDetective/Inspector/03_NUCS/bed_allNucs/'
TPS = ['T5', 'T10', 'T15', 'T20', 'T25', 'T30', 'T35', 'T40']
L = chrom_sizes()
rng = np.random.default_rng(2026)

# ---------------- signals ----------------
raw = {tp: {} for tp in TPS}
for tp in TPS:
    bw = pyBigWig.open(NP + f'MNase-seq_{tp}_qNorm.bw')
    for c in CHROMS:
        raw[tp][c] = np.nan_to_num(np.array(bw.values(c, 0, L[c]), dtype=np.float32))
    bw.close()
def smooth(v, w=21):
    cs = np.concatenate([[0], np.cumsum(v, dtype=np.float64)]); h = w // 2
    lo = np.clip(np.arange(len(v)) - h, 0, len(v)); hi = np.clip(np.arange(len(v)) + h + 1, 0, len(v))
    return ((cs[hi] - cs[lo]) / (hi - lo)).astype(np.float32)
SIG = {'mean8': {c: smooth(np.mean([raw[tp][c] for tp in TPS], 0)) for c in CHROMS}}
for tp in TPS: SIG[tp] = {c: smooth(raw[tp][c]) for c in CHROMS}
ZERO = {c: smooth((np.mean([raw[tp][c] for tp in TPS], 0) == 0).astype(np.float32)) for c in CHROMS}
del raw
MED = {k: float(np.median(np.concatenate([v[c] for c in CHROMS]))) for k, v in SIG.items()}

def call_ndr(sig, thr, minw=50, gap=20):
    out = {}
    for c in CHROMS:
        b = np.concatenate([[False], sig[c] < thr, [False]]).astype(np.int8)
        d = np.diff(b); s = np.where(d == 1)[0]; e = np.where(d == -1)[0]
        if len(s):
            keep_s = [s[0]]; keep_e = [e[0]]
            for a, z in zip(s[1:], e[1:]):
                if a - keep_e[-1] <= gap: keep_e[-1] = z
                else: keep_s.append(a); keep_e.append(z)
            s, e = np.array(keep_s), np.array(keep_e); k = (e - s) >= minw; s, e = s[k], e[k]
        out[c] = (s, e)
    return out

def dyad_gap_ndr(tp, mingap=197):
    d = pd.read_csv(ND + f'MNase-seq_{tp}.bed', sep='\t', header=None)
    d['dy'] = d[3].str.split('_').str[1].astype(int); out = {}
    for c in CHROMS:
        dy = np.sort(d.dy.values[d[0].values == c]); g = np.diff(dy); k = g > mingap
        out[c] = (dy[:-1][k] + 73, dy[1:][k] - 73)
    return out

# ---------------- genes / anchors ----------------
G = load_genes(); bt = pd.read_csv(TRUTH, sep='\t')[['gene_id', 'Consensus_5end__TSS']]
base = G[G.ref_TSS.notna()].merge(bt, on='gene_id', how='left').reset_index(drop=True)
sg = np.where(base.strand.values == '+', 1, -1); atg = base.atg.values.astype(int)
ANCH = {'refined': base.ref_TSS.values, 'PlasmoDB_v71': base.v71_TSS.values, 'experimental': base.Consensus_5end__TSS.values,
        'ATG-341': (atg - 341 * sg).astype(float), 'ATG-506': (atg - 506 * sg).astype(float)}
dpos = (base.ref_TSS.values - base.v71_TSS.values) * sg   # >0: v71 upstream of refined (v71 longer)
SETS = {'calibration_Consensus5end': base.Consensus_5end__TSS.notna().values,
        'discordant_trimmed(v71 >300bp upstream)': dpos > 300, 'discordant_extended(refined >300bp upstream)': dpos < -300,
        'all_refined': np.ones(len(base), bool)}
for g in ['confirmed', 'truncated-only', 'HC model, unconfirmed', 'no HC model']: SETS[f'all|{g}'] = base.lr_group.values == g
# A/T runs
seq = {}; nm = None; buf = []
for l in open(FASTA):
    if l[0] == '>':
        if nm: seq[nm] = ''.join(buf)
        nm = l[1:].split()[0]; buf = []
    else: buf.append(l.strip().upper())
seq[nm] = ''.join(buf)
RUN = {}
for c in CHROMS:
    m = np.zeros(len(seq[c]), bool)
    for mm in re.finditer(r'A{8,}|T{8,}', seq[c]): m[mm.start():mm.end()] = True
    RUN[c] = m
def near_run(X):
    return np.array([np.isfinite(x) and RUN[c][max(0, int(x) - 20):int(x) + 21].any() for c, x in zip(base.chrom.values, X)])

# ---------------- per-gene metrics ----------------
def metrics(ndr, X, win=(0, 150)):
    n = len(base); d = np.full(n, np.nan); inN = np.zeros(n, bool); upN = np.zeros(n, bool); has = np.zeros(n, bool)
    for i in range(n):
        x = X[i]
        if not np.isfinite(x): continue
        c = base.chrom.values[i]; s, e = ndr[c]; a = atg[i]
        if sg[i] == 1:
            lo, hi = a - 3000, a
            j0, j1 = np.searchsorted(e, lo), np.searchsorted(e, hi, side='right'); ss, ee = s[j0:j1], e[j0:j1]; edge = ee
        else:
            lo, hi = a, a + 3000
            j0, j1 = np.searchsorted(s, lo), np.searchsorted(s, hi, side='right'); ss, ee = s[j0:j1], e[j0:j1]; edge = ss - 1
        if len(ss) == 0: continue
        has[i] = True
        dd = (x - edge) * sg[i]; k = np.argmin(np.abs(dd)); d[i] = dd[k]
        inN[i] = np.any((ss <= x) & (x < ee))
        if sg[i] == 1: upN[i] = np.any((ss > x) & (ee <= a))
        else: upN[i] = np.any((ss > a) & (ee <= x))
    hit = np.isfinite(d) & (d >= win[0]) & (d <= win[1])
    return dict(d=d, hit=hit, inN=inN, upN=upN, has=has)

def boot_ci(v, nb=2000):
    v = np.asarray(v, float); idx = rng.integers(0, len(v), (nb, len(v))); b = v[idx].mean(1); return np.percentile(b, [2.5, 97.5])

def summarize(cfg, R, extra_ctrl=None):
    rows = []; prs = []
    for sname, m0 in SETS.items():
        for an in ANCH:
            m = m0 & np.isfinite(ANCH[an])
            if m.sum() < 20: continue
            r = R[an]; lo, hi = boot_ci(r['hit'][m])
            rows.append(dict(config=cfg, gene_set=sname, anchor=an, n=int(m.sum()), frac_has_NDR=r['has'][m].mean(), frac_edge_0_150=r['hit'][m].mean(), CI95_lo=lo, CI95_hi=hi,
                             frac_in_NDR=r['inN'][m].mean(), frac_UTR_contains_NDR=r['upN'][m].mean(), median_d=np.nanmedian(r['d'][m]), median_abs_d=np.nanmedian(np.abs(r['d'][m]))))
        if extra_ctrl is not None:
            m = m0 & np.isfinite(ANCH['refined']); vals = extra_ctrl[:, m].mean(1)
            rows.append(dict(config=cfg, gene_set=sname, anchor='random_shift(refined±U500,20x)', n=int(m.sum()), frac_edge_0_150=vals.mean(), CI95_lo=np.percentile(vals, 2.5), CI95_hi=np.percentile(vals, 97.5)))
        if extra_ctrl is not None:
            for a1 in ['experimental', 'refined', 'PlasmoDB_v71']:
                m = m0 & np.isfinite(ANCH[a1])
                if m.sum() < 20: continue
                dv = R[a1]['hit'][m].astype(float) - extra_ctrl[:, m].mean(0); lo, hi = boot_ci(dv)
                prs.append(dict(config=cfg, gene_set=sname, comparison=f'{a1} vs random_shift', n=int(m.sum()), frac_1=R[a1]['hit'][m].mean(), frac_2=extra_ctrl[:, m].mean(), diff=dv.mean(), diff_CI95_lo=lo, diff_CI95_hi=hi, mcnemar_p=np.nan))
        for a1, a2 in [('refined', 'PlasmoDB_v71'), ('refined', 'experimental'), ('experimental', 'ATG-341'), ('experimental', 'ATG-506'), ('refined', 'ATG-341')]:
            m = m0 & np.isfinite(ANCH[a1]) & np.isfinite(ANCH[a2])
            if m.sum() < 20: continue
            h1, h2 = R[a1]['hit'][m], R[a2]['hit'][m]; n10 = int((h1 & ~h2).sum()); n01 = int((~h1 & h2).sum())
            lo, hi = boot_ci(h1.astype(float) - h2.astype(float))
            d1, d2 = np.abs(R[a1]['d'][m]), np.abs(R[a2]['d'][m]); ok = np.isfinite(d1) & np.isfinite(d2)
            wp = wilcoxon(d1[ok], d2[ok]).pvalue if ok.sum() > 10 and np.any(d1[ok] != d2[ok]) else np.nan
            prs.append(dict(config=cfg, gene_set=sname, comparison=f'{a1} vs {a2}', n=int(m.sum()), frac_1=h1.mean(), frac_2=h2.mean(), diff=h1.mean() - h2.mean(), diff_CI95_lo=lo, diff_CI95_hi=hi,
                            n_only1=n10, n_only2=n01, mcnemar_p=binomtest(n10, n10 + n01).pvalue if n10 + n01 else 1.0,
                            n_wilcoxon=int(ok.sum()), median_abs_d_1=np.median(d1[ok]), median_abs_d_2=np.median(d2[ok]), wilcoxon_absd_p=wp))
    return rows, prs

def random_ctrl(ndr, reps=20):
    out = []
    for _ in range(reps):
        X = base.ref_TSS.values + rng.uniform(-500, 500, len(base))
        X = np.where(sg == 1, np.clip(X, atg - 3000, atg - 1), np.clip(X, atg + 1, atg + 3000))
        out.append(metrics(ndr, X)['hit'])
    return np.array(out, float)

ALLR, ALLP = [], []; PG = None; ATSTRAT = []
configs = [('PRIMARY mean8 thr0.5 w50', 'mean8', 0.5, 50, (0, 150))]
configs += [(f'timepoint {tp} thr0.5 w50', tp, 0.5, 50, (0, 150)) for tp in TPS]
configs += [('mean8 thr0.4 w50', 'mean8', 0.4, 50, (0, 150)), ('mean8 thr0.6 w50', 'mean8', 0.6, 50, (0, 150)),
            ('mean8 thr0.5 w30', 'mean8', 0.5, 30, (0, 150)), ('mean8 thr0.5 w80', 'mean8', 0.5, 80, (0, 150)),
            ('mean8 window[0,100]', 'mean8', 0.5, 50, (0, 100)), ('mean8 window[-50,200]', 'mean8', 0.5, 50, (-50, 200))]
NDRS = {}
for cfg, sig, f, w, win in configs:
    ndr = call_ndr(SIG[sig], f * MED[sig], minw=w); NDRS[cfg] = ndr
    R = {an: metrics(ndr, X, win) for an, X in ANCH.items()}
    ctrl = random_ctrl(ndr) if cfg.startswith('PRIMARY') or cfg.startswith('timepoint T10') or cfg.startswith('timepoint T30') or cfg.startswith('timepoint T40') else None
    r, p = summarize(cfg, R, ctrl); ALLR += r; ALLP += p
    if cfg.startswith('PRIMARY'):
        PRIM = R
        PG = pd.DataFrame({'gene_id': base.gene_id, 'lr_group': base.lr_group, 'strand': base.strand, 'atg': atg})
        for an in ANCH:
            PG[f'{an}_anchor'] = ANCH[an]; PG[f'{an}_d_to_edge'] = R[an]['d']; PG[f'{an}_edge_0_150'] = R[an]['hit']; PG[f'{an}_in_NDR'] = R[an]['inN']; PG[f'{an}_UTR_contains_NDR'] = R[an]['upN']
            PG[f'{an}_near_AT_run'] = near_run(ANCH[an])
        nndr = sum(len(v[0]) for v in ndr.values()); wid = np.concatenate([v[1] - v[0] for v in ndr.values()])
        zf = []
        for c in CHROMS:
            s, e = ndr[c]; cz = np.concatenate([[0], np.cumsum(ZERO[c] > 0.5)]); zf.append((cz[e] - cz[s]) / (e - s))
        zf = np.concatenate(zf)
        NDRSTAT = dict(median_signal=MED['mean8'], threshold=0.5 * MED['mean8'], n_NDR_genome=nndr, median_width=float(np.median(wid)), frac_NDR_mostly_zero_coverage=float(np.mean(zf > 0.5)))
    print(cfg, flush=True)
# dyad-gap alternative
for tp in ['T10', 'T30', 'T40']:
    ndr = dyad_gap_ndr(tp); R = {an: metrics(ndr, X) for an, X in ANCH.items()}
    r, p = summarize(f'dyad-gap NDR {tp} (gap>197)', R, random_ctrl(ndr)); ALLR += r; ALLP += p; print('dyad', tp, flush=True)
# A/T strata on primary (each anchor stratified by its own context; and gene-level both/neither for refined vs v71)
nr = {an: near_run(X) for an, X in ANCH.items()}
for sname in ['calibration_Consensus5end', 'all_refined', 'discordant_trimmed(v71 >300bp upstream)', 'discordant_extended(refined >300bp upstream)']:
    for an in ANCH:
        for lab, flag in [('anchor near A/T run', True), ('anchor not near A/T run', False)]:
            m = SETS[sname] & np.isfinite(ANCH[an]) & (nr[an] == flag)
            if m.sum() < 20: continue
            lo, hi = boot_ci(PRIM[an]['hit'][m])
            ALLR.append(dict(config=f'PRIMARY | {lab}', gene_set=sname, anchor=an, n=int(m.sum()), frac_has_NDR=PRIM[an]['has'][m].mean(), frac_edge_0_150=PRIM[an]['hit'][m].mean(), CI95_lo=lo, CI95_hi=hi,
                             frac_in_NDR=PRIM[an]['inN'][m].mean(), frac_UTR_contains_NDR=PRIM[an]['upN'][m].mean(), median_d=np.nanmedian(PRIM[an]['d'][m]), median_abs_d=np.nanmedian(np.abs(PRIM[an]['d'][m]))))
    for lab, m in [('neither refined nor v71 near A/T run', ~nr['refined'] & ~nr['PlasmoDB_v71']), ('both near A/T run', nr['refined'] & nr['PlasmoDB_v71'])]:
        m = SETS[sname] & m
        h1, h2 = PRIM['refined']['hit'][m], PRIM['PlasmoDB_v71']['hit'][m]; n10 = int((h1 & ~h2).sum()); n01 = int((~h1 & h2).sum()); lo, hi = boot_ci(h1.astype(float) - h2.astype(float))
        ALLP.append(dict(config=f'PRIMARY | {lab}', gene_set=sname, comparison='refined vs PlasmoDB_v71', n=int(m.sum()), frac_1=h1.mean(), frac_2=h2.mean(), diff=h1.mean() - h2.mean(), diff_CI95_lo=lo, diff_CI95_hi=hi,
                         n_only1=n10, n_only2=n01, mcnemar_p=binomtest(n10, n10 + n01).pvalue if n10 + n01 else 1.0))
RES = pd.DataFrame(ALLR); PAIR = pd.DataFrame(ALLP); ok_ = PAIR.mcnemar_p.notna(); PAIR.loc[ok_, 'BH_q'] = multipletests(PAIR.mcnemar_p[ok_], method='fdr_bh')[1]
RES.to_csv('ndr_edge_fractions.tsv', sep='\t', index=False, float_format='%.4g')
PAIR.to_csv('ndr_edge_paired_tests.tsv', sep='\t', index=False, float_format='%.4g')
PG.to_csv('ndr_edge_per_gene_primary.tsv.gz', sep='\t', index=False, float_format='%.0f')
# metagene of time-averaged smoothed occupancy for figure
OFF = np.arange(-1000, 1001); MG = {'offset_bp': OFF}
for sname in ['calibration_Consensus5end', 'discordant_trimmed(v71 >300bp upstream)', 'discordant_extended(refined >300bp upstream)']:
    for an in ANCH:
        m = SETS[sname] & np.isfinite(ANCH[an]); acc = np.zeros(len(OFF)); k = 0
        for c, s, x in zip(base.chrom.values[m], sg[m], ANCH[an][m]):
            x = int(x); v = SIG['mean8'][c]
            if x - 1000 < 0 or x + 1001 > len(v): continue
            w_ = v[x - 1000:x + 1001]; acc += (w_ if s == 1 else w_[::-1]); k += 1
        MG[f'{sname}|{an}|n={k}'] = acc / k / MED['mean8']
pd.DataFrame(MG).to_csv('ndr_edge_metagene.tsv.gz', sep='\t', index=False, float_format='%.4g')
pickle.dump(dict(NDRSTAT=NDRSTAT, MED=MED), open('ndr_edge_stats.pkl', 'wb'))
print(NDRSTAT)
