"""Session-derived computation of tss_decomposition_controls.tsv.

Recovered from the record of the interactive analysis session in which tss_decomposition_controls.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: tss_stage_usage.tsv, report_uorf_translation.md, cage_rep_window_counts.tsv, daft_gene_counts.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import RESULTS_DIR, result  # paths: config/paths.yaml
import pandas as pd
import numpy as np
import re
from scipy import stats, optimize, special
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

# Load gene boundaries
B = RESULTS_DIR + '/'
gb = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep='\t').set_index('gene_id')

# Load stage alt TSS grade AB
AB = pd.read_csv(result('stage_alt_tss_gradeAB.tsv'), sep='\t')

# Load all candidates
CA = pd.read_csv(result('stage_alt_tss_all_candidates.tsv'), sep='\t', low_memory=False)

# Load CAGE window counts
cw = pd.read_csv(result('cage_rep_window_counts.tsv'), sep='\t')

# Load DAFT gene counts
dg = pd.read_csv(result('daft_gene_counts.tsv'), sep='\t', comment='#')

bcols = [c_ for c_ in dg.columns if c_.endswith('.bam')]
meta = pd.DataFrame({
    'col': bcols,
    'run': [re.search(r'(ERR\d+)_t(\d+)', c_).group(1) for c_ in bcols],
    'h': [int(re.search(r'_t(\d+)\.bam', c_).group(1)) for c_ in bcols]
})

DG = dg.set_index('Geneid')[bcols]
DG.columns = [f'{r}_{h}h' for r, h in zip(meta.run, meta.h)]
meta['lib'] = DG.columns

libsize = dict(zip(meta.lib, DG.sum(0)))
CPM = DG / DG.sum(0) * 1e6
cpm_h = {h: CPM[[l for l, hh in zip(meta.lib, meta.h) if hh == h]].mean(1) for h in [0, 8, 16, 24, 32, 40, 48]}

STG = ['ring', 'troph', 'schizont']
h2stage = {0: 'ring', 8: 'ring', 16: 'troph', 24: 'troph', 32: 'schizont', 40: 'schizont'}
cagehpi2daft = {2: 0, 10: 8, 18: 16, 26: 24, 34: 32, 42: 40}

# Build S (gradeAB pairs + controls)
ABp = CA[CA.pair_id.isin(set(AB.alt_pos.astype(str).radd(AB.gene_id + '|')))].copy()
CT = CA[CA.is_control].copy()
S = pd.concat([ABp, CT], ignore_index=True)

sg = np.where(S.strand == '+', 1, -1)
S['distal_pos'] = np.where(S.d_alt_vs_main < 0, S.alt_pos, S.main_TSS).astype(int)
S['proximal_pos'] = np.where(S.d_alt_vs_main < 0, S.main_TSS, S.alt_pos).astype(int)
S['distal_is_alt'] = S.d_alt_vs_main < 0

# Compute windows (simplified - use w column from CA which already has it)
# The windows are stored in S via CA which already has 'w'

# Gene window bounds - needed for prop_models normalization
stops = {k: d for k, d in gb.reset_index().groupby(['chrom', 'strand'])}

def genewin(r):
    s = 1 if r.strand == '+' else -1
    far = max(abs(r.distal_pos - r.atg) + 300, 1000)
    a = r.atg - s * far
    b = r.atg + s * 100
    d = stops[(r.chrom, r.strand)]
    d = d[d.gene_id != r.gene_id]
    if s > 0:
        up = d[(d.stop < r.atg) & (d.stop > a)]
        if len(up): a = max(a, int(up.stop.max()) + 1)
        return int(a), int(b)
    else:
        up = d[(d.stop > r.atg) & (d.stop < a)]
        if len(up): a = min(a, int(up.stop.min()) - 1)
        return int(b), int(a)

gw = S.apply(genewin, axis=1)
S['gwin_lo'] = [x[0] for x in gw]
S['gwin_hi'] = [x[1] for x in gw]

# Helper functions
def nd_hess(f, x, h=1e-4):
    n = len(x)
    H = np.zeros((n, n))
    for i in range(n):
        for j in range(i, n):
            ei = np.zeros(n); ej = np.zeros(n); ei[i] = h; ej[j] = h
            H[i, j] = H[j, i] = (f(x + ei + ej) - f(x + ei - ej) - f(x - ei + ej) + f(x - ei - ej)) / (4 * h * h)
    return H

def bb_nll(params, X, k, n):
    b = params[:-1]; rho = special.expit(params[-1])
    mu = special.expit(X @ b); mu = np.clip(mu, 1e-9, 1 - 1e-9)
    a = mu * (1 - rho) / rho; bb = (1 - mu) * (1 - rho) / rho
    return -np.sum(special.gammaln(n + 1) - special.gammaln(k + 1) - special.gammaln(n - k + 1) +
                   special.betaln(k + a, n - k + bb) - special.betaln(a, bb))

def bb_fit(X, k, n):
    p0 = np.zeros(X.shape[1] + 1)
    p0[0] = special.logit(np.clip(k.sum() / max(n.sum(), 1), 0.01, 0.99)); p0[-1] = -3
    best = None
    for start in [-3, -1, -5]:
        p0[-1] = start
        r = optimize.minimize(bb_nll, p0, args=(X, k, n), method='BFGS')
        if best is None or r.fun < best.fun: best = r
    return best

def prop_models(k, n, stage, cov):
    k = np.asarray(k, float); n = np.asarray(n, float); ok = n > 0
    k, n, stage, cov = k[ok], n[ok], stage[ok], cov[ok]
    pool = {s: (k[stage == s].sum(), n[stage == s].sum()) for s in STG}
    ps = {s: (a / b if b >= 10 else np.nan) for s, (a, b) in pool.items()}
    valid = [s for s in STG if np.isfinite(ps[s])]
    res = dict(n_total=n.sum(), **{f'prop_{s}': ps[s] for s in STG}, **{f'n_{s}': pool[s][1] for s in STG})
    if len(valid) < 2 or n.sum() < 20 or len(set(stage)) < 2:
        return res
    hi = max(valid, key=lambda s: ps[s]); lo = min(valid, key=lambda s: ps[s])
    res.update(dprop=ps[hi] - ps[lo], prop_hi_stage=hi, prop_lo_stage=lo)
    keep = np.isin(stage, valid); k, n, stage, cov = k[keep], n[keep], stage[keep], cov[keep]
    lv = valid; D = np.column_stack([np.ones(len(k))] + [(stage == l).astype(float) for l in lv[1:]])
    Xf = np.column_stack([D, cov]); Xc = np.column_stack([np.ones(len(k)), cov])
    rf = bb_fit(Xf, k, n); rc_ = bb_fit(Xc, k, n); rs = bb_fit(D, k, n); r0 = bb_fit(np.ones((len(k), 1)), k, n)
    dfree = len(lv) - 1
    res['bb_p_stage_adjDAFT'] = stats.chi2.sf(max(2 * (rc_.fun - rf.fun), 0), dfree)
    res['bb_p_stage_unadj'] = stats.chi2.sf(max(2 * (r0.fun - rs.fun), 0), dfree)
    res['bb_rho'] = special.expit(rf.x[-1])
    bvec = np.r_[0, rf.x[1:len(lv)]]
    try:
        H = nd_hess(lambda q: bb_nll(q, Xf, k, n), rf.x); V = np.linalg.inv(H)
        Vs = np.zeros((len(lv), len(lv))); Vs[1:, 1:] = V[1:len(lv), 1:len(lv)]
    except Exception:
        Vs = np.full((len(lv), len(lv)), np.nan)
    i, j = lv.index(hi), lv.index(lo); d = bvec[i] - bvec[j]
    se = np.sqrt(max(Vs[i, i] + Vs[j, j] - 2 * Vs[i, j], 0)) if np.all(np.isfinite(Vs)) else np.nan
    res.update(OR_hi_vs_lo_adj=np.exp(d), OR_lo95=np.exp(d - 1.96 * se), OR_hi95=np.exp(d + 1.96 * se))
    try:
        y = np.column_stack([k, n - k])
        g1 = sm.GLM(y, Xf, family=sm.families.Binomial()).fit(scale='X2')
        g0 = sm.GLM(y, Xc, family=sm.families.Binomial()).fit(scale='X2')
        Fst = (g0.deviance - g1.deviance) / dfree / g1.scale
        res['qb_p_stage_adjDAFT'] = stats.f.sf(Fst, dfree, max(g1.df_resid, 1))
        res['qb_scale'] = g1.scale
    except Exception:
        pass
    return res

def nb_nll(par, X, y, off):
    b = par[:-1]; a = np.exp(par[-1]); mu = np.exp(np.clip(X @ b + off, -50, 50)); r = 1 / a
    return -np.sum(special.gammaln(y + r) - special.gammaln(r) - special.gammaln(y + 1) +
                   r * np.log(r / (r + mu)) + y * np.log(mu / (r + mu) + 1e-300))

def nb_fit2(X, y, off):
    b0 = np.zeros(X.shape[1]); b0[0] = np.log((y.sum() + 0.5) / np.exp(off).sum())
    best = None
    for la in [-3, -1, 1]:
        r = optimize.minimize(nb_nll, np.r_[b0, la], args=(X, y, off), method='BFGS')
        if best is None or r.fun < best.fun: best = r
    return best

def nb_stage_full(y, stage, offset):
    y = np.asarray(y, float)
    if y.sum() < 10: return None
    D = np.column_stack([np.ones(len(y))] + [(stage == l).astype(float) for l in STG[1:]])
    m1 = nb_fit2(D, y, offset); m0 = nb_fit2(D[:, :1], y, offset)
    p = stats.chi2.sf(max(2 * (m0.fun - m1.fun), 0), 2)
    H = nd_hess(lambda q: nb_nll(q, D, y, offset), m1.x, h=1e-3)
    try: V_ = np.linalg.inv(H)
    except Exception: V_ = np.full((4, 4), np.nan)
    b = np.r_[0, m1.x[1:3]]; V = np.zeros((3, 3)); V[1:, 1:] = V_[1:3, 1:3]
    return dict(p=p, b=b, V=V)

def contrast(fit, hi, lo):
    i, j = STG.index(hi), STG.index(lo); d = fit['b'][i] - fit['b'][j]
    V = fit['V']; se = np.sqrt(max(V[i, i] + V[j, j] - 2 * V[i, j], 0))
    return d, se

# Process CAGE window counts
cw['st'] = cw.stage.replace({'trophozoite': 'troph'})
cw['dafth'] = cw.hpi.map(cagehpi2daft)

# Build P (pair table)
P = S[['pair_id', 'gene_id', 'chrom', 'strand', 'atg', 'main_TSS', 'alt_pos',
        'distal_pos', 'proximal_pos', 'distal_is_alt', 'd_alt_vs_main', 'is_control']].copy()
P['set'] = np.where(P.is_control, 'control_singleTSS', 'gradeAB')
P = P.merge(AB[['gene_id', 'alt_pos', 'grade', 'desc']], on=['gene_id', 'alt_pos'], how='left')

# Compute CAGE proportion models
RC = {}
for pid, d in cw.groupby('pair_id'):
    gid = pid.split('|')[0].replace('CTRL_', '')
    cov = np.log2(np.array([cpm_h[h][gid] for h in d.dafth]) + 0.5)
    r = prop_models(d.distal.values, (d.distal + d.proximal).values, d.st.values, cov)
    r['cage_gene_total_by_stage'] = ';'.join(f"{s}:{d[d.st == s].gene_total.sum():.0f}" for s in STG)
    RC[pid] = r
RC_df = pd.DataFrame(RC).T.add_prefix('CAGE_')
RC_df.index.name = 'pair_id'

P = P.join(RC_df, on='pair_id')

# Compute NB fits for expression
dl = meta[meta.h.isin(h2stage)]
st_d = np.array([h2stage[h] for h in dl.h])
off_d = np.log(dl.lib.map(libsize).values.astype(float))

genes = sorted(set(S.gene_id))
NBD = {g_: nb_stage_full(DG.loc[g_, dl.lib].values, st_d, off_d) for g_ in genes}

# Load TSO window counts (u5)
u5 = pd.read_csv('s2/u5_tso_window_counts.tsv', sep='\t')

def run_u5(u5):
    u5 = u5.copy()
    u5['h'] = u5.title.str.extract(r'_(\d+)h_').astype(int)
    u5 = u5[u5.h.isin(h2stage)]
    u5['st'] = u5.h.map(h2stage)
    R = {}
    for pid, d in u5.groupby('pair_id'):
        gid = pid.split('|')[0].replace('CTRL_', '')
        cov = np.log2(np.array([cpm_h[h][gid] for h in d.h]) + 0.5)
        R[pid] = prop_models(d.distal.values, (d.distal + d.proximal).values, d.st.values, cov)
    return pd.DataFrame(R).T.add_prefix('U5_'), {}

RU, NB5 = run_u5(u5)
P2 = P.join(RU, on='pair_id')

# Classify function
def classify(T):
    T = T.copy()
    for D in ['CAGE', 'U5']:
        for c_ in ([f'{D}_bb_p_stage_adjDAFT', f'{D}_bb_p_stage_unadj', f'{D}_qb_p_stage_adjDAFT',
                    f'{D}_dprop', f'{D}_OR_hi_vs_lo_adj', f'{D}_OR_lo95', f'{D}_OR_hi95'] +
                   [f'{D}_prop_{s}' for s in STG]):
            T[c_] = pd.to_numeric(T[c_], errors='coerce')
        m = T[f'{D}_bb_p_stage_adjDAFT'].notna()
        T.loc[m, f'{D}_q'] = multipletests(T.loc[m, f'{D}_bb_p_stage_adjDAFT'], method='fdr_bh')[1]
        m2 = T[f'{D}_qb_p_stage_adjDAFT'].notna()
        T.loc[m2, f'{D}_qb_q'] = multipletests(T.loc[m2, f'{D}_qb_p_stage_adjDAFT'], method='fdr_bh')[1]
        T[f'{D}_changed'] = (T[f'{D}_q'] < 0.1) & (T[f'{D}_dprop'] >= 0.2)
        T[f'{D}_unchanged'] = T[f'{D}_q'].notna() & ((T[f'{D}_dprop'] < 0.15) | (T[f'{D}_q'] >= 0.1))

    def dirc(r):
        if not (r.CAGE_changed and r.U5_changed): return np.nan
        hC, lC = r.CAGE_prop_hi_stage, r.CAGE_prop_lo_stage
        h5, l5 = r.U5_prop_hi_stage, r.U5_prop_lo_stage
        a = r[f'U5_prop_{hC}'] - r[f'U5_prop_{lC}']
        b = r[f'CAGE_prop_{h5}'] - r[f'CAGE_prop_{l5}']
        return bool(np.nan_to_num(a) > 0 and np.nan_to_num(b) > 0)

    T['prop_direction_consistent'] = T.apply(dirc, axis=1)
    T['prop_status'] = np.select(
        [T.prop_direction_consistent == True,
         (T.CAGE_unchanged | T.CAGE_q.isna()) & (T.U5_unchanged | T.U5_q.isna()) & (T.CAGE_q.notna() | T.U5_q.notna())],
        ['changed', 'unchanged'], 'ambiguous')

    pv = []; l2 = []; lo_ = []; hi_ = []; cs = []
    for r in T.itertuples():
        f = NBD.get(r.gene_id)
        if f is None:
            pv.append(np.nan); l2.append(np.nan); lo_.append(np.nan); hi_.append(np.nan); cs.append('')
            continue
        if r.prop_status == 'changed':
            h, l = r.CAGE_prop_hi_stage, r.CAGE_prop_lo_stage; tag = 'switch_contrast'
        else:
            b = f['b']; h, l = STG[int(np.argmax(b))], STG[int(np.argmin(b))]; tag = 'max_min'
        d, se = contrast(f, h, l)
        pv.append(f['p']); l2.append(d / np.log(2))
        lo_.append((d - 1.96 * se) / np.log(2)); hi_.append((d + 1.96 * se) / np.log(2))
        cs.append(f'{h}_vs_{l} ({tag})')

    T['DAFT_p_stage'] = pv; T['DAFT_log2FC'] = l2
    T['DAFT_log2FC_lo95'] = lo_; T['DAFT_log2FC_hi95'] = hi_; T['DAFT_contrast'] = cs
    m = T.DAFT_p_stage.notna()
    T.loc[m, 'DAFT_q'] = multipletests(T.loc[m, 'DAFT_p_stage'], method='fdr_bh')[1]
    T['expr_changed'] = (T.DAFT_q < 0.1) & (T.DAFT_log2FC.abs() >= 1)
    T['expr_unchanged'] = T.DAFT_q.notna() & ((T.DAFT_log2FC.abs() < 1) | (T.DAFT_q >= 0.1))
    T['case'] = np.select(
        [(T.prop_status == 'unchanged') & T.expr_changed,
         (T.prop_status == 'changed') & T.expr_unchanged,
         (T.prop_status == 'changed') & T.expr_changed,
         (T.prop_status == 'unchanged') & T.expr_unchanged],
        ['1_expression_only', '2_TSS_switch_only', '3_both', '4_neither'], '5_unclassified')
    return T

CTT = classify(P2[P2.set == 'control_singleTSS'])

ctl_cols = ['pair_id', 'gene_id', 'desc', 'distal_pos', 'proximal_pos',
            'CAGE_n_total', 'CAGE_dprop', 'CAGE_q', 'U5_n_total', 'U5_dprop', 'U5_q',
            'prop_status', 'DAFT_log2FC', 'DAFT_q', 'case']
CTT[[c_ for c_ in ctl_cols if c_ in CTT.columns]].round(5).to_csv('tss_decomposition_controls.tsv', sep='\t', index=False)
