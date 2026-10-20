from pfconfig import result  # paths: config/paths.yaml
import sys, os; sys.path.append(os.getcwd())
import numpy as np, pandas as pd, pyBigWig, pickle
from common import *
NP = 'ext/nuc/Pf_nucleosome_annotation_of_nucDetective/nucleosome_profiles/'
ND = 'ext/nuc/output_nucDetective/Inspector/03_NUCS/bed_allNucs/'
TPS = ['T5', 'T10', 'T15', 'T20', 'T25', 'T30', 'T35', 'T40']
L = chrom_sizes()
PAIRS = result('divergent_pairs_coexpression_classes.tsv')
D = pd.read_csv(PAIRS, sep='\t')
# The input TSS_*_final / TSS_*_refinedRule columns were saved rounded to 4 significant digits (upstream float_format='%.4g').
# Restore exact coordinates from the gene table: final_TSS reproduces TSS_dist_final for 994/994 pairs and rounds to the stored value.
_G = load_genes().set_index('gene_id')
for side, gcol in [('left', 'gene_left'), ('right', 'gene_right')]:
    ex = _G.final_TSS.reindex(D[gcol]).values
    assert np.all(np.isnan(ex) == D[f'TSS_{side}_final'].isna().values)
    assert all(float(f'{a:.4g}') == b for a, b in zip(ex[~np.isnan(ex)], D[f'TSS_{side}_final'].dropna()))
    D[f'TSS_{side}_final'] = ex
    D[f'TSS_{side}_anchoringRule'] = _G.ref_TSS.reindex(D[gcol]).values   # lr_TSS (confirmed) / chromfree_model_TSS (others)
assert ((D.TSS_right_final - D.TSS_left_final)[D.TSS_dist_final.notna()] == D.TSS_dist_final.dropna()).all()
D['TSS_dist_anchoringRule'] = D.TSS_right_anchoringRule - D.TSS_left_anchoringRule

def smooth(v, w=21):
    cs = np.concatenate([[0], np.cumsum(v, dtype=np.float64)]); h = w // 2
    lo = np.clip(np.arange(len(v)) - h, 0, len(v)); hi = np.clip(np.arange(len(v)) + h + 1, 0, len(v))
    return ((cs[hi] - cs[lo]) / (hi - lo)).astype(np.float32)
raw = {tp: {} for tp in TPS}
for tp in TPS:
    bw = pyBigWig.open(NP + f'MNase-seq_{tp}_qNorm.bw')
    for c in CHROMS: raw[tp][c] = np.nan_to_num(np.array(bw.values(c, 0, L[c]), dtype=np.float32))
    bw.close()
SIG = {'mean8': {c: smooth(np.mean([raw[tp][c] for tp in TPS], 0)) for c in CHROMS}}
for tp in TPS: SIG[tp] = {c: smooth(raw[tp][c]) for c in CHROMS}
del raw
MED = {k: float(np.median(np.concatenate([v[c] for c in CHROMS]))) for k, v in SIG.items()}

def call_ndr(sig, thr, minw=50, gap=20):          # identical to ndr_edge.py
    out = {}
    for c in CHROMS:
        b = np.concatenate([[False], sig[c] < thr, [False]]).astype(np.int8)
        d = np.diff(b); s = np.where(d == 1)[0]; e = np.where(d == -1)[0]
        if len(s):
            ks, ke = [s[0]], [e[0]]
            for a, z in zip(s[1:], e[1:]):
                if a - ke[-1] <= gap: ke[-1] = z
                else: ks.append(a); ke.append(z)
            s, e = np.array(ks), np.array(ke); k = (e - s) >= minw; s, e = s[k], e[k]
        out[c] = (s, e)
    return out

def dyad_gap_ndr(tp, mingap=197):                   # identical to ndr_edge.py
    d = pd.read_csv(ND + f'MNase-seq_{tp}.bed', sep='\t', header=None)
    d['dy'] = d[3].str.split('_').str[1].astype(int); out = {}
    for c in CHROMS:
        dy = np.sort(d.dy.values[d[0].values == c]); g = np.diff(dy); k = g > mingap
        out[c] = (dy[:-1][k] + 73, dy[1:][k] - 73)
    return out

def hits(ndr, c, a, b):
    """indices of NDRs [s,e) intersecting closed window [a,b]"""
    s, e = ndr[c]; i1 = np.searchsorted(s, b, side='right'); i0 = np.searchsorted(e, a, side='right')
    return set(range(i0, i1))

def classify(ndr, tl, tr, c, single_between=False, W=300):
    if not (np.isfinite(tl) and np.isfinite(tr)): return np.nan
    tl, tr = int(tl), int(tr)
    hl = hits(ndr, c, tl, tl + W); hr = hits(ndr, c, tr - W, tr)
    if not hl and not hr: return 'none'
    if not hl or not hr: return 'one_side_only'
    if hl & hr: return 'shared'
    if single_between and tl < tr and len(hits(ndr, c, tl, tr)) == 1: return 'shared'
    return 'separate'

def anystage(cols):
    X = pd.concat(cols, axis=1)
    def f(x):
        if x.isna().all(): return np.nan
        for k in ['shared', 'separate', 'one_side_only', 'none']:
            if (x == k).any(): return k
    return X.apply(f, axis=1)

out = D[['gene_left', 'gene_right', 'chrom', 'TSS_left_final', 'TSS_right_final', 'TSS_dist_final', 'TSS_left_anchoringRule', 'TSS_right_anchoringRule', 'TSS_dist_anchoringRule', 'TSS_dist_refinedRule',
         'pearson_DAFT_log2', 'min_DAFT_max', 'cat', 'phase_diff_h', 'ndr']].rename(columns={'ndr': 'ndr_ATAC_refinedRule_anyStage'}).copy()
NDR_P = call_ndr(SIG['mean8'], 0.5 * MED['mean8'])
run = lambda ndr, lc='TSS_left_final', rc='TSS_right_final', sb=False: [classify(ndr, a, b, c, sb) for a, b, c in zip(D[lc], D[rc], D.chrom)]
out['mnase_P'] = run(NDR_P)
out['mnase_P_anchoringRule'] = run(NDR_P, 'TSS_left_anchoringRule', 'TSS_right_anchoringRule')
out['mnase_S1_singleNDRbetween'] = run(NDR_P, sb=True)
for tp in TPS:
    out[f'mnase_S3_{tp}'] = run(call_ndr(SIG[tp], 0.5 * MED[tp]))
    out[f'mnase_S2_dyadgap_{tp}'] = run(dyad_gap_ndr(tp))
out['mnase_S3_anyTP'] = anystage([out[f'mnase_S3_{tp}'] for tp in TPS])
out['mnase_S2_dyadgap_anyTP'] = anystage([out[f'mnase_S2_dyadgap_{tp}'] for tp in TPS])
# mean occupancy between TSSs and NDR count between TSSs
def between(c, a, b, fn):
    if not (np.isfinite(a) and np.isfinite(b)) or b - a < 1: return np.nan
    return fn(c, int(a), int(b))
out['mnase_mean_occ_interTSS'] = [between(c, a, b, lambda c, a, b: float(SIG['mean8'][c][a:b + 1].mean() / MED['mean8'])) for a, b, c in zip(D.TSS_left_final, D.TSS_right_final, D.chrom)]
out['n_NDR_interTSS'] = [between(c, a, b, lambda c, a, b: len(hits(NDR_P, c, a, b))) for a, b, c in zip(D.TSS_left_final, D.TSS_right_final, D.chrom)]
# nucleosome counts between TSSs (reference map dyads = interval midpoints)
for nm, f in [('allNucs', 'NucPosRef_allNucs.bed'), ('top20', 'NucPosRef_top20.bed')]:
    rf = pd.read_csv('ext/nuc/Pf_nucleosome_annotation_of_nucDetective/' + f, sep='\t', header=None); rf['dy'] = (rf[1] + rf[2]) // 2
    DY = {c: np.sort(rf.dy.values[rf[0].values == c]) for c in CHROMS}
    out[f'n_nuc_interTSS_{nm}'] = [between(c, a, b, lambda c, a, b: int(np.searchsorted(DY[c], b, side='right') - np.searchsorted(DY[c], a))) for a, b, c in zip(D.TSS_left_final, D.TSS_right_final, D.chrom)]
out['nuc_per_kb_interTSS_allNucs'] = out.n_nuc_interTSS_allNucs / (out.TSS_dist_final / 1000)

def gap_sep(ndr, tl, tr, c, W=300):
    if not (np.isfinite(tl) and np.isfinite(tr)): return np.nan
    tl, tr = int(tl), int(tr); hl = hits(ndr, c, tl, tl + W); hr = hits(ndr, c, tr - W, tr)
    if not hl or not hr or (hl & hr): return np.nan
    s, e = ndr[c]; return float(min(s[j] - e[i] for i in hl for j in hr if s[j] >= e[i]) if any(s[j] >= e[i] for i in hl for j in hr) else np.nan)
out['gap_between_separate_NDRs_bp'] = [gap_sep(NDR_P, a, b, c) for a, b, c in zip(D.TSS_left_final, D.TSS_right_final, D.chrom)]
out['mnase_S1b_posthoc_gap_lt147'] = np.where((out.mnase_P == 'separate') & (out.gap_between_separate_NDRs_bp < 147), 'shared', out.mnase_P)
out.loc[out.mnase_P.isna(), 'mnase_S1b_posthoc_gap_lt147'] = np.nan
out.to_csv('divergent_pairs_mnase_ndr.tsv', sep='\t', index=False)
pickle.dump(dict(MED=MED, nNDR=sum(len(v[0]) for v in NDR_P.values()), med_width=float(np.median(np.concatenate([v[1] - v[0] for v in NDR_P.values()])))), open('divergent_mnase_stats.pkl', 'wb'))
print('done', out.mnase_P.value_counts(dropna=False).to_dict())
