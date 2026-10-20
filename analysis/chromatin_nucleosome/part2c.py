import sys, os; sys.path.append(os.getcwd())
import numpy as np, pandas as pd, pickle
from statsmodels.stats.multitest import multipletests
from common import *
base = pickle.load(open('part2_base.pkl', 'rb'))['base']
ND = 'ext/nuc/output_nucDetective/Inspector/03_NUCS/'
TPS = ['T5', 'T10', 'T15', 'T20', 'T25', 'T30', 'T35', 'T40']
rng = np.random.default_rng(11)
def load(path, ref=False):
    d = pd.read_csv(path, sep='\t', header=None)
    dy = ((d[1] + d[2]) // 2).values if ref else d[3].str.split('_').str[1].astype(int).values
    return {c: np.sort(dy[d[0].values == c]) for c in d[0].unique()}
sets = {}
for tp in TPS:
    sets[(tp, 'all_calls')] = load(f'{ND}bed_allNucs/MNase-seq_{tp}.bed'); sets[(tp, 'top20_well_positioned')] = load(f'{ND}bed_top20/MNase-seq_{tp}_top20.bed')
sets[('reference_map', 'all_calls')] = load(f'{ND}referenceMap/NucPosRef_allNucs.bed', True); sets[('reference_map', 'top20_well_positioned')] = load(f'{ND}referenceMap/NucPosRef_top20.bed', True)
gsize = sum(chrom_sizes()[c] for c in CHROMS)
ANCH = {'refined': base.ref_TSS.values.astype(int), 'PlasmoDB_v71': base.v71_TSS.values.astype(int)}
W = 1000; K = np.exp(-0.5 * (np.arange(-40, 41) / 10) ** 2); K /= K.sum()

def offsets(cs, X):
    """per gene: sorted array of dyad offsets within +-W (transcription orientation)"""
    out = []
    for c, s, x in zip(base.chrom.values, base.strand.values, X):
        dy = cs.get(c, np.array([], int)); i0, i1 = np.searchsorted(dy, x - W), np.searchsorted(dy, x + W, side='right')
        o = (dy[i0:i1] - x) * (1 if s == '+' else -1); out.append(o)
    return out

def flat(offl):
    gi = np.concatenate([np.full(len(o), i) for i, o in enumerate(offl)]); of = np.concatenate(offl).astype(int)
    return gi, of

def density(F, wts):
    gi, of = F; h = np.bincount(of + W, weights=wts[gi], minlength=2 * W + 1)
    return np.convolve(h / wts.sum(), K, mode='same')  # dyads per bp per gene

OFF = np.arange(-W, W + 1)
def lm(d, bgd):
    def at(lo, hi, f):
        m = (OFF >= lo) & (OFF <= hi); j = np.where(m)[0][f(d[m])]; return int(OFF[j]), d[j]
    p1, h1 = at(-60, 300, np.argmax); tr, ht = at(-400, p1 - 60, np.argmin); m1, hm1 = at(min(p1 - 600, tr - 200), tr - 30, np.argmax)
    half = bgd + (h1 - bgd) / 2; j = np.where(OFF == p1)[0][0]; l = r = j
    while l > 0 and d[l] > half: l -= 1
    while r < len(d) - 1 and d[r] > half: r += 1
    p2, h2 = at(p1 + 120, p1 + 260, np.argmax)
    return dict(plus1_mode=p1, plus1_height_over_bg=h1 / bgd, plus1_FWHM_above_bg=int(OFF[r] - OFF[l]), NDR_trough=tr, NDR_trough_over_bg=ht / bgd,
                minus1_mode=m1, minus1_height_over_bg=hm1 / bgd, NDR_dyad_to_dyad=p1 - m1, NDR_len_edge_to_edge=p1 - m1 - 147, plus2_mode=p2, plus1_to_plus2=p2 - p1, plus2_height_over_bg=h2 / bgd)

def w20(F, wts, mode):
    gi, of = F; hit = np.bincount(gi, weights=(np.abs(of - mode) <= 20).astype(float), minlength=len(wts)) > 0
    return (hit * wts).sum() / wts.sum()

rows = []
for (tp, kind), cs in sets.items():
    nd = sum(len(v) for c, v in cs.items() if c in CHROMS); bgd = nd / gsize
    O = {an: flat(offsets(cs, X)) for an, X in ANCH.items()}
    for grp in GROUPS:
        for sub in ['paired_all', 'v71_has_UTR']:
            m = np.ones(len(base), bool) if grp == 'all' else base.lr_group.values == grp
            if sub == 'v71_has_UTR': m &= base.v71_has5.values
            idx = np.where(m)[0]; wt0 = m.astype(float)
            r = dict(timepoint=tp, call_set=kind, gene_group=grp, gene_subset=sub, n_genes=len(idx), genome_dyad_density_per_kb=bgd * 1000,
                     expected_w20_random=1 - np.exp(-41 * bgd))
            res = {}
            for an in O:
                L = lm(density(O[an], wt0), bgd); L['frac_genes_dyad_within20_of_plus1_mode'] = w20(O[an], wt0, L['plus1_mode']); res[an] = L
                for k, v in L.items(): r[f'{"ref" if an == "refined" else "v71"}_{k}'] = v
            if grp == 'all' or tp in ('T10', 'T30', 'T40', 'reference_map'):
                nb = 200 if grp == 'all' else 100; dh = []; dw = []
                for _ in range(nb):
                    ii = np.bincount(rng.choice(idx, len(idx)), minlength=len(base)).astype(float)
                    la = lm(density(O['refined'], ii), bgd); lb = lm(density(O['PlasmoDB_v71'], ii), bgd)
                    dh.append(la['plus1_height_over_bg'] - lb['plus1_height_over_bg'])
                    dw.append(w20(O['refined'], ii, la['plus1_mode']) - w20(O['PlasmoDB_v71'], ii, lb['plus1_mode']))
                dh, dw = np.array(dh), np.array(dw)
                r.update(diff_plus1_height=res['refined']['plus1_height_over_bg'] - res['PlasmoDB_v71']['plus1_height_over_bg'], diff_plus1_height_CI95_lo=np.percentile(dh, 2.5), diff_plus1_height_CI95_hi=np.percentile(dh, 97.5),
                         p_height_boot=min(1, 2 * min((dh <= 0).mean(), (dh >= 0).mean()) + 1 / nb),
                         diff_w20=res['refined']['frac_genes_dyad_within20_of_plus1_mode'] - res['PlasmoDB_v71']['frac_genes_dyad_within20_of_plus1_mode'],
                         diff_w20_CI95_lo=np.percentile(dw, 2.5), diff_w20_CI95_hi=np.percentile(dw, 97.5), p_w20_boot=min(1, 2 * min((dw <= 0).mean(), (dw >= 0).mean()) + 1 / nb))
            rows.append(r)
    print(tp, kind, flush=True)
D = pd.DataFrame(rows)
for c in ['p_height_boot', 'p_w20_boot']:
    ok = D[c].notna(); D.loc[ok, c.replace('p_', 'BH_q_')] = multipletests(D.loc[ok, c], method='fdr_bh')[1]
D.to_csv('nucleosome_dyad_density_summary.tsv', sep='\t', index=False, float_format='%.4g')
print('done')
