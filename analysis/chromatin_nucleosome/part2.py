import sys, os; sys.path.append(os.getcwd())
import numpy as np, pandas as pd, pyBigWig, pickle
from statsmodels.stats.multitest import multipletests
from common import *
G = load_genes()
ND = 'ext/nuc/output_nucDetective/Inspector/03_NUCS/'
TPS = ['T5', 'T10', 'T15', 'T20', 'T25', 'T30', 'T35', 'T40']
rng = np.random.default_rng(7)

def load_calls(path):
    d = pd.read_csv(path, sep='\t', header=None, names=['chrom', 's', 'e', 'name', 'fuzz', 'x'])
    d['dyad'] = d.name.str.split('_').str[1].astype(int)
    out = {}
    for c, q in d.groupby('chrom'):
        q = q.sort_values('dyad'); out[c] = (q.dyad.values, q.fuzz.values)
    return out
calls = {tp: load_calls(f'{ND}bed_allNucs/MNase-seq_{tp}.bed') for tp in TPS}
ref = pd.read_csv(f'{ND}referenceMap/NucPosRef_allNucs.bed', sep='\t', header=None, names=['chrom', 's', 'e', 'a', 'b', 'c'])
ref['dyad'] = ((ref.s + ref.e) // 2)
calls['reference_map'] = {c: (np.sort(q.dyad.values), np.zeros(len(q))) for c, q in ref.groupby('chrom')}

base = G[G.ref_TSS.notna()].reset_index(drop=True)
s1 = pd.read_excel('ext/nuc/supp/pcbi.1014557.s005.xlsx', sheet_name='Plus1_Nucleosome_Annotation')
s1 = s1.rename(columns={'gene': 'gene_id', 'start': 'holz_plus1'})[['gene_id', 'holz_plus1']].drop_duplicates('gene_id')
base = base.merge(s1, on='gene_id', how='left')
sgn = np.where(base.strand.values == '+', 1, -1)
ANCH = {'refined': base.ref_TSS.values.astype(int), 'PlasmoDB_v71': base.v71_TSS.values.astype(int)}

def nuc_features(cs, chrom, strand, x, rule):
    """+1: dyad within [-100,+300] of anchor (transcription orientation); rule 'best' = lowest DANPOS fuzziness (Holzinger rule),
    'nearest' = smallest |offset|. -1: nearest dyad upstream of +1 dyad with spacing >= 120 bp, within 1 kb."""
    if chrom not in cs: return (np.nan,) * 3
    dy, fz = cs[chrom]
    lo, hi = (x - 100, x + 300) if strand == '+' else (x - 300, x + 100)
    i0, i1 = np.searchsorted(dy, lo), np.searchsorted(dy, hi, side='right')
    if i1 <= i0: return (np.nan,) * 3
    cand = np.arange(i0, i1); off = (dy[cand] - x) * (1 if strand == '+' else -1)
    if rule == 'best' and np.any(fz[cand] > 0): k = cand[np.argmin(fz[cand])]
    else: k = cand[np.argmin(np.abs(off))]
    p1 = dy[k]
    if strand == '+':
        up = np.where((dy < p1 - 120) & (dy >= p1 - 1000))[0]; m1 = dy[up[-1]] if len(up) else np.nan
    else:
        up = np.where((dy > p1 + 120) & (dy <= p1 + 1000))[0]; m1 = dy[up[0]] if len(up) else np.nan
    o1 = (p1 - x) * (1 if strand == '+' else -1); om1 = (m1 - x) * (1 if strand == '+' else -1)
    return o1, om1, abs(p1 - m1) - 147 if np.isfinite(m1) else np.nan

def mode_kde(v, bw=10):
    h = np.bincount(np.clip(np.round(v).astype(int) + 100, 0, 400), minlength=401).astype(float)
    k = np.exp(-0.5 * (np.arange(-40, 41) / bw) ** 2); dens = np.convolve(h, k, mode='same')
    return int(np.argmax(dens)) - 100

rows = []; per_gene = []
for tp in TPS + ['reference_map']:
    cs = calls[tp]
    for rule in (['best', 'nearest'] if tp != 'reference_map' else ['nearest']):
        F = {an: np.array([nuc_features(cs, c, s, x, rule) for c, s, x in zip(base.chrom, base.strand, X)]) for an, X in ANCH.items()}
        if rule == 'best' or tp == 'reference_map':
            for an in F:
                per_gene.append(pd.DataFrame({'gene_id': base.gene_id, 'timepoint': tp, 'rule': rule, 'anchor': an, 'plus1_offset': F[an][:, 0], 'minus1_offset': F[an][:, 1], 'NDR_len': F[an][:, 2]}))
        for grp in GROUPS:
            for sub in ['paired_all', 'v71_has_UTR']:
                m = np.ones(len(base), bool) if grp == 'all' else base.lr_group.values == grp
                if sub == 'v71_has_UTR': m &= base.v71_has5.values
                ok = m & np.isfinite(F['refined'][:, 0]) & np.isfinite(F['PlasmoDB_v71'][:, 0])
                if ok.sum() < 30: continue
                res = {}
                for an in F:
                    o = F[an][ok, 0]; mo = mode_kde(o)
                    res[an] = dict(mode=mo, med=np.median(o), sd=np.std(o, ddof=1), mad=np.median(np.abs(o - np.median(o))), w20=np.mean(np.abs(o - mo) <= 20),
                                   m1=np.nanmedian(F[an][ok, 1]), ndr=np.nanmedian(F[an][ok, 2]), n_m1=int(np.isfinite(F[an][ok, 1]).sum()))
                # paired bootstrap of SD and within-20 differences (mode re-estimated per replicate)
                oa, ob = F['refined'][ok, 0], F['PlasmoDB_v71'][ok, 0]; nb = 300; dsd = []; dw = []
                for _ in range(nb):
                    ii = rng.integers(0, ok.sum(), ok.sum()); a, b = oa[ii], ob[ii]
                    dsd.append(np.std(a, ddof=1) - np.std(b, ddof=1)); dw.append(np.mean(np.abs(a - mode_kde(a)) <= 20) - np.mean(np.abs(b - mode_kde(b)) <= 20))
                dsd, dw = np.array(dsd), np.array(dw)
                r = dict(timepoint=tp, plus1_rule=rule, gene_group=grp, gene_subset=sub, n_genes_both_plus1=int(ok.sum()), n_genes_in_stratum=int(m.sum()))
                for an, k in [('refined', 'ref'), ('PlasmoDB_v71', 'v71')]:
                    for f, v in res[an].items(): r[f'{k}_{f}'] = v
                r.update(diff_sd_ref_minus_v71=res['refined']['sd'] - res['PlasmoDB_v71']['sd'], diff_sd_CI95_lo=np.percentile(dsd, 2.5), diff_sd_CI95_hi=np.percentile(dsd, 97.5), p_sd_boot=min(1, 2 * min((dsd >= 0).mean(), (dsd <= 0).mean()) + 1 / nb),
                         diff_w20_ref_minus_v71=res['refined']['w20'] - res['PlasmoDB_v71']['w20'], diff_w20_CI95_lo=np.percentile(dw, 2.5), diff_w20_CI95_hi=np.percentile(dw, 97.5), p_w20_boot=min(1, 2 * min((dw >= 0).mean(), (dw <= 0).mean()) + 1 / nb))
                rows.append(r)
    print(tp, flush=True)
T = pd.DataFrame(rows)
for c in ['p_sd_boot', 'p_w20_boot']: T[c.replace('p_', 'BH_q_')] = multipletests(T[c], method='fdr_bh')[1]
T.to_csv('nucleosome_tss_positions_summary.tsv', sep='\t', index=False, float_format='%.4g')
PG = pd.concat(per_gene); PG.to_csv('nucleosome_per_gene_offsets.tsv.gz', sep='\t', index=False, float_format='%.0f')
pickle.dump(dict(base=base), open('part2_base.pkl', 'wb'))
print('done')
