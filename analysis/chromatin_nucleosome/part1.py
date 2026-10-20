import numpy as np, pandas as pd
from scipy.stats import binomtest, wilcoxon
from statsmodels.stats.multitest import multipletests
import sys, os; sys.path.append(os.getcwd()); from common import *

L = chrom_sizes(); G = load_genes()
bt = pd.read_csv(TRUTH, sep='\t')

def load_seg(path, statefun=int):
    arr = {c: np.zeros(L[c] // 200 + 1, dtype=np.int16) for c in CHROMS}
    for line in open(path):
        if line.startswith('track'): continue
        f = line.split('\t')
        if f[0] not in arr: continue
        s, e = int(f[1]), int(f[2]); assert s % 200 == 0
        arr[f[0]][s // 200:(e + 199) // 200] = statefun(f[3])
    return arr

brown = {st: load_seg(f'ext/brown/chrom-states-IGV/{st}_chrom-states.bed') for st in ['10hpi', '30hpi', '40hpi']}
ours = {st: load_seg(f'bench/bench_pkg/chromhmm_10state/{st}_10_segments.bed', lambda x: int(x.strip().lstrip('E'))) for st in ['ring', 'troph', 'schizont']}
stage_lab = {'10hpi': 'ring_10hpi', '30hpi': 'troph_30hpi', '40hpi': 'schizont_40hpi'}
TSSSTATES = [2, 5, 6, 8]; ACTIVE = [2, 5, 6]
NAMES = {1: 'OpenChromatin', 2: 'Active', 3: 'CodingStart', 4: 'GeneBody', 5: 'WeakActive', 6: 'RingActive', 7: 'HP1', 8: 'H2A.Z', 9: 'H3K4me1', 10: 'H3K9me3', 11: 'Quiescent'}
OFFS = np.arange(-1500, 1501, 50); i0 = int(np.where(OFFS == 0)[0][0])

def states_at(arr, chrom, x, strand, offs=OFFS):
    out = np.zeros((len(x), len(offs)), dtype=np.int16)
    for i, (c, p, s) in enumerate(zip(chrom, x, strand)):
        pos = p + offs if s == '+' else p - offs
        a = arr[c]; pos = np.clip(pos, 0, (len(a) - 1) * 200)
        out[i] = a[(pos // 200).astype(int)]
    return out

# ---------------- (a) composition ----------------
WIN = {'point0': (0, 0), 'core_-100_+100': (-100, 100), 'up_-1000_-200': (-1000, -200), 'down_+200_+1000': (200, 1000)}
rows = []; S = {}; BASE = {}
for end in ['TSS', 'TTS']:
    base = G[G[f'ref_{end}'].notna()].reset_index(drop=True); BASE[end] = base
    anchors = {'refined': base[f'ref_{end}'].values.astype(int), 'PlasmoDB_v71': base[f'v71_{end}'].values.astype(int)}
    for an, x in anchors.items():
        for st in brown:
            M = states_at(brown[st], base.chrom.values, x, base.strand.values); S[(end, an, st)] = M
            for grp in GROUPS:
                for sub in ['paired_all', 'v71_has_UTR']:
                    msk = np.ones(len(base), bool) if grp == 'all' else (base.lr_group.values == grp)
                    if sub == 'v71_has_UTR': msk &= base['v71_has5' if end == 'TSS' else 'v71_has3'].values
                    for w, (a, b) in WIN.items():
                        v = M[msk][:, (OFFS >= a) & (OFFS <= b)].ravel(); cnt = np.bincount(v, minlength=12)[1:]
                        r = dict(end=end, anchor=an, anchor_rule=(f'lr_{end} (confirmed) / chromfree_model_{end} (others with refined {end})' if an == 'refined' else f'annot_{end.lower()} (v71)'),
                                 stage=stage_lab[st], gene_group=grp, gene_subset=sub, window=w, n_genes=int(msk.sum()))
                        for k in range(11): r[f'S{k+1}_{NAMES[k+1]}'] = cnt[k] / cnt.sum()
                        r['TSSassoc_2_5_6_8'] = cnt[[1, 4, 5, 7]].sum() / cnt.sum(); r['active_2_5_6'] = cnt[[1, 4, 5]].sum() / cnt.sum()
                        rows.append(r)
for st in brown:
    c = np.bincount(np.concatenate([a[:-1] for a in brown[st].values()]), minlength=12)[1:]
    r = dict(end='genome', anchor='genome_background', anchor_rule='all 200-bp bins chr1-14', stage=stage_lab[st], gene_group='all', gene_subset='genome', window='genome', n_genes=np.nan)
    for k in range(11): r[f'S{k+1}_{NAMES[k+1]}'] = c[k] / c.sum()
    r['TSSassoc_2_5_6_8'] = c[[1, 4, 5, 7]].sum() / c.sum(); r['active_2_5_6'] = c[[1, 4, 5]].sum() / c.sum(); rows.append(r)
comp = pd.DataFrame(rows)
comp.insert(0, 'source', 'Brown2026_BMCGenomics_ChromHMM11 (github seqcode/pfal-chrom-states-trackhub, chrom-states-IGV, Pf3D7_v3 names)')
comp.to_csv('hmm_state_tss_tts_composition.tsv', sep='\t', index=False, float_format='%.4f')

# ---------------- (b) occupancy tests ----------------
rng = np.random.default_rng(0)
def paired_bin(a, b, nb=2000):
    n = len(a); n10 = int((a & ~b).sum()); n01 = int((~a & b).sum())
    p = binomtest(n10, n10 + n01, 0.5).pvalue if n10 + n01 else 1.0
    idx = rng.integers(0, n, (nb, n)); d = a[idx].mean(1) - b[idx].mean(1)
    return a.mean(), b.mean(), a.mean() - b.mean(), np.percentile(d, 2.5), np.percentile(d, 97.5), n10, n01, p
base = BASE['TSS']; sgn = np.where(base.strand == '+', 1, -1)
upm = (OFFS >= -600) & (OFFS <= -100); dnm = (OFFS >= 100) & (OFFS <= 600)
strata = [('gene_group', x) for x in GROUPS] + [('DAFT_quintile', q) for q in sorted(base.Qd.dropna().unique())]
def mask(kind, val, sub):
    m = np.ones(len(base), bool) if val == 'all' else (base.lr_group.values == val if kind == 'gene_group' else base.Qd.values == val)
    if sub == 'v71_has_UTR': m &= base.v71_has5.values
    return m
out = []
for metric, sset in [('in_TSSassoc_state(2/5/6/8)_at_TSS_bp', TSSSTATES), ('in_active_state(2/5/6)_at_TSS_bp', ACTIVE)]:
    for stg in list(brown) + ['any_stage']:
        if stg == 'any_stage':
            A = np.any([np.isin(S[('TSS', 'refined', s)][:, i0], sset) for s in brown], 0); B = np.any([np.isin(S[('TSS', 'PlasmoDB_v71', s)][:, i0], sset) for s in brown], 0)
        else:
            A = np.isin(S[('TSS', 'refined', stg)][:, i0], sset); B = np.isin(S[('TSS', 'PlasmoDB_v71', stg)][:, i0], sset)
        for sub in ['paired_all', 'v71_has_UTR']:
            for kind, val in strata:
                m = mask(kind, val, sub)
                fa, fb, d, lo, hi, n10, n01, p = paired_bin(A[m], B[m])
                out.append(dict(metric=metric, stage=stage_lab.get(stg, stg), stratum_type=kind, stratum=val, gene_subset=sub, n=int(m.sum()), refined=fa, v71=fb, diff=d, diff_CI95_lo=lo, diff_CI95_hi=hi, n_refined_only=n10, n_v71_only=n01, test='McNemar exact', p=p))

def predict_tss(arr, states, chrom, strand, atg, rule, W=3000):
    a = arr[chrom]
    if strand == '+':
        lo = max(atg - W, 0); bins = np.arange(atg // 200, lo // 200 - 1, -1)
    else:
        hi = min(atg + W, (len(a) - 1) * 200); bins = np.arange(atg // 200, hi // 200 + 1)
    ins = np.isin(a[bins], states)
    if not ins.any(): return np.nan
    k = int(np.argmax(ins)); j = k
    while j + 1 < len(bins) and ins[j + 1]: j += 1
    if strand == '+': prox = min(bins[k] * 200 + 200, atg); dist = bins[j] * 200
    else: prox = max(bins[k] * 200, atg); dist = bins[j] * 200 + 200
    return {'proximal_edge': prox, 'distal_edge': dist, 'midpoint': (prox + dist) / 2}[rule]

for st in brown:
    T = {an: np.isin(S[('TSS', an, st)], TSSSTATES) for an in ['refined', 'PlasmoDB_v71']}
    con = {an: T[an][:, upm].mean(1) - T[an][:, dnm].mean(1) for an in T}
    edge = np.array([predict_tss(brown[st], TSSSTATES, c, s, a, 'proximal_edge') for c, s, a in zip(base.chrom, base.strand, base.atg)])
    eo = {an: (base[col].values - edge) * sgn for an, col in [('refined', 'ref_TSS'), ('PlasmoDB_v71', 'v71_TSS')]}
    for sub in ['paired_all', 'v71_has_UTR']:
        for kind, val in strata:
            m = mask(kind, val, sub); a, b = con['refined'][m], con['PlasmoDB_v71'][m]
            idx = rng.integers(0, len(a), (2000, len(a))); bs = (a - b)[idx].mean(1)
            out.append(dict(metric='boundary_contrast: TSSassoc fraction in [-600,-100] minus [+100,+600] (per gene, mean)', stage=stage_lab[st], stratum_type=kind, stratum=val, gene_subset=sub, n=int(m.sum()),
                            refined=a.mean(), v71=b.mean(), diff=(a - b).mean(), diff_CI95_lo=np.percentile(bs, 2.5), diff_CI95_hi=np.percentile(bs, 97.5), test='Wilcoxon signed-rank', p=wilcoxon(a, b).pvalue))
            ea, eb = eo['refined'][m], eo['PlasmoDB_v71'][m]; ok = np.isfinite(ea) & np.isfinite(eb)
            fa, fb, d, lo, hi, n10, n01, p = paired_bin(np.abs(ea[ok]) <= 200, np.abs(eb[ok]) <= 200)
            out.append(dict(metric='TSS within 200 bp of ATG-proximal edge of nearest TSSassoc segment (<=3 kb upstream of ATG)', stage=stage_lab[st], stratum_type=kind, stratum=val, gene_subset=sub, n=int(ok.sum()),
                            refined=fa, v71=fb, diff=d, diff_CI95_lo=lo, diff_CI95_hi=hi, n_refined_only=n10, n_v71_only=n01, test='McNemar exact', p=p,
                            median_offset_to_edge_refined_bp=np.median(ea[ok]), median_offset_to_edge_v71_bp=np.median(eb[ok])))
OCC = pd.DataFrame(out); OCC['BH_q'] = multipletests(OCC.p, method='fdr_bh')[1]

# calibration against experimental 5'-end TSSs (Consensus_5end) on the same genes
TR = bt[['gene_id', 'Consensus_5end__TSS']].dropna().merge(G[['gene_id', 'chrom', 'strand', 'atg', 'lr_group', 'ref_TSS', 'v71_TSS']], on='gene_id')
cal = []
for st in brown:
    ok = TR.ref_TSS.notna().values
    for nm, col in [('Consensus_5end_truth', 'Consensus_5end__TSS'), ('refined', 'ref_TSS'), ('PlasmoDB_v71', 'v71_TSS')]:
        Mx = states_at(brown[st], TR.chrom.values[ok], TR[col].values[ok].astype(int), TR.strand.values[ok]); T = np.isin(Mx, TSSSTATES); A = np.isin(Mx, ACTIVE)
        for metric, val in [('in_TSSassoc_state(2/5/6/8)_at_TSS_bp', T[:, i0].mean()), ('in_active_state(2/5/6)_at_TSS_bp', A[:, i0].mean()),
                            ('boundary_contrast: TSSassoc fraction in [-600,-100] minus [+100,+600] (per gene, mean)', (T[:, upm].mean(1) - T[:, dnm].mean(1)).mean())]:
            cal.append(dict(metric=metric, stage=stage_lab[st], stratum_type='calibration_Consensus_5end_genes', stratum=nm, gene_subset='Consensus_5end genes with refined TSS', n=int(ok.sum()), refined=val if nm == 'refined' else np.nan,
                            v71=val if nm == 'PlasmoDB_v71' else np.nan, truth_value=val if nm.startswith('Consensus') else np.nan))
OCC = pd.concat([OCC, pd.DataFrame(cal)], ignore_index=True)
OCC.to_csv('hmm_state_tss_occupancy_tests.tsv', sep='\t', index=False, float_format='%.5g')
import pickle; pickle.dump(dict(brown=brown, ours=ours), open('segs.pkl', 'wb'))
print('done', comp.shape, OCC.shape)
