from pfconfig import result  # paths: config/paths.yaml
import sys, os; sys.path.append(os.getcwd())
import pickle, numpy as np, pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests
from common import *
G = load_genes(); bt = pd.read_csv(TRUTH, sep='\t')
mp = pd.read_csv(result('benchmark_method_predictions.tsv'), sep='\t')
d = pickle.load(open('segs.pkl', 'rb')); brown, ours = d['brown'], d['ours']

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

TR = bt[['gene_id', 'Consensus_5end__TSS']].dropna().merge(G[['gene_id', 'chrom', 'strand', 'atg', 'lr_group']], on='gene_id') \
    .merge(mp[['gene_id', 'ChromHMM__TSS', 'Ours_chromatin__TSS', 'Ours_nochromatin__TSS', 'PlasmoDB_v71__TSS']], on='gene_id')
assert len(TR) == 1649
truth = TR.Consensus_5end__TSS.values
preds = {}; meta = {}
cfgs = [('Brown2026_11state', 'TSSassoc 2/5/6/8', brown, {'10hpi': 'ring_10hpi', '30hpi': 'troph_30hpi', '40hpi': 'schizont_40hpi'}, [2, 5, 6, 8]),
        ('Brown2026_11state', 'active 2/5/6', brown, {'10hpi': 'ring_10hpi', '30hpi': 'troph_30hpi', '40hpi': 'schizont_40hpi'}, [2, 5, 6]),
        ('inhouse_ChromHMM_10state', 'promoter-like E8/E9/E10 (H2A.Z/H3K9ac/ATAC-high emissions)', ours, {'ring': 'ring', 'troph': 'troph', 'schizont': 'schizont'}, [8, 9, 10])]
for seg, sdesc, segs, smap, sset in cfgs:
    for rule in ['proximal_edge', 'midpoint', 'distal_edge']:
        SC = []
        for st in segs:
            v = np.array([predict_tss(segs[st], sset, c, s, a, rule) for c, s, a in zip(TR.chrom, TR.strand, TR.atg)])
            k = f'{seg}|{sdesc}|{smap[st]}|{rule}'; preds[k] = v; meta[k] = (seg, sdesc, smap[st], rule); SC.append(v)
        SC = np.vstack(SC)
        with np.errstate(all='ignore'):
            k = f'{seg}|{sdesc}|median_of_3_stages|{rule}'; preds[k] = np.nanmedian(SC, 0); meta[k] = (seg, sdesc, 'median_of_3_stages', rule)
        up = np.where(TR.strand.values == '+', TR.atg.values - SC, SC - TR.atg.values); up = np.where(np.isnan(up), np.inf, up)
        k = f'{seg}|{sdesc}|stage_call_nearest_ATG|{rule}'; preds[k] = SC[np.argmin(up, 0), np.arange(SC.shape[1])]; meta[k] = (seg, sdesc, 'stage_call_nearest_ATG', rule)
for k, col, desc in [('Ours_chromatin', 'Ours_chromatin__TSS', 'held-out CNN with chromatin (3-fold chromosome CV)'),
                     ('Ours_nochromatin', 'Ours_nochromatin__TSS', 'held-out CNN without chromatin'),
                     ('inhouse_ChromHMM_10state_original_baseline', 'ChromHMM__TSS', 'original benchmark ChromHMM call (benchmark_method_predictions.tsv)'),
                     ('PlasmoDB_v71', 'PlasmoDB_v71__TSS', 'annotated TSS')]:
    preds[k] = TR[col].values; meta[k] = (k, desc, '', '')
ours_hit = np.abs(preds['Ours_chromatin'] - truth) <= 200; old_hit = np.abs(preds['inhouse_ChromHMM_10state_original_baseline'] - truth) <= 200
rng = np.random.default_rng(1); n = len(TR); BI = rng.integers(0, n, (2000, n))
rows = []
for grp in ['all', 'confirmed', 'truncated-only', 'HC model, unconfirmed', 'no HC model']:
    m = np.ones(n, bool) if grp == 'all' else TR.lr_group.values == grp
    bi = BI if grp == 'all' else rng.integers(0, m.sum(), (2000, m.sum()))
    for k, p in preds.items():
        if grp != 'all' and not (k.startswith('Ours_chromatin') or k.startswith('inhouse_ChromHMM_10state_original') or k == 'PlasmoDB_v71' or 'median_of_3_stages' in k or 'stage_call_nearest_ATG' in k):
            continue
        pp = p[m]; dd = np.abs(pp - truth[m]); hit = dd <= 200
        r = dict(gene_group=grp, method=k, segmentation=meta[k][0], states_or_description=meta[k][1], stage=meta[k][2], rule=meta[k][3], n_genes=int(m.sum()), coverage=np.isfinite(pp).mean())
        for t in (50, 100, 200, 500): r[f'acc@{t}'] = np.mean(dd <= t)
        bs = hit[bi].mean(1); r['acc@200_CI95_lo'], r['acc@200_CI95_hi'] = np.percentile(bs, [2.5, 97.5])
        r['median_abs_err_bp'] = np.nanmedian(dd)
        for ref, rh in [('vs_Ours_chromatin', ours_hit[m]), ('vs_inhouse_ChromHMM_baseline', old_hit[m])]:
            n10 = int((hit & ~rh).sum()); n01 = int((~hit & rh).sum())
            r[f'{ref}_diff_acc@200'] = hit.mean() - rh.mean(); r[f'{ref}_mcnemar_p'] = binomtest(n10, n10 + n01).pvalue if n10 + n01 else 1.0
        rows.append(r)
B = pd.DataFrame(rows)
for c in ['vs_Ours_chromatin_mcnemar_p', 'vs_inhouse_ChromHMM_baseline_mcnemar_p']: B[c.replace('_p', '_BH_q')] = multipletests(B[c], method='fdr_bh')[1]
prim = 'Brown2026_11state|TSSassoc 2/5/6/8|median_of_3_stages|proximal_edge'
B['note'] = ''
B.loc[B.method == prim, 'note'] = 'PRIMARY rule: ATG-proximal edge of nearest TSS-associated (2/5/6/8) segment within 3 kb upstream of ATG; median of 10/30/40 hpi calls'
allb = B[(B.gene_group == 'all') & B.segmentation.eq('Brown2026_11state')]
best = allb.sort_values('acc@200').method.iloc[-1]
B.loc[B.method == best, 'note'] = f'best of {len(allb)} Brown rule/state-set/stage combinations, selected post hoc (optimistic)'
B.to_csv('hmm_state_tss_predictor_benchmark.tsv', sep='\t', index=False, float_format='%.4f')
print(B[B.gene_group == 'all'].sort_values('acc@200', ascending=False)[['method', 'coverage', 'acc@200', 'acc@200_CI95_lo', 'acc@200_CI95_hi', 'median_abs_err_bp']].head(12).round(3).to_string())
print(best)
