"""Recall of v3 truth ends (test chr 1-5, one-to-one, 200 bp) stratified by expression (DAFT-seq, max-stage mean over the CDS) per model,
with calls fixed by validation-tuned threshold/NMS (as in report_fixed)."""
from pfconfig import LAB_DIR  # paths: config/paths.yaml
import os, sys, glob, numpy as np, pandas as pd
sys.path[:0] = ['.', 'boundary_model']
import peak_eval as PE
from bm.evaluate import match
RES = 20
POS = pd.read_csv('data/ml3_end_positives.tsv', sep='\t'); POS = POS[POS.chrom.str[6:8].astype(int) <= 5]
z = np.load('data/daft_stage20.npz')
g = pd.read_csv(glob.glob(f'{LAB_DIR}/PlasmoDB/PlasmoDB-68/*.gff')[0], sep='\t', comment='#', header=None, low_memory=False, names=['chrom','src','type','start','end','score','strand','phase','attr'])
cds = g[g.type == 'CDS'].copy(); cds['gene_id'] = cds.attr.str.extract(r'Parent=([^;.]+)')[0]
env = cds.groupby('gene_id').agg(chrom=('chrom', 'first'), strand=('strand', 'first'), lo=('start', 'min'), hi=('end', 'max'))
def expr(r):
    if r.chrom not in {k.split('__')[1] for k in z.files}: return np.nan
    k = 0 if r.strand == '+' else 1
    return max(float(np.expm1(z[f'{s}__{r.chrom}'][k, r.lo // RES:r.hi // RES + 1].astype(np.float32)).mean()) for s in ('ring', 'troph', 'schizont'))
env['daft'] = [expr(r) for r in env.itertuples()]
POS = POS.merge(env[['daft']], left_on='locus', right_index=True, how='left')
POS['expr_q'] = pd.qcut(POS.daft.rank(method='first'), 4, labels=['Q1 low', 'Q2', 'Q3', 'Q4 high'])
def calls_fixed(names):
    PR = PE.load(names, 'pred'); LG = PE.load(names, 'logit'); va = [c for c in PR if PE.num(c) in (6, 7)]; te = [c for c in PR if PE.num(c) <= 5]; out = {}
    for ch in (0, 1):
        mk = lambda chroms, t, d: {(c, k): PE.call(PR[c][k, ch], LG[c][k, ch][:PR[c].shape[-1]], 'logit_nms', t, d) for c in chroms for k in (0, 1)}
        _, thr, d = max((PE.score(mk(va, t, d), va, ch, 'one_to_one', 200)[2], t, d) for t in np.round(np.arange(0.1, 0.96, 0.05), 2) for d in (60, 100, 160, 300))
        out[ch] = mk(te, thr, d)
    return out
rows = []
for lab, nm in [(a.split('::')[0], a.split('::')[1]) for a in sys.argv[2:]]:
    C = calls_fixed(nm)
    for ch, et in ((0, 'TSS'), (1, 'TTS')):
        P = POS[POS.end_type == et].copy(); P['hit'] = False
        for (c, s), q in P.groupby(['chrom', 'strand']):
            cc = C[ch][(c, 0 if s == '+' else 1)] * RES + RES // 2; _, ti, _ = match(cc, q.pos.values, 200); P.loc[q.index[ti], 'hit'] = True
        for qn, q in P.groupby('expr_q', observed=True): rows.append(dict(model=lab, end_type=et, expr_q=qn, n=len(q), recall200=q.hit.mean(), daft_median=q.daft.median()))
    print(lab, flush=True)
pd.DataFrame(rows).to_csv(sys.argv[1], sep='\t', index=False)
