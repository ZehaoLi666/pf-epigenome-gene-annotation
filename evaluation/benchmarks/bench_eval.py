"""Gene-level TSS/TTS benchmark. Inputs: methods table (gene_id + <method>__TSS / <method>__TTS columns, NaN = no call) and truth table
(gene_id + <truth>__TSS / <truth>__TTS). Metrics per method x truth x gene group:
  n_truth: genes with a truth end; coverage: fraction of those with a call; acc@d: fraction within d bp (no call counts as miss);
  median_err: median |call - truth| over genes with a call."""
import sys, numpy as np, pandas as pd
TOL = (50, 100, 200, 500)
def evaluate(M, T, groups):
    meth = sorted({c.split('__')[0] for c in M if '__' in c}); tru = sorted({c.split('__')[0] for c in T if '__' in c})
    X = M.merge(T, on='gene_id', how='inner').merge(groups, on='gene_id', how='left'); rows = []
    for end in ('TSS', 'TTS'):
        for t in tru:
            tc = f'{t}__{end}'
            if tc not in X: continue
            for g, q0 in [('all', X)] + list(X.groupby('group')):
                q = q0[q0[tc].notna()]
                if len(q) < 10: continue
                for m in meth:
                    mc = f'{m}__{end}'
                    if mc not in q: continue
                    d = (q[mc] - q[tc]).abs(); r = dict(end=end, truth=t, group=g, method=m, n_truth=len(q), coverage=q[mc].notna().mean(), median_err=d.median())
                    for tol in TOL: r[f'acc@{tol}'] = (d <= tol).mean()
                    rows.append(r)
    return pd.DataFrame(rows)
if __name__ == '__main__':
    M = pd.read_csv(sys.argv[1], sep='\t'); T = pd.read_csv(sys.argv[2], sep='\t'); G = pd.read_csv(sys.argv[3], sep='\t')
    evaluate(M, T, G).to_csv(sys.argv[4], sep='\t', index=False)
