"""Session-derived computation of benchmark_results_all.tsv.

Recovered from the record of the interactive analysis session in which benchmark_results_all.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: (none recorded)
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
# Verbatim code of the producing analysis cell.
METH = [('Ours_chromatin', 'Ours: CNN, chromatin + DNA + RNA', '#1f5fa8'), ('Ours_nochromatin', 'Ours: CNN, no chromatin', '#7fb2e5'),
        ('PlasmoDB_v71', 'PlasmoDB v71 UTRs (= Chappell 2020 coverage)', '#d08c2e'), ('FixedOffset', 'Fixed UTR length (median)', '#bdbdbd'),
        ('StringTie_short', 'StringTie, stranded short reads', '#8c6bb1'), ('StringTie_long', 'StringTie, long reads', '#c2a5cf'),
        ('ChromatinPeak', 'Chromatin peak (H2A.Z+H3K9ac+ATAC)', '#74c476'), ('ChromHMM', 'ChromHMM promoter states', '#238b45')]
TRUTHS = [('Consensus_5end', 'Consensus of 3 cap/5′-end datasets\n(≥2 agree within 100 bp)'), ('Adjalley2016_capseq', 'Adjalley 2016 cap-seq'), ('Chappell2020_5UTRseq', 'Chappell 2020 5′UTR-seq'),
          ('Shaw2021_5CAPture', 'Shaw 2021 5CAPture'), ('LongRead', 'Long-read reference (training truth)')]
rng = np.random.default_rng(1)
X = Mx.merge(TRb, on='gene_id').merge(GRP, on='gene_id')
def boot(end, truth, m, tol=200, nb=500, sub=None):
    q = X[X[f'{truth}__{end}'].notna()] if sub is None else X[X[f'{truth}__{end}'].notna() & sub]
    h = ((q[f'{m}__{end}'] - q[f'{truth}__{end}']).abs() <= tol).values.astype(float)
    bs = [h[rng.integers(0, len(h), len(h))].mean() for _ in range(nb)]; return h.mean(), np.percentile(bs, 2.5), np.percentile(bs, 97.5), len(h)
CI = pd.DataFrame([dict(end=e, truth=t, method=m, **dict(zip(['acc200', 'lo', 'hi', 'n'], boot(e, t, m))))
                   for e, ts in (('TSS', [t for t, _ in TRUTHS]), ('TTS', ['LongRead'])) for t in ts for m, _, _ in METH if f'{m}__{e}' in X])
print(CI[CI.truth == 'Consensus_5end'].round(3).to_string(index=False))
R0.to_csv('benchmark_results_all.tsv', sep='\t', index=False); CI.to_csv('benchmark_acc200_bootstrap.tsv', sep='\t', index=False)
