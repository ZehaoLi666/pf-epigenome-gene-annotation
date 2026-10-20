import sys, numpy as np
sys.path[:0] = ['.', 'boundary_model']
import pfprep as P
ends = P.transcript_ends(); refs = {s: P.stage_reference(ends, s) for s in P.STAGES}; excl = P.exclusion_regions()
TR = P.META['chrom_tracks']
Xc = P.normalise_stage({s: P.load_raw('chrom', s) for s in P.STAGES}, TR)
Xr = P.normalise_stage({s: P.load_raw('rna', s) for s in P.STAGES}, P.META['rna_tracks'], quantile_tracks=())
for s in P.STAGES:
    np.savez(f'data/X_chrom_{s}.npz', **Xc[s]); np.savez(f'data/X_rna_{s}.npz', **Xr[s]); refs[s].to_csv(f'data/ref_{s}.tsv', sep='\t', index=False)
excl.to_csv('data/exclusion.tsv', sep='\t', index=False); ends.to_csv('data/transcript_ends.tsv.gz', sep='\t', index=False)
print({s: (refs[s].tier == 'gold').sum() for s in P.STAGES}, len(ends))
