"""Labels v3, step 1: gene set and same-strand gaps (CAGE not used anywhere).
Confirmation per stage = detected & struct_ok & sj_<stage> & tes_ok (the stage-confirmation rule of the long-read reference minus the CAGE term).
Full-CDS model = model span covers the v68 CDS envelope of its locus. Gene set = loci with >=1 confirmed full-CDS mRNA model.
Gaps: as in the v2 flank analysis (nearest same-strand long-read unit on each side, <= 10 kb)."""
from pfconfig import LAB_DIR  # paths: config/paths.yaml
import glob, numpy as np, pandas as pd, sys
TSC, LAB = sys.argv[1], LAB_DIR
S = ['ring', 'troph', 'schizont']
tsc = pd.read_csv(TSC, sep='\t', low_memory=False)
for s in S: tsc[f'conf3_{s}'] = (tsc[f'detect_{s}'] == 'detected') & tsc.struct_ok & tsc[f'sj_{s}'] & tsc.tes_ok
tsc['conf3_any'] = tsc[[f'conf3_{s}' for s in S]].any(axis=1)
g = pd.read_csv(glob.glob(f'{LAB}/PlasmoDB/PlasmoDB-68/*.gff')[0], sep='\t', comment='#', header=None, low_memory=False, names=['chrom','src','type','start','end','score','strand','phase','attr'])
cds = g[g.type == 'CDS'].copy(); cds['gene'] = cds.attr.str.extract(r'Parent=([^;.]+)')[0]
env = cds.groupby('gene').agg(cds_lo=('start', 'min'), cds_hi=('end', 'max')); env['cds_lo'] -= 1
T = tsc.join(env, on='locus'); mr = T.unit_class == 'mRNA'
T['full_cds'] = mr & T.cds_lo.notna() & (T.start <= T.cds_lo) & (T.end >= T.cds_hi)
T['trunc5'] = mr & T.cds_lo.notna() & ~T.full_cds & np.where(T.strand == '+', T.start > T.cds_lo, T.end < T.cds_hi)
T['trunc3'] = mr & T.cds_lo.notna() & ~T.full_cds & ~T.trunc5
T.to_csv('v3/models_v3.tsv.gz', sep='\t', index=False)
F = T[T.full_cds & T.conf3_any]
foc = F.groupby('locus').agg(chrom=('chrom', 'first'), strand=('strand', 'first'), start=('start', 'min'), end=('end', 'max'), n_models=('transcript_id', 'size'))
grp = lambda u: 'mRNA' if u == 'mRNA' else ('lncRNA' if u.startswith('lncRNA') else ('novel ORF' if u.startswith('novel') else 'other ncRNA'))
tsc['g'] = tsc.unit_class.map(grp)
U = tsc.groupby(['locus', 'strand']).agg(chrom=('chrom', 'first'), start=('start', 'min'), end=('end', 'max'), g=('g', 'first')).reset_index().set_index('locus')
U.loc[foc.index, ['start', 'end']] = foc[['start', 'end']].values; U['confirmed_gene'] = U.index.isin(foc.index); U = U.reset_index()
Ui = {k: v.sort_values('start') for k, v in U.groupby(['chrom', 'strand'])}; rows = []
for loc, gg in foc.iterrows():
    q = Ui[(gg.chrom, gg.strand)]; q = q[q.locus != loc]
    left = q[q.start < gg.start]; left = left.loc[[left.end.idxmax()]] if len(left) else left
    right = q[q.end > gg.end]; right = right.loc[[right.start.idxmin()]] if len(right) else right
    for side_g, nb in (('L', left), ('R', right)):
        side = ('5p' if side_g == 'L' else '3p') if gg.strand == '+' else ('3p' if side_g == 'L' else '5p')
        row = dict(gene=loc, side=side, chrom=gg.chrom, strand=gg.strand)
        if len(nb) == 0: row.update(status='no neighbour'); rows.append(row); continue
        nb = nb.iloc[0]; a, b = (nb.end, gg.start) if side_g == 'L' else (gg.end, nb.start); gap = b - a
        row.update(neighbour=nb.locus, neighbour_class=nb.g, a=int(a), b=int(b), gap=int(gap), status='overlap' if gap <= 0 else ('far (>10 kb)' if gap > 10000 else 'gap'))
        row.update(tts_pos=int(a - 1), tss_pos=int(b)) if gg.strand == '+' else row.update(tts_pos=int(b), tss_pos=int(a - 1))
        rows.append(row)
GP = pd.DataFrame(rows); gg_ = GP[GP.status == 'gap']
UG = gg_.drop_duplicates(['chrom', 'strand', 'a', 'b'])[['chrom', 'strand', 'a', 'b', 'tts_pos', 'tss_pos']].reset_index(drop=True)
UG.insert(0, 'gap_id', [f'g3_{i:05d}' for i in range(len(UG))]); UG[['a', 'b', 'tts_pos', 'tss_pos']] = UG[['a', 'b', 'tts_pos', 'tss_pos']].astype(int)
GP = GP.merge(UG, on=['chrom', 'strand', 'a', 'b', 'tts_pos', 'tss_pos'], how='left')
UG.to_csv('v3/gaps_v3.tsv', sep='\t', index=False); GP.to_csv('v3/gene_side_gaps_v3.tsv', sep='\t', index=False)
print('genes', len(foc), 'full-CDS models', len(F), 'sides', GP.status.value_counts().to_dict(), 'unique gaps', len(UG))
