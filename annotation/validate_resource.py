# Validate refined_utr_annotation.gff3 (run after build_resource.py; uses gff, TXI, CDSB from it)
import glob, json
F = pd.read_csv('fig6/refined_utr_annotation.gff3', sep='\t', comment='#', header=None,
                names="chrom src type start end score strand phase attr".split())
F['ID'] = F.attr.str.extract(r'ID=([^;]+)'); F['Parent'] = F.attr.str.extract(r'Parent=([^;]+)')
fai = pd.read_csv(glob.glob(f"{IGV}/*.fai")[0], sep='\t', header=None, index_col=0)[1]
V = {'n_genes': int((F.type == 'protein_coding_gene').sum()), 'features': {k: int(v) for k, v in F.type.value_counts().items()}}
V['start_gt_end'] = int((F.start > F.end).sum())
V['outside_chromosome'] = int(((F.start < 1) | (F.end > F.chrom.map(fai))).sum())
V['duplicate_IDs'] = int(F.ID.duplicated().sum())
mr = F[F.type == 'mRNA'].set_index('ID'); gn = F[F.type == 'protein_coding_gene'].set_index('ID')
ch = F[F.Parent.notna() & (F.type != 'mRNA')]
V['child_strand_or_chrom_mismatch'] = int(((ch.strand.values != mr.loc[ch.Parent, 'strand'].values) |
                                           (ch.chrom.values != mr.loc[ch.Parent, 'chrom'].values)).sum())
V['gene_vs_mRNA_span_mismatch'] = int(((gn.loc[mr.Parent, 'start'].values != mr.start.values) |
                                       (gn.loc[mr.Parent, 'end'].values != mr.end.values)).sum())
V['CDS_phase_invalid'] = int((~F[F.type == 'CDS'].phase.astype(str).isin(['0', '1', '2'])).sum())
cnt = dict(utr5_not_adjacent=0, utr3_not_adjacent=0, utr5_adjacent_via_v71_intron=0, utr3_adjacent_via_v71_intron=0,
           utr_overlaps_own_CDS=0, exon_vs_UTR_CDS_mismatch=0, mRNA_span_vs_exon_mismatch=0, CDS_differs_from_PlasmoDB68=0,
           introns_not_in_PlasmoDB68=0)
def cov(df): return set().union(*[set(range(a, b + 1)) for a, b in zip(df.start, df.end)]) if len(df) else set()
for tx, d in ch.groupby('Parent'):
    s = mr.loc[tx, 'strand']; c = d[d.type == 'CDS']; cs, ce = c.start.min(), c.end.max()
    if sorted(c[['start', 'end']].values.tolist()) != sorted([[a, b] for a, b, _ in CDSB[tx]]): cnt['CDS_differs_from_PlasmoDB68'] += 1
    ivs = set(TXI.get(tx, []))
    for ty, k in (('five_prime_UTR', 'utr5'), ('three_prime_UTR', 'utr3')):
        u = d[d.type == ty]
        if not len(u): continue
        before = (ty == 'five_prime_UTR') == (s == '+')          # UTR lies at lower coordinates than the CDS
        gap = (u.end.max() + 1, cs - 1) if before else (ce + 1, u.start.min() - 1)
        if gap[1] < gap[0]: pass
        elif gap in ivs: cnt[f'{k}_adjacent_via_v71_intron'] += 1
        else: cnt[f'{k}_not_adjacent'] += 1
        if any((a <= ce2) and (b >= cs2) for a, b in zip(u.start, u.end) for cs2, ce2 in zip(c.start, c.end)): cnt['utr_overlaps_own_CDS'] += 1
    e = d[d.type == 'exon'].sort_values('start'); parts = d[d.type != 'exon']
    ce_ = cov(e)
    if ce_ != cov(parts) or len(ce_) != int((parts.end - parts.start + 1).sum()): cnt['exon_vs_UTR_CDS_mismatch'] += 1
    if e.start.min() != mr.loc[tx, 'start'] or e.end.max() != mr.loc[tx, 'end']: cnt['mRNA_span_vs_exon_mismatch'] += 1
    myiv = [(int(e.end.iloc[i]) + 1, int(e.start.iloc[i + 1]) - 1) for i in range(len(e) - 1)]
    cnt['introns_not_in_PlasmoDB68'] += sum(iv not in ivs for iv in myiv)
V.update(cnt)
# QC (not errors): refined UTRs that overlap another gene's CDS (PlasmoDB-68), same or opposite strand
oc = gff[gff.type == 'CDS'][['chrom', 'start', 'end', 'strand', 'gene']]
U = F[F.type.isin(['five_prime_UTR', 'three_prime_UTR']) & F.attr.str.contains('evidence=(?:model|long-read)')].copy()
U['gene'] = U.Parent.str.replace(r'\.\d+$', '', regex=True)
ov = {('five_prime_UTR', 'same'): set(), ('five_prime_UTR', 'opposite'): set(), ('three_prime_UTR', 'same'): set(), ('three_prime_UTR', 'opposite'): set()}
for chrom, uu in U.groupby('chrom'):
    cc = oc[oc.chrom == chrom]; cs_ = cc.start.values; ce_ = cc.end.values
    for r in uu.itertuples():
        hit = cc[(cs_ <= r.end) & (ce_ >= r.start) & (cc.gene.values != r.gene)]
        for st in hit.strand.unique(): ov[(r.type, 'same' if st == r.strand else 'opposite')].add(r.gene)
V['QC_refined_UTR_overlapping_other_gene_CDS'] = {f'{a}_{b}_strand': len(v) for (a, b), v in ov.items()}
json.dump(V, open('fig6/refined_utr_annotation_validation.json', 'w'), indent=1)
