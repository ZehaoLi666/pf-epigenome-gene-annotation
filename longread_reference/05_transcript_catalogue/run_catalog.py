# expects in namespace: TC, GN, Mi, G, ovl_cds, cage_at, indep_units, antisense_units, other_gene_overlap, lnc_ov, flip, np, pd
def chimeric(g, a, b):
    lo, hi = min(g.cds5, g.cds3), max(g.cds5, g.cds3) + 1
    pa, pb = [], []
    for x, y in zip(a, b):           # exon parts outside the gene's own CDS span
        if x < lo: pa.append(x); pb.append(min(y, lo))
        if y > hi: pa.append(max(x, hi)); pb.append(y)
    return ovl_cds(g.chrom, g.strand, pa, pb) if pa else False
rows = []
for g in GN.itertuples(index=False):
    r = dict(transcript=g.transcript, gene=g.gene, chrom=g.chrom, strand=g.strand, cds5=g.cds5, cds3=g.cds3, L5=g.L5 if g.has5 else 0, L3=g.L3 if g.has3 else 0)
    for side, L in ((5, r['L5']), (3, r['L3'])):
        cl, nds, rd, nch, ext = TC.ends_side(Mi, g, side, chimeric=chimeric); P = f'p{side}_'
        sup = [x for x in cl if TC.supported(x)]
        r[P + 'n_ds'] = nds; r[P + 'reads'] = rd; r[P + 'chimeric_reads'] = nch; r[P + 'n_clusters_supported'] = len(sup)
        a, b = TC.interval(g, side, 0, max(L, 1))
        r[P + 'other_gene_same_strand'] = other_gene_overlap(g.chrom, g.strand, a, b, g.gene) if L > 0 else ''
        r[P + 'other_gene_opposite'] = other_gene_overlap(g.chrom, flip[g.strand], a, b, g.gene) if L > 0 else ''
        r[P + 'lnc_sense'] = lnc_ov(g.chrom, g.strand, a, b) if L > 0 else False
        r[P + 'lnc_antisense'] = lnc_ov(g.chrom, flip[g.strand], a, b) if L > 0 else False
        if nds < 2 or not sup: r[P + 'cat'] = 'Insufficient long-read models'; continue
        major = max(sup, key=lambda x: x['reads']); distal = max(sup, key=lambda x: x['peak'])
        ud = ext; r[P + 'major_end'] = major['peak']; r[P + 'distal_end'] = ud; r[P + 'major_share'] = major['reads'] / rd
        alt = TC.alt_ends(cl, rd); r[P + 'n_alt'] = len(alt)
        if side == 5:
            cg = [cage_at(g.chrom, g.strand, (g.cds5 - x['peak']) if g.strand == '+' else (g.cds5 + x['peak'])) for x in alt]
            r['p5_n_alt_cage'] = int(sum(v >= 1 for v in cg)); r['p5_alt_positions'] = ','.join(str(x['peak']) for x in alt)
        if L == 0:
            r[P + 'cat'] = 'No annotated UTR; reads end at CDS (≤50 bp)' if major['peak'] <= 50 else 'No annotated UTR; long reads define a UTR'
            continue
        tol = max(200, 0.25 * L)
        r[P + 'cat'] = 'Long reads extend beyond annotated end' if ud > L + 200 else ('Annotated end supported' if ud >= L - tol else 'Annotated UTR longer than any supported end')
        r[P + 'excess_bp'] = max(0, L - ud)
        U = indep_units(g.chrom, g.strand, a, b, g.gene); okU, ndU, rdU = TC.unit_support(U, G)
        AS = antisense_units(g.chrom, g.strand, a, b, g); okA = bool(len(AS) and AS.dataset.nunique() >= 2)
        r[P + 'sense_unit'] = okU; r[P + 'sense_unit_ds'] = ndU; r[P + 'sense_unit_reads'] = rdU; r[P + 'antisense_unit'] = okA
        if okU:
            far = [((g.cds5 - x.start) if g.strand == '+' else (x.end - 1 - g.cds5)) if side == 5 else ((x.end - 1 - g.cds3) if g.strand == '+' else (g.cds3 - x.start)) for x in U.itertuples()]
            near = [((g.cds5 - x.end + 1) if g.strand == '+' else (x.start - g.cds5)) if side == 5 else ((x.start - g.cds3) if g.strand == '+' else (g.cds3 - x.end + 1)) for x in U.itertuples()]
            r[P + 'unit_far'] = int(np.median(far)); r[P + 'unit_near'] = int(np.median(near))
    rows.append(r)
C = pd.DataFrame(rows)
