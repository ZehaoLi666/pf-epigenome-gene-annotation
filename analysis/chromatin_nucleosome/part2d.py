import sys, os; sys.path.append(os.getcwd())
exec(open('part2c.py').read().split("rows = []")[0])
bt = pd.read_csv(TRUTH, sep='\t')[['gene_id', 'Consensus_5end__TSS', 'Adjalley2016_capseq__TSS']]
b2 = base.merge(bt, on='gene_id', how='left')
base_all = base
res = []
for gsub, mask in [('Consensus_5end genes', b2.Consensus_5end__TSS.notna().values), ('v71 UTR >300bp longer than refined', ((b2.utr5_annot - (b2.ref_TSS - b2.atg).abs()) > 300).values), ('v71 within 100bp of refined', ((b2.v71_TSS - b2.ref_TSS).abs() <= 100).values)]:
    for anc, col in [('Consensus_5end', 'Consensus_5end__TSS'), ('refined', 'ref_TSS'), ('PlasmoDB_v71', 'v71_TSS'), ('ATG', 'atg')]:
        m = mask & b2[col].notna().values
        if gsub != 'Consensus_5end genes' and anc == 'Consensus_5end': continue
        if gsub == 'Consensus_5end genes': m = mask & b2.Consensus_5end__TSS.notna().values
        X = np.where(np.isnan(b2[col].values.astype(float)), 0, b2[col].values.astype(float)).astype(int)
        for (tp, kind) in [('T10', 'all_calls'), ('T30', 'all_calls'), ('T40', 'all_calls'), ('reference_map', 'all_calls'), ('reference_map', 'top20_well_positioned')]:
            cs = sets[(tp, kind)]; bgd = sum(len(v) for c, v in cs.items() if c in CHROMS) / gsize
            F = flat(offsets(cs, X)); wt = m.astype(float)
            L = lm(density(F, wt), bgd); L['frac_genes_dyad_within20_of_plus1_mode'] = w20(F, wt, L['plus1_mode'])
            res.append(dict(gene_subset=gsub, anchor=anc, timepoint=tp, call_set=kind, n_genes=int(m.sum()), expected_w20_random=1 - np.exp(-41 * bgd), **L))
R = pd.DataFrame(res); R.to_csv('nucleosome_dyad_density_calibration.tsv', sep='\t', index=False, float_format='%.4g')
print(R[R.call_set=='all_calls'][['gene_subset','anchor','timepoint','n_genes','plus1_mode','plus1_height_over_bg','frac_genes_dyad_within20_of_plus1_mode','NDR_trough','NDR_trough_over_bg','plus1_to_plus2']].round(3).to_string())
