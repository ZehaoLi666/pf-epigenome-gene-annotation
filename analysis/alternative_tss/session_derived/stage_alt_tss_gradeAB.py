"""Session-derived computation of stage_alt_tss_gradeAB.tsv.

Recovered from the record of the interactive analysis session in which stage_alt_tss_gradeAB.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: alt_tss_classified.tsv, gene_level_isoform_summary.tsv, daft_pairs_counts.tsv.gz, cage_pair_counts.tsv, pairs_positions.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
# Verbatim code of the producing analysis cell.

keep=['pair_id','is_control','gene_id','desc','chrom','strand','atg','main_TSS','main_src','lr_group','alt_pos','sources','n_merged','class','alt_rel_main','d_alt_vs_main','d_alt_vs_ATG','w',
 'alt_within100_model_call','contamination_flag','upstream_sameStrand_3end_genes']
for nm,labs in [('u5',hrs),('cage',['ring','troph','schizont']),('lr',['ring','troph','schizont'])]:
    keep+= [f'{nm}_n_tp_tested',f'{nm}_G',f'{nm}_phi',f'{nm}_p_quasi',f'{nm}_q',f'{nm}_dPSI',f'{nm}_peak',f'{nm}_trough',f'{nm}_sig']+[f'{nm}_PSI_{l}' for l in labs]+[f'{nm}_peak_stage']
keep+=[f'u5_alt_{t}' for t in TPs]+[f'u5_main_{t}' for t in TPs]+[f'cage_alt_{s}' for s in ['ring','trophozoite','schizont']]+[f'cage_main_{s}' for s in ['ring','trophozoite','schizont']]
keep+=[f'lrS_alt_{s}' for s in ST]+[f'lrS_main_{s}' for s in ST]
keep+=['daft_U_total','daft_S_total','daft_len_U','daft_len_S']+[f'daft_ratio_{h}h' for h in hours]+['daft_G','daft_df','daft_phi','daft_p_quasi','daft_q','daft_log2FC','daft_lane_r','daft_up_peak','daft_alt_peak','daft_sig']
keep+=[c_ for c_ in P2.columns if c_.startswith('conc_')]+['grade','grade_evidence']
OUT=P2[keep].copy()
OUT.rename(columns={f'u5_alt_{t}':f'u5_alt_pseudocount_{h}' for t,h in zip(TPs,hrs)}|{f'u5_main_{t}':f'u5_main_pseudocount_{h}' for t,h in zip(TPs,hrs)},inplace=True)
OUT.to_csv('stage_alt_tss_all_candidates.tsv',sep='\t',index=False)
abcols=['grade','grade_evidence','gene_id','desc','families','chrom','strand','atg','main_TSS','main_src','alt_pos','sources','class','alt_rel_main','d_alt_vs_main','d_alt_vs_ATG',
 'consensus_peak_stage','u5_peak','u5_dPSI','u5_q','cage_peak','cage_dPSI','cage_q','lr_peak','lr_dPSI','lr_q','daft_alt_peak','daft_log2FC','daft_q','evidence_detail','expr_peak',
 'alt_within100_model_call','contamination_flag','upstream_sameStrand_3end_genes',
 'alt_utr5_len','main_utr5_len','delta_utr5_len','alt_n_uAUG','main_n_uAUG','delta_n_uAUG','alt_n_uORF_contained','main_n_uORF_contained','alt_n_oORF_overlap_CDS','main_n_oORF_overlap_CDS','alt_inframe_Nterm_ext','main_inframe_Nterm_ext']
ABo=AB[abcols].sort_values(['grade','gene_id']); ABo['desc']=ABo.desc.str.replace('%2C',',')
ABo.to_csv('stage_alt_tss_gradeAB.tsv',sep='\t',index=False)
ENR['name']=ENR.name.astype(str); ENR.to_csv('stage_alt_tss_enrichment.tsv',sep='\t',index=False)
CMP.to_csv('stage_alt_tss_vs_prior33.tsv',sep='\t',index=False)
print(OUT.shape, ABo.shape, ABo.gene_id.nunique(), ABo.alt_inframe_Nterm_ext.sum(), ABo.main_inframe_Nterm_ext.sum())
