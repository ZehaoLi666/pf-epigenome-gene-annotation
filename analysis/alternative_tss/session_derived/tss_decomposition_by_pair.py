"""Session-derived computation of tss_decomposition_by_pair.tsv.

Recovered from the record of the interactive analysis session in which tss_decomposition_by_pair.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: tss_stage_usage.tsv, alt_tss_classified.tsv, gene_level_isoform_summary.tsv, report_uorf_translation.md, ribo_region_counts2.tsv, ribo_region_counts.tsv.gz, cage_rep_window_counts.tsv, daft_gene_counts.tsv, cds6_cage_profile.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
# Verbatim code of the producing analysis cell.

hst={'0h':'ring','8h':'ring','16h':'troph','24h':'troph','32h':'schizont','40h':'schizont'}
cm=CPM[[c_ for c_ in CPM.columns if c_.split('_')[1] in hst]]
stg_cpm=cm.T.groupby(cm.columns.map(lambda c_: hst[c_.split('_')[1]])).mean().T
stg_cpm.columns=[f'DAFT_CPM_{s}' for s in stg_cpm.columns]
OUT=ABT.merge(F,on='pair_id',how='left').merge(stg_cpm,left_on='gene_id',right_index=True,how='left')
extra=[c_ for c_ in sel.columns if c_ not in OUT.columns or c_=='pair_id']
OUT=OUT.merge(sel[extra],on='pair_id',how='left')
OUT['contamination_flag']=OUT.gene_id.map(contam).fillna('')
OUT['case_label_cn']=OUT.case.map({'1_expression_only':'1 仅表达变化','2_TSS_switch_only':'2 仅TSS比例切换','3_both':'3 两者皆变','4_neither':'4 均不变','5_unclassified':'5 不可分类'})
lead=['pair_id','gene_id','desc','grade','case','case_label_cn','contamination_flag','chrom','strand','atg','main_TSS','alt_pos','distal_pos','proximal_pos','distal_is_alt','prop_status','prop_direction_consistent']
OUT=OUT[lead+[c_ for c_ in OUT.columns if c_ not in lead and c_!='dominant_at_hi_stage']]
OUT=OUT.round(5); OUT.to_csv('tss_decomposition_by_pair.tsv',sep='\t',index=False)
ctl_cols=['pair_id','gene_id','desc','distal_pos','proximal_pos','CAGE_n_total','CAGE_dprop','CAGE_q','U5_n_total','U5_dprop','U5_q','prop_status','DAFT_log2FC','DAFT_q','case']
CTT[[c_ for c_ in ctl_cols if c_ in CTT.columns]].round(5).to_csv('tss_decomposition_controls.tsv',sep='\t',index=False)
c6o=C6.copy(); c6o['desc']=c6o.desc.str.replace('%2C',',')
c6o=c6o.merge(pv.reset_index().assign(gene_id=lambda d:d.pair_id.str.replace('CDSTSS_',''))[['gene_id','Nterm_before_altATG','after_altATG','Nterm_over_downstream']].rename(columns={'Nterm_before_altATG':'ribo_Psite_density_Nterm_before_altATG','after_altATG':'ribo_Psite_density_after_altATG','Nterm_over_downstream':'ribo_Nterm_over_downstream_ratio'}),on='gene_id',how='left')
c6o.round(4).to_csv('cds_internal_tss_calls.tsv',sep='\t',index=False)
print(OUT.shape, len(CTT), c6o.shape)
