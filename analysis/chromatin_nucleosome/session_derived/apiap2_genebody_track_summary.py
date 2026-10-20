"""Session-derived computation of apiap2_genebody_track_summary.tsv.

Recovered from the record of the interactive analysis session in which apiap2_genebody_track_summary.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: apiap2_tss_profiles.tsv, apiap2_tss_aggregate_profiles_long.tsv.gz, H6.tsv, ref_tracks.tsv, S.tsv, gc_curves.tsv, controls.tsv, fc_bundle.tar
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
# Verbatim code of the producing analysis cell.

g29=S[S.track_class=='genebody_or_TSS_depleted']
print('T1',g29.assigned_causes.str.contains('T1').sum(),'T2',g29.assigned_causes.str.contains('T2').sum(),'B1',g29.assigned_causes.str.contains('B1').sum(),'B2',g29.assigned_causes.str.contains('B2').sum(),'B3',g29.assigned_causes.str.contains('B3').sum(),'S1',g29.assigned_causes.str.contains('S1').sum())
print('GB expr|GC,HP1 median/max',g29.body_expr_given_GC_HP1.median().round(3),g29.body_expr_given_GC_HP1.max().round(3),'GC|expr median',g29.body_GC_given_expr_HP1.median().round(3))
print('fPbelow GB range',g29.frac_promoter_bins_IP_below_input.min().round(2),g29.frac_promoter_bins_IP_below_input.max().round(2),'median',g29.frac_promoter_bins_IP_below_input.median().round(2),'fG',g29.frac_genic_bins_IP_below_input.min().round(2),g29.frac_genic_bins_IP_below_input.max().round(2))
print('GB GCcorr genic range',g29.GCcorr_genic.min().round(3),g29.GCcorr_genic.max().round(3),'P-B range',g29.GCcorr_P_minus_B.min().round(2),g29.GCcorr_P_minus_B.max().round(2))
hf=['AP2-HC','AP2-G5','AP2-O5','AP2-12','AP2-G2','AP2-11A','AP2-exp']
allB1=S.GCcorr_HP1_minus_rest>max(0.25,GFP_HP1)
print('tracks with HP1 excess (all 39):',S.index[allB1].tolist())
print(pd.crosstab(S.factor.isin(hf),allB1), fe(pd.crosstab(S.factor.isin(hf),allB1).values))
print('IP GC slope GB median',g29.IP_GC_top_vs_bottom_decile_log2.median().round(2),'input',g29.input_GC_top_vs_bottom_decile_log2.median().round(2),'promoter IP',S.loc[S.track_class!='genebody_or_TSS_depleted','IP_GC_top_vs_bottom_decile_log2'].median().round(2))
print('wilcoxon IP vs input GC slope GB', __import__('scipy.stats',fromlist=['wilcoxon']).wilcoxon(g29.IP_GC_top_vs_bottom_decile_log2,g29.input_GC_top_vs_bottom_decile_log2))
print('GB tracks GC rho intergenic median',g29.bins_rho_GC_intergenic.median().round(2))
print(ve[:2], spearmanr(pc1,GCbody)[0])
S.to_csv('apiap2_genebody_track_summary.tsv',sep='\t',float_format='%.4g')
