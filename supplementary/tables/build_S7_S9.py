from pfconfig import FIGURES_DIR, result  # paths: config/paths.yaml
import pandas as pd, numpy as np, re, urllib.parse
import os
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'common.py')).read())
FIG = FIGURES_DIR
rd = lambda v: pd.read_csv(result(v), sep='\t', compression='gzip' if result(v).endswith('.gz') else 'infer', low_memory=False)
STG = {'ring': 'ring', 'troph': 'trophozoite', 'schizont': 'schizont', 'trophozoite': 'trophozoite'}
def unq(d):
    for c in ('desc',):
        if c in d: d[c] = d[c].map(lambda s: urllib.parse.unquote(s) if isinstance(s, str) else s)
    return d

# ======================= S7 =======================
s7 = []
V1 = 'model_tss_calls_validation.tsv'; C = rd(V1)
dC = {'gene_id': 'PlasmoDB gene ID.', 'chrom': 'Chromosome.', 'strand': 'Strand.', 'pos': 'Model TSS call (0-based, 20-bp grid).', 'role': 'main = highest-probability call (refined/model TSS); alt = other call in the gene window (alternative TSS).',
      'atg': 'ATG (0-based).', 'stop': 'Last base of stop codon (0-based).', 'rel_atg': 'Strand-aware position relative to the ATG (negative = upstream).', 'cls': 'upstream (of ATG) / CDS / downstream.',
      'cage_ring': 'Adjalley 2016 CAGE CPM summed within ±50 bp of the call, ring.', 'cage_trophozoite': 'Same, trophozoite.', 'cage_schizont': 'Same, schizont.', 'cage_flank_perbp': 'Mean CAGE CPM per bp in flanks 100–500 bp on both sides (all stages).',
      'cage_win': 'Sum of the three stage CAGE windows.', 'cage_enrich': '(cage_win / 101 bp) / (cage_flank_perbp + 0.01).', 'u5_hit': 'A Chappell 2020 5′UTR-seq tag-cluster dominant position lies within 50 bp (same strand).', 'u5_score': 'Highest 5′UTR-seq cluster score among hits.',
      'cap5_hit': 'A Shaw 2021 5CAPture cluster dominant position lies within 50 bp.', 'cap5_score': 'Highest 5CAPture cluster score among hits.', 'cage_hit': 'CAGE support: cage_win ≥5 CPM and cage_enrich ≥2.',
      'n_sup': 'Number of independent 5′-end datasets supporting the call (CAGE, 5′UTR-seq, 5CAPture; 0–3).', 'sup1': 'n_sup ≥1.', 'sup2': 'n_sup ≥2 ("validated").'}
s7.append(dict(name='Model_TSS_calls_validation', df=C, desc=dC, source=f'model_tss_calls_validation.tsv ({V1})', note=f'{len(C)} calls (main {int((C.role=="main").sum())}, alt {int((C.role=="alt").sum())})'))
V2 = 'model_validated_alt_tss_per_gene.tsv'; PG = unq(rd(V2))
assert (PG.n_alt_sup2 >= 1).sum() == 609 and len(PG) == 4181
dP = {'gene_id': 'PlasmoDB gene ID (genes with a main model TSS call).', 'desc': 'Product description.', 'lr_group': 'Long-read status.', 'n_alt_calls': 'Alternative TSS calls in the gene window.', 'n_alt_sup1': 'Alternative calls supported by ≥1 independent 5′ dataset.',
      'n_alt_sup2': 'Alternative calls supported by ≥2 datasets (validated alternative TSSs; 609 genes have ≥1).', 'main_sup2': 'Main TSS call supported by ≥2 datasets.', 'main_cage': 'CAGE CPM (cage_win) at the main call.', 'DAFT_max': 'Maximum DAFT-seq CPM across time points.',
      'up_ig': 'Upstream intergenic length (bp) to the neighbouring gene.', 'cds_len': 'CDS length (as in the gene table; unit to verify).'}
s7.append(dict(name='Validated_altTSS_per_gene', df=PG, desc=dP, source=f'model_validated_alt_tss_per_gene.tsv ({V2})', note='609/4,181 genes (14.6%) with ≥1 alternative TSS supported by ≥2 datasets'))
V3 = 'model_alt_tss_family_enrichment.tsv'; FE = rd(V3); V4 = 'model_alt_tss_go_enrichment.tsv'; GO = rd(V4)
dE = {'family': 'Gene family / class.', 'n_genes': 'Genes in the family (with a main call).', 'n_alt_total': 'Validated alternative TSSs in the family.', 'RR_adj': 'Rate ratio of validated alternative TSSs per gene vs other genes, adjusted for expression and window length (model in source analysis).',
      'lo': '95% CI lower.', 'hi': '95% CI upper.', 'p': 'p value.', 'q': 'BH q value.', 'frac_ge1': 'Fraction of family genes with ≥1 validated alternative TSS.', 'tier': 'Support tier used (≥1 or ≥2 datasets).',
      'go_id': 'GO term.', 'n': 'Genes annotated with the term.', 'name': 'GO term name.', 'ns': 'GO namespace.'}
s7.append(dict(name='altTSS_family_enrichment', df=FE, desc=dE, source=f'model_alt_tss_family_enrichment.tsv ({V3})', note='Fig S9'))
s7.append(dict(name='altTSS_GO_enrichment', df=GO, desc=dE, source=f'model_alt_tss_go_enrichment.tsv ({V4})', note='Fig S9'))
V5 = 'stage_alt_tss_gradeAB.tsv'; GA = unq(rd(V5))
assert GA.grade.value_counts().to_dict() == {'B': 41, 'A': 24} and GA.gene_id.nunique() == 57
dG = {'grade': 'A = 5′UTR-seq + CAGE agree (or CAGE + long reads); B = CAGE + DAFT-seq, 5′UTR-seq + DAFT-seq or long reads + DAFT-seq.', 'grade_evidence': 'Sources that are significant and agree in stage direction.',
      'gene_id': 'PlasmoDB gene ID.', 'desc': 'Product description.', 'families': 'Family / feature flags.', 'chrom': 'Chromosome.', 'strand': 'Strand.', 'atg': 'ATG (0-based).', 'main_TSS': 'Main TSS (0-based; long-read TSS for reference genes, else model TSS).',
      'main_src': 'Source of main_TSS.', 'alt_pos': 'Alternative TSS (0-based).', 'sources': 'Sources that proposed the alternative TSS candidate (long-read cluster / model alt call / ...).', 'class': 'Position class of the alternative TSS (upstream_of_ATG / in CDS).',
      'alt_rel_main': 'Alternative upstream or downstream of the main TSS.', 'd_alt_vs_main': 'Strand-aware distance alt − main (bp; negative = upstream).', 'd_alt_vs_ATG': 'Strand-aware distance alt − ATG (bp).', 'consensus_peak_stage': 'Stage of peak alternative-TSS usage agreed by the significant sources.',
      'u5_peak': '5′UTR-seq time point with the highest alternative PSI.', 'u5_dPSI': '5′UTR-seq ΔPSI (max − min PSI of the alternative site across time points).', 'u5_q': '5′UTR-seq G-test BH q.', 'cage_peak': 'CAGE stage with highest alternative PSI.', 'cage_dPSI': 'CAGE ΔPSI.', 'cage_q': 'CAGE quasi-likelihood test BH q.',
      'lr_peak': 'Long-read stage with highest PSI.', 'lr_dPSI': 'Long-read ΔPSI.', 'lr_q': 'Long-read G-test BH q.', 'daft_alt_peak': 'DAFT-seq time point where the unique/shared coverage ratio peaks.', 'daft_log2FC': 'DAFT-seq log2 fold change of the ratio.', 'daft_q': 'DAFT-seq BH q.',
      'evidence_detail': 'Per-source summary string.', 'expr_peak': 'DAFT-seq time point of peak gene expression.', 'alt_within100_model_call': 'Alternative TSS within 100 bp of a model TSS call.', 'contamination_flag': 'Upstream same-strand gene 3′ end overlaps the window (possible read-through).',
      'upstream_sameStrand_3end_genes': 'Upstream same-strand genes whose 3′ ends fall in the window.', 'alt_utr5_len': '5′UTR length from the alternative TSS (bp).', 'main_utr5_len': '5′UTR length from the main TSS.', 'delta_utr5_len': 'alt − main 5′UTR length.',
      'alt_n_uAUG': 'uAUGs in the alternative 5′UTR (genomic sequence; unspliced).', 'main_n_uAUG': 'uAUGs in the main 5′UTR.', 'delta_n_uAUG': 'Difference.', 'alt_n_uORF_contained': 'uORFs ending before the ATG (alternative UTR).', 'main_n_uORF_contained': 'Same, main UTR.',
      'alt_n_oORF_overlap_CDS': 'Out-of-frame ORFs overlapping the CDS start (alternative UTR).', 'main_n_oORF_overlap_CDS': 'Same, main UTR.', 'alt_inframe_Nterm_ext': 'In-frame upstream ATG without stop (potential N-terminal extension), alternative UTR.', 'main_inframe_Nterm_ext': 'Same, main UTR.'}
s7.append(dict(name='Stage_specific_GradeAB', df=GA, desc=dG, source=f'stage_alt_tss_gradeAB.tsv ({V5})', note='65 pairs / 57 genes (A 24, B 41)'))
V6 = 'tss_decomposition_by_pair.tsv'; DC = unq(rd(V6))
assert DC.case.value_counts().to_dict() == {'5_unclassified': 33, '2_TSS_switch_only': 16, '3_both': 9, '4_neither': 6, '1_expression_only': 1}
def ddesc(c):
    fixed = {'pair_id': 'gene_id|alternative TSS position.', 'gene_id': 'PlasmoDB gene ID.', 'desc': 'Product description.', 'grade': 'Grade from Stage_specific_GradeAB.', 'case': 'Decomposition class: 1 expression only; 2 TSS-usage switch only; 3 both; 4 neither; 5 unclassified.',
             'case_label_cn': 'Class label (Chinese).', 'contamination_flag': 'flag (upstream gene overlap) or exclude (MSP4 read-through; MSP5 pair excluded from case 3 in the manuscript).', 'chrom': 'Chromosome.', 'strand': 'Strand.', 'atg': 'ATG (0-based).',
             'main_TSS': 'Main TSS (0-based).', 'alt_pos': 'Alternative TSS (0-based).', 'distal_pos': 'More upstream TSS of the pair (0-based).', 'proximal_pos': 'More downstream TSS (0-based).', 'distal_is_alt': 'True if the distal TSS is the alternative one.',
             'prop_status': 'Distal-share status: changed / unchanged / ambiguous (both CAGE and 5′UTR-seq must change, q <0.1, Δ ≥0.2, same direction).', 'prop_direction_consistent': 'CAGE and 5′UTR-seq change in the same direction.', 'd_alt_vs_main': 'Strand-aware distance alt − main (bp).',
             'is_control': 'Single-TSS pseudo-pair control.', 'set': 'gradeAB or control.', 'DAFT_p_stage': 'NB-GLM likelihood-ratio p for stage effect on gene expression (DAFT-seq).', 'DAFT_log2FC': 'DAFT-seq log2 fold change for the contrast in DAFT_contrast.', 'DAFT_log2FC_lo95': 'Wald 95% CI lower.',
             'DAFT_log2FC_hi95': 'Wald 95% CI upper.', 'DAFT_contrast': 'Stage contrast (highest vs lowest distal-share stage if changed, else max − min).', 'DAFT_q': 'BH q.', 'expr_changed': 'DAFT q <0.1 and |log2FC| ≥1.', 'expr_unchanged': '|log2FC| <1 or q ≥0.1.',
             'min_expr_change_contrast_with_switch': 'Post hoc: stage pair with both Δ ≥0.2 and the smallest |DAFT log2FC|.', 'DAFT_CPM_ring': 'DAFT-seq CPM, ring.', 'DAFT_CPM_troph': 'DAFT-seq CPM, trophozoite.', 'DAFT_CPM_schizont': 'DAFT-seq CPM, schizont.',
             'dominant_iso_distal_high_stage': 'Dominant isoform (distal/proximal) at the high-distal stage.', 'dominant_iso_distal_low_stage': 'Dominant isoform at the low-distal stage.', 'dominant_isoform_changes': 'Dominant isoform differs between the two stages.',
             'TE_hi': 'Polysome TE (log2, Bunnik 2013) at the high-distal stage.', 'TE_lo': 'TE at the low-distal stage.', 'dTE': 'TE_hi − TE_lo.', 'dTE_resid': 'dTE residual after covariates (CDS length, GC, expression).', 'dTE_z': 'dTE z-score.',
             'TE_pass_both': 'TE available at both stages.', 'TE_contrast': 'Stages used for TE (0 h ring TE considered unreliable).', 'dTE_percentile_vs_singleTSS_controls': 'Percentile of dTE among single-TSS controls.', 'dTE_pct_expr_matched': 'Percentile of dTE among expression-matched background genes (manuscript: median 0.11, n = 16).',
             'n_matched_controls': 'Number of matched controls.', 'n_matched_bg': 'Number of expression-matched background genes.', 'ribo_distalUTR_log2ratio_hi_vs_lo': 'Caro 2014 Ribo-seq: log2 ratio of distal-unique leader footprint density (relative to CDS), high- vs low-distal stage.'}
    if c in fixed: return fixed[c]
    m = re.match(r'(CAGE|U5)_(.*)', c)
    if m:
        src = 'CAGE (Adjalley 2016)' if m.group(1) == 'CAGE' else '5′UTR-seq (Chappell 2020)'; k = m.group(2)
        d = {'n_total': 'total tag count in the two windows', 'dprop': 'Δ distal share (highest − lowest stage)', 'prop_hi_stage': 'stage with highest distal share', 'prop_lo_stage': 'stage with lowest distal share',
             'bb_p_stage_adjDAFT': 'β-binomial LRT p for stage, adjusted for log2 DAFT CPM', 'bb_p_stage_unadj': 'β-binomial LRT p, unadjusted', 'bb_rho': 'β-binomial overdispersion', 'OR_hi_vs_lo_adj': 'adjusted odds ratio hi vs lo stage (unstable; reference only)',
             'OR_lo95': 'OR 95% CI lower', 'OR_hi95': 'OR 95% CI upper', 'qb_p_stage_adjDAFT': 'quasi-binomial F-test p (robustness)', 'qb_scale': 'quasi-binomial Pearson scale', 'q': 'BH q of the β-binomial p', 'qb_q': 'BH q of the quasi-binomial p',
             'changed': 'q <0.1 and Δ ≥0.2', 'unchanged': 'Δ <0.15 or q ≥0.1', 'cage_gene_total_by_stage': 'CAGE gene total per stage'}
        mm = re.match(r'(prop|n)_(ring|troph|schizont)$', k)
        if mm: return f'{src}: {"distal share distal/(distal+proximal)" if mm.group(1) == "prop" else "tag count in both windows"} in {STG[mm.group(2)]}.'
        return f'{src}: {d.get(k, k)}.'
    m = re.match(r'(distal|proximal)_(.*)', c)
    if m:
        d = {'utr5_len': '5′UTR length (bp)', 'n_uAUG': 'uAUG count', 'n_uORF_contained': 'uORFs ending before the ATG', 'n_uORF_overlap_CDS': 'uORFs overlapping the CDS', 'n_Nterm_ext': 'in-frame N-terminal extension ATGs', 'n_uAUG_kozak_m3R_p4G': 'uAUGs with R at −3 and G at +4',
             'ATG_m3': 'base at −3 of the main ATG', 'ATG_p4': 'base at +4', 'ATG_kozak_m3R': 'main ATG has R at −3', 'ATG_kozak_p4G': 'main ATG has G at +4', 'utr5_AT': 'AT fraction of the 5′UTR'}
        return f'{m.group(1)} isoform: {d.get(m.group(2), m.group(2))}.'
    m = re.match(r'ribo_(.*)_(ring|troph|schizont)$', c)
    if m:
        d = {'distal_unique_psites': 'P-sites in the distal-unique leader', 'distal_unique_density_rel_CDS': 'distal-unique leader footprint density relative to CDS', 'uORF_distalOnly_psites': 'P-sites in distal-only uORFs',
             'uORF_distalOnly_density_rel_CDS': 'distal-only uORF density relative to CDS', 'proximal_utr_psites': 'P-sites in the proximal 5′UTR', 'proximal_utr_density_rel_CDS': 'proximal 5′UTR density relative to CDS', 'CDS_psites': 'CDS P-sites'}
        return f'Caro 2014 Ribo-seq (W2), {STG[m.group(2)]}: {d.get(m.group(1), m.group(1))}.'
    return None
s7.append(dict(name='Decomposition_by_pair', df=DC, desc={c: ddesc(c) for c in DC.columns if ddesc(c)}, source=f'tss_decomposition_by_pair.tsv ({V6})', note='65 pairs: case 2 = 16, case 3 = 9 (8 after excluding MSP5), case 1 = 1, case 4 = 6, case 5 = 33'))
V7 = 'tss_decomposition_controls.tsv'; CT = unq(rd(V7))
assert len(CT) == 888 and (CT.prop_status == 'changed').sum() == 0
s7.append(dict(name='Decomposition_controls', df=CT, desc={c: ddesc(c) for c in CT.columns if ddesc(c)}, source=f'tss_decomposition_controls.tsv ({V7})', note='888 single-TSS pseudo-pairs; 0 classed as distal-share change'))
F5J = pd.read_csv(f'{FIG}/Figure5/data/Figure5_J_uAUG_pairs.tsv', sep='\t'); F5I = pd.read_csv(f'{FIG}/Figure5/data/Figure5_I_focal_genes.tsv', sep='\t')
s7.append(dict(name='Fig5I_focal_genes', df=F5I, desc={**{c: ddesc(c) for c in F5I.columns if ddesc(c)}, 'gene_name': 'Gene name.'}, source='Figure5/data/Figure5_I_focal_genes.tsv', note='PfA-M1, LDH, eIF3i'))
s7.append(dict(name='Fig5J_uAUG_pairs', df=F5J, desc={c: ddesc(c) for c in F5J.columns if ddesc(c)}, source='Figure5/data/Figure5_J_uAUG_pairs.tsv', note='24 switch pairs (cases 2 + 3, MSP5 excluded)'))
V8 = 'uorf_census.tsv'; U = rd(V8)
def udesc(c):
    fixed = {'gene_id': 'PlasmoDB gene ID.', 'chrom': 'Chromosome.', 'strand': 'Strand.', 'atg': 'ATG (0-based).', 'lr_group': 'Long-read status.', 'final_src': 'Refined TSS source.', 'final_TSS': 'Refined TSS (0-based).', 'annot_tss': 'v71 TSS (0-based).',
             'utr5_change_vs_v71': '5′ change class vs v71 (Table S4).', 'refined_utr5_source': 'Source of the refined 5′UTR in the refined annotation (long-read / model / v71 kept).', 'n_uAUG_shared': 'uAUGs present in both refined and v71 UTRs.',
             'n_uAUG_gained': 'uAUGs only in the refined UTR.', 'n_uAUG_lost': 'uAUGs only in the v71 UTR.', 'n_uAUG_lost_trimmed': 'Lost uAUGs upstream of the refined TSS (removed by trimming).', 'n_uORF_shared': 'uORFs in both.', 'n_uORF_gained': 'uORFs only refined.',
             'n_uORF_lost': 'uORFs only v71.', 'n_uORF_lost_trimmed': 'Lost uORFs upstream of the refined TSS.', 'n_oORF_gained': 'Overlapping ORFs gained.', 'n_oORF_lost': 'Overlapping ORFs lost.', 'ref_uORF_status_vs_v71': 'uORF_both / uORF_gained / uORF_lost / uORF_neither.'}
    if c in fixed: return fixed[c]
    m = re.match(r'(ref|v71)_(.*)', c)
    if not m: return None
    w = 'refined 5′UTR' if m.group(1) == 'ref' else 'v71 5′UTR'
    d = {'utr5_len': 'length (nt, spliced)', 'utr5_AT': 'AT fraction', 'n_uAUG': 'uAUG count', 'n_uORF': 'uORFs (≥2 codons)', 'n_uORF_contained': 'uORFs ending before the ATG', 'n_oORF_overlapping': 'uORFs crossing the ATG out of frame',
         'n_startstop': 'start-stop uORFs (AUG immediately followed by stop)', 'n_Nterm_ext_inframe': 'in-frame uAUGs without stop (N-terminal extension)', 'max_uORF_codons': 'longest uORF (codons)', 'median_uORF_codons': 'median uORF length (codons)',
         'n_uAUG_m3A': 'uAUGs with A at −3', 'n_uAUG_m3R_p4G': 'uAUGs with R at −3 and G at +4', 'max_uAUG_kozak_pwm': 'best uAUG Kozak PWM score', 'n_uORF_strongctx': 'uORFs with strong context (PWM ≥ 0.387)', 'dist_uAUG_nearest_ATG': 'distance of the nearest uAUG to the ATG (nt)', 'oORF_any': 'any overlapping ORF'}
    return f'{w}: {d.get(m.group(2), m.group(2))}.'
s7.append(dict(name='uORF_census', df=U, desc={c: udesc(c) for c in U.columns if udesc(c)}, source=f'uorf_census.tsv ({V8})', note='5,285 genes; refined vs v71 5′UTR'))
V9 = 'translation_efficiency.tsv'; TE = rd(V9)
def tdesc(c):
    fixed = {'gene_id': 'PlasmoDB gene ID.', 'lr_group': 'Long-read status.', 'final_src': 'Refined TSS source.', 'cds_len': 'CDS length (nt).', 'cds_GC': 'CDS GC fraction.', 'upstream_intergenic': 'Upstream intergenic length (bp).', 'DAFT_max': 'Maximum DAFT-seq CPM.',
             'Qd': 'DAFT-seq expression quintile.', 'stage_quality_note': 'Library-quality note (0 h polysome library is low quality; 0 h TE unreliable).'}
    if c in fixed: return fixed[c]
    m = re.match(r'(.*)_(0h|18h|36h|36h_r1|36h_r2)$', c)
    if not m: return None
    d = {'poly_reads': 'polysome CDS fragment count', 'ss_reads': 'steady-state mRNA CDS fragment count', 'poly_CPM': 'polysome CDS CPM', 'ss_CPM': 'steady-state CDS CPM', 'log2TE_raw': 'log2[(polysome CPM + 0.5 read)/(steady-state CPM + 0.5 read)]',
         'log2TE': 'log2 TE centred on the stage median', 'TE_pass': 'steady-state ≥20 fragments', 'cds_pairs_poly': 'polysome CDS read pairs', 'cds_pairs_ss': 'steady-state CDS read pairs', 'proxy_log2TE_polyBW_over_DAFT': 'proxy TE: polysome bigWig CDS mean / DAFT (sensitivity analysis)',
         'DAFT_nearest': 'DAFT CPM at the nearest time point'}
    return f'{d.get(m.group(1), m.group(1))}, {m.group(2)} (Bunnik 2013).'
s7.append(dict(name='Translation_efficiency', df=TE, desc={c: tdesc(c) for c in TE.columns if tdesc(c)}, source=f'translation_efficiency.tsv ({V9})', note='Bunnik 2013 polysome vs matched steady-state mRNA'))
V10 = 'altTSS_distal_share_vs_TE.tsv'; DS = rd(V10)
s7.append(dict(name='Distal_share_vs_TE', df=DS, desc={'tp': 'Time point.', 'n': 'Genes.', 'rho_distal_TE': 'Spearman ρ between distal-isoform share and TE.', 'p_rho': 'p.', 'beta_distal_adj': 'Adjusted regression coefficient of distal share on TE.', 'se': 'SE.', 'p_beta': 'p.',
    'beta_wlogUTR': 'Coefficient with log UTR length weighting (see source analysis).', 'p_wlog': 'p.', 'R2_weighted': 'Model R² (weighted).', 'R2_mainonly': 'Model R² with main isoform only.'}, source=f'altTSS_distal_share_vs_TE.tsv ({V10})', note='genome-wide; column definitions partly to verify'))
intro7 = ['Alternative TSSs: (i) model-called alternative TSSs and their validation by independent 5′-end data (609 genes, 14.6%, with ≥1 alternative TSS supported by ≥2 datasets; Fig S9); (ii) stage-specific Grade A/B TSS pairs (65 pairs, 57 genes); (iii) decomposition of each pair into TSS-usage switch vs expression change, with 888 single-TSS controls (Figure 5H–J); (iv) uAUG/uORF census and translation efficiency.',
          'Decomposition criteria: distal share = distal/(distal + proximal); changed = CAGE and 5′UTR-seq both q <0.1, Δ ≥0.2 and same direction; expression changed = DAFT-seq q <0.1 and |log2FC| ≥1.',
          'Not included: uaug_catalog.tsv.gz (per-uAUG catalogue, uaug_catalog.tsv.gz) and te_association_models.tsv (te_association_models.tsv); not cited in the manuscript, available on request.']
p7, m7 = write_table('TableS7_alternative_TSS.xlsx', 'Table S7. Alternative and stage-specific TSSs, uORFs and translation efficiency', intro7, s7)
print(p7, m7)

# ======================= S8 =======================
s8 = []
VA = 'bidirectional_mnase_pairs.tsv'; BM = rd(VA); VB = 'divergent_pairs_mnase_ndr.tsv'; DM = rd(VB); VC = 'divergent_pairs_coexpression_classes.tsv'; CX = rd(VC)
assert len(BM) == 994 and len(DM) == 1356 and len(CX) == 1356
NDR = 'MNase NDR class (Kensche 2016 via nucDetective; 8-time-point mean, 21-bp smoothing; NDR = <0.5 × genome median, ≥50 bp; within 300 bp upstream of each TSS): shared / separate / one_side_only / none.'
ATC = 'ATAC (Toenhake 2018) NDR class: 101-bp smoothing, genome-wide 90th percentile, 300-bp upstream windows; shared / separate / one_side_only / none.'
def pdesc(c):
    fixed = {'gene_left': 'Left gene of the pair (lower coordinate; minus strand for divergent pairs).', 'gene_right': 'Right gene (plus strand for divergent pairs).', 'chrom': 'Chromosome.', 'orientation': 'Pair orientation.', 'gene_gap_bp': 'Gap between the two genes (bp).',
             'TSS_left_final': 'Refined TSS of the left gene (0-based; long-read or model).', 'TSS_right_final': 'Refined TSS of the right gene.', 'TSS_dist_final': 'Distance between refined TSSs (bp).', 'TSS_left_v71': 'v71 TSS, left gene.', 'TSS_right_v71': 'v71 TSS, right gene.',
             'dist_v71': 'Distance between v71 TSSs (bp).', 'TSS_dist_v71': 'Distance between v71 TSSs (bp).', 'r': 'Pearson r of log2 DAFT-seq expression over 7 time points.', 'pearson_DAFT_log2': 'Pearson r of log2 DAFT-seq expression over 7 time points.',
             'min_DAFT_max': 'Smaller of the two genes\' maximum DAFT CPM.', 'cat': 'Co-expression class (anti / uncorrelated −0.3–0.3 / moderate 0.3–0.7 / strong).', 'mnase_P': NDR + ' Primary rule.', 'ndr_ATAC_refinedRule_anyStage': ATC + ' Any stage, anchoring rule TSSs.',
             'n_NDR_interTSS': 'Number of MNase NDRs between the two TSSs.', 'n_nuc_interTSS_allNucs': 'Reference-map nucleosome dyads between the two TSSs.', 'n_nuc_interTSS_top20': 'Same, top-20% best-positioned nucleosomes.', 'mn_occ_rel': 'Mean Kensche MNase occupancy in ±150 bp around the midpoint between the two TSSs, divided by the genome median.',
             'mnase_mean_occ_interTSS': 'Mean MNase occupancy between the TSSs.', 'atac_mean': 'Mean ATAC signal (model ATAC tracks) in ±150 bp around the midpoint between the two TSSs.', 'nuc_per_kb_interTSS_allNucs': 'Nucleosomes per kb between the TSSs.', 'gap_between_separate_NDRs_bp': 'Gap between two separate NDRs (bp).',
             'mnase_S1_singleNDRbetween': 'Secondary rule S1 (one NDR between TSSs; preregistered but uninformative).', 'mnase_S1b_posthoc_gap_lt147': 'Post hoc rule S1b: separate NDRs with gap <147 bp counted as shared.', 'mnase_S3_anyTP': 'NDR class from any single time point.', 'mnase_S2_dyadgap_anyTP': 'Dyad-gap NDR class, any time point.',
             'mnase_P_anchoringRule': 'Primary MNase class with anchoring-rule TSSs (reference genes: long-read TSS; others: chromatin-free model TSS).', 'TSS_left_anchoringRule': 'Left TSS under the anchoring rule.', 'TSS_right_anchoringRule': 'Right TSS under the anchoring rule.',
             'TSS_dist_anchoringRule': 'TSS distance under the anchoring rule.', 'TSS_dist_refinedRule': 'TSS distance (bp) with refined-rule anchors (anchors used for chromatin read-outs: long-read TSS for reference genes, otherwise chromatin-free model TSS; to verify). Source of the manuscript numbers: 940 pairs with both distances, median 469 bp (v71) vs 920 bp, overlapping/abutting 105 vs 11. Saved upstream with 4 significant digits.', 'phase_diff_h': 'Difference in peak-expression phase (h).', 'dist': 'TSS distance used in the regression (bp).',
             'atac_finalAnchors_anyStage': ATC + ' Refined (final) TSS anchors.', 'atac_anchoringRule_anyStage': ATC + ' Anchoring-rule TSSs.', 'atg_left': 'ATG of left gene.', 'atg_right': 'ATG of right gene.', 'stop_left': 'Stop of left gene.', 'stop_right': 'Stop of right gene.',
             'strand_left': 'Strand, left gene.', 'strand_right': 'Strand, right gene.', 'src_left': 'Refined TSS source, left.', 'src_right': 'Refined TSS source, right.', 'v71utr_left': 'Left gene has a v71 5′UTR.', 'v71utr_right': 'Right gene has a v71 5′UTR.', 'ATG_dist': 'Distance between ATGs (bp).',
             'CDS_gap': 'Gap between CDSs (bp).', 'DAFT_max_left': 'Max DAFT CPM, left.', 'DAFT_max_right': 'Max DAFT CPM, right.', 'NDR_refinedRule_anyStage': ATC + ' Any stage.', 'NDR_v71_anyStage': ATC + ' Using v71 TSSs.', 'phase_L': 'Peak phase (h), left.', 'phase_R': 'Peak phase (h), right.',
             'min_amp': 'Smaller expression amplitude.', 'min_expr': 'Smaller expression level.', 'any_het': 'Either gene in heterochromatin (HP1 domain).', 'ndr': ATC + ' (column used in the original analysis).',
             'd_left_refined': 'Strand-aware distance (bp) from the nearest MNase NDR downstream edge (exit; NDRs within 3 kb upstream of the ATG) to the refined TSS (positive = TSS downstream of the exit).', 'd_left_v71': 'Same with the v71 TSS.', 'd_right_refined': 'As d_left_refined for the right gene.', 'd_right_v71': 'Same, v71.', 'hit_left_refined': 'Left refined TSS lies 0–150 bp downstream of an NDR exit.',
             'hit_left_v71': 'Same, v71.', 'hit_right_refined': 'Right refined TSS lies 0–150 bp downstream of an NDR exit.', 'hit_right_v71': 'Same, v71.', 'spans_left_refined': 'An entire MNase NDR lies inside the left gene\'s 5′UTR (TSS..ATG), refined TSS.', 'spans_left_v71': 'Same, v71.', 'spans_right_refined': 'Same for the right gene.', 'spans_right_v71': 'Same, v71.'}
    if c in fixed: return fixed[c]
    m = re.match(r'mnase_(S3|S2_dyadgap)_T(\d+)$', c)
    if m: return f'MNase NDR class, {"single time point" if m.group(1) == "S3" else "dyad-gap rule"}, T{m.group(2)}.'
    m = re.match(r'NDR_(refinedRule|v71)_(Ring|Trophozoite|Schizont)$', c)
    if m: return f'ATAC NDR class, {m.group(2)}, {"refined" if m.group(1) == "refinedRule" else "v71"} TSSs.'
    m = re.match(r'(ATAC|MNase)_meanStages_interTSS_(refinedRule|v71)$', c)
    if m: return f'Mean {m.group(1)} signal (model input tracks, mean of 3 stages; to verify) between the two TSSs, {m.group(2)} anchors (saved with 4 significant digits).'
    return None
for nm, df, src, nt in (('MNase_NDR_994_pairs', BM, f'bidirectional_mnase_pairs.tsv ({VA})', 'Figure 5F; 994 divergent pairs with MNase classification'),
                        ('Divergent_MNase_vs_ATAC', DM, f'divergent_pairs_mnase_ndr.tsv ({VB})', '1,356 divergent pairs; all MNase rules, time points and ATAC comparison'),
                        ('Divergent_coexpression', CX, f'divergent_pairs_coexpression_classes.tsv ({VC})', '1,356 divergent pairs; spacing (refined vs v71), co-expression r, ATAC NDR')):
    s8.append(dict(name=nm, df=df, desc={c: pdesc(c) for c in df.columns if pdesc(c)}, source=src, note=nt))
SMN = rd('divergent_mnase_ndr_summary.tsv'); SCX = rd('coexpression_summary.tsv')
dsm = {c: 'Summary statistic (see "table"/"analysis" columns; definitions in report_divergent_mnase_ndr.md report_divergent_mnase_ndr.md).' for c in SMN.columns}
dsm.update({'class': 'NDR definition evaluated.', 'table': 'Sub-table of the summary.', 'beta_shared': 'Fisher-z regression coefficient for "shared" (adjusted for log10 distance and log10 min_DAFT_max, HC3).', 'beta_CI_lo': '95% CI lower.', 'beta_CI_hi': '95% CI upper.', 'kappa': 'Cohen κ vs ATAC.'})
s8.append(dict(name='MNase_NDR_summary', df=SMN, desc=dsm, source='divergent_mnase_ndr_summary.tsv (divergent_mnase_ndr_summary.tsv)', note='κ, regressions, class fractions'))
s8.append(dict(name='Coexpression_summary', df=SCX, desc={'table': 'Sub-table.', 'expr_filter_min_DAFTmax': 'Minimum DAFT max CPM of both genes.', 'group': 'Pair class (divergent / tandem / convergent / random).', 'n': 'Pairs.', 'mean_r': 'Mean Pearson r.', 'ci_lo': '95% CI lower.', 'ci_hi': '95% CI upper.', 'median_r': 'Median r.'},
    source='coexpression_summary.tsv (coexpression_summary.tsv)', note='earlier summary with DAFT expression filters (min DAFT max ≥1 or ≥5); its medians (e.g. divergent 0.55) differ from Figure 5E / manuscript (0.50) - use Neighbour_coexpr_Fig5E for the paper'))
F5E = pd.read_csv(f'{FIG}/Figure5/data/Figure5_E_neighbour_coexpression.tsv', sep='\t'); F5G = pd.read_csv(f'{FIG}/Figure5/data/Figure5_G_pairs_by_distance_and_NDR.tsv', sep='\t')
med = F5E.groupby('group').pearson_r.median().round(2).to_dict()
assert med == {'convergent': 0.14, 'divergent': 0.5, 'random_interchromosomal': 0.09, 'tandem': 0.22}, med
s8.insert(0, dict(name='Neighbour_coexpr_Fig5E', df=F5E, desc={'group': 'Pair orientation: divergent (head-to-head) / tandem / convergent (tail-to-tail) / random_interchromosomal.', 'gene_left': 'Left gene.', 'gene_right': 'Right gene.', 'pearson_r': 'Pearson r of log2 DAFT-seq expression over 7 time points.'},
    source='Figure5/data/Figure5_E_neighbour_coexpression.tsv', note='medians: head-to-head 0.50, tandem 0.22, tail-to-tail 0.14, random 0.09 (manuscript)'))
s8.insert(1, dict(name='Pairs_by_distance_NDR_Fig5G', df=F5G, desc={'gene_left': 'Left gene.', 'gene_right': 'Right gene.', 'TSS_distance_bp': 'Distance between refined TSSs (bp).', 'distance_bin': 'Distance bin.', 'n_NDR_interTSS': 'MNase NDRs between the TSSs.', 'ndr_class': '1 NDR / ≥2 NDRs / none.', 'pearson_r': 'Co-expression r.'},
    source='Figure5/data/Figure5_G_pairs_by_distance_and_NDR.tsv', note='Figure 5G; 300–800 bp: Mann–Whitney p = 0.77 (Figure5_statistics.tsv)'))
intro8 = ['Head-to-head (divergent) gene pairs: spacing between refined and v71 TSSs, MNase- and ATAC-defined NDRs between the TSSs, nucleosome counts, and DAFT-seq co-expression (Figure 5E–G, Fig S8).',
          'Exact TSS coordinates were checked against Table S4 (all non-missing TSS_left/right_final values match). TSS_dist_refinedRule and the ATAC/MNase mean columns were saved with 4 significant digits upstream.']
p8, m8 = write_table('TableS8_head_to_head_pairs.xlsx', 'Table S8. Head-to-head gene pairs', intro8, s8)
print(p8, m8)

# ======================= S9 =======================
VS = 'apiap2_genebody_track_summary.tsv'; AP = rd(VS)
assert len(AP) == 39 and AP.track_class.value_counts().to_dict() == {'genebody_or_TSS_depleted': 29, 'promoter_type_diffuse': 7, 'promoter_type_position_constrained': 3}
dA = {'track': 'Track (factor_stage).', 'factor': 'ApiAP2 factor.', 'stage': 'Stage.', 'gene_id': 'PlasmoDB ID of the factor (AP2-11B = PF3D7_1115500 inferred from naming, not confirmed).', 'track_class': 'Aggregate-profile class: genebody_or_TSS_depleted (29) / promoter_type_diffuse (7) / promoter_type_position_constrained (3). Promoter-type tracks (10) are used in Figure 5B.',
      'AP2_mRNA_CPM_stage': 'DAFT-seq CPM of the factor in the ChIP stage.', 'AP2_mRNA_pct_stage': 'Genome-wide percentile of that expression.', 'AP2_mRNA_rel_to_own_max': 'Expression relative to the factor\'s own IDC maximum.', 'frac_bp_nonzero': 'Fraction of genome with non-zero clipped signal.',
      'max_10bp': 'Maximum 10-bp bin value.', 'top1pct_bins_signal_share': 'Share of clipped signal in the top 1% bins.', 'IP_reads_total': 'IP fragments (500-bp bins, featureCounts).', 'input_reads_total': 'Input fragments.', 'IP_top1pct_share': 'Share of IP reads in the top 1% bins.',
      'input_top1pct_share': 'Same for input.', 'rho_IP_vs_input_bins': 'Spearman ρ IP vs input across bins.', 'rep_pearson_log2CPM_bins': 'Replicate Pearson r of log2 CPM.', 'rep_pearson_log2ratio_bins': 'Replicate Pearson r of log2(IP/input).',
      'IP_GC_top_vs_bottom_decile_log2': 'IP log2 CPM difference between highest and lowest GC deciles (GC slope).', 'input_GC_top_vs_bottom_decile_log2': 'Same for input.', 'rho_log2ratio_vs_GC': 'Spearman ρ of log2(IP/input) with 500-bp GC content.',
      'log2ratio_promoter': 'log2(IP/input) in promoters (TSS−500..TSS).', 'log2ratio_body': 'log2(IP/input) in gene bodies (ATG..stop).', 'log2ratio_interg': 'log2(IP/input) in intergenic bins (≥500 bp from CDS).', 'ratio_class_promoter_minus_body': 'Promoter − body log2 ratio.',
      'frac_promoter_bins_IP_below_input': 'Fraction of promoter bins with IP CPM < input CPM.', 'frac_genic_bins_IP_below_input': 'Same for genic bins.', 'GCcorr_promoter': 'GC-corrected log2 ratio, promoters.', 'GCcorr_genic': 'GC-corrected, genic.', 'GCcorr_P_minus_B': 'GC-corrected promoter − body.',
      'GCcorr_HP1_minus_rest': 'GC-corrected excess in HP1-domain bins.', 'rho_LR_vs_GFPmock': 'Spearman ρ of the genome-wide log2-ratio pattern with the GFP mock IP.', 'rho_LR_vs_IgG': 'Same with IgG (Pol II series).', 'rho_body_expr_stage': 'ρ of body signal with expression in the ChIP stage.',
      'rho_body_PolIIS2P': 'ρ of body signal with Pol II Ser2P.', 'rho_body_H3': 'ρ with H3.', 'rho_body_GC': 'ρ with gene-body GC.', 'pc1_loading': 'Loading on PC1 of gene-body signal across tracks.', 'body_GC_given_expr_HP1': 'Partial correlation body signal–GC given expression and HP1.',
      'body_expr_given_GC_HP1': 'Partial correlation body signal–expression given GC and HP1.', 'body_HP1_given_GC_expr': 'Partial correlation body signal–HP1 given GC and expression.', 'bins_rho_GC_all': 'ρ log2 ratio vs GC, all bins.', 'bins_rho_GC_intergenic': 'Same, intergenic bins.',
      'HP1bin_enrichment': 'Enrichment in HP1-domain bins.', 'frac_signal_in_HP1_bins': 'Fraction of signal in HP1-domain bins.', 'summit_region': 'Region class of the top summits.', 'lit_motif': 'Literature motif tested.', 'OR_lit': 'Odds ratio of the literature motif at summits vs GC-matched controls.',
      'q_lit_BH': 'BH q.', 'top_kmer': 'Top enriched 6-mer at summits.', 'top_kmer_OR': 'Its odds ratio.', 'top_kmer_q': 'Its q.', 'IP_minus_input_top1pct': 'IP − input top-1% read share (percentage points).', 'shang2022_AP2HF': 'Factor is an AP2-HF (heterochromatin-associated) factor in Shang 2022.',
      'assigned_causes': 'Assigned causes: T1 GC bias; T2 low enrichment / mock-like; B1 heterochromatin; B2 sparse true promoter peaks; B3 possible in-body motif binding; S1 low AP2 mRNA.', 'raw_body_LR_rho_expr': 'Uncorrected body log2 ratio vs expression ρ.', 'raw_body_LR_rho_GC': 'Uncorrected body log2 ratio vs GC ρ.',
      'raw_body_LR_expr_given_GC_HP1': 'Uncorrected partial r with expression.', 'raw_body_LR_GC_given_expr_HP1': 'Uncorrected partial r with GC.', 'top50_frac_in_promoter1kb': 'Fraction of the 50 strongest peaks within 1 kb upstream of a TSS (genome expectation 0.204).',
      'top50_frac_in_body': 'Fraction in gene bodies.', 'top50_height_over_median_peak': 'Height of top-50 peaks over the median peak.', 'binom_p_promoter': 'Binomial p for promoter enrichment of top-50 peaks.', 'q_promoter_BH': 'BH q.', 'primary_explanation': 'Primary explanation of the track\'s profile.'}
intro9 = ['ApiAP2 ChIP-seq track summary (Shang 2022; 13 factors × 3 stages = 39 tracks): profile class and QC metrics explaining why 29 tracks show gene-body signal (GC bias of IP libraries; heterochromatin binding in some AP2-HF factors) and the selection of the 10 promoter-type tracks used in Figure 5B (Fig S6).',
          'Source analysis: report_apiap2_genebody.md (report_apiap2_genebody.md). AP2-06A, AP2-G4 and AP2-O2 (single-stage tracks) are not included.']
p9, m9 = write_table('TableS9_ApiAP2_tracks.xlsx', 'Table S9. ApiAP2 ChIP-seq track classification and QC', intro9, [dict(name='ApiAP2_tracks', df=AP, desc=dA, source=f'apiap2_genebody_track_summary.tsv ({VS})', note='39 tracks')], coord=False)
print(p9, m9)
