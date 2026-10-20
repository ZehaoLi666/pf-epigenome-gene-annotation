from pfconfig import FIGURES_DIR, result  # paths: config/paths.yaml
import pandas as pd, numpy as np, urllib.parse
import os
exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'common.py')).read())
FIG = FIGURES_DIR
V_G = 'gene_boundaries_v3_all_genes.tsv'; V_M = 'models_v3.tsv.gz'; V_VR = 'v71_vs_reference_ends.tsv'
G = pd.read_csv(result(V_G), sep='\t', low_memory=False)
G['description'] = G.desc.map(lambda s: urllib.parse.unquote(s) if isinstance(s, str) else s)
M = pd.read_csv(result(V_M), sep='\t', compression='gzip', low_memory=False)
VR = pd.read_csv(result(V_VR), sep='\t')
NUC = M.chrom.str.match(r'Pf3D7_\d\d_v3$')
REFM = M[NUC & M.conf3_any & M.full_cds & (M.cgroup == 'mRNA')].copy()
assert len(REFM) == 9142 and REFM.locus.nunique() == 2389
S = ['ring', 'troph', 'schizont']

# ======================= S2 =======================
REFM['TSS'] = np.where(REFM.strand == '+', REFM.start, REFM.end - 1); REFM['TTS'] = np.where(REFM.strand == '+', REFM.end - 1, REFM.start)
models = REFM[['transcript_id', 'locus', 'chrom', 'start', 'end', 'strand', 'TSS', 'TTS', 'n_exons', 'tier', 'structural_category', 'struct_ok', 'tes_ok']
              + [f'{p}_{s}' for s in S for p in ('FL', 'detect', 'sj', 'conf3')] + ['cds_lo', 'cds_hi']].rename(
    columns={'locus': 'gene_id', 'start': 'start_0based_BED', 'end': 'end_BED_exclusive', 'TSS': 'TSS_0based', 'TTS': 'TTS_0based', 'cds_lo': 'CDS_start_0based_BED', 'cds_hi': 'CDS_end_BED_exclusive'})
for c in ['CDS_start_0based_BED', 'CDS_end_BED_exclusive']: models[c] = models[c].astype('Int64')
models = models.sort_values(['chrom', 'start_0based_BED', 'transcript_id']).reset_index(drop=True)
cnt = REFM.groupby('locus').agg(n_reference_models=('transcript_id', 'size'), **{f'n_models_confirmed_{s}': (f'conf3_{s}', 'sum') for s in S})
genes = G[G.lr_group == 'confirmed'][['gene_id', 'chrom', 'strand', 'description', 'atg', 'stop', 'lr_TSS', 'lr_TTS', 'annot_tss', 'annot_tts']].merge(
    VR, on='gene_id', how='left').merge(cnt, left_on='gene_id', right_index=True, how='left')
genes = genes.rename(columns={'lr_TSS': 'ref_primary_TSS', 'lr_TTS': 'ref_primary_TTS', 'annot_tss': 'v71_TSS', 'annot_tts': 'v71_TTS', 'L5_v71': 'v71_5UTR_len', 'L3_v71': 'v71_3UTR_len',
                              'ref_TSS_dists': 'ref_5end_cluster_dists_from_ATG', 'ref_TTS_dists': 'ref_3end_cluster_dists_from_stop',
                              'c5_100': 'class_5prime_100bp', 'c3_100': 'class_3prime_100bp', 'c5_200': 'class_5prime_200bp', 'c3_200': 'class_3prime_200bp'})
for c in ['ref_primary_TSS', 'ref_primary_TTS']: genes[c] = genes[c].astype('Int64')
assert len(genes) == 2389 and genes.n_reference_models.sum() == 9142
bad = lambda s: ~s.isin(['agree', 'not assessable'])
both = ((genes.class_5prime_100bp == 'agree') & (genes.class_3prime_100bp == 'agree')).sum()
dis100 = (bad(genes.class_5prime_100bp) | bad(genes.class_3prime_100bp)).sum(); dis200 = (bad(genes.class_5prime_200bp) | bad(genes.class_3prime_200bp)).sum()
assert (both, dis100, dis200) == (572, 1612, 1183), (both, dis100, dis200)
summ = []
for end, c1, c2 in (('5′', 'class_5prime_100bp', 'class_5prime_200bp'), ('3′', 'class_3prime_100bp', 'class_3prime_200bp')):
    for tol, c in ((100, c1), (200, c2)):
        for k, v in genes[c].value_counts().items(): summ.append(dict(end=end, tolerance_bp=tol, v71_class=k, n_genes=int(v)))
summ += [dict(end='both', tolerance_bp=100, v71_class='both ends agree', n_genes=int(both)), dict(end='either', tolerance_bp=100, v71_class='≥1 end disagrees (assessable)', n_genes=int(dis100)),
         dict(end='either', tolerance_bp=200, v71_class='≥1 end disagrees (assessable)', n_genes=int(dis200))]
summ = pd.DataFrame(summ)
funnel = pd.read_csv(f"{FIG}/Figure1/data/Figure1_A_reference_funnel.tsv", sep='\t')
cm = pd.read_csv(f"{FIG}/Figure1/data/Figure1_G_class_matrix.tsv", sep='\t')
comp = pd.read_csv(f"{FIG}/Figure1/data/Figure1_C_model_completeness.tsv", sep='\t')
CLS = ('v71 UTR end vs all long-read end clusters of the gene (primary and alternative), tolerance as in the column name: agree = within tolerance of any cluster; over-long = v71 end beyond the farthest cluster; '
       'too short = v71 end inside the nearest cluster; no annotated UTR = v71 has no UTR on this side; between ends = between two clusters but within tolerance of neither; not assessable = no qualifying long-read end (Figure 1F rules).')
dg = {'gene_id': 'PlasmoDB gene ID.', 'chrom': 'Chromosome.', 'strand': 'Strand.', 'description': 'PlasmoDB product description (v68 GFF).', 'atg': 'ATG position (0-based).', 'stop': 'Last base of stop codon (0-based).',
      'ref_primary_TSS': 'Primary long-read TSS of the reference (0-based; highest-support 5′ end cluster).', 'ref_primary_TTS': 'Primary long-read TTS (0-based).',
      'v71_TSS': 'PlasmoDB v71 transcript start (0-based; = ATG when v71 has no 5′UTR).', 'v71_TTS': 'PlasmoDB v71 transcript end (0-based; = stop when no 3′UTR).',
      'class_5prime_100bp': CLS, 'class_3prime_100bp': CLS, 'class_5prime_200bp': CLS, 'class_3prime_200bp': CLS,
      'excess5_content': 'For over-long v71 5′UTRs: what lies in the excess region beyond the farthest reference end (catalogued lncRNA / antisense unit only / same-strand long-read unit / no other unit detected).',
      'excess3_content': 'As excess5_content for the 3′ side.', 'excess5_bp': 'Length (bp) of v71 5′UTR beyond the farthest long-read 5′ end.', 'excess3_bp': 'Length (bp) of v71 3′UTR beyond the farthest long-read 3′ end.',
      'v71_5UTR_len': 'v71 5′UTR length (bp; 0 = none).', 'v71_3UTR_len': 'v71 3′UTR length (bp; 0 = none).',
      'ref_5end_cluster_dists_from_ATG': 'Distances (bp, upstream of ATG) of all long-read 5′ end clusters of the gene, ";"-separated (Figure 1H).', 'ref_3end_cluster_dists_from_stop': 'Distances (bp, downstream of stop) of all long-read 3′ end clusters, ";"-separated (Figure 1I).',
      'n_reference_models': 'Number of reference (stage-confirmed, full-CDS) models for the gene.'}
for s in S: dg[f'n_models_confirmed_{s}'] = f'Number of the gene\'s reference models confirmed in the {s} stage.'
dm = {'transcript_id': 'Long-read model ID (TAMA merge).', 'gene_id': 'PlasmoDB gene ID of the locus.', 'chrom': 'Chromosome.', 'start_0based_BED': 'Model start, BED-style (0-based, inclusive).',
      'end_BED_exclusive': 'Model end, BED-style (exclusive; = 1-based end).', 'strand': 'Strand.', 'TSS_0based': 'Model 5′ end (0-based).', 'TTS_0based': 'Model 3′ end (0-based).', 'n_exons': 'Number of exons.',
      'tier': 'Evidence tier in the merged catalogue (Gold / Silver).', 'structural_category': 'SQANTI3 structural category vs PlasmoDB v48 annotation.', 'struct_ok': 'SQANTI3 structure passes (not fusion, intron-retention fragment or RT-switching artefact).',
      'tes_ok': '3′ end supported by ≥2 reads and not flagged as internal priming.', 'CDS_start_0based_BED': 'Annotated CDS envelope start (BED-style, 0-based).', 'CDS_end_BED_exclusive': 'Annotated CDS envelope end (BED-style, exclusive).'}
for s in S:
    dm[f'FL_{s}'] = f'Full-length PacBio Iso-Seq reads (GSE186109) supporting the model in {s}.'
    dm[f'detect_{s}'] = f'Detection call in {s}: detected / not_detected_powered / not_assessable.'
    dm[f'sj_{s}'] = f'Every intron supported by ≥2 short-read splice reads from the {s} stage (True for single-exon models).'
    dm[f'conf3_{s}'] = f'Model stage-confirmed in {s} (detected + struct_ok + sj + tes_ok; CAGE not used).'
ds = {'end': 'Transcript end.', 'tolerance_bp': 'Tolerance for "agree".', 'v71_class': 'Comparison class (see Reference_genes).', 'n_genes': 'Number of reference genes.'}
dfun = {c: c for c in funnel.columns}; dfun.update({'step': 'Construction step (1 merge, 2 stage confirmation, 3 complete gene models).', 'step_label': 'Step label.', 'row_label': 'Filter row.', 'n_loci_or_genes': 'Loci or genes remaining.', 'unit': 'Unit of the count.', 'n_models': 'Models remaining.', 'source': 'Origin of the number.'})
dcm = {c: 'Number of reference genes with this 3′ class (column) and 5′ class (row), 100-bp tolerance.' for c in cm.columns}; dcm[cm.columns[0]] = '5′ class (rows); columns are 3′ classes.'
dco = {'cls': 'Model class among stage-confirmed mRNA models.', 'n_models': 'Number of models.', 'n_genes_total': 'Genes with ≥1 stage-confirmed mRNA model.', 'label': 'Label in Figure 1C.', 'percent': 'Percent of models.'}
intro2 = ['Long-read reference: 2,389 nuclear protein-coding genes with ≥1 stage-confirmed model covering the whole annotated CDS (9,142 models), and the comparison of each gene\'s PlasmoDB v71 UTR ends with the long-read end clusters (Figure 1F–I).',
          'Checks reproduced from the saved data: 9,142 models / 2,389 genes; both ends agree with v71 at 100 bp: 572 genes; ≥1 assessable end disagrees: 1,612 at 100 bp and 1,183 at 200 bp (manuscript Results §1).',
          'Full exon structures of the models are in the reference GTF/BED (Data S file, prepared separately); this table gives model spans, ends and confirmation flags.']
p2, m2 = write_table('TableS2_longread_reference.xlsx', 'Table S2. Long-read reference gene set and comparison with PlasmoDB v71', intro2, [
    dict(name='Reference_genes', df=genes, desc=dg, source=f'gene_boundaries_v3_all_genes.tsv ({V_G}); v71_vs_reference_ends.tsv ({V_VR}); models_v3.tsv.gz ({V_M})', note='2,389 genes'),
    dict(name='Reference_models', df=models, desc=dm, source=f'models_v3.tsv.gz ({V_M}); filter: nuclear, cgroup = mRNA, conf3_any, full_cds', note='9,142 models'),
    dict(name='v71_class_summary', df=summ, desc=ds, source='computed from Reference_genes (counts only)', note='Figure 1G numbers'),
    dict(name='v71_class_matrix_100bp', df=cm, desc=dcm, source='Figure1/data/Figure1_G_class_matrix.tsv', note='Figure 1G'),
    dict(name='Reference_funnel', df=funnel, desc=dfun, source='Figure1/data/Figure1_A_reference_funnel.tsv', note='Figure 1A'),
    dict(name='Model_completeness', df=comp, desc=dco, source='Figure1/data/Figure1_C_model_completeness.tsv', note='Figure 1C')])
print(p2, m2)

# ======================= S3 =======================
V_E = 'ml3_end_positives.tsv'; V_F = 'flank_gaps_v3_longread_tiers.tsv'
E = pd.read_csv(result(V_E), sep='\t'); FG = pd.read_csv(result(V_F), sep='\t')
LC = pd.read_csv(f"{FIG}/Figure2/data/Figure2_B_label_counts.tsv", sep='\t')
assert len(E) == 10225 and (E.end_type == 'TSS').sum() == 6789 and (E.end_type == 'TTS').sum() == 3436
assert FG.tier.value_counts().to_dict() == {'A': 1225, 'spanned': 550, 'B': 489, 'C': 86}
E = E.rename(columns={'locus': 'gene_id', 'pos': 'pos_0based'})
FG2 = FG.rename(columns={'tts_pos': 'tts_pos_0based', 'tss_pos': 'tss_pos_0based'})
FG2['label_weight'] = FG2.tier.map({'A': 1.0, 'B': 0.6, 'C': 0.3, 'spanned': 0.0})
de = {'chrom': 'Chromosome.', 'strand': 'Strand of the labelled end.', 'end_type': 'TSS or TTS.', 'pos_0based': 'Label position (0-based; read-support-weighted mean of the end cluster, ≤50-bp clustering).',
      'weight': 'Training weight of the positive window (±100 bp). TSS terminal: max(0.3, read share) × (1.0 with cap-G, 0.5 without); TTS terminal: max(0.3, read share); truncated-model ends 1.0; ends inside the CDS 0.5.',
      'gene_id': 'PlasmoDB gene ID of the locus.', 'kind': 'terminal = end cluster of full-CDS reference models; truncated-model end = supported end of a truncated confirmed model outside the CDS; internal = supported end inside the CDS envelope.',
      'end_evidence': 'TSS: cap-G support (≥2 reads with a non-templated 5′ G and ≥20% of reads within ±20 bp). TTS: ≥2 poly(A)-tailed reads and ≥5% of reads, downstream A <60%. (Terminal TTS labels with False are still positive.)'}
dfg = {'gap_id': 'Same-strand gap ID.', 'chrom': 'Chromosome.', 'strand': 'Strand.', 'a': 'Gap start (0-based; from gaps_v3.tsv).', 'b': 'Gap end (0-based; from gaps_v3.tsv; half-open convention to verify).',
       'tts_pos_0based': 'TTS of the upstream long-read unit bounding the gap.', 'tss_pos_0based': 'TSS of the downstream long-read unit bounding the gap.',
       'left': 'Long reads in the unit on one side of the gap (summed over datasets).', 'right': 'Long reads in the unit on the other side.', 'span': 'Reads spanning the gap.',
       'n3': 'Read 3′ ends at the gap TTS.', 'n3_tail': 'Of n3, reads with a poly(A) tail.', 'n3_arich': 'Of n3, ends with A-rich downstream sequence.', 'n5_G1': 'Read 5′ ends at the gap TSS with a non-templated 5′ G.',
       'n3_drna': 'Read 3′ ends at the TTS from direct-RNA datasets.', 'downA_tts': 'Fraction A in the 20 nt downstream of the TTS.', 'no_span': 'True if span = 0, or span ≤2 and <2% of the weaker side.',
       'tts_ok': 'TTS evidenced (≥1 direct-RNA end, or ≥1 tailed read, or ≥2 non-A-rich ends with downstream A <60%).', 'tss_ok': 'TSS evidenced (≥1 cap-G read).', 'enough': '≥5 reads on both sides.',
       'tier': 'spanned (masked, weight 0) / A (both ends evidenced, weight 1.0) / B (one end, 0.6) / C (<5 reads on a side or no end evidence, 0.3).', 'label_weight': 'Flank label weight implied by tier (added for convenience).'}
dlc = {'bar': 'Label class.', 'bar_order': 'Plot order.', 'segment': 'Sub-category.', 'segment_order': 'Plot order.', 'n': 'Count.', 'bar_total': 'Total of the class.'}
intro3 = ['Training labels of the final (v3, CAGE-free) boundary model, built per strand from the long-read reference: 6,789 TSS and 3,436 TTS positive ends and 2,350 same-strand gaps (1,800 flank labels in tiers A/B/C; 550 gaps spanned by reads, masked).',
          'Gene-body labels are the union of full-CDS reference models (Table S2, Reference_models) plus bodies of truncated-only confirmed loci; they are not listed separately.',
          'Masked (weight 0): ±300 bp around every other model end, terminal TTS failing the A-content/flag test, unsupported truncated ends, and gaps spanned by reads. Positive windows are ±100 bp around each end; labels are Gaussian-softened (σ = 40 bp) at 20-bp resolution for training.',
          'Rule source: v3/labels_v3.py in v3_model_code_weights.tar.gz (v3_model_code_weights.tar.gz); final_model_labels_inputs.md (final_model_labels_inputs.md).']
p3, m3 = write_table('TableS3_training_labels.xlsx', 'Table S3. Training labels (TSS, TTS, flank) with weights and supporting evidence', intro3, [
    dict(name='End_labels', df=E, desc=de, source=f'ml3_end_positives.tsv ({V_E}) = Figure2/data/Figure2_B_source_end_positives.tsv', note='10,225 ends'),
    dict(name='Flank_gaps', df=FG2, desc=dfg, source=f'flank_gaps_v3_longread_tiers.tsv ({V_F}) = Figure2/data/Figure2_B_source_flank_gaps.tsv', note='2,350 gaps'),
    dict(name='Label_counts', df=LC, desc=dlc, source='Figure2/data/Figure2_B_label_counts.tsv', note='Figure 2B')])
print(p3, m3)

# ======================= S4 =======================
def chg(fin, ann):
    if pd.isna(fin): return np.nan
    if ann <= 0: return 'UTR added'
    d = fin - ann; return 'agrees' if abs(d) <= 100 else ('trimmed' if d < 0 else 'extended')
G['c5'] = [chg(f, a) for f, a in zip(G.utr5_final, G.utr5_annot.fillna(0))]; G['c3'] = [chg(f, a) for f, a in zip(G.utr3_final, G.utr3_annot.fillna(0))]
G['src5'] = np.where(G.lr_group == 'confirmed', 'long-read', np.where(G.final_TSS.notna(), 'model', 'none'))
G['src3'] = np.where(G.lr_group == 'confirmed', 'long-read', np.where(G.final_TTS.notna(), 'model', 'none'))
vc5 = G.c5.value_counts().to_dict(); vc3 = G.c3.value_counts().to_dict()
assert vc5 == {'trimmed': 1923, 'agrees': 1642, 'extended': 421, 'UTR added': 286} and vc3 == {'agrees': 1691, 'trimmed': 1243, 'extended': 1138, 'UTR added': 430}
assert G.c5.notna().sum() == 4272 and G.c3.notna().sum() == 4502
X = G.merge(VR[['gene_id', 'ref_TSS_dists', 'ref_TTS_dists']], on='gene_id', how='left')
per = pd.DataFrame({'gene_id': X.gene_id, 'chrom': X.chrom, 'strand': X.strand, 'description': X.description, 'atg': X.atg, 'stop': X.stop, 'longread_group': X.lr_group,
    'refined_TSS': X.final_TSS, 'TSS_source': X.src5, 'refined_5UTR_len': X.utr5_final, 'change_5prime_vs_v71': X.c5,
    'refined_TTS': X.final_TTS, 'TTS_source': X.src3, 'refined_3UTR_len': X.utr3_final, 'change_3prime_vs_v71': X.c3,
    'v71_TSS': X.annot_tss, 'v71_5UTR_len': X.utr5_annot, 'v71_TTS': X.annot_tts, 'v71_3UTR_len': X.utr3_annot,
    'model_TSS': X.model_TSS, 'model_TSS_prob': X.model_TSS_prob, 'model_alt_TSS': X.model_alt_TSS, 'n_model_TSS_calls': X.n_TSS_calls,
    'model_TTS': X.model_TTS, 'model_TTS_prob': X.model_TTS_prob, 'model_alt_TTS': X.model_alt_TTS, 'n_model_TTS_calls': X.n_TTS_calls,
    'chromfree_model_TSS': X.chromfree_model_TSS, 'chromfree_model_TTS': X.chromfree_model_TTS,
    'longread_primary_TSS': X.lr_TSS, 'longread_primary_TTS': X.lr_TTS, 'longread_5end_cluster_dists_from_ATG': X.ref_TSS_dists, 'longread_3end_cluster_dists_from_stop': X.ref_TTS_dists,
    'CAGE_cpm_at_model_TSS': X.CAGE_cpm_at_model_TSS, 'CAGE_pct_model_TSS': X.CAGE_pct_model_TSS, 'CAGE_pct_v71_TSS': X.CAGE_pct_v68_TSS, 'CAGE_pct_distance_matched_ctrl': X.CAGE_pct_distance_matched_ctrl,
    'dist_model_TSS_to_CAGE_peak': X.dist_model_TSS_to_CAGE_peak})
for c in ['refined_TSS', 'refined_TTS', 'refined_5UTR_len', 'refined_3UTR_len', 'model_TSS', 'model_TTS', 'chromfree_model_TSS', 'chromfree_model_TTS', 'longread_primary_TSS', 'longread_primary_TTS', 'dist_model_TSS_to_CAGE_peak']:
    per[c] = per[c].round().astype('Int64')
per.loc[per.refined_TSS.isna(), 'TSS_source'] = 'none (v71 structure kept)'; per.loc[per.refined_TTS.isna(), 'TTS_source'] = 'none (v71 structure kept)'
per = per.sort_values(['chrom', 'atg']).reset_index(drop=True)
CH = ['agrees', 'trimmed', 'extended', 'UTR added']; rows = []
for end, cc, sc in (('5′ end', 'c5', 'src5'), ('3′ end', 'c3', 'src3')):
    for src in ('long-read', 'model', 'all refined'):
        q = G[G[sc].isin(['long-read', 'model'])] if src == 'all refined' else G[G[sc] == src]; vc = q[cc].value_counts()
        rows.append(dict(end=end, source=src, n=len(q), **{f'n_{c}': int(vc.get(c, 0)) for c in CH}, **{f'frac_{c}': vc.get(c, 0) / len(q) for c in CH}))
    rows.append(dict(end=end, source='none (no refined end)', n=int((G[sc] == 'none').sum())))
SUM = pd.DataFrame(rows)
f4d = pd.read_csv(f"{FIG}/Figure4/data/Figure4_D_change_vs_v71.tsv", sep='\t')
for _, r in f4d.dropna().iterrows():
    q = SUM[(SUM.end == r['end']) & (SUM.source == r['source'])].iloc[0]
    assert q['n'] == r['n'] and all(abs(q[f'frac_{c}'] - r[c]) < 1e-8 for c in CH)
dp = {'gene_id': 'PlasmoDB gene ID (5,285 nuclear protein-coding genes).', 'chrom': 'Chromosome.', 'strand': 'Strand.', 'description': 'PlasmoDB product description.', 'atg': 'ATG (0-based).', 'stop': 'Last base of stop codon (0-based).',
      'longread_group': 'Long-read status: confirmed (reference gene, Table S2) / truncated-only / HC model, unconfirmed / no HC model.',
      'refined_TSS': 'Refined TSS (0-based): long-read primary TSS for reference genes, otherwise the model TSS (empty = no refined 5′ end; v71 kept).', 'TSS_source': 'Source of the refined TSS: long-read / model / none.',
      'refined_5UTR_len': '|ATG − refined TSS| (bp).', 'change_5prime_vs_v71': 'Refined 5′ end vs v71 (Figure 4D rule): agrees = |difference in UTR length| ≤100 bp; trimmed = refined UTR >100 bp shorter; extended = >100 bp longer; UTR added = v71 has no 5′UTR.',
      'refined_TTS': 'Refined TTS (0-based), as refined_TSS.', 'TTS_source': 'Source of the refined TTS.', 'refined_3UTR_len': '|refined TTS − stop| (bp).', 'change_3prime_vs_v71': 'As change_5prime_vs_v71 for the 3′UTR.',
      'v71_TSS': 'PlasmoDB v71 transcript start (0-based; = ATG if no 5′UTR).', 'v71_5UTR_len': 'v71 5′UTR length (bp; 0 = none).', 'v71_TTS': 'v71 transcript end (0-based; = stop if no 3′UTR).', 'v71_3UTR_len': 'v71 3′UTR length (bp; 0 = none).',
      'model_TSS': 'Highest-probability model TSS call within 3 kb upstream of the ATG (not crossing the upstream same-strand CDS), from the out-of-fold full model (V3_top8; each chromosome predicted by the fold that did not train on it).',
      'model_TSS_prob': 'Model probability (sigmoid output) at model_TSS.', 'model_alt_TSS': 'Other TSS calls in the window (alternative TSSs), comma-separated, in decreasing probability.', 'n_model_TSS_calls': 'Number of TSS calls in the window.',
      'model_TTS': 'Highest-probability model TTS call within 2 kb downstream of the stop codon.', 'model_TTS_prob': 'Model probability at model_TTS.', 'model_alt_TTS': 'Other TTS calls in the window (alternative TTSs).', 'n_model_TTS_calls': 'Number of TTS calls in the window.',
      'chromfree_model_TSS': 'TSS from the model trained without chromatin inputs (used as anchor for chromatin read-outs in non-reference genes).', 'chromfree_model_TTS': 'TTS from the chromatin-free model.',
      'longread_primary_TSS': 'Long-read primary TSS (reference genes only).', 'longread_primary_TTS': 'Long-read primary TTS (reference genes only).',
      'longread_5end_cluster_dists_from_ATG': 'Alternative long-read 5′ ends: distances (bp upstream of ATG) of all long-read 5′ end clusters, ";"-separated (reference genes).', 'longread_3end_cluster_dists_from_stop': 'Distances (bp downstream of stop) of all long-read 3′ end clusters (reference genes).',
      'CAGE_cpm_at_model_TSS': 'Adjalley 2016 CAGE CPM at model_TSS (validation only).', 'CAGE_pct_model_TSS': 'Percentile of CAGE signal at model_TSS within the gene\'s upstream window.', 'CAGE_pct_v71_TSS': 'Same percentile at the v71 (= v68) TSS.',
      'CAGE_pct_distance_matched_ctrl': 'Same percentile at a distance-matched control position.', 'dist_model_TSS_to_CAGE_peak': 'Distance (bp) from model_TSS to the strongest CAGE position in the window.'}
dsum = {'end': 'Transcript end.', 'source': 'Source of refined end.', 'n': 'Genes.'}
for c in CH: dsum[f'n_{c}'] = f'Genes classed "{c}".'; dsum[f'frac_{c}'] = f'Fraction classed "{c}" (equals Figure4_D_change_vs_v71.tsv).'
intro4 = ['Per-gene refined TSS/TTS for all 5,285 nuclear protein-coding genes, their source (long-read reference or model), model probabilities, alternative ends and the change relative to PlasmoDB v71.',
          'Reproduces Figure 4D exactly: 5′ end 1,642 agrees / 1,923 trimmed / 421 extended / 286 UTR added of 4,272 genes; 3′ end 1,691 / 1,243 / 1,138 / 430 of 4,502 genes; 1,013 genes without a refined 5′ end and 783 without a refined 3′ end keep the v71 structure (Figure4D_summary sheet; fractions checked against Figure4/data/Figure4_D_change_vs_v71.tsv).',
          'Model calls: gene_boundaries.py (smoothed-logit local maxima, P ≥ threshold, NMS) on genome-wide out-of-fold predictions (folds: chr 1–5, 6–10, 11–14).']
p4, m4 = write_table('TableS4_refined_gene_boundaries.xlsx', 'Table S4. Refined TSS and TTS per gene', intro4, [
    dict(name='Refined_boundaries', df=per, desc=dp, source=f'gene_boundaries_v3_all_genes.tsv ({V_G}); v71_vs_reference_ends.tsv ({V_VR}); change classes recomputed with the Figure 4D rule from Figure4/code/prep_data.py', note='5,285 genes'),
    dict(name='Figure4D_summary', df=SUM, desc=dsum, source='computed from Refined_boundaries; checked against Figure4/data/Figure4_D_change_vs_v71.tsv (Figure4_D_change_vs_v71.tsv)', note='Figure 4D')])
print(p4, m4)
