"""Session-derived computation of gene_expression_strata_calls.tsv.gz.

Recovered from the record of the interactive analysis session in which gene_expression_strata_calls.tsv.gz was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: stage_TPM_9stages.tsv, gene_boundaries_v3_all_genes.tsv, benchmark_method_predictions.tsv, benchmark_truth_sets.tsv, siegel2014_pas_per_gene.tsv, ml3_end_positives.tsv
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import result  # paths: config/paths.yaml
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

B = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep='\t')
tpm = pd.read_csv(result('stage_TPM_9stages.tsv'), sep='\t')
TR = pd.read_csv(result('benchmark_truth_sets.tsv'), sep='\t')
SG = pd.read_csv(result('siegel2014_pas_per_gene.tsv'), sep='\t')
MX = pd.read_csv(result('benchmark_method_predictions.tsv'), sep='\t')
P3 = pd.read_csv(result('ml3_end_positives.tsv'), sep='\t')

T = tpm.set_index('Geneid')
D = B[['gene_id', 'chrom', 'strand', 'atg', 'stop', 'lr_group', 'lr_TSS', 'lr_TTS']].copy()
D['IDC_TPM'] = D.gene_id.map(T[['Ring', 'Trophozoite', 'Schizont']].max(axis=1))
D['tau'] = D.gene_id.map(T.tau)
D['peak_stage'] = D.gene_id.map(T.peak_stage)
D['tau_IDC'] = D.gene_id.map(T[['Ring', 'Trophozoite', 'Schizont']].pipe(lambda x: (1 - x.div(x.max(axis=1), axis=0)).sum(axis=1) / 2))
D['train'] = D.gene_id.isin(set(P3.locus))
D['confirmed'] = D.lr_TSS.notna()
D = D.merge(MX[['gene_id', 'Ours_chromatin__TSS', 'Ours_chromatin__TTS', 'Ours_nochromatin__TSS', 'Ours_nochromatin__TTS', 'PlasmoDB_v71__TSS', 'PlasmoDB_v71__TTS', 'FixedOffset__TSS', 'FixedOffset__TTS']], on='gene_id', how='left')
D = D.merge(TR[['gene_id', 'Consensus_5end__TSS']], on='gene_id', how='left').merge(SG[['gene_id', 'PAS_major']], on='gene_id', how='left')
D['logTPM'] = np.log10(D.IDC_TPM.fillna(0) + 1)
D['Qall'] = pd.qcut(D.IDC_TPM.rank(method='first'), 5, labels=['Q1 (lowest)', 'Q2', 'Q3', 'Q4', 'Q5 (highest)'])

SEX = {'Gametocyte', 'Gii', 'GV', 'Ookinete'}
D['stage_class'] = np.select([D.peak_stage.isin(SEX), D.tau_IDC >= D.tau_IDC.quantile(2 / 3)], ['Peak in sexual/mosquito stage', 'Blood-stage specific (top tertile τ)'], 'Broadly expressed in blood stages')
for c_ in ('TSS', 'TTS'):
    D[f'OursFB__{c_}'] = D[f'Ours_chromatin__{c_}'].fillna(D[f'FixedOffset__{c_}'])

D.drop(columns=['logTPM']).to_csv('gene_expression_strata_calls.tsv.gz', sep='\t', index=False)
