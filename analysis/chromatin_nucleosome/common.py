from pfconfig import LAB_DIR, RESULTS_DIR, result  # paths: config/paths.yaml
import numpy as np, pandas as pd
ART = RESULTS_DIR + '/'
GB = result('gene_boundaries_v3_all_genes.tsv')
EXPR = result('gene_daft_expression_truth.tsv.gz')
TRUTH = result('benchmark_truth_sets.tsv')
FASTA = f"{LAB_DIR}/IGV_longread_transcript_models/PlasmoDB-68_Pfalciparum3D7_Genome.fasta"
CHROMS = [f'Pf3D7_{i:02d}_v3' for i in range(1, 15)]
GROUPS = ['all', 'confirmed', 'truncated-only', 'HC model, unconfirmed', 'no HC model']

def chrom_sizes():
    L = {}; n = None
    for l in open(FASTA):
        if l[0] == '>':
            n = l[1:].split()[0]; L[n] = 0
        else:
            L[n] += len(l.strip())
    return L

def load_genes():
    g = pd.read_csv(GB, sep='\t')
    ex = pd.read_csv(EXPR, sep='\t')
    G = g.merge(ex[['gene_id', 'DAFT_max', 'Qd']], on='gene_id', how='left')
    conf = G.lr_group == 'confirmed'
    # anchoring rule: chromatin is the read-out -> confirmed genes on long-read ends, others on chromatin-free model
    G['ref_TSS'] = np.where(conf, G.lr_TSS, np.where(G.final_TSS.notna(), G.chromfree_model_TSS, np.nan))
    G['ref_TTS'] = np.where(conf, G.lr_TTS, np.where(G.final_TTS.notna(), G.chromfree_model_TTS, np.nan))
    G['v71_TSS'] = G.annot_tss.astype(float); G['v71_TTS'] = G.annot_tts.astype(float)
    G['v71_has5'] = G.utr5_annot > 0; G['v71_has3'] = G.utr3_annot > 0
    return G
