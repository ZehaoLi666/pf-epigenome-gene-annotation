"""Session-derived computation of genes_without_utr_prediction.tsv.

Recovered from the record of the interactive analysis session in which genes_without_utr_prediction.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: stage_TPM_9stages.tsv, gene_heterochromatin_groups.tsv, PlasmoDB-71_Pfalciparum3D7.gff, gene_boundaries_v3_all_genes.tsv, oof_top8_oof_prob.npz, daft_stage20.npz, oof_top8_genes.tsv.gz
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import glob, numpy as np, pandas as pd, pyBigWig
from urllib.parse import unquote

# Load GFF files
cols = ['chrom','src','type','start','end','score','strand','phase','attr']
def rd(p):
    g = pd.read_csv(p, sep='\t', comment='#', header=None, names=cols, low_memory=False)
    return g

g68 = rd(glob.glob(f'{LAB_DIR}/PlasmoDB/PlasmoDB-68/*.gff')[0])
g71 = rd(result('PlasmoDB-71_Pfalciparum3D7.gff'))
hdr71 = [l for l in open(result('PlasmoDB-71_Pfalciparum3D7.gff')) if l.startswith('##')][:30]

# Build name/desc tables from both GFF versions
def gdesc(g):
    q = g[g.type == 'protein_coding_gene']
    return pd.DataFrame({'gene_id': q.attr.str.extract(r'ID=([^;]+)')[0].values,
                         'desc': q.attr.str.extract(r'description=([^;]+)')[0].fillna('').map(unquote).values,
                         'name': q.attr.str.extract(r'Name=([^;]+)')[0].fillna('').values})

D = gdesc(g68).merge(gdesc(g71), on='gene_id', suffixes=('_68', '_71'))

# Load main gene boundaries table
B = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep='\t').merge(D[['gene_id', 'desc_71']], on='gene_id', how='left')

# Load oof_top8 model output
M8 = pd.read_csv(result('oof_top8_genes.tsv.gz'), sep='\t')[['gene_id', 'tss_win_lo', 'tss_win_hi', 'tts_win_lo', 'tts_win_hi', 'max_TSS_prob_window', 'max_TTS_prob_window']]
B = B.merge(M8, on='gene_id')

# Load TPM and heterochromatin tables
tpm = pd.read_csv(result('stage_TPM_9stages.tsv'), sep='\t')
het = pd.read_csv(result('gene_heterochromatin_groups.tsv'), sep='\t')

# Load DAFT-seq stage data
z = np.load(result('daft_stage20.npz'))
DK = {k: z[k] for k in z.files}

def daft(r):
    k = 0 if r.strand == '+' else 1; lo, hi = sorted((r.atg, r.stop))
    v = [float(np.expm1(DK[f'{s}__{r.chrom}'][k, lo // 20:hi // 20 + 1].astype(np.float32)).mean()) for s in ('ring', 'troph', 'schizont') if f'{s}__{r.chrom}' in DK]
    return max(v) if v else np.nan

B['daft_sense'] = [daft(r) for r in B.itertuples()]

T = tpm.set_index('Geneid')
B = B.join(T[['Ring', 'Trophozoite', 'Schizont', 'Gametocyte', 'ET', 'LT', 'Gii', 'GV', 'Ookinete', 'peak_stage']], on='gene_id')
B['IDC_TPM'] = B[['Ring', 'Trophozoite', 'Schizont']].max(axis=1)
B['sex_TPM'] = B[['Gametocyte', 'Gii', 'GV', 'Ookinete']].max(axis=1)

B = B.merge(het[['gene_id', 'family', 'HP1_domain', 'group']].rename(columns={'group': 'het_group'}), on='gene_id', how='left')

B['tss_win_len'] = B.tss_win_hi - B.tss_win_lo
B['tts_win_len'] = B.tts_win_hi - B.tts_win_lo

L = {l.split()[1]: l.split()[3] for l in hdr71 if 'sequence-region' in l}
B['dist_chr_end'] = [min(min(r.atg, r.stop), int(L[r.chrom]) - max(r.atg, r.stop)) for r in B.itertuples()]

B['expr_idc'] = (B.IDC_TPM.fillna(0) >= 10) & (B.daft_sense.fillna(0) >= 1)
B['expr_sex'] = B.sex_TPM.fillna(0) >= 10

SEX = {'Gametocyte', 'Gii', 'GV', 'Ookinete'}

# Load oof_top8 probability array
PRo = {k: v for k, v in np.load(result('oof_top8_oof_prob.npz')).items()}

# CAGE bigwig files
CB = {(s, st): pyBigWig.open(f"{LAB_DIR}/lncRNA project tracks/cage_tss/CTSS_{s}_{st}_CPM.bw")
      for s in ('ring', 'trophozoite', 'schizont') for st in ('pos', 'neg')}

def track(c, strand, a, b):
    k = 'pos' if strand == '+' else 'neg'
    return np.max([np.nan_to_num(np.array(CB[(s, k)].values(c, max(0, a), b), float)) for s in ('ring', 'trophozoite', 'schizont')], axis=0)

rng = np.random.default_rng(0)

def rescue(r):
    k = 0 if r.strand == '+' else 1; lo, hi = int(r.tss_win_lo), int(r.tss_win_hi)
    if hi - lo < 300: return (np.nan, np.nan, np.nan)
    p = PRo[r.chrom][k, 0, lo // 20:hi // 20 + 1].astype(np.float32); pos = lo + int(np.argmax(p)) * 20 + 10
    v = track(r.chrom, r.strand, lo - 100, hi + 100); sm = np.convolve(v, np.ones(201), 'valid')
    if sm.max() <= 0: return (np.nan, np.nan, 0.0)
    pc = lambda q: float((sm < sm[min(len(sm) - 1, max(0, q - lo))]).mean())
    return (pc(pos), pc(int(rng.integers(lo, hi))), float(sm.max()))

lab3 = '3 Not expressed in blood or sexual stages (incl. liver/sporozoite genes: LISP1/2, LSA1)'

def reason2(r):
    if r.het_group in ('var/rifin/stevor/Pfmc-2TM', 'other HP1+H3K9me3-high'): return '1 Heterochromatic / clonally variant'
    if not r.expr_idc:
        return '2 Not expressed in blood stages; expressed in gametocyte/mosquito stages' if r.expr_sex else lab3
    if r.tss_win_len < 500: return '4 upstream/downstream'
    if r.peak_stage in SEX: return '5a Expressed at low level in blood stages; peak in gametocyte/mosquito stages'
    if r.daft_sense < 5: return '5b Weakly expressed in blood stages (DAFT-seq <5 CPM)'
    return '5c Expressed in blood stages; signal below threshold'

B['no_TSS_reason'] = B.apply(lambda r: reason2(r) if r.final_src == 'none' else None, axis=1)

RS = {}
for cat in ['5c', '5b', '5a', '2', '1']:
    q = B[B.no_TSS_reason.str.startswith(cat, na=False)]
    rr = np.array([rescue(r) for r in q.itertuples()])
    RS[cat] = rr

def reason_tts(r):
    if not pd.isna(r.final_TTS): return None
    if r.het_group in ('var/rifin/stevor/Pfmc-2TM', 'other HP1+H3K9me3-high'): return '1 Heterochromatic / clonally variant'
    if not r.expr_idc:
        return '2 Not expressed in blood stages; expressed in gametocyte/mosquito stages' if r.expr_sex else lab3
    if r.tts_win_len < 300: return '4 upstream/downstream'
    if r.peak_stage in SEX: return '5a Expressed at low level in blood stages; peak in gametocyte/mosquito stages'
    if r.daft_sense < 5: return '5b Weakly expressed in blood stages (DAFT-seq <5 CPM)'
    return '5c Expressed in blood stages; TTS signal below threshold'

B['no_TTS_reason'] = B.apply(reason_tts, axis=1)

nm = D.set_index('gene_id').name_71
B['name71'] = B.gene_id.map(nm)

outc = ['gene_id', 'name71', 'desc_71', 'chrom', 'strand', 'atg', 'stop', 'lr_group', 'final_src', 'no_TSS_reason', 'no_TTS_reason',
        'IDC_TPM', 'Ring', 'Trophozoite', 'Schizont', 'Gametocyte', 'Gii', 'GV', 'Ookinete', 'peak_stage',
        'daft_sense', 'het_group', 'tss_win_len', 'tts_win_len', 'max_TSS_prob_window', 'max_TTS_prob_window']

B[B.no_TSS_reason.notna() | B.no_TTS_reason.notna()][outc].rename(
    columns={'name71': 'name_v71', 'desc_71': 'description_v71', 'daft_sense': 'DAFTseq_sense_CPM_max_stage'}
).to_csv('genes_without_utr_prediction.tsv', sep='\t', index=False)
