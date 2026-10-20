"""Derive every Figure 4 source table and write it to ../data/.

The computations are copied verbatim from make_fig4.py ,
the script that produced the composite Figure 4. Only the
plotting was removed and the outputs are written at full precision.

Inputs: intermediate result tables read from RESULTS_DIR (config/paths.yaml), or from a directory
given by --inputs DIR or the FIG4_INPUTS environment variable.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import os, sys, argparse
from pathlib import Path
import numpy as np, pandas as pd
from scipy import stats

INPUTS = {'benchmark_method_predictions.tsv': 'benchmark_method_predictions.tsv',
          'benchmark_truth_sets.tsv': 'benchmark_truth_sets.tsv',
          'gene_boundaries_v3_all_genes.tsv': 'gene_boundaries_v3_all_genes.tsv',
          'siegel2014_pas_per_gene.tsv': 'siegel2014_pas_per_gene.tsv',
          'chappell2020_5UTRseq_tagclusters.tsv': 'chappell2020_5UTRseq_tagclusters.tsv',
          'shaw2021_5CAPture_5CN_clusters.tsv': 'shaw2021_5CAPture_5CN_clusters.tsv',
          'daft_stratified_reference_vs_model.tsv': 'daft_stratified_reference_vs_model.tsv'}
# Adjalley 2016 cap-seq CTSS bigWigs (CPM, stranded, ring / trophozoite / schizont) on the lab machine
CAGE = os.environ.get('FIG4_CAGE', f'{LAB_DIR}/lncRNA project tracks/cage_tss/')
# Tables saved with the approved composite (used only to check that this script reproduces them)
APPROVED = {'A_tss_consensus_acc200': 'Figure4_A_tss_consensus_acc200.tsv', 'B_tts_siegel_acc200': 'Figure4_B_tts_siegel_acc200.tsv',
            'C_refined_vs_v71_signal': 'Figure4_C_refined_vs_v71_signal.tsv', 'C_profiles': 'Figure4_C_profiles.tsv',
            'D_change_vs_v71': 'Figure4_D_change_vs_v71.tsv', 'E_utr_length': 'Figure4_E_utr_length.tsv',
            'F_accuracy_by_expression': 'Figure4_F_accuracy_by_expression.tsv'}

HERE = Path(__file__).resolve().parent if '__file__' in globals() else Path.cwd()
OUT = HERE.parent / 'data'
ap = argparse.ArgumentParser(); ap.add_argument('--inputs', default=os.environ.get('FIG4_INPUTS'))
args, _ = ap.parse_known_args(sys.argv[1:] if __name__ == '__main__' else [])
def P(name):
    return str(Path(args.inputs) / name) if args.inputs else result(INPUTS[name])
def rd(name): return pd.read_csv(P(name), sep='\t')

import pyBigWig  # only needed here (panel C); the panel scripts need pandas/numpy/matplotlib/scipy only
COL = dict(agree='#4a78c2', ol_lnc='#e08214', short='#7b4fa0', noutr='#2a9d8f')

# ================================================================== verbatim from make_fig4.py (computation)
M = rd('benchmark_method_predictions.tsv'); T = rd('benchmark_truth_sets.tsv'); G = rd('gene_boundaries_v3_all_genes.tsv')
SP = rd('siegel2014_pas_per_gene.tsv'); CH = rd('chappell2020_5UTRseq_tagclusters.tsv'); SH = rd('shaw2021_5CAPture_5CN_clusters.tsv')
DSF = rd('daft_stratified_reference_vs_model.tsv')
X = M.merge(T, on='gene_id').merge(G[['gene_id', 'chrom', 'strand', 'atg', 'stop', 'lr_group', 'final_src', 'final_TSS', 'final_TTS',
                                      'annot_tss', 'annot_tts', 'utr5_annot', 'utr3_annot', 'utr5_final', 'utr3_final']], on='gene_id', how='left')
X = X.merge(SP[['gene_id', 'n_clusters', 'PAS_major']].rename(columns={'PAS_major': 'PAS_major_sp'}), on='gene_id', how='left')
X['v71all__TSS'] = X.annot_tss; X['v71all__TTS'] = X.annot_tts   # v71 as annotated (no UTR -> ATG / stop codon)
rng = np.random.default_rng(0); NB = 2000
LAB = {'Ours_chromatin': 'Our model', 'Ours_nochromatin': 'Our model, no chromatin', 'FixedOffset': 'Fixed UTR length', 'v71all': 'PlasmoDB v71',
       'ChromHMM': 'ChromHMM promoter', 'StringTie_short': 'StringTie, short reads', 'StringTie_long': 'StringTie, long reads'}

def acc_rows(df, meths, tcol, end):
    out = []
    for m in meths:
        h = ((df[f'{m}__{end}'] - df[tcol]).abs() <= 200).values.astype(float); bs = h[rng.integers(0, len(h), (NB, len(h)))].mean(1)
        out.append(dict(key=m, method=LAB[m], n=len(h), acc200=h.mean(), ci_lo=np.percentile(bs, 2.5), ci_hi=np.percentile(bs, 97.5),
                        median_abs_err=(df[f'{m}__{end}'] - df[tcol]).abs().median(), coverage=df[f'{m}__{end}'].notna().mean()))
    return pd.DataFrame(out)
Q = X[X.Consensus_5end__TSS.notna()]
TA = acc_rows(Q, ['Ours_chromatin', 'Ours_nochromatin', 'FixedOffset', 'v71all', 'ChromHMM', 'StringTie_short', 'StringTie_long'], 'Consensus_5end__TSS', 'TSS')
S = X[X.n_clusters > 0]
TB = acc_rows(S, ['Ours_chromatin', 'Ours_nochromatin', 'FixedOffset', 'v71all', 'StringTie_short', 'StringTie_long'], 'PAS_major_sp', 'TTS')
assert len(Q) == 1649 and len(S) == 3439
assert abs(TA.set_index('key').loc['Ours_chromatin', 'acc200'] - 0.7538) < 1e-3

# ---- C: independent signal centred on refined vs v71 ends
W = 1000; OFF = np.arange(-W, W + 1)
CL = pyBigWig.open(CAGE + 'CTSS_ring_pos_CPM.bw').chroms()
def dense_bw(tag):
    out = {}
    for c, L in CL.items():
        a = np.zeros(L, np.float32)
        for st in ('ring', 'trophozoite', 'schizont'):
            b = pyBigWig.open(CAGE + f'CTSS_{st}_{tag}_CPM.bw')
            if c in b.chroms(): a += np.nan_to_num(np.array(b.values(c, 0, L), np.float32))
        out[c] = a
    return out
def dense_pts(chrom, strand, pos, w):
    d = {s: {c: np.zeros(L, np.float32) for c, L in CL.items()} for s in '+-'}
    for c, s, p, v in zip(chrom, strand, pos, w):
        if c in CL and 0 <= p < CL[c]: d[s][c][int(p)] += v
    return d
SIG = {'Adjalley 2016': {'+': dense_bw('pos'), '-': dense_bw('neg')},
       'Chappell 2020': dense_pts(CH.Chromosome, CH.Strand, CH.Dominant_tag, CH.Quantile_TPM),
       'Shaw 2021': dense_pts(SH.chrom, SH.strand, SH.dominant, SH.score)}
ANC = {'Refined': ('final_TSS', 'final_TTS'), 'v71': ('v71all__TSS', 'v71all__TTS')}
ALAB = {'Refined': 'Refined annotation', 'v71': 'PlasmoDB v71'}
def oriented(arr, a, s):
    a = int(a); lo, hi = a - W, a + W + 1
    v = np.pad(arr[max(lo, 0):hi], (max(0, -lo), max(0, hi - len(arr))))
    return v if s == '+' else v[::-1]
PROF, NPROF, NSET = {}, {}, {}
for ds in SIG:
    PROF[ds], NPROF[ds] = {}, {}
    for k, (ct, _) in ANC.items():
        q = X[X[ct].notna()]; NSET[(k, 'TSS')] = len(q); Pl = []
        for r in q.itertuples():
            v = oriented(SIG[ds][r.strand][r.chrom], getattr(r, ct), r.strand); t = v.sum()
            if t > 0: Pl.append(v / t)
        PROF[ds][k] = np.mean(Pl, 0); NPROF[ds][k] = len(Pl)
PROF['Siegel 2014'], NPROF['Siegel 2014'] = {}, {}
for k, (_, ct) in ANC.items():
    q = X[X[ct].notna() & (X.n_clusters > 0)]; NSET[(k, 'TTS')] = int(X[ct].notna().sum())
    d = (q.PAS_major_sp - q[ct]).values * np.where(q.strand == '+', 1, -1)
    PROF['Siegel 2014'][k] = np.bincount(np.clip(d + W, 0, 2 * W).astype(int)[np.abs(d) <= W], minlength=2 * W + 1) / len(q)
    NPROF['Siegel 2014'][k] = len(q)
DSETS = ['Adjalley 2016', 'Chappell 2020', 'Shaw 2021', 'Siegel 2014']
TC = pd.DataFrame([dict(dataset=ds, anchor=ALAB[k], n_genes=NPROF[ds][k], frac_pm100=PROF[ds][k][np.abs(OFF) <= 100].sum(),
                        frac_pm50=PROF[ds][k][np.abs(OFF) <= 50].sum(), peak_offset_bp=int(OFF[np.argmax(PROF[ds][k])])) for ds in DSETS for k in ANC])

# ---- D: change of refined ends relative to v71
def chg(fin, ann):
    if pd.isna(fin): return np.nan
    if ann <= 0: return 'UTR added'
    d = fin - ann; return 'agrees' if abs(d) <= 100 else ('trimmed' if d < 0 else 'extended')
CHG = ['agrees', 'trimmed', 'extended', 'UTR added']
G['c5'] = [chg(f, a) for f, a in zip(G.utr5_final, G.utr5_annot.fillna(0))]; G['c3'] = [chg(f, a) for f, a in zip(G.utr3_final, G.utr3_annot.fillna(0))]
G['src5'] = np.where(G.lr_group == 'confirmed', 'long-read', np.where(G.final_TSS.notna(), 'model', 'none'))
G['src3'] = np.where(G.lr_group == 'confirmed', 'long-read', np.where(G.final_TTS.notna(), 'model', 'none'))
rowsD = []
for end, cc, sc in (('5′ end', 'c5', 'src5'), ('3′ end', 'c3', 'src3')):
    for src in ('long-read', 'model'):
        q = G[G[sc] == src]; vc = q[cc].value_counts()
        rowsD.append(dict(end=end, source=src, n=len(q), **{c: vc.get(c, 0) / len(q) for c in CHG}))
    rowsD.append(dict(end=end, source='none', n=int((G[sc] == 'none').sum())))
TD = pd.DataFrame(rowsD)
_d = TD.set_index(['end', 'source'])
assert abs(_d.loc[('5′ end', 'long-read'), 'trimmed'] - 0.483) < 0.002 and abs(_d.loc[('5′ end', 'model'), 'agrees'] - 0.426) < 0.002

# ---- E: UTR length distributions (refined: all genes with a refined end; v71: genes with a UTR)
UL = {('5′', 'Refined'): G.utr5_final.dropna(), ('5′', 'v71'): G.utr5_annot[G.utr5_annot > 0],
      ('3′', 'Refined'): G.utr3_final.dropna(), ('3′', 'v71'): G.utr3_annot[G.utr3_annot > 0]}
TE = pd.DataFrame([dict(utr=u, set=ALAB[k], n=len(v), median=v.median(), q25=v.quantile(.25), q75=v.quantile(.75), n_zero=int((v <= 0).sum()))
                   for (u, k), v in UL.items()])
p5 = stats.mannwhitneyu(UL[('5′', 'Refined')], UL[('5′', 'v71')]).pvalue; p3 = stats.mannwhitneyu(UL[('3′', 'Refined')], UL[('3′', 'v71')]).pvalue
# ================================================================== end of verbatim block

# ---- write per-panel source tables
OUT.mkdir(exist_ok=True)
FF = '%.10g'
def w(df, name): df.to_csv(OUT / f'Figure4_{name}.tsv', sep='\t', index=False, float_format=FF); return df
w(TA, 'A_tss_consensus_acc200')
w(TB, 'B_tts_siegel_acc200')
w(TC, 'C_refined_vs_v71_signal')
w(pd.DataFrame({'offset_bp': OFF, **{f'{ds}|{ALAB[k]}': PROF[ds][k] for ds in DSETS for k in ANC}}), 'C_profiles')
w(pd.DataFrame([dict(anchor=ALAB[k], end=e, n_genes=NSET[(k, e)]) for k in ANC for e in ('TSS', 'TTS')]), 'C_anchor_gene_counts')
w(TD, 'D_change_vs_v71')
TE['mannwhitney_p_refined_vs_v71'] = TE.utr.map({'5′': p5, '3′': p3})
w(TE, 'E_utr_length')
w(pd.concat([pd.DataFrame(dict(utr=u, set=ALAB[k], gene_id=G.loc[v.index, 'gene_id'].values, length_bp=v.values)) for (u, k), v in UL.items()]),
  'E_utr_lengths_per_gene')
w(DSF, 'F_accuracy_by_expression')

# ---- check against the tables saved with the original composite (only if they are present in RESULTS_DIR)
if all(os.path.exists(result(v)) for v in APPROVED.values()):
    for nm, vid in APPROVED.items():
        a = pd.read_csv(result(vid), sep='\t'); b = pd.read_csv(OUT / f'Figure4_{nm}.tsv', sep='\t')[a.columns]
        num = a.select_dtypes('number').columns
        assert a.shape == b.shape and (a.drop(columns=num).astype(str).values == b.drop(columns=num).astype(str).values).all(), nm
        assert np.allclose(a[num].values.astype(float), b[num].values.astype(float), rtol=1e-4, atol=1e-9, equal_nan=True), nm
    print('all tables match the approved Figure 4 source tables')
print('wrote', sorted(p.name for p in OUT.glob('Figure4_*.tsv')))
