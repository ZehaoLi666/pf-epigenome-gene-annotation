"""Session-derived computation of bidirectional_pairs.tsv.

Recovered from the record of the interactive analysis session in which bidirectional_pairs.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: (none recorded)
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import pandas as pd
import numpy as np
import pyBigWig
import os
from scipy import ndimage

# Load gene boundaries
g = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep='\t')

# Load genome
seq = {}
name = None; buf = []
for l in open(f'{LAB_DIR}/IGV_longread_transcript_models/PlasmoDB-68_Pfalciparum3D7_Genome.fasta'):
    if l.startswith('>'):
        if name: seq[name] = ''.join(buf).upper()
        name = l[1:].split()[0].split('|')[0].strip(); buf = []
    else: buf.append(l.strip())
seq[name] = ''.join(buf).upper()

comp = str.maketrans('ACGTN', 'TGCAN')
rc = lambda s: s.translate(comp)[::-1]

g['sgn'] = np.where(g.strand == '+', 1, -1)
ref = g.lr_group == 'confirmed'
g['ref_TSS'] = np.where(ref, g.lr_TSS, g.chromfree_model_TSS)
g['anc_src'] = np.where(ref, 'long-read', np.where(g.chromfree_model_TSS.notna(), 'chromfree_model', 'none'))
g['v71_has_utr'] = g.utr5_annot > 0
g['utr5_ref'] = (g.atg - g.ref_TSS) * g.sgn
g['utr5_v71'] = (g.atg - g.annot_tss) * g.sgn

TR = f"{LAB_DIR}/lncRNA project tracks/"
nuc = [c for c in seq if c in g.chrom.unique()]
off = {}; o = 0
for c in seq: off[c] = o; o += len(seq[c])
GL = o

def load(path):
    b = pyBigWig.open(path); arr = np.zeros(GL, np.float32)
    for c in seq:
        if c in b.chroms():
            v = np.array(b.values(c, 0, b.chroms()[c]), dtype=np.float32); v[np.isnan(v)] = 0
            arr[off[c]:off[c]+len(v)] = v
    b.close(); return arr

def smooth(a, w=101):
    cs = np.concatenate([[0], np.cumsum(a, dtype=np.float64)]); h = w // 2
    i = np.arange(len(a)); lo = np.clip(i - h, 0, len(a)); hi = np.clip(i + h + 1, 0, len(a))
    return ((cs[hi] - cs[lo]) / (hi - lo)).astype(np.float32)

nucmask = np.zeros(GL, bool)
for c in nuc: nucmask[off[c]:off[c]+len(seq[c])] = True

# Load GFF
gff = pd.read_csv(f'{LAB_DIR}/PlasmoDB/PlasmoDB-68/PlasmoDB-68_Pfalciparum3D7.gff', sep='\t', comment='#', header=None,
                  names=['chrom', 'src', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attr'], low_memory=False)
genes_all = gff[gff.type.str.endswith('gene')].copy()
genes_all['gid'] = genes_all.attr.str.extract(r'ID=([^;]+)')
genes_all = genes_all[genes_all.chrom.isin(nuc)]

gi = g.set_index('gene_id')
pairs = []
for c, sub in genes_all.sort_values(['chrom', 'start']).groupby('chrom'):
    sub = sub.reset_index(drop=True)
    for i in range(len(sub) - 1):
        a, b = sub.iloc[i], sub.iloc[i+1]
        if a.gid in gi.index and b.gid in gi.index:
            ori = {'-+': 'divergent', '++': 'tandem', '--': 'tandem', '+-': 'convergent'}[a.strand + b.strand]
            pairs.append((a.gid, b.gid, c, ori, b.start - a.end))
P = pd.DataFrame(pairs, columns=['gene_left', 'gene_right', 'chrom', 'orientation', 'gene_gap_bp'])

# Load expression
ex = pd.read_csv(result('daft_gene_expression_v68.tsv.gz'), sep='\t').set_index('gene_id')
tp = ['h0', 'h8', 'h16', 'h24', 'h32', 'h40', 'h48']
LX = np.log2(ex[tp] + 1); dmax = ex[tp].max(1)

gi = g.set_index('gene_id')
def col(ids, c): return gi.loc[ids, c].values
L_, R_ = P.gene_left.values, P.gene_right.values
for nm, c in [('refinedRule', 'ref_TSS'), ('final', 'final_TSS'), ('v71', 'annot_tss')]:
    P[f'TSS_left_{nm}'] = col(L_, c); P[f'TSS_right_{nm}'] = col(R_, c)
P['atg_left'] = col(L_, 'atg'); P['atg_right'] = col(R_, 'atg')
P['stop_left'] = col(L_, 'stop'); P['stop_right'] = col(R_, 'stop')
P['strand_left'] = col(L_, 'strand'); P['strand_right'] = col(R_, 'strand')
P['src_left'] = col(L_, 'anc_src'); P['src_right'] = col(R_, 'anc_src')
P['v71utr_left'] = col(L_, 'v71_has_utr'); P['v71utr_right'] = col(R_, 'v71_has_utr')
dv = P.orientation == 'divergent'
for nm in ['refinedRule', 'final', 'v71']:
    P[f'TSS_dist_{nm}'] = np.where(dv, P[f'TSS_right_{nm}'] - P[f'TSS_left_{nm}'], np.nan)
P.loc[dv & ~(P.v71utr_left & P.v71utr_right), 'TSS_dist_v71'] = np.nan
P['ATG_dist'] = np.where(dv, P.atg_right - P.atg_left, np.nan)
def cds_bounds(a, s, st): return (np.minimum(a, s), np.maximum(a, s))
lo_l, hi_l = cds_bounds(P.atg_left, P.stop_left, None); lo_r, hi_r = cds_bounds(P.atg_right, P.stop_right, None)
P['CDS_gap'] = lo_r - hi_l

rng = np.random.default_rng(1)
r = []
for a, b in zip(L_, R_):
    if a in LX.index and b in LX.index: r.append(np.corrcoef(LX.loc[a], LX.loc[b])[0, 1])
    else: r.append(np.nan)
P['pearson_DAFT_log2'] = r
P['DAFT_max_left'] = dmax.reindex(L_).values; P['DAFT_max_right'] = dmax.reindex(R_).values
P['min_DAFT_max'] = P[['DAFT_max_left', 'DAFT_max_right']].min(1)

# NDR classification
acc = {}
for st in ['Ring', 'Trophozoite', 'Schizont']:
    sm = smooth(load(TR + 'accessibility_nucleosome/ATAC_%s_minusInput.bw' % st), 101)
    thr = np.percentile(sm[nucmask], 90)
    lab, _ = ndimage.label(sm > thr); acc[st] = (sm, lab, thr)
mn = {st: smooth(load(TR + 'accessibility_nucleosome/MNase_%s_CPM.bw' % st), 51) for st in ['Ring', 'Trophozoite', 'Schizont']}

def ndr_class(c, tl, tr, st, W=300):
    sm, lab, thr = acc[st]; o = off[c]; tl = int(tl); tr = int(tr)
    wl = lab[o+tl:o+tl+W+1]; wr = lab[o+tr-W:o+tr+1]
    sl = set(wl[wl > 0]); sr = set(wr[wr > 0])
    if not sl and not sr: return 'none'
    if not sl or not sr: return 'one_side_only'
    return 'shared' if sl & sr else 'separate'

def mean_sig(arr, c, a, b):
    a, b = int(min(a, b)), int(max(a, b)); o = off[c]; return float(arr[o+a:o+b+1].mean()) if b > a else np.nan

Dv = P[dv].copy()
for anc in ['refinedRule', 'v71']:
    ok = Dv[f'TSS_dist_{anc}'].notna()
    for st in ['Ring', 'Trophozoite', 'Schizont']:
        Dv[f'NDR_{anc}_{st}'] = [ndr_class(r.chrom, r[f'TSS_left_{anc}'], r[f'TSS_right_{anc}'], st) if k else np.nan for (_, r), k in zip(Dv.iterrows(), ok)]
    Dv[f'ATAC_meanStages_interTSS_{anc}'] = [np.mean([mean_sig(acc[st][0], r.chrom, r[f'TSS_left_{anc}'], r[f'TSS_right_{anc}']) for st in acc]) if k and r[f'TSS_dist_{anc}'] > 50 else np.nan for (_, r), k in zip(Dv.iterrows(), ok)]
    Dv[f'MNase_meanStages_interTSS_{anc}'] = [np.mean([mean_sig(mn[st], r.chrom, r[f'TSS_left_{anc}'], r[f'TSS_right_{anc}']) for st in mn]) if k and r[f'TSS_dist_{anc}'] > 50 else np.nan for (_, r), k in zip(Dv.iterrows(), ok)]
    cols = [f'NDR_{anc}_{st}' for st in acc]
    Dv[f'NDR_{anc}_anyStage'] = Dv[cols].apply(lambda x: np.nan if x.isna().all() else ('shared' if (x == 'shared').any() else ('separate' if (x == 'separate').any() else ('one_side_only' if (x == 'one_side_only').any() else 'none'))), axis=1)

ndrcols = [c for c in Dv.columns if c.startswith('NDR_') or 'interTSS' in c]
BP = P.merge(Dv[['gene_left', 'gene_right'] + ndrcols], on=['gene_left', 'gene_right'], how='left')

os.makedirs('out', exist_ok=True)
BP.to_csv('out/bidirectional_pairs.tsv', sep='\t', index=False, float_format='%.4g')
