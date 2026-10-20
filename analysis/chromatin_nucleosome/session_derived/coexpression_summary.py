"""Session-derived computation of coexpression_summary.tsv.

Recovered from the record of the interactive analysis session in which coexpression_summary.tsv was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: (none recorded)
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import LAB_DIR, result  # paths: config/paths.yaml
import pandas as pd, numpy as np, pyBigWig, os, scipy, sys
from scipy.stats import mannwhitneyu, spearmanr, binomtest, wilcoxon
from scipy import ndimage

g = pd.read_csv(result('gene_boundaries_v3_all_genes.tsv'), sep='\t')

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
    cs = np.concatenate([[0], np.cumsum(a, dtype=np.float64)]); h = w//2
    i = np.arange(len(a)); lo = np.clip(i-h, 0, len(a)); hi = np.clip(i+h+1, 0, len(a))
    return ((cs[hi]-cs[lo])/(hi-lo)).astype(np.float32)

nucmask = np.zeros(GL, bool)
for c in nuc: nucmask[off[c]:off[c]+len(seq[c])] = True

g['sgn'] = np.where(g.strand == '+', 1, -1)
ref = g.lr_group == 'confirmed'
g['ref_TSS'] = np.where(ref, g.lr_TSS, g.chromfree_model_TSS)
g['anc_src'] = np.where(ref, 'long-read', np.where(g.chromfree_model_TSS.notna(), 'chromfree_model', 'none'))
g['v71_has_utr'] = g.utr5_annot > 0
g['utr5_ref'] = (g.atg - g.ref_TSS) * g.sgn
g['utr5_v71'] = (g.atg - g.annot_tss) * g.sgn

S = g[(g.anc_src != 'none') & g.v71_has_utr].reset_index(drop=True)
sg = S.sgn.values; so = S.chrom.map(off).values
anch = {'refined': (so + S.ref_TSS.values).astype(int), 'v71': (so + S.annot_tss.values).astype(int), 'ATG': (so + S.atg.values).astype(int)}
utr = {'refined': S.utr5_ref.values, 'v71': S.utr5_v71.values, 'ATG': np.zeros(len(S))}
rel = np.arange(-1000, 1001)

def peakstats(prof):
    p = smooth(prof.astype(np.float32), 51)
    base = np.r_[p[:200], p[-200:]].mean()
    i = int(np.argmax(p)); h = p[i] - base
    half = base + h/2; l = i
    while l > 0 and p[l] > half: l -= 1
    r = i
    while r < len(p)-1 and p[r] > half: r += 1
    return rel[i], h, p[i]/base if base > 0 else np.nan, r-l

rng = np.random.default_rng(1)

TR = f"{LAB_DIR}/lncRNA project tracks/"
files = sorted(os.listdir(TR+'apiap2_tf'))

def gather(arr, anchor_gpos, sgn, flank=1000):
    rel_arr = np.arange(-flank, flank+1)
    idx = anchor_gpos[:, None] + sgn[:, None] * rel_arr[None, :]
    idx = np.clip(idx, 0, GL-1)
    return arr[idx]

def boundary(prof):
    p = smooth(prof.astype(np.float32), 51)
    up = p[rel <= -100]; i_max = int(np.argmax(np.where(rel <= 200, p, -np.inf)))
    top = p[i_max]; low = p[(rel >= 300) & (rel <= 1000)].mean()
    half = low + (top-low)/2; j = i_max
    while j < len(p)-1 and p[j] > half: j += 1
    hi90 = low + 0.9*(top-low); lo10 = low + 0.1*(top-low)
    a_ = i_max
    while a_ < len(p)-1 and p[a_] > hi90: a_ += 1
    b_ = a_
    while b_ < len(p)-1 and p[b_] > lo10: b_ += 1
    drop = p[(rel >= -200) & (rel < 0)].mean()/max(p[(rel > 0) & (rel <= 200)].mean(), 1e-6)
    return rel[j], b_-a_, drop

relw = np.arange(-2500, 0)
rows = []; profs = {}; summits = {}
tracks = [('apiap2_tf/'+f, f.split('_')[0], f.split('_')[1]) for f in files] + \
         [('accessibility_nucleosome/ATAC_%s_minusInput.bw' % s, 'ATAC', s) for s in ['Ring', 'Trophozoite', 'Schizont']] + \
         [('accessibility_nucleosome/MNase_%s_CPM.bw' % s, 'MNase', s) for s in ['Ring', 'Trophozoite', 'Schizont']]
lr = (S.anc_src == 'long-read').values
for path, fac, stage in tracks:
    a = load(TR+path); sm = smooth(a, 101); thr = np.percentile(sm[nucmask], 99)
    idx = np.clip(anch['ATG'][:, None]+sg[:, None]*relw[None, :], 0, GL-1); M = sm[idx]
    am_ = M.argmax(1); summ = relw[am_]; edge = (am_ < 50) | (am_ > len(relw)-50)
    bound = (M.max(1) > thr) & ~edge
    summits[(fac, stage)] = (bound.copy(), summ.copy())
    Xr = gather(a, anch['refined'], sg); pr = smooth(Xr.mean(0).astype(np.float32), 51)
    cls_ratio = pr[(rel >= -1000) & (rel <= -100)].mean()/max(pr[(rel >= 100) & (rel <= 1000)].mean(), 1e-6)
    for an in ['refined', 'v71', 'ATG']:
        X = Xr if an == 'refined' else gather(a, anch[an], sg)
        sets = [('all', np.ones(len(S), bool)), ('bound', bound)]
        if an == 'refined': sets += [('all_longread_anchor', lr), ('all_chromfree_anchor', ~lr)]
        for gs, mask in sets:
            if mask.sum() < 10: continue
            Xm = X[mask]; prof = Xm.mean(0); profs[(fac, stage, an, gs)] = prof
            pos, h, fold, fw = peakstats(prof); bpos, bw_, drop = boundary(prof)
            bs = []
            for _ in range(200):
                bi = rng.integers(0, len(Xm), len(Xm)); pb = Xm[bi].mean(0)
                bs.append(peakstats(pb)[:1]+boundary(pb))
            bs = np.array(bs, float)
            ci = lambda k: (np.percentile(bs[:, k], 2.5), np.percentile(bs[:, k], 97.5))
            d = dict(factor=fac, stage=stage, profile_class_ratio_up_vs_down_refined=float(cls_ratio), anchor=an, gene_set=gs, n_genes=int(mask.sum()), n_bound=int(bound.sum()),
                     agg_peak_pos=int(pos), agg_peak_pos_ci=f'{ci(0)[0]:.0f};{ci(0)[1]:.0f}', agg_peak_height_over_flank=float(h), agg_peak_fold_over_flank=float(fold), agg_fwhm_bp=int(fw),
                     signal_at_anchor=float(smooth(prof.astype(np.float32), 51)[1000]),
                     downstream_halfmax_edge=int(bpos), edge_ci=f'{ci(1)[0]:.0f};{ci(1)[1]:.0f}', edge_width_90to10=int(bw_), width_ci=f'{ci(2)[0]:.0f};{ci(2)[1]:.0f}',
                     drop_ratio_up200_vs_down200=float(drop), drop_ci=f'{ci(3)[0]:.2f};{ci(3)[1]:.2f}')
            if gs == 'bound':
                off_ = summ[mask]+utr[an][mask]
                d.update(summit_offset_median=float(np.median(off_)), summit_offset_q25=float(np.percentile(off_, 25)), summit_offset_q75=float(np.percentile(off_, 75)),
                         frac_summit_in_m600_0=float(np.mean((off_ >= -600) & (off_ <= 0))), frac_summit_downstream_of_anchor=float(np.mean(off_ > 0)))
                if an != 'ATG':
                    u = utr[an][mask]; s_ = summ[mask]
                    obs = np.mean((off_ >= -600) & (off_ <= 0))
                    null = np.array([np.mean(((s_+(pu := rng.permutation(u))) >= -600) & ((s_+pu) <= 0)) for _ in range(1000)])
                    sl = np.polyfit(u, s_, 1)[0]; bsl = [np.polyfit(u[bi], s_[bi], 1)[0] for bi in [rng.integers(0, len(u), len(u)) for _ in range(500)]]
                    d.update(perm_null_frac_in_m600_0=float(null.mean()), perm_p_window=float((1+(null >= obs).sum())/1001),
                             slope_summitATG_on_utr5=float(sl), slope_ci=f'{np.percentile(bsl,2.5):.2f};{np.percentile(bsl,97.5):.2f}',
                             spearman_summitATG_vs_utr5=float(spearmanr(u, s_)[0]))
            rows.append(d)
    print(fac, stage, int(bound.sum()), round(cls_ratio, 2), flush=True)
ap = pd.DataFrame(rows)

def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p); r = np.empty(n); r[o] = np.minimum.accumulate((p[o]*n/np.arange(1, n+1))[::-1])[::-1]; return np.minimum(r, 1)

from scipy.stats import pearsonr
mc = []
for (fac, stage), (bound, summ) in summits.items():
    s_ = summ[bound]; dr = (s_+utr['refined'][bound]) > 0; dv = (s_+utr['v71'][bound]) > 0
    b10 = int((dr & ~dv).sum()); b01 = int((~dr & dv).sum())
    p = binomtest(b10, b10+b01, 0.5).pvalue if b10+b01 > 0 else 1
    mc.append(dict(factor=fac, stage=stage, n_bound=int(bound.sum()), frac_summit_in_5UTR_refined=dr.mean(), frac_summit_in_5UTR_v71=dv.mean(), discord_ref_only=b10, discord_v71_only=b01, mcnemar_p=p))
mc = pd.DataFrame(mc)
mc['mcnemar_q_BH'] = bh(mc.mcnemar_p)

def binned(a, w=100):
    n = (GL//w)*w; return a[:n].reshape(-1, w).mean(1)
nb_ = binned(nucmask.astype(np.float32)) > 0.99
atac = {s: binned(load(TR+'accessibility_nucleosome/ATAC_%s_minusInput.bw' % s)) for s in ['Ring', 'Trophozoite', 'Schizont']}
cor = {}
for f in files:
    fac, stage = f.split('_')[0], f.split('_')[1]
    v = binned(load(TR+'apiap2_tf/'+f)); cor[(fac, stage)] = pearsonr(np.log1p(np.clip(v[nb_], 0, None)), np.log1p(np.clip(atac[stage][nb_], 0, None)))[0]
mc['pearson_log_signal_vs_ATAC_sameStage_100bp'] = [cor.get((r.factor, r.stage), np.nan) for r in mc.itertuples()]

gff = pd.read_csv(f'{LAB_DIR}/PlasmoDB/PlasmoDB-68/PlasmoDB-68_Pfalciparum3D7.gff', sep='\t', comment='#', header=None,
                  names=['chrom', 'src', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attr'], low_memory=False)
genes_all = gff[gff.type.str.endswith('gene')].copy()
genes_all['gid'] = genes_all.attr.str.extract(r'ID=([^;]+)')
genes_all = genes_all[genes_all.chrom.isin(nuc)]
gi = g.set_index('gene_id')
pairs = []
for c, sub in genes_all.sort_values(['chrom', 'start']).groupby('chrom'):
    sub = sub.reset_index(drop=True)
    for i in range(len(sub)-1):
        a, b = sub.iloc[i], sub.iloc[i+1]
        if a.gid in gi.index and b.gid in gi.index:
            ori = {'-+': 'divergent', '++': 'tandem', '--': 'tandem', '+-': 'convergent'}[a.strand+b.strand]
            pairs.append((a.gid, b.gid, c, ori, b.start-a.end))
P = pd.DataFrame(pairs, columns=['gene_left', 'gene_right', 'chrom', 'orientation', 'gene_gap_bp'])

ex = pd.read_csv(result('daft_gene_expression_v68.tsv.gz'), sep='\t').set_index('gene_id')
tp = ['h0', 'h8', 'h16', 'h24', 'h32', 'h40', 'h48']
LX = np.log2(ex[tp]+1); dmax = ex[tp].max(1)
gi = g.set_index('gene_id')
def col(ids, c): return gi.loc[ids, c].values
L_ = P.gene_left.values; R_ = P.gene_right.values
for nm, c in [('refinedRule', 'ref_TSS'), ('final', 'final_TSS'), ('v71', 'annot_tss')]:
    P[f'TSS_left_{nm}'] = col(L_, c); P[f'TSS_right_{nm}'] = col(R_, c)
P['atg_left'] = col(L_, 'atg'); P['atg_right'] = col(R_, 'atg'); P['stop_left'] = col(L_, 'stop'); P['stop_right'] = col(R_, 'stop')
P['strand_left'] = col(L_, 'strand'); P['strand_right'] = col(R_, 'strand')
P['src_left'] = col(L_, 'anc_src'); P['src_right'] = col(R_, 'anc_src')
P['v71utr_left'] = col(L_, 'v71_has_utr'); P['v71utr_right'] = col(R_, 'v71_has_utr')
dv = P.orientation == 'divergent'
for nm in ['refinedRule', 'final', 'v71']:
    P[f'TSS_dist_{nm}'] = np.where(dv, P[f'TSS_right_{nm}']-P[f'TSS_left_{nm}'], np.nan)
P.loc[dv & ~(P.v71utr_left & P.v71utr_right), 'TSS_dist_v71'] = np.nan
P['ATG_dist'] = np.where(dv, P.atg_right-P.atg_left, np.nan)
def cds_bounds(a, s, st): return (np.minimum(a, s), np.maximum(a, s))
lo_l, hi_l = cds_bounds(P.atg_left, P.stop_left, None); lo_r, hi_r = cds_bounds(P.atg_right, P.stop_right, None)
P['CDS_gap'] = lo_r-hi_l
r_list = []
for a, b in zip(L_, R_):
    if a in LX.index and b in LX.index: r_list.append(np.corrcoef(LX.loc[a], LX.loc[b])[0, 1])
    else: r_list.append(np.nan)
P['pearson_DAFT_log2'] = r_list
P['DAFT_max_left'] = dmax.reindex(L_).values; P['DAFT_max_right'] = dmax.reindex(R_).values
P['min_DAFT_max'] = P[['DAFT_max_left', 'DAFT_max_right']].min(1)

rng2 = np.random.default_rng(7)
ids = np.array([i for i in g.gene_id if i in LX.index])
chr_ = gi.loc[ids, 'chrom'].values
ra = rng2.integers(0, len(ids), 20000); rb = rng2.integers(0, len(ids), 20000)
keep = chr_[ra] != chr_[rb]; ra, rb = ra[keep][:10000], rb[keep][:10000]
X = LX.loc[ids].values; Xc = X-X.mean(1, keepdims=True); nrm = np.linalg.norm(Xc, axis=1)
rr = (Xc[ra]*Xc[rb]).sum(1)/(nrm[ra]*nrm[rb]); dm = dmax.loc[ids].values
R = pd.DataFrame(dict(gene_left=ids[ra], gene_right=ids[rb], orientation='random_interchrom', pearson_DAFT_log2=rr, min_DAFT_max=np.minimum(dm[ra], dm[rb])))
ALL = pd.concat([P[['gene_left', 'gene_right', 'orientation', 'pearson_DAFT_log2', 'min_DAFT_max']], R])

def fz(x): return np.arctanh(np.clip(x, -0.999, 0.999))
def boot_mean(x, n=2000):
    x = np.asarray(x); b = [np.mean(x[rng.integers(0, len(x), len(x))]) for _ in range(n)]; return np.mean(x), np.percentile(b, 2.5), np.percentile(b, 97.5)

D = P[dv & (P.min_DAFT_max >= 1) & P.pearson_DAFT_log2.notna()].copy()
both = D.TSS_dist_refinedRule.notna() & D.TSS_dist_v71.notna()
edges = [-1e9, 0, 250, 500, 1000, 2000, 1e9]; labs = ['<=0', '1-250', '251-500', '501-1000', '1001-2000', '>2000']
BIN_rows = []
for nm in ['TSS_dist_refinedRule', 'TSS_dist_v71', 'ATG_dist']:
    Dd = D[both]
    bb_ = pd.cut(Dd[nm], edges, labels=labs)
    for lab, s in Dd.groupby(bb_, observed=False):
        if len(s) < 5: BIN_rows.append(dict(distance=nm, bin=lab, n=len(s))); continue
        m, lo, hi = boot_mean(s.pearson_DAFT_log2, 1000); BIN_rows.append(dict(distance=nm, bin=lab, n=len(s), mean_r=m, ci_lo=lo, ci_hi=hi, median_r=s.pearson_DAFT_log2.median()))
BIN = pd.DataFrame(BIN_rows)

T = P[(P.orientation == 'tandem') & (P.min_DAFT_max >= 1) & P.pearson_DAFT_log2.notna()]

acc = {}
for st in ['Ring', 'Trophozoite', 'Schizont']:
    sm = smooth(load(TR+'accessibility_nucleosome/ATAC_%s_minusInput.bw' % st), 101)
    thr = np.percentile(sm[nucmask], 90)
    lab, _ = ndimage.label(sm > thr); acc[st] = (sm, lab, thr)
mn = {st: smooth(load(TR+'accessibility_nucleosome/MNase_%s_CPM.bw' % st), 51) for st in ['Ring', 'Trophozoite', 'Schizont']}

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

E = Dv[(Dv.min_DAFT_max >= 1) & Dv.pearson_DAFT_log2.notna() & Dv.TSS_dist_refinedRule.notna()].copy()

os.makedirs('out', exist_ok=True)

cs = []
for thr in [1, 5]:
    sub = ALL[(ALL.min_DAFT_max >= thr) & ALL.pearson_DAFT_log2.notna()]
    for o, s in sub.groupby('orientation'):
        m, lo, hi = boot_mean(s.pearson_DAFT_log2); cs.append(dict(table='pair_class', expr_filter_min_DAFTmax=thr, group=o, n=len(s), mean_r=m, ci_lo=lo, ci_hi=hi, median_r=s.pearson_DAFT_log2.median()))
for r in BIN.itertuples(): cs.append(dict(table='divergent_by_'+r.distance, expr_filter_min_DAFTmax=1, group=r.bin, n=r.n, mean_r=r.mean_r, ci_lo=r.ci_lo, ci_hi=r.ci_hi, median_r=r.median_r))
for lab_, s in T.groupby(pd.cut(T.CDS_gap, [-1e9, 500, 1000, 2000, 4000, 1e9]), observed=False):
    m, lo, hi = boot_mean(s.pearson_DAFT_log2, 1000); cs.append(dict(table='tandem_by_CDS_gap', expr_filter_min_DAFTmax=1, group=str(lab_), n=len(s), mean_r=m, ci_lo=lo, ci_hi=hi, median_r=s.pearson_DAFT_log2.median()))
for lab_, s in E.groupby([pd.cut(E.TSS_dist_refinedRule, edges, labels=labs), 'NDR_refinedRule_anyStage'], observed=True):
    if len(s) >= 5:
        m, lo, hi = boot_mean(s.pearson_DAFT_log2, 1000); cs.append(dict(table='divergent_by_refinedTSSdist_x_NDR', expr_filter_min_DAFTmax=1, group=f'{lab_[0]}|{lab_[1]}', n=len(s), mean_r=m, ci_lo=lo, ci_hi=hi, median_r=s.pearson_DAFT_log2.median()))
pd.DataFrame(cs).to_csv('coexpression_summary.tsv', sep='\t', index=False, float_format='%.4g')
