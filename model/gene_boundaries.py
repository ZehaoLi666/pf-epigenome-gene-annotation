"""Gene-level TSS/TTS/UTR assignment from genome-wide out-of-fold predictions of the selected model.
Usage: python gene_boundaries.py <foldA names> <foldB names> <foldC names> <thr_TSS> <thr_TTS> <nms_bp> <out_prefix>
Each fold's ensemble (names joined by '+') supplies predictions only for its held-out chromosomes (A: 1-5, B: 6-10, C: 11-14).
Per v68 protein-coding gene (CDS envelope from the GFF):
  TSS window: from the ATG up to 3 kb upstream, truncated at the nearest same-strand CDS end upstream (+50 bp).
  TTS window: from the stop codon up to 2 kb downstream, truncated at the nearest same-strand CDS start downstream (-50 bp).
  Calls: smoothed-logit local maxima with P >= thr, NMS nms_bp (same caller as the evaluation). Predicted TSS/TTS = highest-P call in the window;
  all calls in the window are listed as alternatives. 5'UTR = |ATG - TSS|, 3'UTR = |TTS - stop|."""
from pfconfig import LAB_DIR  # paths: config/paths.yaml
import sys, glob, numpy as np, pandas as pd
sys.path[:0] = ['.', 'boundary_model']
import peak_eval as PE
RES = 20; LAB = LAB_DIR
names = sys.argv[1:4]; thr = {0: float(sys.argv[4]), 1: float(sys.argv[5])}; D = dict(zip((0, 1), map(int, sys.argv[6].split(',')))); out = sys.argv[7]
rngs = [(1, 5), (6, 10), (11, 14)]; LG = {}
for nm, (lo, hi) in zip(names, rngs):
    for c, v in PE.load(nm, 'logitall').items():
        if lo <= int(c[6:8]) <= hi: LG[c] = v
assert len(LG) == 14, len(LG)
PR = {c: 1 / (1 + np.exp(-v)) for c, v in LG.items()}
np.savez_compressed(f'{out}_oof_prob.npz', **{c: v.astype(np.float16) for c, v in PR.items()})
CALL = {}
for c in PR:
    for k in (0, 1):
        for ch in (0, 1):
            idx = PE.call(PR[c][k, ch], LG[c][k, ch], 'logit_nms', thr[ch], D[ch]); CALL[(c, k, ch)] = (idx * RES + RES // 2, PR[c][k, ch, idx])
g = pd.read_csv(glob.glob(f'{LAB}/PlasmoDB/PlasmoDB-68/*.gff')[0], sep='\t', comment='#', header=None, low_memory=False, names=['chrom','src','type','start','end','score','strand','phase','attr'])
g = g[g.chrom.isin(PR)]
gene = g[g.type.isin(['protein_coding_gene'])].copy(); gene['gene_id'] = gene.attr.str.extract(r'ID=([^;]+)')[0]; gene['desc'] = gene.attr.str.extract(r'description=([^;]+)')[0]
cds = g[g.type == 'CDS'].copy(); cds['gene_id'] = cds.attr.str.extract(r'Parent=([^;.]+)')[0]
env = cds.groupby('gene_id').agg(cds_lo=('start', 'min'), cds_hi=('end', 'max')); env['cds_lo'] -= 1
G = gene[['gene_id', 'chrom', 'strand', 'start', 'end', 'desc']].merge(env, left_on='gene_id', right_index=True)
G['start'] -= 1; G['atg'] = np.where(G.strand == '+', G.cds_lo, G.cds_hi - 1); G['stop'] = np.where(G.strand == '+', G.cds_hi - 1, G.cds_lo)
G['annot_tss'] = np.where(G.strand == '+', G.start, G.end - 1); G['annot_tts'] = np.where(G.strand == '+', G.end - 1, G.start)
rows = []
for (c, s), q in G.groupby(['chrom', 'strand']):
    q = q.sort_values('cds_lo'); lo_, hi_ = q.cds_lo.values, q.cds_hi.values; k = 0 if s == '+' else 1
    for i, r in enumerate(q.itertuples()):
        prev_hi = hi_[i - 1] if i > 0 else -10**9; next_lo = lo_[i + 1] if i + 1 < len(q) else 10**12
        if s == '+': w5 = (max(r.cds_lo - 3000, prev_hi + 50), r.cds_lo + 1); w3 = (r.cds_hi - 1, min(r.cds_hi + 2000, next_lo - 50))
        else: w5 = (r.cds_hi - 1, min(r.cds_hi + 3000, next_lo - 50)); w3 = (max(r.cds_lo - 2000, prev_hi + 50), r.cds_lo + 1)
        d = dict(gene_id=r.gene_id, chrom=c, strand=s, atg=int(r.atg), stop=int(r.stop), annot_tss=int(r.annot_tss), annot_tts=int(r.annot_tts), desc=r.desc,
                 tss_win_lo=int(w5[0]), tss_win_hi=int(w5[1]), tts_win_lo=int(w3[0]), tts_win_hi=int(w3[1]))
        for ch, w, et in ((0, w5, 'TSS'), (1, w3, 'TTS')):
            p, sc = CALL[(c, k, ch)]; m = (p >= w[0]) & (p <= w[1]); p, sc = p[m], sc[m]
            d[f'n_{et}_calls'] = int(len(p))
            if len(p):
                j = int(np.argmax(sc)); d[f'pred_{et}'] = int(p[j]); d[f'pred_{et}_prob'] = float(sc[j])
                d[f'alt_{et}'] = ','.join(str(int(x)) for x in p[np.argsort(-sc)][1:]); d[f'max_{et}_prob_window'] = float(PR[c][k, ch, max(0, w[0] // RES):w[1] // RES + 1].max())
            else: d[f'max_{et}_prob_window'] = float(PR[c][k, ch, max(0, w[0] // RES):w[1] // RES + 1].max()) if w[1] > w[0] else np.nan
        rows.append(d)
R = pd.DataFrame(rows)
R['utr5_pred'] = (R.atg - R.pred_TSS).abs(); R['utr3_pred'] = (R.pred_TTS - R.stop).abs()
R['utr5_annot'] = (R.atg - R.annot_tss).abs(); R['utr3_annot'] = (R.annot_tts - R.stop).abs()
R.to_csv(f'{out}_genes.tsv.gz', sep='\t', index=False); print(len(R), R.pred_TSS.notna().mean().round(3), R.pred_TTS.notna().mean().round(3))
