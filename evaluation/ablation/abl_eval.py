"""Ablation evaluation on held-out chr 1-5. For each 3-model ensemble:
 (1) peak_eval.report_fixed against long-read ends (calls fixed once; thr / NMS tuned on chr 6-7);
 (2) with the same calls: per-gene main TSS / TTS (highest-probability call in the gene's search window) -> abl_gene_calls.tsv
     (scored locally against the TSS consensus set and Siegel 2014 poly(A) sites);
 (3) gene-free long-read lncRNA 5' ends with a TSS call <=200 bp (same strand), and strand mirroring
     (long-read TSSs with a TSS call <=200 bp on the OPPOSITE strand)."""
import sys, os, numpy as np, pandas as pd
sys.path[:0] = ['.', 'boundary_model']
import peak_eval as PE
RES = 20
W = pd.read_csv('hpcc/gene_windows_truth.tsv', sep='\t'); W = W[W.chrom.str[6:8].astype(int) <= 5]
LN = pd.read_csv('hpcc/genefree_lncRNA_tss.tsv', sep='\t'); LN = LN[LN.chrom.str[6:8].astype(int) <= 5]
cfg = [l.split('\t')[0] for l in open('hpcc/configs.tsv').read().strip().split('\n')]
ENS = [(n, f'{n}+{n}_s1+{n}_drop') for n in cfg] + [('V3_top8_original', 'V3_top8+V3_top8_s1+V3_top8_drop'), ('V3_nochrom_original', 'V3_nochrom+V3_nochrom_s1+V3_nochrom_drop')]
done = lambda nm: all(os.path.exists(f'exp_ml/{k}_{x}.npz') for x in nm.split('+') for k in ('pred', 'logit'))
pts, genes, extra = [], [], []
for lab, nm in ENS:
    if not done(nm): print('missing', lab, flush=True); continue
    rows = PE.report_fixed(lab, nm, 'logit_nms'); pts += rows
    PR = PE.load(nm, 'pred'); LG = PE.load(nm, 'logit')
    par = {r['end_type']: (r['thr'], r['nms_bp']) for r in rows if r['scheme'] == 'one_to_one'}
    C = {}
    for ch, et in ((0, 'TSS'), (1, 'TTS')):
        t, d = par[et]
        for c in PR:
            if int(c[6:8]) > 5: continue
            for k in (0, 1):
                idx = PE.call(PR[c][k, ch], LG[c][k, ch][:PR[c].shape[-1]], 'logit_nms', t, d); C[(c, k, ch)] = (idx * RES + RES // 2, PR[c][k, ch, idx])
    for r in W.itertuples():
        k = 0 if r.strand == '+' else 1; d = dict(model=lab, gene_id=r.gene_id)
        for ch, et, lo, hi in ((0, 'TSS', r.tss_win_lo, r.tss_win_hi), (1, 'TTS', r.tts_win_lo, r.tts_win_hi)):
            p, s = C[(r.chrom, k, ch)]; m = (p >= lo) & (p <= hi)
            d[et] = int(p[m][np.argmax(s[m])]) if m.any() else np.nan
        genes.append(d)
    def near(c, k, x, tol=200):
        p = C[(c, k, 0)][0]; return bool(len(p) and np.abs(p - x).min() <= tol)
    lnc = np.mean([near(r.chrom, 0 if r.strand == '+' else 1, r.tss) for r in LN.itertuples()])
    q = W[W.lr_TSS.notna()]
    mir = np.mean([near(r.chrom, 1 if r.strand == '+' else 0, r.lr_TSS) for r in q.itertuples()])
    sense = np.mean([near(r.chrom, 0 if r.strand == '+' else 1, r.lr_TSS) for r in q.itertuples()])
    extra.append(dict(model=lab, lncRNA_TSS_called=lnc, n_lncRNA=len(LN), mirror_TSS_opposite=mir, sense_TSS_called=sense, n_lr=len(q),
                      thr_TSS=par['TSS'][0], nms_TSS=par['TSS'][1], thr_TTS=par['TTS'][0], nms_TTS=par['TTS'][1]))
    print('evaluated', lab, flush=True)
pd.DataFrame(pts).to_csv('abl_point_metrics.tsv', sep='\t', index=False); pd.DataFrame(genes).to_csv('abl_gene_calls.tsv', sep='\t', index=False)
pd.DataFrame(extra).to_csv('abl_extra.tsv', sep='\t', index=False); print('EVAL_DONE')
