from pfconfig import LAB_DIR  # paths: config/paths.yaml
import sys, os; sys.path.append(os.getcwd())
import numpy as np, pandas as pd, pyBigWig, pickle
from common import *
base = pickle.load(open('part2_base.pkl', 'rb'))['base'] if os.path.exists('part2_base.pkl') else None
G = load_genes()
if base is None:
    base = G[G.ref_TSS.notna()].reset_index(drop=True)
LT = f'{LAB_DIR}/lncRNA project tracks/'
NP = 'ext/nuc/Pf_nucleosome_annotation_of_nucDetective/nucleosome_profiles/'
tracks = [(f'nucDetective_MNase_{tp}_qNorm (Kensche2016 PE, mono-nuc fragments)', NP + f'MNase-seq_{tp}_qNorm.bw', 'nuc') for tp in ['T5', 'T10', 'T15', 'T20', 'T25', 'T30', 'T35', 'T40']]
tracks += [(f'inhouse_MNase_{s}_CPM (Bunnik2014, reads not extended)', LT + f'accessibility_nucleosome/MNase_{s}_CPM.bw', 'other') for s in ['Ring', 'Trophozoite', 'Schizont']]
tracks += [(f'inhouse_H3_{s}_minusInput', LT + f'histone_marks/H3_{s}_minusInput.bw', 'other') for s in ['Ring', 'Trophozoite', 'Schizont']]
tracks += [(f'inhouse_ATAC_{s}_minusInput', LT + f'accessibility_nucleosome/ATAC_{s}_minusInput.bw', 'other') for s in ['Ring', 'Trophozoite', 'Schizont']]
W = 1000; OFF = np.arange(-W, W + 1)
ANCH = {'refined': base.ref_TSS.values.astype(int), 'PlasmoDB_v71': base.v71_TSS.values.astype(int)}
rng = np.random.default_rng(3)

def profile_matrix(bw, X):
    M = np.full((len(base), len(OFF)), np.nan, np.float32); L = bw.chroms()
    for i, (c, s, x) in enumerate(zip(base.chrom.values, base.strand.values, X)):
        lo, hi = x - W, x + W + 1
        if lo < 0 or hi > L[c]: continue
        v = np.array(bw.values(c, lo, hi), dtype=np.float32)
        M[i] = v if s == '+' else v[::-1]
    return M

def smooth(p, w=21):
    k = np.ones(w) / w; return np.convolve(np.nan_to_num(p), k, mode='same')

def landmarks(p, kind):
    ps = smooth(p)
    def at(lo, hi, f):
        m = (OFF >= lo) & (OFF <= hi); j = np.where(m)[0][f(ps[m])]; return OFF[j], ps[j]
    if kind == 'nuc':
        p1, h1 = at(0, 300, np.argmax); tr, ht = at(-300, p1 - 60, np.argmin); m1, hm1 = at(-600, tr - 60, np.argmax)
        half = ht + (h1 - ht) / 2; j = np.where(OFF == p1)[0][0]; l = j; r = j
        while l > 0 and ps[l] > half: l -= 1
        while r < len(ps) - 1 and ps[r] > half: r += 1
        return dict(plus1_peak=p1, NDR_trough=tr, minus1_peak=m1, plus1_over_trough=h1 / ht if ht > 0 else np.nan, plus1_FWHM=OFF[r] - OFF[l], minus1_to_plus1=p1 - m1)
    mx, _ = at(-500, 500, np.argmax); mn, _ = at(-500, 500, np.argmin)
    return {'max_pos[-500,+500]': mx, 'min_pos[-500,+500]': mn}

rows = []; prof_out = {}
for name, path, kind in tracks:
    bw = pyBigWig.open(path)
    for an, X in ANCH.items():
        M = profile_matrix(bw, X)
        for grp in GROUPS:
            m = np.ones(len(base), bool) if grp == 'all' else base.lr_group.values == grp
            m &= ~np.isnan(M).all(1)
            Mm = M[m]; p = np.nanmean(Mm, 0); lm = landmarks(p, kind)
            bs = [landmarks(np.nanmean(Mm[rng.integers(0, len(Mm), len(Mm))], 0), kind) for _ in range(100 if grp == 'all' else 0)]
            r = dict(track=name, anchor=an, gene_group=grp, n_genes=int(m.sum()))
            for k, v in lm.items():
                r[k] = v
                if bs: vv = np.array([b[k] for b in bs], float); r[k + '_CI95_lo'], r[k + '_CI95_hi'] = np.nanpercentile(vv, [2.5, 97.5])
            rows.append(r)
            if grp == 'all': prof_out[(name, an)] = p
    bw.close(); print(name, flush=True)
R = pd.DataFrame(rows); R.to_csv('nucleosome_meta_profile_landmarks.tsv', sep='\t', index=False, float_format='%.4g')
P = pd.DataFrame({f'{k[0]}|{k[1]}': v for k, v in prof_out.items()}); P.insert(0, 'offset_bp', OFF)
P.to_csv('nucleosome_meta_profiles_all_genes.tsv.gz', sep='\t', index=False, float_format='%.4g')
print('done')
