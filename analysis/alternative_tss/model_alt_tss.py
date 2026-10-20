"""Validation of model-predicted (held-out) TSS calls with three independent 5'-end datasets,
and gene-level alternative-TSS analysis based only on validated model calls."""
from pfconfig import LAB_DIR  # paths: config/paths.yaml
import numpy as np, pandas as pd, pyBigWig, os

CAGE_DIR = f"{LAB_DIR}/lncRNA project tracks/cage_tss"
STAGES = ["ring", "trophozoite", "schizont"]
W = 50            # half-window for support around a call (bp)
FL_IN, FL_OUT = 100, 500   # flank for local background
MIN_SEP = 150     # alt call must be >= this from the main call

def load_calls(gb):
    rows = []
    for r in gb.itertuples():
        if pd.isna(r.model_TSS):
            continue
        rows.append((r.gene_id, r.chrom, r.strand, int(r.model_TSS), "main"))
        if isinstance(r.model_alt_TSS, str):
            for p in r.model_alt_TSS.split(","):
                p = int(float(p))
                if abs(p - int(r.model_TSS)) >= MIN_SEP:
                    rows.append((r.gene_id, r.chrom, r.strand, p, "alt"))
    C = pd.DataFrame(rows, columns=["gene_id", "chrom", "strand", "pos", "role"])
    a = gb.set_index("gene_id")
    C["atg"] = C.gene_id.map(a.atg); C["stop"] = C.gene_id.map(a.stop)
    sgn = np.where(C.strand == "+", 1, -1)
    C["rel_atg"] = (C.pos - C.atg) * sgn          # negative = upstream of ATG
    cds_len = (C.stop - C.atg) * sgn
    C["cls"] = np.where(C.rel_atg < 0, "upstream", np.where(C.rel_atg <= cds_len, "CDS", "downstream"))
    return C

def make_decoys(C, rng):
    """Shift each call by +-(300..1000) bp, keeping the same class relative to ATG and
    >= MIN_SEP from every model call of the same gene."""
    calls_by_gene = C.groupby("gene_id").pos.apply(np.array).to_dict()
    out = []
    for r in C.itertuples():
        for _ in range(50):
            d = int(rng.integers(300, 1001)) * (1 if rng.random() < 0.5 else -1)
            p = r.pos + d
            sg = 1 if r.strand == "+" else -1
            rel = (p - r.atg) * sg; cl = (r.stop - r.atg) * sg
            cls = "upstream" if rel < 0 else ("CDS" if rel <= cl else "downstream")
            if cls != r.cls or np.min(np.abs(calls_by_gene[r.gene_id] - p)) < MIN_SEP:
                continue
            out.append((r.gene_id, r.chrom, r.strand, p, "decoy_" + r.role, r.atg, r.stop, rel, cls)); break
    return pd.DataFrame(out, columns=list(C.columns))

def cage_metrics(C):
    bws = {(s, st): pyBigWig.open(f"{CAGE_DIR}/CTSS_{s}_{'pos' if st=='+' else 'neg'}_CPM.bw") for s in STAGES for st in "+-"}
    def win_sum(bw, ch, a, b):
        a = max(0, a); b = min(bw.chroms(ch), b)
        if b <= a: return 0.0
        iv = bw.intervals(ch, a, b) or []
        return float(sum((min(e, b) - max(s0, a)) * v for s0, e, v in iv))
    res = {f"cage_{s}": [] for s in STAGES}; res["cage_flank_perbp"] = []
    for r in C.itertuples():
        tot_fl = 0.0
        for s in STAGES:
            bw = bws[(s, r.strand)]
            res[f"cage_{s}"].append(win_sum(bw, r.chrom, r.pos - W, r.pos + W + 1))
            tot_fl += win_sum(bw, r.chrom, r.pos - FL_OUT, r.pos - FL_IN) + win_sum(bw, r.chrom, r.pos + FL_IN + 1, r.pos + FL_OUT + 1)
        res["cage_flank_perbp"].append(tot_fl / (2 * (FL_OUT - FL_IN)))
    for k, v in res.items(): C[k] = v
    C["cage_win"] = C[[f"cage_{s}" for s in STAGES]].sum(axis=1)
    C["cage_enrich"] = (C.cage_win / (2 * W + 1)) / (C.cage_flank_perbp + 0.01)
    return C

def cluster_hits(C, clusters, pos_col, name, score_col):
    """Nearest same-strand cluster (dominant position) within +-W bp."""
    C[f"{name}_hit"] = False; C[f"{name}_score"] = 0.0
    for (ch, st), sub in C.groupby(["chrom", "strand"]):
        cl = clusters[(clusters.chrom == ch) & (clusters.strand == st)]
        if cl.empty: continue
        cp = np.sort(cl[pos_col].values); order = np.argsort(cl[pos_col].values); sc = cl[score_col].values[order]
        i = np.searchsorted(cp, sub.pos.values)
        best = np.zeros(len(sub)); hit = np.zeros(len(sub), bool)
        for off in (-1, 0):
            j = np.clip(i + off, 0, len(cp) - 1)
            ok = np.abs(cp[j] - sub.pos.values) <= W
            hit |= ok; best = np.where(ok, np.maximum(best, sc[j]), best)
        # also scan neighbours for multiple clusters in window (take max score)
        C.loc[sub.index, f"{name}_hit"] = hit; C.loc[sub.index, f"{name}_score"] = best
    return C
