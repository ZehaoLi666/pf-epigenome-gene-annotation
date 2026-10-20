"""Assemble per-(stage|chrom) input arrays for the segmentation ladder."""
import os, numpy as np, pfprep as P
from bm.tracks import RES
S = P.STAGES
def dna_onehot():
    f = "data/dna_onehot_20bp.npz"
    if os.path.exists(f): return dict(np.load(f))
    seq, name = {}, None
    for line in open(f"{P.LAB}/PlasmoDB/PlasmoDB-48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta"):
        if line.startswith(">"): name = line[1:].split()[0].split("|")[0].strip(); seq[name] = []
        else: seq[name].append(line.strip().upper())
    out = {}
    for c in P.CHROMS:
        s = np.frombuffer("".join(seq[c]).encode(), np.uint8); nb = P.CS[c] // RES; s = s[:nb * RES].reshape(nb, RES)
        out[c] = np.stack([(s == ord(b)).mean(1) for b in "ACGT"]).astype(np.float16)
    np.savez_compressed(f, **out); return out
def rna():
    m = np.load("data/mrna_20bp.npz")
    raw = {s: {c: np.log1p(m[f"{s}__{c}"].astype(np.float32)) for c in P.CHROMS} for s in S}
    Xr = {s: dict(np.load(f"data/X_rna_{s}.npz")) for s in S}  # z-scored GRO, polysome
    out = {}
    for s in S:
        tr = np.concatenate([raw[s][c] for c in P.TRAIN]); mu, sd = tr.mean(), tr.std() + 1e-6
        for c in P.CHROMS: out[f"{s}|{c}"] = np.stack([(raw[s][c] - mu) / sd, Xr[s][c][0]]).astype(np.float16)
    return out  # rows: mRNA, GRO
def glm():
    f = "data/caduceus_pca16_int8.npz"
    if not os.path.exists(f): return None
    z = np.load(f); sc = z["_scale"]; out = {}
    for c in P.CHROMS:
        a, b = z[f"{c}__fwd"].astype(np.float32) * sc[:, None], z[f"{c}__rc"].astype(np.float32) * sc[:, None]
        n = min(a.shape[1], P.CS[c] // RES); out[c] = (a[:, :n], b[:, :n])
    sd = np.concatenate([out[c][0] for c in P.TRAIN], 1).std(1, keepdims=True) + 1e-6
    return {c: np.clip(np.concatenate([a / sd, b / sd]), -8, 8).astype(np.float16) for c, (a, b) in out.items()}
def assemble_pool(spec, tracks=None):
    """POOL mode: one key per chromosome ('pool|chrom'); chromatin of ring, troph, schizont stacked (3 x 25 rows); RNA 3 x 2 rows."""
    SR = P.STAGES_REAL
    ti = list(range(25)) if tracks is None else [P.META["chrom_tracks"].index(t) for t in tracks]
    Xc = {s: dict(np.load(f"data/X_chrom_{s}.npz")) for s in SR} if "chrom" in spec else None
    R = None
    rt = (["mrna", "gro"] if "rna" in spec else []) + [t for t in ("mrna", "gro", "poly") if t in spec and "rna" not in spec]
    if rt:
        m = np.load("data/mrna_20bp.npz") if "mrna" in rt else None; Xr = {s: dict(np.load(f"data/X_rna_{s}.npz")) for s in SR}; R = {c: [] for c in P.CHROMS}
        for s in SR:
            if m is not None:
                raw = {c: np.log1p(m[f"{s}__{c}"].astype(np.float32)) for c in P.CHROMS}; tr = np.concatenate([raw[c] for c in P.TRAIN]); mu, sd = tr.mean(), tr.std() + 1e-6
            for c in P.CHROMS:
                for t in rt: R[c].append((raw[c] - mu) / sd if t == "mrna" else Xr[s][c][0 if t == "gro" else 1])
        R = {c: np.stack([v[: min(len(x) for x in R[c])] for v in R[c]]).astype(np.float16) for c in P.CHROMS}
    OH = dna_onehot() if "onehot" in spec else None
    G = glm() if "glm" in spec else None
    X = {}
    for c in P.CHROMS:
        parts = []
        if Xc is not None: parts += [Xc[s][c][ti].astype(np.float16) for s in SR]
        if R is not None: parts.append(R[c])
        if OH is not None: parts.append(OH[c])
        if G is not None: parts.append(G[c])
        n = min(p.shape[1] for p in parts); X[f"pool|{c}"] = np.concatenate([p[:, :n] for p in parts])
    n_unstr = (3 * len(ti) if Xc is not None else 0) + (3 * len(rt) if R is not None else 0)
    return X, n_unstr, (16 if G is not None else 0), OH is not None
def assemble(spec, tracks=None):
    if P.STAGES == ["pool"]: return assemble_pool(spec, tracks)
    """spec: set of {'chrom','rna','onehot','glm'} -> X dict, n_unstr, k_glm, has_onehot"""
    Xc = {f"{s}|{c}": v for s in S for c, v in np.load(f"data/X_chrom_{s}.npz").items()} if "chrom" in spec else None
    ti = list(range(25)) if tracks is None else [P.META["chrom_tracks"].index(t) for t in tracks]
    if Xc is not None: Xc = {k: v[ti] for k, v in Xc.items()}
    R = rna() if "rna" in spec else None
    OH = dna_onehot() if "onehot" in spec else None
    G = glm() if "glm" in spec else None
    if "glm" in spec: assert G is not None, "gLM embeddings missing"
    X = {}
    for s in S:
        for c in P.CHROMS:
            k = f"{s}|{c}"; parts = []
            if Xc is not None: parts.append(Xc[k].astype(np.float16))
            if R is not None: parts.append(R[k])
            if OH is not None: parts.append(OH[c])
            if G is not None: parts.append(G[c])
            n = min(p.shape[1] for p in parts); X[k] = np.concatenate([p[:, :n] for p in parts])
    n_unstr = (len(ti) if Xc is not None else 0) + (2 if R is not None else 0)
    return X, n_unstr, (16 if G is not None else 0), OH is not None
