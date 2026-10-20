"""Caduceus-ph per-bin embeddings for the Pf 3D7 v48 nuclear genome, both strands.
Output: $HPCC_BIGDATA/boundary_seg/glm/caduceus_pca16_20bp.npz with keys
  <chrom>__fwd, <chrom>__rc : float16 (16, n_bins); rc = embedding of the reverse-complement
  strand, re-indexed to forward coordinates. PCA fitted on chromosomes 8-14 only."""
from pfconfig import HPCC_BIGDATA, HPCC_HOME, HPCC_LAB  # paths: config/paths.yaml
import os, sys, time, numpy as np, torch
os.environ.setdefault("HF_HOME", f"{HPCC_HOME}/.cache/huggingface"); os.environ["HF_HUB_OFFLINE"] = "1"; os.environ["TRANSFORMERS_OFFLINE"] = "1"
from transformers import AutoTokenizer, AutoModelForMaskedLM
FA = f"{HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta"
OUT = f"{HPCC_BIGDATA}/boundary_seg/glm"; os.makedirs(OUT, exist_ok=True)
RES, WIN, STEP, BATCH, NPC = 20, 4096, 2048, 16, 16
ck = "kuleshov-group/caduceus-ph_seqlen-131k_d_model-256_n_layer-16"
tok = AutoTokenizer.from_pretrained(ck, trust_remote_code=True, local_files_only=True)
model = AutoModelForMaskedLM.from_pretrained(ck, trust_remote_code=True, local_files_only=True).cuda().eval().to(torch.bfloat16)
T0 = time.time(); log = lambda *a: print(f"[{time.time()-T0:6.0f}s]", *a, flush=True)
seqs = {}; name = None
for line in open(FA):
    if line.startswith(">"): name = line[1:].split()[0].split("|")[0].strip(); seqs[name] = []
    else: seqs[name].append(line.strip().upper())
seqs = {k: "".join(v) for k, v in seqs.items() if k.startswith("Pf3D7_") and k[6:8].isdigit()}
log("chroms", len(seqs), sum(map(len, seqs.values())))
COMP = str.maketrans("ACGTN", "TGCAN")
def embed_track(s):
    """per-base hidden states pooled to RES bins, windows WIN/STEP keeping centre."""
    L = len(s); nb = L // RES; acc = np.zeros((nb * RES, 256), np.float32); cnt = np.zeros(nb * RES, np.float32)
    starts = list(range(0, max(1, L - WIN + 1), STEP))
    if starts[-1] + WIN < L: starts.append(max(0, L - WIN))
    for i in range(0, len(starts), BATCH):
        bs = starts[i:i + BATCH]
        ids = tok([s[a:a + WIN] for a in bs], add_special_tokens=False, return_tensors="pt", padding=True)["input_ids"].cuda()
        with torch.no_grad():
            h = model(input_ids=ids, output_hidden_states=True).hidden_states[-1].float().cpu().numpy()
        for j, a in enumerate(bs):
            w = min(WIN, L - a); lo = 0 if a == 0 else WIN // 4; hi = w if a + WIN >= L else 3 * WIN // 4
            e = min(a + hi, nb * RES)
            if a + lo >= e: continue
            acc[a + lo:e] += h[j, lo:e - a]; cnt[a + lo:e] += 1
    acc /= np.maximum(cnt, 1)[:, None]
    return acc.reshape(nb, RES, 256).mean(1).astype(np.float16)
raw = {}
for c, s in seqs.items():
    f = embed_track(s)
    r = embed_track(s.translate(COMP)[::-1])  # RC strand, in RC coordinates
    L = len(s); nb = L // RES
    # map RC bins back to forward coordinates: RC base k <-> forward base L-1-k
    rb = np.zeros_like(f); off = L - nb * RES  # leading bases trimmed in forward coords
    rc_full = np.repeat(r, RES, axis=0)[::-1]  # per-base, forward order, length nb*RES aligned to [L-nb*RES, L)
    rc_fwd = np.zeros((L, 256), np.float16); rc_fwd[off:] = rc_full
    rb = rc_fwd[:nb * RES].reshape(nb, RES, 256).astype(np.float32).mean(1).astype(np.float16)
    raw[c] = (f, rb); log(c, f.shape)
train = [c for c in raw if int(c[6:8]) >= 8]
rng = np.random.default_rng(0)
X = np.concatenate([np.concatenate([raw[c][0], raw[c][1]]) for c in train]).astype(np.float32)
X = X[rng.choice(len(X), min(400000, len(X)), replace=False)]
mu = X.mean(0); U, S, Vt = np.linalg.svd(X - mu, full_matrices=False); P = Vt[:NPC]
ev = (S[:NPC] ** 2).sum() / (S ** 2).sum()
out = {}
for c, (f, r) in raw.items():
    out[f"{c}__fwd"] = ((f.astype(np.float32) - mu) @ P.T).T.astype(np.float16)
    out[f"{c}__rc"] = ((r.astype(np.float32) - mu) @ P.T).T.astype(np.float16)
np.savez_compressed(f"{OUT}/caduceus_pca16_20bp.npz", **out, _pca_mean=mu.astype(np.float32), _pca_components=P.astype(np.float32), _explained=np.float32(ev))
log("done explained_var", float(ev))
