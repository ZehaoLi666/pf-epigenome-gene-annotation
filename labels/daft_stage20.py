"""Stage-binned stranded DAFT-seq input: 20-bp bins, RNA strand (row 0 = +, row 1 = -), mean CPM of the stage's time points, log1p, float16.
ring = 8,16 h; troph = 24,32 h; schizont = 40 h (Chappell 2020, ERP001570)."""
from pfconfig import HPCC_BIGDATA  # paths: config/paths.yaml
import numpy as np
W = f"{HPCC_BIGDATA}/flank_coverage_2026"; z = np.load(f"{W}/daft_cov10.npz")
ST = {"ring": (8, 16), "troph": (24, 32), "schizont": (40,)}
chroms = sorted({k.split("|")[1] for k in z.files if not k.startswith("pool")})
LEN = {c: int(l) for c, l in (line.split() for line in open(f"{W}/chrom_sizes.txt"))} if False else None
out = {}
for c in chroms:
    for s, tps in ST.items():
        v = np.mean([z[f"{t}|{c}"] for t in tps], axis=0)          # (2, n10)
        n20 = v.shape[1] // 2; v = v[:, :2 * n20].reshape(2, n20, 2).mean(2)
        out[f"{s}__{c}"] = np.log1p(v).astype(np.float16)
np.savez_compressed(f"{W}/daft_stage20.npz", **out); print("ok", len(out))
