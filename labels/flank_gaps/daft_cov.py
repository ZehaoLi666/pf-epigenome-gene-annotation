"""Strand-specific DAFT-seq coverage (Chappell 2020, ERP001570). Reverse-stranded library: RNA strand = opposite of read-1 strand.
Per time point (2 lanes pooled), per chromosome: 10-bp bin mean coverage, normalised to per-million read pairs (CPM).
Output: daft_cov10.npz with keys '<tp>|<chrom>' and 'pool|<chrom>' (mean of the 7 time-point CPM tracks), each float32 (2, nbins) [+ strand, - strand]."""
from pfconfig import HPCC_BIGDATA  # paths: config/paths.yaml
import sys, glob, os, collections, numpy as np, pysam
W = f"{HPCC_BIGDATA}/flank_coverage_2026"; RES = 10
B = collections.defaultdict(list)
for f in sorted(glob.glob(f"{W}/bam/*.bam")): B[int(f.split("_t")[-1][:-4])].append(f)
TPS = sorted(B); print("time points", TPS, {t: len(v) for t, v in B.items()}, flush=True)
h = pysam.AlignmentFile(B[TPS[0]][0]); CH = {c: L for c, L in zip(h.references, h.lengths) if c.startswith("Pf3D7_") and c[6:8].isdigit()}
npairs = collections.Counter(); OUT = {}
for c, L in CH.items():
    nb = (L + RES - 1) // RES; P = {}
    for t in TPS:
        d = np.zeros((2, L + 1), np.int32)
        for f in B[t]:
            bam = pysam.AlignmentFile(f); sp = [[], []]; ep = [[], []]
            for r in bam.fetch(c):
                if r.is_secondary or r.is_supplementary or not r.is_proper_pair: continue
                if r.is_read1:
                    npairs[t] += 1; plus = r.is_reverse
                else: plus = not r.is_reverse
                k = 0 if plus else 1
                for b0, b1 in r.get_blocks(): sp[k].append(b0); ep[k].append(b1)
            for k in (0, 1):
                if sp[k]: np.add.at(d[k], np.asarray(sp[k]), 1); np.add.at(d[k], np.asarray(ep[k]), -1)
        cov = np.cumsum(d[:, :L], axis=1).astype(np.float32)
        pad = nb * RES - L
        if pad: cov = np.concatenate([cov, np.zeros((2, pad), np.float32)], axis=1)
        P[t] = cov.reshape(2, nb, RES).mean(2)
    for t in TPS: OUT[f"{t}|{c}"] = P[t]
    print(c, flush=True)
for k in list(OUT):
    t = int(k.split("|")[0]); OUT[k] = (OUT[k] / (npairs[t] / 1e6)).astype(np.float32)
for c in CH: OUT[f"pool|{c}"] = np.mean([OUT[f"{t}|{c}"] for t in TPS], axis=0).astype(np.float32)
np.savez_compressed(f"{W}/daft_cov10.npz", **OUT)
open(f"{W}/daft_pairs.tsv", "w").write("".join(f"{t}\t{npairs[t]}\n" for t in TPS))
try:
    import pyBigWig
    for k, nm in ((0, "plus"), (1, "minus")):
        bw = pyBigWig.open(f"{W}/DAFTseq_pool_{nm}_CPM_10bp.bw", "w"); bw.addHeader([(c, CH[c]) for c in CH])
        for c in CH:
            v = OUT[f"pool|{c}"][k]; st = np.arange(len(v)) * RES; en = np.minimum(st + RES, CH[c])
            bw.addEntries([c] * len(v), st.tolist(), ends=en.tolist(), values=v.astype(float).tolist())
        bw.close()
except ImportError: print("no pyBigWig; bigWigs skipped")
print("ALLDONE", dict(npairs), flush=True)
