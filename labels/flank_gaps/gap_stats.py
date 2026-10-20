"""Per same-strand gap between a long-read confirmed gene and its neighbouring long-read unit: evidence that the gap separates two
transcription units (spec UTR_BOUNDARY_MODEL_DESIGN.md v2, section 2).
gaps.tsv: gap_id chrom strand a b tts_pos tss_pos   ([a,b) = gap; tts_pos = 3' end of the upstream unit, tss_pos = 5' end of the downstream unit)
Per long-read dataset (sense reads only, primary, MAPQ>=1, no intron > 5 kb; GSE186109 PacBio isoforms weighted by n_FL):
  left / right : reads with an aligned block within 300 bp inside the left / right unit edge
  span         : reads with aligned blocks on BOTH sides (cover or splice across the gap) -> same molecule
  n3, n3_tail, n3_arich : read 3' ends within +-50 bp of tts_pos; with non-templated A tail in the 3' clip; genomic A-rich downstream
  n5, n5_G1    : read 5' ends within +-50 bp of tss_pos; with >=1 non-templated G in the 5' clip (template-switching cDNA)
DAFT-seq (stranded, pooled CPM, 10-bp bins): mean over gap, minimum 100-bp window in gap, mean of 300 bp inside each unit edge, and
antisense-strand gap mean; plus the same per time point."""
from pfconfig import HPCC_BIGDATA, HPCC_LAB  # paths: config/paths.yaml
import sys, os, collections, numpy as np, pysam
W = f"{HPCC_BIGDATA}/flank_coverage_2026"; R = f"{HPCC_BIGDATA}/longread_reprocess_2026"; G = f"{HPCC_BIGDATA}/GSE186109_longread/bam"
FA = pysam.FastaFile(f"{HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta")
DS = collections.OrderedDict()
DS["LeeVV2021_dRNA"] = [f"{R}/ont/LeeVV_2021_dRNA.sorted.bam"]
DS["Hoshizaki2022_dRNA"] = [f"{R}/ont/Hoshizaki_2022_dRNA.sorted.bam"]
DS["Yang2021_IsoSeq"] = [f"{R}/pacbio/CRR2165{i}_flnc_all.sorted.bam" for i in range(16, 21)]
for st in ("ring", "troph", "schizont"):
    DS[f"GSE186109_PacBio_{st}"] = [f"{G}/{st}_{k}/{t}.sorted.bam" for k in ("total", "nuclear") for t in ("hq", "lq")]
DS["GSE186109_ONT_cDNA"] = [f"{G}/ont_mixed_oriented/oriented.sorted.bam"]
for k in list(DS): DS[k] = [p for p in DS[k] if os.path.exists(p)]; assert DS[k], k
CDNA = {"Yang2021_IsoSeq", "GSE186109_PacBio_ring", "GSE186109_PacBio_troph", "GSE186109_PacBio_schizont", "GSE186109_ONT_cDNA"}
W8 = {}
for i, line in enumerate(open(f"{HPCC_BIGDATA}/GSE186109_longread/isoform_fl_counts.tsv")):
    if i == 0: continue
    f = line.rstrip("\n").split("\t"); W8[f[2]] = max(int(f[4]), 1)
comp = str.maketrans("ACGTN", "TGCAN"); rc = lambda s: s.translate(comp)[::-1]
def ok(r): return not (r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 1) and not any(op == 3 and n > 5000 for op, n in r.cigartuples)
def read_ends(r):
    q = r.query_sequence or ""; cig = r.cigartuples
    lc = cig[0][1] if cig[0][0] == 4 else 0; rcl = cig[-1][1] if cig[-1][0] == 4 else 0
    if not r.is_reverse: return r.reference_start, r.reference_end - 1, q[:lc], (q[len(q) - rcl:] if rcl else "")
    return r.reference_end - 1, r.reference_start, (rc(q[len(q) - rcl:]) if rcl else ""), rc(q[:lc])
def g1(c5): return bool(c5) and len(c5) - len(c5.rstrip("G")) >= 1
def tail(c3): t = c3[:15]; return len(t) >= 8 and t.count("A") >= 0.8 * len(t)
def arich(c, s, p3):
    w = FA.fetch(c, p3 + 1, p3 + 21).upper() if s == "+" else rc(FA.fetch(c, max(0, p3 - 20), p3).upper())
    return len(w) > 0 and w.count("A") >= 0.6 * len(w)
COV = np.load(f"{W}/daft_cov10.npz"); TPS = sorted({k.split("|")[0] for k in COV.files} - {"pool"}, key=int)
def covstats(key, c, s, a, b):
    v = COV[f"{key}|{c}"]; k = 0 if s == "+" else 1; ia, ib = a // 10, max(a // 10 + 1, (b + 9) // 10)
    g = v[k, ia:ib]; win = np.convolve(g, np.ones(10) / 10, "valid") if len(g) >= 10 else np.array([g.mean()])
    left = v[k, max(0, ia - 30):ia].mean() if ia > 0 else np.nan; right = v[k, ib:ib + 30].mean()
    return g.mean(), win.min(), left, right, v[1 - k, ia:ib].mean()
gaps = [l.rstrip("\n").split("\t") for l in open(sys.argv[1])]; hdr = gaps[0]; gaps = [dict(zip(hdr, x)) for x in gaps[1:]]
if len(sys.argv) > 3: gaps = gaps[:int(sys.argv[3])]
out = open(sys.argv[2], "w")
cols = ["left", "right", "span", "n3", "n3_tail", "n3_arich", "n5", "n5_G1"]
out.write("\t".join(["gap_id", "dataset"] + cols) + "\n")
cv = open(sys.argv[2].replace(".tsv", "_daft.tsv"), "w")
cv.write("\t".join(["gap_id", "track", "gap_mean", "gap_min100", "left_edge", "right_edge", "antisense_gap_mean"]) + "\n")
for i, g in enumerate(gaps):
    c, s = g["chrom"], g["strand"]; a, b, p3, p5 = int(g["a"]), int(g["b"]), int(g["tts_pos"]), int(g["tss_pos"])
    for key in ["pool"] + TPS:
        cv.write("\t".join([g["gap_id"], key] + [f"{x:.4f}" for x in covstats(key, c, s, a, b)]) + "\n")
    for ds, paths in DS.items():
        K = collections.Counter(); wt = ds.startswith("GSE186109_PacBio")
        for p in paths:
            bam = pysam.AlignmentFile(p)
            for r in bam.fetch(c, max(0, a - 350), b + 350):
                if not ok(r) or (r.is_reverse != (s == "-")): continue
                w = W8.get(r.query_name, 1) if wt else 1
                bl = r.get_blocks()
                L = any(x0 < a and x1 > a - 300 for x0, x1 in bl); Rr = any(x0 < b + 300 and x1 > b for x0, x1 in bl)
                K["left"] += w * L; K["right"] += w * Rr; K["span"] += w * (L and Rr)
                e5, e3, c5, c3 = read_ends(r)
                if abs(e3 - p3) <= 50:
                    K["n3"] += w; K["n3_tail"] += w * tail(c3); K["n3_arich"] += w * arich(c, s, e3)
                if abs(e5 - p5) <= 50:
                    K["n5"] += w; K["n5_G1"] += w * (g1(c5) and ds in CDNA)
        out.write("\t".join([g["gap_id"], ds] + [str(K[k]) for k in cols]) + "\n")
    if i % 200 == 0: print(i, flush=True); out.flush(); cv.flush()
out.close(); cv.close(); print("ALLDONE", flush=True)
