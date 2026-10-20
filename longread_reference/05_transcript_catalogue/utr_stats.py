"""Per annotated UTR (PlasmoDB v68 mRNAs; v48 and v68 nuclear sequences are identical) x long-read dataset: read-level
connectivity statistics, for validating over-long 5'/3' UTRs.

Coordinates are in transcription direction relative to the CDS:
  5' UTR: u = distance UPSTREAM of the first CDS base (u > 0 inside the 5' UTR; u < 0 inside the CDS). L5 = u of annotated TSS.
  3' UTR: v = distance DOWNSTREAM of the last CDS base (v > 0 inside the 3' UTR). L3 = v of annotated transcript end.
Only primary alignments with MAPQ >= 1 on the gene strand ("sense") are classified; opposite-strand reads are counted as antisense.
All datasets are strand-informative (ONT direct RNA = RNA sense; Iso-Seq FLNC and GSE186109 reads oriented 5'->3').

5' UTR read classes (sense reads whose span overlaps [TSS-200, CDS start+300]):
  anchored      : read 3' end >= 100 bp into the CDS (u3 <= -100)  -> informative about how far upstream the same molecule reaches
     reach_cds  : ... and read 5' end upstream of the CDS start (u5 > 0)
     reach_fXX  : ... and u5 >= XX% of L5 ; beyond_tss : u5 > L5 + 50
  utr_only      : read 5' end AND 3' end upstream of the CDS start (u3 >= 50): candidate independent upstream unit
     ..._Arich  : genomic 20 nt downstream of the read 3' end >= 60% A (possible internal priming; relevant for cDNA only)
     ..._tail   : non-templated A tail retained in a 3' soft clip (only seen in ONT cDNA)
  touch         : 5' end in UTR, 3' end within the first 100 bp of the CDS (ambiguous; counted separately)
  cds_start     : read 5' end inside the CDS (u5 <= 0): 5'-truncated relative to the UTR, uninformative except as anchor
5' end cap signal (G1): >= 1 non-templated G in the 5' soft clip (template-switching cDNA only).
Histograms (100-bp bins of u from -300 to 6000; '>' last bin): 5' ends of anchored reads, 5' and 3' ends of utr_only reads.

3' UTR read classes (sense reads overlapping [CDS end-300, TES+200]):
  cds_connected : read 5' end >= 100 bp upstream of the CDS end (v5 <= -100) and 3' end beyond the CDS end (v3 > 0)
     reach_fXX / beyond_tes on v3 (fraction of L3)
  utr3_only     : read 5' end >= 50 bp downstream of the CDS end (v5 >= 50); with G1 cap on its 5' end -> possible independent start
"""
from pfconfig import HPCC_BIGDATA, HPCC_LAB  # paths: config/paths.yaml
import sys, os, collections, numpy as np, pysam
GFF, OUT = sys.argv[1], sys.argv[2]
FA = pysam.FastaFile(f"{HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta")
R = f"{HPCC_BIGDATA}/longread_reprocess_2026"; G = f"{HPCC_BIGDATA}/GSE186109_longread/bam"
DS = collections.OrderedDict()
DS["LeeVV2021_dRNA"] = [f"{R}/ont/LeeVV_2021_dRNA.sorted.bam"]
DS["Hoshizaki2022_dRNA"] = [f"{R}/ont/Hoshizaki_2022_dRNA.sorted.bam"]
DS["Yang2021_IsoSeq"] = [f"{R}/pacbio/CRR2165{i}_flnc_all.sorted.bam" for i in range(16, 21)]
for st in ("ring", "troph", "schizont"):
    DS[f"GSE186109_PacBio_{st}"] = [f"{G}/{st}_{k}/{t}.sorted.bam" for k in ("total", "nuclear") for t in ("hq", "lq")]
DS["GSE186109_PacBio_schiztroph_cyto"] = [f"{G}/schiztroph_cytoplasmic/{t}.sorted.bam" for t in ("hq", "lq")]
DS["GSE186109_ONT_cDNA"] = [f"{G}/ont_mixed_oriented/oriented.sorted.bam"]
for k in list(DS):
    DS[k] = [p for p in DS[k] if os.path.exists(p)]
    assert DS[k], f"no BAM for {k}"
# GSE186109 PacBio records are Iso-Seq clustered isoforms: weight each by its number of full-length reads (n_fl)
W8 = {}
for i, line in enumerate(open(f"{HPCC_BIGDATA}/GSE186109_longread/isoform_fl_counts.tsv")):
    if i == 0: continue
    f = line.rstrip("\n").split("\t"); W8[f[2]] = max(int(f[4]), 1)
MAXINTRON = 5000
MINID = float(os.environ.get("MINID", "0"))     # e.g. 0.85 -> identity = 1 - de (gap-compressed divergence, minimap2)
MINCOV = float(os.environ.get("MINCOV", "0"))   # e.g. 0.99 -> aligned fraction of the read (TAMA -c)
def ok_aln(r, lo, hi):
    """primary, MAPQ>=1, no intron > MAXINTRON, and an aligned block overlapping [lo, hi)."""
    if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < 1: return False
    if any(op == 3 and n > MAXINTRON for op, n in r.cigartuples): return False
    if MINID > 0 and r.has_tag("de") and 1 - r.get_tag("de") < MINID: return False
    if MINCOV > 0:
        ql = r.infer_read_length() or 0
        if ql == 0 or r.query_alignment_length / ql < MINCOV: return False
    return any(b0 < hi and b1 > lo for b0, b1 in r.get_blocks())
comp = str.maketrans("ACGTN", "TGCAN"); rc = lambda s: s.translate(comp)[::-1]
NB = 63; BIN = lambda d: str(min(max(int((d + 300) // 100), 0), NB))   # bin 0 = [-300,-200) ... bin 63 = >= 6000

# ---------------- transcripts from GFF ----------------
tx = {}; cds = collections.defaultdict(list); utr5 = collections.defaultdict(list); utr3 = collections.defaultdict(list)
for line in open(GFF):
    if line.startswith("#") or "\t" not in line: continue
    f = line.rstrip("\n").split("\t")
    if not (f[0].startswith("Pf3D7_") and f[0][6:8].isdigit()): continue
    a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x); s, e = int(f[3]) - 1, int(f[4])   # 0-based half-open
    if f[2] == "mRNA": tx[a["ID"]] = dict(chrom=f[0], strand=f[6], start=s, end=e, gene=a.get("Parent", ""))
    elif f[2] == "CDS": cds[a["Parent"]].append((s, e))
    elif f[2] == "five_prime_UTR": utr5[a["Parent"]].append((s, e))
    elif f[2] == "three_prime_UTR": utr3[a["Parent"]].append((s, e))
rows = []
for t, d in tx.items():
    if t not in cds: continue
    cs, ce = min(x[0] for x in cds[t]), max(x[1] for x in cds[t])
    if d["strand"] == "+": cds5, cds3, L5, L3 = cs, ce - 1, cs - d["start"], d["end"] - ce
    else: cds5, cds3, L5, L3 = ce - 1, cs, d["end"] - ce, cs - d["start"]
    rows.append(dict(transcript=t, gene=d["gene"], chrom=d["chrom"], strand=d["strand"], tx_start=d["start"], tx_end=d["end"],
                     cds5=cds5, cds3=cds3, L5=L5, L3=L3, has5=t in utr5, has3=t in utr3,
                     utr5_exons=len(utr5.get(t, [])), utr3_exons=len(utr3.get(t, []))))
if len(sys.argv) > 3: rows = rows[:int(sys.argv[3])]
print("transcripts", len(rows), "with 5'UTR", sum(r["has5"] for r in rows), "with 3'UTR", sum(r["has3"] for r in rows), flush=True)

def read_ends(r):
    """sense 5' and 3' end (0-based genomic), soft clips in sense orientation."""
    q = r.query_sequence or ""; cig = r.cigartuples
    lc = cig[0][1] if cig[0][0] == 4 else 0; rcl = cig[-1][1] if cig[-1][0] == 4 else 0
    if not r.is_reverse:
        return r.reference_start, r.reference_end - 1, q[:lc], (q[len(q) - rcl:] if rcl else "")
    return r.reference_end - 1, r.reference_start, (rc(q[len(q) - rcl:]) if rcl else ""), rc(q[:lc])
def g1(clip5): return len(clip5) - len(clip5.rstrip("G")) >= 1 if clip5 else False
def tail(clip3): t = clip3[:15]; return len(t) >= 8 and t.count("A") >= 0.8 * len(t)
def arich(c, s, p3):
    if s == "+": w = FA.fetch(c, p3 + 1, p3 + 21).upper()
    else: w = rc(FA.fetch(c, max(0, p3 - 20), p3).upper())
    return len(w) > 0 and w.count("A") >= 0.6 * len(w)

FR = (25, 50, 75, 100)
out = open(OUT, "w"); hdr = None
for ds, paths in DS.items():
    bams = [pysam.AlignmentFile(p) for p in paths]; WEIGHT = ds.startswith("GSE186109_PacBio")
    for n_done, T in enumerate(rows):
        c, s = T["chrom"], T["strand"]; sg = 1 if s == "+" else -1
        u = lambda x: (T["cds5"] - x) * sg          # upstream distance from CDS start
        v = lambda x: (x - T["cds3"]) * sg          # downstream distance from CDS end
        K = collections.Counter(); H = {k: collections.Counter() for k in ("anch5", "uo5", "uo3", "c3")}
        # ---- 5' UTR window ----
        if s == "+": lo, hi = T["tx_start"] - 200, T["cds5"] + 300
        else: lo, hi = T["cds5"] - 300, T["tx_end"] + 200
        seen = set()
        for b in bams:
            for r in b.fetch(c, max(0, lo), hi):
                if r.query_name in seen or not ok_aln(r, max(0, lo), hi): continue
                seen.add(r.query_name); w = W8.get(r.query_name, 1) if WEIGHT else 1
                rs = "-" if r.is_reverse else "+"
                if rs != s:
                    # antisense: overlaps the 5' UTR proper (not the CDS)
                    if T["L5"] > 0 and r.reference_end > (T["tx_start"] if s == "+" else T["cds5"] + 1) and r.reference_start < (T["cds5"] if s == "+" else T["tx_end"]):
                        K["as5_reads"] += w; p5, _, cl5, _ = read_ends(r)
                        K["as5_G1"] += w * g1(cl5)
                    continue
                p5, p3, cl5, cl3 = read_ends(r); u5, u3 = u(p5), u(p3)
                if u3 <= -100:
                    K["anchored"] += w; H["anch5"][BIN(u5)] += w
                    if u5 > 0:
                        K["reach_cds"] += w
                        for f in FR:
                            if T["L5"] > 0 and u5 >= f / 100 * T["L5"]: K[f"reach_f{f}"] += w
                        if u5 > T["L5"] + 50: K["beyond_tss"] += w
                        K["anchored_G1"] += w * g1(cl5)
                elif u5 > 0 and u3 >= 50:
                    K["utr_only"] += w; H["uo5"][BIN(u5)] += w; H["uo3"][BIN(u3)] += w
                    K["utr_only_Arich"] += w * arich(c, s, p3); K["utr_only_tail"] += w * tail(cl3); K["utr_only_G1"] += w * g1(cl5)
                    K["utr_only_in_utr"] += w * (u5 <= T["L5"] + 50)
                elif u5 > 0: K["touch"] += w
                else: K["cds_start"] += w
        # ---- 3' UTR window ----
        if s == "+": lo, hi = T["cds3"] - 300, T["tx_end"] + 200
        else: lo, hi = T["tx_start"] - 200, T["cds3"] + 300
        seen = set()
        for b in bams:
            for r in b.fetch(c, max(0, lo), hi):
                if r.query_name in seen or ("-" if r.is_reverse else "+") != s or not ok_aln(r, max(0, lo), hi): continue
                seen.add(r.query_name); w = W8.get(r.query_name, 1) if WEIGHT else 1
                p5, p3, cl5, cl3 = read_ends(r); v5, v3 = v(p5), v(p3)
                if v5 <= -100 and v3 > 0:
                    K["cds_connected3"] += w; H["c3"][BIN(v3)] += w
                    for f in FR:
                        if T["L3"] > 0 and v3 >= f / 100 * T["L3"]: K[f"reach3_f{f}"] += w
                    if v3 > T["L3"] + 50: K["beyond_tes"] += w
                    K["cds_connected3_Arich"] += w * arich(c, s, p3); K["cds_connected3_tail"] += w * tail(cl3)
                elif v5 >= 50:
                    K["utr3_only"] += w; K["utr3_only_G1"] += w * g1(cl5)
                elif v3 <= 0: K["ends_in_cds3"] += w
        cols = ["anchored", "reach_cds"] + [f"reach_f{f}" for f in FR] + ["beyond_tss", "anchored_G1", "utr_only", "utr_only_in_utr", "utr_only_Arich",
                "utr_only_tail", "utr_only_G1", "touch", "cds_start", "as5_reads", "as5_G1", "cds_connected3"] + [f"reach3_f{f}" for f in FR] + \
               ["beyond_tes", "cds_connected3_Arich", "cds_connected3_tail", "utr3_only", "utr3_only_G1", "ends_in_cds3"]
        hc = ["hist_anchored_5end", "hist_utronly_5end", "hist_utronly_3end", "hist_cds3_connected_3end"]
        if hdr is None:
            hdr = ["dataset"] + list(T) + cols + hc; out.write("\t".join(hdr) + "\n")
        hs = [",".join(f"{k}:{n}" for k, n in sorted(H[x].items(), key=lambda kv: int(kv[0]))) for x in ("anch5", "uo5", "uo3", "c3")]
        out.write("\t".join([ds] + [str(T[k]) for k in T] + [str(K[k]) for k in cols] + hs) + "\n")
        if n_done % 1000 == 0: print(ds, n_done, flush=True)
    out.flush(); print("done", ds, flush=True)
out.close(); print("ALLDONE")
