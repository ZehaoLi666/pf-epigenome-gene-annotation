"""Session-derived computation of hc_model_evidence.tsv.gz.

Recovered from the record of the interactive analysis session in which hc_model_evidence.tsv.gz was produced
(it was not saved as a stand-alone script). It documents the computation; it may rely on objects defined
earlier in that session, so check it before re-running.
Inputs recorded for this table: (none recorded)
Files named in result(...) are read from RESULTS_DIR (config/paths.yaml); code files named there
(e.g. common.py, tracks.py, cv_score.py) are the repository versions.
"""
from pfconfig import HPCC_BIGDATA, HPCC_LAB  # paths: config/paths.yaml
import os, re, gzip, glob, bisect, collections, math, csv, json, random
import numpy as np
from scipy.stats import chi2_contingency, fisher_exact
import pysam

P   = f"{HPCC_BIGDATA}/GSE186109_longread"
COL = P + "/tama/collapse_up"
MRG = P + "/tama/merge"
OUT = P + "/reference_v1"; os.makedirs(OUT, exist_ok=True)
GEN = f"{HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta"
GTF = f"{HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7.gtf"
CAGE_RAW = f"{HPCC_BIGDATA}/cage_tss_2026/raw"
CAGE_MAN = f"{HPCC_BIGDATA}/cage_tss_2026/manifest.tsv"
FILELIST = MRG + "/merge_filelist.txt"
PFX = "MASTER"
W_TSS, W_TES, W_CAGE = 20, 50, 50
THREE_P = {"LeeVV2021", "Marcus", "Mengquan"}
STAGES = ["ring", "troph", "schizont"]

# ---------------- GSE186109 library table ----------------
srr_lib = {}
for row in csv.DictReader(open(P + "/manifest_isoform_tier.csv")):
    srr_lib[row["srr"]] = (row["sample"], row["stage"], row["fraction"], row["tier"])
srr_lib["SRR18188665"] = ("ont_mixed", "mixed", "total", "ont")
srr_lib["SRR18188666"] = ("ont_mixed", "mixed", "total", "ont")
flw = {}
for row in csv.DictReader(open(P + "/isoform_fl_counts.tsv"), delimiter="\t"):
    flw[row["read"]] = max(1, int(row["n_fl"] or 0))
def lib_of(read):
    srr = read.split(".")[0]
    return srr_lib.get(srr)

# ---------------- sources and source transcripts ----------------
sources = []
for ln in open(FILELIST):
    f = ln.rstrip("\n").split("\t")
    if len(f) >= 4: sources.append((f[3], f[0], f[1]))
src_names = [s[0] for s in sources]
def dataset_of(src):
    return "GSE186109" if src.startswith("GSE_") else src

# merge.txt: source transcript coords -> merged transcript
src2m = {}
for ln in open(f"{MRG}/{PFX}_merge.txt"):
    f = ln.rstrip("\n").split("\t")
    mid, stx = f[3].split(";")
    src2m[stx] = mid

# read -> (source, merged tx, read start, read end, strand, blocks)
reads_by_m = collections.defaultdict(list)
for src, bedpath, cap in sources:
    trb = bedpath.replace(".bed", "_trans_read.bed")
    for ln in open(trb):
        f = ln.rstrip("\n").split("\t")
        tid, rid = f[3].split(";", 1)
        m = src2m.get(f"{src}_{tid}")
        if m is None: continue
        reads_by_m[m].append((src, rid, int(f[1]), int(f[2]), f[5]))

# ---------------- master models ----------------
models = []
for ln in open(f"{MRG}/{PFX}.bed"):
    f = ln.rstrip("\n").split("\t")
    gid, tid = f[3].split(";")
    s, e = int(f[1]), int(f[2])
    sz = [int(x) for x in f[10].rstrip(",").split(",")]
    st = [int(x) for x in f[11].rstrip(",").split(",")]
    ex = [(s + a, s + a + b) for a, b in zip(st, sz)]
    introns = [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]
    models.append(dict(gene=gid, tx=tid, chrom=f[0], start=s, end=e, strand=f[5],
                       exons=ex, introns=introns, bed=f))
trep = {}
for row in csv.DictReader(open(f"{MRG}/{PFX}_trans_report.txt"), delimiter="\t"):
    trep[row["transcript_id"]] = row

# ---------------- v48 annotation ----------------
v48_ex = collections.defaultdict(list); v48_int = set(); v48_gene = {}
tx_ex = collections.defaultdict(list)
for ln in open(GTF):
    if ln.startswith("#"): continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] != "exon": continue
    gid = re.search('gene_id "([^"]+)"', f[8]).group(1)
    tid = re.search('transcript_id "([^"]+)"', f[8]).group(1)
    s, e = int(f[3]) - 1, int(f[4])
    tx_ex[(f[0], f[6], gid, tid)].append((s, e))
    g = v48_gene.get(gid)
    v48_gene[gid] = [f[0], min(s, g[1]) if g else s, max(e, g[2]) if g else e, f[6]]
for (ch, sd, gid, tid), ex in tx_ex.items():
    ex.sort()
    for i in range(len(ex) - 1): v48_int.add((ch, sd, ex[i][1], ex[i + 1][0]))
gidx = collections.defaultdict(list)
for gid, (ch, s, e, sd) in v48_gene.items(): gidx[ch].append((s, e, sd, gid))
for ch in gidx: gidx[ch].sort()
gstarts = {ch: [x[0] for x in L] for ch, L in gidx.items()}
def v48_overlap(ch, s, e, sd):
    L = gidx.get(ch, []); i = bisect.bisect_right(gstarts.get(ch, []), e)
    same, anti = [], []
    for j in range(max(0, i - 80), i):
        a, b, gs, gid = L[j]
        if b > s and a < e:
            ov = min(b, e) - max(a, s)
            (same if gs == sd else anti).append((ov, gid))
    same.sort(reverse=True); anti.sort(reverse=True)
    n_multi = sum(1 for ov, g in same if ov >= 0.2 * (e - s) or ov >= 0.5 * (v48_gene[g][2] - v48_gene[g][1]))
    return (same[0][1] if same else ""), (anti[0][1] if anti else ""), n_multi

fa = pysam.FastaFile(GEN)
def down_A_frac(ch, sd, tes, n=20):
    L = fa.get_reference_length(ch)
    if sd == "+":
        s = fa.fetch(ch, tes, min(L, tes + n)).upper(); return s.count("A") / max(len(s), 1)
    s = fa.fetch(ch, max(0, tes - n), tes).upper(); return s.count("T") / max(len(s), 1)

# ---------------- CAGE (per stage, per replicate) ----------------
cage = {}  # (stage, rep, strand, chrom) -> (pos array, cumsum)
man = list(csv.DictReader(open(CAGE_MAN), delimiter="\t"))
acc = collections.defaultdict(lambda: collections.defaultdict(float))
cage_tot = collections.Counter()
for row in man:
    st = {"ring": "ring", "trophozoite": "troph", "schizont": "schizont"}[row["stage"]]
    sd = "+" if row["strand"] == "POS" else "-"
    with gzip.open(os.path.join(CAGE_RAW, row["file"]), "rt") as fh:
        for ln in fh:
            f = ln.split("\t")
            v = abs(float(f[3]))
            if v <= 0: continue
            acc[(st, row["rep"], sd, f[0])][int(f[1])] += v
            cage_tot[(st, row["rep"])] += v
for k, d in acc.items():
    pos = np.array(sorted(d)); val = np.array([d[p] for p in pos]); cage[k] = (pos, np.concatenate([[0], np.cumsum(val)]))
del acc
def cage_sum(stage, rep, sd, ch, x, w=W_CAGE):
    t = cage.get((stage, rep, sd, ch))
    if t is None: return 0.0
    pos, cs = t
    i = np.searchsorted(pos, x - w, "left"); j = np.searchsorted(pos, x + w, "right")
    return float(cs[j] - cs[i])
def cage_cpm(stage, sd, ch, x):
    return {rep: 1e6 * cage_sum(stage, rep, sd, ch, x) / max(cage_tot[(stage, rep)], 1) for rep in ("A", "B")}

# ---------------- per-model evidence ----------------
rows = []; stage_counts = {}
tot_fl = collections.Counter()
for m in reads_by_m.values():
    for src, rid, a, b, sd in m:
        lib = lib_of(rid)
        if lib and lib[1] in STAGES: tot_fl[lib[1]] += flw.get(rid, 1)
for M in models:
    tid, sd, ch = M["tx"], M["strand"], M["chrom"]
    tss = M["start"] if sd == "+" else M["end"] - 1
    tes = M["end"] - 1 if sd == "+" else M["start"]
    R = reads_by_m.get(tid, [])
    cap_libs, cap_fl, nuc_fl = set(), 0, 0
    tes_ds, span = set(), False
    n_reads = 0; n_fl = 0; tes_reads = 0
    stc = collections.Counter(); stn = collections.Counter()
    ds_all = set()
    for src, rid, a, b, rsd in R:
        r5 = a if sd == "+" else b - 1; r3 = b - 1 if sd == "+" else a
        ds = dataset_of(src); ds_all.add(ds)
        lib = lib_of(rid) if ds == "GSE186109" else None
        w = flw.get(rid, 1) if lib else 1
        if lib:
            if lib[1] in STAGES: stc[lib[1]] += w
            if lib[1] in STAGES and lib[2] == "nuclear": stn[lib[1]] += w
            if abs(r5 - tss) <= W_TSS:
                cap_libs.add((lib[0], lib[2])); cap_fl += w
                if lib[2] == "nuclear": nuc_fl += w
        n_reads += 1; n_fl += w
        if abs(r3 - tes) <= W_TES: tes_ds.add(ds); tes_reads += 1
        if abs(r5 - tss) <= W_TES and abs(r3 - tes) <= W_TES: span = True
    cg = {s: cage_cpm(s, sd, ch, tss) for s in STAGES}
    cage_any = max(v for s in STAGES for v in cg[s].values())
    cage_reps = sum(1 for s in STAGES for v in cg[s].values() if v >= 1.0)
    af = down_A_frac(ch, sd, tes)
    ip = af >= 0.8
    three_src = {dataset_of(x) for x in trep.get(tid, {}).get("sources", "").split(",") if x} & THREE_P
    nint = len(M["introns"])
    int_v48 = sum((ch, sd, a, b) in v48_int for a, b in M["introns"])
    int_multi = nint if len(ds_all) >= 2 else 0
    same_g, anti_g, n_multi = v48_overlap(ch, M["start"], M["end"], sd)
    tss_strong = (len(cap_libs) >= 2) or (len(cap_libs) >= 1 and cage_any >= 1.0)
    tss_weak = len(cap_libs) >= 1 or cage_any >= 1.0
    tes_strong = (len(tes_ds) >= 2) and not ip
    tes_weak = (len(tes_ds) >= 1) and not ip
    junc_ok = (nint == 0) or (int_v48 == nint) or (int_multi == nint)
    if tss_strong and tes_strong and junc_ok and span: tier = "Gold"
    elif tss_weak and tes_weak and junc_ok and span: tier = "Silver"
    else: tier = "Partial"
    reason = []
    if not tss_weak: reason.append("no_5p_evidence")
    if not tes_weak: reason.append("3p_IP_or_unsupported" if ip else "no_3p_evidence")
    if not span: reason.append("no_spanning_read")
    if not junc_ok: reason.append("junction_single_dataset_unannotated")
    rows.append(dict(transcript_id=tid, gene_id=M["gene"], chrom=ch, start=M["start"], end=M["end"], strand=sd,
        n_exons=len(M["exons"]), length=sum(b - a for a, b in M["exons"]),
        sources=trep.get(tid, {}).get("sources", ""), datasets=",".join(sorted(ds_all)),
        n_3p_datasets=len(three_src), tss=tss, tss_cap_libraries=len(cap_libs), tss_cap_FL=cap_fl, tss_nuclear_FL=nuc_fl,
        cage_ring_cpm=round(np.mean(list(cg["ring"].values())), 3), cage_troph_cpm=round(np.mean(list(cg["troph"].values())), 3),
        cage_schizont_cpm=round(np.mean(list(cg["schizont"].values())), 3), cage_reps_ge1cpm=cage_reps,
        tes=tes, tes_datasets=len(tes_ds), tes_down20_A=round(af, 2), tes_internal_priming_flag=int(ip),
        spanning_read=int(span), n_introns=nint, introns_in_v48=int_v48,
        v48_gene_same_strand=same_g, v48_gene_antisense=anti_g, n_v48_genes_spanned=n_multi,
        n_reads=n_reads, n_reads_FLweighted=n_fl, tes_reads=tes_reads,
        FL_ring=stc["ring"], FL_troph=stc["troph"], FL_schizont=stc["schizont"],
        FLnuc_ring=stn["ring"], FLnuc_troph=stn["troph"], FLnuc_schizont=stn["schizont"],
        tier=tier, tier_reason=";".join(reason) if reason else "all_criteria_met"))

# structural class vs v48
for r in rows:
    if r["v48_gene_same_strand"]: r["class"] = "v48_gene_overlap"
    elif r["v48_gene_antisense"]: r["class"] = "antisense_to_v48"
    else: r["class"] = "intergenic_novel"

# ---------------- locus key (TAMA merge genes chain neighbouring Pf genes into mega-loci) ----------------
mix = {M["tx"]: M for M in models}
nov = collections.defaultdict(list)
for r in rows:
    if r["v48_gene_same_strand"]: r["locus"] = r["v48_gene_same_strand"]
    else: nov[r["gene_id"]].append(r)
for g, L in nov.items():
    L.sort(key=lambda r: r["start"]); comp = 0; cur_end = -1
    for r in L:
        if r["start"] >= cur_end: comp += 1
        cur_end = max(cur_end, r["end"]); r["locus"] = f"novel_{g}_{comp}"

# ---------------- decay candidates (alternative internal TSS) ----------------
by_gene = collections.defaultdict(list)
for r in rows: by_gene[r["locus"]].append(r)
for g, L in by_gene.items():
    main = max(L, key=lambda r: (r["tss_cap_FL"], r["FL_ring"] + r["FL_troph"] + r["FL_schizont"]))
    for r in L:
        internal = (r["strand"] == "+" and r["tss"] > main["tss"] + W_CAGE) or (r["strand"] == "-" and r["tss"] < main["tss"] - W_CAGE)
        r["tss_role"] = "primary" if r is main else ("internal_alt" if internal else "upstream_alt")
        r["decay_candidate"] = int(r["tss_role"] == "internal_alt" and r["tss_cap_libraries"] <= 1
                                   and r["cage_reps_ge1cpm"] == 0 and r["tss_nuclear_FL"] == 0)
        if r["decay_candidate"] and r["tier"] != "Partial":
            r["tier"] = "Partial"; r["tier_reason"] = (r["tier_reason"].replace("all_criteria_met", "") + ";decay_candidate_5p").strip(";")

# ---------------- detectability (lambda) per stage ----------------
for r in rows:
    for s in STAGES:
        other = [x for x in STAGES if x != s]
        k_o = sum(r["FL_" + x] for x in other); n_o = sum(tot_fl[x] for x in other)
        lam = tot_fl[s] * k_o / n_o if n_o else 0.0
        r["lambda_" + s] = round(lam, 2)
        if r["FL_" + s] > 0: r["detect_" + s] = "detected"
        elif lam >= 3: r["detect_" + s] = "not_detected_powered"
        else: r["detect_" + s] = "not_assessable"

for r in rows:
    r["high_confidence"] = int(r["tier"] in ("Gold", "Silver") and not r["decay_candidate"] and
                               (r["n_reads_FLweighted"] >= 2 or len(r["datasets"].split(",")) >= 2))

cols = list(rows[0].keys())
with open(f"{OUT}/master_model_evidence.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows: fh.write("\t".join(str(r[c]) for c in cols) + "\n")

# HC-only evidence table
with open(f"{OUT}/hc_model_evidence.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        if r["high_confidence"]: fh.write("\t".join(str(r[c]) for c in cols) + "\n")

import shutil
shutil.copy(f"{OUT}/hc_model_evidence.tsv", "hc_model_evidence.tsv")
import gzip
with open("hc_model_evidence.tsv", "rb") as f_in, gzip.open("hc_model_evidence.tsv.gz", "wb") as f_out:
    f_out.write(f_in.read())
