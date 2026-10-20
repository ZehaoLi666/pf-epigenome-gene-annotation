"""Transcript classes + transcription-context classification of every validated TSS/TTS.
Inputs: HC reference (reference_v1), v48 GTF/GFF, short-read junctions + TPM, optional SQANTI3 classification."""
from pfconfig import HPCC_BIGDATA, HPCC_LAB  # paths: config/paths.yaml
import os, re, bisect, collections, json, sys
import numpy as np, pandas as pd, pysam

P    = f"{HPCC_BIGDATA}/GSE186109_longread"
REF  = P + "/reference_v1"
SR   = P + "/shortread"
OUT  = P + "/contexts"; os.makedirs(OUT + "/bed", exist_ok=True)
GEN  = f"{HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta"
GTF  = f"{HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7.gtf"
GFF  = f"{HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7.gff"
SQ   = P + "/sqanti3/qc/hc_classification.txt"
PROX = 500          # promoter / terminator proximal window (bp)
MIN_ORF_AA = 100
STAGES = ["Ring", "Trophozoite", "Schizont", "Gametocyte", "ET", "LT", "Gii", "GV", "Ookinete"]
fa = pysam.FastaFile(GEN)

# ---------------- v48 genes: type, description, CDS, span ----------------
gtype, gdesc = {}, {}
for ln in open(GFF):
    if ln.startswith("#"): continue
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9: continue
    a = dict(x.split("=", 1) for x in f[8].split(";") if "=" in x)
    if f[2] == "gene":
        gdesc[a.get("ID")] = a.get("description", "").replace("+", " ")
    elif f[2] in ("mRNA", "tRNA", "rRNA", "snoRNA", "snRNA", "ncRNA"):
        par = a.get("Parent", "")
        gtype[par] = "mRNA" if f[2] == "mRNA" else ("ncRNA_" + f[2])
v48 = {}
for ln in open(GTF):
    f = ln.rstrip("\n").split("\t")
    if len(f) < 9 or f[2] not in ("exon", "CDS"): continue
    gid = re.search('gene_id "([^"]+)"', f[8]).group(1)
    s, e = int(f[3]) - 1, int(f[4])
    g = v48.setdefault(gid, dict(chrom=f[0], strand=f[6], start=s, end=e, cds_s=None, cds_e=None))
    g["start"] = min(g["start"], s); g["end"] = max(g["end"], e)
    if f[2] == "CDS":
        g["cds_s"] = s if g["cds_s"] is None else min(g["cds_s"], s)
        g["cds_e"] = e if g["cds_e"] is None else max(g["cds_e"], e)

# ---------------- HC models, exons, junction support ----------------
ev = pd.read_csv(REF + "/hc_model_evidence.tsv", sep="\t", low_memory=False)
exons = collections.defaultdict(list)
for ln in open(REF + "/hc_reference.gtf"):
    f = ln.rstrip("\n").split("\t")
    if f[2] != "exon": continue
    tid = re.search('transcript_id "([^"]+)"', f[8]).group(1)
    exons[tid].append((int(f[3]) - 1, int(f[4])))
sj = pd.read_csv(SR + "/junctions/pooled_filtered.SJ.out.tab", sep="\t", header=None,
                 names=["chrom", "start", "end", "strand", "motif", "annot", "uniq", "multi", "ovh"])
sjset = set(zip(sj.chrom, sj.start - 1, sj.end))           # 0-based intron start, end-exclusive
sjc = dict(zip(zip(sj.chrom, sj.start - 1, sj.end), sj.uniq))
jrows = []
for r in ev.itertuples():
    ex = sorted(exons[r.transcript_id]); intr = [(ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)]
    sup = [(r.chrom, a, b) in sjset for a, b in intr]
    jrows.append((r.transcript_id, len(intr), sum(sup), min([sjc.get((r.chrom, a, b), 0) for a, b in intr], default=0)))
jt = pd.DataFrame(jrows, columns=["transcript_id", "n_introns", "n_introns_sr_supported", "min_sr_reads_per_intron"])
jt["all_introns_sr_supported"] = ((jt.n_introns == 0) | (jt.n_introns == jt.n_introns_sr_supported)).astype(int)
ev = ev.merge(jt.drop(columns="n_introns"), on="transcript_id")
jt.to_csv(OUT + "/hc_model_sr_junction_support.tsv", sep="\t", index=False)

# ---------------- ORF (longest ATG..stop on the transcript strand) ----------------
comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
stops = {"TAA", "TAG", "TGA"}
def longest_orf(ch, sd, ex):
    seq = "".join(fa.fetch(ch, a, b) for a, b in sorted(ex)).upper()
    if sd == "-": seq = seq.translate(comp)[::-1]
    best = 0
    for fr in range(3):
        start = None
        for i in range(fr, len(seq) - 2, 3):
            c = seq[i:i + 3]
            if start is None and c == "ATG": start = i
            elif start is not None and c in stops:
                best = max(best, (i - start) // 3); start = None
    return best

# ---------------- optional SQANTI3 ----------------
sq = None
if os.path.exists(SQ):
    sq = pd.read_csv(SQ, sep="\t", low_memory=False)
    keep = [c for c in ["isoform", "structural_category", "subcategory", "coding", "ORF_length", "RTS_stage",
                        "all_canonical", "min_sample_cov", "within_CAGE_peak", "dist_to_CAGE_peak", "perc_A_downstream_TTS",
                        "polyA_motif_found"] if c in sq.columns]
    sq = sq[keep].rename(columns={"isoform": "transcript_id"})
    ev = ev.merge(sq, on="transcript_id", how="left")

# ---------------- loci (units) ----------------
gb = pd.read_csv(REF + "/gene_boundaries_by_stage.tsv", sep="\t", low_memory=False).set_index("v48_gene")
tpm = pd.read_csv(SR + "/counts/locus_TPM_9stages.tsv", sep="\t", index_col=0)
units = {}
for loc, L in ev.groupby("locus"):
    ch, sd = L.chrom.iloc[0], L.strand.iloc[0]
    if loc in v48 and loc in gb.index and not pd.isna(gb.loc[loc, "dominant_tss"]):
        tss, tes = int(gb.loc[loc, "dominant_tss"]), int(gb.loc[loc, "dominant_tes"])
    else:
        a = L.sort_values(["tss_cap_FL", "n_reads_FLweighted"], ascending=False).iloc[0]
        b = L.sort_values(["tes_reads", "n_reads_FLweighted"], ascending=False).iloc[0]
        tss, tes = int(a.tss), int(b.tes)
    if loc in v48:
        cls = gtype.get(loc, "mRNA")
    else:
        orf = max(longest_orf(ch, sd, exons[t]) for t in L.transcript_id)
        base = "lncRNA_antisense" if (L["class"] == "antisense_to_v48").any() else "lncRNA_intergenic"
        cls = base if orf < MIN_ORF_AA else "novel_ORF_" + base.split("_")[1]
    lo, hi = min(tss, tes), max(tss, tes) + 1
    if hi - lo < 50: lo, hi = min(L.start), max(L.end)
    u = dict(locus=loc, chrom=ch, strand=sd, tss=tss, tes=tes, lo=lo, hi=hi, cls=cls, lr_validated=1,
             readthrough=int((L.n_v48_genes_spanned >= 2).any()), n_models=len(L),
             cds_s=v48.get(loc, {}).get("cds_s"), cds_e=v48.get(loc, {}).get("cds_e"),
             max_orf_aa=(None if loc in v48 else orf) if loc not in v48 else None,
             desc=gdesc.get(loc, ""))
    units[loc] = u
for gid, g in v48.items():                                   # v48 genes without HC models (annotation only)
    if gid in units: continue
    tss = g["start"] if g["strand"] == "+" else g["end"] - 1
    tes = g["end"] - 1 if g["strand"] == "+" else g["start"]
    units[gid] = dict(locus=gid, chrom=g["chrom"], strand=g["strand"], tss=tss, tes=tes, lo=g["start"], hi=g["end"],
                      cls=gtype.get(gid, "mRNA"), lr_validated=0, readthrough=0, n_models=0,
                      cds_s=g["cds_s"], cds_e=g["cds_e"], max_orf_aa=None, desc=gdesc.get(gid, ""))
byc = collections.defaultdict(list)
for u in units.values(): byc[u["chrom"]].append(u)
for ch in byc: byc[ch].sort(key=lambda u: u["lo"])
starts = {ch: [u["lo"] for u in L] for ch, L in byc.items()}
maxlen = max(u["hi"] - u["lo"] for u in units.values())

def overlapping(ch, p):
    L = byc.get(ch, []); i = bisect.bisect_right(starts.get(ch, []), p)
    out = []
    j = i - 1
    while j >= 0 and L[j]["lo"] > p - maxlen:
        if L[j]["lo"] <= p < L[j]["hi"]: out.append(L[j])
        j -= 1
    return out

def region_in(u, p):
    """position of p inside host u, in host orientation"""
    if u["cls"] == "mRNA" and u["cds_s"] is not None:
        cs, ce = u["cds_s"], u["cds_e"]
        if cs <= p < ce: return "CDS"
        before = p < cs
        return ("5'UTR" if before else "3'UTR") if u["strand"] == "+" else ("3'UTR" if before else "5'UTR")
    L = u["hi"] - u["lo"]; f = (p - u["lo"]) / max(L, 1)
    if u["strand"] == "-": f = 1 - f
    return "5'end" if f < 0.2 else ("3'end" if f > 0.8 else "body")

CLS_RANK = {"mRNA": 0, "lncRNA_antisense": 1, "lncRNA_intergenic": 1, "novel_ORF_antisense": 1, "novel_ORF_intergenic": 1}
def host_label(c):
    return "mRNA" if c == "mRNA" else ("lncRNA" if c.startswith(("lncRNA", "novel_ORF")) else "ncRNA")

def neighbour(ch, p, direction):
    """nearest unit entirely on one side of p: direction=-1 -> units ending before p; +1 -> starting after p"""
    L = byc.get(ch, []); best = None; bd = None
    i = bisect.bisect_right(starts.get(ch, []), p)
    if direction > 0:
        for u in L[i:i + 60]:
            if u["lo"] > p:
                d = u["lo"] - p
                if bd is None or d < bd: best, bd = u, d
    else:
        for u in L[max(0, i - 200):i]:
            if u["hi"] <= p:
                d = p - u["hi"] + 1
                if bd is None or d < bd: best, bd = u, d
    return best, bd

def dist_bin(d):
    return "≤500bp" if d <= 500 else ("0.5-1kb" if d <= 1000 else ("1-2kb" if d <= 2000 else ">2kb"))

def classify(owner, p, end):
    sd, ch = owner["strand"], owner["chrom"]
    hosts = [u for u in overlapping(ch, p) if u["locus"] != owner["locus"]]
    if hosts:
        hosts.sort(key=lambda u: (CLS_RANK.get(u["cls"], 2), -min(p - u["lo"], u["hi"] - p)))
        h = hosts[0]; rel = "sense" if h["strand"] == sd else "antisense"
        reg = region_in(h, p)
        return dict(context_group="embedded", relation=rel, host=h["locus"], host_class=h["cls"],
                    host_region=reg, n_hosts=len(hosts), arrangement="", neighbour="", neighbour_dist=np.nan,
                    context=f"in {rel} {host_label(h['cls'])} {reg}")
    # proximal: within PROX bp upstream of another unit's TSS or downstream of its TES
    prox = []
    for u in byc.get(ch, []):
        pass
    cand = []
    L = byc.get(ch, []); i = bisect.bisect_right(starts.get(ch, []), p)
    for u in L[max(0, i - 200): i + 60]:
        if u["locus"] == owner["locus"]: continue
        up = (u["tss"] - p) if u["strand"] == "+" else (p - u["tss"])      # >0: p upstream of u TSS
        dn = (p - u["tes"]) if u["strand"] == "+" else (u["tes"] - p)      # >0: p downstream of u TES
        if 0 < up <= PROX: cand.append((up, u, "promoter-proximal"))
        if 0 < dn <= PROX: cand.append((dn, u, "terminator-proximal"))
    # intergenic arrangement relative to nearest neighbour on the relevant side
    if end == "TSS":
        side = -1 if sd == "+" else 1           # upstream side
        nb, d = neighbour(ch, p, side)
        arr = "" if nb is None else ("divergent" if nb["strand"] != sd else "tandem(after upstream TES)")
    else:
        side = 1 if sd == "+" else -1           # downstream side
        nb, d = neighbour(ch, p, side)
        arr = "" if nb is None else ("convergent" if nb["strand"] != sd else "tandem(before downstream TSS)")
    base = dict(arrangement=arr, neighbour=nb["locus"] if nb else "", neighbour_dist=d if nb else np.nan,
                n_hosts=0, host_region="")
    if cand:
        cand.sort(key=lambda x: x[0]); d2, u, kind = cand[0]
        rel = "sense" if u["strand"] == sd else "antisense"
        base.update(context_group="proximal", relation=rel, host=u["locus"], host_class=u["cls"],
                    context=f"{kind} ({rel} {host_label(u['cls'])}, ≤{PROX}bp)")
        return base
    base.update(context_group="intergenic", relation="", host="", host_class="",
                context=f"intergenic {arr} {dist_bin(d)}" if nb else "intergenic (no neighbour)")
    return base

# ---------------- ends ----------------
tsscl = pd.read_csv(REF + "/tss_clusters.tsv", sep="\t", low_memory=False)
cagec = [c for c in tsscl.columns if c.startswith("cage_")]
rows = []
for loc, u in units.items():
    for end in ("TSS", "TTS"):
        p = u["tss"] if end == "TSS" else u["tes"]
        c = classify(u, p, end)
        rows.append(dict(end_id=f"{loc}|{end}|dominant", locus=loc, end_type=end, end_role="dominant",
                         owner_class=u["cls"], lr_validated=u["lr_validated"], readthrough=u["readthrough"],
                         chrom=u["chrom"], pos=p, strand=u["strand"], **c))
# alternative TSS clusters (not decay, >=2 FL or CAGE >=1 CPM in both reps of any stage)
for r in tsscl.itertuples():
    u = units.get(r.locus)
    if u is None or r.decay_candidate == 1 or abs(r.tss - u["tss"]) <= 20: continue
    fl = r.FL_ring + r.FL_troph + r.FL_schizont
    cage_ok = any(getattr(r, f"cage_{s}_A") >= 1 and getattr(r, f"cage_{s}_B") >= 1 for s in ("ring", "troph", "schizont"))
    if fl < 2 and not cage_ok: continue
    p = int(r.tss)
    inside_own = u["lo"] <= p < u["hi"]
    c = classify(u, p, "TSS")
    own = region_in(u, p) if inside_own else "upstream of dominant TSS"
    rows.append(dict(end_id=f"{r.locus}|TSS|{r.tss_cluster}", locus=r.locus, end_type="TSS", end_role="alternative",
                     owner_class=u["cls"], lr_validated=1, readthrough=u["readthrough"],
                     chrom=u["chrom"], pos=p, strand=u["strand"], own_locus_position=own, **c))
E = pd.DataFrame(rows)
E["own_locus_position"] = E["own_locus_position"].fillna("")
E["situation"] = E.end_type + ":" + E.owner_class.map(host_label) + " | " + E.context
E.loc[E.end_role == "alternative", "situation"] = "alt TSS (" + E.own_locus_position + ") | " + E.context

# ---------------- attach evidence / expression ----------------
t = tpm.reindex(E.locus)
for s in STAGES: E[f"TPM_{s}"] = t[s].values
E["tau"] = t["tau"].values; E["peak_stage"] = t["peak_stage"].values; E["max_TPM"] = t["max_TPM"].values
ht = tpm.reindex(E.host.replace("", np.nan))
E["host_max_TPM"] = ht["max_TPM"].values; E["host_peak_stage"] = ht["peak_stage"].values
# nearest TSS cluster CAGE/FL for TSS ends
tc = {k: g.sort_values("tss") for k, g in tsscl.groupby("locus")}
cage_cols = {s: [f"cage_{s}_A", f"cage_{s}_B"] for s in ("ring", "troph", "schizont")}
for s in ("ring", "troph", "schizont"):
    E[f"CAGE_{s}_cpm"] = np.nan; E[f"LR_FL_{s}"] = np.nan
for i, r in E[E.end_type == "TSS"].iterrows():
    g = tc.get(r.locus)
    if g is None: continue
    j = (g.tss - r.pos).abs().idxmin()
    if abs(g.loc[j, "tss"] - r.pos) > 20: continue
    for s in ("ring", "troph", "schizont"):
        E.at[i, f"CAGE_{s}_cpm"] = round(g.loc[j, cage_cols[s]].mean(), 3); E.at[i, f"LR_FL_{s}"] = g.loc[j, f"FL_{s}"]
# locus-level SQANTI / junction summaries
agg = ev.groupby("locus").agg(n_hc_models=("transcript_id", "size"),
                              frac_models_all_introns_sr=("all_introns_sr_supported", "mean"),
                              n_gold=("tier", lambda x: (x == "Gold").sum()))
if sq is not None and "structural_category" in ev.columns:
    agg = agg.join(ev.groupby("locus").structural_category.agg(lambda x: x.value_counts().index[0]).rename("main_sqanti_category"))
    if "coding" in ev.columns:
        agg = agg.join(ev.groupby("locus").coding.agg(lambda x: (x == "coding").mean()).rename("frac_sqanti_coding"))
E = E.merge(agg, left_on="locus", right_index=True, how="left")
E["description"] = E.locus.map(lambda x: units[x]["desc"])
E.to_csv(OUT + "/end_contexts.tsv", sep="\t", index=False)
U = pd.DataFrame(units.values()).drop(columns=["lo", "hi"]).rename(columns={"cls": "transcript_class"})
U.to_csv(OUT + "/loci_classes.tsv", sep="\t", index=False)

# ---------------- BED per situation (validated ends only) ----------------
def safe(s): return re.sub(r"[^A-Za-z0-9]+", "_", s).strip("_")[:90]
V = E[E.lr_validated == 1]
for sit, g in V.groupby("situation"):
    if len(g) < 20: continue
    with open(f"{OUT}/bed/{safe(sit)}.bed", "w") as fh:
        for r in g.sort_values(["chrom", "pos"]).itertuples():
            fh.write(f"{r.chrom}\t{r.pos}\t{r.pos+1}\t{r.end_id}\t0\t{r.strand}\n")
summ = dict(n_units=len(units), unit_classes=U.transcript_class.value_counts().to_dict(),
            lr_validated_units=int(U.lr_validated.sum()),
            n_ends=len(E), n_ends_validated=len(V),
            by_situation=V.situation.value_counts().to_dict(),
            by_context_group=V.groupby(["end_type", "end_role", "context_group"]).size().rename("n").reset_index().values.tolist(),
            hc_models_all_introns_sr=float(jt.all_introns_sr_supported.mean()),
            hc_multiexon_all_introns_sr=float(jt[jt.n_introns > 0].all_introns_sr_supported.mean()),
            sqanti_available=sq is not None)
json.dump(summ, open(OUT + "/summary.json", "w"), indent=1, default=str)
print(json.dumps({k: v for k, v in summ.items() if k != "by_situation"}, indent=1, default=str))
print("\n".join(f"{v}\t{k}" for k, v in list(summ["by_situation"].items())[:60]))
