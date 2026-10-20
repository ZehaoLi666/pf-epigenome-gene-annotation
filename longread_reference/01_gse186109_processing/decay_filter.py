from pfconfig import HPCC_BIGDATA, HPCC_LAB  # paths: config/paths.yaml

import csv, re, os, bisect
from collections import defaultdict

PROJ = f"{HPCC_BIGDATA}/GSE186109_longread"
GTF = f"{HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7.gtf"
ATLAS = os.path.join(PROJ, "readend_atlas.tsv")

genes = []
with open(GTF) as f:
    for line in f:
        f_ = line.rstrip("\n").split("\t")
        if len(f_) < 9 or f_[2] != "gene":
            continue
        chrom, start, end, strand, attr = f_[0], int(f_[3]), int(f_[4]), f_[6], f_[8]
        if not (chrom.startswith("Pf3D7_") and "_v3" in chrom):
            continue
        m = re.search(r'gene_id "([^"]+)"', attr)
        genes.append({"gene_id": m.group(1) if m else "NA", "chrom": chrom, "start": start, "end": end, "strand": strand})

atlas = defaultdict(list)
with open(ATLAS) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        key = (row["chrom"], row["strand"], row["endtype"])
        atlas[key].append({"mode_pos": int(row["mode_pos"]), "n_reads": int(row["n_reads"]),
                            "n_distinct_sampletier": int(row["n_distinct_sampletier"])})
for k in atlas:
    atlas[k].sort(key=lambda c: c["mode_pos"])

MAXD_PRIMARY = 1000   # window to find the primary (annotation-matching) TSS
WINDOW_ALT = None      # search alt TSS across the gene body (start..end, extended -1kb upstream)

def clusters_in_range(chrom, strand, endtype, lo, hi):
    cl_list = atlas.get((chrom, strand, endtype), [])
    positions = [c["mode_pos"] for c in cl_list]
    i0 = bisect.bisect_left(positions, lo)
    i1 = bisect.bisect_right(positions, hi)
    return cl_list[i0:i1]

primary_alt_counts = {"n_genes_with_alt": 0, "alt_total": 0, "alt_reproducible": 0, "alt_single_library": 0}
gene_rows = []

for g in genes:
    chrom, strand = g["chrom"], g["strand"]
    if strand == "+":
        tss_anchor = g["start"]
        search_lo, search_hi = g["start"] - 1000, g["start"] + int(0.5*(g["end"]-g["start"]))
    else:
        tss_anchor = g["end"]
        search_lo, search_hi = g["end"] - int(0.5*(g["end"]-g["start"])), g["end"] + 1000

    cands = clusters_in_range(chrom, strand, "5p", search_lo, search_hi)
    if not cands:
        continue
    # primary = closest to tss_anchor
    primary = min(cands, key=lambda c: abs(c["mode_pos"] - tss_anchor))
    alts = [c for c in cands if c is not primary]
    if alts:
        primary_alt_counts["n_genes_with_alt"] += 1
        primary_alt_counts["alt_total"] += len(alts)
        for a in alts:
            if a["n_distinct_sampletier"] >= 2:
                primary_alt_counts["alt_reproducible"] += 1
                cls = "likely_genuine_alt_TSS"
            else:
                primary_alt_counts["alt_single_library"] += 1
                cls = "likely_decay_fragment"
            gene_rows.append({
                "gene_id": g["gene_id"], "chrom": chrom, "strand": strand,
                "primary_pos": primary["mode_pos"], "alt_pos": a["mode_pos"],
                "alt_offset_from_primary": a["mode_pos"] - primary["mode_pos"],
                "alt_n_reads": a["n_reads"], "alt_n_distinct_sampletier": a["n_distinct_sampletier"],
                "classification": cls
            })

print("genes with >=1 alternative 5' cluster (beyond primary):", primary_alt_counts["n_genes_with_alt"])
print("total alternative-5'-cluster instances:", primary_alt_counts["alt_total"])
print(f"  classified likely_genuine (>=2 independent sample/tier): {primary_alt_counts['alt_reproducible']} "
      f"({100*primary_alt_counts['alt_reproducible']/max(primary_alt_counts['alt_total'],1):.1f}%)")
print(f"  classified likely_decay_fragment (single sample/tier only): {primary_alt_counts['alt_single_library']} "
      f"({100*primary_alt_counts['alt_single_library']/max(primary_alt_counts['alt_total'],1):.1f}%)")

out = os.path.join(PROJ, "alt_tss_classification.tsv")
with open(out, "w") as f:
    if gene_rows:
        hdr = list(gene_rows[0].keys())
        f.write("\t".join(hdr) + "\n")
        for row in gene_rows:
            f.write("\t".join(str(row[h]) for h in hdr) + "\n")
print("wrote", out, "rows:", len(gene_rows))
