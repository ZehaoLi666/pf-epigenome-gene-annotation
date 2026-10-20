from pfconfig import HPCC_BIGDATA  # paths: config/paths.yaml

import csv, glob, os
from collections import defaultdict

PROJ = f"{HPCC_BIGDATA}/GSE186109_longread"
manifest = os.path.join(PROJ, "bed12_manifest.csv")

rows = []
with open(manifest) as f:
    r = csv.DictReader(f)
    for row in r:
        rows.append(row)

# ends[(chrom,strand,endtype)] = list of (pos, sample, tier, platform, read_len)
ends = defaultdict(list)

for row in rows:
    sample, tier, platform, bedpath = row["sample"], row["tier"], row["platform"], row["bed12"]
    with open(bedpath) as bf:
        for line in bf:
            f_ = line.rstrip("\n").split("\t")
            chrom, start, end, name, score, strand = f_[0], int(f_[1]), int(f_[2]), f_[3], f_[4], f_[5]
            read_len = end - start
            if strand == "+":
                five, three = start, end - 1
            else:
                five, three = end - 1, start
            ends[(chrom, strand, "5p")].append((five, sample, tier, platform))
            ends[(chrom, strand, "3p")].append((three, sample, tier, platform))

print("total keys:", len(ends))
total_pts = sum(len(v) for v in ends.values())
print("total end points:", total_pts)

# greedy clustering within WINDOW bp
WINDOW = 20
clusters = []  # chrom, strand, endtype, cl_min, cl_max, mode_pos, n_reads, n_distinct_samples, samples_set

for (chrom, strand, endtype), pts in ends.items():
    pts.sort(key=lambda x: x[0])
    i = 0
    n = len(pts)
    while i < n:
        j = i
        cl_positions = [pts[i][0]]
        cl_samples = set()
        cl_samples.add((pts[i][1], pts[i][2]))
        while j + 1 < n and pts[j+1][0] - cl_positions[-1] <= WINDOW:
            j += 1
            cl_positions.append(pts[j][0])
            cl_samples.add((pts[j][1], pts[j][2]))
        # mode = most frequent position
        from collections import Counter
        cnt = Counter(cl_positions)
        mode_pos, mode_count = cnt.most_common(1)[0]
        clusters.append({
            "chrom": chrom, "strand": strand, "endtype": endtype,
            "cl_min": cl_positions[0], "cl_max": cl_positions[-1],
            "mode_pos": mode_pos, "n_reads": len(cl_positions),
            "n_distinct_sampletier": len(cl_samples),
            "samples": ";".join(sorted(f"{s}:{t}" for s,t in cl_samples))
        })
        i = j + 1

print("total clusters:", len(clusters))
out = os.path.join(PROJ, "readend_atlas.tsv")
with open(out, "w") as f:
    hdr = ["chrom","strand","endtype","cl_min","cl_max","mode_pos","n_reads","n_distinct_sampletier","samples"]
    f.write("\t".join(hdr) + "\n")
    for c in clusters:
        f.write("\t".join(str(c[h]) for h in hdr) + "\n")
print("wrote", out)
