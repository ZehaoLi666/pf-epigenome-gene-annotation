from pfconfig import HPCC_BIGDATA  # paths: config/paths.yaml

import csv, statistics
PROJ = f"{HPCC_BIGDATA}/GSE186109_longread"
path = PROJ + "/gene_vs_readend_v48.tsv"
tss_offsets = []
tts_offsets = []
tss_multi = 0
tts_multi = 0
n = 0
with open(path) as f:
    r = csv.DictReader(f, delimiter="\t")
    for row in r:
        n += 1
        if row["tss_offset"] != "":
            tss_offsets.append(abs(int(row["tss_offset"])))
            if int(row["tss_nsamples"]) >= 2:
                tss_multi += 1
        if row["tts_offset"] != "":
            tts_offsets.append(abs(int(row["tts_offset"])))
            if int(row["tts_nsamples"]) >= 2:
                tts_multi += 1

def summ(vals, label):
    vals = sorted(vals)
    print(f"{label}: n={len(vals)} median={statistics.median(vals):.0f} "
          f"mean={statistics.mean(vals):.1f} "
          f"p25={vals[len(vals)//4]} p75={vals[3*len(vals)//4]} "
          f"within50bp={sum(1 for v in vals if v<=50)} "
          f"within200bp={sum(1 for v in vals if v<=200)}")

summ(tss_offsets, "5' offset (bp)")
summ(tts_offsets, "3' offset (bp)")
print(f"5' supported by >=2 distinct sample/tier combos: {tss_multi}/{len(tss_offsets)} ({100*tss_multi/len(tss_offsets):.1f}%)")
print(f"3' supported by >=2 distinct sample/tier combos: {tts_multi}/{len(tts_offsets)} ({100*tts_multi/len(tts_offsets):.1f}%)")
print(f"total genes: {n}")
