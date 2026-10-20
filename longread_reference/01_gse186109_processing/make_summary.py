from pfconfig import HPCC_BIGDATA  # paths: config/paths.yaml

import csv, json, statistics
PROJ = f"{HPCC_BIGDATA}/GSE186109_longread"

def load_offsets(path, tss_col, tts_col):
    tss, tts = [], []
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            if row[tss_col] not in ("", None):
                tss.append(int(row[tss_col]))
            if row[tts_col] not in ("", None):
                tts.append(int(row[tts_col]))
    return tss, tts

v48_tss, v48_tts = load_offsets(PROJ+"/gene_vs_readend_v48.tsv", "tss_offset", "tts_offset")
ref_tss, ref_tts = load_offsets(PROJ+"/manuscript_ref_vs_readend.tsv", "tss_offset_new", "tts_offset_new")

def hist(vals, bins):
    absvals = [abs(v) for v in vals]
    counts = [0]*(len(bins))
    for v in absvals:
        for i, b in enumerate(bins):
            if v <= b:
                counts[i] += 1
                break
        else:
            counts[-1] += 1  # shouldn't happen since last bin is inf
    return counts

BINS = [10,25,50,100,200,500,1000]
BIN_LABELS = ["<=10","11-25","26-50","51-100","101-200","201-500","501-1000"]

out = {
    "v48_n_genes_total": 5670,
    "v48_tss_n_matched": len(v48_tss), "v48_tts_n_matched": len(v48_tts),
    "v48_tss_hist": hist(v48_tss, BINS), "v48_tts_hist": hist(v48_tts, BINS),
    "v48_tss_median": statistics.median([abs(v) for v in v48_tss]),
    "v48_tts_median": statistics.median([abs(v) for v in v48_tts]),
    "ref_n_genes_total": 2313,
    "ref_tss_n_matched": len(ref_tss), "ref_tts_n_matched": len(ref_tts),
    "ref_tss_hist": hist(ref_tss, BINS), "ref_tts_hist": hist(ref_tts, BINS),
    "bin_labels": BIN_LABELS,
    "decay_classification": {"likely_genuine_alt_TSS": 33717, "likely_decay_fragment": 62780},
    "single_source_promotion": {"n_single_source": 1019, "gaining_tss": 1013, "gaining_tts": 1011},
}
with open(PROJ+"/figure_summary.json", "w") as f:
    json.dump(out, f, indent=2)
print(json.dumps(out, indent=2))
