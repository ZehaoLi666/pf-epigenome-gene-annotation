"""Export TAMA collapsed transcript models (span, strand, block structure, supporting read count) for the four datasets."""
from pfconfig import HPCC_BIGDATA, HPCC_LAB  # paths: config/paths.yaml
import collections, glob, gzip
B = f"{HPCC_LAB}/SRA_toolkit/Long_read_RNA_datasets"; G = f"{HPCC_BIGDATA}/GSE186109_longread/tama"
SRC = {"LeeVV2021_dRNA": [(f"{B}/LeeVV_2021_nanopore/LeeVV_ONT.bed", f"{B}/LeeVV_2021_nanopore/LeeVV_ONT_trans_read.bed")],
       "Hoshizaki2022_dRNA": [(f"{B}/Marcus_Nanopore/Marcus_ONT.bed", f"{B}/Marcus_Nanopore/Marcus_ONT_trans_read.bed")],
       "Yang2021_IsoSeq": [(f"{B}/Mengquan_long_read/Mengquan_PacBio.bed", f"{B}/Mengquan_long_read/Mengquan_PacBio_trans_read.bed")]}
for st in ("ring", "troph", "schizont"):
    SRC[f"GSE186109_PacBio_{st}"] = [(f"{G}/collapse_up/GSE_{st}.bed", f"{G}/collapse_up/GSE_{st}_trans_read.bed")]
SRC["GSE186109_ONT_cDNA"] = [(b, b.replace(".bed", "_trans_read.bed")) for b in sorted(glob.glob(f"{G}/collapse_ont_chr/c*.bed")) if not b.endswith("_trans_read.bed")]
out = gzip.open("tama_models.tsv.gz", "wt"); out.write("dataset\tmodel\tchrom\tstrand\tstart\tend\tn_reads\tblock_sizes\tblock_starts\n")
for ds, pairs in SRC.items():
    n = 0
    for bed, tr in pairs:
        cnt = collections.Counter()
        try:
            for line in open(tr): cnt[line.split("\t")[3].split(";")[0]] += 1
        except FileNotFoundError: pass
        for line in open(bed):
            f = line.rstrip("\n").split("\t")
            if len(f) < 12: continue
            tid = f[3].split(";")[1] if ";" in f[3] else f[3]
            out.write(f"{ds}\t{tid}\t{f[0]}\t{f[5]}\t{f[1]}\t{f[2]}\t{cnt.get(tid, 0)}\t{f[10]}\t{f[11]}\n"); n += 1
    print(ds, len(pairs), n, flush=True)
out.close()
