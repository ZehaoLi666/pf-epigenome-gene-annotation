#!/usr/bin/env python3
import os, shutil, subprocess, shlex

FEATURE_FILES = {
    # "ATAC-seq": ".../ATAC-seq.txt",
    # "H3R17me2": ".../H3R17me2.txt",
    # ...
    "MNase":  "/rhome/zli529/lab/chip-seq/chipseq/MNase-seq_results/merge_bam/MNseq_sub.txt",
    "H3K9ac": "/rhome/zli529/lab/chip-seq/chipseq/H3K9ac_results/merge_bam/H3K9ac.sub.txt",
    "H3K36me3": "/rhome/zli529/lab/chip-seq/chipseq/H3K36me3_results/H3K36me3.Nor.txt",
}

FAI   = "/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta.fai"
OUTDIR= "/rhome/zli529/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw"
os.makedirs(OUTDIR, exist_ok=True)
chrom_sizes = os.path.join(OUTDIR, "genome.chrom.sizes")

def run(cmd):
    print("[cmd]", cmd)
    subprocess.run(cmd, shell=True, check=True)

# --- build chrom.sizes from .fai ---
with open(FAI) as fin, open(chrom_sizes, "w") as fout:
    for line in fin:
        toks = line.rstrip("\n").split("\t")
        if len(toks) >= 2:
            fout.write(toks[0] + "\t" + toks[1] + "\n")

# --- check UCSC tools ---
if not shutil.which("bedGraphToBigWig"):
    raise SystemExit("ERROR: bedGraphToBigWig not found in PATH. Try: conda install -c bioconda ucsc-bedgraphtobigwig")

# --- conversion + cleaning awk ---
# Notes:
#  * loads chrom sizes (ok[], len[])
#  * skips header/blank/NA/NaN/inf/non-numeric
#  * renames mito to v48 name (Pf3D7_MIT_v3 -> Pf_M76611; also M76611 -> Pf_M76611)
#  * enforces 4 fields and clamps to chrom length
AWK_CLEAN = r"""
awk -v OFS='\t' '
FNR==NR {len[$1]=$2; ok[$1]=1; next}
NR==1   { if($2 !~ /^[0-9]+$/ || $3 !~ /^-?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) next }
{
  chr=$1; pos=$2; val=$3;
  # rename mito to v48 alias
  if(chr=="Pf3D7_MIT_v3") chr="Pf_M76611";
  if(chr=="M76611")       chr="Pf_M76611";
  # keep only chroms present in chrom.sizes
  if(!ok[chr]) next;
  # numeric checks
  if(pos !~ /^[0-9]+$/) next;
  if(val=="NA" || val=="NaN" || val=="nan" || val=="inf" || val=="-inf") next;
  if(val !~ /^-?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) next;
  # bedGraph [start,end) (0-based start)
  start=pos-1; if(start<0) start=0; end=pos;
  if(end>len[chr]) end=len[chr];
  if(start<end) print chr, start, end, val;
}'
"""

for name, txt in FEATURE_FILES.items():
    if not os.path.exists(txt) or os.stat(txt).st_size == 0:
        print(f"[warn] Skipping {name} (missing/empty): {txt}")
        continue

    bg  = os.path.join(OUTDIR, f"{name}.bedgraph")
    bgs = os.path.join(OUTDIR, f"{name}.sorted.bedgraph")
    bw  = os.path.join(OUTDIR, f"{name}.bw")

    # 1) txt -> cleaned bedGraph (4 columns guaranteed)
    run(f"{AWK_CLEAN} {shlex.quote(chrom_sizes)} {shlex.quote(txt)} > {shlex.quote(bg)}")

    # 2) sort
    run(f"LC_ALL=C sort -k1,1 -k2,2n {shlex.quote(bg)} -o {shlex.quote(bgs)}")

    # 2.5) validate: every line must have 4 fields
    run(f"awk 'NF!=4{{bad++}} END{{if(bad>0){{print \"[error] bad lines:\",bad > \"/dev/stderr\"; exit 1}}}}' {shlex.quote(bgs)}")

    # 3) convert to bigWig
    run(f"bedGraphToBigWig {shlex.quote(bgs)} {shlex.quote(chrom_sizes)} {shlex.quote(bw)}")

    # 4) cleanup (optional)
    try:
        os.remove(bg)
        os.remove(bgs)
    except OSError:
        pass

    print(f"[✓] Wrote {bw}")

print(f"All done. Output dir: {OUTDIR}")
