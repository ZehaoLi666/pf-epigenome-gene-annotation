#!/bin/bash
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
set -uo pipefail
PKG=/opt/linux/rocky/8.x/x86_64/pkgs
SAMTOOLS=$PKG/samtools/1.22.1/bin/samtools
BEDTOOLS=$PKG/bedtools/2.30.0/bin/bedtools
export LD_LIBRARY_PATH=$PKG/htslib/1.22.1/lib:$PKG/libdeflate/1.17/lib:${LD_LIBRARY_PATH:-}
PROJ=${PF_HPCC_BIGDATA}/GSE186109_longread
OUT=$PROJ/bed12
mkdir -p $OUT

echo "sample,tier,platform,bed12" > $PROJ/bed12_manifest.csv
while IFS=$'\t' read -r sample platform tier fastq; do
  bam=$PROJ/bam/$sample/${tier}.sorted.bam
  bed=$OUT/${sample}__${tier}.bed12
  # -F 0x900 excludes secondary(0x100) + supplementary(0x800); primary mapped reads only
  "$SAMTOOLS" view -b -F 0x900 "$bam" | "$BEDTOOLS" bamtobed -bed12 -i - > "$bed"
  n=$(wc -l < "$bed")
  echo "$sample,$tier,$platform,$bed" >> $PROJ/bed12_manifest.csv
  echo "$sample/$tier: $n primary alignments -> $bed"
done < $PROJ/scripts/align_tasks.tsv
