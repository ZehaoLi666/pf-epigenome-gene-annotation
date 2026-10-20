#!/bin/bash
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
set -uo pipefail
PROJ=${PF_HPCC_BIGDATA}/GSE186109_longread
FQ=$PROJ/fastq
BY=$PROJ/by_sample
mkdir -p $BY
MANIFEST=$PROJ/manifest_isoform_tier.csv
echo "sample,stage,fraction,platform,tier,srr,n_reads,path" > $MANIFEST

declare -A SAMPLES=(
  [schizont_total]="SRR16480181:hq SRR16480182:lq"
  [ring_total]="SRR16480185:hq SRR16480186:lq"
  [troph_total]="SRR16480190:hq SRR16480191:lq"
  [ring_nuclear]="SRR16480195:hq SRR16480196:lq"
  [troph_nuclear]="SRR16480199:hq SRR16480200:lq"
  [schizont_nuclear]="SRR16480203:hq SRR16480204:lq"
  [schiztroph_cytoplasmic]="SRR16480206:hq SRR16480207:lq"
  [ont_mixed]="SRR18188665:main"
)

for sample in "${!SAMPLES[@]}"; do
  mkdir -p $BY/$sample
  for pair in ${SAMPLES[$sample]}; do
    srr="${pair%%:*}"; tier="${pair##*:}"
    src="$FQ/$srr.fastq"
    dst="$BY/$sample/${tier}.fastq"
    if [ -s "$src" ]; then
      ln -sf "$src" "$dst"
      n=$(( $(wc -l < "$src") / 4 ))
      plat="pacbio"; [[ "$sample" == "ont_mixed" ]] && plat="ont"
      stage=$(echo $sample | sed 's/_.*//')
      frac=$(echo $sample | sed 's/^[a-z]*_//')
      echo "$sample,$stage,$frac,$plat,$tier,$srr,$n,$dst" >> $MANIFEST
    else
      echo "$sample,,,,${tier},$srr,MISSING,$src" >> $MANIFEST
    fi
  done
done
cat $MANIFEST
