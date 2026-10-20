#!/bin/bash
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
P=${PF_HPCC_BIGDATA}/GSE186109_longread
OUT=$P/isoform_fl_counts.tsv
echo -e "sample\ttier\tread\tcluster\tn_fl\tn_partial" > $OUT
tail -n +2 $P/manifest_isoform_tier.csv | while IFS=, read sample stage frac plat tier srr n path; do
  [ "$plat" == "pacbio" ] || continue
  awk -v s=$sample -v t=$tier 'NR%4==1{split($2,a,"/"); fl=a[2]; sub(/^f/,"",fl); split(fl,b,"p"); print s"\t"t"\t"substr($1,2)"\t"a[1]"\t"b[1]"\t"b[2]}' $path >> $OUT
done
awk -F'\t' 'NR>1{n[$1"/"$2]++; f[$1"/"$2]+=$5; p[$1"/"$2]+=$6} END{for(k in n) print k, n[k], f[k], p[k]}' $OUT | sort
echo FL_DONE
