#!/bin/bash
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
# StringTie baselines for the boundary benchmark (de novo, no reference guide). Outputs in ${PF_HPCC_BIGDATA}/benchmark_2026/stringtie
set -e
module load stringtie/3.0.3 2>/dev/null || true
PKG=/opt/linux/rocky/8.x/x86_64/pkgs; export LD_LIBRARY_PATH=$PKG/htslib/1.19.1/lib:$PKG/libdeflate/1.10/lib:$PKG/isa-l/8b7c1b8/lib:$LD_LIBRARY_PATH; ST=$PKG/samtools/1.19.2/bin/samtools
O=${PF_HPCC_BIGDATA}/benchmark_2026/stringtie; mkdir -p $O/short $O/long $O/mix
D=${PF_HPCC_BIGDATA}/flank_coverage_2026/bam; G=${PF_HPCC_BIGDATA}/GSE186109_longread/bam; R=${PF_HPCC_BIGDATA}/longread_reprocess_2026
stringtie --version
# 1 short reads: stranded DAFT-seq (reverse-stranded library -> --rf)
for b in $D/*.bam; do n=$(basename $b .bam); [ -s $O/short/$n.gtf ] || stringtie $b --rf -p 2 -m 150 -o $O/short/$n.gtf -l DAFT_$n; done
ls $O/short/*.gtf > $O/short/list.txt; stringtie --merge -p 2 -o $O/ST_short_merged.gtf -l STs $O/short/list.txt
# 2 long reads
LR="$G/ring_total/hq.sorted.bam $G/ring_nuclear/hq.sorted.bam $G/troph_total/hq.sorted.bam $G/troph_nuclear/hq.sorted.bam $G/schizont_total/hq.sorted.bam $G/schizont_nuclear/hq.sorted.bam $G/ont_mixed_oriented/oriented.sorted.bam $R/ont/LeeVV_2021_dRNA.sorted.bam $R/ont/Hoshizaki_2022_dRNA.sorted.bam $(ls $R/pacbio/*_flnc_all.sorted.bam)"
for b in $LR; do [ -s $b ] || continue; n=$(echo $b | sed 's#.*/bam/##; s#.*/ont/##; s#.*/pacbio/##; s#/#_#g; s#.sorted.bam##'); [ -s $O/long/$n.gtf ] || stringtie -L $b -p 2 -o $O/long/$n.gtf -l LR_$n; done
ls $O/long/*.gtf > $O/long/list.txt; stringtie --merge -p 2 -o $O/ST_long_merged.gtf -l STl $O/long/list.txt
# 3 mixed: merged DAFT + merged long reads
[ -s $O/mix/daft_all.bam ] || $ST merge -@2 -f $O/mix/daft_all.bam $D/*.bam
[ -s $O/mix/long_all.bam ] || $ST merge -@2 -f $O/mix/long_all.bam $LR
$ST index $O/mix/daft_all.bam; $ST index $O/mix/long_all.bam
[ -s $O/ST_mix.gtf ] || stringtie --mix --rf -p 2 -o $O/ST_mix.gtf -l STm $O/mix/daft_all.bam $O/mix/long_all.bam
rm -f $O/mix/*.bam $O/mix/*.bai
ls -la $O/*.gtf; echo ALLDONE
