#!/bin/bash
#SBATCH --job-name=lrcov2
#SBATCH --account=lerochlab
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=6:00:00
#SBATCH --output=logs/lrcov_%j.out
#SBATCH --error=logs/lrcov_%j.err
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
set -eu
PKG=/opt/linux/rocky/8.x/x86_64/pkgs
export LD_LIBRARY_PATH=$PKG/htslib/1.22.1/lib:$PKG/libdeflate/1.17/lib:${LD_LIBRARY_PATH:-}
ST=$PKG/samtools/1.22.1/bin/samtools; BT=$PKG/bedtools/2.30.0/bin/bedtools
$ST --version >/dev/null; $BT --version >/dev/null
W=${PF_HPCC_BIGDATA}/GSE186109_longread/lr_coverage; cd $W
export TMPDIR=$W/tmp; mkdir -p tmp out2; : > out2/nlines.tsv
G=$W/in/genome.fai
GS=${PF_HPCC_BIGDATA}/GSE186109_longread/bam
LR=${PF_HPCC_LAB}/SRA_toolkit/Long_read_RNA_datasets
declare -A B
B[Shaw2022_PacBio]="$(ls $GS/{ring,troph,schizont}_{total,nuclear}/{hq,lq}.sorted.bam $GS/schiztroph_cytoplasmic/{hq,lq}.sorted.bam | tr '\n' ' ')"
B[Shaw2022_ONT]="$GS/ont_mixed_oriented/oriented.sorted.bam"
B[Lee2021]="$LR/LeeVV_2021_nanopore/SRR11094274.no_long_introns.sorted.bam"
B[Hoshizaki2022]="$LR/Marcus_Nanopore/ERR9746080.no_long_introns.sorted.bam"
B[Yang2021]="$LR/Mengquan_long_read/merged.no_long_introns.sorted.bam"
export LC_ALL=C
for s in Shaw2022_PacBio Shaw2022_ONT Lee2021 Hoshizaki2022 Yang2021; do
  echo "== $s $(date)"
  : > tmp/$s.b12
  for b in ${B[$s]}; do $ST view -u -F 0x904 $b | $BT bamtobed -bed12 -i stdin >> tmp/$s.b12; done
  sort -k1,1 -k2,2n -S 8G -T tmp tmp/$s.b12 > tmp/$s.sorted.b12; rm tmp/$s.b12
  echo -e "$s\t$(wc -l < tmp/$s.sorted.b12)" >> out2/nlines.tsv
  # span of reads whose largest intron <= 2 kb
  awk 'BEGIN{OFS="\t"}{n=$10; split($11,sz,","); split($12,st,","); mg=0; for(i=2;i<=n;i++){g=st[i]-(st[i-1]+sz[i-1]); if(g>mg)mg=g}; if(mg<=2000) print $1,$2,$3,$4,$5,$6}' tmp/$s.sorted.b12 > tmp/$s.span.bed
  $BT intersect -a in/cds_env.bed -b tmp/$s.span.bed -s -f 1.0 -c > out2/$s.fullspan_counts.tsv
  rm tmp/$s.span.bed
done
cat tmp/Shaw2022_PacBio.sorted.b12 tmp/Shaw2022_ONT.sorted.b12 tmp/Lee2021.sorted.b12 tmp/Hoshizaki2022.sorted.b12 tmp/Yang2021.sorted.b12 > tmp/allu.b12
echo -e "ALL_unsorted\t$(wc -l < tmp/allu.b12)" >> out2/nlines.tsv
sort -k1,1 -k2,2n -S 16G -T tmp tmp/allu.b12 > tmp/all.b12; rm tmp/allu.b12
echo -e "ALL_sorted\t$(wc -l < tmp/all.b12)" >> out2/nlines.tsv
$BT coverage -a in/cds_blocks.bed -b tmp/all.b12 -s -split > out2/ALL.cds_blocks_cov.tsv
$BT coverage -a in/cds_env.bed -b tmp/all.b12 -s -split -counts > out2/ALL.gene_counts.tsv
$BT genomecov -i tmp/all.b12 -g $G -split -bg > tmp/all.bg
awk '$4>=1' tmp/all.bg | $BT merge -i stdin > out2/ALL.cov1.bed
awk '$4>=5' tmp/all.bg | $BT merge -i stdin > out2/ALL.cov5.bed
rm -f tmp/all.bg tmp/all.b12 tmp/*.sorted.b12
tar czf out2.tgz out2; md5sum out2.tgz > out2.md5
echo DONE
