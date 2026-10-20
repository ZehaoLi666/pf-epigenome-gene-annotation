#!/bin/bash
#SBATCH --job-name=tama_col_up
#SBATCH --account=lerochlab
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=2
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --array=1-7
#SBATCH -o ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/tamaup_%A_%a.out
#SBATCH -e ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/tamaup_%A_%a.err
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
set -uo pipefail
PROJ=${PF_HPCC_BIGDATA}/GSE186109_longread
PKG=/opt/linux/rocky/8.x/x86_64/pkgs
SAMTOOLS=$PKG/samtools/1.22.1/bin/samtools
export LD_LIBRARY_PATH=$PKG/htslib/1.22.1/lib:$PKG/libdeflate/1.17/lib:${LD_LIBRARY_PATH:-}
export TMPDIR=$PROJ/tmp
PY=${PF_HPCC_BIGDATA}/.conda/envs/tama_py27/bin/python
$PY -c "import Bio, pysam" || { echo "tama env broken"; exit 1; }
"$SAMTOOLS" --version >/dev/null || { echo "samtools libs"; exit 1; }
TAMA=$PROJ/tama/tama_upstream/tama_collapse.py
GEN=${PF_HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta
LR=${PF_HPCC_LAB}/SRA_toolkit/Long_read_RNA_datasets
case $SLURM_ARRAY_TASK_ID in
 1) ST=GSE_ring;     CAP=capped; BAMS="$PROJ/bam/ring_total/hq.sorted.bam $PROJ/bam/ring_total/lq.sorted.bam $PROJ/bam/ring_nuclear/hq.sorted.bam $PROJ/bam/ring_nuclear/lq.sorted.bam";;
 2) ST=GSE_troph;    CAP=capped; BAMS="$PROJ/bam/troph_total/hq.sorted.bam $PROJ/bam/troph_total/lq.sorted.bam $PROJ/bam/troph_nuclear/hq.sorted.bam $PROJ/bam/troph_nuclear/lq.sorted.bam";;
 3) ST=GSE_schizont; CAP=capped; BAMS="$PROJ/bam/schizont_total/hq.sorted.bam $PROJ/bam/schizont_total/lq.sorted.bam $PROJ/bam/schizont_nuclear/hq.sorted.bam $PROJ/bam/schizont_nuclear/lq.sorted.bam";;
 4) ST=GSE_schiztroph_cyto; CAP=capped; BAMS="$PROJ/bam/schiztroph_cytoplasmic/hq.sorted.bam $PROJ/bam/schiztroph_cytoplasmic/lq.sorted.bam";;
 5) ST=LeeVV2021; CAP=no_cap; BAMS="$LR/LeeVV_2021_nanopore/SRR11094274.no_long_introns.sorted.bam";;
 6) ST=Marcus;    CAP=no_cap; BAMS="$LR/Marcus_Nanopore/ERR9746080.no_long_introns.sorted.bam";;
 7) ST=Mengquan;  CAP=no_cap; BAMS="$LR/Mengquan_long_read/merged.no_long_introns.sorted.bam";;
esac
OUT=$PROJ/tama/collapse_up; mkdir -p $OUT; cd $OUT
set -- $BAMS
if [ $# -gt 1 ]; then "$SAMTOOLS" merge -f -@ 2 -o $OUT/$ST.merged.bam $BAMS; IN=$OUT/$ST.merged.bam; else IN=$1; fi
"$SAMTOOLS" view -h -F 0x904 $IN | "$SAMTOOLS" sort -@ 2 -m 2G -T $OUT/$ST.tmp -O SAM -o $OUT/$ST.pri.sam -
rm -f $OUT/$ST.merged.bam
echo "[$(date)] $ST input_primary=$(grep -vc '^@' $OUT/$ST.pri.sam)"
$PY $TAMA -s $OUT/$ST.pri.sam -b SAM -f $GEN -p $ST -x $CAP -e common_ends -a 20 -z 20 -c 99 -i 85 > $OUT/$ST.collapse.log 2>&1
RC=$?
echo "[$(date)] $ST rc=$RC tx=$(wc -l < $OUT/$ST.bed 2>/dev/null) accepted=$(awk '$3==\"accepted\"' $OUT/${ST}_read.txt | wc -l)"
[ $RC -eq 0 ] && rm -f $OUT/$ST.pri.sam
awk '{n++; if($6=="+")p++; split($11,s,","); split($12,st,","); k=$10; if(st[k]+s[k]==$3-$2 && st[1]==0)g++} END{print "QC", n, "plus="p/n, "valid="g/n}' $OUT/$ST.bed
