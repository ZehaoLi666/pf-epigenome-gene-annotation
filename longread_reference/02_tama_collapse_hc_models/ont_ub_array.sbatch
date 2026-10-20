#!/bin/bash
#SBATCH --job-name=gse_ont_ub
#SBATCH --account=lerochlab
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=30G
#SBATCH --time=06:00:00
#SBATCH --array=1-2
#SBATCH -o ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/ontub_%A_%a.out
#SBATCH -e ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/ontub_%A_%a.err
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
set -uo pipefail
PROJ=${PF_HPCC_BIGDATA}/GSE186109_longread
PKG=/opt/linux/rocky/8.x/x86_64/pkgs
MM2=$PKG/minimap2/2.30/bin/minimap2
SAMTOOLS=$PKG/samtools/1.22.1/bin/samtools
export LD_LIBRARY_PATH=$PKG/htslib/1.22.1/lib:$PKG/libdeflate/1.17/lib:${LD_LIBRARY_PATH:-}
export TMPDIR=$PROJ/tmp
"$SAMTOOLS" --version >/dev/null || { echo "samtools libs"; exit 1; }
GENOME=${PF_HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta
case $SLURM_ARRAY_TASK_ID in
 1) SRR=SRR18188665; TAG=run1;;
 2) SRR=SRR18188666; TAG=run2;;
esac
OUT=$PROJ/bam/ont_mixed_ub; mkdir -p $OUT
FQ=$PROJ/fastq/$SRR.fastq
# -ub: infer transcript strand from splice motifs (ts tag), soft-clips kept (-Y) for poly(A)/poly(T) orientation
"$MM2" -ax splice -ub -k14 -Y --MD -t 8 "$GENOME" "$FQ" 2> $OUT/$TAG.mm2.log | "$SAMTOOLS" sort -@ 4 -m 2G -T $OUT/$TAG.tmp -o $OUT/$TAG.sorted.bam -
"$SAMTOOLS" index $OUT/$TAG.sorted.bam
NFQ=$(( $(wc -l < "$FQ") / 4 ))
PRI=$("$SAMTOOLS" view -c -F 0x904 $OUT/$TAG.sorted.bam)
UNM=$("$SAMTOOLS" view -c -f 4 $OUT/$TAG.sorted.bam)
echo "READS_IN_FASTQ=$NFQ PRIMARY_MAPPED=$PRI UNMAPPED=$UNM SUM=$((PRI+UNM))" > $OUT/$TAG.align_stats.txt
"$SAMTOOLS" quickcheck $OUT/$TAG.sorted.bam && echo "DONE $TAG"
