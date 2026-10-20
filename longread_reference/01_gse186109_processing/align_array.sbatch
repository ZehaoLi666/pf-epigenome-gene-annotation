#!/bin/bash
#SBATCH --job-name=gse186109_align
#SBATCH --account=lerochlab
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --array=1-15
#SBATCH -o ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/align_%A_%a.out
#SBATCH -e ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/align_%A_%a.err
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths

set -uo pipefail
PROJ=${PF_HPCC_BIGDATA}/GSE186109_longread
PKG=/opt/linux/rocky/8.x/x86_64/pkgs
MM2=$PKG/minimap2/2.30/bin/minimap2
SAMTOOLS=$PKG/samtools/1.22.1/bin/samtools
export LD_LIBRARY_PATH=$PKG/htslib/1.22.1/lib:$PKG/libdeflate/1.17/lib:$LD_LIBRARY_PATH

GENOME=${PF_HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta
TASKFILE=$PROJ/scripts/align_tasks.tsv
OUTDIR=$PROJ/bam
mkdir -p $OUTDIR

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $TASKFILE)
SAMPLE=$(echo "$LINE" | cut -f1)
PLATFORM=$(echo "$LINE" | cut -f2)
TIER=$(echo "$LINE" | cut -f3)
FASTQ=$(echo "$LINE" | cut -f4)

echo "=== $SAMPLE / $TIER / $PLATFORM : $FASTQ ==="
"$SAMTOOLS" --version | head -1
"$MM2" --version

mkdir -p $OUTDIR/$SAMPLE
OUT=$OUTDIR/$SAMPLE/${TIER}.sorted.bam

if [ "$PLATFORM" == "ont" ]; then
  PRESET="-ax splice -uf -k14"
else
  PRESET="-ax splice -uf -C5"
fi

"$MM2" $PRESET -t ${SLURM_CPUS_PER_TASK} --MD "$GENOME" "$FASTQ" 2> $OUTDIR/$SAMPLE/${TIER}.mm2.log \
  | "$SAMTOOLS" sort -@ 2 -m 2G -T $OUTDIR/$SAMPLE/${TIER}.tmpsort -o "$OUT" -
"$SAMTOOLS" index "$OUT"

TOTAL=$("$SAMTOOLS" view -c "$FASTQ" 2>/dev/null || echo NA)
MAPPED=$("$SAMTOOLS" view -c -F 4 "$OUT")
UNMAPPED=$("$SAMTOOLS" view -c -f 4 "$OUT")
NREADS_FQ=$(( $(wc -l < "$FASTQ") / 4 ))
echo "READS_IN_FASTQ=$NREADS_FQ MAPPED=$MAPPED UNMAPPED=$UNMAPPED SUM=$((MAPPED+UNMAPPED))" > $OUTDIR/$SAMPLE/${TIER}.align_stats.txt
echo "DONE $SAMPLE $TIER"
