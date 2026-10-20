#!/bin/bash
#SBATCH --job-name=daft_map
#SBATCH --account=lerochlab
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --array=1-14
#SBATCH --output=logs/daft_%A_%a.out
#SBATCH --error=logs/daft_%A_%a.err
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
# Chappell et al. 2020 BMC Genomics, ENA ERP001570: 3D7 DAFT-seq (directional, amplification-free), 7 time points x 2 lanes, paired-end.
# Download -> HISAT2 (v48 index; v48 = v68 nuclear sequence) -> properly paired, MAPQ>=10, sorted BAM. FASTQs are kept (raw data).
set -euo pipefail
module load hisat2/2.2.1 samtools
W=${PF_HPCC_BIGDATA}/flank_coverage_2026
IDX=${PF_HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome
read RUN TP U1 U2 < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $W/scripts/daft_runs.tsv)
mkdir -p $W/fastq $W/bam $W/stats
cd $W/fastq
for u in $U1 $U2; do f=$(basename $u); [ -s $f ] || curl -sS --retry 5 -o $f.part $u && mv -f $f.part $f 2>/dev/null || true; done
[ -s ${RUN}_1.fastq.gz ] && [ -s ${RUN}_2.fastq.gz ] || { echo "ERROR download $RUN"; exit 1; }
gzip -t ${RUN}_1.fastq.gz && gzip -t ${RUN}_2.fastq.gz
hisat2 -p 8 -x $IDX -1 ${RUN}_1.fastq.gz -2 ${RUN}_2.fastq.gz --max-intronlen 5000 --no-unal --new-summary --summary-file $W/stats/${RUN}.hisat2.txt \
  | samtools view -@ 2 -b -q 10 -f 2 - | samtools sort -@ 4 -m 2G -T $TMPDIR/${RUN} -o $W/bam/${RUN}_t${TP}.bam -
samtools index $W/bam/${RUN}_t${TP}.bam
samtools flagstat $W/bam/${RUN}_t${TP}.bam > $W/stats/${RUN}.flagstat.txt
echo "DONE $RUN t$TP"
