#!/bin/bash
#SBATCH --job-name=riboseq_map --account=lerochlab --partition=epyc --cpus-per-task=8 --mem=24G --time=16:00:00 --array=1-8
#SBATCH --output=logs/map_%A_%a.out
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
# Caro et al. 2014 eLife ribosome footprints (GSE58402 / PRJNA252511; W2 strain; stages 1-4 = ring, early troph, late troph, schizont; 2 sequencing runs).
# Download (FASTQ kept as raw data) -> trim 3' linker CTGTAGGCACCATCAAT, keep 25-35 nt inserts -> bowtie2 end-to-end to 3D7 v48 (= v68 nuclear sequence)
# -> MAPQ>=10 -> counts of read 5' ends by (chrom, read strand, 5' position, read length): sites/<run>.tsv.gz
set -euo pipefail
W=${PF_HPCC_BIGDATA}/riboseq_caro2014; PKG=/opt/linux/rocky/8.x/x86_64/pkgs
export LD_LIBRARY_PATH=$PKG/htslib/1.19.1/lib:$PKG/libdeflate/1.10/lib:$PKG/isa-l/8b7c1b8/lib:${LD_LIBRARY_PATH:-}; ST=$PKG/samtools/1.19.2/bin/samtools
BT=${PF_HPCC_BIGDATA}/.conda/envs/cleavage/bin/bowtie2; PY=${PF_HPCC_BIGDATA}/.conda/envs/cleavage/bin/python3
IDX=${PF_HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/bowtie2_index/Pfalciparum_v48
read RUN STG URL < <(sed -n "${SLURM_ARRAY_TASK_ID}p" $W/scripts/runs.tsv)
mkdir -p $W/fastq $W/sites $W/stats; cd $W/fastq
FQD=${PF_HPCC_BIGDATA}/SRA_toolkit/sratoolkit.3.0.10-ubuntu64/bin/fasterq-dump
if [ ! -s $RUN.sra ]; then curl -sS --retry 5 -o $RUN.sra.part https://sra-pub-run-odp.s3.amazonaws.com/sra/$RUN/$RUN && mv $RUN.sra.part $RUN.sra; fi
$FQD --stdout --threads 4 -t $TMPDIR $RUN.sra 2> $W/stats/$RUN.fqd.txt \
 | $PY $W/scripts/trim.py > $TMPDIR/$RUN.trim.fq 2> $W/stats/$RUN.trim.txt
$BT -p 8 -x $IDX -U $TMPDIR/$RUN.trim.fq --end-to-end --no-unal 2> $W/stats/$RUN.bt2.txt | $ST view -q 10 -F 0x904 - \
 | awk 'BEGIN{OFS="\t"} {L=length($10); if (and($2,16)) {s="-"; cig=$6; n=0; m=0; while (match(cig,/[0-9]+[MDN=X]/)) {n+=substr(cig,RSTART,RLENGTH-1); cig=substr(cig,RSTART+RLENGTH)}; p=$4-1+n-1} else {s="+"; p=$4-1}; print $3,s,p,L}' \
 | sort -S 4G -T $TMPDIR | uniq -c | awk 'BEGIN{OFS="\t"}{print $2,$3,$4,$5,$1}' | gzip -c > $W/sites/$RUN.tsv.gz
rm -f $TMPDIR/$RUN.trim.fq
echo "DONE $RUN $STG $(zcat $W/sites/$RUN.tsv.gz | wc -l)"
