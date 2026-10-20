#!/bin/bash
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
export PATH=/opt/linux/rocky/8.x/x86_64/pkgs/slurm/24.11.1/bin:$PATH
JOBID=29010567
PKG=/opt/linux/rocky/8.x/x86_64/pkgs
SEQKIT=$PKG/seqkit/2.4.0/bin/seqkit
OUT=${PF_HPCC_BIGDATA}/GSE186109_longread/logs/qc.log
cd ${PF_HPCC_BIGDATA}/GSE186109_longread/fastq
echo "---HEADERS---" > $OUT
for f in SRR16480181 SRR16480182 SRR16480183 SRR16480194 SRR16480205 SRR18188665 SRR18188666; do
  echo "$f:" >> $OUT
  head -1 $f.fastq >> $OUT
done
echo "---STATS---" >> $OUT
srun --jobid=$JOBID --overlap $SEQKIT stats -a *.fastq >> $OUT 2>> $OUT
echo "QC_DONE $(date)" >> $OUT
touch ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/QC_DONE
