#!/bin/bash
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
set -uo pipefail
cd ${PF_HPCC_BIGDATA}/GSE186109_longread
pigz -f -p 4 fastq/SRR16480183.fastq
if gzip -t fastq/SRR16480183.fastq.gz; then
  echo "GZIP_OK_183 $(date)" >> logs/download.log
  rm -rf sra/SRR16480183
else
  echo "GZIP_FAIL_183 $(date)" >> logs/download.log
fi
bash scripts/download2.sh
