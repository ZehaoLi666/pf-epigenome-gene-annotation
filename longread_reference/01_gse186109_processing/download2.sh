#!/bin/bash
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
set -uo pipefail
PROJ=${PF_HPCC_BIGDATA}/GSE186109_longread
mkdir -p $PROJ/sra $PROJ/fastq $PROJ/logs
export TMPDIR=$PROJ/tmp
mkdir -p $TMPDIR
PREFETCH=/opt/linux/rocky/8.x/x86_64/pkgs/sratoolkit/3.2.0/bin/prefetch
FQDUMP=/opt/linux/rocky/8.x/x86_64/pkgs/sratoolkit/3.2.0/bin/fasterq-dump
RUNS="SRR16480181 SRR16480182 SRR16480184 SRR16480185 SRR16480186 SRR16480187 SRR16480188 SRR16480189 SRR16480190 SRR16480191 SRR16480192 SRR16480193 SRR16480194 SRR16480195 SRR16480196 SRR16480197 SRR16480198 SRR16480199 SRR16480200 SRR16480201 SRR16480202 SRR16480203 SRR16480204 SRR16480205 SRR16480206 SRR16480207 SRR16480208 SRR18188665 SRR18188666"
echo "START $(date)" >> $PROJ/logs/download.log
for r in $RUNS; do
  echo "=== $r $(date) ===" >> $PROJ/logs/download.log
  if [ -f "$PROJ/fastq/$r.fastq.gz" ]; then
    echo "SKIP $r already done" >> $PROJ/logs/download.log
    continue
  fi
  $PREFETCH -O $PROJ/sra $r >> $PROJ/logs/download.log 2>&1
  $FQDUMP -O $PROJ/fastq -t $TMPDIR -e 4 $PROJ/sra/$r >> $PROJ/logs/download.log 2>&1
  if [ -s "$PROJ/fastq/$r.fastq" ]; then
    pigz -f -p 4 $PROJ/fastq/$r.fastq >> $PROJ/logs/download.log 2>&1
    if gzip -t "$PROJ/fastq/$r.fastq.gz" 2>>$PROJ/logs/download.log; then
      echo "GZIP_OK $r" >> $PROJ/logs/download.log
      rm -rf $PROJ/sra/$r
    else
      echo "GZIP_FAIL $r -- keeping sra for retry" >> $PROJ/logs/download.log
    fi
  else
    echo "DUMP_FAIL $r -- no fastq produced" >> $PROJ/logs/download.log
  fi
done
echo "ALL_DONE $(date)" >> $PROJ/logs/download.log
touch $PROJ/logs/DOWNLOAD_DONE
