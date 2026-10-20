#!/bin/bash
#SBATCH --account=lerochlab --partition=epyc --cpus-per-task=8 --mem=32G --time=12:00:00
#SBATCH --job-name=abl_train --array=1-13 --output=logs/train_%A_%a.out
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
source ${PF_HPCC_BIGDATA}/ablation_2026/work/hpcc/env.sh
T8=H2BZ,MNase,H3K4me3,PolII_S2P,H3K9ac,H2AZ,ATAC,H3K27ac
IFS=$'\t' read -r NAME SPEC TRK < <(sed -n "${SLURM_ARRAY_TASK_ID}p" hpcc/configs.tsv)
TA=""; [ "$TRK" = "T8" ] && TA="--tracks $T8"
echo "task $SLURM_ARRAY_TASK_ID $NAME spec=$SPEC tracks=$TRK host=$(hostname)"
r(){ n=$1; shift
  if [ ! -s exp_ml/$n.pt ]; then $PY exp_ml.py --name $n --spec $SPEC $TA --threads 8 "$@" > logs/fit_$n.txt 2>&1 || { echo "FAIL fit $n"; tail -5 logs/fit_$n.txt; return 1; }; fi
  [ -s exp_ml/logit_$n.npz ] || $PY predict_logit.py $n 8 2>&1 | grep -v -i warn | tail -1
  echo "done $n $(date +%H:%M)"; }
r ${NAME}; r ${NAME}_s1 --seed 1; r ${NAME}_drop --dropout 0.3 --wd 0.01
echo TASK_DONE $NAME
