#!/bin/bash
#SBATCH --account=lerochlab --partition=epyc --cpus-per-task=4 --mem=48G --time=08:00:00
#SBATCH --job-name=abl_eval --output=logs/eval_%j.out
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
source ${PF_HPCC_BIGDATA}/ablation_2026/work/hpcc/env.sh
for n in V3_nochrom V3_nochrom_s1 V3_nochrom_drop; do
  [ -s exp_ml/logit_$n.npz ] || $PY predict_logit.py $n 4 2>&1 | grep -v -i warn | tail -1
  [ -s exp_ml/pred_$n.npz ] || $PY hpcc/pred_prob.py $n 4 2>&1 | grep -v -i warn | tail -1; done
$PY hpcc/abl_eval.py 2>&1 | grep -v -i warn
cd ${PF_HPCC_BIGDATA}/ablation_2026/work && tar -czf abl_results.tgz abl_point_metrics.tsv abl_gene_calls.tsv abl_extra.tsv exp_ml/*_hist.json && md5sum abl_results.tgz
