#!/bin/bash
#SBATCH --account=lerochlab --partition=epyc --cpus-per-task=8 --mem=32G --time=04:00:00
#SBATCH --job-name=abl_repro --output=logs/repro_%j.out
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
source ${PF_HPCC_BIGDATA}/ablation_2026/work/hpcc/env.sh
for n in V3_top8 V3_top8_s1 V3_top8_drop; do $PY predict_logit.py $n 8 2>&1 | grep -v -i warn | tail -1; done
$PY - <<'PYEOF' 2>&1 | grep -v -i warn
import pandas as pd, peak_eval as PE
R = pd.DataFrame(PE.report_fixed('V3 top-8 (re-predicted on HPCC)', 'V3_top8+V3_top8_s1+V3_top8_drop', 'logit_nms'))
R.to_csv('repro_check.tsv', sep='\t', index=False)
print(R.pivot_table(index='end_type', columns='scheme', values='F1@200').round(3).to_string())
PYEOF
echo REPRO_DONE
