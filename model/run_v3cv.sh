#!/bin/bash
cd "$(dirname "$0")"
export PYTHONPATH=.:boundary_model POOL=1 ML_LABELS=data/Lml3soft40_pool.npz
NEW=chrom,onehot,glm,gro,poly,strand,daft; T8=H2BZ,MNase,H3K4me3,PolII_S2P,H3K9ac,H2AZ,ATAC,H3K27ac
r(){ python exp_ml.py --name $1 ${@:2} --threads 4 > exp_ml/log_$1.txt 2>&1 && PRED_ALL=1 python predict_logit.py $1 4 > /dev/null 2>&1; echo "done $1"; }
for f in "B --test 6,7,8,9,10 --val 11" "C --test 11,12,13,14 --val 10"; do set -- $f; F=$1; shift
  r V3cv${F}_top8 --spec $NEW --tracks $T8 "$@"; r V3cv${F}_top8_s1 --spec $NEW --tracks $T8 --seed 1 "$@"
  r V3cv${F}_nochrom --spec onehot,glm,gro,poly,strand,daft "$@"; r V3cv${F}_nochrom_s1 --spec onehot,glm,gro,poly,strand,daft --seed 1 "$@"; done
echo ALLDONE
