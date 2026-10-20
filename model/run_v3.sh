#!/bin/bash
cd "$(dirname "$0")"
export PYTHONPATH=.:boundary_model POOL=1 ML_LABELS=data/Lml3soft40_pool.npz
NEW=chrom,onehot,glm,gro,poly,strand,daft
r(){ python exp_ml.py --name $1 ${@:2} --threads 5 > exp_ml/log_$1.txt 2>&1 && python predict_logit.py $1 5 > /dev/null 2>&1; echo "done $1"; }
three(){ n=$1; shift; r ${n} "$@"; r ${n}_s1 "$@" --seed 1; r ${n}_drop "$@" --dropout 0.3 --wd 0.01; }
( three V3_all --spec $NEW; three V3_noAP2 --spec $NEW --tracks H2AZ,H2BZ,H3,H3K27ac,H3K4me1,H3K4me3,H3K9ac,H3K9me3,ATAC,MNase,HP1,PolII_S2P; three V3_old_inputs --spec chrom,onehot,glm,rna,strand ) &
( three V3_top8 --spec $NEW --tracks H2BZ,MNase,H3K4me3,PolII_S2P,H3K9ac,H2AZ,ATAC,H3K27ac; three V3_nochrom --spec onehot,glm,gro,poly,strand,daft ) &
wait; echo ALLDONE
