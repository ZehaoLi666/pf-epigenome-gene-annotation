#!/bin/bash
#SBATCH --job-name=tama_merge_ev
#SBATCH --account=lerochlab
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=2
#SBATCH --mem=60G
#SBATCH --time=24:00:00
#SBATCH -o ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/merge_ev_%j.out
#SBATCH -e ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/merge_ev_%j.err
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
set -uo pipefail
PROJ=${PF_HPCC_BIGDATA}/GSE186109_longread
C=$PROJ/tama/collapse_up
M=$PROJ/tama/merge; mkdir -p $M; cd $M

# ---- concatenate per-chromosome ONT collapses with globally renumbered gene ids
O=$PROJ/tama/collapse_ont_chr
: > $C/GSE_ont_mixed.bed; : > $C/GSE_ont_mixed_trans_read.bed; OFF=0
for i in $(seq 1 16); do
  T=c$i; [ -f $O/$T.EMPTY ] && continue
  [ -f $O/$T.DONE ] || { echo "ONT chunk $T not done"; exit 1; }
  awk -v off=$OFF 'BEGIN{OFS="\t"}{split($4,a,";"); g=substr(a[1],2)+off; split(a[2],b,"."); $4="G"g";G"g"."b[2]; print}' $O/$T.bed >> $C/GSE_ont_mixed.bed
  awk -v off=$OFF 'BEGIN{OFS="\t"}{split($4,a,";"); split(a[1],b,"."); g=substr(b[1],2)+off; $4="G"g"."b[2]";"a[2]; print}' $O/${T}_trans_read.bed >> $C/GSE_ont_mixed_trans_read.bed
  MX=$(awk '{split($4,a,";"); g=substr(a[1],2)+0; if(g>m)m=g}END{print m+0}' $O/$T.bed)
  OFF=$((OFF+MX))
done
echo "ONT concat tx=$(wc -l < $C/GSE_ont_mixed.bed) reads=$(wc -l < $C/GSE_ont_mixed_trans_read.bed) genes_offset=$OFF"
for s in GSE_ring GSE_troph GSE_schizont GSE_schiztroph_cyto GSE_ont_mixed LeeVV2021_relaxed Marcus_relaxed Mengquan_relaxed; do
  [ -s $C/$s.bed ] && [ -s $C/${s}_trans_read.bed ] || { echo "MISSING $s"; exit 1; }
done
printf '%s\t%s\t%s\t%s\n' \
 $C/GSE_ring.bed capped 1,1,2 GSE_ring \
 $C/GSE_troph.bed capped 1,1,2 GSE_troph \
 $C/GSE_schizont.bed capped 1,1,2 GSE_schizont \
 $C/GSE_schiztroph_cyto.bed capped 1,1,2 GSE_schiztroph_cyto \
 $C/GSE_ont_mixed.bed capped 2,2,3 GSE_ont_mixed \
 $C/LeeVV2021_relaxed.bed no_cap 3,2,1 LeeVV2021 \
 $C/Marcus_relaxed.bed no_cap 3,2,1 Marcus \
 $C/Mengquan_relaxed.bed no_cap 3,2,1 Mengquan > $M/merge_filelist.txt
cat $M/merge_filelist.txt
PY=${PF_HPCC_BIGDATA}/.conda/envs/tama_py27/bin/python
$PY $PROJ/tama/tama_upstream/tama_merge.py -f $M/merge_filelist.txt -p MASTER -a 20 -z 100 -m 10 -d merge_dup > $M/merge.log 2>&1
echo "merge rc=$? tx=$(wc -l < $M/MASTER.bed) genes=$(cut -f4 $M/MASTER.bed | cut -d';' -f1 | sort -u | wc -l)"
tail -2 $M/merge.log
awk '{n++; if($6=="+")p++; split($11,s,","); split($12,st,","); k=$10; if(st[k]+s[k]==$3-$2 && st[1]==0)g++} END{print "QC", n, "plus="p/n, "valid="g/n}' $M/MASTER.bed
source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda activate ${PF_HPCC_BIGDATA}/.conda/envs/cleavage
python $PROJ/scripts/evidence_tiers.py > $PROJ/logs/evidence.log 2>&1
echo "evidence rc=$?"; tail -40 $PROJ/logs/evidence.log
