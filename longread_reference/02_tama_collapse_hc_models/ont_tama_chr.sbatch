#!/bin/bash
#SBATCH --job-name=ont_tama_chr
#SBATCH --account=lerochlab
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH --array=1-16
#SBATCH -o ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/ontchr_%A_%a.out
#SBATCH -e ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/ontchr_%A_%a.err
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
set -uo pipefail
PROJ=${PF_HPCC_BIGDATA}/GSE186109_longread
PKG=/opt/linux/rocky/8.x/x86_64/pkgs
SAMTOOLS=$PKG/samtools/1.22.1/bin/samtools
export LD_LIBRARY_PATH=$PKG/htslib/1.22.1/lib:$PKG/libdeflate/1.17/lib:${LD_LIBRARY_PATH:-}
export TMPDIR=$PROJ/tmp
PY=${PF_HPCC_BIGDATA}/.conda/envs/tama_py27/bin/python
GEN=${PF_HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta
BAM=$PROJ/bam/ont_mixed_oriented/oriented.sorted.bam
CHR=$("$SAMTOOLS" idxstats $BAM | awk -v i=$SLURM_ARRAY_TASK_ID 'NR==i{print $1}')
[ -n "$CHR" ] || exit 0
O=$PROJ/tama/collapse_ont_chr; mkdir -p $O; cd $O
T=c${SLURM_ARRAY_TASK_ID}
"$SAMTOOLS" view -h -F 0x904 $BAM $CHR > $O/$T.sam
N=$(grep -vc '^@' $O/$T.sam)
echo "$CHR $T reads=$N"
if [ "$N" -eq 0 ]; then rm -f $O/$T.sam; touch $O/$T.EMPTY; exit 0; fi
$PY $PROJ/tama/tama_upstream/tama_collapse.py -s $O/$T.sam -b SAM -f $GEN -p $T -x capped -e common_ends -a 20 -z 20 -c 95 -i 75 > $O/$T.log 2>&1
RC=$?; echo "rc=$RC tx=$(wc -l < $O/$T.bed)"
[ $RC -eq 0 ] && rm -f $O/$T.sam && touch $O/$T.DONE
