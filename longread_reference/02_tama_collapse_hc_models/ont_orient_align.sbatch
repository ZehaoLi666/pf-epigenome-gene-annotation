#!/bin/bash
#SBATCH --job-name=ont_orient
#SBATCH --account=lerochlab
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=8
#SBATCH --mem=40G
#SBATCH --time=12:00:00
#SBATCH -o ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/ontorient_%j.out
#SBATCH -e ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/ontorient_%j.err
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
set -uo pipefail
PROJ=${PF_HPCC_BIGDATA}/GSE186109_longread
PKG=/opt/linux/rocky/8.x/x86_64/pkgs
MM2=$PKG/minimap2/2.30/bin/minimap2
SAMTOOLS=$PKG/samtools/1.22.1/bin/samtools
export LD_LIBRARY_PATH=$PKG/htslib/1.22.1/lib:$PKG/libdeflate/1.17/lib:${LD_LIBRARY_PATH:-}
export TMPDIR=$PROJ/tmp
source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda activate ${PF_HPCC_BIGDATA}/.conda/envs/cleavage
python $PROJ/scripts/ont_orient.py || exit 1
cat $PROJ/diag_ont_orientation.tsv
GENOME=${PF_HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta
OUT=$PROJ/bam/ont_mixed_oriented; mkdir -p $OUT
"$MM2" -ax splice -uf -k14 --MD -t 8 "$GENOME" $PROJ/by_sample/ont_mixed/oriented.fastq 2> $OUT/oriented.mm2.log | "$SAMTOOLS" sort -@ 4 -m 2G -T $OUT/tmp -o $OUT/oriented.sorted.bam -
"$SAMTOOLS" index $OUT/oriented.sorted.bam
NFQ=$(( $(wc -l < $PROJ/by_sample/ont_mixed/oriented.fastq) / 4 ))
echo "ORIENTED_FASTQ=$NFQ PRIMARY_MAPPED=$("$SAMTOOLS" view -c -F 0x904 $OUT/oriented.sorted.bam)" | tee $OUT/align_stats.txt
# TAMA capped collapse of oriented ONT (upstream TAMA, py2.7)
PY=${PF_HPCC_BIGDATA}/.conda/envs/tama_py27/bin/python
C=$PROJ/tama/collapse_up; cd $C
"$SAMTOOLS" view -h -F 0x904 $OUT/oriented.sorted.bam > $C/GSE_ont_mixed.pri.sam
$PY $PROJ/tama/tama_upstream/tama_collapse.py -s $C/GSE_ont_mixed.pri.sam -b SAM -f $GENOME -p GSE_ont_mixed -x capped -e common_ends -a 20 -z 20 -c 95 -i 75 > $C/GSE_ont_mixed.collapse.log 2>&1
echo "TAMA rc=$? tx=$(wc -l < $C/GSE_ont_mixed.bed)"
rm -f $C/GSE_ont_mixed.pri.sam
