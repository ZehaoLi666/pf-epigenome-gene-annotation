#!/bin/bash
#SBATCH --job-name=ont_c14_cap
#SBATCH --account=lerochlab
#SBATCH --partition=epyc
#SBATCH --cpus-per-task=1
#SBATCH --mem=40G
#SBATCH --time=24:00:00
#SBATCH -o ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/ontc14_%j.out
#SBATCH -e ${PF_HPCC_BIGDATA}/GSE186109_longread/logs/ontc14_%j.err
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
set -uo pipefail
PROJ=${PF_HPCC_BIGDATA}/GSE186109_longread
PKG=/opt/linux/rocky/8.x/x86_64/pkgs
SAMTOOLS=$PKG/samtools/1.22.1/bin/samtools
export LD_LIBRARY_PATH=$PKG/htslib/1.22.1/lib:$PKG/libdeflate/1.17/lib:${LD_LIBRARY_PATH:-}
PY=${PF_HPCC_BIGDATA}/.conda/envs/tama_py27/bin/python
GEN=${PF_HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta
O=$PROJ/tama/collapse_ont_chr; cd $O; rm -f c14.DONE
# cap read density: at most 1500 reads per (strand, 5-kb start bin), deterministic hash-based selection
"$SAMTOOLS" view -h -F 0x904 $PROJ/bam/ont_mixed_oriented/oriented.sorted.bam Pf3D7_14_v3 | awk 'BEGIN{OFS="\t"} /^@/{print; next} {k=and($2,16)"_"int($4/5000); n[k]++; if(n[k]<=1500) print; else d++} END{print "dropped="d+0 > "/dev/stderr"}' > c14.sam
echo "kept=$(grep -vc '^@' c14.sam)"
$PY $PROJ/tama/tama_upstream/tama_collapse.py -s c14.sam -b SAM -f $GEN -p c14 -x capped -e common_ends -a 20 -z 20 -c 95 -i 75 > c14.log 2>&1
RC=$?; echo "rc=$RC tx=$(wc -l < c14.bed)"
[ $RC -eq 0 ] && rm -f c14.sam && touch c14.DONE
