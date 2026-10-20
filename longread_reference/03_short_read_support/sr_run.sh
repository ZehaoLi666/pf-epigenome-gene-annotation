#!/bin/bash
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
# runs on interactive node via nohup setsid srun --overlap; writes DONE flag
set -uo pipefail
PROJ=${PF_HPCC_BIGDATA}/GSE186109_longread
B=$PROJ/shortread; mkdir -p $B/junctions $B/counts $B/logs
source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda activate ${PF_HPCC_BIGDATA}/.conda/envs/cleavage
P=${PF_HPCC_LAB}/SRA_toolkit/RNAseq_7stages/processed_data
# 1) junctions, 2 at a time
for pair in "Ring Trophozoite" "Schizont Gametocyte" "ET LT" "Gii GV" "Ookinete"; do
  for s in $pair; do [ -s $B/junctions/$s.SJ.out.tab ] || python $PROJ/scripts/sr_junctions.py $s > $B/logs/junc_$s.log 2>&1 & done; wait
done
# 2) combined annotation: HC loci + v48 genes lacking an HC locus
V48=${PF_HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7.gtf
HC=$PROJ/reference_v1/hc_reference.gtf
cut -f9 $HC | grep -o 'gene_id "[^"]*"' | sort -u > $B/counts/hc_loci.txt
awk -F'\t' 'NR==FNR{split($0,a,"\""); h[a[2]]=1; next} $3=="exon"{match($9,/gene_id "[^"]*"/); g=substr($9,RSTART+9,RLENGTH-10); if(!(g in h)) print}' $B/counts/hc_loci.txt $V48 > $B/counts/v48_uncovered_exons.gtf
awk -F'\t' '$3=="exon"' $HC | cat - $B/counts/v48_uncovered_exons.gtf > $B/counts/combined_loci.gtf
echo "combined exons $(wc -l < $B/counts/combined_loci.gtf) uncovered_v48_genes $(cut -f9 $B/counts/v48_uncovered_exons.gtf | grep -o 'gene_id "[^"]*"' | sort -u | wc -l)"
module load subread
featureCounts -v 2>&1 | head -2
featureCounts -a $B/counts/combined_loci.gtf -t exon -g gene_id -s 0 -T 2 -p --countReadPairs -B -C -o $B/counts/featurecounts_pe.txt $P/Ring/Ring_sorted.bam $P/Trophozoite/Trophozoite_sorted.bam $P/Schizont/Schizont_sorted.bam $P/Gametocyte/Gametocyte_sorted.bam > $B/logs/fc_pe.log 2>&1
echo "fc_pe rc=$?"
featureCounts -a $B/counts/combined_loci.gtf -t exon -g gene_id -s 0 -T 2 -o $B/counts/featurecounts_se.txt $P/ET/ET_sorted.bam $P/LT/LT_sorted.bam $P/Gii/Gii_sorted.bam $P/GV/GV_sorted.bam $P/Ookinete/Ookinete_sorted.bam > $B/logs/fc_se.log 2>&1
echo "fc_se rc=$?"
python $PROJ/scripts/sr_pool_tpm.py
echo "pool rc=$?"
touch $B/SR_DONE
