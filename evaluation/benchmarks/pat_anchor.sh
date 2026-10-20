#!/bin/bash
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
# PAT-seq cleavage sites from poly(T)-anchored reads only. Read structure: [adapter 0-50 nt][T>=15][insert, antisense to RNA].
# Insert 5' end (first aligned base) = poly(A) junction = cleavage site. Output: anchored/<SRR>.sites.tsv (chrom, RNA strand, pos 0-based, MAPQ).
set -eo pipefail
source /opt/linux/rhel/8.x/x86_64/pkgs/miniconda3/py39_4.12.0/etc/profile.d/conda.sh
conda activate ${PF_HPCC_BIGDATA}/.conda/envs/cleavage
O=${PF_HPCC_BIGDATA}/benchmark_2026/patseq_anchor; mkdir -p $O; cd $O
T=${PF_HPCC_BIGDATA}/lncRNA_longread_discovery/teseq/fastq
IDX=${PF_HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/bowtie2_index/Pfalciparum_v48
for S in SRR39981131 SRR39981130 SRR39981126 SRR39981125 SRR39981124 SRR39981123 SRR39981122 SRR39981121 SRR39981120; do
  [ -s $S.sites.tsv ] && continue
  python - $T/$S.fastq $S.anchor.fq <<'PY'
import sys, re
fi, fo = open(sys.argv[1]), open(sys.argv[2], 'w'); pat = re.compile('T{15,}'); n = k = 0
while True:
    h = fi.readline()
    if not h: break
    s = fi.readline().rstrip(); fi.readline(); q = fi.readline().rstrip(); n += 1
    m = pat.search(s, 0, 70)
    if not m or m.start() > 55: continue
    j = m.end()
    if len(s) - j < 25: continue
    fo.write(f"@{h[1:].split()[0]}\n{s[j:]}\n+\n{q[j:]}\n"); k += 1
print(sys.argv[1].split('/')[-1], 'reads', n, 'anchored', k, file=sys.stderr)
PY
  bowtie2 -x $IDX -U $S.anchor.fq -p 2 --no-unal --end-to-end 2> $S.bt2.log | python -c "
import sys, pysam
f = pysam.AlignmentFile('-', 'r'); o = open('$S.sites.tsv', 'w')
for r in f:
    if r.is_unmapped or r.is_secondary or r.is_supplementary: continue
    # read is antisense to RNA: forward read -> RNA on '-', site = alignment start; reverse read -> RNA on '+', site = alignment end - 1
    if r.is_reverse: o.write(f'{r.reference_name}\t+\t{r.reference_end - 1}\t{r.mapping_quality}\n')
    else: o.write(f'{r.reference_name}\t-\t{r.reference_start}\t{r.mapping_quality}\n')
"
  echo "$S $(grep overall $S.bt2.log) sites=$(wc -l < $S.sites.tsv)"
  rm -f $S.anchor.fq
done
cat *.sites.tsv | wc -l; echo ALLDONE
