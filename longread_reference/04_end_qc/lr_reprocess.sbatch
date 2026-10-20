#!/bin/bash
#SBATCH --job-name=lr_reprocess
#SBATCH --account=lerochlab
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --array=1-7
#SBATCH --output=logs/lr_%A_%a.out
#SBATCH --error=logs/lr_%A_%a.err
source "${PF_REPO:?set PF_REPO to the repository path}/config/hpcc.env"   # site paths
# Re-processing of three published long-read datasets for UTR validation (new folder; originals untouched).
#  1  LeeVV 2021 (SRR11094274)  ONT direct RNA  -> minimap2 -ax splice -uf -k14 -G 5000
#  2  Hoshizaki 2022 (ERR9746080) ONT direct RNA -> same
#  3-7 Yang 2021 (CRR216516-20) PacBio Iso-Seq: isoseq refine WITHOUT --require-polya (all FLNC), then minimap2 splice:hq
set -euo pipefail
W=${PF_HPCC_BIGDATA}/longread_reprocess_2026
B=${PF_HPCC_LAB}/SRA_toolkit/Long_read_RNA_datasets
GENOME=${PF_HPCC_LAB}/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta
PKG=/opt/linux/rocky/8.x/x86_64/pkgs
MM2=$PKG/minimap2/2.24/bin/minimap2
SAMTOOLS=$PKG/samtools/1.19.2/bin/samtools
E=${PF_HPCC_BIGDATA}/.conda/envs/cleavage/bin
export LD_LIBRARY_PATH=$PKG/htslib/1.19.1/lib:$PKG/libdeflate/1.10/lib:$PKG/isa-l/8b7c1b8/lib:${LD_LIBRARY_PATH:-}
export TMPDIR=$W/tmp/$SLURM_ARRAY_TASK_ID; mkdir -p $TMPDIR $W/ont $W/pacbio $W/stats
"$MM2" --version >/dev/null 2>&1 || { echo "ERROR minimap2"; exit 1; }
"$SAMTOOLS" --version >/dev/null 2>&1 || { echo "ERROR samtools libs"; exit 1; }
[ -s "$GENOME" ] || { echo "ERROR genome"; exit 1; }
T=$SLURM_ARRAY_TASK_ID
map_sort () {  # $1 preset args, $2 fastq, $3 out prefix
  "$MM2" -t 8 $1 --secondary=no "$GENOME" "$2" 2> "$3.mm2.log" | "$SAMTOOLS" sort -@ 4 -m 2G -T "$TMPDIR/s" -o "$3.sorted.bam" -
  "$SAMTOOLS" index "$3.sorted.bam"
  "$SAMTOOLS" flagstat "$3.sorted.bam" > "$W/stats/$(basename $3).flagstat.txt"
}
case $T in
  1|2)
    if [ $T = 1 ]; then NAME=LeeVV_2021_dRNA; FQ=$B/LeeVV_2021_nanopore/SRR11094274.fastq; else NAME=Hoshizaki_2022_dRNA; FQ=$B/Marcus_Nanopore/ERR9746080.fastq; fi
    [ -s "$FQ" ] || { echo "ERROR missing $FQ"; exit 1; }
    echo "input reads $NAME: $(( $(wc -l < $FQ) / 4 ))"
    map_sort "-ax splice -uf -k14 -G 5000" "$FQ" "$W/ont/$NAME"
    ;;
  3|4|5|6|7)
    ID=CRR2165$(printf '%02d' $((T + 13)))
    S=$B/Mengquan_long_read/Second_analysis; FL=$S/$ID.demux.primerIIA_5p--PrimerIIA_3p.bam; P=$B/Mengquan_long_read/primers.fasta
    [ -s "$FL" ] && [ -s "$P" ] || { echo "ERROR missing $FL or $P"; exit 1; }
    "$E/isoseq" refine "$FL" "$P" "$W/pacbio/$ID.flnc_all.bam" -j 8 --log-level INFO
    "$SAMTOOLS" fastq "$W/pacbio/$ID.flnc_all.bam" > "$TMPDIR/$ID.flnc_all.fastq"
    echo "flnc reads $ID: $(( $(wc -l < $TMPDIR/$ID.flnc_all.fastq) / 4 ))"
    map_sort "-ax splice:hq -uf -G 5000" "$TMPDIR/$ID.flnc_all.fastq" "$W/pacbio/${ID}_flnc_all"
    rm -f "$TMPDIR/$ID.flnc_all.fastq"
    ;;
esac
rm -rf "$TMPDIR"
echo "DONE task $T"
