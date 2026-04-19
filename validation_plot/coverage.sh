
#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BW_DIR="${BW_DIR:-${SCRIPT_DIR}/../window1000_bin100/coverage_pattern/new_coverage/bw}"

# protein-coding gene TSS 
computeMatrix reference-point --referencePoint TSS -R predicted_UTR5_with_v68_replacements.bed intergenic.bed\
    -S "${BW_DIR}/Cage.bw" \
    -b 500 -a 500 --binSize 10 --skipZeros --missingDataAsZero -p max -o predicted_Cage_TSS.gz
plotProfile -m predicted_Cage_TSS.gz --plotType se -T "Signal around TSS (±500bp)" -out predicted_TSS_Cage_profile.png --outFileNameData predicted_TSS_Cage.tsv\
       --regionsLabel "Refined_genes" "intergenic" --samplesLabel Cage-seq --dpi 300 



# protein-coding gene TES 
computeMatrix reference-point --referencePoint TES -R predicted_UTR5_with_v68_replacements.bed intergenic.bed\
    -S "${BW_DIR}/PolyA.bw" \
    -b 500 -a 500 --binSize 10 --skipZeros --missingDataAsZero -p max -o predicted_PolyA_TES.gz
plotProfile -m predicted_PolyA_TES.gz --plotType se -T "Signal around TES (±500bp)" -out predicted_TES_PolyA_profile.png --outFileNameData predicted_TES_PolyA.tsv\
       --regionsLabel "Refined_genes" "intergenic" --samplesLabel PolyA-seq --dpi 300 


# lncRNA TSS 
computeMatrix reference-point --referencePoint TSS -R novel_units.bed intergenic.bed\
    -S "${BW_DIR}/GRO.bw" \
    -b 250 -a 250 --binSize 10 --skipZeros --missingDataAsZero -p max -o lncRNA_Cage_TSS.gz
plotProfile -m lncRNA_Cage_TSS.gz --plotType se -T "Signal around TSS (±500bp)" -out lncRNA_TSS_Cage_profile.png --outFileNameData lncRNA_TSS_Cage.tsv\
       --regionsLabel "lncRNA" "intergenic" --samplesLabel Cage-seq --dpi 300 

# lncRNA TES 
computeMatrix reference-point --referencePoint TES -R novel_units.bed intergenic.bed\
    -S "${BW_DIR}/PolyA.bw" \
    -b 250 -a 250 --binSize 10 --skipZeros --missingDataAsZero -p max -o lncRNA_PolyA_TES.gz
plotProfile -m lncRNA_PolyA_TES.gz --plotType se -T "Signal around TES (±500bp)" -out lncRNA_TES_PolyA_profile.png --outFileNameData lncRNA_TES_PolyA.tsv\
       --regionsLabel "lncRNA" "intergenic" --samplesLabel PolyA-seq --dpi 300 






