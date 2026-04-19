# regular RNA-seq TSS coverage 
computeMatrix reference-point --referencePoint TSS -R predicted_genes_TSS.dedup.bed v68_TSS.bed -S /rhome/zli529/lab/SRA_toolkit/RNAseq_7stages/processed_data/LT/RNA.bw -b 500 -a 500 --binSize 10 --skipZeros --missingDataAsZero -p max -o predicted_RNAmatrix_TSS.gz
plotProfile -m predicted_RNAmatrix_TSS.gz --plotType se -T "Signal around TSS (±500bp)" -out predicted_TSS_RNA_profile.png --outFileNameData predicted_TSS_RNAprofile_values.tsv

# GRO-seq TSS coverage 
bamCoverage -b /rhome/zli529/lab/SRA_toolkit/GRO_seq/aligned/SRR4019517.sorted.q10.pp.bam \
  -o GRO.bw \
  --normalizeUsing CPM \
  --binSize 10 \
  --minMappingQuality 10 \
  --ignoreDuplicates \
  -p max
computeMatrix reference-point --referencePoint TSS -R predicted_genes_TSS.dedup.bed v68_TSS.bed -S GRO.bw -b 500 -a 500 --binSize 10 --skipZeros --missingDataAsZero -p max -o predicted_GROmatrix_TSS.gz
plotProfile -m predicted_GROmatrix_TSS.gz --plotType se -T "Signal around TSS (±500bp)" -out predicted_TSS_GRO_profile.png --outFileNameData  predicted_TSS_GROprofile_values.tsv

# cage-seq TSS coverage 
bamCoverage -b /rhome/zli529/lab/SRA_toolkit/cage_seq_data/paired-end/aligned_hisat2_pe/SRR2031965.sorted.bam \
  -o Cage.bw \
  --normalizeUsing CPM \
  --binSize 10 \
  --minMappingQuality 10 \
  --ignoreDuplicates \
  -p max

computeMatrix reference-point --referencePoint TSS -R predicted_genes_TSS.dedup.bed v68_TSS.bed -S Cage.bw -b 500 -a 500 --binSize 10 --skipZeros --missingDataAsZero -p max -o predicted_Cagematrix_TSS.gz
plotProfile -m predicted_Cagematrix_TSS.gz --plotType se -T "Signal around TSS (±500bp)" -out predicted_TSS_Cage_profile.png --outFileNameData  predicted_TSS_Cageprofile_values.tsv




# regular RNA-seq TES coverage 
computeMatrix reference-point --referencePoint TES -R predicted_genes_TES.bed v68_TES.bed -S /rhome/zli529/lab/SRA_toolkit/RNAseq_7stages/processed_data/LT/RNA.bw -b 500 -a 500 --binSize 10 --skipZeros --missingDataAsZero -p max -o predicted_matrix_TES.gz
plotProfile -m predicted_matrix_TES.gz --plotType se -T "Signal around TES (±500bp)" -out predicted_TES_RNA_profile.png --outFileNameData predicted_TES_RNAprofile_values.tsv

# regular GRO-seq TES coverage 
computeMatrix reference-point --referencePoint TES -R predicted_genes_TES.bed v68_TES.bed -S GRO.bw -b 500 -a 500 --binSize 10 --skipZeros --missingDataAsZero -p max -o predicted_GROmatrix_TES.gz
plotProfile -m predicted_GROmatrix_TES.gz --plotType se -T "Signal around TES (±500bp)" -out predicted_TES_GRO_profile.png --outFileNameData predicted_TES_GROprofile_values.tsv


# regular PolyA-seq TES coverage 
bamCoverage -b /rhome/zli529/lab/SRA_toolkit/Poly-A-seq/Bunnik_2013/aligned/SRR836071.sorted.q10.pp.bam \
  -o PolyA.bw \
  --normalizeUsing CPM \
  --binSize 10 \
  --minMappingQuality 10 \
  --ignoreDuplicates \
  -p max

computeMatrix reference-point --referencePoint TES -R predicted_genes_TES.bed v68_TES.bed -S PolyA.bw  -b 500 -a 500 --binSize 10 --skipZeros --missingDataAsZero -p max -o predicted_PolyAmatrix_TES.gz
plotProfile -m predicted_PolyAmatrix_TES.gz --plotType se -T "Signal around TES (±500bp)" -out predicted_TES_PolyA_profile.png --outFileNameData predicted_TES_PolyAprofile_values.tsv



bamCoverage -b /rhome/zli529/lab/SRA_toolkit/Long_read_RNA_datasets/LeeVV_2021_long_read/merged.no_long_introns.sorted.bam -o Long_reads_coverage/LeeVV_2021.bw --normalizeUsing CPM --binSize 10 --minMappingQuality 10 --ignoreDuplicates -p max

computeMatrix reference-point --referencePoint TSS -R ../predicted_genes_TSS.dedup.bed ../v68_TSS.bed -S LeeVV_2021.bw  -b 500 -a 500 --binSize 10 --skipZeros --missingDataAsZero -p max -o LeeVV_2021_TSS_Matrix.gz
plotProfile -m LeeVV_2021_TSS_Matrix.gz --plotType se -T "Signal around TSS (±500bp)" -out predicted_TSS_LeeVV_2021.png --outFileNameData predicted_TSS_LeeVV_2021.tsv



computeMatrix scale-regions -S ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/ATAC-seq.bw \
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H2A.Zac.bw \
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H2A.Z.bw \
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H2B.Z.bw \
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K18ac.bw\
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K9ac.bw\
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K27ac.bw\
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K4me3.bw\
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K27me.bw\
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K4me1.bw\
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K4me2.bw\
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K4me.bw\
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3R17me2.bw\
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/MNase.bw\
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K18me.bw\
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K36me3.bw\
                                ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3.bw\
                            -R ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/Ref_gene.bed ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/intergenic.bed \
                            --beforeRegionStartLength  1000 --regionBodyLength 1500 --afterRegionStartLength 1000 --binSize 50 --skipZeros --missingDataAsZero -p 32 \
                            -o all_Epi_gene_matrix.gz 

plotProfile -m all_Epi_gene_matrix.gz -out all_Epi_gene.png --numPlotsPerRow 4  --plotTitle "Epigenetic patterns around genes" \
            --yMin 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \
            --yMax 8 20 20 12 20 1.5 10 6 2 6 5 4 5 6 5 5 5\
            --regionsLabel "genes" "intergenic" \
            --yAxisLabel "normalized coverage"\
            --legendLocation upper-right \
            --dpi 300 

computeMatrix scale-regions -S ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/ATAC-seq.bw \
                               ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/Cage.bw\
                               /rhome/zli529/lab/SRA_toolkit/RNAseq_7stages/processed_data/LT/RNA.bw\
                            -R ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/Ref_gene.bed ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/Pf3D7_v68_Refgenes.bed\
                            --beforeRegionStartLength  1000 --regionBodyLength 1500 --afterRegionStartLength 1000 --binSize 50 --skipZeros --missingDataAsZero -p 32 \
                            -o ATAC_RNA_Refgene_matrix.gz

plotProfile -m ATAC_RNA_Refgene_matrix.gz -out ATAC_RNA_v68_genes.png   --plotTitle "" \
            --yAxisLabel "normalized coverage"\
            --legendLocation upper-right \
            --dpi 300\
            --perGroup\
            --yMin 0 0\
            --yMax 8 8\
            --samplesLabel ATAC-seq Cage-seq RNA-seq\
            --regionsLabel  refined-genes PlasmoDB-v68-genes


plotProfile -m ATAC_RNA_Refgene_matrix.gz \
      --perGroup \
      --kmeans 2 \
      --plotType heatmap \
      -out ATAC_RNA_v68_genes_heatmap.png\
      --samplesLabel ATAC-seq Cage-seq RNA-seq\
      --dpi 300\



bamCoverage -b /rhome/zli529/lab/SRA_toolkit/GRO_seq/aligned/SRR4019517.sorted.q10.pp.bam \
  -o GRO.bw \
  --normalizeUsing CPM \
  --binSize 10 \
  --minMappingQuality 10 \
  --ignoreDuplicates \
  -p max