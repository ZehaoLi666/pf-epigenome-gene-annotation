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
                            -R LncRNA_all.bed ~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/intergenic.bed \
                            --beforeRegionStartLength  500 --regionBodyLength 500 --afterRegionStartLength 500 --binSize 50 --skipZeros --missingDataAsZero -p 32 \
                            -o all_lncRNA_Epi_matrix.gz 


plotProfile -m all_lncRNA_Epi_matrix.gz -out all_lncRNA_Epi.png --numPlotsPerRow 4  --plotTitle "Epigenetic patterns around lncRNAs" \
            --yMin 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 \
            --yMax 8 25 25 20 20 1.5 15 6 2 6 5 4 5 6 5 5 5\
            --regionsLabel "lncRNAs" "intergenic" \
            --yAxisLabel "normalized coverage"\
            --legendLocation upper-right \
            --dpi 300 