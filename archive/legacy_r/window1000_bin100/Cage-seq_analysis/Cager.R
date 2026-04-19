setwd("/rhome/zli529/lab/LncRNA_chip_prediction/Final/window1000_bin100/Cage-seq_analysis/") 

# BiocManager::install(c("CAGEr","GenomicRanges","rtracklayer"))
library(CAGEr)
library(GenomicRanges)
library(rtracklayer)
#library(Gviz)
##library(CAGEfightR)
#library(GenomicInteractions)
#library(Gviz)

library(BSgenome.Pfalciparum.PlasmoDB.v24)

library(CAGEr)

# Path to your BAM files
bam_dir <- "~/lab/SRA_toolkit/cage_seq_data/paired-end/aligned_hisat2_pe"

bam_files <- list.files(bam_dir, pattern="*.sorted.bam$", full.names=TRUE)

# Create a sample table
# Assign meaningful sample names instead of SRR IDs if you like
sample.info <- data.frame(
  name = sub(".sorted.bam","", basename(bam_files)),
  bamFiles = bam_files,
  stringsAsFactors = FALSE
)

# Initialize CAGEset object
cage <- CAGEexp(
  genomeName     = "BSgenome.Pfalciparum.PlasmoDB.v24",   # use the genome you have installed
  inputFiles     = "/rhome/zli529/lab/SRA_toolkit/cage_seq_data/paired-end/aligned_hisat2_pe/SRR2031965.sorted.bam",
  inputFilesType = "bam",
  sampleLabels   = "SRR2031965"
)

# Perform standard steps
ce <- getCTSS(cage, useMulticore = TRUE, nrCores = 4)                               # extract CTSS from BAMs

## Sanity check: you should now see a 'counts' assay
assays(ce)           # expect includes "counts"
sum(assay(ce, "counts"))  # > 0

## 3) Normalize to TPM (the '.simpleTpm' path your error mentioned)
ce <- normalizeTagCount(ce, method = "simpleTpm")

## Now downstream calls will work, e.g.:
tpm <- CTSStagCountTPM(ce)   


clusterCTSS(cage, threshold=1, thresholdIsTpm=TRUE,
            nrPassThreshold=1, method="distclu", maxDist=20, removeSingletons=TRUE)

exportCTSStoBedGraph(cage, values="normalized", format="bedGraph", oneFile=FALSE)


