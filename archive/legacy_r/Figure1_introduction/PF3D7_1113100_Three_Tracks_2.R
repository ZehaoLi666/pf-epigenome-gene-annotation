# --- 1. Setup and Libraries ---
setwd("/rhome/zli529/lab/LncRNA_chip_prediction/Final/Figure1_introduction/")
library(Gviz)
library(rtracklayer)
library(GenomicRanges)
options(ucscChromosomeNames = FALSE)

# --- 2. Define Parameters and Paths ---
# Coordinates
chr_target <- "Pf3D7_11_v3"
from_pos   <- 509191
to_pos     <- 516497

# File Paths
#RNA_BED  <- "/rhome/zli529/lab/LncRNA_chip_prediction/NEW/GRO-RNA_coverage/filtered_bed_files/filtered_GSE85478_GametocyteV_scaleFactor-0.201.bed"
RNA_BED  <- "/rhome/zli529/lab/LncRNA_chip_prediction/NEW/GRO-RNA_coverage/filtered_bed_files/filtered_GSE85478_GROSeq_LT_scaleFactor-0.202.bed"

# *** UPDATED ONT-SEQ PATH ***
ONT_BAM  <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/Figure1_introduction/minimap2_ONT.sorted.bam"
#ONT_BAM  <- "/rhome/zli529/lab/SRA_toolkit/Long_read_RNA_datasets/LeeVV_2021_nanopore/SRR11094274.no_long_introns.sorted.bam"

GFF_V68  <- "/rhome/zli529/lab/PlasmoDB_Genome/PlasmnDB_v68/PlasmoDB-68_Pfalciparum3D7.gff"
#GENE_ID  <- "PF3D7_1476600"
GENE_ID  <- "PF3D7_1113100"


# --- 3. Track 1: RNA-seq Coverage (Blue) ---
# Custom reader for your 3-column format: [Chrom] [Pos] [Score]
rna_df <- read.table(RNA_BED, header = FALSE, stringsAsFactors = FALSE)

# Create GRanges assuming Col2 is both Start and End (single base position)
rna_gr <- GRanges(
  seqnames = rna_df[, 1],
  ranges   = IRanges(start = rna_df[, 2], end = rna_df[, 2]),
  score    = as.numeric(rna_df[, 3])
)

# Subset to the specific region immediately to reduce memory usage
rna_gr <- subsetByOverlaps(rna_gr, GRanges(chr_target, IRanges(from_pos, to_pos)))

track_rnaseq <- DataTrack(
  range = rna_gr,
  genome = "Pf3D7",
  name = " ",
  type = "polygon",      
  chromosome = chr_target,
  col = "grey30",       
  fill.mountain = c("grey30", "grey30"),
  background.title = "white",
  col.title = "black",
  col.axis = "black",
  cex.title = 2,
  rotation.title = 0,
  margin = 15 ,
  title.width = 2.5 ,
  # --- FONT SIZE 30 SETUP ---
  fontsize = 20,             # Sets the base font size to 30 (matches ggplot base_size)
  cex.title = 1.0,           # 1.0 means exactly 30pt. Increase to 1.2 for "extra large"
  cex.axis = 0.8,
  fontface.title = 1,
  # --- UPDATED Y-AXIS 0-100 ---
  ylim = c(0, 100),
  yTicksAt = c(0, 25, 75, 100),
)


### ONT-seq track 

# 1. Coverage Track (The histogram at the top)
track_ont_cov <- AlignmentsTrack(
  ONT_BAM,
  isPaired = FALSE,
  chromosome = chr_target, # Explicitly set chromosome
  name = " ",
  type = "coverage",
  fill.coverage = "skyblue3",
  col.coverage = NA,
  # --- Add these for the title formatting ---
  background.title = "transparent", # Or "transparent" to match the RNA-seq look
  col.title = "black",           # Changes text from white to black
  cex.title = 2,
  title.width = 2.5,
  rotation.title = 0,
  margin = 15  ,
  # --- FONT SIZE 30 SETUP ---
  fontsize = 20,             # Sets the base font size to 30 (matches ggplot base_size)
  cex.title = 1.0,           # 1.0 means exactly 30pt. Increase to 1.2 for "extra large"
  cex.axis = 0.8,
  fontface.title = 1,
  # --- UPDATED Y-AXIS 0-100 ---
  ylim = c(0, 100),
  yTicksAt = c(0, 25, 75, 100),
  
)

# 2. Reads Track (The individual thick reads)
track_ont_reads <- AlignmentsTrack(
  ONT_BAM,
  isPaired = FALSE,
  chromosome = chr_target, # Explicitly set chromosome
  name = " ",
  type = "pileup",         # Only shows reads, no histogram
  stacking = "full",
  min.height = 50,         # This makes them thick!
  #max.height = 30,
  col.reads = "white",     # Outline of the read
  fill = "#A0A0A0"   ,      # Color of the read body
  max.draw.reads = 40,
  showIndels = TRUE,
  showScrollbars = FALSE,
  # --- Add these for the title formatting ---
  background.title = "transparent", # Or "transparent" to match the RNA-seq look
  col.title = "black",           # Changes text from white to black
  cex.title = 2,
  title.width = 2.5,
  rotation.title = 0,
  margin = 15 
)



### Track 3: Gene Model (Helper Function) --- 

make_gene_blocks_track <- function(gff_path, gene_id, name_for_track, chr_hint = NULL) {
  gff <- rtracklayer::import(gff_path)
  
  # Find gene
  gene_like <- gff[gff$type %in% c("gene","protein_coding_gene","lncRNA_gene","RNA_gene","pseudogene")]
  gene_hits <- gene_like[!is.na(mcols(gene_like)$ID) & as.character(mcols(gene_like)$ID) == gene_id]
  
  if (length(gene_hits) == 0) {
    id_cols <- intersect(c("Name","gene_id","locus_tag"), colnames(mcols(gene_like)))
    if (length(id_cols)) {
      keep <- Reduce(`|`, lapply(id_cols, function(col) {
        x <- as.character(mcols(gene_like)[[col]]); !is.na(x) & x == gene_id
      }))
      gene_hits <- gene_like[keep]
    }
  }
  if (length(gene_hits) == 0) stop("Gene ", gene_id, " not found in ", gff_path)
  
  gene <- gene_hits[1]
  chr2 <- as.character(seqnames(gene))
  if (!is.null(chr_hint)) chr2 <- chr_hint 
  gsta2 <- start(gene); gend2 <- end(gene); gstr2 <- as.character(strand(gene))
  
  # Find CDS
  cds <- gff[gff$type == "CDS" & as.character(seqnames(gff)) == chr2]
  par <- mcols(cds)$Parent
  if (inherits(par, "CharacterList")) par <- vapply(par, function(x) x[[1]], "")
  cds_hits <- cds[ !is.na(par) & grepl(gene_id, par) ] 
  
  has_cds <- length(cds_hits) > 0
  
  if (has_cds) {
    cds_start <- min(start(cds_hits))
    cds_end   <- max(end(cds_hits))
    
    if (gstr2 == "+") {
      seg <- data.frame(
        start = c(gsta2,         cds_start,   cds_end + 1),
        end   = c(cds_start - 1, cds_end,     gend2),
        id    = c("5'UTR",       "CDS",       "3'UTR"),
        feature = c("five_prime_UTR", "CDS", "three_prime_UTR")
      )
    } else {
      seg <- data.frame(
        start = c(gsta2,         cds_start,   cds_end + 1),
        end   = c(cds_start - 1, cds_end,     gend2),
        id    = c("3'UTR",       "CDS",       "5'UTR"),
        feature = c("three_prime_UTR", "CDS", "five_prime_UTR")
      )
    }
    seg <- seg[seg$end >= seg$start, ]
  } else {
    seg <- data.frame(start = gsta2, end = gend2, id = "gene", feature = "gene")
  }
  
  tr <- AnnotationTrack(
    start = seg$start, end = seg$end, chromosome = chr2, strand = TRUE,
    id = seg$id, feature = seg$feature,
    name = name_for_track, stacking = "dense", shape = "box"
  )
  
  fill_map <- c(three_prime_UTR="#7570b3", CDS="grey30", five_prime_UTR="#1b9e77", gene="#888888")
  
  displayPars(tr) <- list(
    featureAnnotation = "id", showFeatureId = TRUE, cex = 1.2,
    fill = fill_map,
    thinBoxFeature = c("five_prime_UTR","three_prime_UTR"),
    thinBoxHeight  = 0.5,
    col = NA, arrowHeadSize = 0, add53 = FALSE,
    background.title = "white", col.title = "black", rotation.title = 90
  )
  tr
}

track_gene <- make_gene_blocks_track(GFF_V68, GENE_ID, name_for_track = paste0(GENE_ID, " (v68)"), chr_hint = chr_target)





# --- 6. Plotting ---
png("PF3D7_1113100_Three_Tracks_2.png", width = 20, height =8, units = "in", res = 600)

tracks_to_plot <- list(track_rnaseq, track_ont_cov,track_ont_reads, track_gene)
track_sizes <- c(1, 1,2, 0.2)

plotTracks(
  tracks_to_plot,
  chromosome = chr_target,
  from = from_pos,
  to = to_pos,
  sizes = track_sizes,
  title.width = 2,
  col.axis = "black",
  col.border.title = "transparent",
  showScrollbars = FALSE
)

dev.off()

