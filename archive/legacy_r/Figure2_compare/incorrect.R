setwd("/bigdata/lerochlab/zli529/LncRNA_chip_prediction/Final/Figure2_compare")


# --- deps ---
library(Gviz)
library(rtracklayer)
library(GenomicRanges)
options(ucscChromosomeNames = FALSE)

# ---------- function: build ONE BigWig track ----------
make_bw_track <- function(bw_path, chr, from, to,
                          name = NULL,
                          color = "#E69F00",
                          label_cex = 1.0,      # y-axis tick size
                          title_cex = 0.8,      # track name size
                          rotate_title = 90,    # 0 or 90
                          ...) {                # <--- 1. Added "..." to accept extra args
  
  if (is.null(name)) name <- tools::file_path_sans_ext(basename(bw_path))
  
  region <- GenomicRanges::GRanges(chr, IRanges::IRanges(from, to))
  bw_gr  <- rtracklayer::import(bw_path, which = region, as = "GRanges")
  
  if (length(bw_gr) == 0) stop("No data in region for: ", bw_path)
  
  # <--- 2. Fix: Force negative values to 0 to prevent plotting below x-axis
  bw_gr$score[bw_gr$score < 0] <- 0
  
  # <--- 3. Logic: Check if user passed 'ylim' in '...'. If not, calculate it.
  args <- list(...)
  if ("ylim" %in% names(args)) {
    final_ylim <- args$ylim
  } else {
    ymax <- max(1, stats::quantile(bw_gr$score, 0.99, na.rm = TRUE))
    final_ylim <- c(0, ymax)
  }
  
  # Create the track, passing "..." to Gviz
  tr <- Gviz::DataTrack(range = bw_gr, type = "histogram",
                        name = name, baseline = 0, 
                        ylim = final_ylim, 
                        ...) # <--- Pass extra args (like showTitle) here
  
  # Setup display parameters
  # Note: yTicks are only calculated if using the auto-generated ylim to avoid conflicts
  if (!"ylim" %in% names(args)) {
    y_ticks <- pretty(final_ylim, n = 5)
  } else {
    y_ticks <- pretty(args$ylim, n = 5)
  }
  
  Gviz::displayPars(tr) <- modifyList(Gviz::displayPars(tr), list(
    fill.histogram = color, col.histogram = NA,
    cex.axis = label_cex, yTicksAt = y_ticks, yTickLabels = y_ticks,
    rotation.title = rotate_title, cex.title = title_cex, fontcolor.title = "black"
  ))
  
  tr
}

# If you already have chr/from/to from your gene ± pad, reuse them.
# If not, set them manually like this:
chr  <- "Pf3D7_14_v3"   # <-- put your seqname
pad <- 500L
from <- 1482677 - pad      # start
to   <- 1606261 + pad      # end



bw_file <- c("/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Krishanpal_Chip/substract_readcounts/by_antibody/bigwig/H3K14ac.bw" ,
             "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Krishanpal_Chip/substract_readcounts/by_antibody/bigwig/H3K27ac.bw", 
             "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Krishanpal_Chip/substract_readcounts/by_antibody/bigwig/H3K9ac.bw",
             "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Krishanpal_Chip/substract_readcounts/by_antibody/bigwig/H3K4me3.bw",
             "/rhome/zli529/lab/SRA_toolkit/Chip-seq_datasets/Krishanpal_Chip/substract_readcounts/by_antibody/bigwig/H3K9me3.bw"
             
)



H3K14ac <- make_bw_track(bw_file[1], chr, from, to, name = " ", title_cex = 0.1, color = "#E69F00",  label_cex = 0.5)
H3K27ac     <- make_bw_track(bw_file[2], chr, from, to, name = " ",  title_cex = 0.1, color = "#56B4E9", label_cex = 0.5)
H3K9ac     <- make_bw_track(bw_file[3], chr, from, to,name = " ",  title_cex = 0.1, color = "#009E73",  label_cex = 0.5)
H3K4me3     <- make_bw_track(bw_file[4], chr, from, to, name = " ", title_cex = 0.1, color = "#F0E442",  label_cex = 0.5)
H3K9me3     <- make_bw_track(bw_file[5], chr, from, to, name = " ", title_cex = 0.1, color = "#0072B2",  label_cex = 0.5)


displayPars(H3K14ac) <- list(col.baseline = "black")
displayPars(H3K27ac) <- list(col.baseline = "black")
displayPars(H3K9ac) <- list(col.baseline = "black")
displayPars(H3K4me3) <- list(col.baseline = "black")
displayPars(H3K9me3) <- list(col.baseline = "black")


png("incorrect.png", width = 10, height = 3, units = "in", res = 600)
plotTracks(list(H3K14ac, H3K27ac, H3K9ac, H3K4me3,H3K9me3 ), 
           chromosome = chr, from = from, to = to,
           title.width = 0.5, background.title = "white",
           col.axis = "black", col.title = "black")

dev.off()




library(Gviz)
library(rtracklayer)
library(GenomicRanges)

GFF_V48 <- "/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7.gff"

region <- GRanges(chr, IRanges(from, to))
gff_gr <- import(GFF_V48, which = region)          # GFF3
genes  <- gff_gr[gff_gr$type == "gene"]

# Use gene Name (if present) otherwise ID
lbl <- if ("Name" %in% colnames(mcols(genes))) genes$Name else genes$ID

geneTrack <- AnnotationTrack(genes,
                             name = " ",
                             group = lbl,
                             id    = lbl,
                             stacking = "dense",
                             shape = "box",
                             fill  = "#666666", col = NA,
                             cex.title = 0.6, cex.group = 0.5, fontsize = 8
)


png("incorrect_v48.png", width = 12, height = 4, units = "in", res = 600)
plotTracks(
  list(  H3K14ac, H3K27ac, H3K9ac,H3K4me3, H3K9me3,geneTrack),
  chromosome = chr, from = from, to = to,
  sizes = c( 0.8, 0.8, 0.8 ,0.8,0.8,0.2),      # relative heights
  background.title = "white", col.axis = "black", col.title = "black",
  title.width = 0.6
)
dev.off()


