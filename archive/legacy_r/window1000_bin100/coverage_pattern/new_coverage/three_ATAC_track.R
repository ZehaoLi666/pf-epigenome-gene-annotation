setwd("/rhome/zli529/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage")


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
                          rotate_title = 90) {  # 0 or 90
  if (is.null(name)) name <- tools::file_path_sans_ext(basename(bw_path))
  region <- GenomicRanges::GRanges(chr, IRanges::IRanges(from, to))
  bw_gr  <- rtracklayer::import(bw_path, which = region, as = "GRanges")
  if (length(bw_gr) == 0) stop("No data in region for: ", bw_path)
  
  ymax <- max(1, stats::quantile(bw_gr$score, 0.99, na.rm = TRUE))
  tr <- Gviz::DataTrack(range = bw_gr, type = "histogram",
                        name = name, baseline = 0, ylim = c(0, ymax))
  
  y_ticks <- pretty(c(0, ymax), n = 5)
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

bw_file <- c("~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/ATAC-seq.bw" ,
             "/rhome/zli529/lab/SRA_toolkit/ATAC-seq_datasets/Parul_2025/aligned/ATAC_sub.bw",
             "/rhome/zli529/lab/SRA_toolkit/ATAC-seq_datasets/Ruiz_2018/bigwig/ATAC_sub.bw" 
             
)

# promoter 
Toenhake_2018 <- make_bw_track(bw_file[1], chr, from, to,
                          name = "Toenhake_2018", title_cex=0.2, color = "#E69F00", label_cex = 0.2)

Parul_2025 <- make_bw_track(bw_file[2], chr, from, to,
                             name = "Parul_2025",title_cex=0.2,  color = "#56B4E9", label_cex = 0.2)
Ruiz_2018 <- make_bw_track(bw_file[3], chr, from, to,
                           name = "Ruiz_2018", title_cex=0.2, color = "#009E73", label_cex = 0.2)
displayPars(Toenhake_2018) <- list(col.baseline = "black")
displayPars(Parul_2025)    <- list(col.baseline = "black")
displayPars(Ruiz_2018)     <- list(col.baseline = "black")


# Quick preview plot of this single track:
png("three_ATAC_track.png", width = 10, height = 2, units = "in", res = 600)
plotTracks(list(Toenhake_2018, Parul_2025, Ruiz_2018), chromosome = chr, from = from, to = to,
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
                             name = "Genes (v48)",
                             group = lbl,
                             id    = lbl,
                             stacking = "dense",
                             shape = "box",
                             fill  = "#666666", col = NA,
                             cex.title = 0.6, cex.group = 0.5, fontsize = 8
)


png("ATAC_plus_genes_v48.png", width = 12, height = 3.2, units = "in", res = 600)
plotTracks(
  list(  Toenhake_2018, Parul_2025, Ruiz_2018, geneTrack),
  chromosome = chr, from = from, to = to,
  sizes = c( 0.2, 0.8, 0.8, 0.8),      # relative heights
  background.title = "white", col.axis = "black", col.title = "black",
  title.width = 0.6
)
dev.off()

