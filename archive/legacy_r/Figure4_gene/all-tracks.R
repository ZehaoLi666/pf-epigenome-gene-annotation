setwd("/rhome/zli529/lab/LncRNA_chip_prediction/Final/Figure4_gene")


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
chr  <- "Pf3D7_11_v3"   # <-- put your seqname
pad <- 500L
from <- 509191 - pad      # start
to   <- 516497 + pad      # end

bw_file <- c("~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/ATAC-seq.bw" ,
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H2A.Zac.bw", 
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H2A.Z.bw",
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H2B.Z.bw",
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K18ac.bw",
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K27ac.bw",
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K4me3.bw",
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K9ac.bw", 
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3.bw",
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K18me.bw",
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K27me.bw",
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K36me3.bw",
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K4me1.bw",
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K4me2.bw",
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3K4me.bw",
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/H3R17me2.bw",
             "~/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/MNase.bw")

# promoter 
ATAC <- make_bw_track(bw_file[1], chr, from, to, name = " ", title_cex = 0.1, color = "#E69F00",  label_cex = 0.5)
H2A.Zac     <- make_bw_track(bw_file[2], chr, from, to, name = " ",  title_cex = 0.1, color = "#56B4E9", label_cex = 0.5)
H2A.Z     <- make_bw_track(bw_file[3], chr, from, to,name = " ",  title_cex = 0.1, color = "#009E73",  label_cex = 0.5)
H2B.Z     <- make_bw_track(bw_file[4], chr, from, to, name = " ", title_cex = 0.1, color = "#F0E442",  label_cex = 0.5)
H3K18ac     <- make_bw_track(bw_file[5], chr, from, to, name = " ", title_cex = 0.1, color = "#0072B2",  label_cex = 0.5)
H3K27ac     <- make_bw_track(bw_file[6], chr, from, to,name = " ",  title_cex = 0.1,  color = "#D55E00",  label_cex = 0.5)
H3K4me3     <- make_bw_track(bw_file[7], chr, from, to,name = " ", title_cex = 0.1,color = "#CC79A7", label_cex = 0.5)
H3K9ac     <- make_bw_track(bw_file[8], chr, from, to,name = " ", title_cex = 0.1,color = "#1f77b4", label_cex = 0.5)


# gene body
H3     <- make_bw_track(bw_file[9], chr, from, to, name = " ",  title_cex = 0.1,color = "#aec7e8", label_cex = 0.5)
H3K18me     <- make_bw_track(bw_file[10], chr, from, to,name = " ", title_cex = 0.1, color = "#ff7f0e", label_cex = 0.5)
H3K27me     <- make_bw_track(bw_file[11], chr, from, to, name = " ",  title_cex = 0.1, color = "#ffbb78",  label_cex = 0.5)
H3K36me3     <- make_bw_track(bw_file[12], chr, from, to,name = " ", title_cex = 0.1, color = "#2ca02c", label_cex = 0.5)
H3K4me1     <- make_bw_track(bw_file[13], chr, from, to,name = " ",title_cex = 0.1, color = "#98df8a", label_cex = 0.5)
H3K4me2     <- make_bw_track(bw_file[14], chr, from, to,name = " ", title_cex = 0.1, color = "#e377c2", label_cex = 0.5)
H3K4me     <- make_bw_track(bw_file[15], chr, from, to,name = " ", title_cex = 0.1, color = "#17becf", label_cex = 0.5)
H3R17me2     <- make_bw_track(bw_file[2], chr, from, to,name = " ",  title_cex = 0.1,  color = "#d62728",  label_cex = 0.5)
MNase     <- make_bw_track(bw_file[16], chr, from, to, name = " ", title_cex = 0.1, color = "#9467bd",  label_cex = 0.5)

displayPars(ATAC) <- list(col.baseline = "black")
displayPars(H2A.Zac)     <- list(col.baseline = "black")
displayPars(H2A.Z) <- list(col.baseline = "black")
displayPars(H2B.Z)     <- list(col.baseline = "black")
displayPars(H3K18ac) <- list(col.baseline = "black")
displayPars(H3K27ac)     <- list(col.baseline = "black")
displayPars(H3K4me3) <- list(col.baseline = "black")
displayPars(H3K9ac)     <- list(col.baseline = "black")
displayPars(H3) <- list(col.baseline = "black")
displayPars(H3K18me)     <- list(col.baseline = "black")
displayPars(H3K27me) <- list(col.baseline = "black")
displayPars(H3K36me3)     <- list(col.baseline = "black")
displayPars(H3K4me1)     <- list(col.baseline = "black")
displayPars(H3K4me2)     <- list(col.baseline = "black")
displayPars(H3K4me)     <- list(col.baseline = "black")
displayPars(H3R17me2)     <- list(col.baseline = "black")
displayPars(MNase)     <- list(col.baseline = "black")


# Quick preview plot of this single track:
png("epi_track.png", width = 10, height = 12, units = "in", res = 800)
plotTracks(list(ATAC, H2A.Zac, H2A.Z, H2B.Z,H3K18ac,H3K27ac,H3K4me3,H3K9ac,
                H3,H3K18me, H3K27me, H3K36me3, H3K4me1, H3K4me2, H3K4me,H3R17me2,MNase ), chromosome = chr, from = from, to = to,
           title.width = 0.5, background.title = "white",
           col.axis = "black", col.title = "black")

dev.off()







### gene model track ###
library(Gviz)
library(rtracklayer)
library(GenomicRanges)
options(ucscChromosomeNames = FALSE)

# ---- Inputs (edit if paths differ) ----
gene_id <- "PF3D7_1113100"
pad     <- 1000L

GFF_V68 <- "/rhome/zli529/lab/PlasmoDB_Genome/PlasmnDB_v68/PlasmoDB-68_Pfalciparum3D7.gff"
GFF_V48 <- "/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7.gff"



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
    featureAnnotation = "id", showFeatureId = TRUE, cex = 1,
    fill = fill_map,
    thinBoxFeature = c("five_prime_UTR","three_prime_UTR"),
    thinBoxHeight  = 0.5,
    col = NA, arrowHeadSize = 0, add53 = FALSE,
    background.title = "white", col.title = "black", rotation.title = 90
  )
  tr
}

track_gene <- make_gene_blocks_track(GFF_V68, GENE_ID, name_for_track = paste0(GENE_ID, " (v68)"), chr_hint = chr_target)

# ---- Build gene model tracks ----
track_v68 <- make_gene_blocks_track(GFF_V68, gene_id, name_for_track = "PlasmoDB_v68", chr_hint = chr)
track_v48 <- make_gene_blocks_track(GFF_V48, gene_id, name_for_track = "PlasmoDB_v48", chr_hint = chr)


# Quick preview plot of this single track:
png("gene_v48_v68_track.png", width = 12, height = 1, units = "in", res = 800)
plotTracks(list(track_v48, track_v68 ), chromosome = chr, from = from, to = to,
           title.width = 0.5, background.title = "white",
           col.axis = "black", col.title = "black")

dev.off()





### add ONT-seq data ### 
ONT_BAM <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/Figure4_gene/Pf3D7_11_v3_subset.bam"  # <-- edit; must have .bai
# 2) ONT per-read alignments (individual reads)
# (Make sure the BAM is indexed; if not: system(sprintf("samtools index %s", shQuote(ONT_BAM))))
ont_track <- AlignmentsTrack(
  ONT_BAM,
  isPaired = FALSE,
  name = " ",
  stacking = "full", # compact rows
  type = "pileup",
  min.height = 30,lwd = 1.5,
  showIndels = TRUE,
  add53 = FALSE ,
  stackHeight = 10 ,
  col = NA,
)

displayPars(ont_track) <- list(
  cex.title = 1,   
  fill.coverage = "#A0A0A0", # smaller title
  col.coverage    = NA ,
  fontcolor.title = "black",
  list(min.height = 3)
)

png("ONT_track.png", width = 12, height = 2, units = "in", res = 600)

tracks <- list(ont_track)
plotTracks(
  tracks,
  chromosome = chr, from = from, to = to,
  
  background.title = "white", col.axis = "black", col.title = "black",
  title.width = 0.5        # uncomment if labels are clipped
  # , reverseStrand = (gstr == "-")
)
dev.off()



### predicted gene model track## 
PRED_BED <- "predicted_gene.bed"



# ---- Helper: build a predicted-gene track from BED ----
make_predicted_track_from_bed <- function(bed_path, gene_id, name_for_track = NULL, chr_hint = NULL) {
  stopifnot(file.exists(bed_path))
  bed <- rtracklayer::import(bed_path, format = "BED")
  if (!"name" %in% colnames(mcols(bed))) mcols(bed)$name <- NA_character_
  
  # Prefer exact name match; else allow partial
  hit <- bed[!is.na(mcols(bed)$name) & as.character(mcols(bed)$name) == gene_id]
  if (length(hit) == 0) {
    hit <- bed[!is.na(mcols(bed)$name) & grepl(paste0("\\b", gene_id, "\\b"), as.character(mcols(bed)$name))]
  }
  if (length(hit) == 0) stop("Predicted gene ", gene_id, " not found in ", bed_path)
  
  hit <- hit[1]
  chr2 <- as.character(seqnames(hit))
  if (!is.null(chr_hint)) chr2 <- chr_hint
  
  # If BED12, expand exon blocks; else single block
  blk_sizes  <- mcols(hit)$blockSizes
  blk_starts <- mcols(hit)$blockStarts
  
  if (!is.null(blk_sizes) && !is.na(blk_sizes) && nchar(blk_sizes) > 0) {
    sizes  <- as.integer(strsplit(as.character(blk_sizes), ",")[[1]])
    starts <- as.integer(strsplit(as.character(blk_starts), ",")[[1]])
    starts_abs <- start(hit) + starts
    ends_abs   <- starts_abs + sizes - 1L
    seg <- data.frame(start = starts_abs, end = ends_abs,
                      id = paste0("exon", seq_along(starts_abs)),
                      feature = "exon", stringsAsFactors = FALSE)
  } else {
    seg <- data.frame(start = start(hit), end = end(hit),
                      id = "Refined gene model", feature = "Gene", stringsAsFactors = FALSE)
  }
  
  if (is.null(name_for_track)) name_for_track <- " "
  
  tr <- Gviz::AnnotationTrack(
    start = seg$start, end = seg$end, chromosome = chr2, strand = as.character(strand(hit)),
    id = seg$id, feature = seg$feature,
    name = name_for_track, stacking = "dense", shape = "box"
  )
  
  # style: distinct color so it pops vs v48/v68
  fill_map <- c(exon = "#e377c2", gene = "#e377c2")  # magenta family
  Gviz::displayPars(tr) <- list(
    featureAnnotation = "id", showFeatureId = TRUE, cex = 1,
    fill = fill_map, thinBoxFeature = "exon", thinBoxHeight = 1,
    col = NA, arrowHeadSize = 0, add53 = FALSE,
    rotation.title = 360, cex.title = 1, title.just = c(1, 1)
  )
  tr
}

# ---- Build the predicted track and add it to your stack ----
track_pred <- make_predicted_track_from_bed(PRED_BED, gene_id,
                                            name_for_track = " ",
                                            chr_hint = chr)


# Quick preview plot of this single track:
png("track_pred.png", width = 12, height = 0.5, units = "in", res = 800)
plotTracks(list(track_pred ), chromosome = chr, from = from, to = to,
           title.width = 0.5, background.title = "white",
           col.axis = "black", col.title = "black")

dev.off()





################ all tracks ##########################

# Quick preview plot of this single track:
png("ALl_track.png", width = 12, height = 12, units = "in", res = 600)
plotTracks(list(ATAC, H2A.Zac, H2A.Z, H2B.Z,H3K18ac,H3K27ac,H3K4me3,H3K9ac,
                H3,H3K18me, H3K27me, H3K36me3, H3K4me1, H3K4me2, H3K4me,H3R17me2,MNase,
                track_pred,track_v48, track_v68, ont_track  ), 
           sizes = c( 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6,0.6, 0.6,0.6,
                      0.3,0.3,0.3,2.5),
           chromosome = chr, from = from, to = to,
           title.width = 0.5, background.title = "white",
           col.axis = "black", col.title = "black")

dev.off()

