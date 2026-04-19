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
chr  <- "Pf3D7_08_v3"   # <-- put your seqname
pad <- 500L
from <- 945892 - pad      # start
to   <- 946599 + pad      # end

 

# promoter 
trk_ATAC <- make_bw_track(bw_file[1], chr, from, to,
                          name = "ATAC-seq", title_cex=0.2, color = "#E69F00", label_cex = 0.2)

trk_H2A.Zac <- make_bw_track(bw_file[2], chr, from, to,
                          name = "H2A.Zac",title_cex=0.2,  color = "#56B4E9", label_cex = 0.2)
trk_H2A.Z <- make_bw_track(bw_file[3], chr, from, to,
                          name = "H2A.Z", title_cex=0.2, color = "#009E73", label_cex = 0.2)
trk_H2B.Z <- make_bw_track(bw_file[4], chr, from, to,
                          name = "H2B.Z.bw", title_cex=0.2, color = "#F0E442", label_cex = 0.2)
trk_H3K18ac <- make_bw_track(bw_file[5], chr, from, to,
                          name = "H3K18ac", title_cex=0.2, color = "#0072B2", label_cex = 0.2)
trk_H3K27ac <- make_bw_track(bw_file[6], chr, from, to,
                          name = "H3K27ac", title_cex=0.2, color = "#D55E00", label_cex = 0.2)
trk_H3K4me3 <- make_bw_track(bw_file[7], chr, from, to,
                          name = "H3K4me3", title_cex=0.2, color = "#CC79A7", label_cex = 0.2)
trk_H3K9ac <- make_bw_track(bw_file[8], chr, from, to,
                          name = "H3K9ac", title_cex=0.2, color = "#1f77b4", label_cex = 0.2)

# gene body
trk_H3 <- make_bw_track(bw_file[9], chr, from, to,
                          name = "H3", title_cex=0.2, color = "#aec7e8", label_cex = 0.2)
trk_H3K18me <- make_bw_track(bw_file[10], chr, from, to,
                          name = "H3K18me", title_cex=0.2, color = "#ff7f0e", label_cex = 0.2)
trk_H3K27me <- make_bw_track(bw_file[11], chr, from, to,
                          name = "H3K27me", title_cex=0.2, color = "#ffbb78", label_cex = 0.2)
trk_H3K36me3 <- make_bw_track(bw_file[12], chr, from, to,
                          name = "H3K36me3", title_cex=0.2, color = "#2ca02c", label_cex = 0.2)
trk_H3K4me1 <- make_bw_track(bw_file[13], chr, from, to,
                          name = "H3K4me1", title_cex=0.2, color = "#98df8a", label_cex = 0.2)
trk_H3K4me2 <- make_bw_track(bw_file[14], chr, from, to,
                          name = "H3K4me2", title_cex=0.2, color = "#e377c2", label_cex = 0.2)
trk_H3K4me <- make_bw_track(bw_file[15], chr, from, to,
                          name = "H3K4me", title_cex=0.2, color = "#17becf", label_cex = 0.2)
trk_H3R17me2 <- make_bw_track(bw_file[16], chr, from, to,
                          name = "H3R17me2", title_cex=0.2, color = "#d62728", label_cex = 0.2)
trk_MNase <- make_bw_track(bw_file[17], chr, from, to,
                          name = "MNase", title_cex=0.2, color = "#9467bd", label_cex = 0.2)




# Quick preview plot of this single track:
png("novel_gene_track.png", width = 10, height = 3, units = "in", res = 300)
plotTracks(list(trk_ATAC, trk_H2A.Zac, trk_H2A.Z, trk_H2B.Z,trk_H3K18ac,trk_H3K27ac,trk_H3K4me3,trk_H3K9ac,
                trk_H3,trk_H3K18me, trk_H3K27me, trk_H3K36me3, trk_H3K4me1, trk_H3K4me2, trk_H3K4me,trk_H3R17me2,trk_MNase ), chromosome = chr, from = from, to = to,
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



# ---- Helper: robustly find the gene in a GFF by ID/Name/Parent ----
find_gene_in_gff <- function(gff_path, gene_id) {
  gff <- rtracklayer::import(gff_path)
  gene_like <- gff[gff$type %in% c("gene","protein_coding_gene","lncRNA_gene","RNA_gene","pseudogene")]
  hits <- gene_like[!is.na(mcols(gene_like)$ID) & as.character(mcols(gene_like)$ID) == gene_id]
  
  if (length(hits) == 0) {
    id_cols <- intersect(c("Name","gene_id","locus_tag"), colnames(mcols(gene_like)))
    if (length(id_cols)) {
      keep <- Reduce(`|`, lapply(id_cols, function(col) {
        x <- as.character(mcols(gene_like)[[col]]); !is.na(x) & x == gene_id
      }))
      hits <- gene_like[keep]
    }
  }
  if (length(hits) == 0) {
    tx <- gff[gff$type %in% c("mRNA","transcript","ncRNA","pseudogenic_transcript")]
    if ("Parent" %in% colnames(mcols(tx))) {
      p <- as.character(mcols(tx)$Parent)
      tx_hits <- tx[!is.na(p) & grepl(paste0("\\b", gene_id, "\\b"), p)]
      if (length(tx_hits) > 0) hits <- range(tx_hits)
    }
  }
  if (length(hits) == 0) stop("Gene ", gene_id, " not found in ", gff_path)
  hits[1]
}

# ---- Get region from v68 (for coordinates and seqname) ----
gene_v68 <- find_gene_in_gff(GFF_V68, gene_id)
chr  <- as.character(seqnames(gene_v68))
gsta <- start(gene_v68); gend <- end(gene_v68); gstr <- as.character(strand(gene_v68))
from <- max(1, gsta - pad); to <- gend + pad
region <- GRanges(chr, IRanges(from, to), strand = strand(gene_v68))


# ---- Helper: make UTR/CDS block track for a gene from a GFF ----
make_gene_blocks_track <- function(gff_path, gene_id, name_for_track, chr_hint = NULL) {
  gff <- rtracklayer::import(gff_path)
  
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
  if (length(gene_hits) == 0) {
    tx <- gff[gff$type %in% c("mRNA","transcript","ncRNA","pseudogenic_transcript")]
    if ("Parent" %in% colnames(mcols(tx))) {
      p <- as.character(mcols(tx)$Parent)
      tx_hits <- tx[!is.na(p) & grepl(paste0("\\b", gene_id, "\\b"), p)]
      if (length(tx_hits) > 0) gene_hits <- range(tx_hits)
    }
  }
  if (length(gene_hits) == 0) stop("Gene ", gene_id, " not found in ", gff_path)
  
  gene <- gene_hits[1]
  chr2 <- as.character(seqnames(gene))
  if (!is.null(chr_hint)) chr2 <- chr_hint  # keep on plotting chromosome if labels differ
  gsta2 <- start(gene); gend2 <- end(gene); gstr2 <- as.character(strand(gene))
  
  mrna <- gff[gff$type %in% c("mRNA","transcript","ncRNA","pseudogenic_transcript") &
                as.character(seqnames(gff)) == chr2]
  # transcript -> gene mapping
  pid <- mcols(mrna)$Parent
  if (is(pid, "CharacterList")) {
    pid <- vapply(pid, function(x) if (length(x)) x[[1]] else NA_character_, "")
  }
  gid <- mcols(mrna)$gene_id
  tx2gene <- ifelse(!is.na(pid), pid, gid)
  names(tx2gene) <- as.character(mcols(mrna)$ID)
  
  cds <- gff[gff$type == "CDS" & as.character(seqnames(gff)) == chr2]
  par <- mcols(cds)$Parent
  if (inherits(par, "CharacterList")) {
    par <- vapply(par, function(x) if (length(x)) x[[1]] else NA_character_, "")
  }
  cds <- cds[!is.na(par) & !is.na(tx2gene[par]) & tx2gene[par] == gene_id]
  has_cds <- length(cds) > 0
  if (has_cds) {
    cds_start <- min(start(cds)); cds_end <- max(end(cds))
  } else {
    cds_start <- cds_end <- NA_integer_
  }
  
  if (has_cds) {
    if (gstr2 == "+") {
      seg <- data.frame(
        start = c(gsta2,        cds_start,  cds_end + 1),
        end   = c(cds_start - 1, cds_end,    gend2),
        id    = c("5'UTR","CDS","3'UTR"),
        feature = c("five_prime_UTR","CDS","three_prime_UTR")
      )
    } else {
      seg <- data.frame(
        start = c(gsta2,        cds_start,  cds_end + 1),
        end   = c(cds_start - 1, cds_end,    gend2),
        id    = c("3'UTR","CDS","5'UTR"),
        feature = c("three_prime_UTR","CDS","five_prime_UTR")
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
  # color map per version
  if (grepl("v68", name_for_track, ignore.case = TRUE)) {
    fill_map <- c(three_prime_UTR="#7570b3", CDS="grey30", five_prime_UTR="#1b9e77", gene="#888888")
  } else {
    fill_map <- c(three_prime_UTR="#7570b3", CDS="grey30", five_prime_UTR="#1b9e77", gene="#888888")
  }
  displayPars(tr) <- list(
    featureAnnotation = "id", showFeatureId = TRUE, cex = 0.3,
    fill = fill_map,
    thinBoxFeature = c("five_prime_UTR","three_prime_UTR"),
    thinBoxHeight  = 1,
    col = NA, arrowHeadSize = 0, add53 = FALSE,
    cex.title = 0.5, fontcolor.title = "black"
  )
  tr
}

# ---- Build gene model tracks ----
track_v68 <- make_gene_blocks_track(GFF_V68, gene_id, name_for_track = paste0(gene_id, "_v68"), chr_hint = chr)
track_v48 <- make_gene_blocks_track(GFF_V48, gene_id, sprintf("%s_%s", gene_id, "v48"), chr_hint = chr)


# rotate the track names 90 degrees and slightly shrink the title text

Gviz::displayPars(track_v68) <- modifyList(
  Gviz::displayPars(track_v68),
  list( rotation.title = 360, cex.title= 0.3,title.just = c(1, 1), cex= 0.3, sizes_list=1,               
    fill  = c(three_prime_UTR="#7570b3", CDS="#ff7f0e",five_prime_UTR="#1b9e77",gene="#888888"),
    col = NA  ))

Gviz::displayPars(track_v48) <- modifyList(
  Gviz::displayPars(track_v48),
  list(rotation.title = 360,cex.title= 0.3,title.just= c(1,1),cex= 0.3,sizes_list=1, 
    fill= c(three_prime_UTR="#7570b3",CDS="#ff7f0e", five_prime_UTR="#1b9e77",gene="#888888"),
    col = NA))

png("ATAC_seq_track.png", width = 12, height = 6, units = "in", res = 300)

tracks <- list(trk_ATAC, trk_H2A.Zac, trk_H2A.Z, trk_H2B.Z, trk_H3K18ac, trk_H3K27ac,
               trk_H3K4me3, trk_H3K9ac, trk_H3, trk_H3K18me, trk_H3K27me, trk_H3K36me3,
               trk_H3K4me1, trk_H3K4me2, trk_H3K4me, trk_H3R17me2, trk_MNase, track_v68, track_v48)
sizes <- c(rep(0.5, 17), 0.3, 0.3)
plotTracks(
  tracks,
  chromosome = chr, from = from, to = to,
  sizes = sizes,
  background.title = "white", col.axis = "black", col.title = "black",
  title.width = 0.5        # uncomment if labels are clipped
  # , reverseStrand = (gstr == "-")
)
dev.off()



### add ONT-seq data ### 
ONT_BAM <- "/rhome/zli529/lab/SRA_toolkit/Long_read_RNA_datasets/LeeVV_2021_nanopore/SRR11094274.no_long_introns.sorted.bam"  # <-- edit; must have .bai
# 2) ONT per-read alignments (individual reads)
# (Make sure the BAM is indexed; if not: system(sprintf("samtools index %s", shQuote(ONT_BAM))))
ont_track <- AlignmentsTrack(
  ONT_BAM,
  isPaired = FALSE,
  name = "ONT-seq (coverage, reads)",
  stacking = "squish", # compact rows
  type = "pileup",
  showIndels = TRUE,
  add53 = FALSE                # no 5'/3' markers here
)

displayPars(ont_track) <- list(
  cex.title = 0.2,   
  fill.coverage = "#A0A0A0", # smaller title
  col.coverage    = NA ,
  fontcolor.title = "black"
)

png("ATAC_seq_track.png", width = 12, height = 6, units = "in", res = 300)

tracks <- list(trk_ATAC, trk_H2A.Zac, trk_H2A.Z, trk_H2B.Z, trk_H3K18ac, trk_H3K27ac,
               trk_H3K4me3, trk_H3K9ac, trk_H3, trk_H3K18me, trk_H3K27me, trk_H3K36me3,
               trk_H3K4me1, trk_H3K4me2, trk_H3K4me, trk_H3R17me2, trk_MNase, track_v68, track_v48, ont_track)
sizes <- c(rep(0.5, 17), 0.3, 0.3,4)
plotTracks(
  tracks,
  chromosome = chr, from = from, to = to,
  sizes = sizes,
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
                      id = "gene", feature = "gene", stringsAsFactors = FALSE)
  }
  
  if (is.null(name_for_track)) name_for_track <- paste0(gene_id, "_pred")
  
  tr <- Gviz::AnnotationTrack(
    start = seg$start, end = seg$end, chromosome = chr2, strand = as.character(strand(hit)),
    id = seg$id, feature = seg$feature,
    name = name_for_track, stacking = "dense", shape = "box"
  )
  
  # style: distinct color so it pops vs v48/v68
  fill_map <- c(exon = "#e377c2", gene = "#e377c2")  # magenta family
  Gviz::displayPars(tr) <- list(
    featureAnnotation = "id", showFeatureId = TRUE, cex = 0.3,
    fill = fill_map, thinBoxFeature = "exon", thinBoxHeight = 1,
    col = NA, arrowHeadSize = 0, add53 = FALSE,
    rotation.title = 360, cex.title = 0.3, title.just = c(1, 1)
  )
  tr
}

# ---- Build the predicted track and add it to your stack ----
track_pred <- make_predicted_track_from_bed(PRED_BED, gene_id,
                                            name_for_track = paste0(gene_id, "_pred"),
                                            chr_hint = chr)

# Insert it near the other gene models (before v68/v48) or wherever you prefer
tracks <- list(
  trk_ATAC, trk_H2A.Zac, trk_H2A.Z, trk_H2B.Z, trk_H3K18ac, trk_H3K27ac,
  trk_H3K4me3, trk_H3K9ac, trk_H3, trk_H3K18me, trk_H3K27me, trk_H3K36me3,
  trk_H3K4me1, trk_H3K4me2, trk_H3K4me, trk_H3R17me2, trk_MNase,
  track_pred, track_v68, track_v48, ont_track
)

# Give the predicted model a compact row like the others
sizes <- c(rep(0.5, 17), 0.3, 0.3,0.3, 3)  # add one more 0.3 for the new track

png("ATAC_seq_track.png", width = 12, height = 6, units = "in", res = 300)
plotTracks(
  tracks,
  chromosome = chr, from = from, to = to,
  sizes = sizes,
  background.title = "white", col.axis = "black", col.title = "black", 
  title.width = 0.5,
  col.baseline = "black"
  # , reverseStrand = (gstr == "-")
)
dev.off()
