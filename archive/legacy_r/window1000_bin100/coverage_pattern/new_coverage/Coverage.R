# install.packages("ggcoverage")           # if needed
# BiocManager::install(c("rtracklayer","GenomicRanges"))  # if needed
setwd("/rhome/zli529/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/PF3D7_1113100")
library(ggcoverage)
library(rtracklayer)
library(GenomicRanges)
library(Gviz)
library(rtracklayer)
library(GenomicRanges)

gene_id <- "PF3D7_1113100"
gff <- import(GFF)

# allow multiple gene-like types in v68
gene_like <- gff[gff$type %in% c("gene","protein_coding_gene","lncRNA_gene","RNA_gene","pseudogene")]

# match by any of the usual ID fields
id_cols <- intersect(c("ID","Name","gene_id","locus_tag"), colnames(mcols(gene_like)))
keep <- Reduce(`|`, lapply(id_cols, function(col) {
  x <- as.character(mcols(gene_like)[[col]])
  !is.na(x) & x == gene_id
}))
gene_hits <- gene_like[keep]

# fallback: try Parent field (some GFFs anchor transcripts to gene ID)
if (length(gene_hits) == 0 && "Parent" %in% colnames(mcols(gff))) {
  tx <- gff[gff$type %in% c("mRNA","transcript","ncRNA","pseudogenic_transcript")]
  p <- as.character(mcols(tx)$Parent)
  tx_hits <- tx[!is.na(p) & grepl(paste0("\\b", gene_id, "\\b"), p)]
  if (length(tx_hits) > 0) gene_hits <- range(tx_hits)
}

stopifnot(length(gene_hits) >= 1)

# build padded region for plotting
pad <- 1000L
region <- GRanges(
  seqnames = seqnames(gene_hits)[1],
  ranges   = IRanges(start(gene_hits)[1] - pad, end(gene_hits)[1] + pad),
  strand   = strand(gene_hits)[1]
)

region


RNA_BW  <- "~/lab/SRA_toolkit/RNAseq_7stages/processed_data/LT/RNA.bw"
ONT_BAM <- "/rhome/zli529/lab/SRA_toolkit/Long_read_RNA_datasets/LeeVV_2021_nanopore/SRR11094274.no_long_introns.sorted.bam"  # needs .bai

gff <- import(GFF)

stopifnot(length(gene_hits) >= 1)
gene <- gene_hits[1]
pad <- 1000L
chr <- as.character(seqnames(gene))
reg <- GRanges(chr, IRanges(start(gene)-pad, end(gene)+pad), strand = strand(gene))

# RNA-seq coverage from bigWig
rna_track <- DataTrack(RNA_BW, type = "histogram",
                       name = "RNA-seq (coverage)",
                       baseline = 0)
# ONT per-read alignments
ont_track <- AlignmentsTrack(ONT_BAM, isPaired = FALSE,
                             name = "ONT-seq (reads)")

# Gene model from GFF
options(ucscChromosomeNames = FALSE) 
chr <- as.character(seqnames(region))
gff <- import(GFF)
tx <- gff[gff$type %in% c("mRNA") &
            as.character(seqnames(gff)) == chr]
gene_track <- GeneRegionTrack(tx, chromosome = chr,
                              name = "v68 genes", showId = TRUE,
                              transcriptAnnotation = "gene_id",
                              collapseTranscripts = "meta",
                              ucscChromosomeNames=FALSE)

pdf("PF3D7_1113100_RNA_ONT_Gviz.pdf", width = 12, height = 6)
plotTracks(list(rna_track, ont_track, gene_track),
           chromosome = chr, from = start(reg), to = end(reg),
           background.title = "white", col.axis = "black", col.title = "black")
dev.off()





library(Gviz)
library(rtracklayer)
library(GenomicRanges)
library(S4Vectors)
options(ucscChromosomeNames = FALSE)

# --- inputs ---
GFF     <- "/rhome/zli529/lab/PlasmoDB_Genome/PlasmnDB_v68/PlasmoDB-68_Pfalciparum3D7.gff"
gene_id <- "PF3D7_1113100"
pad     <- 1000L


## find the gene record
gff   <- import(GFF)
gene  <- gff[gff$type %in% c("gene","protein_coding_gene","lncRNA_gene","RNA_gene","pseudogene") &
               mcols(gff)$ID == gene_id][1]
chr   <- as.character(seqnames(gene)); gstr <- as.character(strand(gene))
gsta  <- start(gene); gend <- end(gene)

## get union CDS span for this gene (min start, max end across tx)
mrna <- gff[gff$type %in% c("mRNA","transcript","ncRNA","pseudogenic_transcript")]
tx2gene <- setNames(as.character(ifelse(!is.na(mcols(mrna)$Parent), mcols(mrna)$Parent, mcols(mrna)$gene_id)),
                    as.character(mcols(mrna)$ID))

cds <- gff[gff$type=="CDS" & as.character(seqnames(gff))==chr]
par <- mcols(cds)$Parent
if (inherits(par,"CharacterList")) par <- vapply(par, function(x) if (length(x)) x[[1]] else NA_character_, "")
cds <- cds[!is.na(par) & tx2gene[par]==gene_id]
cds_start <- min(start(cds)); cds_end <- max(end(cds))

## define contiguous 3'UTR, CDS, 5'UTR blocks along the genome
if (gstr == "+") {
  seg <- data.frame(
    start = c(cds_end+1, cds_start, gsta),
    end   = c(gend,      cds_end,   cds_start-1),
    id    = c("3'UTR","CDS","5'UTR"),
    feature = c("three_prime_UTR","CDS","five_prime_UTR")
  )
} else { # minus strand (PF3D7_1113100 is '-'): genomic L→R is 3'UTR → CDS → 5'UTR
  seg <- data.frame(
    start = c(gsta,      cds_start, cds_end+1),
    end   = c(cds_start-1, cds_end,  gend),
    id    = c("3'UTR","CDS","5'UTR"),
    feature = c("three_prime_UTR","CDS","five_prime_UTR")
  )
}
seg <- seg[seg$end >= seg$start, ]  # drop empty ends if UTR missing

## Track 1: labeled blocks for 3'UTR / CDS / 5'UTR (one line)
boxes <- AnnotationTrack(
  start = seg$start, end = seg$end, chromosome = chr, strand = TRUE,
  id = seg$id, feature = seg$feature,
  name = "PF3D7_1113100", stacking = "dense",
  shape = "box"
)
displayPars(boxes) <- list(
  featureAnnotation = "id", showFeatureId = TRUE, cex = 0.9,
  fill = c(three_prime_UTR="#7570b3", CDS="grey30", five_prime_UTR="#1b9e77"),
  thinBoxFeature = c("five_prime_UTR","three_prime_UTR"),
  thinBoxHeight  = 1,  
  col = NA,
  arrowHeadSize  = 0, 
  add53 = FALSE 
)



## Overlay both so they render as a single panel/line
ov <- OverlayTrack(trackList = list(boxes), name = "PF3D7_1113100")

## Plot (only this panel; no nearby genes)
plotTracks(
  list(ov),
  chromosome = chr, from = gsta - pad, to = gend + pad,
  # uncomment to always display 5'→3' left-to-right:
  # reverseStrand = (gstr == "-"),
  background.title = "white", col.axis = "black", col.title = "black"
)


RNA_BW  <- "~/lab/SRA_toolkit/RNAseq_7stages/processed_data/LT/RNA.bw"          # <-- edit
ONT_BAM <- "/rhome/zli529/lab/SRA_toolkit/Long_read_RNA_datasets/LeeVV_2021_nanopore/SRR11094274.no_long_introns.sorted.bam"  # <-- edit; must have .bai

# Region to plot (uses your gene ± pad)
region <- GRanges(chr, IRanges(gsta - pad, gend + pad), strand = strand(gene))
from <- start(region); to <- end(region)

# 1) RNA-seq coverage track (histogram)
bw_gr <- rtracklayer::import(RNA_BW, which = region, as = "GRanges")
ymax  <- max(1, stats::quantile(bw_gr$score, 0.99, na.rm = TRUE))  # trim spikes
rna_track <- DataTrack(
  range = bw_gr,
  type  = "histogram",
  name  = "RNA-seq (coverage)",
  baseline = 0,
  ylim  = c(0, ymax)
)

# 2) ONT per-read alignments (individual reads)
# (Make sure the BAM is indexed; if not: system(sprintf("samtools index %s", shQuote(ONT_BAM))))
ont_track <- AlignmentsTrack(
  ONT_BAM,
  isPaired = FALSE,
  name = "ONT-seq (coverage, reads)",
  stacking = "squish",         # compact rows
  showIndels = TRUE,
  add53 = FALSE                # no 5'/3' markers here
)

displayPars(ont_track) <- list(
  cex.title = 0.7,   
  fill.coverage = "#A0A0A0", # smaller title
  col.coverage    = NA ,
  fontcolor.title = "black"
)



displayPars(rna_track) <- list(
  fill.area  = "#4E79A7",
  col.histogram  = NA ,
  cex.title = 0.7,        # smaller title
  fontcolor.title = "black"
)


displayPars(ov) <- list(
  rotation.title = 180, 
  cex.title = 0.5,
  fontcolor.title = "black"                 # no label
)


# 3) Plot all three: RNA coverage, ONT reads, gene model (your 'ov')
png("PF3D7_1113100_RNA_ONT_geneModel.png", width = 2400, height = 1200, res = 300)
plotTracks(
  list(rna_track, ont_track, ov),
  chromosome = chr,
  from = from, to = to,
  sizes = c(0.6, 1.0, 0.08),
  background.title = "white",
  col.axis = "black", col.title = "black"
)
dev.off()





# ---- Inputs (edit if paths differ) ----
gene_id <- "PF3D7_1113100"
pad     <- 1000L

GFF_V68 <- "/rhome/zli529/lab/PlasmoDB_Genome/PlasmnDB_v68/PlasmoDB-68_Pfalciparum3D7.gff"
GFF_V48 <- "/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7.gff"

RNA_BW  <- "/rhome/zli529/lab/SRA_toolkit/RNAseq_7stages/processed_data/LT/RNA.bw"
ONT_BAM <- "/rhome/zli529/lab/SRA_toolkit/Long_read_RNA_datasets/LeeVV_2021_nanopore/SRR11094274.no_long_introns.sorted.bam"

stopifnot(file.exists(GFF_V68), file.exists(GFF_V48), file.exists(RNA_BW), file.exists(ONT_BAM))

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

# ---- Build RNA-seq coverage track (histogram) from BigWig ----
bw_gr <- rtracklayer::import(RNA_BW, which = region, as = "GRanges")
ymax  <- max(1, stats::quantile(bw_gr$score, 0.99, na.rm = TRUE))  # trim spikes
rna_track <- DataTrack(
  range = bw_gr, type = "histogram", name = "RNA-seq (coverage)",
  baseline = 0, ylim = c(0, ymax)
)
displayPars(rna_track) <- list(fill.area = "#4E79A7", col.histogram = NA,
                               cex.title = 0.7, fontcolor.title = "black")

# ---- Build ONT alignments track (coverage + reads) ----
# Ensure BAM is indexed: create SRR11094274.no_long_introns.sorted.bam.bai if missing
ont_track <- AlignmentsTrack(ONT_BAM, isPaired = FALSE, name = "ONT-seq (reads)",
                             stacking = "squish", showIndels = TRUE, add53 = FALSE)
displayPars(ont_track) <- list(cex.title = 0.7, fontcolor.title = "black")

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
    fill_map <- c(three_prime_UTR="#e7298a", CDS="#555555", five_prime_UTR="#66a61e", gene="#888888")
  }
  displayPars(tr) <- list(
    featureAnnotation = "id", showFeatureId = TRUE, cex = 0.8,
    fill = fill_map,
    thinBoxFeature = c("five_prime_UTR","three_prime_UTR"),
    thinBoxHeight  = 1,
    col = NA, arrowHeadSize = 0, add53 = FALSE,
    cex.title = 0.6, fontcolor.title = "black"
  )
  tr
}

# ---- Build gene model tracks ----
track_v68 <- make_gene_blocks_track(GFF_V68, gene_id, name_for_track = "v68 gene model", chr_hint = chr)
track_v48 <- make_gene_blocks_track(GFF_V48, gene_id, name_for_track = "v48 gene model", chr_hint = chr)


# rotate the track names 90 degrees and slightly shrink the title text
displayPars(track_v68); dp$rotation.title <- 90; dp$cex.title <- 0.6; displayPars(track_v68) <- dp
displayPars(track_v48); dp$rotation.title <- 90; dp$cex.title <- 0.6; displayPars(track_v48) <- dp

displayPars(ont_track) <- list(
  cex.title = 0.7,   
  fill.coverage = "#A0A0A0", # smaller title
  col.coverage    = NA ,
  fontcolor.title = "black"
)



displayPars(rna_track) <- list(
  fill.area  = "#4E79A7",
  col.histogram  = NA ,
  cex.title = 0.7,        # smaller title
  fontcolor.title = "black"
)



png(sprintf("%s_RNA_ONT_v68_v48.png", gene_id), width = 2400, height = 1400, res = 300)
plotTracks(
  list(rna_track, ont_track, track_v68, track_v48),
  chromosome = chr, from = from, to = to,
  sizes = c(0.60, 1.00, 0.12, 0.12),
  background.title = "white", col.axis = "black", col.title = "black"
  # , title.width = 1.6        # uncomment if labels are clipped
  # , reverseStrand = (gstr == "-")
)
dev.off()




# ---- Plot everything ----
png(sprintf("%s_RNA_ONT_v68_v48.png", gene_id), width = 2400, height = 1400, res = 300)
plotTracks(
  list(rna_track, ont_track, track_v68, track_v48),
  chromosome = chr, from = from, to = to,
  sizes = c(0.60, 1.00, 0.12, 0.12),
  background.title = "white", col.axis = "black", col.title = "black"
  # , reverseStrand = (gstr == "-")  # uncomment if you prefer 5'→3' left-to-right
)
dev.off()





