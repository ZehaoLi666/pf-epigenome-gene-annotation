setwd("/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content")

library(Rsamtools)
library(Biostrings)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(readr)
library(GenomeInfoDb)


# ---------- Paths ----------
REF_FA   <- "/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta"

NOVEL_BED <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content/novel_units.bed"
UTR5_refined  <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content/predicted_UTR5_from_CDS.bed"
UTR3_refined  <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content/predicted_UTR3_from_CDS.bed"
UTR5_v68 <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content/v68_gene_UTR5_longest.bed"
UTR3_v68 <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content/v68_gene_UTR3_longest.bed"
CDS_v48 <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content/CDS.bed"
intergenic <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content/intergenic.bed"


# If you don't already have a CDS BED, point to the v48 GFF and we’ll read CDS from there:
V48_GFF   <- "/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7.gff"
OUTDIR    <- "~/lab/LncRNA_chip_prediction/Final/GC_content_R"
dir.create(path.expand(OUTDIR), showWarnings = FALSE, recursive = TRUE)

# ---------- FASTA indexing ----------
faf <- FaFile(REF_FA)
if (!file.exists(paste0(REF_FA, ".fai"))) {
  message("Indexing FASTA...")
  indexFa(REF_FA)
}
open(faf)
on.exit(close(faf), add = TRUE)
fa_seqlens <- seqlengths(scanFaIndex(REF_FA))

# ---------- Helpers ----------
read_bed_gr <- function(bed_path) {
  gr <- import(con = bed_path, format = "BED")
  # enforce genome seqlevels present in FASTA
  gr <- keepSeqlevels(gr, intersect(seqlevels(gr), names(fa_seqlens)), pruning.mode = "coarse")
  # drop any zero/negative widths just in case
  gr <- gr[width(gr) > 0]
  gr
}

gc_for_gr <- function(gr, label) {
  if (length(gr) == 0) {
    warning(sprintf("No intervals for %s; returning empty tibble.", label))
    return(tibble(Category = character(), GC = numeric(), width = integer()))
  }
  # get sequences and compute GC fraction
  seqs <- getSeq(faf, gr)
  gc   <- as.numeric(letterFrequency(seqs, letters = "GC", as.prob = TRUE))
  tibble(
    Category = label,
    GC       = gc,
    width    = width(gr)
  )
}

# ---------- Load regions ----------
gr_novel <- read_bed_gr(NOVEL_BED)
gr_5utr  <- read_bed_gr(UTR5_BED)
gr_3utr  <- read_bed_gr(UTR3_BED)

# CDS from GFF (v48)
gff <- import(V48_GFF)
gr_cds <- gff[gff$type == "CDS"]
# Harmonize seqlevels with FASTA
seqlevels(gr_cds) <- intersect(seqlevels(gr_cds), names(fa_seqlens))
gr_cds <- keepSeqlevels(gr_cds, names(fa_seqlens), pruning.mode = "coarse")
gr_cds <- gr_cds[width(gr_cds) > 0]

# ---------- Compute GC ----------
df_gc <- bind_rows(
  gc_for_gr(gr_novel, "Novel units"),
  gc_for_gr(gr_5utr,  "Predicted 5′UTR"),
  gc_for_gr(gr_3utr,  "Predicted 3′UTR"),
  gc_for_gr(gr_cds,   "v48 CDS")
)

# Sanity filters (keep 0..1)
df_gc <- df_gc %>% filter(!is.na(GC), GC >= 0, GC <= 1)

# Save table
csv_path <- file.path(path.expand(OUTDIR), "Pf_GC.gc_by_region.csv")
write_csv(df_gc, csv_path)
message("Wrote: ", csv_path)

# ---------- Plots ----------
# Violin + boxplot
p_violin <- ggplot(df_gc, aes(x = Category, y = GC)) +
  geom_violin(trim = FALSE, scale = "width") +
  geom_boxplot(width = 0.15, outlier.alpha = 0.3) +
  coord_cartesian(ylim = c(0,1)) +
  labs(title = "GC content distribution by region class",
       y = "GC fraction", x = NULL) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

png_violin <- file.path(path.expand(OUTDIR), "Pf_GC.violin.png")
ggsave(png_violin, p_violin, width = 8, height = 5, dpi = 300)
message("Wrote: ", png_violin)

# Overlaid density (nice for AT-rich genomes)
p_dens <- ggplot(df_gc, aes(x = GC, color = Category)) +
  geom_density(adjust = 1, linewidth = 1) +
  coord_cartesian(xlim = c(0,1)) +
  labs(title = "GC content distributions",
       x = "GC fraction", y = "Density") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

png_dens <- file.path(path.expand(OUTDIR), "Pf_GC.density.png")
ggsave(png_dens, p_dens, width = 8, height = 5, dpi = 300)
message("Wrote: ", png_dens)

# ---------- Quick summary in console ----------
summary_tbl <- df_gc %>%
  group_by(Category) %>%
  summarize(n = n(),
            mean_GC = mean(GC),
            sd_GC   = sd(GC),
            q25 = quantile(GC, 0.25),
            median = median(GC),
            q75 = quantile(GC, 0.75)) %>%
  arrange(Category)

print(summary_tbl)














# ---------- Inputs ----------
REF_FA <- "/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta"

NOVEL_BED    <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content/novel_units.bed"
UTR5_refined <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content/predicted_UTR5_from_CDS.bed"
UTR3_refined <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content/predicted_UTR3_from_CDS.bed"
UTR5_v68     <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content/v68_gene_UTR5_longest.bed"
UTR3_v68     <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content/v68_gene_UTR3_longest.bed"
CDS_v48      <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content/CDS.fixed.bed"
intergenic   <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content/intergenic.bed"

OUTDIR <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/GC_content/plots"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

# ---------- FASTA indexing ----------
faf <- FaFile(REF_FA)
if (!file.exists(paste0(REF_FA, ".fai"))) indexFa(REF_FA)
open(faf); on.exit(close(faf), add = TRUE)
fa_idx <- scanFaIndex(REF_FA)
fa_seqlens <- seqlengths(fa_idx)
fa_seqs <- names(fa_seqlens)

# ---------- Helpers ----------
read_bed_gr <- function(bed_path) {
  gr <- import(con = bed_path, format = "BED")
  # enforce genome seqlevels present in FASTA
  gr <- keepSeqlevels(gr, intersect(seqlevels(gr), names(fa_seqlens)), pruning.mode = "coarse")
  # drop any zero/negative widths just in case
  gr <- gr[width(gr) > 0]
  gr
}

gc_for_gr <- function(gr, label) {
  if (length(gr) == 0L) return(tibble(Category = character(), GC = numeric(), width = integer()))
  seqs <- getSeq(faf, gr)
  tibble(
    Category = label,
    GC       = as.numeric(letterFrequency(seqs, letters = "GC", as.prob = TRUE)),
    width    = width(gr)
  ) |> filter(!is.na(GC), GC >= 0, GC <= 1)
}

# ---------- Load regions ----------
gr_list <- list(
  "predicted lncRNAs"          = read_bed_gr(NOVEL_BED),
  "5′UTR (refined)" = read_bed_gr(UTR5_refined),
  "3′UTR (refined)" = read_bed_gr(UTR3_refined),
  "v68 5′UTR"  = read_bed_gr(UTR5_v68),
  "v68 3′UTR"  = read_bed_gr(UTR3_v68),
  "CDS"              = read_bed_gr(CDS_v48),
  "Intergenic"           = read_bed_gr(intergenic)
)

# ---------- Compute GC ----------
df_gc <- bind_rows(lapply(names(gr_list), function(nm) gc_for_gr(gr_list[[nm]], nm)))

# Save table
csv_path <- file.path(OUTDIR, "GC_by_region.csv")
write_csv(df_gc, csv_path)
message("Wrote: ", csv_path)

# Quick summary
summary_tbl <- df_gc |>
  group_by(Category) |>
  summarize(n = n(),
            mean_GC = mean(GC),
            sd_GC = sd(GC),
            q25 = quantile(GC, 0.25),
            median = median(GC),
            q75 = quantile(GC, 0.75),
            .groups = "drop") |>
  arrange(Category)
print(summary_tbl, n = Inf)

# ---------- GC Density Plot ----------
# Order: intergenic/novel/UTRs/CDS for readability
cat_order <- c("Intergenic","predicted lncRNAs",
               "5′UTR (refined)","v68 5′UTR",
               "3′UTR (refined)","v68 3′UTR",
               "CDS")
df_gc$Category <- factor(df_gc$Category, levels = intersect(cat_order, unique(df_gc$Category)))

p_dens <- ggplot(df_gc, aes(x = GC, color = Category)) +
  geom_density(linewidth = 1) +                       # <- size instead of linewidth
  coord_cartesian(xlim = c(0, 1)) +              # <- xlim instead of xmin/xmax
  labs(title = "GC content distributions by region",
       x = "GC fraction", y = "Density") +
  theme_bw(base_size = 12) +
  theme(legend.position = "right")

png_dens <- file.path(OUTDIR, "GC_density.png")
ggsave(png_dens, p_dens, width = 8, height = 5, dpi = 300)
message("Wrote: ", png_dens)
