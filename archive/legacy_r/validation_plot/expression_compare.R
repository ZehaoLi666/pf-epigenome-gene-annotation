# -------- Inputs --------
setwd("/rhome/zli529/lab/LncRNA_chip_prediction/Final/validation_plot")


# ---- Inputs ----
novel_bed      <- "predicted_UTR3_from_CDS.bed"  # or "/mnt/data/novel_units.bed"
rna_bw         <- "/rhome/zli529/lab/LncRNA_chip_prediction/Final/window1000_bin100/coverage_pattern/new_coverage/bw/GRO.bw"
fai_path       <- "/rhome/zli529/lab/PlasmoDB_Genome/PlasmoDB_v48/PlasmoDB-48_Pfalciparum3D7_Genome.fasta.fai"
intergenic_bed <- "intergenic.bed"   # if missing, script will auto-build intergenic tiles

out_dir <- "plots"; dir.create(out_dir, showWarnings = FALSE)

stopifnot(file.exists(novel_bed), file.exists(rna_bw), file.exists(fai_path))

# ---- Load novel lncRNAs ----
gr_novel <- rtracklayer::import(novel_bed, format = "BED")

# ---- Genome sizes from .fai ----
fai <- read.table(fai_path, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
chrom_sizes <- setNames(as.numeric(fai[[2]]), fai[[1]])

# Harmonize seqlevels
keep <- intersect(names(chrom_sizes), as.character(GenomeInfoDb::seqnames(gr_novel)))
gr_novel    <- gr_novel[as.character(GenomeInfoDb::seqnames(gr_novel)) %in% keep]
chrom_sizes <- chrom_sizes[keep]
GenomeInfoDb::seqlevels(gr_novel, pruning.mode = "coarse") <- names(chrom_sizes)

# ---- Import RNA bigWig as RleList ----
bw_rle <- rtracklayer::import(rna_bw, as = "RleList")
common_seq <- intersect(names(bw_rle), GenomeInfoDb::seqlevels(gr_novel))
bw_rle   <- bw_rle[common_seq]
gr_novel <- gr_novel[as.character(GenomeInfoDb::seqnames(gr_novel)) %in% common_seq]
GenomeInfoDb::seqlevels(gr_novel, pruning.mode = "coarse") <- common_seq

# ---- Helper: mean signal over regions ----
region_means <- function(rle_list, gr) {
  chr <- as.character(GenomeInfoDb::seqnames(gr))
  idx_split <- split(seq_along(gr), chr)
  res <- numeric(length(gr))
  for (nm in names(idx_split)) {
    i <- idx_split[[nm]]
    if (!nm %in% names(rle_list)) { res[i] <- NA_real_; next }
    v <- IRanges::Views(rle_list[[nm]], IRanges::ranges(gr[i]))
    res[i] <- IRanges::viewMeans(v)
  }
  res
}

# ---- Novel expression ----
novel_expr <- region_means(bw_rle, gr_novel)

# ---- Intergenic set (use BED if present, else auto-build) ----
if (!is.na(intergenic_bed) && file.exists(intergenic_bed)) {
  gr_inter <- rtracklayer::import(intergenic_bed)
  gr_inter <- gr_inter[as.character(GenomeInfoDb::seqnames(gr_inter)) %in% common_seq]
  GenomeInfoDb::seqlevels(gr_inter, pruning.mode = "coarse") <- common_seq
  if (length(gr_inter) > length(gr_novel)) {
    set.seed(1); gr_inter <- gr_inter[sample(seq_along(gr_inter), length(gr_novel))]
  }
} else {
  med_len <- as.integer(stats::median(GenomicRanges::width(gr_novel)))
  tiles <- GenomicRanges::tileGenome(seqlengths = chrom_sizes[common_seq],
                                     tilewidth   = max(200L, med_len),
                                     cut.last.tile.in.chrom = TRUE)
  tiles <- tiles[common_seq]
  tiles <- S4Vectors::unlist(tiles, use.names = FALSE)
  tiles <- tiles[GenomicRanges::countOverlaps(tiles, gr_novel, ignore.strand = TRUE) == 0]
  tbl_novel <- table(as.character(GenomeInfoDb::seqnames(gr_novel)))
  picks <- logical(length(tiles)); picks[] <- FALSE
  set.seed(1)
  for (nm in names(tbl_novel)) {
    k <- as.integer(tbl_novel[[nm]])
    cand_idx <- which(as.character(GenomeInfoDb::seqnames(tiles)) == nm)
    if (length(cand_idx) == 0) next
    picks[sample(cand_idx, min(k, length(cand_idx)))] <- TRUE
  }
  gr_inter <- tiles[picks]
}

# ---- Intergenic expression ----
inter_expr <- region_means(bw_rle, gr_inter)

# ---- Data frame & drop NAs ----
df <- dplyr::bind_rows(
  dplyr::tibble(group = "Predicted 5'UTR", expr = novel_expr),
  dplyr::tibble(group = "Intergenic",  expr = inter_expr)
) %>% dplyr::filter(is.finite(expr))

# ---- OUTLIER REMOVAL (Tukey 1.5×IQR per group) ----
df_bounds <- df %>%
  dplyr::group_by(group) %>%
  dplyr::summarize(
    q1  = stats::quantile(expr, 0.25, na.rm = TRUE),
    q3  = stats::quantile(expr, 0.75, na.rm = TRUE),
    iqr = IQR(expr, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(lo = q1 - 1.5 * iqr,
                hi = q3 + 1.5 * iqr)

df2 <- df %>% dplyr::left_join(df_bounds, by = "group") %>%
  dplyr::mutate(is_outlier = expr < lo | expr > hi)

removed_counts <- df2 %>% dplyr::group_by(group) %>%
  dplyr::summarize(n_before = dplyr::n(),
                   n_removed = sum(is_outlier),
                   n_after = n_before - n_removed,
                   .groups = "drop")
print(removed_counts)

df_clean <- df2 %>% dplyr::filter(!is_outlier) %>% dplyr::select(group, expr)

# ---- Stats (filtered) + Wilcoxon ----
stats_clean <- df_clean %>%
  dplyr::group_by(group) %>%
  dplyr::summarize(n = dplyr::n(),
                   mean = mean(expr),
                   median = median(expr),
                   sd = sd(expr),
                   .groups = "drop")
wil_clean <- wilcox.test(expr ~ group, data = df_clean, exact = FALSE)

readr::write_tsv(stats_clean, file.path(out_dir, "novel_vs_intergenic_RNA_summary_no_outliers.tsv"))
readr::write_tsv(removed_counts, file.path(out_dir, "novel_vs_intergenic_outliers_removed.tsv"))


# ============================================================
# (2) Bar plot of mean with 95% CI (Seaborn-like)
# ============================================================
# --- Summary for bar plot (mean ± 95% CI) ---
df_bar <- df_clean %>%
  dplyr::group_by(group) %>%
  dplyr::summarize(
    n    = dplyr::n(),
    mean = mean(expr),
    sd   = sd(expr),
    se   = sd / sqrt(n),
    ci   = qt(0.975, df = n - 1) * se,  # half-width of 95% CI
    .groups = "drop"
  )

# Ensure the order is Intergenic -> Novel gene
df_bar$group <- factor(df_bar$group, levels = c("Intergenic", "Predicted 5'UTR"))

# Build the bar plot (labels use "Novel gene")
p2 <- ggplot(df_bar, aes(x = group, y = mean, fill = group)) +
  geom_col(width = 0.6, alpha = 0.9, color = NA) +
  geom_errorbar(aes(ymin = mean - ci, ymax = mean + ci),
                width = 0.15, linewidth = 0.9) +
  labs(x = NULL, y = "Mean RNA signal",
       title = "Cage-seq expression (mean ± 95% CI)",
       subtitle = sprintf("n = %d vs %d",
                          df_bar$n[df_bar$group == "Intergenic"],
                          df_bar$n[df_bar$group == "Predicted 5'UTR"])) +
  scale_fill_manual(values = c("Intergenic" = "#ff7f0e", "Predicted 5'UTR" = "#1f77b4")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        plot.title = element_text(face = "bold"))

# --- Add Wilcoxon p-value bracket ---
if (!exists("wil_clean")) {
  wil_clean <- wilcox.test(expr ~ group, data = df_clean, exact = FALSE)
}
pval  <- wil_clean$p.value
stars <- if (pval < 1e-3) "***" else if (pval < 1e-2) "**" else if (pval < 0.05) "*" else "ns"

x1 <- which(levels(df_bar$group) == "Intergenic")
x2 <- which(levels(df_bar$group) == "Predicted 5'UTR")

y_top  <- max(df_bar$mean + df_bar$ci, na.rm = TRUE)
span   <- diff(range(df_clean$expr, na.rm = TRUE)); if (span <= 0) span <- max(y_top, 1e-6)
y_bar  <- y_top + 0.06 * span
y_tick <- y_top + 0.04 * span

p2 <- p2 +
  annotate("segment", x = x1, xend = x2, y = y_bar,  yend = y_bar,  linewidth = 0.9) +
  annotate("segment", x = x1, xend = x1, y = y_tick, yend = y_bar,  linewidth = 0.9) +
  annotate("segment", x = x2, xend = x2, y = y_tick, yend = y_bar,  linewidth = 0.9) +
  annotate("text", x = (x1 + x2) / 2, y = y_bar + 0.03 * span,
           label = sprintf("p = %.3g  %s", pval, stars), size = 4) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.20)))

print(p2)
ggsave(file.path(out_dir, "Predicted 5'UTR_Cage.png"), p2,
       width = 7.2, height = 4.8, dpi = 300)
