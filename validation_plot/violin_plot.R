# ---------- Packages ----------
library(readr)
library(dplyr)
library(ggplot2)
library(scales)

get_script_dir <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 0) {
    return(getwd())
  }
  dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/", mustWork = FALSE))
}

parse_args <- function(args) {
  parsed <- list()
  idx <- 1
  while (idx <= length(args)) {
    arg <- args[[idx]]
    if (!startsWith(arg, "--")) {
      stop(sprintf("Unexpected positional argument: %s", arg), call. = FALSE)
    }
    key <- sub("^--", "", arg)
    if (grepl("=", key, fixed = TRUE)) {
      parts <- strsplit(key, "=", fixed = TRUE)[[1]]
      parsed[[parts[1]]] <- parts[2]
    } else if (idx < length(args) && !startsWith(args[[idx + 1]], "--")) {
      parsed[[key]] <- args[[idx + 1]]
      idx <- idx + 1
    } else {
      parsed[[key]] <- TRUE
    }
    idx <- idx + 1
  }
  parsed
}

script_dir <- get_script_dir()
args <- parse_args(commandArgs(trailingOnly = TRUE))

resolve_path <- function(cli_key, env_key, default_path) {
  value <- args[[cli_key]]
  if (is.null(value) || identical(value, TRUE) || identical(value, "")) {
    env_value <- Sys.getenv(env_key, unset = "")
    value <- if (nzchar(env_value)) env_value else default_path
  }
  normalizePath(path.expand(value), winslash = "/", mustWork = FALSE)
}

is_true <- function(value, default = FALSE) {
  if (is.null(value)) {
    return(default)
  }
  if (is.logical(value)) {
    return(value)
  }
  tolower(as.character(value)) %in% c("1", "true", "yes", "y")
}

pred5_bed <- resolve_path(
  "pred5-bed",
  "PRED5_BED",
  file.path(script_dir, "predicted_UTR5_with_v68_replacements.bed")
)
pred3_bed <- resolve_path(
  "pred3-bed",
  "PRED3_BED",
  file.path(script_dir, "predicted_UTR3_from_CDS.bed")
)
v68_5_bed <- resolve_path(
  "v68-5-bed",
  "V68_5_BED",
  file.path(script_dir, "v68_gene_UTR5_longest.bed")
)
v68_3_bed <- resolve_path(
  "v68-3-bed",
  "V68_3_BED",
  file.path(script_dir, "v68_gene_UTR3_longest.bed")
)
out_png <- resolve_path(
  "out",
  "VIOLIN_OUT",
  file.path(script_dir, "UTR_length_violin_predicted_vs_v68.png")
)

use_log10_y <- !is_true(args[["linear-scale"]], default = FALSE)
min_len_bp <- 1

required_inputs <- c(pred5_bed, pred3_bed, v68_5_bed, v68_3_bed)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0) {
  stop(
    sprintf("Missing input file(s): %s", paste(missing_inputs, collapse = ", ")),
    call. = FALSE
  )
}

# ---------- Helpers ----------
read_bed6 <- function(path) {
  # BED6: chrom, start, end, name, score, strand (0-based, half-open)
  read_tsv(
    path, col_names = c("chrom","start","end","name","score","strand"),
    col_types = "ciicic", progress = FALSE
  ) %>%
    mutate(length = end - start) %>%
    filter(length >= min_len_bp)
}

# ---------- Load & label ----------
pred5 <- read_bed6(pred5_bed) %>%
  mutate(Source = "Predicted", UTR = "5′UTR")
pred3 <- read_bed6(pred3_bed) %>%
  mutate(Source = "Predicted", UTR = "3′UTR")

v685 <- read_bed6(v68_5_bed) %>%
  mutate(Source = "v68", UTR = "5′UTR")
v683 <- read_bed6(v68_3_bed) %>%
  mutate(Source = "v68", UTR = "3′UTR")

df <- bind_rows(pred5, pred3, v685, v683) %>%
  mutate(
    # nice ordering
    UTR = factor(UTR, levels = c("5′UTR","3′UTR")),
    Source = factor(Source, levels = c("Predicted","v68"))
  )

# ---------- Quick summary (prints to console) ----------
summ <- df %>%
  group_by(UTR, Source) %>%
  summarise(
    n = n(),
    median_bp = median(length),
    mean_bp = mean(length),
    p25 = quantile(length, 0.25),
    p75 = quantile(length, 0.75),
    .groups = "drop"
  )
print(summ)

# ---------- Plot ----------
p <- ggplot(df, aes(x = Source, y = length, fill = Source)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.9) +
  stat_summary(fun = median, geom = "point", size = 2.2, shape = 23, fill = "white") +
  facet_wrap(~ UTR, nrow = 1, scales = "free_x") +
  labs(
    title = NULL,
    x = NULL,
    y = "UTR length (bp)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "gray95")
  )

if (use_log10_y) {
  p <- p + scale_y_log10(labels = comma) +
    annotation_logticks(sides = "l")
} else {
  p <- p + scale_y_continuous(labels = comma)
}

print(p)

# Optionally save
ggsave(out_png, p, width = 9, height = 4, dpi = 300)
