setwd("/rhome/zli529/lab/LncRNA_chip_prediction/Final/Figure1_introduction")

suppressPackageStartupMessages({
  library(rtracklayer)
  library(GenomicRanges)
  library(dplyr)
  library(ggplot2)
})

get_tx_utr_cds <- function(gff_file) {
  gr <- rtracklayer::import(gff_file)
  
  # Pull Parent/ID as *plain character vectors* (handles CharacterList)
  parent_chr <- if ("Parent" %in% names(mcols(gr))) as.character(mcols(gr)[["Parent"]]) else rep(NA_character_, length(gr))
  id_chr     <- if ("ID"     %in% names(mcols(gr))) as.character(mcols(gr)[["ID"]])     else rep(NA_character_, length(gr))
  
  # Clean: take first parent if comma-separated
  parent_chr <- ifelse(is.na(parent_chr), NA_character_, sub(",.*$", "", parent_chr))
  
  tx_id <- parent_chr
  tx_id[is.na(tx_id)] <- id_chr[is.na(tx_id)]
  mcols(gr)$tx_id <- tx_id
  
  sum_len <- function(type_vec) {
    gr2 <- gr[mcols(gr)$type %in% type_vec]
    if (length(gr2) == 0) return(tibble::tibble(tx_id = character(), len = integer()))
    tibble::tibble(tx_id = as.character(mcols(gr2)$tx_id), w = width(gr2)) |>
      dplyr::filter(!is.na(tx_id)) |>
      dplyr::group_by(tx_id) |>
      dplyr::summarise(len = sum(w, na.rm = TRUE), .groups = "drop")
  }
  
  cds  <- sum_len("CDS") |> dplyr::rename(cds_len = len)
  
  utr5 <- sum_len(c("five_prime_UTR", "5UTR")) |> dplyr::rename(utr5_len = len)
  utr3 <- sum_len(c("three_prime_UTR", "3UTR")) |> dplyr::rename(utr3_len = len)
  
  cds |>
    dplyr::left_join(utr5, by = "tx_id") |>
    dplyr::left_join(utr3, by = "tx_id") |>
    dplyr::mutate(
      utr5_len = dplyr::coalesce(utr5_len, 0L),
      utr3_len = dplyr::coalesce(utr3_len, 0L)
    )
}

plot_utr_vs_cds <- function(df, x_max = 5000,
                            title = "UTR length vs gene length") {
  
  df_f <- df |>
    dplyr::filter(!is.na(cds_len), cds_len >= 0, cds_len <= x_max)
  
  long_df <- dplyr::bind_rows(
    df_f |> dplyr::transmute(tx_id, cds_len, UTR_len = utr3_len, UTR_type = "3'UTR"),
    df_f |> dplyr::transmute(tx_id, cds_len, UTR_len = utr5_len, UTR_type = "5'UTR")
  )
  
  ggplot(long_df, aes(x = cds_len, y = UTR_len, color = UTR_type, shape = UTR_type)) +
    geom_point(alpha = 0.35, size = 2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_x_continuous(limits = c(0, x_max)) +
    scale_y_continuous(limits = c(0, x_max)) +
    scale_shape_manual(values = c("3'UTR" = 20, "5'UTR" = 20)) +  # square vs circle
    labs(title = title, x = "CDS length (bp)", y = "UTR length (bp)") +
    theme_classic(base_size = 20) +
    theme(
      legend.position = "top",
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )
}

# Run
df68 <- get_tx_utr_cds("/rhome/zli529/lab/PlasmoDB_Genome/PlasmnDB_v68/PlasmoDB-68_Pfalciparum3D7.gff")
p <- plot_utr_vs_cds(df68, x_max = 5000)
print(p)

ggsave("UTR_vs_CDS_0_5000.png", p, width = 10, height = 6, dpi = 300)







plot_utr5_vs_cds <- function(df, x_max = 5000,
                            title = "  ") {
  df_f <- df |>
    dplyr::filter(!is.na(cds_len), cds_len >= 0, cds_len <= x_max)
  
  ggplot2::ggplot(df_f, ggplot2::aes(x = cds_len, y = utr5_len)) +
    ggplot2::geom_point(alpha = 0.35, size = 1.8) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::scale_x_continuous(limits = c(0, x_max)) +
    ggplot2::scale_y_continuous(limits = c(0, x_max)) +
    ggplot2::labs(title = title, x = "CDS length (bp)", y = "5'UTR length (bp)") +
    ggplot2::theme_classic(base_size = 20) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}

plot_utr3_vs_cds <- function(df, x_max = 5000,
                             title = " ") {
  df_f <- df |>
    dplyr::filter(!is.na(cds_len), cds_len >= 0, cds_len <= x_max)
  
  ggplot2::ggplot(df_f, ggplot2::aes(x = cds_len, y = utr3_len)) +
    ggplot2::geom_point(alpha = 0.35, size = 1.8) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggplot2::scale_x_continuous(limits = c(0, x_max)) +
    ggplot2::scale_y_continuous(limits = c(0, x_max)) +
    ggplot2::labs(title = title, x = "CDS length (bp)", y = "3'UTR length (bp)") +
    ggplot2::theme_classic(base_size = 20) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
}





p5 <- plot_utr5_vs_cds(df68, x_max = 5000)
p3 <- plot_utr3_vs_cds(df68, x_max = 5000)

print(p5)
print(p3)

ggplot2::ggsave("UTR5_vs_CDS_0_5000.png", p5, width = 5, height = 6, dpi = 600)
ggplot2::ggsave("UTR3_vs_CDS_0_5000.png", p3, width = 5, height = 6, dpi = 600)

