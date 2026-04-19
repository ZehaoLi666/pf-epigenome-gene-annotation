# -------- Inputs --------
setwd("/rhome/zli529/lab/LncRNA_chip_prediction/Final/Figure4_gene/")


library(ggplot2)
library(dplyr)

# ============================================================
# 1. Load Data (Based on your provided summary files)
# ============================================================

# Source 1: 5'UTR Data
df_5utr <- data.frame(
  group = c("Intergenic", "Predicted 5'UTR"),
  n = c(825, 4531),
  mean = c(6.848613205540564, 27.889952342213416),
  sd = c(13.581916153824668, 43.0689537346879)
)

# Source 2: 3'UTR Data
df_3utr <- data.frame(
  group = c("Intergenic", "Predicted 3'UTR"),
  n = c(860, 4651),
  mean = c(4.303131267440052, 58.538529752261006),
  sd = c(6.137556334189176, 78.86401690854953)
)

# ============================================================
# 2. Define Plotting Function
# ============================================================

plot_genomic_summary <- function(df_input, target_label, plot_title, y_label, output_filename) {
  
  # --- Calculate SE and CI ---
  df_bar <- df_input %>%
    mutate(
      se = sd / sqrt(n),
      ci = qt(0.975, df = n - 1) * se
    )
  
  df_bar$group <- factor(df_bar$group, levels = c("Intergenic", target_label))
  
  # --- Welch's T-test ---
  g1 <- df_bar[df_bar$group == "Intergenic", ]
  g2 <- df_bar[df_bar$group == target_label, ]
  
  t_stat <- (g1$mean - g2$mean) / sqrt((g1$sd^2 / g1$n) + (g2$sd^2 / g2$n))
  num <- ((g1$sd^2 / g1$n) + (g2$sd^2 / g2$n))^2
  den <- ((g1$sd^2 / g1$n)^2 / (g1$n - 1)) + ((g2$sd^2 / g2$n)^2 / (g2$n - 1))
  pval <- 2 * pt(-abs(t_stat), df = num / den)
  
  stars <- if (pval < 1e-3) "***" else if (pval < 1e-2) "**" else if (pval < 0.05) "*" else "ns"
  
  # --- Plotting Setup ---
  my_colors <- c("#ff7f0e", "#1f77b4")
  names(my_colors) <- c("Intergenic", target_label)
  
  y_top  <- max(df_bar$mean + df_bar$ci, na.rm = TRUE)
  y_bar  <- y_top * 1.10
  y_tick <- y_top * 1.05
  y_text <- y_top * 1.15
  
  p <- ggplot(df_bar, aes(x = group, y = mean, fill = group)) +
    geom_col(width = 0.6, alpha = 0.9, color = NA) +
    geom_errorbar(aes(ymin = pmax(0, mean - ci), ymax = mean + ci),
                  width = 0.15, linewidth = 0.9) +
    # --- UPDATED LABELS HERE ---
    labs(x = NULL, 
         y = y_label,          # Set Y-axis to TPM
         title = plot_title,   # Custom Title
         subtitle = sprintf(" ",
                            df_bar$n[df_bar$group == "Intergenic"],
                            df_bar$n[df_bar$group == target_label])) +
    scale_fill_manual(values = my_colors) +
    theme_minimal(base_size = 20) +
    theme(
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      plot.title = element_text(),
      # Large Bold X-axis labels
      axis.text.x = element_text(size = 20, color = "black"), 
      axis.title.y = element_text(size = 16, margin = margin(r = 10))
    ) +
    
    # Statistical Bracket
    annotate("segment", x = 1, xend = 2, y = y_bar,  yend = y_bar,  linewidth = 0.9) +
    annotate("segment", x = 1, xend = 1, y = y_tick, yend = y_bar,  linewidth = 0.9) +
    annotate("segment", x = 2, xend = 2, y = y_tick, yend = y_bar,  linewidth = 0.9) +
    annotate("text", x = 1.5, y = y_text,
             label = sprintf("p = %.3g %s", pval, stars), size = 5) +
    
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.25)))
  
  print(p)
  ggsave(output_filename, p, width = 7.2, height = 4.8, dpi = 300)
}

# ============================================================
# 3. Execute with Specific Titles
# ============================================================

# 1. Plot for 5'UTR (Cage-seq)
plot_genomic_summary(
  df_input = df_5utr, 
  target_label = "Predicted 5'UTR", 
  plot_title = "Cage-seq coverage",   # Specific title for 5'UTR
  y_label = "TPM",                      # Y-axis label
  output_filename = "Predicted_5UTR_Cage.png"
)

# 2. Plot for 3'UTR (PolyA-seq)
plot_genomic_summary(
  df_input = df_3utr, 
  target_label = "Predicted 3'UTR", 
  plot_title = "PolyA-seq coverage",  # Specific title for 3'UTR
  y_label = "TPM",                      # Y-axis label
  output_filename = "Predicted_3UTR_PolyA.png"
)

