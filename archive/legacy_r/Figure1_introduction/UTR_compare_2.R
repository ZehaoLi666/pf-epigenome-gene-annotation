setwd("/rhome/zli529/lab/LncRNA_chip_prediction/Final/Figure1_introduction")

# Load libraries
library(dplyr)
library(ggplot2)
library(tidyr)

# Assuming df68 is already loaded
plot_utr_vs_cds <- function(df, x_max = 5000, title = "UTR length vs CDS length") {
  
  # 1. Data Preparation
  long_df <- df |>
    dplyr::filter(!is.na(cds_len), cds_len >= 0, cds_len <= x_max) |>
    dplyr::rename(UTR3 = utr3_len, UTR5 = utr5_len) |>
    tidyr::pivot_longer(cols = c(UTR3, UTR5), 
                        names_to = "UTR_type", 
                        values_to = "UTR_len") |>
    dplyr::mutate(UTR_type = factor(UTR_type, levels = c("UTR3", "UTR5"))) |>
    dplyr::filter(!is.na(UTR_type))
  
  # 2. Define Explicit Mappings
  my_shapes <- c("UTR3" = 16, "UTR5" = 17)
  my_colors <- c("UTR3" = "#2980B9", "UTR5" = "#27AE60")
  
  # 3. Create Plot with Refined (Thinner) Elements
  p_out <- ggplot(long_df, aes(x = cds_len, y = UTR_len, color = UTR_type, shape = UTR_type)) +
    geom_point(alpha = 0.4, size = 2) + # Slightly smaller points
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) + # Thinner reference line
    
    scale_x_continuous(limits = c(0, x_max)) +
    scale_y_continuous(limits = c(0, x_max)) +
    
    scale_shape_manual(values = my_shapes) + 
    scale_color_manual(values = my_colors) +
    
    labs(title = title, 
         x = "CDS length (bp)", 
         y = "UTR length (bp)") +
    
    # --- THINNER STYLE ADJUSTMENTS ---
    theme_classic(base_size = 18) + # Reduced base size from 30 to 18
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b = 15)),
      axis.title = element_text(face = "plain", size = 16), # Thinner font face
      axis.text = element_text(color = "black", size = 12),
      legend.text = element_text(size = 14),
      legend.position = "top",
      legend.title = element_blank(),
      # Reduced linewidths for a thinner appearance
      axis.line = element_line(linewidth = 0.6), # Reduced from 2 to 0.6
      axis.ticks = element_line(linewidth = 0.6), # Reduced from 2 to 0.6
      axis.ticks.length = unit(0.2, "cm")
    )
  
  return(p_out)
}

# Run the plot
p <- plot_utr_vs_cds(df68, x_max = 5000)

# Save the plot
ggsave("UTR_vs_CDS_Publication_Thin.png", 
       plot = p, 
       width = 8, # Slightly smaller width for better proportion with thinner lines
       height = 4, 
       units = "in", 
       dpi = 600)

