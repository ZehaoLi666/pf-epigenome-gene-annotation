# ---------- Packages ----------

setwd("/bigdata/lerochlab/zli529/LncRNA_chip_prediction/Final/Figure4_gene")


# 1. READ DATA (Make sure these run successfully first!)
pred5 <- read_bed6(pred5_bed) %>% mutate(Source = "Predicted", UTR = "5′UTR")
pred3 <- read_bed6(pred3_bed) %>% mutate(Source = "Predicted", UTR = "3′UTR")
v685  <- read_bed6(v68_5_bed) %>% mutate(Source = "v68", UTR = "5′UTR")
v683  <- read_bed6(v68_3_bed) %>% mutate(Source = "v68", UTR = "3′UTR")

# 2. COMBINE (Renamed 'df' to 'plot_data' to avoid error)
plot_data <- bind_rows(pred5, pred3, v685, v683) %>%
  mutate(
    UTR = factor(UTR, levels = c("5′UTR","3′UTR")),
    Source = factor(Source, levels = c("Predicted","v68"))
  )

# 3. SUMMARY (Use 'plot_data' here)
summ <- plot_data %>%
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

# 4. PLOT (Use 'plot_data' here too)
p <- ggplot(plot_data, aes(x = Source, y = length, fill = Source)) +
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

# 5. SAVE PLOT
ggsave(
  filename = file.path("/bigdata/lerochlab/zli529/LncRNA_chip_prediction/Final/Figure4_gene", "Figure4E_UTR_length_comparison.png"),
  plot = p,
  width = 10,
  height = 4,
  dpi = 800
)





