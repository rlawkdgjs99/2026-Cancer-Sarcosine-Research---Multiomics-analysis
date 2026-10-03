###############################################################################
# PRJEB10878_CRC: Lollipop plot of differential bacterial abundance
#
# Input  : results_bacteria/diff_abundance_species.csv
#          (already prevalence-filtered at >=10% in analysis_healthy_vs_cancer.R)
# Output : results_bacteria/lollipop_diff_abundance.png
#
# Selection : top 15 cancer-enriched + top 15 healthy-enriched (by adj. p-value)
#             using same significance criteria as the volcano plot:
#             p_adj < 0.05 AND |log2FC| > 1.
###############################################################################

library(ggplot2)
library(dplyr)

# PNG device wrapper: cairo unavailable on this machine (libSM.6.dylib missing),
# so force quartz, which is supported (capabilities()['aqua'] == TRUE).
.ggsave <- function(filename, plot, width, height, dpi, ...) {
  grDevices::png(filename = filename, width = width, height = height,
                 units = "in", res = dpi, type = "quartz", bg = "white")
  print(plot)
  grDevices::dev.off()
  invisible(filename)
}

project_dir <- "PRJEB10878_CRC"   # relative to the analysis-folder root (working dir)
setwd(project_dir)

DATASET      <- "PRJEB10878_CRC"
CANCER_LABEL <- "CRC"
N_PER_SIDE   <- 15
P_ADJ_CUT    <- 0.05
LFC_CUT      <- 1

format_pval <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.0001) return("p<0.0001")
  if (p < 0.001)  return(sprintf("p=%.4f", p))
  return(sprintf("p=%.3f", p))
}

cat("=== ", DATASET, ": Lollipop differential abundance ===\n", sep = "")

df <- read.csv("results_bacteria/diff_abundance_species.csv", stringsAsFactors = FALSE)
cat("Total species in CSV (after upstream >=10% prevalence filter): ", nrow(df), "\n", sep = "")

df$Significance <- "Not Significant"
df$Significance[df$p_adj < P_ADJ_CUT & df$log2FC >  LFC_CUT] <- "Enriched in Cancer"
df$Significance[df$p_adj < P_ADJ_CUT & df$log2FC < -LFC_CUT] <- "Enriched in Healthy"

cancer_pool  <- df %>% filter(Significance == "Enriched in Cancer")  %>% arrange(p_adj)
healthy_pool <- df %>% filter(Significance == "Enriched in Healthy") %>% arrange(p_adj)
cat("Significant cancer-enriched : ", nrow(cancer_pool),  "\n", sep = "")
cat("Significant healthy-enriched: ", nrow(healthy_pool), "\n", sep = "")

top_cancer  <- head(cancer_pool,  N_PER_SIDE)
top_healthy <- head(healthy_pool, N_PER_SIDE)
plot_df <- bind_rows(top_cancer, top_healthy)

if (nrow(plot_df) == 0) {
  stop("No species pass p_adj<0.05 & |log2FC|>1; nothing to plot.")
}

cat("Plotting: ", nrow(top_cancer), " cancer-enriched + ",
    nrow(top_healthy), " healthy-enriched = ", nrow(plot_df), " species\n", sep = "")

plot_df <- plot_df %>% arrange(log2FC)
plot_df$Species_display <- gsub("_", " ", plot_df$Species)
plot_df$Species_display <- factor(plot_df$Species_display, levels = plot_df$Species_display)
plot_df$pval_text <- vapply(plot_df$p_adj, format_pval, character(1))

x_max <- max(abs(plot_df$log2FC), na.rm = TRUE) * 1.30
direction_colors <- c("Enriched in Cancer"  = "#B2182B",
                      "Enriched in Healthy" = "#1B7837")

p_lolli <- ggplot(plot_df, aes(x = log2FC, y = Species_display, color = Significance)) +
  geom_vline(xintercept = 0, linetype = "solid",  color = "gray60", linewidth = 0.4) +
  geom_vline(xintercept = c(-LFC_CUT, LFC_CUT),
             linetype = "dashed", color = "gray70", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = log2FC,
                   y = Species_display, yend = Species_display),
               linewidth = 0.6) +
  geom_point(size = 3) +
  geom_text(aes(x = ifelse(log2FC > 0, log2FC + x_max * 0.03, log2FC - x_max * 0.03),
                label = pval_text),
            hjust = ifelse(plot_df$log2FC > 0, 0, 1),
            size = 2.7, color = "gray25", fontface = "italic",
            show.legend = FALSE) +
  scale_color_manual(values = direction_colors, name = "Direction") +
  scale_x_continuous(limits = c(-x_max, x_max), expand = c(0, 0)) +
  labs(title = paste0(DATASET, ": Differential Bacterial Abundance"),
       subtitle = paste0("Healthy vs ", CANCER_LABEL,
                         " | top ", nrow(top_cancer), " cancer-enriched + top ",
                         nrow(top_healthy), " healthy-enriched (by adj. p-value)"),
       x = expression(log[2]~"Fold Change (Cancer / Healthy)"),
       y = NULL,
       caption = paste0("Wilcoxon rank-sum test, BH-adjusted | criteria: p_adj < ",
                        P_ADJ_CUT, " & |log2FC| > ", LFC_CUT,
                        " | dashed lines: log2FC = +/-", LFC_CUT,
                        "\nSpecies pre-filtered at >=10% prevalence in analysis_healthy_vs_cancer.R")) +
  theme_minimal(base_size = 11) +
  theme(plot.title    = element_text(face = "bold", size = 13, hjust = 0),
        plot.subtitle = element_text(size = 10, color = "gray40"),
        plot.caption  = element_text(size  = 8,  color = "gray50", hjust = 0),
        axis.text.y   = element_text(face = "italic", size = 9),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 9, face = "bold"),
        plot.margin = margin(10, 18, 10, 10))

fig_h <- max(5, nrow(plot_df) * 0.32 + 2)
out_path <- "results_bacteria/lollipop_diff_abundance.png"
.ggsave(out_path, p_lolli, width = 10, height = fig_h, dpi = 200)
cat("Saved: ", file.path(project_dir, out_path), "\n", sep = "")
