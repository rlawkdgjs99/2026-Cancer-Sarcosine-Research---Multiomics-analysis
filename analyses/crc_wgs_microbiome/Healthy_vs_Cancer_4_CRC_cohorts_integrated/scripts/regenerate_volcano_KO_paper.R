###############################################################################
# Regenerate the KEGG KO differential-abundance volcano plot for paper use:
# larger fonts throughout. Faithfully reproduces Part 2 of
# integrated_analysis_pooled.R (same significance rule), but re-renders from the
# precomputed CSV instead of reloading the KO JSON files.
#
# Significance (verbatim from Part 2): p_adj < 0.05 & |log2FC| > 0.5
# Labels: top 20 KOs by p_adj (regardless of log2FC).
#
# Output: results_integrated/kegg/volcano_KO_pooled.png  (paper style)
###############################################################################

options(bitmapType = "quartz")
library(ggplot2)
library(dplyr)
library(ggrepel)

setwd("Healthy_vs_Cancer_4_CRC_cohorts_integrated")   # relative to the analysis-folder root (working dir)

# ---- CONFIG (fonts / sizes) -------------------------------------------------
FS_TITLE    <- 26
FS_SUBTITLE <- 20
FS_AXIS_T   <- 22   # axis titles
FS_AXIS_X   <- 18   # axis tick text
FS_LEGEND_T <- 20
FS_LEGEND_X <- 18
LABEL_SIZE  <- 6.5  # enzyme labels (only significant sarcosine KOs)
POINT_SIZE  <- 1.6
N_SAMPLES   <- 1647 # = nrow(ko_mat) in Part 2 (pooled KO-matrix sample count)

OUT <- "results_integrated/kegg/volcano_KO_pooled.png"

# ---- Data (precomputed) -----------------------------------------------------
ko_diff <- read.csv("results_integrated/kegg/diff_KO_abundance_pooled.csv",
                    stringsAsFactors = FALSE)

# Significance classification — identical to Part 2
ko_diff$Significance <- "Not Significant"
ko_diff$Significance[ko_diff$p_adj < 0.05 & ko_diff$log2FC > 0.5]  <- "Enriched in Cancer"
ko_diff$Significance[ko_diff$p_adj < 0.05 & ko_diff$log2FC < -0.5] <- "Enriched in Healthy"
ko_diff <- ko_diff %>% arrange(p_adj)

# Label ONLY significant sarcosine-related KOs (Part 3 panel, p_adj < 0.05).
# Note: only KOs passing prevalence >=10% appear on the volcano at all, so
# low-prevalence sarcosine KOs (soxD/soxG/PIPOX/DMGDH/GNMT/SOX-mono) are absent.
sarc_panel <- data.frame(
  KO    = c("K00301","K00302","K00303","K00304","K00305","K00306","K00315","K00552","K08688"),
  Short = c("SOX(mono)","soxA","soxB","soxD","soxG","PIPOX","DMGDH","GNMT","Creatinase"),
  RoleF = c("Degradation","Degradation","Degradation","Degradation","Degradation",
            "Degradation","Production","Production","Production"),
  stringsAsFactors = FALSE)

sarc_lab <- ko_diff %>%
  dplyr::inner_join(sarc_panel, by = "KO") %>%
  filter(p_adj < 0.05) %>%
  mutate(fillcol = ifelse(log2FC > 0, "#B2182B", "#1B7837"),
         lab = paste0(KO, " — ", Short, "\n(", RoleF, ")"),
         nx  = ifelse(log2FC > 0, 1.0, -1.0))

cat("KOs plotted:", nrow(ko_diff),
    "| Cancer:", sum(ko_diff$Significance == "Enriched in Cancer"),
    "| Healthy:", sum(ko_diff$Significance == "Enriched in Healthy"), "\n")
cat("Sarcosine KOs labeled:", paste0(sarc_lab$KO, "/", sarc_lab$Short, collapse = ", "), "\n")

# ---- Plot -------------------------------------------------------------------
p <- ggplot(ko_diff, aes(x = log2FC, y = -log10(p_adj))) +
  geom_point(aes(color = Significance), alpha = 0.6, size = POINT_SIZE) +
  scale_color_manual(values = c("Enriched in Cancer" = "#B2182B",
                                "Enriched in Healthy" = "#1B7837",
                                "Not Significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  # label only the significant sarcosine KOs — segment + text, no marker
  ggrepel::geom_text_repel(data = sarc_lab, aes(label = lab),
             color = "black", fontface = "bold", size = LABEL_SIZE,
             bg.color = "white", bg.r = 0.12, lineheight = 0.9,
             box.padding = 1.0, point.padding = 0.6,
             segment.color = "grey20", segment.size = 0.6,
             min.segment.length = 0, max.overlaps = Inf,
             nudge_x = sarc_lab$nx, nudge_y = 4, show.legend = FALSE) +
  labs(title = paste0("Integrated 4 CRC Cohorts: Differential KO Abundance (n=", N_SAMPLES, ")"),
       subtitle = "Healthy vs Colorectal Neoplasms (Pooled) — significant sarcosine KOs highlighted",
       x = "log2FC (Cancer/Healthy)", y = "-log10(adj. p-value)") +
  theme_bw() +
  theme(plot.title    = element_text(size = FS_TITLE, face = "bold"),
        plot.subtitle = element_text(size = FS_SUBTITLE),
        axis.title    = element_text(size = FS_AXIS_T),
        axis.text     = element_text(size = FS_AXIS_X),
        legend.title  = element_text(size = FS_LEGEND_T),
        legend.text   = element_text(size = FS_LEGEND_X),
        legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

ggsave(OUT, p, width = 12, height = 8.5, dpi = 300, bg = "white",
       device = grDevices::png, type = "quartz")
cat("Saved:", OUT, "\n")
