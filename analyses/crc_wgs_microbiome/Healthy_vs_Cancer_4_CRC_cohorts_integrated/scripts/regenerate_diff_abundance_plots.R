###############################################################################
# [SUPERSEDED 2026-05-21]  Do NOT run this script.
#
# Fully superseded by the per-folder regeneration scripts:
#   PRJEB10878_CRC/scripts/regenerate_plots.R
#   PRJEB27928_CRC/scripts/regenerate_plots.R
#   PRJEB6070_CRC_AdenomatousPolyps/scripts/regenerate_plots.R
#   PRJNA429097_CRC/scripts/regenerate_plots.R
#   Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/regenerate_plots.R
# which re-render the lollipop and sarcosine-bacteria figures at dissertation
# font sizes (6.5-in canvas, 300 dpi). Running THIS script would overwrite
# those figures with the older 14-in / smaller-font versions. Kept for
# reference only.
###############################################################################
###############################################################################
# Regenerate Differential-Abundance Plots with Larger Fonts
#
# Purpose: re-render two figure types at larger, more readable font sizes:
#   (1) Bacterial differential-abundance lollipop plots
#   (2) Sarcosine-associated-bacteria differential-abundance plots
# for the pooled analysis AND each of the 4 individual cohorts.
#
# IMPORTANT - this is a RE-RENDER ONLY:
#   - It reads the already-saved differential-abundance CSVs and re-plots them.
#     No data is reloaded and NO statistic is recomputed. The plotted values,
#     species, ordering, significance and selection are bit-for-bit identical to
#     the original figures. The changes are font sizes (titles, subtitles,
#     captions, axis text incl. species names, legends, in-plot p-value/role
#     labels) and a proportionally wider canvas so the larger text does not clip.
#   - Original analysis scripts are NOT modified. Output PNGs overwrite the
#     existing figures at the same paths (same pattern as
#     regenerate_diversity_plots.R).
#   - Font sizes follow regenerate_diversity_plots.R so the figure set stays
#     visually consistent.
#
# Inputs  (per cohort + pooled): results_*/diff_abundance_species*.csv,
#         results_*/sarcosine_bacteria_diff_abundance*.csv,
#         results_integrated/supplementary/cohort_sample_summary.csv
# Outputs (overwritten in place): 5 lollipop PNGs + 5 sarcosine-bacteria PNGs
###############################################################################

library(ggplot2)
library(dplyr)

# quartz PNG backend (cairo unavailable on this machine; same as the other scripts)
options(bitmapType = "quartz")

.ggsave <- function(filename, plot, width, height, dpi = 200) {
  grDevices::png(filename = filename, width = width, height = height,
                 units = "in", res = dpi, type = "quartz", bg = "white")
  print(plot)
  grDevices::dev.off()
  invisible(filename)
}

# --- Current project root ----------------------------------------------------
ROOT       <- "."   # analysis-folder root (working dir)
setwd(ROOT)
INTEGRATED <- "Healthy_vs_Cancer_4_CRC_cohorts_integrated"

# --- Larger font sizes (consistent with regenerate_diversity_plots.R) --------
# Original sizes are shown in comments for reference.
FS_TITLE   <- 26    # plot title           (was 13; raised to 20 per user request)
FS_SUBT    <- 16    # subtitle             (was 10)
FS_CAPTION <- 13    # caption              (was 8)
FS_AXIS_T  <- 17    # axis title           (was ~11, theme_minimal default)
FS_AXIS_X  <- 16    # x-axis tick text     (was ~9,  theme_minimal default)
FS_SPECIES <- 16    # y-axis tick text = species names (was 9)
FS_LEG_T   <- 16    # legend title         (was 9)
FS_LEG_X   <- 14    # legend text          (was ~9, theme_minimal default)
FS_PVAL    <- 5.2   # in-plot p-value labels (was 2.7 lollipop / 2.8 sarcosine)
FS_ROLE    <- 4.6   # in-plot Deg/Prod role labels (was 2.5)

# --- Shared constants (verbatim from the original scripts) -------------------
P_ADJ_CUT  <- 0.05
LFC_CUT    <- 1
N_PER_SIDE <- 15

format_pval <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.0001) return("p<0.0001")
  if (p < 0.001)  return(sprintf("p=%.4f", p))
  return(sprintf("p=%.3f", p))
}

###############################################################################
# FIGURE 1: Bacterial differential-abundance lollipop
#   Re-render of lollipop_diff_abundance(_pooled).R - verbatim plot logic,
#   larger fonts only.
###############################################################################
make_lollipop <- function(csv_path, out_path, dataset_label, cancer_label,
                           prefilter_note) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)

  df$Significance <- "Not Significant"
  df$Significance[df$p_adj < P_ADJ_CUT & df$log2FC >  LFC_CUT] <- "Enriched in Cancer"
  df$Significance[df$p_adj < P_ADJ_CUT & df$log2FC < -LFC_CUT] <- "Enriched in Healthy"

  cancer_pool  <- df %>% filter(Significance == "Enriched in Cancer")  %>% arrange(p_adj)
  healthy_pool <- df %>% filter(Significance == "Enriched in Healthy") %>% arrange(p_adj)

  top_cancer  <- head(cancer_pool,  N_PER_SIDE)
  top_healthy <- head(healthy_pool, N_PER_SIDE)
  plot_df <- bind_rows(top_cancer, top_healthy)

  if (nrow(plot_df) == 0) {
    cat("  [SKIP]", out_path, "- no species pass thresholds.\n")
    return(invisible(NULL))
  }

  plot_df <- plot_df %>% arrange(log2FC)
  plot_df$Species_display <- gsub("_", " ", plot_df$Species)
  plot_df$Species_display <- factor(plot_df$Species_display, levels = plot_df$Species_display)
  plot_df$pval_text <- vapply(plot_df$p_adj, format_pval, character(1))

  x_max <- max(abs(plot_df$log2FC), na.rm = TRUE) * 1.40
  direction_colors <- c("Enriched in Cancer"  = "#B2182B",
                        "Enriched in Healthy" = "#1B7837")

  p_lolli <- ggplot(plot_df, aes(x = log2FC, y = Species_display, color = Significance)) +
    geom_vline(xintercept = 0, linetype = "solid", color = "gray60", linewidth = 0.4) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT),
               linetype = "dashed", color = "gray70", linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = log2FC,
                     y = Species_display, yend = Species_display),
                 linewidth = 0.6) +
    geom_point(size = 3) +
    geom_text(aes(x = ifelse(log2FC > 0, log2FC + x_max * 0.03, log2FC - x_max * 0.03),
                  label = pval_text),
              hjust = ifelse(plot_df$log2FC > 0, 0, 1),
              size = FS_PVAL, color = "gray25", fontface = "italic",
              show.legend = FALSE) +
    scale_color_manual(values = direction_colors, name = "Direction") +
    scale_x_continuous(limits = c(-x_max, x_max), expand = c(0, 0)) +
    labs(title = paste0(dataset_label, ": Differential Bacterial Abundance"),
         subtitle = paste0("Healthy vs ", cancer_label,
                           " | top ", nrow(top_cancer), " cancer-enriched + top ",
                           nrow(top_healthy), " healthy-enriched (by adj. p-value)"),
         x = expression(log[2]~"Fold Change (Cancer / Healthy)"),
         y = NULL,
         caption = paste0("Wilcoxon rank-sum test, BH-adjusted | criteria: p_adj < ",
                          P_ADJ_CUT, " & |log2FC| > ", LFC_CUT,
                          " | dashed lines: log2FC = +/-", LFC_CUT,
                          "\nSpecies pre-filtered at >=10% prevalence in ", prefilter_note)) +
    theme_minimal(base_size = 16) +
    theme(plot.title    = element_text(face = "bold", size = FS_TITLE, hjust = 0),
          plot.subtitle = element_text(size = FS_SUBT, color = "gray40"),
          plot.caption  = element_text(size = FS_CAPTION, color = "gray50", hjust = 0),
          axis.title.x  = element_text(size = FS_AXIS_T),
          axis.text.x   = element_text(size = FS_AXIS_X),
          axis.text.y   = element_text(face = "italic", size = FS_SPECIES),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(size = FS_LEG_T, face = "bold"),
          legend.text  = element_text(size = FS_LEG_X),
          plot.margin = margin(10, 18, 10, 10))

  fig_h <- max(6, nrow(plot_df) * 0.34 + 2.5)
  .ggsave(out_path, p_lolli, width = 14, height = fig_h, dpi = 200)
  cat("  [SAVED]", out_path, "(", nrow(plot_df), "species )\n")
}

###############################################################################
# FIGURE 2: Sarcosine-associated bacteria differential abundance
#   Re-render of Part 6 of analysis_healthy_vs_cancer.R /
#   integrated_analysis_pooled.R - verbatim plot logic, larger fonts only.
###############################################################################
make_sarc_bact <- function(csv_path, out_path, title_txt, subtitle_txt, caption_txt) {
  d <- read.csv(csv_path, stringsAsFactors = FALSE)

  d$Species_display <- gsub("_", " ", d$Species)
  d <- d %>% arrange(log2FC)
  d$Species_display <- factor(d$Species_display, levels = d$Species_display)

  d$p_text <- sapply(seq_len(nrow(d)), function(i) {
    p <- d$p_adj[i]
    if (is.na(p) || p >= 0.05) return("")
    format_pval(p)
  })
  d$Role_label <- ifelse(d$Role == "Degradation", "Deg",
                  ifelse(d$Role == "Production", "Prod", "Deg & Prod"))

  x_max <- max(abs(d$log2FC), na.rm = TRUE) * 1.35
  role_colors <- c("Deg" = "#2196F3", "Prod" = "#FF9800", "Deg & Prod" = "#9C27B0")

  p <- ggplot(d, aes(x = log2FC, y = Species_display)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    geom_col(aes(fill = Direction), width = 0.7, alpha = 0.85) +
    geom_text(data = d[d$p_text != "", ],
      aes(x = ifelse(log2FC > 0, log2FC + x_max * 0.03, log2FC - x_max * 0.03),
          label = p_text),
      hjust = ifelse(d$log2FC[d$p_text != ""] > 0, 0, 1),
      size = FS_PVAL, color = "gray20", fontface = "italic") +
    geom_point(aes(x = x_max * 0.92, color = Role_label), size = 3, shape = 15) +
    geom_text(aes(x = x_max * 0.97, label = Role_label, color = Role_label),
              size = FS_ROLE, hjust = 0, fontface = "bold") +
    scale_fill_manual(values = c("Cancer-enriched"  = "#B2182B",
                                 "Healthy-enriched" = "#1B7837"), name = "Direction") +
    scale_color_manual(values = role_colors, name = "Sarcosine\nPathway") +
    scale_x_continuous(limits = c(-x_max, x_max * 1.15), expand = c(0, 0)) +
    labs(title = title_txt, subtitle = subtitle_txt,
         x = expression(log[2]~"Fold Change (Cancer / Healthy)"), y = NULL,
         caption = caption_txt) +
    theme_minimal(base_size = 16) +
    theme(plot.title    = element_text(face = "bold", size = FS_TITLE, hjust = 0),
          plot.subtitle = element_text(size = FS_SUBT, color = "gray40"),
          plot.caption  = element_text(size = FS_CAPTION, color = "gray50", hjust = 0),
          axis.title.x  = element_text(size = FS_AXIS_T),
          axis.text.x   = element_text(size = FS_AXIS_X),
          axis.text.y   = element_text(face = "italic", size = FS_SPECIES),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          legend.position = "bottom", legend.box = "horizontal",
          legend.title = element_text(size = FS_LEG_T, face = "bold"),
          legend.text  = element_text(size = FS_LEG_X),
          plot.margin = margin(10, 15, 10, 10)) +
    guides(fill = guide_legend(order = 1), color = guide_legend(order = 2))

  fig_h <- max(7, nrow(d) * 0.38 + 2.5)
  .ggsave(out_path, p, width = 14, height = fig_h, dpi = 200)
  cat("  [SAVED]", out_path, "(", nrow(d), "species )\n")
}

###############################################################################
# Cohort configuration (labels verbatim from the original per-cohort scripts)
###############################################################################
cohorts <- list(
  list(key = "PRJEB10878",  dir = "PRJEB10878_CRC",                  cancer = "CRC"),
  list(key = "PRJEB27928",  dir = "PRJEB27928_CRC",                  cancer = "CRC"),
  list(key = "PRJEB6070",   dir = "PRJEB6070_CRC_AdenomatousPolyps", cancer = "CRC"),
  list(key = "PRJNA429097", dir = "PRJNA429097_CRC",                 cancer = "CRC")
)

# Sample counts for the sarcosine-bacteria subtitles (metadata counts, exactly
# as the original scripts compute them).
samp <- read.csv(file.path(INTEGRATED,
                 "results_integrated/supplementary/cohort_sample_summary.csv"),
                 stringsAsFactors = FALSE)

cat("=== Regenerating differential-abundance plots with larger fonts ===\n\n")

# --- Pooled ------------------------------------------------------------------
cat("--- Pooled (4 CRC cohorts integrated) ---\n")
pooled_h <- sum(samp$Healthy)
pooled_c <- sum(samp$Cancer)

make_lollipop(
  csv_path = file.path(INTEGRATED, "results_integrated/bacteria/diff_abundance_species_pooled.csv"),
  out_path = file.path(INTEGRATED, "results_integrated/bacteria/lollipop_diff_abundance_pooled.png"),
  dataset_label  = "Pooled (4 CRC cohorts integrated)",
  cancer_label   = "Colorectal Neoplasms (Pooled)",
  prefilter_note = "integrated_analysis_pooled.R")

make_sarc_bact(
  csv_path = file.path(INTEGRATED, "results_integrated/sarcosine/sarcosine_bacteria_diff_abundance_pooled.csv"),
  out_path = file.path(INTEGRATED, "results_integrated/sarcosine/sarcosine_bacteria_diff_abundance_pooled.png"),
  title_txt    = "Sarcosine-Associated Bacteria: Differential Abundance (Pooled)",
  subtitle_txt = paste0("Integrated 4 CRC Cohorts | Healthy (n=", pooled_h,
                        ") vs Cancer (n=", pooled_c, ")"),
  caption_txt  = "Top species by |rho| from sarcosine KO correlation\nBars: log2FC | p-values: BH-adjusted")

# --- Per cohort --------------------------------------------------------------
for (co in cohorts) {
  cat("--- ", co$dir, " ---\n", sep = "")
  row <- samp[samp$Cohort == co$key, ]
  n_h <- row$Healthy
  n_c <- row$Cancer

  make_lollipop(
    csv_path = file.path(co$dir, "results_bacteria/diff_abundance_species.csv"),
    out_path = file.path(co$dir, "results_bacteria/lollipop_diff_abundance.png"),
    dataset_label  = co$dir,
    cancer_label   = co$cancer,
    prefilter_note = "analysis_healthy_vs_cancer.R")

  make_sarc_bact(
    csv_path = file.path(co$dir, "results_sarcosine/sarcosine_bacteria_diff_abundance.csv"),
    out_path = file.path(co$dir, "results_sarcosine/sarcosine_bacteria_diff_abundance.png"),
    title_txt    = "Sarcosine-Associated Bacteria: Differential Abundance",
    subtitle_txt = paste0(co$dir, " | Healthy (n=", n_h, ") vs Cancer (n=", n_c, ")"),
    caption_txt  = "Top species by |rho| from sarcosine KO correlation analysis\nBars: log2FC of mean abundance | p-values: BH-adjusted")
}

cat("\n=== Done. 10 figures re-rendered with larger fonts. ===\n")
