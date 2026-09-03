###############################################################################
# Per-cohort: regenerate dissertation figures with enlarged fonts
# (the cohort this copy operates on is set in the Config block below)
# =============================================================================
# Re-renders this cohort's three figure types at dissertation-ready sizes:
#     (1) Bacterial differential-abundance lollipop
#     (2) Sarcosine-associated-bacteria differential abundance
#     (3) Sarcosine Production-vs-Degradation pathway balance
#
# WHAT CHANGES vs the previous figures: FONT SIZES only, plus a canvas authored
#   at the dissertation text-column width (6.5 in) at 300 dpi. Long titles /
#   subtitles / captions are wrapped onto multiple lines so they fit the 6.5-in
#   width (ggplot does not auto-wrap) — wording is unchanged. Plotted values,
#   species, ordering, colours, significance thresholds and layout are NOT
#   changed.
#
# HOW each figure is produced:
#   Figures 1-2  pure re-render from the saved result CSVs (no recomputation).
#   Figure 3     no per-sample CSV was ever saved for the pathway plot, so it is
#                recomputed from KO_relative_abundance.tsv exactly as in
#                analysis_healthy_vs_cancer.R Part 5 (rowSums of degradation /
#                production KOs + Wilcoxon). The computation is deterministic;
#                the numbers are identical to the original figure.
#
# Supersedes, for this cohort's 3 figures:
#   - regenerate_diff_abundance_plots.R   (lollipop + sarcosine-bacteria)
#   - analysis_healthy_vs_cancer.R Part 5 (pathway balance)
# It does NOT modify any data file or statistical result.
#
# Run:  Rscript regenerate_plots.R   (COHORT_DIR is absolute; run from anywhere)
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(jsonlite)
})

# === Config — the ONLY lines that differ between the 4 per-cohort copies =====
COHORT_DIR <- "PRJEB6070_CRC_AdenomatousPolyps"   # relative to the analysis-folder root (working dir)
DATASET    <- "PRJEB6070_CRC_AdenomatousPolyps"   # label used verbatim in figure titles/subtitles
# =============================================================================

setwd(COHORT_DIR)

# --- PNG backend: cairo is unavailable on this machine; quartz works ---------
options(bitmapType = "quartz")
save_png <- function(plot_obj, filename, width_in, height_in, dpi = 300) {
  grDevices::png(filename = filename, width = width_in, height = height_in,
                 units = "in", res = dpi, type = "quartz", bg = "white")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot_obj)
  invisible(filename)
}

# --- Figure spec (approved): dissertation, full text-column width ------------
# The canvas is authored at the final on-page width (6.5 in) so a point size set
# here renders as that same point size in the dissertation (no down-scaling).
# 300 dpi for print sharpness.
FIG_W         <- 6.5     # canvas width (in) = dissertation text-column width
FIG_W_SARC    <- 9.0     # sarcosine-bacteria figure: wider canvas — its role
                         #   tags + bidirectional p-value labels need > 6.5 in
FIG_DPI       <- 300     # print resolution
FIG_H_PATHWAY <- 4.8     # pathway-balance canvas height (in); 3-facet boxplot
# Font sizes (points; = on-page points because canvas width = display width)
FS_TITLE   <- 22   # plot title (bold) — pathway & sarcosine-bacteria figures
FS_TITLE_LOLLI <- 22   # lollipop plot title (bold) — enlarged per request
FS_SUBT    <- 11   # subtitle
FS_CAPTION <- 8.5  # caption
FS_AXIS_T  <- 12   # axis titles
FS_AXIS_X  <- 10   # axis tick labels
FS_SPECIES <- 10   # y-axis species / category names
FS_STRIP   <- 11   # facet strip labels
FS_LEG_T   <- 10   # legend title
FS_LEG_X   <- 9    # legend text
FS_PVAL    <- 3.3  # geom_text p-value labels      (~9.4 pt on the page)
FS_ROLE    <- 3.0  # geom_text Deg/Prod role tags  (~8.5 pt on the page)

# --- Sarcosine-bacteria figure: authored larger for its wider canvas --------
# Figures 1 & 3 are authored at 6.5 in, so the point sizes above are also their
# on-page sizes. The sarcosine-bacteria figure is authored at 9 in (FIG_W_SARC)
# yet is still placed in the dissertation at the 6.5-in text-column width, so on
# the page it is scaled DOWN by 6.5/9. Authoring every absolute dimension of
# that one figure larger by SARC_SCALE = 9/6.5 cancels that down-scaling, so its
# on-page typography matches the other two figures. Wording/values unchanged.
SARC_SCALE   <- FIG_W_SARC / FIG_W       # 1.3846 ; cancels the 6.5/9 page scaling
FS_TITLE_S   <- FS_TITLE   * SARC_SCALE
FS_SUBT_S    <- FS_SUBT    * SARC_SCALE
FS_CAPTION_S <- FS_CAPTION * SARC_SCALE
FS_AXIS_T_S  <- FS_AXIS_T  * SARC_SCALE
FS_AXIS_X_S  <- FS_AXIS_X  * SARC_SCALE
FS_SPECIES_S <- FS_SPECIES * SARC_SCALE
FS_LEG_T_S   <- FS_LEG_T   * SARC_SCALE
FS_LEG_X_S   <- FS_LEG_X   * SARC_SCALE
FS_PVAL_S    <- FS_PVAL    * SARC_SCALE
FS_ROLE_S    <- FS_ROLE    * SARC_SCALE

# --- Thresholds (verbatim from the original scripts) -------------------------
P_ADJ_CUT  <- 0.05   # BH-adjusted significance cutoff
LFC_CUT    <- 1      # |log2FC| cutoff for "enriched"
N_PER_SIDE <- 15     # lollipop: top N cancer- + top N healthy-enriched

# --- Sarcosine KO panel (KO + pathway role; verbatim from the main script) ---
sarcosine_kos <- data.frame(
  KO   = c("K00301","K00302","K00303","K00304","K00305","K00306",
           "K00315","K00552","K08688"),
  Role = c("Degradation","Degradation","Degradation","Degradation",
           "Degradation","Degradation",
           "Production","Production","Production"),
  stringsAsFactors = FALSE
)

format_pval <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.0001) return("p<0.0001")
  if (p < 0.001)  return(sprintf("p=%.4f", p))
  return(sprintf("p=%.3f", p))
}

# Wrap a long string to <= `width` characters per line. ggplot does not
# auto-wrap titles/subtitles/captions, so without this they overflow (clip) the
# 6.5-in canvas. Wrapping inserts line breaks only — the wording is unchanged.
wrap_txt <- function(s, width) paste(strwrap(s, width = width), collapse = "\n")
WRAP_TITLE   <- 25   # max chars/line for the 14-pt bold title
WRAP_SUBT    <- 42   # max chars/line for the 11-pt subtitle
WRAP_CAPTION <- 56   # max chars/line for the 8.5-pt caption

cat("=== ", DATASET, " : regenerating figures (dissertation fonts) ===\n", sep = "")

# === Load metadata (Health + cancer, WGS only) ===============================
# Same parsing as analysis_healthy_vs_cancer.R: line 1 = project id, line 2 =
# header, data from line 3. list.files(...)[1] matches the original script.
cat("--- Loading metadata ---\n")
meta_file  <- list.files(pattern = "^selected_project_.*\\.txt$")[1]
meta_lines <- readLines(meta_file)
header     <- strsplit(meta_lines[2], "\t")[[1]]
data_lines <- meta_lines[3:length(meta_lines)]
data_lines <- data_lines[data_lines != ""]
data_list  <- strsplit(data_lines, "\t")
data_list  <- lapply(data_list, function(x) {
  if (length(x) >= length(header)) x[1:length(header)]
  else c(x, rep(NA, length(header) - length(x)))
})
meta <- as.data.frame(do.call(rbind, data_list), stringsAsFactors = FALSE)
colnames(meta) <- gsub(" ", ".", header)

all_pheno    <- unique(meta$Phenotype.name)
cancer_pheno <- all_pheno[!all_pheno %in% c("Health", "Adenomatous Polyps")]
if (length(cancer_pheno) == 0) stop("No cancer phenotype found in metadata.")
cancer_label <- cancer_pheno[1]

meta_hc <- meta %>%
  filter(Phenotype.name %in% c("Health", cancer_label)) %>%
  mutate(Group = ifelse(Phenotype.name == "Health", "Healthy", "Cancer"))
if ("Assay.type" %in% colnames(meta_hc)) {
  meta_hc <- meta_hc %>% filter(Assay.type == "WGS")
}
# Sample-inclusion invariants: Healthy vs CRC, WGS only.
stopifnot(!any(meta_hc$Phenotype.name == "Adenomatous Polyps"))
stopifnot(all(meta_hc$Assay.type == "WGS" | is.na(meta_hc$Assay.type)))
n_healthy <- sum(meta_hc$Group == "Healthy")
n_cancer  <- sum(meta_hc$Group == "Cancer")
cat("  Healthy = ", n_healthy, ", Cancer = ", n_cancer,
    "  (cancer phenotype: ", cancer_label, ")\n", sep = "")

# === FIGURE 1 : bacterial differential-abundance lollipop ====================
# Pure re-render of results_bacteria/diff_abundance_species.csv.
cat("--- Figure 1: lollipop (re-render from CSV) ---\n")
df <- read.csv("results_bacteria/diff_abundance_species.csv", stringsAsFactors = FALSE)

df$Significance <- "Not Significant"
df$Significance[df$p_adj < P_ADJ_CUT & df$log2FC >  LFC_CUT] <- "Enriched in Cancer"
df$Significance[df$p_adj < P_ADJ_CUT & df$log2FC < -LFC_CUT] <- "Enriched in Healthy"

cancer_pool  <- df %>% filter(Significance == "Enriched in Cancer")  %>% arrange(p_adj)
healthy_pool <- df %>% filter(Significance == "Enriched in Healthy") %>% arrange(p_adj)
top_cancer   <- head(cancer_pool,  N_PER_SIDE)
top_healthy  <- head(healthy_pool, N_PER_SIDE)
plot_df      <- bind_rows(top_cancer, top_healthy)

if (nrow(plot_df) == 0) {
  cat("  [SKIP] lollipop: no species pass p_adj<0.05 & |log2FC|>1.\n")
} else {
  plot_df <- plot_df %>% arrange(log2FC)
  plot_df$Species_display <- gsub("_", " ", plot_df$Species)
  plot_df$Species_display <- factor(plot_df$Species_display,
                                    levels = plot_df$Species_display)

  # x_max: small headroom beyond the longest lollipop so the end points sit
  # inside the panel (scale_x expand = c(0,0) leaves no automatic margin).
  x_max <- max(abs(plot_df$log2FC), na.rm = TRUE) * 1.15
  direction_colors <- c("Enriched in Cancer"  = "#B2182B",
                        "Enriched in Healthy" = "#1B7837")

  p_lolli <- ggplot(plot_df,
                    aes(x = log2FC, y = Species_display, color = Significance)) +
    geom_vline(xintercept = 0, linetype = "solid",
               color = "gray60", linewidth = 0.4) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT),
               linetype = "dashed", color = "gray70", linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = log2FC,
                     y = Species_display, yend = Species_display),
                 linewidth = 0.6) +
    geom_point(size = 3) +
    scale_color_manual(values = direction_colors, name = "Direction") +
    scale_x_continuous(limits = c(-x_max, x_max), expand = c(0, 0)) +
    labs(title = paste0("Differential Bacterial Abundance (", sub("_.*", "", DATASET), ")"),
         subtitle = wrap_txt(paste0("Top ", nrow(top_cancer), " cancer- + top ", nrow(top_healthy), " healthy-enriched"),
                             WRAP_SUBT),
         x = expression(log[2]~"Fold Change (Cancer / Healthy)"), y = NULL,
         caption = wrap_txt(paste0("Wilcoxon rank-sum test, BH-adjusted | ",
                                   "criteria: p_adj < ", P_ADJ_CUT,
                                   " & |log2FC| > ", LFC_CUT,
                                   " | dashed lines: log2FC = +/-", LFC_CUT,
                                   " | species pre-filtered at >=10% ",
                                   "prevalence in analysis_healthy_vs_cancer.R"),
                            WRAP_CAPTION)) +
    theme_minimal(base_size = 11) +
    theme(plot.title    = element_text(face = "bold", size = FS_TITLE_LOLLI, hjust = 0), plot.title.position = "plot",
          plot.subtitle = element_text(size = FS_SUBT, color = "gray40"),
          plot.caption  = element_text(size = FS_CAPTION, color = "gray50",
                                       hjust = 0),
          axis.title.x  = element_text(size = FS_AXIS_T),
          axis.text.x   = element_text(size = FS_AXIS_X),
          axis.text.y   = element_text(face = "italic", size = FS_SPECIES),
          panel.grid.major.y = element_blank(),
          panel.grid.minor   = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(size = FS_LEG_T, face = "bold"),
          legend.text  = element_text(size = FS_LEG_X),
          plot.margin  = margin(10, 18, 10, 10))

  # Height: ~0.20 in per species row + 3.5 in for the (multi-line) title,
  # axes, legend and caption; tuned for the 6.5-in canvas and 10-pt labels.
  fig_h <- max(4.5, nrow(plot_df) * 0.20 + 3.5)
  save_png(p_lolli, "results_bacteria/lollipop_diff_abundance.png",
           width_in = FIG_W + 1.5, height_in = fig_h, dpi = FIG_DPI)   # +1.5in so the title incl. cohort ID fits
  cat("  saved results_bacteria/lollipop_diff_abundance.png  (",
      nrow(plot_df), " species, ", FIG_W, " x ", round(fig_h, 2), " in )\n",
      sep = "")
}

# === FIGURE 2 : sarcosine-associated bacteria differential abundance =========
# Pure re-render of results_sarcosine/sarcosine_bacteria_diff_abundance.csv.
cat("--- Figure 2: sarcosine-bacteria DA (re-render from CSV) ---\n")
d <- read.csv("results_sarcosine/sarcosine_bacteria_diff_abundance.csv",
              stringsAsFactors = FALSE)
d$Species_display <- gsub("_", " ", d$Species)
d <- d %>% arrange(log2FC)
d$Species_display <- factor(d$Species_display, levels = d$Species_display)
d$p_text <- vapply(seq_len(nrow(d)), function(i) {
  p <- d$p_adj[i]
  if (is.na(p) || p >= 0.05) return("")
  format_pval(p)
}, character(1))
d$Role_label <- ifelse(d$Role == "Degradation", "Deg",
                ifelse(d$Role == "Production",  "Prod", "Deg & Prod"))

x_max <- max(abs(d$log2FC), na.rm = TRUE) * 1.35
role_colors <- c("Deg" = "#2196F3", "Prod" = "#FF9800", "Deg & Prod" = "#9C27B0")

p_sarc <- ggplot(d, aes(x = log2FC, y = Species_display)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray50", linewidth = 0.5 * SARC_SCALE) +
  geom_col(aes(fill = Direction), width = 0.7, alpha = 0.85) +
  geom_text(data = d[d$p_text != "", ],
            aes(x = ifelse(log2FC > 0, log2FC + x_max * 0.03,
                                        log2FC - x_max * 0.03),
                label = p_text),
            hjust = ifelse(d$log2FC[d$p_text != ""] > 0, 0, 1),
            size = FS_PVAL_S, color = "gray20", fontface = "italic") +
  # Role-tag strip: placed well right of the bars so the SARC_SCALE-enlarged
  # p-value labels on the longest bars clear it. Positions are multiples of
  # x_max (= 1.35 x max|log2FC|), so the layout is identical across cohorts.
  geom_point(aes(x = x_max * 1.53, color = Role_label),
             size = 3 * SARC_SCALE, shape = 15) +
  geom_text(aes(x = x_max * 1.60, label = Role_label, color = Role_label),
            size = FS_ROLE_S, hjust = 0, fontface = "bold") +
  scale_fill_manual(values = c("Cancer-enriched"  = "#B2182B",
                               "Healthy-enriched" = "#1B7837"),
                    name = "Direction") +
  scale_color_manual(values = role_colors, name = "Sarcosine\nPathway") +
  # x-limits widened for the SARC_SCALE-enlarged labels: left holds the
  # healthy-enriched p-value labels; right holds the longest cancer bars'
  # p-value labels plus the role-tag strip.
  scale_x_continuous(limits = c(-x_max * 1.25, x_max * 2.60), expand = c(0, 0)) +
  labs(title = wrap_txt(paste0("Sarcosine-associated gut bacteria in ", sub("_.*", "", DATASET)),
                        WRAP_TITLE),
       subtitle = wrap_txt(paste0("Healthy n=", n_healthy, " vs Cancer n=", n_cancer),
                           WRAP_SUBT),
       x = expression(log[2]~"Fold Change (Cancer / Healthy)"), y = NULL,
       caption = wrap_txt(paste0("Top species by |rho| from sarcosine KO ",
                                 "correlation analysis | Bars: log2FC of mean ",
                                 "abundance | p-values: BH-adjusted"),
                          WRAP_CAPTION)) +
  theme_minimal(base_size = 11 * SARC_SCALE) +
  theme(plot.title    = element_text(face = "bold", size = FS_TITLE_S, hjust = 0), plot.title.position = "plot",
        plot.subtitle = element_text(size = FS_SUBT_S, color = "gray40"),
        plot.caption  = element_text(size = FS_CAPTION_S, color = "gray50",
                                     hjust = 0),
        axis.title.x  = element_text(size = FS_AXIS_T_S),
        axis.text.x   = element_text(size = FS_AXIS_X_S),
        axis.text.y   = element_text(face = "italic", size = FS_SPECIES_S),
        panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        legend.position = "bottom", legend.box = "vertical",
        legend.title = element_text(size = FS_LEG_T_S, face = "bold"),
        legend.text  = element_text(size = FS_LEG_X_S),
        plot.margin  = margin(10, 15, 10, 10)) +
  guides(fill = guide_legend(order = 1), color = guide_legend(order = 2))

fig_h <- SARC_SCALE * max(4.5, nrow(d) * 0.20 + 3.9)
save_png(p_sarc, "results_sarcosine/sarcosine_bacteria_diff_abundance.png",
         width_in = FIG_W_SARC, height_in = fig_h, dpi = FIG_DPI)
cat("  saved results_sarcosine/sarcosine_bacteria_diff_abundance.png  (",
    nrow(d), " species, ", FIG_W_SARC, " x ", round(fig_h, 2), " in )\n", sep = "")

# === FIGURE 3 : sarcosine Production-vs-Degradation pathway balance ==========
# No per-sample CSV was saved for this figure, so it is recomputed from the KO
# table exactly as analysis_healthy_vs_cancer.R Part 5: rowSums of degradation /
# production KOs, log2(Prod/Deg) ratio, and Wilcoxon tests. The computation is
# deterministic — the numbers are identical to the original figure.
cat("--- Figure 3: pathway balance (recompute from KO data) ---\n")

# KO_relative_abundance.tsv is JSON despite the .tsv extension.
ko_raw  <- fromJSON("KO_relative_abundance.tsv")
ko_raw  <- ko_raw %>% filter(run_id %in% meta_hc$Run.ID)
ko_wide <- ko_raw %>%
  select(ko, run_id, abundance) %>%
  pivot_wider(names_from = ko, values_from = abundance, values_fill = 0)
ko_mat  <- as.data.frame(ko_wide)
rownames(ko_mat) <- ko_mat$run_id
ko_mat$run_id <- NULL
ko_mat  <- as.matrix(ko_mat)
ko_samples <- intersect(rownames(ko_mat), meta_hc$Run.ID)
ko_mat  <- ko_mat[ko_samples, , drop = FALSE]
meta_ko <- meta_hc %>% filter(Run.ID %in% ko_samples)
meta_ko <- meta_ko[match(ko_samples, meta_ko$Run.ID), ]
# Healthy-first group order so the pathway boxplot has Healthy on the left.
meta_ko$Group <- factor(meta_ko$Group, levels = c("Healthy", "Cancer"))
rm(ko_raw, ko_wide); gc(verbose = FALSE)
cat("  KO matrix: ", nrow(ko_mat), " samples x ", ncol(ko_mat), " KOs\n", sep = "")

deg_kos  <- intersect(sarcosine_kos$KO[sarcosine_kos$Role == "Degradation"],
                      colnames(ko_mat))
prod_kos <- intersect(sarcosine_kos$KO[sarcosine_kos$Role == "Production"],
                      colnames(ko_mat))
cat("  Degradation KOs: ", paste(deg_kos,  collapse = ","), "\n", sep = "")
cat("  Production  KOs: ", paste(prod_kos, collapse = ","), "\n", sep = "")

if (length(deg_kos) > 0 && length(prod_kos) > 0) {
  deg_sum    <- rowSums(ko_mat[, deg_kos,  drop = FALSE])
  prod_sum   <- rowSums(ko_mat[, prod_kos, drop = FALSE])
  pseudo     <- 1e-8   # pseudocount for the log-ratio (verbatim from Part 5)
  log2_ratio <- log2((prod_sum + pseudo) / (deg_sum + pseudo))

  ratio_df <- data.frame(
    Run.ID              = rownames(ko_mat),
    Group               = meta_ko$Group,
    Degradation_sum     = deg_sum,
    Production_sum      = prod_sum,
    Log2_Prod_Deg_Ratio = log2_ratio,
    stringsAsFactors    = FALSE)

  ratio_test <- tryCatch(wilcox.test(Log2_Prod_Deg_Ratio ~ Group, data = ratio_df),
                         error = function(e) list(p.value = NA))
  deg_test   <- tryCatch(wilcox.test(Degradation_sum ~ Group, data = ratio_df),
                         error = function(e) list(p.value = NA))
  prod_test  <- tryCatch(wilcox.test(Production_sum ~ Group, data = ratio_df),
                         error = function(e) list(p.value = NA))
  cat("  Wilcoxon p — Degradation: ", format_pval(deg_test$p.value),
      " | Production: ", format_pval(prod_test$p.value),
      " | log2(Prod/Deg): ", format_pval(ratio_test$p.value), "\n", sep = "")

  metric_levels <- c(paste0("Degradation\n(sum of ", length(deg_kos), " KOs)"),
                     paste0("Production\n(sum of ", length(prod_kos), " KOs)"),
                     "log2(Prod/Deg)\nRatio")

  plot_data <- ratio_df %>%
    select(Group, Degradation_sum, Production_sum, Log2_Prod_Deg_Ratio) %>%
    pivot_longer(cols = c(Degradation_sum, Production_sum, Log2_Prod_Deg_Ratio),
                 names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(Metric,
             levels = c("Degradation_sum","Production_sum","Log2_Prod_Deg_Ratio"),
             labels = metric_levels))

  pval_labels <- data.frame(
    Metric    = factor(metric_levels, levels = metric_levels),
    pval_text = c(format_pval(deg_test$p.value),
                  format_pval(prod_test$p.value),
                  format_pval(ratio_test$p.value)),
    stringsAsFactors = FALSE)
  ypos <- plot_data %>% group_by(Metric) %>%
    summarise(ymax = max(Value, na.rm = TRUE), .groups = "drop")
  pval_labels <- merge(pval_labels, ypos, by = "Metric")

  p_ratio <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 0.6, alpha = 0.3) +
    facet_wrap(~Metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
    geom_text(data = pval_labels, aes(x = 1.5, y = ymax * 1.1, label = pval_text),
              inherit.aes = FALSE, size = FS_PVAL, fontface = "italic") +
    labs(title = paste0("Sarcosine pathway balance (",
                        sub("_.*", "", DATASET), ")"),
         subtitle = paste0("Deg KOs: ", paste(deg_kos, collapse = ","),
                           "\nProd KOs: ", paste(prod_kos, collapse = ",")),
         x = "", y = "Abundance / Ratio") +
    theme_bw() +
    theme(plot.title    = element_text(size = FS_TITLE, face = "bold"), plot.title.position = "plot",
          plot.subtitle = element_text(size = FS_SUBT),
          strip.text    = element_text(size = FS_STRIP),
          axis.title    = element_text(size = FS_AXIS_T),
          axis.text     = element_text(size = FS_AXIS_X),
          legend.title  = element_text(size = FS_LEG_T),
          legend.text   = element_text(size = FS_LEG_X),
          legend.position = "bottom")

  set.seed(42)  # reproducible geom_jitter scatter (cosmetic; no data/stat effect)
  save_png(p_ratio, "results_sarcosine/sarcosine_prod_vs_deg.png",
           width_in = FIG_W, height_in = FIG_H_PATHWAY, dpi = FIG_DPI)
  cat("  saved results_sarcosine/sarcosine_prod_vs_deg.png  (",
      FIG_W, " x ", FIG_H_PATHWAY, " in )\n", sep = "")
} else {
  cat("  [SKIP] pathway: need both degradation and production KOs.\n")
}

cat("=== ", DATASET, " : done — 3 figures regenerated ===\n", sep = "")
print(sessionInfo())
