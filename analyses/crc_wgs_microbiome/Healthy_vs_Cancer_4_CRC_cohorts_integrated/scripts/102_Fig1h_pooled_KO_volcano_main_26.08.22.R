#!/usr/bin/env Rscript

# Main Figure 1h: pooled CRC WGS KEGG-ortholog volcano plot.
#
# Design: visual re-render only. No statistical test is recomputed here.
# Input: the saved 7,106-KO pooled differential-abundance results used for
# Supplementary Figure 3f (relative-abundance profiles; >=10% prevalence).
# Statistics retained: two-sided Wilcoxon rank-sum test, BH correction.
# Display threshold retained: q < 0.05 and |log2FC| > 0.5.
#
# Scientific caveat: the current pooled result is indexed by sequencing run.
# PRJEB6070 and PRJEB27928 contain multiple runs per biosample, so a
# biosample-level sensitivity analysis is required before final submission.

set.seed(42)

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

required_packages <- c("ragg", "digest")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

root <- here::here()
analysis_rel <- file.path(
  "공공_Metabolomics&Metagenomics_분석모음",
  "HGMT_CRC_WGS-Healthy_vs_Cancer",
  "Healthy_vs_Cancer_4_CRC_cohorts_integrated"
)
panel_source_rel <- file.path(
  "Figure_Panel_Source_Data_26.07.23", "SupFig2f.csv"
)
raw_stats_rel <- file.path(
  analysis_rel, "results_integrated", "kegg", "diff_KO_abundance_pooled.csv"
)
output_rel <- file.path(
  "공공_Metabolomics&Metagenomics_분석모음",
  "HGMT_CRC_WGS-Healthy_vs_Cancer",
  "Manuscript_Final_Panels_26.08.22"
)
script_rel <- file.path(
  analysis_rel, "scripts", "102_Fig1h_pooled_KO_volcano_main_26.08.22.R"
)

panel_source_path <- file.path(root, panel_source_rel)
raw_stats_path <- file.path(root, raw_stats_rel)
output_dir <- file.path(root, output_rel)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(panel_source_path, raw_stats_path)) {
  if (!file.exists(path)) stop("Missing required input: ", path)
}

dat <- read.csv(panel_source_path, check.names = FALSE, stringsAsFactors = FALSE)
raw_stats <- read.csv(raw_stats_path, check.names = FALSE, stringsAsFactors = FALSE)

required_cols <- c("KO", "log2FC", "p_adj", "neg_log10_p_adj", "Significance")
if (!all(required_cols %in% names(dat))) {
  stop("Panel source is missing: ", paste(setdiff(required_cols, names(dat)), collapse = ", "))
}
if (!all(c("KO", "log2FC", "p_adj") %in% names(raw_stats))) {
  stop("Raw statistics table is missing required columns.")
}
stopifnot(
  nrow(dat) == 7106L,
  nrow(raw_stats) == 7106L,
  !anyDuplicated(dat$KO),
  !anyDuplicated(raw_stats$KO),
  !anyNA(dat[, required_cols]),
  all(is.finite(dat$log2FC)),
  all(is.finite(dat$p_adj)),
  all(dat$p_adj > 0 & dat$p_adj <= 1)
)

comparison <- inner_join(
  dat %>% select(KO, log2FC, p_adj),
  raw_stats %>% select(KO, log2FC, p_adj),
  by = "KO", suffix = c("_panel", "_raw")
)
stopifnot(
  nrow(comparison) == 7106L,
  max(abs(comparison$log2FC_panel - comparison$log2FC_raw)) < 1e-12,
  max(abs(comparison$p_adj_panel - comparison$p_adj_raw)) < 1e-12
)

expected_significance <- case_when(
  dat$p_adj < 0.05 & dat$log2FC < -0.5 ~ "Enriched in Healthy",
  dat$p_adj < 0.05 & dat$log2FC >  0.5 ~ "Enriched in Cancer",
  TRUE ~ "Not Significant"
)
stopifnot(identical(dat$Significance, expected_significance))

direction_levels <- c("Healthy", "Cancer", "Not significant")
dat <- dat %>%
  mutate(
    Direction = recode(
      Significance,
      "Enriched in Healthy" = "Healthy",
      "Enriched in Cancer" = "Cancer",
      "Not Significant" = "Not significant"
    ),
    Direction = factor(Direction, levels = direction_levels)
  )

target_labels <- tibble::tribble(
  ~KO,      ~Label,                         ~NudgeX, ~NudgeY,
  "K00303", "soxB (K00303)\nDegradation",    -1.35,    2.8,
  "K08688", "Creatinase (K08688)\nProduction",  1.30,    5.0
)
targets <- dat %>%
  inner_join(target_labels, by = "KO") %>%
  arrange(KO)

stopifnot(
  identical(targets$KO, sort(target_labels$KO)),
  targets$Significance[targets$KO == "K00303"] == "Enriched in Healthy",
  targets$Significance[targets$KO == "K08688"] == "Enriched in Cancer"
)

colors <- c(
  Healthy = "#1B9E8F",
  Cancer = "#C43C3C",
  "Not significant" = "#A7A9AC"
)

p <- ggplot() +
  geom_point(
    data = filter(dat, Direction == "Not significant"),
    aes(log2FC, neg_log10_p_adj, color = Direction),
    size = 1.65, alpha = 0.42
  ) +
  geom_point(
    data = filter(dat, Direction != "Not significant"),
    aes(log2FC, neg_log10_p_adj, color = Direction),
    size = 1.8, alpha = 0.72
  ) +
  geom_vline(
    xintercept = c(-0.5, 0.5), linetype = "dashed",
    color = "grey52", linewidth = 0.55
  ) +
  geom_hline(
    yintercept = -log10(0.05), linetype = "dashed",
    color = "grey52", linewidth = 0.55
  ) +
  geom_point(
    data = targets,
    aes(log2FC, neg_log10_p_adj, fill = Direction),
    shape = 21, color = "black", stroke = 0.75, size = 3.2,
    show.legend = FALSE
  ) +
  ggrepel::geom_text_repel(
    data = targets,
    aes(log2FC, neg_log10_p_adj, label = Label),
    family = "Arial", fontface = "bold", color = "black",
    size = 4.6, lineheight = 0.92,
    bg.color = "white", bg.r = 0.10,
    box.padding = 0.65, point.padding = 0.55,
    min.segment.length = 0, max.overlaps = Inf,
    segment.color = "grey20", segment.size = 0.55,
    nudge_x = targets$NudgeX, nudge_y = targets$NudgeY,
    seed = 42, show.legend = FALSE
  ) +
  scale_color_manual(values = colors, drop = FALSE) +
  scale_fill_manual(values = colors, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.08))) +
  scale_y_continuous(expand = expansion(mult = c(0.015, 0.08))) +
  labs(
    x = expression(log[2]~fold~change~("Cancer/Healthy")),
    y = expression(-log[10]~adjusted~italic(P)),
    color = NULL
  ) +
  guides(color = guide_legend(
    nrow = 1, byrow = TRUE,
    override.aes = list(size = 3.1, alpha = 1)
  )) +
  theme_classic(base_family = "Arial", base_size = 15) +
  theme(
    axis.title = element_text(size = 17, color = "black"),
    axis.text = element_text(size = 14.5, color = "black"),
    axis.line = element_line(linewidth = 0.65, color = "black"),
    axis.ticks = element_line(linewidth = 0.55, color = "black"),
    panel.grid.major.y = element_line(color = "grey91", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.text = element_text(size = 14),
    legend.key.width = grid::unit(0.20, "in"),
    legend.spacing.x = grid::unit(0.06, "in"),
    plot.margin = margin(8, 10, 3, 8)
  )

png_path <- file.path(output_dir, "Fig1h_pooled_KO_volcano.png")
pdf_path <- file.path(output_dir, "Fig1h_pooled_KO_volcano.pdf")
source_path <- file.path(output_dir, "Fig1h_pooled_KO_volcano_source_data.csv")
stats_path <- file.path(output_dir, "Fig1h_pooled_KO_volcano_statistics.tsv")
readme_path <- file.path(output_dir, "README.txt")
session_path <- file.path(output_dir, "sessionInfo.txt")

ragg::agg_png(
  png_path, width = 7.2, height = 5.3, units = "in", res = 600,
  background = "white", scaling = 1
)
print(p)
invisible(dev.off())

grDevices::quartz(
  type = "pdf", file = pdf_path, width = 7.2, height = 5.3,
  family = "Arial"
)
print(p)
invisible(dev.off())

write.csv(
  dat %>% select(KO, log2FC, p_adj, neg_log10_p_adj, Significance),
  source_path, row.names = FALSE, quote = TRUE
)

stats_summary <- bind_rows(
  targets %>%
    transmute(
      item = KO,
      value = sprintf("log2FC=%.6f; p_adj=%.8g; %s", log2FC, p_adj, Significance)
    ),
  tibble(
    item = c("KOs plotted", "Healthy-enriched", "Cancer-enriched", "Not significant"),
    value = as.character(c(
      nrow(dat),
      sum(dat$Direction == "Healthy"),
      sum(dat$Direction == "Cancer"),
      sum(dat$Direction == "Not significant")
    ))
  )
)
write.table(
  stats_summary, stats_path, sep = "\t", row.names = FALSE, quote = FALSE
)

input_sha <- digest::digest(file = panel_source_path, algo = "sha256")
raw_sha <- digest::digest(file = raw_stats_path, algo = "sha256")
readme <- c(
  "Main Figure 1h — pooled CRC WGS KEGG-ortholog volcano plot",
  "Generated: 2026-08-22",
  "",
  "This is a visual re-render of the statistical result used for Supplementary Figure 3f.",
  "No differential-abundance test was recomputed and no value was changed.",
  "",
  paste0("Script: ", script_rel),
  paste0("Panel source: ", panel_source_rel),
  paste0("Panel-source SHA-256: ", input_sha),
  paste0("Raw statistics: ", raw_stats_rel),
  paste0("Raw-statistics SHA-256: ", raw_sha),
  "Rows: 7,106 unique KOs; no missing values.",
  "Input value space: pooled KO relative-abundance differential statistics.",
  "Filter: KOs present in >=10% of samples (upstream).",
  "Statistics: two-sided Wilcoxon rank-sum test with Benjamini-Hochberg correction (upstream).",
  "Display threshold: q < 0.05 and |log2FC| > 0.5.",
  "Labels: soxB/K00303 (degradation) and creatinase/K08688 (production).",
  "",
  "Important limitation: the current pooled WGS result is indexed by sequencing run.",
  "PRJEB6070 and PRJEB27928 contain multiple runs per biosample; a biosample-level",
  "sensitivity analysis is required before the panel is treated as submission-final."
)
writeLines(readme, readme_path, useBytes = TRUE)
capture.output(sessionInfo(), file = session_path)

outputs <- c(png_path, pdf_path, source_path, stats_path, readme_path, session_path)
stopifnot(all(file.exists(outputs)), all(file.info(outputs)$size > 0))

cat("Saved Main Figure 1h files to:\n", output_dir, "\n", sep = "")
cat("KOs plotted:", nrow(dat), "\n")
cat("Significant by display rule:", sum(dat$Direction != "Not significant"), "\n")
print(targets %>% select(KO, log2FC, p_adj, Significance))
