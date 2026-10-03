#!/usr/bin/env Rscript

# Candidate main Figure 2g: NSCLC sarcosine pathway-score meta-analysis.
#
# Scientific scope
# - Biological comparison: ICI responders (R) versus non-responders (NR).
# - Unit of observation upstream: one patient/sample within each cohort.
# - Design: three independent NSCLC discovery cohorts; the Korean validation
#   cohort is deliberately excluded from this discovery meta-analysis.
# - Displayed statistics: the saved DerSimonian-Laird random-effects pooled
#   standardized mean differences, 95% CIs, nominal P values, and I2 values.
# - This script is display-only. It does not read raw abundance data, recompute
#   pathway scores, estimate effect sizes, fit a model, or adjust P values.
# - The historical code calls the cohort effects "Hedges' g", but its saved
#   effect-size calculation does not apply the small-sample J correction.
#   Therefore this candidate uses the accurate generic label "standardized mean
#   difference" and preserves the saved numerical values exactly.

set.seed(42)

suppressPackageStartupMessages({
  library(digest)
  library(ggplot2)
  library(here)
  library(ragg)
})

script_rel <- paste0(
  "공공_Metabolomics&Metagenomics_분석모음/",
  "HGMT_NSCLC_ICI_RvsNR_WGS/pooled_analysis/R_scripts/",
  "101_candidate_Fig2g_pathway_score_meta_26.08.22.R"
)
here::i_am(script_rel)

input_csv <- here::here(
  "공공_Metabolomics&Metagenomics_분석모음",
  "HGMT_NSCLC_ICI_RvsNR_WGS",
  "pooled_analysis",
  "results",
  "pooled_sarcosine_score_SMD_meta.csv"
)
output_dir <- here::here(
  "공공_Metabolomics&Metagenomics_분석모음",
  "HGMT_NSCLC_ICI_RvsNR_WGS",
  "pooled_analysis",
  "results",
  "Fig2g_pathway_score_meta_CANDIDATE_26.08.22"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_csv))
input_sha256 <- digest::digest(file = input_csv, algo = "sha256", serialize = FALSE)
expected_sha256 <- "7ce88475fcd722f8f65256f326a8fb293ffeb5c0278f632c2076f1d4b4594052"
stopifnot(identical(input_sha256, expected_sha256))

d <- read.csv(input_csv, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("score", "pooled_SMD", "ci_lb", "ci_ub", "p", "I2")
stopifnot(
  identical(names(d), required_columns),
  nrow(d) == 3L,
  !anyNA(d),
  identical(d$score, c("degradation", "production", "prod_deg_log2ratio"))
)

# Fixed-value guard: stop if the selected source table changes unexpectedly.
expected_values <- data.frame(
  score = c("degradation", "production", "prod_deg_log2ratio"),
  pooled_SMD = c(0.145183202020446, 0.0334695965445164, -0.104581323812219),
  ci_lb = c(0.0082018038688885, -0.189108931331644, -0.243244446901541),
  ci_ub = c(0.282164600172003, 0.256048124420677, 0.0340817992771028),
  p = c(0.0377722706208384, 0.768204831566314, 0.139346921918441),
  I2 = c(0, 53.4725853642873, 0),
  stringsAsFactors = FALSE
)
stopifnot(
  identical(d$score, expected_values$score),
  max(abs(as.matrix(d[-1]) - as.matrix(expected_values[-1]))) < 1e-12,
  all(d$ci_lb <= d$pooled_SMD),
  all(d$pooled_SMD <= d$ci_ub),
  all(d$p >= 0 & d$p <= 1),
  all(d$I2 >= 0 & d$I2 <= 100)
)

COL_R <- "#2E5F8A"
COL_NR <- "#C47B3B"
COL_TEXT <- "#202020"
COL_ZERO <- "#8A8A8A"

d$display <- c("Degradation", "Production", "Production/degradation")
d$y <- c(3, 2, 1)
d$direction <- ifelse(d$pooled_SMD >= 0, "R", "NR")
d$estimate_label <- sprintf(
  "%.2f [%.2f, %.2f]   P = %.3f   I² = %.0f%%",
  d$pooled_SMD, d$ci_lb, d$ci_ub, d$p, d$I2
)

pt_to_mm <- function(x) x / ggplot2::.pt
one_pt_mm <- 0.3528

p <- ggplot(d, aes(x = pooled_SMD, y = y, colour = direction)) +
  geom_vline(
    xintercept = 0,
    colour = COL_ZERO,
    linetype = "22",
    linewidth = one_pt_mm
  ) +
  geom_segment(
    aes(x = ci_lb, xend = ci_ub, yend = y),
    linewidth = 0.46,
    lineend = "round",
    show.legend = FALSE
  ) +
  geom_point(size = 2.15, shape = 16, show.legend = FALSE) +
  geom_text(
    aes(x = 0.315, label = estimate_label),
    hjust = 0,
    colour = COL_TEXT,
    family = "Arial",
    size = pt_to_mm(5.7),
    show.legend = FALSE
  ) +
  annotate(
    "text",
    x = -0.17,
    y = 3.47,
    label = "3 cohorts • DL",
    hjust = 0.5,
    colour = "#555555",
    family = "Arial",
    size = pt_to_mm(5.5)
  ) +
  annotate(
    "text",
    x = 0.315,
    y = 3.47,
    label = "SMD [95% CI]     P     I²",
    hjust = 0,
    colour = "#555555",
    family = "Arial",
    fontface = "bold",
    size = pt_to_mm(5.5)
  ) +
  annotate(
    "text",
    x = -0.20,
    y = 0.54,
    label = "NR higher",
    hjust = 0.5,
    colour = COL_NR,
    family = "Arial",
    size = pt_to_mm(5.5)
  ) +
  annotate(
    "text",
    x = 0.20,
    y = 0.54,
    label = "R higher",
    hjust = 0.5,
    colour = COL_R,
    family = "Arial",
    size = pt_to_mm(5.5)
  ) +
  scale_colour_manual(values = c(R = COL_R, NR = COL_NR)) +
  scale_x_continuous(
    breaks = c(-0.2, 0, 0.2),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    breaks = d$y,
    labels = d$display,
    limits = c(0.42, 3.55),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(xlim = c(-0.30, 0.30), clip = "off") +
  labs(x = "Standardized mean difference (R − NR)", y = NULL) +
  theme_classic(base_family = "Arial", base_size = 6) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    axis.title.x = element_text(size = 7, colour = "black", margin = margin(t = 3)),
    axis.text.x = element_text(size = 6, colour = "black"),
    axis.text.y = element_text(size = 6.2, colour = "black", margin = margin(r = 4)),
    axis.line.x = element_line(linewidth = one_pt_mm, colour = "black"),
    axis.line.y = element_blank(),
    axis.ticks.x = element_line(linewidth = one_pt_mm, colour = "black"),
    axis.ticks.y = element_blank(),
    axis.ticks.length = grid::unit(1.4, "pt"),
    legend.position = "none",
    plot.margin = margin(t = 2, r = 118, b = 2, l = 2, unit = "pt")
  )

png_path <- file.path(output_dir, "Fig2g_pathway_score_meta_candidate.png")
ggsave(
  filename = png_path,
  plot = p,
  width = 3.61,
  height = 1.87,
  units = "in",
  dpi = 600,
  bg = "white",
  device = ragg::agg_png
)
stopifnot(file.exists(png_path), file.info(png_path)$size > 5000)

write.csv(
  d[c(
    "score", "display", "pooled_SMD", "ci_lb", "ci_ub", "p", "I2",
    "direction", "estimate_label"
  )],
  file.path(output_dir, "Fig2g_pathway_score_meta_plotted_values.csv"),
  row.names = FALSE
)

input_info <- file.info(input_csv)
write.table(
  data.frame(
    file = normalizePath(input_csv),
    bytes = input_info$size,
    sha256 = input_sha256,
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "Fig2g_pathway_score_meta_input_checksum.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

writeLines(
  c(
    "Candidate main Figure 2g: NSCLC sarcosine pathway-score meta-analysis",
    "",
    "Status: candidate render only; not installed in the authoritative publication folder.",
    "Input: pooled_sarcosine_score_SMD_meta.csv (saved statistical results).",
    "Scope: three NSCLC discovery cohorts only; PRJEB26531 validation is excluded.",
    "Display: saved DerSimonian-Laird random-effects SMD, 95% CI, nominal P, and I2.",
    "No upstream data, pathway score, effect size, model, P value, or I2 was recomputed.",
    "The generic label standardized mean difference is used because the historical",
    "per-cohort formula does not apply the Hedges small-sample J correction.",
    "Final-size render: 3.61 x 1.87 inches at 600 dpi; place without scaling.",
    paste0("Input SHA-256: ", input_sha256)
  ),
  file.path(output_dir, "README.txt")
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

message("Created candidate: ", normalizePath(png_path))
