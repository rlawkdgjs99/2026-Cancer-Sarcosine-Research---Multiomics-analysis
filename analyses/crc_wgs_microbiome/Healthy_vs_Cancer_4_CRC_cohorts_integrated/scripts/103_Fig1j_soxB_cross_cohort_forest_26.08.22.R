#!/usr/bin/env Rscript

# Proposed Figure 1j: cross-cohort reproducibility of soxB (K00303) depletion.
#
# Biological question
#   Is the lower pooled abundance of soxB in CRC reproduced independently in
#   each of the four CRC WGS cohorts?
#
# Statistical design
#   - Unit of observation: one WGS stool sample.
#   - Design: two independent groups (Healthy versus CRC), stratified by cohort.
#   - Input: HGMT-exported KO relative abundance; values are continuous,
#     non-negative, and may contain genuine zeros.
#   - Effect: r_rb = P(Healthy > CRC) - P(Healthy < CRC); ties contribute zero.
#     Positive values therefore mean that soxB is higher in Healthy controls.
#   - Uncertainty: group-stratified percentile bootstrap, 5,000 repetitions,
#     seed 42. Wilcoxon tests are reported as supporting statistics; because
#     this is one pre-specified KO, no genome-wide feature selection is done.
#   - This is a cohort-stratified effect-size forest plot. It is not labelled a
#     formal meta-analysis because no pooled meta-analytic estimate is fitted.
#
# Workflow
#   1. Inspect and validate the exported sample-level table.
#   2. Calculate r_rb independently within each cohort.
#   3. Cross-check r_rb by direct pairwise comparison and Wilcoxon W.
#   4. Bootstrap within Healthy and CRC strata, then export data and artwork.

suppressPackageStartupMessages({
  library(digest)
  library(ggplot2)
  library(here)
  library(ragg)
  library(svglite)
})

set.seed(
  42L,
  kind = "Mersenne-Twister",
  normal.kind = "Inversion",
  sample.kind = "Rejection"
)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path.")
SCRIPT_PATH <- normalizePath(sub("^--file=", "", script_arg))
BASE_DIR <- normalizePath(file.path(dirname(SCRIPT_PATH), ".."))
setwd(BASE_DIR)
here::i_am("scripts/103_Fig1j_soxB_cross_cohort_forest_26.08.22.R")

SHARED_THEME <- normalizePath(file.path(
  here::here(), "..", "..", "_shared", "theme_nc_26.08.18.R"
))
source(SHARED_THEME)

INPUT_FILE <- here::here(
  "results_integrated", "sarcosine", "sarcosine_KO_per_sample_pooled.csv"
)
OUT_DIR <- normalizePath(file.path(
  here::here(), "..", "Manuscript_Final_Panels_26.08.22"
))

OUT_BASE <- file.path(OUT_DIR, "Fig1j_soxB_cross_cohort_forest")
OUT_PNG <- paste0(OUT_BASE, ".png")
OUT_PDF <- paste0(OUT_BASE, ".pdf")
OUT_SVG <- paste0(OUT_BASE, ".svg")
OUT_STATS <- paste0(OUT_BASE, "_statistics.csv")
OUT_CHECKSUM <- paste0(OUT_BASE, "_input_sha256.tsv")
OUT_SESSION <- paste0(OUT_BASE, "_sessionInfo.txt")
OUT_README <- paste0(OUT_BASE, "_README.txt")

if (!file.exists(INPUT_FILE)) stop("Missing input: ", INPUT_FILE)

dat <- read.csv(
  INPUT_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

required_columns <- c("Run.ID", "K00303", "Group", "Cohort")
missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns)) {
  stop("Input is missing: ", paste(missing_columns, collapse = ", "))
}
dat <- dat[, required_columns]

COHORTS <- c("PRJEB6070", "PRJEB10878", "PRJEB27928", "PRJNA429097")
GROUPS <- c("Healthy", "Cancer")

if (!setequal(unique(dat$Cohort), COHORTS)) stop("Unexpected cohort labels.")
if (!setequal(unique(dat$Group), GROUPS)) stop("Unexpected group labels.")
if (anyNA(dat) || any(!is.finite(dat$K00303))) {
  stop("Missing or non-finite required values detected.")
}
if (any(dat$K00303 < 0)) stop("Negative KO relative abundance detected.")
if (anyDuplicated(dat$Run.ID)) stop("Duplicated Run.ID values detected.")

sample_counts <- with(dat, table(Cohort, Group))
if (any(sample_counts < 3L)) stop("A cohort/group has fewer than three samples.")

rank_biserial <- function(healthy, cancer) {
  n_h <- length(healthy)
  n_c <- length(cancer)
  ranks <- rank(c(healthy, cancer), ties.method = "average")
  u_h <- sum(ranks[seq_len(n_h)]) - n_h * (n_h + 1) / 2
  2 * u_h / (n_h * n_c) - 1
}

rank_biserial_direct <- function(healthy, cancer) {
  mean(outer(healthy, cancer, FUN = function(x, y) sign(x - y)))
}

BOOT_REPS <- 5000L
BOOT_SEED <- 42L
results <- vector("list", length(COHORTS))

for (i in seq_along(COHORTS)) {
  cohort <- COHORTS[i]
  z <- dat[dat$Cohort == cohort, , drop = FALSE]
  healthy <- z$K00303[z$Group == "Healthy"]
  cancer <- z$K00303[z$Group == "Cancer"]

  effect <- rank_biserial(healthy, cancer)
  effect_direct <- rank_biserial_direct(healthy, cancer)
  if (abs(effect - effect_direct) > 1e-10) {
    stop(cohort, ": rank-biserial implementations disagree.")
  }

  wt <- suppressWarnings(wilcox.test(healthy, cancer, exact = FALSE))
  effect_from_w <- 2 * unname(wt$statistic) /
    (length(healthy) * length(cancer)) - 1
  if (!isTRUE(all.equal(effect, effect_from_w, tolerance = 1e-12))) {
    stop(cohort, ": effect size disagrees with Wilcoxon W.")
  }

  boot_effect <- replicate(
    BOOT_REPS,
    rank_biserial(
      sample(healthy, length(healthy), replace = TRUE),
      sample(cancer, length(cancer), replace = TRUE)
    )
  )
  ci <- unname(quantile(
    boot_effect,
    probs = c(0.025, 0.975),
    type = 7,
    names = FALSE
  ))

  results[[i]] <- data.frame(
    Cohort = cohort,
    n_Healthy = length(healthy),
    n_CRC = length(cancer),
    zero_Healthy = sum(healthy == 0),
    zero_CRC = sum(cancer == 0),
    median_Healthy = median(healthy),
    median_CRC = median(cancer),
    rank_biserial_Healthy_vs_CRC = effect,
    ci_95_low = ci[1],
    ci_95_high = ci[2],
    wilcox_p = wt$p.value,
    bootstrap_reps = BOOT_REPS,
    bootstrap_seed = BOOT_SEED,
    stringsAsFactors = FALSE
  )
}

stats <- do.call(rbind, results)
stats$wilcox_BH_q_across_4_cohorts <- p.adjust(stats$wilcox_p, method = "BH")
stats$Direction <- ifelse(
  stats$rank_biserial_Healthy_vs_CRC > 0,
  "Healthy higher",
  ifelse(stats$rank_biserial_Healthy_vs_CRC < 0, "CRC higher", "No direction")
)
stats$CI_status <- ifelse(
  stats$ci_95_low > 0 | stats$ci_95_high < 0,
  "95% CI excludes 0",
  "95% CI overlaps 0"
)

if (nrow(stats) != length(COHORTS) || anyDuplicated(stats$Cohort)) {
  stop("Expected one result per cohort.")
}
if (any(abs(stats$rank_biserial_Healthy_vs_CRC) > 1 + 1e-12)) {
  stop("Effect size outside [-1, 1].")
}
if (any(stats$ci_95_low > stats$ci_95_high)) stop("Reversed CI detected.")

write.csv(stats, OUT_STATS, row.names = FALSE, quote = TRUE)

plot_stats <- stats
plot_stats$Cohort <- factor(plot_stats$Cohort, levels = rev(COHORTS))
plot_stats$y <- as.numeric(plot_stats$Cohort)
plot_stats$Cohort_label <- sprintf(
  "%s\nH/CRC: %d/%d",
  as.character(plot_stats$Cohort),
  plot_stats$n_Healthy,
  plot_stats$n_CRC
)

write.table(
  data.frame(
    File = file.path(
      "results_integrated", "sarcosine",
      "sarcosine_KO_per_sample_pooled.csv"
    ),
    SHA256 = digest::digest(INPUT_FILE, algo = "sha256", file = TRUE),
    stringsAsFactors = FALSE
  ),
  OUT_CHECKSUM,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

row_background <- data.frame(
  ymin = c(0.5, 2.5),
  ymax = c(1.5, 3.5)
)

p <- ggplot(
  plot_stats,
  aes(x = rank_biserial_Healthy_vs_CRC, y = y)
) +
  geom_rect(
    data = row_background,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#F5F5F5",
    colour = NA
  ) +
  annotate(
    "segment",
    x = 0,
    xend = 0,
    y = 0.5,
    yend = 4.38,
    linewidth = NC_AXIS_LW,
    linetype = "22",
    colour = "#606060"
  ) +
  geom_segment(
    aes(x = ci_95_low, xend = ci_95_high, yend = y),
    linewidth = NC_AXIS_LW,
    colour = "#333333",
    lineend = "round"
  ) +
  geom_point(
    shape = 21,
    size = 3.2,
    stroke = NC_AXIS_LW,
    fill = COL_HEALTHY,
    colour = "#202020"
  ) +
  annotate(
    "text",
    x = -0.008,
    y = 4.92,
    label = "soxB",
    family = "Arial",
    fontface = "bold.italic",
    size = mm_text(NC_STRIP_PT),
    colour = "black",
    hjust = 1
  ) +
  annotate(
    "text",
    x = 0.008,
    y = 4.92,
    label = "(K00303)",
    family = "Arial",
    fontface = "bold",
    size = mm_text(NC_STRIP_PT),
    colour = "black",
    hjust = 0
  ) +
  annotate(
    "text",
    x = -0.38,
    y = 4.57,
    label = "CRC higher",
    family = "Arial",
    fontface = "bold",
    size = mm_text(NC_ANNOT_PT),
    colour = COL_CANCER
  ) +
  annotate(
    "text",
    x = 0.38,
    y = 4.57,
    label = "Healthy higher",
    family = "Arial",
    fontface = "bold",
    size = mm_text(NC_ANNOT_PT),
    colour = COL_HEALTHY
  ) +
  scale_x_continuous(
    name = expression("Rank-biserial effect size ("*r[rb]*")"),
    breaks = c(-0.5, 0, 0.5),
    labels = c("-0.5", "0", "0.5"),
    limits = c(-0.55, 0.55),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = plot_stats$y,
    labels = plot_stats$Cohort_label,
    limits = c(0.5, 5.08),
    expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_nc(base_pt = NC_TICK_PT) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(
      family = "Arial",
      size = NC_TICK_PT,
      colour = "black",
      lineheight = 0.88,
      margin = margin(r = 3)
    ),
    axis.title.x = element_text(
      family = "Arial",
      size = NC_TITLE_PT,
      colour = "black",
      margin = margin(t = 3)
    ),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(t = 1, r = 3, b = 1, l = 4)
  )

WIDTH_IN <- 2.80
HEIGHT_IN <- 2.10
DPI <- 600L

ggsave(
  OUT_PNG,
  p,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  dpi = DPI,
  bg = "white",
  device = ragg::agg_png
)
ggsave(
  OUT_PDF,
  p,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  bg = "white",
  device = function(filename, width, height, bg = "white", ...) {
    grDevices::quartz(
      type = "pdf",
      file = filename,
      width = width,
      height = height,
      family = "Arial",
      bg = bg
    )
  }
)
ggsave(
  OUT_SVG,
  p,
  width = WIDTH_IN,
  height = HEIGHT_IN,
  bg = "white",
  device = svglite::svglite
)

expected_files <- c(OUT_PNG, OUT_PDF, OUT_SVG, OUT_STATS, OUT_CHECKSUM)
if (!all(file.exists(expected_files))) stop("One or more outputs were not written.")
if (any(file.info(c(OUT_PNG, OUT_PDF, OUT_SVG))$size <= 5000)) {
  stop("An artwork file is unexpectedly small.")
}
if (any(file.info(c(OUT_STATS, OUT_CHECKSUM))$size <= 100)) {
  stop("A provenance table is unexpectedly small.")
}

png_dim <- dim(png::readPNG(OUT_PNG))
expected_dim <- c(round(HEIGHT_IN * DPI), round(WIDTH_IN * DPI))
if (!identical(as.integer(png_dim[1:2]), as.integer(expected_dim))) {
  stop("Unexpected PNG dimensions: ", paste(png_dim[1:2], collapse = " x "))
}

readme_lines <- c(
  "Proposed Figure 1j: four-cohort soxB (K00303) effect-size forest plot",
  "",
  "Input:",
  "  results_integrated/sarcosine/sarcosine_KO_per_sample_pooled.csv",
  "  HGMT WGS KO relative abundance; one row per Healthy or CRC sample.",
  "",
  "Effect definition:",
  "  r_rb = P(Healthy > CRC) - P(Healthy < CRC); ties contribute zero.",
  "  Positive values mean that soxB is higher in Healthy controls.",
  "",
  "Uncertainty:",
  "  5,000 group-stratified percentile bootstrap repetitions; seed 42.",
  "",
  "Interpretation:",
  "  Cohort-specific estimates only; no pooled meta-analytic estimate is fitted.",
  "  The figure should therefore be described as a cohort-stratified effect-size",
  "  forest plot, not as a formal meta-analysis.",
  "",
  sprintf("Final footprint: %.2f x %.2f inches; PNG: %d dpi.", WIDTH_IN, HEIGHT_IN, DPI),
  "Do not rescale after insertion if exact final type sizes are required.",
  "",
  "Reproduce from Healthy_vs_Cancer_4_CRC_cohorts_integrated:",
  "  Rscript scripts/103_Fig1j_soxB_cross_cohort_forest_26.08.22.R"
)
writeLines(readme_lines, OUT_README)
capture.output(sessionInfo(), file = OUT_SESSION)

cat("Input dimensions: ", nrow(dat), " samples x ", ncol(dat), " selected columns\n", sep = "")
print(sample_counts)
print(stats[, c(
  "Cohort", "n_Healthy", "n_CRC", "rank_biserial_Healthy_vs_CRC",
  "ci_95_low", "ci_95_high", "wilcox_p", "wilcox_BH_q_across_4_cohorts"
)], row.names = FALSE)
cat("Wrote:\n", paste(c(expected_files, OUT_README, OUT_SESSION), collapse = "\n"), "\n")
