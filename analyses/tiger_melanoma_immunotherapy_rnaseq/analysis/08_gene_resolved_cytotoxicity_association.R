#!/usr/bin/env Rscript
# =============================================================================
# Candidate Figure 4 panel: gene-resolved decomposition of the association
# between host sarcosine degradation and tumour cytolytic activity (CYT).
#
# Frozen design (2026-08-24)
#   Biological question:
#     Which host sarcosine-metabolism components account for the already
#     reported association between the degradation score and tumour CYT?
#   Observational unit:
#     One pretreatment melanoma biopsy per patient (n = 73).
#   Effect measure:
#     Two-sided Spearman correlation (rho).
#   Intervals:
#     Patient-level percentile bootstrap 95% CIs (10,000 resamples).
#   Multiplicity:
#     BH correction across the same six exploratory tests defined in
#     102_Fig2d_CYT_association_26.08.23.R. Existing rho/P/q values must be
#     reproduced exactly (numerical tolerance 1e-12) before plotting.
#   Interpretation limit:
#     These are unadjusted, cross-sectional expression associations. They do
#     not establish causality or cell-type specificity.
#
# Inputs are existing, immutable outputs from analysis 102. This script writes
# only to a new dated candidate folder and does not edit the manuscript or PPT.
# =============================================================================

set.seed(4201)

suppressPackageStartupMessages({
  library(ggplot2)
  library(here)
})

# Locate the project from this script, then register a portable here() root.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run this analysis with Rscript.")
script_file <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
setwd(project_dir)
here::i_am(file.path("analysis", basename(script_file)))

input_dir <- here::here(
  "results", "manuscript_figures", "Fig2d_CYT_candidate_26.08.23"
)
source_file <- file.path(input_dir, "SourceData_Fig2d_CYT_PRE73.csv")
reference_file <- file.path(input_dir, "secondary_component_correlations.csv")
out_dir <- here::here(
  "results", "manuscript_figures",
  "Fig4_gene_resolved_CYT_candidate_26.08.24"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_file) || !file.exists(reference_file)) {
  stop("Required analysis-102 input is missing.")
}

source_data <- read.csv(source_file, stringsAsFactors = FALSE, check.names = FALSE)
reference <- read.csv(reference_file, stringsAsFactors = FALSE)

required_source_columns <- c(
  "sample_id", "patient_name", "timepoint",
  "SARDH_log2_FPKM_plus1", "PIPOX_log2_FPKM_plus1",
  "GNMT_log2_FPKM_plus1", "DMGDH_log2_FPKM_plus1",
  "GZMA_TPM", "PRF1_TPM", "Degradation_score", "CYT_score"
)
required_reference_columns <- c(
  "predictor", "outcome", "n", "spearman_rho", "p_nominal", "p_BH"
)
if (!all(required_source_columns %in% names(source_data))) {
  stop("Source-data columns are missing: ", paste(
    setdiff(required_source_columns, names(source_data)), collapse = ", "
  ))
}
if (!all(required_reference_columns %in% names(reference))) {
  stop("Reference-result columns are missing: ", paste(
    setdiff(required_reference_columns, names(reference)), collapse = ", "
  ))
}
stopifnot(
  nrow(source_data) == 73L,
  length(unique(source_data$patient_name)) == 73L,
  !anyDuplicated(source_data$patient_name),
  !anyDuplicated(source_data$sample_id),
  all(source_data$timepoint == "PRE"),
  nrow(reference) == 6L,
  !anyNA(source_data[, required_source_columns]),
  !anyNA(reference[, required_reference_columns])
)

pair_spec <- data.frame(
  predictor = c(
    "SARDH", "PIPOX", "GNMT", "DMGDH",
    "Degradation_score", "Degradation_score"
  ),
  outcome = c(
    "CYT_score", "CYT_score", "CYT_score", "CYT_score",
    "GZMA_TPM", "PRF1_TPM"
  ),
  predictor_column = c(
    "SARDH_log2_FPKM_plus1", "PIPOX_log2_FPKM_plus1",
    "GNMT_log2_FPKM_plus1", "DMGDH_log2_FPKM_plus1",
    "Degradation_score", "Degradation_score"
  ),
  outcome_column = c(
    "CYT_score", "CYT_score", "CYT_score", "CYT_score",
    "GZMA_TPM", "PRF1_TPM"
  ),
  display_label = c(
    "SARDH", "PIPOX", "GNMT", "DMGDH",
    "Score - GZMA", "Score - PRF1"
  ),
  section = c(
    rep("Enzyme expression vs CYT", 4L),
    rep("Degradation score vs CYT genes", 2L)
  ),
  stringsAsFactors = FALSE
)

reference$key <- paste(reference$predictor, reference$outcome, sep = "__")
pair_spec$key <- paste(pair_spec$predictor, pair_spec$outcome, sep = "__")
if (!setequal(reference$key, pair_spec$key)) {
  stop("The six frozen analysis pairs no longer match the reference results.")
}
reference <- reference[match(pair_spec$key, reference$key), , drop = FALSE]
stopifnot(identical(reference$key, pair_spec$key))

bootstrap_spearman <- function(x, y, replicates = 10000L, seed) {
  stopifnot(
    is.numeric(x), is.numeric(y), length(x) == length(y),
    length(x) == 73L, !anyNA(x), !anyNA(y)
  )
  set.seed(seed)
  n <- length(x)
  rho <- replicate(replicates, {
    index <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[index], y[index], method = "spearman"))
  })
  rho <- rho[is.finite(rho)]
  if (length(rho) < 0.99 * replicates) {
    stop("More than 1% of bootstrap correlations were non-finite.")
  }
  c(
    ci_lower = unname(quantile(rho, 0.025, names = FALSE)),
    ci_upper = unname(quantile(rho, 0.975, names = FALSE)),
    finite_replicates = length(rho)
  )
}

recomputed <- lapply(seq_len(nrow(pair_spec)), function(i) {
  x <- source_data[[pair_spec$predictor_column[i]]]
  y <- source_data[[pair_spec$outcome_column[i]]]
  test <- suppressWarnings(cor.test(
    x, y, method = "spearman", exact = FALSE, alternative = "two.sided"
  ))
  ci <- bootstrap_spearman(x, y, replicates = 10000L, seed = 4200L + i)
  data.frame(
    pair_spec[i, c(
      "predictor", "outcome", "display_label", "section"
    )],
    n = length(x),
    spearman_rho = unname(test$estimate),
    ci_lower = ci[["ci_lower"]],
    ci_upper = ci[["ci_upper"]],
    p_nominal = test$p.value,
    finite_bootstrap_replicates = ci[["finite_replicates"]],
    stringsAsFactors = FALSE
  )
})
results <- do.call(rbind, recomputed)
results$p_BH <- p.adjust(results$p_nominal, method = "BH")

# Hard validation gate: plotting stops if the saved analysis-102 values are not
# exactly reproduced to floating-point tolerance.
for (column in c("spearman_rho", "p_nominal", "p_BH")) {
  if (!isTRUE(all.equal(
    results[[column]], reference[[column]], tolerance = 1e-12,
    check.attributes = FALSE
  ))) {
    stop("Recomputed ", column, " does not reproduce analysis 102.")
  }
}
stopifnot(
  all(results$ci_lower >= -1), all(results$ci_upper <= 1),
  all(results$ci_lower < results$spearman_rho),
  all(results$ci_upper > results$spearman_rho),
  all(results$finite_bootstrap_replicates >= 9900L)
)

write.csv(
  results,
  file.path(out_dir, "SourceData_Fig4_gene_resolved_CYT_candidate.csv"),
  row.names = FALSE
)

# Use the shared final-size typography and linewidth system.
workspace_root <- normalizePath(file.path(project_dir, "..", "..", ".."))
theme_candidates <- Sys.glob(file.path(
  workspace_root, "*Metabolomics*", "_shared", "theme_nc_26.08.18.R"
))
if (length(theme_candidates) != 1L) {
  stop("Expected exactly one shared theme; found ", length(theme_candidates))
}
source(theme_candidates)

results$display_label <- factor(
  results$display_label,
  levels = rev(pair_spec$display_label)
)
results$evidence <- factor(
  ifelse(results$p_BH < 0.05, "q < 0.05", "q ≥ 0.05"),
  levels = c("q < 0.05", "q ≥ 0.05")
)

# The interval scale occupies the full panel width. BH significance is encoded
# by filled versus open points; exact q values remain in the source-data table.
# Dimensions match the current Figure 4c PowerPoint slot (1.927 x 2.648 in)
# without post-hoc shrinking.
x_limits <- c(-0.48, 0.82)

panel <- ggplot(results, aes(y = display_label, x = spearman_rho)) +
  annotate(
    "rect", xmin = -Inf, xmax = Inf, ymin = 0.55, ymax = 2.45,
    fill = "#F5F5F5", colour = NA
  ) +
  annotate(
    "segment", x = 0, xend = 0, y = 0.55, yend = 2.25,
    colour = "#8B8B8B", linewidth = pt_lw(1.0), linetype = "22"
  ) +
  annotate(
    "segment", x = 0, xend = 0, y = 2.75, yend = 6.25,
    colour = "#8B8B8B", linewidth = pt_lw(1.0), linetype = "22"
  ) +
  geom_segment(
    aes(x = ci_lower, xend = ci_upper, yend = display_label),
    colour = "#4A4A4A", linewidth = pt_lw(1.0), lineend = "round"
  ) +
  geom_point(
    aes(fill = evidence), shape = 21, colour = "#252525",
    size = 2.25, stroke = pt_lw(0.8)
  ) +
  annotate(
    "text", x = x_limits[1], y = 6.48,
    label = "Enzyme expression vs CYT", hjust = 0,
    family = "Arial", fontface = "bold", size = mm_text(NC_STRIP_PT)
  ) +
  annotate(
    "text", x = x_limits[1], y = 2.48,
    label = "Score vs GZMA/PRF1", hjust = 0,
    family = "Arial", fontface = "bold", size = mm_text(NC_STRIP_PT)
  ) +
  scale_fill_manual(
    values = c("q < 0.05" = "#3B3B3B", "q ≥ 0.05" = "white"),
    breaks = c("q < 0.05", "q ≥ 0.05"),
    labels = c("< 0.05", ">= 0.05"),
    name = "BH q", drop = FALSE
  ) +
  scale_x_continuous(
    limits = x_limits, breaks = c(-0.4, 0, 0.4, 0.8),
    expand = expansion(mult = 0)
  ) +
  scale_y_discrete(expand = expansion(add = c(0.35, 0.75))) +
  labs(x = expression("Spearman " * rho * " (95% CI)"), y = NULL) +
  theme_nc() +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(
      family = "Arial", size = NC_TICK_PT, hjust = 1,
      margin = margin(r = 2.5)
    ),
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.title = element_text(
      family = "Arial", size = 5.5, face = "bold"
    ),
    legend.text = element_text(family = "Arial", size = 5.5),
    legend.key.width = grid::unit(7, "pt"),
    legend.spacing.x = grid::unit(1.5, "pt"),
    legend.box.spacing = grid::unit(0, "pt"),
    panel.grid = element_blank(),
    plot.margin = margin(2, 2, 1, 1)
  ) +
  guides(fill = guide_legend(
    nrow = 1, byrow = TRUE,
    override.aes = list(shape = 21, size = 2.0, stroke = pt_lw(0.8))
  ))

figure_stem <- file.path(out_dir, "Fig4_gene_resolved_CYT_candidate")
ggsave(
  paste0(figure_stem, ".png"), panel,
  width = 1.93, height = 2.65, units = "in", dpi = 600,
  bg = "white", device = ragg::agg_png
)

# Vector PDF for final redraw/placement. Use the same verified Arial Type1 map
# as analysis 102 because cairo_pdf fails on this host's R build.
grDevices::pdfFonts(Arial = grDevices::pdfFonts("ArialMT")[[1]])
ggsave(
  paste0(figure_stem, ".pdf"), panel,
  width = 1.93, height = 2.65, units = "in", bg = "white",
  device = grDevices::pdf, family = "Arial", useDingbats = FALSE
)

writeLines(
  c(
    "Candidate Figure 4 gene-resolved CYT panel",
    "Observational unit: 73 unique pretreatment melanoma biopsies.",
    "Effect: two-sided Spearman rho; intervals: 10,000 patient bootstrap resamples.",
    "Multiplicity: BH across the six frozen exploratory component tests.",
    "All saved analysis-102 rho, nominal P and BH q values reproduced within 1e-12.",
    "Final-size footprint: 1.93 x 2.65 inches; PNG is 1158 x 1590 pixels at 600 dpi.",
    "Interpretation: cross-sectional association only; no causal or cell-type-specific claim."
  ),
  file.path(out_dir, "README_candidate.txt")
)

expected_png <- paste0(figure_stem, ".png")
expected_pdf <- paste0(figure_stem, ".pdf")
stopifnot(
  file.exists(expected_png), file.info(expected_png)$size > 5000,
  file.exists(expected_pdf), file.info(expected_pdf)$size > 5000
)

if (!requireNamespace("png", quietly = TRUE)) {
  stop("Package 'png' is required for the pixel-dimension assertion.")
}
png_info <- png::readPNG(expected_png, info = TRUE)
stopifnot(identical(dim(png_info)[1:2], c(1590L, 1158L)))

cat("PASS: gene-resolved CYT candidate created in\n", out_dir, "\n")
print(results[, c(
  "display_label", "n", "spearman_rho", "ci_lower", "ci_upper", "p_BH"
)], row.names = FALSE, digits = 5)
