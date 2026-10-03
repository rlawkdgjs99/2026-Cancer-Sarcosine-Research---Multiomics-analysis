#!/usr/bin/env Rscript

# Four-cohort Healthy-versus-CRC effect-size forests for five sarcosine KOs.
# This complements, rather than replaces, the pooled five-KO panels.
# r_rb = P(Healthy > CRC) - P(Healthy < CRC); ties contribute zero.

suppressPackageStartupMessages({
  library(digest)
  library(dplyr)
  library(ggplot2)
  library(ragg)
})

options(stringsAsFactors = FALSE)
set.seed(42L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", script_arg))
base_dir <- normalizePath(file.path(dirname(script_path), ".."))
analysis_root <- normalizePath(file.path(base_dir, ".."))
project_root <- normalizePath(file.path(base_dir, "..", ".."))
workspace_root <- normalizePath(file.path(base_dir, "..", "..", ".."))

theme_path <- file.path(project_root, "_shared", "theme_nc_26.08.18.R")
if (!file.exists(theme_path)) stop("Missing shared theme: ", theme_path)
source(theme_path)

input_path <- file.path(base_dir, "results_integrated", "sarcosine", "sarcosine_KO_per_sample_pooled.csv")
pooled_frozen_path <- file.path(base_dir, "results_integrated", "sarcosine", "sarcosine_KO_comparison_filtered_pooled.csv")
if (!file.exists(input_path) || !file.exists(pooled_frozen_path)) stop("Missing pooled source/frozen file")

out_dir <- file.path(base_dir, "results_integrated", "publication_style_plots_REMAINING_5KO_4COHORT_26.09.01")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
data_copy_dir <- file.path(
  analysis_root, "사용데이터_모음",
  "CRC_WGS_remaining_5KO_4cohort_plot_source_data_26.09.01"
)
dir.create(data_copy_dir, recursive = TRUE, showWarnings = FALSE)
transport_dir <- file.path(workspace_root, "Manuscript작업", "FigDesign&Manuscript")
if (!dir.exists(transport_dir)) stop("Missing manuscript transport directory")

ko_info <- tibble::tribble(
  ~KO,      ~Display,              ~Slug,
  "K00301", "Sarcosine oxidase", "K00301_Sarcosine_oxidase",
  "K00302", "soxA",               "K00302_soxA",
  "K00305", "soxG",               "K00305_soxG",
  "K00306", "PIPOX",              "K00306_PIPOX",
  "K00315", "DMGDH",              "K00315_DMGDH"
)
cohorts <- c("PRJEB6070", "PRJEB10878", "PRJEB27928", "PRJNA429097")
frozen_paths <- c(
  PRJEB6070 = file.path(analysis_root, "PRJEB6070_CRC_AdenomatousPolyps", "results_sarcosine", "sarcosine_KO_comparison.csv"),
  PRJEB10878 = file.path(analysis_root, "PRJEB10878_CRC", "results_sarcosine", "sarcosine_KO_comparison.csv"),
  PRJEB27928 = file.path(analysis_root, "PRJEB27928_CRC", "results_sarcosine", "sarcosine_KO_comparison.csv"),
  PRJNA429097 = file.path(analysis_root, "PRJNA429097_CRC", "results_sarcosine", "sarcosine_KO_comparison.csv")
)
if (!all(file.exists(frozen_paths))) stop("Missing cohort-level frozen KO summary")

dat <- read.csv(input_path, check.names = FALSE)
pooled_frozen <- read.csv(pooled_frozen_path, check.names = FALSE)
required <- c("Run.ID", ko_info$KO, "Group", "Cohort")
if (!all(required %in% names(dat))) stop("Input lacks required columns")
if (nrow(dat) != 1647L || anyDuplicated(dat$Run.ID)) stop("Expected 1,647 unique Run IDs")
if (!setequal(unique(dat$Group), c("Healthy", "Cancer"))) stop("Unexpected Group labels")
if (!setequal(unique(dat$Cohort), cohorts)) stop("Unexpected cohort labels")
expected_counts <- matrix(
  c(476, 590, 54, 74, 120, 140, 95, 98), nrow = 4, byrow = TRUE,
  dimnames = list(cohorts, c("Healthy", "Cancer"))
)
observed_counts <- with(dat, table(factor(Cohort, levels = cohorts), factor(Group, levels = c("Healthy", "Cancer"))))
if (!all(unname(observed_counts) == unname(expected_counts))) stop("Cohort/group counts changed")

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

bootstrap_ci <- function(healthy, cancer, reps = 5000L, seed = 42L) {
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  boot_effect <- replicate(
    reps,
    rank_biserial(
      sample(healthy, length(healthy), replace = TRUE),
      sample(cancer, length(cancer), replace = TRUE)
    )
  )
  unname(quantile(boot_effect, c(0.025, 0.975), type = 7))
}

standardize_png <- function(path) {
  tmp <- tempfile(fileext = ".png")
  cmd <- sprintf(
    "/usr/bin/sips -s format png -s formatOptions best -s profile '/System/Library/ColorSync/Profiles/sRGB Profile.icc' %s --out %s >/dev/null",
    shQuote(path), shQuote(tmp)
  )
  if (system(cmd) != 0L || !file.exists(tmp)) stop("PNG standardization failed: ", path)
  if (!file.copy(tmp, path, overwrite = TRUE)) stop("Could not replace PNG: ", path)
  unlink(tmp)
  Sys.chmod(path, "0644")
  system2("/usr/bin/xattr", c("-c", shQuote(path)), stdout = FALSE, stderr = FALSE)
  system2("/usr/bin/chflags", c("nohidden", shQuote(path)), stdout = FALSE, stderr = FALSE)
}

copy_exact_visible <- function(source, destination) {
  if (!file.copy(source, destination, overwrite = TRUE)) stop("Copy failed: ", destination)
  if (!identical(digest(source, "sha256", file = TRUE), digest(destination, "sha256", file = TRUE))) {
    stop("Hash mismatch after copy: ", destination)
  }
  Sys.chmod(destination, "0644")
  system2("/usr/bin/xattr", c("-c", shQuote(destination)), stdout = FALSE, stderr = FALSE)
  system2("/usr/bin/chflags", c("nohidden", shQuote(destination)), stdout = FALSE, stderr = FALSE)
}

verify_frozen_row <- function(result_row) {
  cohort <- as.character(result_row$Cohort)
  frozen <- read.csv(frozen_paths[[cohort]], check.names = FALSE)
  frow <- frozen[frozen$KO == result_row$KO, , drop = FALSE]
  if (!result_row$Estimable) {
    if (nrow(frow) > 1L) stop("Duplicate all-zero frozen row: ", result_row$KO, " / ", cohort)
    if (nrow(frow) == 1L) {
      vals <- unlist(frow[, intersect(c("Mean_Healthy", "Mean_Cancer", "Median_Healthy", "Median_Cancer"), names(frow)), drop = FALSE])
      if (any(abs(vals) > 1e-15, na.rm = TRUE)) stop("All-zero result disagrees with frozen table")
    }
    return(invisible(TRUE))
  }
  if (nrow(frow) != 1L) stop("Missing/duplicate frozen row: ", result_row$KO, " / ", cohort)
  checks <- c(
    Mean_Healthy = result_row$mean_Healthy,
    Mean_Cancer = result_row$mean_CRC,
    Median_Healthy = result_row$median_Healthy,
    Median_Cancer = result_row$median_CRC,
    p_value = result_row$wilcox_p
  )
  for (nm in names(checks)) {
    if (!isTRUE(all.equal(as.numeric(checks[[nm]]), as.numeric(frow[[nm]]), tolerance = 1e-10))) {
      stop("Frozen mismatch: ", result_row$KO, " / ", cohort, " / ", nm)
    }
  }
  invisible(TRUE)
}

results_all <- vector("list", nrow(ko_info))
for (i in seq_len(nrow(ko_info))) {
  meta <- ko_info[i, ]
  ko <- meta$KO
  pooled_meta <- pooled_frozen[pooled_frozen$KO == ko, , drop = FALSE]
  if (nrow(pooled_meta) != 1L) stop("Missing pooled prevalence metadata: ", ko)

  rows <- lapply(cohorts, function(cohort) {
    z <- dat[dat$Cohort == cohort, , drop = FALSE]
    healthy <- as.numeric(z[z$Group == "Healthy", ko])
    cancer <- as.numeric(z[z$Group == "Cancer", ko])
    if (any(!is.finite(c(healthy, cancer))) || any(c(healthy, cancer) < 0)) stop("Invalid KO value")
    all_zero <- all(c(healthy, cancer) == 0)
    effect <- if (all_zero) 0 else rank_biserial(healthy, cancer)
    direct <- if (all_zero) 0 else rank_biserial_direct(healthy, cancer)
    if (abs(effect - direct) > 1e-10) stop("Rank-biserial cross-check failed: ", ko, " / ", cohort)
    if (all_zero) {
      ci <- c(0, 0)
      p_value <- NA_real_
    } else {
      wt <- suppressWarnings(wilcox.test(healthy, cancer, exact = FALSE))
      effect_from_w <- 2 * unname(wt$statistic) / (length(healthy) * length(cancer)) - 1
      if (!isTRUE(all.equal(effect, effect_from_w, tolerance = 1e-12))) stop("Wilcoxon effect check failed")
      ci <- bootstrap_ci(healthy, cancer, reps = 5000L, seed = 42L)
      p_value <- wt$p.value
    }
    data.frame(
      KO = ko,
      Enzyme = pooled_meta$Enzyme,
      EC = pooled_meta$EC,
      Role = pooled_meta$Role,
      Cohort = cohort,
      n_Healthy = length(healthy),
      n_CRC = length(cancer),
      zero_Healthy = sum(healthy == 0),
      zero_CRC = sum(cancer == 0),
      prevalence_Healthy_percent = 100 * mean(healthy > 0),
      prevalence_CRC_percent = 100 * mean(cancer > 0),
      overall_pooled_prevalence_percent = pooled_meta$Prev_Overall,
      passes_10_percent_overall_pooled_prevalence = pooled_meta$Pass_Prevalence,
      median_Healthy = median(healthy),
      median_CRC = median(cancer),
      mean_Healthy = mean(healthy),
      mean_CRC = mean(cancer),
      rank_biserial_Healthy_vs_CRC = effect,
      ci_95_low = ci[1],
      ci_95_high = ci[2],
      wilcox_p = p_value,
      Estimable = !all_zero,
      Estimability_note = if (all_zero) "All zero in both groups; Wilcoxon test and effect direction not estimable" else "Estimable",
      bootstrap_reps = if (all_zero) 0L else 5000L,
      bootstrap_seed = 42L,
      stringsAsFactors = FALSE
    )
  }) |> bind_rows()

  estimable <- which(rows$Estimable & is.finite(rows$wilcox_p))
  rows$wilcox_BH_q_across_estimable_cohorts_within_KO <- NA_real_
  rows$wilcox_BH_q_across_estimable_cohorts_within_KO[estimable] <- p.adjust(rows$wilcox_p[estimable], "BH")
  rows$Direction <- ifelse(
    !rows$Estimable, "No variation",
    ifelse(rows$rank_biserial_Healthy_vs_CRC > 0, "Higher in Healthy",
           ifelse(rows$rank_biserial_Healthy_vs_CRC < 0, "Higher in CRC", "No direction"))
  )
  rows$CI_status <- ifelse(
    !rows$Estimable, "Not estimable: all zero",
    ifelse(rows$ci_95_low > 0 | rows$ci_95_high < 0, "Excludes zero", "Includes zero")
  )
  for (j in seq_len(nrow(rows))) verify_frozen_row(rows[j, , drop = FALSE])

  plot_data_path <- file.path(out_dir, paste0(meta$Slug, "_4cohort_forest_plot_data.csv"))
  write.csv(rows, plot_data_path, row.names = FALSE, na = "")
  raw_data <- data.frame(
    Run.ID = dat$Run.ID,
    Cohort = dat$Cohort,
    Group = dat$Group,
    KO = ko,
    KO_relative_abundance = dat[[ko]],
    stringsAsFactors = FALSE
  )
  raw_data_path <- file.path(out_dir, paste0(meta$Slug, "_4cohort_raw_values.csv"))
  write.csv(raw_data, raw_data_path, row.names = FALSE, na = "")

  plot_rows <- rows
  plot_rows$y <- match(plot_rows$Cohort, rev(cohorts))
  plot_rows$Cohort_label <- sprintf(
    "%s\nH/CRC: %d/%d%s", plot_rows$Cohort, plot_rows$n_Healthy, plot_rows$n_CRC,
    ifelse(plot_rows$Estimable, "", " · all zero")
  )
  fill_values <- ifelse(
    plot_rows$passes_10_percent_overall_pooled_prevalence,
    ifelse(plot_rows$Direction == "Higher in Healthy", COL_HEALTHY,
           ifelse(plot_rows$Direction == "Higher in CRC", COL_CANCER, "white")),
    "white"
  )
  outline_values <- ifelse(
    !plot_rows$Estimable, "#929292",
    ifelse(plot_rows$Direction == "Higher in Healthy", COL_HEALTHY, COL_CANCER)
  )

  row_background <- data.frame(ymin = c(0.5, 2.5), ymax = c(1.5, 3.5))
  p <- ggplot(plot_rows, aes(x = rank_biserial_Healthy_vs_CRC, y = y)) +
    geom_rect(
      data = row_background,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "#F5F5F5", colour = NA
    ) +
    annotate("segment", x = 0, xend = 0, y = 0.5, yend = 4.38,
             linewidth = NC_AXIS_LW, linetype = "22", colour = "#606060") +
    geom_segment(
      data = plot_rows |> filter(Estimable),
      aes(x = ci_95_low, xend = ci_95_high, yend = y),
      linewidth = NC_AXIS_LW, colour = "#8F8F8F", lineend = "round"
    ) +
    geom_point(shape = 21, size = 3.2, stroke = 0.85, fill = fill_values, colour = outline_values) +
    geom_text(
      data = plot_rows |> filter(!Estimable),
      aes(x = 0.035, label = "all zero"), hjust = 0,
      family = "Arial", size = 1.75, colour = "#808080"
    ) +
    annotate("text", x = 0, y = 4.98, label = paste0(meta$Display, " (", ko, ")"),
             family = "Arial", fontface = ifelse(meta$Display %in% c("soxA", "soxG"), "bold.italic", "bold"),
             size = ifelse(meta$Display == "Sarcosine oxidase", 2.7, 3.0)) +
    annotate("text", x = -0.38, y = 4.60, label = "CRC higher", family = "Arial",
             fontface = "bold", size = 2.35, colour = COL_CANCER) +
    annotate("text", x = 0.38, y = 4.60, label = "Healthy higher", family = "Arial",
             fontface = "bold", size = 2.35, colour = COL_HEALTHY) +
    scale_x_continuous(
      name = expression("Rank-biserial effect size ("*r[rb]*")"),
      breaks = c(-0.5, 0, 0.5), limits = c(-0.55, 0.55), expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = plot_rows$y, labels = plot_rows$Cohort_label,
      limits = c(0.5, 5.08), expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      y = NULL,
      caption = if (!isTRUE(pooled_meta$Pass_Prevalence)) {
        if (any(!rows$Estimable))
          paste0("Open colored circles: pooled prevalence ", pooled_meta$Prev_Overall, "% (<10%); gray = all zero")
        else paste0("Open circles: pooled prevalence ", pooled_meta$Prev_Overall, "% (<10%)")
      } else if (any(!rows$Estimable)) {
        "Gray open circle: all zero in both groups"
      } else NULL
    ) +
    theme_nc(base_pt = NC_TICK_PT) +
    theme(
      plot.margin = margin(7, 7, 5, 7, unit = "pt"),
      axis.text.y = element_text(size = 5.35, lineheight = 0.9, hjust = 1),
      axis.text.x = element_text(size = 6.0),
      axis.title.x = element_text(size = 6.1, margin = margin(t = 4)),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 4.5, colour = "#666666", hjust = 0.5, margin = margin(t = 3))
    )

  png_path <- file.path(out_dir, paste0(meta$Slug, "_4cohort_forest.png"))
  pdf_path <- file.path(out_dir, paste0(meta$Slug, "_4cohort_forest.pdf"))
  ggsave(png_path, p, width = 2.80, height = 2.25, units = "in", dpi = 600, device = ragg::agg_png, bg = "white")
  ggsave(pdf_path, p, width = 2.80, height = 2.25, units = "in", device = cairo_pdf, bg = "white")
  standardize_png(png_path)

  transport_path <- file.path(transport_dir, paste0("CRC_WGS_4COHORT_", meta$Slug, "_forest.png"))
  copy_exact_visible(png_path, transport_path)
  copy_exact_visible(
    plot_data_path,
    file.path(data_copy_dir, paste0("CRC_WGS_4COHORT_", meta$Slug, "_forest_plot_data.csv"))
  )
  copy_exact_visible(
    raw_data_path,
    file.path(data_copy_dir, paste0("CRC_WGS_4COHORT_", meta$Slug, "_raw_values.csv"))
  )
  results_all[[i]] <- rows
}

combined <- bind_rows(results_all)
combined_path <- file.path(out_dir, "CRC_WGS_remaining_5KO_4cohort_forest_plot_data_combined.csv")
write.csv(combined, combined_path, row.names = FALSE, na = "")
copy_exact_visible(
  combined_path,
  file.path(data_copy_dir, "CRC_WGS_remaining_5KO_4cohort_forest_plot_data_combined.csv")
)

input_files <- c(input_path, pooled_frozen_path, unname(frozen_paths))
input_manifest <- data.frame(
  File = normalizePath(input_files),
  Bytes = file.info(input_files)$size,
  SHA256 = vapply(input_files, digest, character(1), algo = "sha256", file = TRUE),
  stringsAsFactors = FALSE
)
write.table(input_manifest, file.path(out_dir, "INPUT_SHA256.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

readme <- c(
  "# Remaining five KO four-cohort Healthy-versus-CRC effect-size forests",
  "",
  "Generated: 2026-09-01",
  "",
  "These figures show separate effects for PRJEB6070, PRJEB10878, PRJEB27928 and PRJNA429097.",
  "They complement the pooled five-KO panels and do not calculate a pooled/meta-analysis estimate.",
  "",
  "- Effect: rank-biserial r_rb = P(Healthy > CRC) - P(Healthy < CRC).",
  "- Interval: 95% group-stratified percentile-bootstrap CI; 5,000 resamples; seed 42.",
  "- P: two-sided Wilcoxon rank-sum test within each cohort.",
  "- q: Benjamini-Hochberg correction across the estimable cohort-specific tests within each KO.",
  "- All-zero cohort/group comparisons are explicitly non-estimable and retain P/q as NA.",
  "- Open colored circles denote an KO with pooled prevalence below 10%; gray open circles denote all-zero comparisons.",
  "",
  "Each KO has one four-row plot/statistics CSV and one 1,647-row raw-value CSV.",
  "Means, medians and raw Wilcoxon P values were required to match the frozen cohort-level KO summaries where estimable."
)
writeLines(readme, file.path(out_dir, "README.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"), useBytes = TRUE)

output_files <- list.files(out_dir, full.names = TRUE)
output_files <- output_files[basename(output_files) != "OUTPUT_SHA256.tsv"]
output_manifest <- data.frame(
  File = basename(output_files),
  Bytes = file.info(output_files)$size,
  SHA256 = vapply(output_files, digest, character(1), algo = "sha256", file = TRUE),
  stringsAsFactors = FALSE
)
write.table(output_manifest, file.path(out_dir, "OUTPUT_SHA256.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

message("Created five four-cohort KO panels: ", out_dir)
message("Copied per-KO plot and raw-value CSVs: ", data_copy_dir)
