#!/usr/bin/env Rscript

# Pooled four-cohort Healthy-versus-CRC effect plots for the five sarcosine KOs
# not already shown as dedicated soxB (K00303) and creatinase (K08688) panels.
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
frozen_path <- file.path(base_dir, "results_integrated", "sarcosine", "sarcosine_KO_comparison_filtered_pooled.csv")
if (!file.exists(input_path) || !file.exists(frozen_path)) stop("Missing pooled source file")

out_dir <- file.path(base_dir, "results_integrated", "publication_style_plots_REMAINING_5KO_POOLED_26.09.01")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
data_copy_dir <- file.path(
  analysis_root, "사용데이터_모음",
  "CRC_WGS_remaining_5KO_pooled_plot_source_data_26.09.01"
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
seven_targeted <- c("K00301", "K00302", "K00303", "K00305", "K00306", "K00315", "K08688")
cohorts_expected <- c("PRJEB6070", "PRJEB10878", "PRJEB27928", "PRJNA429097")

dat <- read.csv(input_path, check.names = FALSE)
frozen <- read.csv(frozen_path, check.names = FALSE)
required <- c("Run.ID", seven_targeted, "Group", "Cohort")
if (!all(required %in% names(dat))) stop("Pooled input lacks required columns")
if (nrow(dat) != 1647L || anyDuplicated(dat$Run.ID)) stop("Expected 1,647 unique Run IDs")
if (!setequal(unique(dat$Group), c("Healthy", "Cancer"))) stop("Unexpected Group labels")
if (!setequal(unique(dat$Cohort), cohorts_expected)) stop("Unexpected cohort labels")
if (sum(dat$Group == "Healthy") != 745L || sum(dat$Group == "Cancer") != 902L) stop("Pooled group counts changed")
if (!setequal(frozen$KO, seven_targeted)) stop("Frozen pooled summary is not the seven targeted KOs")

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

# Confirm that the frozen q column is exactly BH across all seven targeted KO P values.
q_recomputed <- p.adjust(frozen$p_value, method = "BH")
if (!isTRUE(all.equal(q_recomputed, frozen$p_adj, tolerance = 1e-12, check.attributes = FALSE))) {
  stop("Frozen pooled q values do not match BH correction across seven targeted KOs")
}

results <- vector("list", nrow(ko_info))
for (i in seq_len(nrow(ko_info))) {
  meta <- ko_info[i, ]
  ko <- meta$KO
  healthy <- as.numeric(dat[dat$Group == "Healthy", ko])
  cancer <- as.numeric(dat[dat$Group == "Cancer", ko])
  if (any(!is.finite(c(healthy, cancer))) || any(c(healthy, cancer) < 0)) stop("Invalid KO values: ", ko)

  effect <- rank_biserial(healthy, cancer)
  effect_direct <- rank_biserial_direct(healthy, cancer)
  if (abs(effect - effect_direct) > 1e-10) stop("Rank-biserial cross-check failed: ", ko)
  wt <- suppressWarnings(wilcox.test(healthy, cancer, exact = FALSE))
  effect_from_w <- 2 * unname(wt$statistic) / (length(healthy) * length(cancer)) - 1
  if (!isTRUE(all.equal(effect, effect_from_w, tolerance = 1e-12))) stop("Wilcoxon effect cross-check failed: ", ko)
  ci <- bootstrap_ci(healthy, cancer, reps = 5000L, seed = 42L)

  frow <- frozen[frozen$KO == ko, , drop = FALSE]
  if (nrow(frow) != 1L) stop("Missing/duplicate frozen row: ", ko)
  comparisons <- c(
    Mean_Healthy = mean(healthy), Mean_Cancer = mean(cancer),
    Median_Healthy = median(healthy), Median_Cancer = median(cancer),
    p_value = wt$p.value
  )
  for (nm in names(comparisons)) {
    if (!isTRUE(all.equal(as.numeric(comparisons[[nm]]), as.numeric(frow[[nm]]), tolerance = 1e-10))) {
      stop("Frozen pooled mismatch: ", ko, " / ", nm)
    }
  }

  stats <- data.frame(
    KO = ko,
    Enzyme = frow$Enzyme,
    EC = frow$EC,
    Role = frow$Role,
    Analysis = "Pooled four CRC WGS cohorts; Healthy vs Cancer",
    Cohorts_pooled = paste(cohorts_expected, collapse = ";"),
    n_Healthy = length(healthy),
    n_CRC = length(cancer),
    zero_Healthy = sum(healthy == 0),
    zero_CRC = sum(cancer == 0),
    prevalence_Healthy_percent = 100 * mean(healthy > 0),
    prevalence_CRC_percent = 100 * mean(cancer > 0),
    overall_prevalence_percent = frow$Prev_Overall,
    passes_10_percent_overall_prevalence = frow$Pass_Prevalence,
    median_Healthy = median(healthy),
    median_CRC = median(cancer),
    mean_Healthy = mean(healthy),
    mean_CRC = mean(cancer),
    rank_biserial_Healthy_vs_CRC = effect,
    ci_95_low = ci[1],
    ci_95_high = ci[2],
    wilcox_p = wt$p.value,
    BH_q_across_7_targeted_KOs = frow$p_adj,
    Direction = ifelse(effect > 0, "Higher in Healthy", ifelse(effect < 0, "Higher in CRC", "No direction")),
    CI_status = ifelse(ci[1] > 0 | ci[2] < 0, "Excludes zero", "Includes zero"),
    bootstrap_reps = 5000L,
    bootstrap_seed = 42L,
    stringsAsFactors = FALSE
  )

  plot_data_path <- file.path(out_dir, paste0(meta$Slug, "_pooled_forest_plot_data.csv"))
  write.csv(stats, plot_data_path, row.names = FALSE, na = "")
  raw_data <- data.frame(
    Run.ID = dat$Run.ID,
    Cohort = dat$Cohort,
    Group = dat$Group,
    KO = ko,
    KO_relative_abundance = dat[[ko]],
    stringsAsFactors = FALSE
  )
  raw_data_path <- file.path(out_dir, paste0(meta$Slug, "_pooled_raw_values.csv"))
  write.csv(raw_data, raw_data_path, row.names = FALSE, na = "")

  point_fill <- if (isTRUE(stats$passes_10_percent_overall_prevalence)) {
    if (stats$Direction == "Higher in Healthy") COL_HEALTHY else COL_CANCER
  } else {
    "white"
  }
  point_outline <- if (stats$Direction == "Higher in Healthy") COL_HEALTHY else COL_CANCER

  p <- ggplot(stats, aes(x = rank_biserial_Healthy_vs_CRC, y = 1)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.55, ymax = 1.45, fill = "#F7F7F7") +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = 0.45, color = "#777777") +
    geom_segment(aes(x = ci_95_low, xend = ci_95_high, yend = 1), color = "#8F8F8F", linewidth = 0.75) +
    geom_point(shape = 21, size = 3.4, stroke = 0.9, fill = point_fill, color = point_outline) +
    annotate("text", x = -0.38, y = 1.77, label = "CRC higher", color = COL_CANCER,
             family = "Arial", fontface = "bold", size = 2.35) +
    annotate("text", x = 0.38, y = 1.77, label = "Healthy higher", color = COL_HEALTHY,
             family = "Arial", fontface = "bold", size = 2.35) +
    annotate("text", x = 0, y = 2.08, label = paste0(meta$Display, " (", ko, ")"),
             family = "Arial", fontface = ifelse(meta$Display %in% c("soxA", "soxG"), "bold.italic", "bold"),
             size = ifelse(meta$Display == "Sarcosine oxidase", 2.7, 3.0)) +
    scale_x_continuous(limits = c(-0.55, 0.55), breaks = c(-0.5, 0, 0.5), expand = c(0, 0)) +
    scale_y_continuous(
      breaks = 1, labels = "Pooled 4 cohorts\nH/CRC: 745/902",
      limits = c(0.45, 2.17), expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x = expression(paste("Rank-biserial effect size (", r[rb], ")")), y = NULL,
      caption = if (!isTRUE(stats$passes_10_percent_overall_prevalence))
        paste0("Open circle: overall prevalence ", format(stats$overall_prevalence_percent, trim = TRUE), "% (<10%)")
      else NULL
    ) +
    theme_nc(base_pt = 7.2) +
    theme(
      plot.margin = margin(6, 7, 5, 7, unit = "pt"),
      axis.text.y = element_text(size = 5.7, lineheight = 0.9, hjust = 1),
      axis.text.x = element_text(size = 6.0),
      axis.title.x = element_text(size = 6.2, margin = margin(t = 4)),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      plot.caption = element_text(size = 4.8, color = "#666666", hjust = 0.5, margin = margin(t = 3))
    )

  png_path <- file.path(out_dir, paste0(meta$Slug, "_pooled_forest.png"))
  pdf_path <- file.path(out_dir, paste0(meta$Slug, "_pooled_forest.pdf"))
  ggsave(png_path, p, width = 2.80, height = 1.55, units = "in", dpi = 600, device = ragg::agg_png, bg = "white")
  ggsave(pdf_path, p, width = 2.80, height = 1.55, units = "in", device = cairo_pdf, bg = "white")
  standardize_png(png_path)

  transport_path <- file.path(transport_dir, paste0("CRC_WGS_POOLED_", meta$Slug, "_forest.png"))
  copy_exact_visible(png_path, transport_path)
  copy_exact_visible(
    plot_data_path,
    file.path(data_copy_dir, paste0("CRC_WGS_POOLED_", meta$Slug, "_forest_plot_data.csv"))
  )
  copy_exact_visible(
    raw_data_path,
    file.path(data_copy_dir, paste0("CRC_WGS_POOLED_", meta$Slug, "_raw_values.csv"))
  )
  results[[i]] <- stats
}

combined <- bind_rows(results)
combined_path <- file.path(out_dir, "CRC_WGS_remaining_5KO_pooled_forest_plot_data_combined.csv")
write.csv(combined, combined_path, row.names = FALSE, na = "")
copy_exact_visible(
  combined_path,
  file.path(data_copy_dir, "CRC_WGS_remaining_5KO_pooled_forest_plot_data_combined.csv")
)

input_manifest <- data.frame(
  File = normalizePath(c(input_path, frozen_path)),
  Bytes = file.info(c(input_path, frozen_path))$size,
  SHA256 = vapply(c(input_path, frozen_path), digest, character(1), algo = "sha256", file = TRUE),
  stringsAsFactors = FALSE
)
write.table(input_manifest, file.path(out_dir, "INPUT_SHA256.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

readme <- c(
  "# Remaining five KO pooled Healthy-versus-CRC effect plots",
  "",
  "Generated: 2026-09-01",
  "",
  "These panels pool all four CRC WGS cohorts at the Run level: 745 Healthy and 902 CRC samples.",
  "They do not show separate cohort rows and are not meta-analysis estimates.",
  "",
  "- Effect: rank-biserial r_rb = P(Healthy > CRC) - P(Healthy < CRC).",
  "- Interval: 95% group-stratified percentile bootstrap CI; 5,000 resamples; seed 42.",
  "- P: two-sided pooled Wilcoxon rank-sum test.",
  "- q: Benjamini-Hochberg correction across all seven targeted sarcosine KOs.",
  "- Open points denote overall prevalence below 10%, matching the low-prevalence caution used in the seven-KO pooled lollipop panel.",
  "",
  "Each KO has two exact user-facing CSVs: one one-row plot/statistics table and one 1,647-row raw-value table.",
  "Means, medians, P values, q values, and prevalence were required to match the frozen pooled KO summary.",
  "Rank-biserial effects were independently cross-checked by rank-sum, direct pairwise signs, and the Wilcoxon W statistic.",
  "PNG copies were standardized to 8-bit sRGB at 600 dpi and are byte-identical to the analysis outputs."
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

message("Created five pooled KO panels: ", out_dir)
message("Copied pooled per-KO plot and raw-value CSVs: ", data_copy_dir)
