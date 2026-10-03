#!/usr/bin/env Rscript

# Re-render the seven saved four-cohort KO forest plots with an explicit,
# common visual key. No effect estimate, confidence interval, P value, or
# prevalence is recomputed; all plotted values come from the frozen plot-data
# CSV package exported on 2026-09-01.

suppressPackageStartupMessages({
  library(digest)
  library(ggplot2)
  library(ragg)
})

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", script_arg))
base_dir <- normalizePath(file.path(dirname(script_path), ".."))
analysis_root <- normalizePath(file.path(base_dir, ".."))

source_dir <- file.path(
  analysis_root, "사용데이터_모음",
  "CRC_WGS_7KO_4cohort_forest_source_data_26.09.01"
)
archive_dir <- file.path(base_dir, "results_integrated", "아카이브")
out_dir <- file.path(
  base_dir, "results_integrated",
  "CRC_WGS_7KO_4cohort_forests_EXPLICIT_LEGEND_26.09.01"
)
backup_dir <- file.path(out_dir, "pre_explicit_legend_archive_originals")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(backup_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(source_dir)) stop("Missing frozen source-data directory: ", source_dir)
if (!dir.exists(archive_dir)) stop("Missing archive directory: ", archive_dir)

COL_CRC <- "#CC3D3D"
COL_HEALTHY <- "#20A394"
COL_ZERO <- "#929292"
COHORTS <- c("PRJEB6070", "PRJEB10878", "PRJEB27928", "PRJNA429097")

ko_info <- data.frame(
  KO = c("K00301", "K00302", "K00303", "K00305", "K00306", "K00315", "K08688"),
  Display = c("Sarcosine oxidase", "soxA", "soxB", "soxG", "PIPOX", "DMGDH", "Creatinase"),
  SourceSlug = c(
    "K00301_Sarcosine_oxidase", "K00302_soxA", "K00303_soxB",
    "K00305_soxG", "K00306_PIPOX", "K00315_DMGDH", "K08688_Creatinase"
  ),
  ArchivePNG = c(
    "K00301_Sarcosine_oxidase_4cohort_forest.png",
    "K00302_soxA_4cohort_forest.png",
    "K00303_SoxB_4cohorts.png",
    "K00305_soxG_4cohort_forest.png",
    "K00306_PIPOX_4cohort_forest.png",
    "K00315_DMGDH_4cohort_forest.png",
    "Creatinase_K08688_cross_cohort_forest.png"
  ),
  PooledPrevalence = c(1.3, 52.4, 92.8, 9.0, 9.8, 3.2, 55.5),
  stringsAsFactors = FALSE
)

required_common <- c(
  "Cohort", "n_Healthy", "n_CRC", "zero_Healthy", "zero_CRC",
  "rank_biserial_Healthy_vs_CRC", "ci_95_low", "ci_95_high", "Direction"
)

standardize_png <- function(path, dpi = 600L) {
  srgb_profile <- "/System/Library/ColorSync/Profiles/sRGB Profile.icc"
  if (!file.exists(srgb_profile)) stop("System sRGB profile not found")
  tmp <- tempfile(fileext = ".png")
  result <- system2(
    "/usr/bin/sips",
    c(
      "-s", "dpiWidth", as.character(dpi),
      "-s", "dpiHeight", as.character(dpi),
      "-e", shQuote(srgb_profile), shQuote(path), "--out", shQuote(tmp)
    ),
    stdout = TRUE, stderr = TRUE
  )
  if (!file.exists(tmp) || file.info(tmp)$size <= 5000) {
    stop("PNG standardization failed: ", paste(result, collapse = "\n"))
  }
  if (!file.copy(tmp, path, overwrite = TRUE)) stop("Could not install standardized PNG")
  unlink(tmp)
}

copy_exact <- function(from, to) {
  if (!file.copy(from, to, overwrite = TRUE)) stop("Copy failed: ", to)
  if (!identical(
    digest(from, algo = "sha256", file = TRUE),
    digest(to, algo = "sha256", file = TRUE)
  )) stop("Hash mismatch after copy: ", to)
  Sys.chmod(to, "0644")
  invisible(system2("/usr/bin/xattr", c("-c", shQuote(to)), stdout = FALSE, stderr = FALSE))
  invisible(system2("/usr/bin/chflags", c("nohidden", shQuote(to)), stdout = FALSE, stderr = FALSE))
}

title_expression <- function(ko) {
  switch(
    ko,
    K00302 = expression(bolditalic(soxA)~bold("(K00302)")),
    K00303 = expression(bolditalic(soxB)~bold("(K00303)")),
    K00305 = expression(bolditalic(soxG)~bold("(K00305)")),
    K00301 = expression(bold("Sarcosine oxidase (K00301)")),
    K00306 = expression(bold("PIPOX (K00306)")),
    K00315 = expression(bold("DMGDH (K00315)")),
    K08688 = expression(bold("Creatinase (K08688)"))
  )
}

rendered <- character(0)
for (i in seq_len(nrow(ko_info))) {
  meta <- ko_info[i, , drop = FALSE]
  source_csv <- file.path(
    source_dir,
    paste0("CRC_WGS_4COHORT_", meta$SourceSlug, "_forest_plot_data.csv")
  )
  if (!file.exists(source_csv)) stop("Missing frozen plot data: ", source_csv)
  dat <- read.csv(source_csv, check.names = FALSE)
  if (!all(required_common %in% names(dat))) stop("Missing required columns: ", source_csv)
  if (nrow(dat) != 4L || !setequal(dat$Cohort, COHORTS)) {
    stop("Expected exactly four named cohort rows: ", source_csv)
  }
  dat <- dat[match(COHORTS, dat$Cohort), , drop = FALSE]
  if (any(!is.finite(dat$rank_biserial_Healthy_vs_CRC))) stop("Non-finite effect: ", source_csv)
  if (any(!is.finite(dat$ci_95_low)) || any(!is.finite(dat$ci_95_high))) {
    stop("Non-finite CI: ", source_csv)
  }

  dat$Estimable <- !(dat$zero_Healthy == dat$n_Healthy & dat$zero_CRC == dat$n_CRC)
  dat$DirectionCanonical <- ifelse(
    !dat$Estimable, "All zero",
    ifelse(dat$rank_biserial_Healthy_vs_CRC > 0, "Healthy higher",
           ifelse(dat$rank_biserial_Healthy_vs_CRC < 0, "CRC higher", "No direction"))
  )
  supplied_direction <- sub("Higher in Healthy", "Healthy higher", dat$Direction, fixed = TRUE)
  supplied_direction <- sub("Higher in CRC", "CRC higher", supplied_direction, fixed = TRUE)
  supplied_direction[supplied_direction %in% c("No variation", "No direction")] <- dat$DirectionCanonical[
    supplied_direction %in% c("No variation", "No direction")
  ]
  estimable_direction_rows <- dat$Estimable & dat$DirectionCanonical != "No direction"
  if (any(supplied_direction[estimable_direction_rows] != dat$DirectionCanonical[estimable_direction_rows])) {
    stop("Direction disagrees with effect sign: ", source_csv)
  }

  dat$PooledPrevalence <- meta$PooledPrevalence
  dat$PassPrevalence <- meta$PooledPrevalence >= 10
  dat$PointOutline <- ifelse(
    !dat$Estimable, COL_ZERO,
    ifelse(dat$DirectionCanonical == "Healthy higher", COL_HEALTHY, COL_CRC)
  )
  dat$PointFill <- ifelse(
    !dat$Estimable | !dat$PassPrevalence, "white", dat$PointOutline
  )
  dat$y <- match(dat$Cohort, rev(COHORTS))
  dat$CohortLabel <- sprintf(
    "%s\nH/CRC: %d/%d%s",
    dat$Cohort, dat$n_Healthy, dat$n_CRC,
    ifelse(dat$Estimable, "", " · all zero")
  )

  row_background <- data.frame(ymin = c(0.5, 2.5), ymax = c(1.5, 3.5))
  legend_text <- paste0(
    "Fill (same for all rows): filled = pooled prevalence ≥10%;\n",
    "open = pooled prevalence <10%\n",
    "Color: red = CRC higher; teal = Healthy higher\n",
    "Gray open = all zero / not estimable\n",
    "Whisker = 95% bootstrap CI; dashed line = no effect"
  )

  p <- ggplot(dat, aes(x = rank_biserial_Healthy_vs_CRC, y = y)) +
    geom_rect(
      data = row_background,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE, fill = "#F5F5F5", colour = NA
    ) +
    annotate(
      "segment", x = 0, xend = 0, y = 0.5, yend = 4.38,
      linewidth = 0.35, linetype = "22", colour = "#606060"
    ) +
    geom_segment(
      data = dat[dat$Estimable, , drop = FALSE],
      aes(x = ci_95_low, xend = ci_95_high, yend = y),
      linewidth = 0.35, colour = "#8F8F8F", lineend = "round"
    ) +
    geom_point(
      shape = 21, size = 3.2, stroke = 0.85,
      fill = dat$PointFill, colour = dat$PointOutline
    ) +
    geom_text(
      data = dat[!dat$Estimable, , drop = FALSE],
      aes(x = 0.035, label = "all zero"), hjust = 0,
      family = "Arial", size = 1.75, colour = "#808080"
    ) +
    annotate(
      "text", x = -0.38, y = 4.58, label = "CRC higher",
      family = "Arial", fontface = "bold", size = 2.35, colour = COL_CRC
    ) +
    annotate(
      "text", x = 0.38, y = 4.58, label = "Healthy higher",
      family = "Arial", fontface = "bold", size = 2.35, colour = COL_HEALTHY
    ) +
    scale_x_continuous(
      name = expression("Rank-biserial effect size ("*r[rb]*")"),
      breaks = c(-0.5, 0, 0.5), limits = c(-0.55, 0.55), expand = c(0, 0)
    ) +
    scale_y_continuous(
      breaks = dat$y, labels = dat$CohortLabel,
      limits = c(0.5, 4.72), expand = c(0, 0)
    ) +
    coord_cartesian(clip = "off") +
    labs(title = title_expression(meta$KO), y = NULL, caption = legend_text) +
    theme_classic(base_family = "Arial", base_size = 6) +
    theme(
      plot.margin = margin(5, 6, 5, 6, unit = "pt"),
      plot.title = element_text(size = 7.2, hjust = 0.5, margin = margin(b = 2)),
      axis.text.y = element_text(size = 5.25, lineheight = 0.9, hjust = 1, colour = "black"),
      axis.text.x = element_text(size = 5.8, colour = "black"),
      axis.title.x = element_text(size = 5.9, margin = margin(t = 3), colour = "black"),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.grid = element_blank(),
      plot.caption = element_text(
        size = 4.10, lineheight = 1.04, colour = "#444444",
        hjust = 0, margin = margin(t = 4)
      )
    )

  out_png <- file.path(out_dir, meta$ArchivePNG)
  out_pdf <- file.path(out_dir, sub("\\.png$", ".pdf", meta$ArchivePNG))
  ggsave(out_png, p, width = 2.80, height = 2.80, units = "in", dpi = 600,
         device = ragg::agg_png, bg = "white")
  ggsave(out_pdf, p, width = 2.80, height = 2.80, units = "in",
         device = cairo_pdf, bg = "white")
  standardize_png(out_png, 600L)

  archive_png <- file.path(archive_dir, meta$ArchivePNG)
  backup_png <- file.path(backup_dir, meta$ArchivePNG)
  if (file.exists(archive_png) && !file.exists(backup_png)) copy_exact(archive_png, backup_png)
  copy_exact(out_png, archive_png)
  rendered <- c(rendered, out_png, out_pdf, archive_png)
}

readme <- c(
  "# Seven CRC WGS KO four-cohort forest plots with explicit visual key",
  "",
  "Generated: 2026-09-01",
  "",
  "No statistic was recomputed. All plotted cohort effects and intervals were read from",
  "the frozen plot-data CSVs in CRC_WGS_7KO_4cohort_forest_source_data_26.09.01.",
  "",
  "Visual key printed inside every PNG/PDF:",
  "- Point fill (KO-wide): filled = overall pooled KO prevalence >=10%; open = <10%.",
  "- Point color: red = CRC higher; teal = Healthy higher; gray = all zero/not estimable.",
  "- Horizontal whisker = 95% group-stratified bootstrap CI; dashed vertical line = no effect.",
  "",
  "The fill state is KO-wide and therefore identical across the four cohort rows.",
  "It does not encode statistical significance.",
  "",
  "The pre-change archive PNGs are retained under pre_explicit_legend_archive_originals/.",
  "The revised PNGs were also copied over the seven historical archive filenames."
)
writeLines(readme, file.path(out_dir, "README.md"), useBytes = TRUE)

manifest_files <- unique(c(rendered, file.path(out_dir, "README.md")))
manifest <- data.frame(
  File = normalizePath(manifest_files),
  Bytes = file.info(manifest_files)$size,
  SHA256 = vapply(manifest_files, digest, character(1), algo = "sha256", file = TRUE),
  stringsAsFactors = FALSE
)
write.table(
  manifest, file.path(out_dir, "OUTPUT_SHA256.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)

message("Rendered and installed seven explicit-legend forest plots: ", out_dir)
