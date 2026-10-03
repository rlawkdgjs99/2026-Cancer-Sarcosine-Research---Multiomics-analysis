#!/usr/bin/env Rscript

# Cross-cohort effect-size forest plot for creatinase (K08688).
# r_rb = P(Healthy > CRC) - P(Healthy < CRC); ties contribute zero.
# The four cohort estimates are not combined as a formal meta-analysis.

suppressPackageStartupMessages({
  library(digest)
  library(ggplot2)
  library(magick)
  library(ragg)
})

set.seed(42L, kind = "Mersenne-Twister", normal.kind = "Inversion",
         sample.kind = "Rejection")

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine script path.")
script_path <- normalizePath(sub("^--file=", "", script_arg))
base_dir <- normalizePath(file.path(dirname(script_path), ".."))
project_root <- normalizePath(file.path(base_dir, "..", ".."))
workspace_root <- normalizePath(file.path(base_dir, "..", "..", ".."))
source(file.path(project_root, "_shared", "theme_nc_26.08.18.R"))

ko_id <- "K08688"
ko_label <- "Creatinase"
cohorts <- c("PRJEB6070", "PRJEB10878", "PRJEB27928", "PRJNA429097")
groups <- c("Healthy", "Cancer")
boot_reps <- 5000L
boot_seed <- 42L

input_file <- file.path(
  base_dir, "results_integrated", "sarcosine",
  "sarcosine_KO_per_sample_pooled.csv"
)
out_dir <- file.path(
  base_dir, "results_integrated",
  "publication_style_plots_CREATINASE_K08688_FOREST_26.09.01"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_base <- file.path(out_dir, "Creatinase_K08688_cross_cohort_forest")
out_png <- paste0(out_base, ".png")
out_pdf <- paste0(out_base, ".pdf")
out_stats <- paste0(out_base, "_statistics.csv")
out_checksum <- paste0(out_base, "_input_sha256.tsv")
out_readme <- paste0(out_base, "_README.txt")
out_session <- paste0(out_base, "_sessionInfo.txt")
transport_png <- file.path(
  workspace_root, "Manuscript작업", "FigDesign&Manuscript",
  "PPT_INSERT_Creatinase_K08688_forest.png"
)

if (!file.exists(input_file)) stop("Missing input: ", input_file)
dat <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("Run.ID", ko_id, "Group", "Cohort")
missing_columns <- setdiff(required, names(dat))
if (length(missing_columns)) {
  stop("Missing columns: ", paste(missing_columns, collapse = ", "))
}
dat <- dat[, required]
value <- dat[[ko_id]]

if (!setequal(unique(dat$Cohort), cohorts)) stop("Unexpected cohort labels.")
if (!setequal(unique(dat$Group), groups)) stop("Unexpected group labels.")
if (anyNA(dat) || any(!is.finite(value))) stop("Missing/non-finite values.")
if (any(value < 0)) stop("Negative KO relative abundance detected.")
if (anyDuplicated(dat$Run.ID)) stop("Duplicated Run.ID values detected.")
sample_counts <- with(dat, table(Cohort, Group))
if (any(sample_counts < 3L)) stop("A cohort/group has fewer than 3 samples.")
if (sum(sample_counts) != nrow(dat)) stop("Sample count mismatch.")

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

results <- vector("list", length(cohorts))
for (i in seq_along(cohorts)) {
  cohort <- cohorts[i]
  z <- dat[dat$Cohort == cohort, , drop = FALSE]
  healthy <- z[[ko_id]][z$Group == "Healthy"]
  cancer <- z[[ko_id]][z$Group == "Cancer"]
  effect <- rank_biserial(healthy, cancer)

  if (abs(effect - rank_biserial_direct(healthy, cancer)) > 1e-10) {
    stop(cohort, ": rank-biserial cross-check failed.")
  }
  wt <- suppressWarnings(wilcox.test(healthy, cancer, exact = FALSE))
  effect_from_w <- 2 * unname(wt$statistic) /
    (length(healthy) * length(cancer)) - 1
  if (!isTRUE(all.equal(effect, effect_from_w, tolerance = 1e-12))) {
    stop(cohort, ": Wilcoxon effect-size cross-check failed.")
  }

  boot_effect <- replicate(
    boot_reps,
    rank_biserial(
      sample(healthy, length(healthy), replace = TRUE),
      sample(cancer, length(cancer), replace = TRUE)
    )
  )
  ci <- unname(quantile(boot_effect, c(0.025, 0.975), type = 7))

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
    bootstrap_reps = boot_reps,
    bootstrap_seed = boot_seed,
    stringsAsFactors = FALSE
  )
}

stats <- do.call(rbind, results)
stats$wilcox_BH_q_across_4_cohorts <- p.adjust(stats$wilcox_p, method = "BH")
stats$Direction <- ifelse(
  stats$rank_biserial_Healthy_vs_CRC > 0, "Healthy higher",
  ifelse(stats$rank_biserial_Healthy_vs_CRC < 0, "CRC higher", "No direction")
)
stats$CI_status <- ifelse(
  stats$ci_95_low > 0 | stats$ci_95_high < 0,
  "95% CI excludes 0", "95% CI overlaps 0"
)

if (nrow(stats) != 4L || anyDuplicated(stats$Cohort)) {
  stop("Expected one result per cohort.")
}
if (any(abs(stats$rank_biserial_Healthy_vs_CRC) > 1 + 1e-12)) {
  stop("Effect size outside [-1, 1].")
}
if (any(stats$ci_95_low > stats$ci_95_high)) stop("Reversed CI detected.")
write.csv(stats, out_stats, row.names = FALSE, quote = TRUE)

plot_stats <- stats
plot_stats$Cohort <- factor(plot_stats$Cohort, levels = rev(cohorts))
plot_stats$y <- as.numeric(plot_stats$Cohort)
plot_stats$Cohort_label <- sprintf(
  "%s\nH/CRC: %d/%d", as.character(plot_stats$Cohort),
  plot_stats$n_Healthy, plot_stats$n_CRC
)

write.table(
  data.frame(
    File = "results_integrated/sarcosine/sarcosine_KO_per_sample_pooled.csv",
    SHA256 = digest::digest(input_file, algo = "sha256", file = TRUE)
  ),
  out_checksum, sep = "\t", row.names = FALSE, quote = FALSE
)

row_background <- data.frame(ymin = c(0.5, 2.5), ymax = c(1.5, 3.5))
p <- ggplot(plot_stats, aes(x = rank_biserial_Healthy_vs_CRC, y = y)) +
  geom_rect(
    data = row_background,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE, fill = "#F5F5F5", colour = NA
  ) +
  annotate(
    "segment", x = 0, xend = 0, y = 0.5, yend = 4.38,
    linewidth = NC_AXIS_LW, linetype = "22", colour = "#606060"
  ) +
  geom_segment(
    aes(x = ci_95_low, xend = ci_95_high, yend = y),
    linewidth = NC_AXIS_LW, colour = "#333333", lineend = "round"
  ) +
  geom_point(
    aes(fill = Direction), shape = 21, size = 3.2,
    stroke = NC_AXIS_LW, colour = "#202020"
  ) +
  scale_fill_manual(
    values = c("Healthy higher" = COL_HEALTHY, "CRC higher" = COL_CANCER),
    guide = "none"
  ) +
  annotate(
    "text", x = -0.008, y = 4.92, label = ko_label,
    family = "Arial", fontface = "bold", size = mm_text(NC_STRIP_PT),
    colour = "black", hjust = 1
  ) +
  annotate(
    "text", x = 0.008, y = 4.92, label = paste0("(", ko_id, ")"),
    family = "Arial", fontface = "bold", size = mm_text(NC_STRIP_PT),
    colour = "black", hjust = 0
  ) +
  annotate(
    "text", x = -0.38, y = 4.57, label = "CRC higher",
    family = "Arial", fontface = "bold", size = mm_text(NC_ANNOT_PT),
    colour = COL_CANCER
  ) +
  annotate(
    "text", x = 0.38, y = 4.57, label = "Healthy higher",
    family = "Arial", fontface = "bold", size = mm_text(NC_ANNOT_PT),
    colour = COL_HEALTHY
  ) +
  scale_x_continuous(
    name = expression("Rank-biserial effect size ("*r[rb]*")"),
    breaks = c(-0.5, 0, 0.5), labels = c("-0.5", "0", "0.5"),
    limits = c(-0.55, 0.55), expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = plot_stats$y, labels = plot_stats$Cohort_label,
    limits = c(0.5, 5.08), expand = c(0, 0)
  ) +
  coord_cartesian(clip = "off") +
  theme_nc(base_pt = NC_TICK_PT) +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_text(
      family = "Arial", size = NC_TICK_PT, colour = "black",
      lineheight = 0.88, margin = margin(r = 3)
    ),
    axis.title.x = element_text(
      family = "Arial", size = NC_TITLE_PT, colour = "black",
      margin = margin(t = 3)
    ),
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    plot.margin = margin(t = 1, r = 3, b = 1, l = 4)
  )

width_in <- 2.80
height_in <- 2.10
dpi <- 600L
ggsave(
  out_png, p, width = width_in, height = height_in, dpi = dpi,
  bg = "white", device = ragg::agg_png
)
ggsave(
  out_pdf, p, width = width_in, height = height_in, bg = "white",
  device = grDevices::cairo_pdf
)

srgb_profile <- "/System/Library/ColorSync/Profiles/sRGB Profile.icc"
if (!file.exists(srgb_profile)) stop("System sRGB profile not found.")
standardize_png <- function(path, dpi) {
  img <- magick::image_read(path)
  img <- magick::image_background(img, "white", flatten = TRUE)
  img <- magick::image_convert(
    img, format = "png", colorspace = "sRGB", depth = 8, type = "TrueColor"
  )
  img <- magick::image_strip(img)
  magick::image_write(img, path = path, format = "png", depth = 8)
  before <- png::readPNG(path, native = FALSE, info = TRUE)
  tmp <- tempfile(fileext = ".png")
  sips_output <- system2(
    "/usr/bin/sips",
    args = c(
      "-s", "dpiWidth", as.character(dpi),
      "-s", "dpiHeight", as.character(dpi),
      "-e", shQuote(srgb_profile), shQuote(path), "--out", shQuote(tmp)
    ),
    stdout = TRUE, stderr = TRUE
  )
  if (!file.exists(tmp) || file.info(tmp)$size <= 5000) {
    stop("sips failed: ", paste(sips_output, collapse = "\n"))
  }
  if (!file.copy(tmp, path, overwrite = TRUE)) stop("Could not install standardized PNG.")
  unlink(tmp)
  after <- png::readPNG(path, native = FALSE, info = TRUE)
  if (!identical(dim(before), dim(after)) ||
      !identical(as.numeric(before), as.numeric(after))) {
    stop("Embedding sRGB/dpi metadata changed decoded PNG pixels.")
  }
}

standardize_png(out_png, dpi)
if (!file.copy(out_png, transport_png, overwrite = TRUE)) {
  stop("Could not create PowerPoint transport PNG.")
}
for (path in c(out_png, out_pdf, out_stats, out_checksum, transport_png)) {
  if (file.exists(path)) {
    invisible(system2("/usr/bin/xattr", c("-c", shQuote(path)), stdout = TRUE, stderr = TRUE))
    invisible(system2("/bin/chmod", c("644", shQuote(path)), stdout = TRUE, stderr = TRUE))
    invisible(system2("/usr/bin/chflags", c("nohidden", shQuote(path)), stdout = TRUE, stderr = TRUE))
  }
}
if (!identical(unname(tools::md5sum(out_png)),
               unname(tools::md5sum(transport_png)))) {
  stop("PowerPoint transport PNG is not byte-identical to analytical PNG.")
}

expected_files <- c(out_png, out_pdf, out_stats, out_checksum, transport_png)
if (!all(file.exists(expected_files))) stop("An output was not written.")
if (any(file.info(c(out_png, out_pdf))$size <= 5000)) {
  stop("Artwork unexpectedly small.")
}
png_dim <- dim(png::readPNG(out_png))
expected_dim <- c(round(height_in * dpi), round(width_in * dpi))
if (!identical(as.integer(png_dim[1:2]), as.integer(expected_dim))) {
  stop("Unexpected PNG dimensions: ", paste(png_dim[1:2], collapse = " x "))
}

writeLines(c(
  "Creatinase (K08688): four-cohort effect-size forest plot",
  "",
  "Input: results_integrated/sarcosine/sarcosine_KO_per_sample_pooled.csv",
  "Effect: r_rb = P(Healthy > CRC) - P(Healthy < CRC); ties contribute zero.",
  "Uncertainty: 5,000 group-stratified percentile bootstrap replicates; seed 42.",
  "Cohort-specific estimates only; no pooled meta-analytic estimate was fitted.",
  sprintf("Final footprint: %.2f x %.2f inches; PNG: %d dpi.",
          width_in, height_in, dpi),
  paste0("PowerPoint transport PNG: ", transport_png),
  "",
  "Reproduce:",
  "Rscript scripts/104_creatinase_K08688_cross_cohort_forest_26.09.01.R"
), out_readme)
capture.output(sessionInfo(), file = out_session)

cat("Input: ", nrow(dat), " samples\n", sep = "")
print(sample_counts)
print(stats[, c(
  "Cohort", "n_Healthy", "n_CRC", "rank_biserial_Healthy_vs_CRC",
  "ci_95_low", "ci_95_high", "wilcox_p", "wilcox_BH_q_across_4_cohorts"
)], row.names = FALSE)
cat("Wrote:\n", paste(c(expected_files, out_readme, out_session), collapse = "\n"), "\n")
