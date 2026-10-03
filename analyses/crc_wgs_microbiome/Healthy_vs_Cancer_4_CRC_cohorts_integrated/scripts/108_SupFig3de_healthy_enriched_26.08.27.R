#!/usr/bin/env Rscript

# Supplementary Figure 3d/e: recurrent Healthy-enriched species in four CRC WGS cohorts
# Date: 2026-08-27
#
# This script does NOT re-fit any statistical model. It reads the four frozen
# per-cohort differential-abundance tables and mirrors the published rule:
# Cancer enriched = BH q < 0.05 and log2FC(Cancer/Healthy) > 1;
# Healthy enriched = BH q < 0.05 and log2FC(Cancer/Healthy) < -1.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggvenn)
  library(ggtext)
  library(scales)
  library(readr)
})

required_packages <- c("digest", "magick", "png")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Run this file with Rscript so its path can be resolved.")
SCRIPT_PATH <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_PATH)
INTEGRATED_DIR <- dirname(SCRIPT_DIR)
CRC_ROOT <- dirname(INTEGRATED_DIR)
ANALYSIS_ROOT <- dirname(CRC_ROOT)

THEME_FILE <- file.path(ANALYSIS_ROOT, "_shared", "theme_nc_26.08.18.R")
if (!file.exists(THEME_FILE)) stop("Shared theme not found: ", THEME_FILE)
source(THEME_FILE)

OUTDIR <- file.path(
  INTEGRATED_DIR, "results_integrated", "SupFig3de_HEALTHY_ENRICHED_26.08.27"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

P_ADJ_CUT <- 0.05
LFC_CUT <- -1
DPI <- 600
SRGB_PROFILE <- "/System/Library/ColorSync/Profiles/sRGB Profile.icc"
if (!file.exists(SRGB_PROFILE)) stop("System sRGB profile not found: ", SRGB_PROFILE)

INPUT_FILES <- c(
  PRJEB10878 = file.path(
    CRC_ROOT, "PRJEB10878_CRC", "results_bacteria", "diff_abundance_species.csv"
  ),
  PRJEB27928 = file.path(
    CRC_ROOT, "PRJEB27928_CRC", "results_bacteria", "diff_abundance_species.csv"
  ),
  PRJEB6070 = file.path(
    CRC_ROOT, "PRJEB6070_CRC_AdenomatousPolyps", "results_bacteria",
    "diff_abundance_species.csv"
  ),
  PRJNA429097 = file.path(
    CRC_ROOT, "PRJNA429097_CRC", "results_bacteria", "diff_abundance_species.csv"
  )
)

EXPECTED_SHA256 <- c(
  PRJEB10878 = "340b744449d63899f02d6830dc66fb572988971efc3ee996b06ff31bed00b207",
  PRJEB27928 = "4fc21ec81d4fa3249ecd01b6c7d5ac6f59c7818ba0881519f4fb1d1f2a385bce",
  PRJEB6070 = "5ab68bb9aaf894542b3b69e06867aab7db0689c2fbd29bcef41329f540ea23c7",
  PRJNA429097 = "85dfa142cbae6e7da0791dafbb4576b009163b4bbf1b0f2314c6da8877d9a2eb"
)
EXPECTED_ROWS <- c(
  PRJEB10878 = 557L, PRJEB27928 = 484L, PRJEB6070 = 378L, PRJNA429097 = 396L
)
EXPECTED_COLUMNS <- c(
  "Species_full", "Species", "Mean_Healthy", "Mean_Cancer",
  "log2FC", "p_value", "p_adj"
)

if (!all(file.exists(INPUT_FILES))) {
  stop("Missing input(s): ", paste(INPUT_FILES[!file.exists(INPUT_FILES)], collapse = "; "))
}

input_sha256 <- vapply(
  INPUT_FILES, digest::digest, character(1),
  algo = "sha256", file = TRUE, serialize = FALSE
)
if (!identical(unname(input_sha256), unname(EXPECTED_SHA256[names(input_sha256)]))) {
  mismatch <- names(input_sha256)[input_sha256 != EXPECTED_SHA256[names(input_sha256)]]
  stop("Input SHA-256 mismatch: ", paste(mismatch, collapse = ", "))
}

cohort_data <- lapply(names(INPUT_FILES), function(cohort) {
  x <- readr::read_csv(INPUT_FILES[[cohort]], show_col_types = FALSE, progress = FALSE)
  if (!identical(names(x), EXPECTED_COLUMNS)) {
    stop(cohort, ": unexpected columns: ", paste(names(x), collapse = ", "))
  }
  if (nrow(x) != EXPECTED_ROWS[[cohort]]) {
    stop(cohort, ": expected ", EXPECTED_ROWS[[cohort]], " rows, found ", nrow(x))
  }
  if (anyDuplicated(x$Species_full)) stop(cohort, ": duplicated Species_full values")
  if (any(is.na(x$Species_full)) || any(x$Species_full == "")) {
    stop(cohort, ": missing Species_full values")
  }
  if (any(!is.finite(x$log2FC[!is.na(x$log2FC)]))) stop(cohort, ": non-finite log2FC")
  if (any(x$p_adj < 0 | x$p_adj > 1, na.rm = TRUE)) stop(cohort, ": invalid p_adj")
  x
})
names(cohort_data) <- names(INPUT_FILES)

healthy_rows <- lapply(cohort_data, function(x) {
  x %>% filter(!is.na(p_adj), !is.na(log2FC), p_adj < P_ADJ_CUT, log2FC < LFC_CUT)
})
healthy_sets <- lapply(healthy_rows, `[[`, "Species_full")

set_sizes <- vapply(healthy_sets, length, integer(1))
union_species <- Reduce(union, healthy_sets)
intersection_all4 <- Reduce(intersect, healthy_sets)

membership <- tibble(Species_full = union_species)
for (cohort in names(cohort_data)) {
  dat <- cohort_data[[cohort]]
  idx <- match(membership$Species_full, dat$Species_full)
  membership[[paste0(cohort, "_Species")]] <- dat$Species[idx]
  membership[[paste0(cohort, "_log2FC")]] <- dat$log2FC[idx]
  membership[[paste0(cohort, "_p_adj")]] <- dat$p_adj[idx]
  membership[[paste0(cohort, "_strict")]] <- membership$Species_full %in% healthy_sets[[cohort]]
}
strict_cols <- paste0(names(cohort_data), "_strict")
membership$n_cohorts_strict <- rowSums(as.data.frame(membership[strict_cols]))
membership$Species <- apply(
  as.data.frame(membership[paste0(names(cohort_data), "_Species")]), 1,
  function(z) z[which(!is.na(z) & z != "")[1]]
)
membership$mean_log2FC <- rowMeans(
  as.data.frame(membership[paste0(names(cohort_data), "_log2FC")]), na.rm = TRUE
)
membership <- membership %>%
  arrange(desc(n_cohorts_strict), mean_log2FC, Species)

membership_distribution <- table(factor(membership$n_cohorts_strict, levels = 1:4))

# Frozen-result assertions: stop if anything upstream changed.
stopifnot(
  identical(set_sizes, c(PRJEB10878 = 6L, PRJEB27928 = 107L,
                         PRJEB6070 = 51L, PRJNA429097 = 41L)),
  length(union_species) == 163L,
  length(intersection_all4) == 0L,
  unname(membership_distribution) == c(129L, 26L, 8L, 0L)
)

targets <- membership %>% filter(n_cohorts_strict == 3L)
stopifnot(nrow(targets) == 8L)

all_long <- bind_rows(lapply(names(cohort_data), function(cohort) {
  cohort_data[[cohort]] %>% mutate(Cohort = cohort)
}))

target_long <- all_long %>%
  filter(Species_full %in% targets$Species_full) %>%
  select(Cohort, Species_full, Species, Mean_Healthy, Mean_Cancer,
         log2FC, p_value, p_adj) %>%
  mutate(
    strict_healthy = !is.na(p_adj) & !is.na(log2FC) &
      p_adj < P_ADJ_CUT & log2FC < LFC_CUT,
    neglog10_q = pmin(-log10(pmax(p_adj, .Machine$double.xmin)), 12)
  )
stopifnot(
  nrow(target_long) == 32L,
  !anyDuplicated(target_long[c("Cohort", "Species_full")]),
  all(table(target_long$Species_full) == 4L),
  all(vapply(split(target_long$strict_healthy, target_long$Species_full), sum, integer(1)) == 3L)
)

# Tables retain every traceable frozen value used in the two plots.
per_cohort_list <- bind_rows(lapply(names(healthy_rows), function(cohort) {
  healthy_rows[[cohort]] %>% mutate(Cohort = cohort, .before = 1)
}))
summary_table <- tibble(
  Metric = c(
    paste0("Strict Healthy-enriched species: ", names(set_sizes)),
    "Union across four cohorts", "Exactly one cohort", "Exactly two cohorts",
    "Exactly three cohorts", "All four cohorts"
  ),
  Value = c(
    unname(set_sizes), length(union_species),
    unname(membership_distribution[1]), unname(membership_distribution[2]),
    unname(membership_distribution[3]), unname(membership_distribution[4])
  )
)

readr::write_csv(per_cohort_list, file.path(OUTDIR, "healthy_enriched_per_cohort.csv"))
readr::write_csv(membership, file.path(OUTDIR, "healthy_enriched_membership_matrix.csv"))
readr::write_csv(targets, file.path(OUTDIR, "healthy_enriched_exactly_3of4_species.csv"))
readr::write_csv(target_long, file.path(OUTDIR, "healthy_enriched_3of4_effects_all_cohorts.csv"))
readr::write_csv(summary_table, file.path(OUTDIR, "analysis_summary.csv"))
readr::write_tsv(
  tibble(Cohort = names(INPUT_FILES), File = unname(INPUT_FILES), SHA256 = input_sha256),
  file.path(OUTDIR, "input_sha256.tsv")
)

# Render and then flatten/rewrite as standard 8-bit sRGB PNG for PowerPoint.
# The explicit sips pass is required because image_strip() removes both the
# embedded colour profile and the physical-resolution metadata.
standardize_png <- function(path) {
  img <- magick::image_read(path)
  img <- magick::image_background(img, "white", flatten = TRUE)
  img <- magick::image_convert(img, format = "png", colorspace = "sRGB", depth = 8)
  img <- magick::image_strip(img)
  magick::image_write(img, path = path, format = "png", depth = 8)

  decoded_before_metadata <- png::readPNG(path, native = FALSE, info = TRUE)
  metadata_tmp <- tempfile(fileext = ".png")
  sips_output <- system2(
    "/usr/bin/sips",
    args = c(
      "-s", "dpiWidth", as.character(DPI),
      "-s", "dpiHeight", as.character(DPI),
      "-e", shQuote(SRGB_PROFILE),
      shQuote(path), "--out", shQuote(metadata_tmp)
    ),
    stdout = TRUE, stderr = TRUE
  )
  if (!file.exists(metadata_tmp) || file.info(metadata_tmp)$size <= 5000) {
    stop("sips failed to create the PowerPoint-safe PNG: ", paste(sips_output, collapse = "\n"))
  }
  if (!file.copy(metadata_tmp, path, overwrite = TRUE)) {
    stop("Could not replace PNG after embedding resolution/profile metadata: ", path)
  }
  unlink(metadata_tmp)

  decoded <- png::readPNG(path, native = FALSE, info = TRUE)
  if (length(dim(decoded)) != 3L || dim(decoded)[3] != 3L) {
    stop("PNG is not standard RGB after conversion: ", path)
  }
  # readPNG(info = TRUE) attaches metadata attributes, so compare only the
  # decoded RGB array dimensions and values.
  if (!identical(dim(decoded_before_metadata), dim(decoded)) ||
      !identical(as.numeric(decoded_before_metadata), as.numeric(decoded))) {
    stop("Embedding PNG metadata changed decoded scientific pixels: ", path)
  }

  sips_metadata <- system2(
    "/usr/bin/sips",
    args = c(
      "-g", "profile", "-g", "dpiWidth", "-g", "dpiHeight", shQuote(path)
    ),
    stdout = TRUE, stderr = TRUE
  )
  if (!any(grepl("profile: sRGB IEC61966-2.1", sips_metadata, fixed = TRUE)) ||
      !any(grepl("dpiWidth: 600.000", sips_metadata, fixed = TRUE)) ||
      !any(grepl("dpiHeight: 600.000", sips_metadata, fixed = TRUE))) {
    stop("PNG profile/resolution metadata validation failed: ",
         paste(sips_metadata, collapse = "\n"))
  }
  invisible(path)
}

save_final_png <- function(plot, path, width_in, height_in) {
  px_w <- as.integer(round(width_in * DPI))
  px_h <- as.integer(round(height_in * DPI))
  grDevices::png(
    filename = path, width = px_w, height = px_h, units = "px", res = DPI,
    type = "cairo-png", bg = "white"
  )
  print(plot)
  grDevices::dev.off()
  standardize_png(path)
  info <- magick::image_info(magick::image_read(path))
  stopifnot(
    file.exists(path), file.info(path)$size > 5000,
    info$width == px_w, info$height == px_h, info$colorspace == "sRGB"
  )
  invisible(path)
}

# Supplementary Figure 3d: exact strict sets, including the true central zero.
venn_colors <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")
p_venn <- ggvenn::ggvenn(
  healthy_sets,
  show_percentage = FALSE, show_stats = "c", show_set_totals = "none",
  fill_color = venn_colors, fill_alpha = 0.34,
  stroke_color = "#252525", stroke_size = 0.42,
  set_name_color = "#111111", set_name_size = 2.10,
  text_color = "#111111", text_size = 1.95, padding = 0.11
) +
  labs(
    title = "Healthy-enriched species",
    subtitle = expression(italic(q) < 0.05 ~ "and" ~ log[2] * "FC" < -1)
  ) +
  theme_void(base_family = "Arial") +
  theme(
    plot.title = element_text(
      family = "Arial", face = "bold", size = 8.5,
      hjust = 0.5, margin = margin(b = 0.6)
    ),
    plot.subtitle = element_text(
      family = "Arial", size = 5.7, colour = "#4A4A4A",
      hjust = 0.5, margin = margin(b = 0.5)
    ),
    plot.margin = margin(1.5, 2.5, 1.5, 2.5)
  )

venn_png <- file.path(OUTDIR, "SupFig3d_healthy_enriched_species_venn.png")
save_final_png(p_venn, venn_png, width_in = 2.852, height_in = 2.561)
ggsave(
  file.path(OUTDIR, "SupFig3d_healthy_enriched_species_venn.pdf"),
  p_venn, width = 2.852, height = 2.561, device = cairo_pdf, bg = "white"
)

# Supplementary Figure 3e: eight species recurring in exactly three cohorts.
species_order <- target_long %>%
  group_by(Species) %>%
  summarise(mean_log2FC = mean(log2FC, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_log2FC, Species) %>%
  pull(Species)

# Use plotmath expressions rather than markdown. This is robust to the ggplot2/
# ggtext version combination on this workstation and preserves scientific italics.
species_plotmath <- c(
  Adlercreutzia_equolifaciens = "italic(Adlercreutzia~equolifaciens)",
  Anaerobutyricum_hallii = "italic(Anaerobutyricum~hallii)",
  Anaerostipes_hadrus = "italic(Anaerostipes~hadrus)",
  Blautia_stercoris = "italic(Blautia~stercoris)",
  Blautia_wexlerae = "italic(Blautia~wexlerae)",
  Clostridium_sp_AF34_13 = "italic(Clostridium)~'sp. AF34-13'",
  Eubacterium_ventriosum = "italic(Eubacterium~ventriosum)",
  Faecalibacillus_intestinalis = "italic(Faecalibacillus~intestinalis)"
)
stopifnot(all(species_order %in% names(species_plotmath)))

target_long <- target_long %>%
  mutate(
    Species_label = factor(Species, levels = rev(species_order)),
    Cohort_label = factor(
      Cohort, levels = names(INPUT_FILES),
      labels = c("PRJEB\n10878", "PRJEB\n27928", "PRJEB\n6070", "PRJNA\n429097")
    )
  )

p_matrix <- ggplot(target_long, aes(x = Cohort_label, y = Species_label)) +
  geom_hline(
    yintercept = seq(1.5, 7.5, by = 1), linewidth = 0.18, colour = "#EEEEEE"
  ) +
  geom_vline(
    xintercept = seq(1.5, 3.5, by = 1), linewidth = 0.18, colour = "#EEEEEE"
  ) +
  geom_point(
    aes(size = neglog10_q, fill = log2FC),
    shape = 21, colour = "#B0B0B0", stroke = 0.22
  ) +
  geom_point(
    data = filter(target_long, strict_healthy),
    aes(size = neglog10_q, fill = log2FC),
    shape = 21, colour = "#202020", stroke = 0.42
  ) +
  scale_fill_gradient2(
    low = COL_HEALTHY, mid = "white", high = COL_CANCER,
    midpoint = 0, limits = c(-4.5, 1.0), oob = scales::squish,
    breaks = c(-4, -2, 0, 1), name = expression(log[2] * "FC")
  ) +
  scale_size_continuous(
    range = c(1.4, 4.0), limits = c(0, 12), breaks = c(2, 6, 10),
    name = expression(-log[10] * "(" * italic(q) * ")")
  ) +
  scale_y_discrete(
    labels = function(x) parse(text = unname(species_plotmath[x]))
  ) +
  labs(
    title = "Recurrent Healthy-enriched species",
    subtitle = "Eight species meeting the strict criterion in 3 of 4 cohorts",
    x = NULL, y = NULL
  ) +
  guides(
    fill = guide_colorbar(
      title.position = "top", title.hjust = 0.5,
      barwidth = unit(40, "pt"), barheight = unit(4.2, "pt"), order = 1
    ),
    size = guide_legend(
      title.position = "top", title.hjust = 0.5,
      direction = "horizontal", nrow = 1, order = 2,
      override.aes = list(fill = "#808080", colour = "#202020", stroke = 0.25)
    )
  ) +
  theme_minimal(base_family = "Arial", base_size = 6) +
  theme(
    plot.title = element_text(
      family = "Arial", face = "bold", size = 8.0,
      hjust = 0, margin = margin(b = 0.3)
    ),
    plot.subtitle = element_text(
      family = "Arial", size = 5.5, colour = "#4A4A4A",
      hjust = 0, margin = margin(b = 1.8)
    ),
    axis.text.x = element_text(
      family = "Arial", size = 5.5, colour = "#111111",
      lineheight = 0.90, margin = margin(t = 1.2)
    ),
    axis.text.y = element_text(
      family = "Arial", size = 5.7, colour = "#111111",
      margin = margin(r = 2.0)
    ),
    panel.grid = element_blank(),
    legend.position = "bottom", legend.box = "horizontal",
    legend.box.spacing = unit(0.5, "pt"), legend.spacing.x = unit(1.0, "pt"),
    legend.margin = margin(t = 0.5, r = 0, b = 0, l = 0),
    legend.title = element_text(family = "Arial", size = 5.3),
    legend.text = element_text(family = "Arial", size = 5.1),
    plot.margin = margin(2, 2.5, 0.5, 1.5)
  )

matrix_png <- file.path(OUTDIR, "SupFig3e_healthy_enriched_3of4_matrix.png")
save_final_png(p_matrix, matrix_png, width_in = 3.372, height_in = 2.090)
ggsave(
  file.path(OUTDIR, "SupFig3e_healthy_enriched_3of4_matrix.pdf"),
  p_matrix, width = 3.372, height = 2.090, device = cairo_pdf, bg = "white"
)

report_species <- gsub("_", " ", targets$Species, fixed = TRUE)
report_species[targets$Species == "Clostridium_sp_AF34_13"] <-
  "Clostridium sp. AF34-13"

report_lines <- c(
  "Supplementary Figure 3d/e Healthy-enriched analysis (2026-08-27)",
  "", "Direction convention: log2FC = Cancer / Healthy.",
  "Strict criterion: BH q < 0.05 and log2FC < -1.",
  "No statistical model was re-fit; values were carried from frozen per-cohort outputs.",
  "", paste0(names(set_sizes), ": ", unname(set_sizes), " strict species"),
  paste0("Union: ", length(union_species)),
  paste0("Exactly 1 cohort: ", membership_distribution[[1]]),
  paste0("Exactly 2 cohorts: ", membership_distribution[[2]]),
  paste0("Exactly 3 cohorts: ", membership_distribution[[3]]),
  paste0("All 4 cohorts: ", membership_distribution[[4]]),
  "", "Exactly 3-of-4 species:",
  paste0("- ", report_species),
  "", "PNG validation:",
  paste0("- SupFig3d: ", round(2.852 * DPI), " x ", round(2.561 * DPI),
         " px, 600 dpi, flattened 8-bit sRGB PNG"),
  paste0("- SupFig3e: ", round(3.372 * DPI), " x ", round(2.090 * DPI),
         " px, 600 dpi, flattened 8-bit sRGB PNG")
)
writeLines(report_lines, file.path(OUTDIR, "VALIDATION_REPORT.txt"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))

manifest_files <- list.files(OUTDIR, full.names = TRUE)
manifest_files <- manifest_files[basename(manifest_files) != "output_sha256.tsv"]
output_sha <- vapply(
  manifest_files, digest::digest, character(1),
  algo = "sha256", file = TRUE, serialize = FALSE
)
readr::write_tsv(
  tibble(File = basename(manifest_files), SHA256 = output_sha),
  file.path(OUTDIR, "output_sha256.tsv")
)

cat("Analysis complete.\n")
print(summary_table)
cat("\n3-of-4 species:\n")
print(targets %>% select(Species, n_cohorts_strict, mean_log2FC))
cat("\nOutputs:\n", OUTDIR, "\n", sep = "")
