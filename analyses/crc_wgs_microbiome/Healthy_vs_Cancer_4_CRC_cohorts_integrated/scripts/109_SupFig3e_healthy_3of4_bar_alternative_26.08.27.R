#!/usr/bin/env Rscript

# Display-only alternative for Supplementary Figure 3e.
# No model is re-fit and no statistic is recomputed. The script reads the
# verified 8-species x 4-cohort table produced by script 108 and presents the
# same 32 effect estimates as small-multiple horizontal bars.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(scales)
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
if (length(file_arg) != 1) stop("Run with Rscript so the script path can be resolved.")
SCRIPT_PATH <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_PATH)
INTEGRATED_DIR <- dirname(SCRIPT_DIR)
CRC_ROOT <- dirname(INTEGRATED_DIR)
ANALYSIS_ROOT <- dirname(CRC_ROOT)

THEME_FILE <- file.path(ANALYSIS_ROOT, "_shared", "theme_nc_26.08.18.R")
if (!file.exists(THEME_FILE)) stop("Shared theme not found: ", THEME_FILE)
source(THEME_FILE)

INPUT_FILE <- file.path(
  INTEGRATED_DIR, "results_integrated", "SupFig3de_HEALTHY_ENRICHED_26.08.27",
  "healthy_enriched_3of4_effects_all_cohorts.csv"
)
EXPECTED_INPUT_SHA256 <-
  "5bef913c8b28a80580292e464365d5891681ee18a6db96d0c8aeab9731dc97cb"
OUTDIR <- file.path(
  INTEGRATED_DIR, "results_integrated", "SupFig3e_HEALTHY_3OF4_BAR_ALT_26.08.27"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

DPI <- 600
WIDTH_IN <- 3.372
HEIGHT_IN <- 2.090
SRGB_PROFILE <- "/System/Library/ColorSync/Profiles/sRGB Profile.icc"

if (!file.exists(INPUT_FILE)) stop("Verified input table not found: ", INPUT_FILE)
if (!file.exists(SRGB_PROFILE)) stop("System sRGB profile not found: ", SRGB_PROFILE)
input_sha <- digest::digest(
  INPUT_FILE, algo = "sha256", file = TRUE, serialize = FALSE
)
if (!identical(input_sha, EXPECTED_INPUT_SHA256)) {
  stop("Input SHA-256 mismatch. Expected ", EXPECTED_INPUT_SHA256, "; found ", input_sha)
}

d <- readr::read_csv(INPUT_FILE, show_col_types = FALSE, progress = FALSE)
required_columns <- c(
  "Cohort", "Species_full", "Species", "log2FC", "p_adj",
  "strict_healthy", "neglog10_q"
)
if (!all(required_columns %in% names(d))) {
  stop("Missing required columns: ", paste(setdiff(required_columns, names(d)), collapse = ", "))
}
if (nrow(d) != 32L || n_distinct(d$Species) != 8L || n_distinct(d$Cohort) != 4L) {
  stop("Expected a complete 8-species x 4-cohort table (32 rows).")
}
if (anyDuplicated(d[c("Cohort", "Species")])) stop("Duplicated cohort-species rows")
if (any(!is.finite(d$log2FC)) || any(!is.finite(d$p_adj)) ||
    any(d$p_adj < 0 | d$p_adj > 1)) {
  stop("Invalid log2FC or BH q values")
}
strict_per_species <- d %>% count(Species, wt = as.integer(strict_healthy))
if (any(strict_per_species$n != 3L)) {
  stop("Every displayed species must meet the strict rule in exactly 3 cohorts")
}

cohort_order <- c("PRJEB10878", "PRJEB27928", "PRJEB6070", "PRJNA429097")
if (!setequal(unique(d$Cohort), cohort_order)) stop("Unexpected cohort identifiers")

species_order <- d %>%
  group_by(Species) %>%
  summarise(mean_log2FC = mean(log2FC), .groups = "drop") %>%
  arrange(mean_log2FC, Species) %>%
  pull(Species)

species_plotmath <- c(
  Adlercreutzia_equolifaciens = "italic(A.~equolifaciens)",
  Anaerobutyricum_hallii = "italic(A.~hallii)",
  Anaerostipes_hadrus = "italic(A.~hadrus)",
  Blautia_stercoris = "italic(B.~stercoris)",
  Blautia_wexlerae = "italic(B.~wexlerae)",
  Clostridium_sp_AF34_13 = "italic(Clostridium)~'sp. AF34-13'",
  Eubacterium_ventriosum = "italic(E.~ventriosum)",
  Faecalibacillus_intestinalis = "italic(F.~intestinalis)"
)
if (!all(species_order %in% names(species_plotmath))) stop("Missing plotmath labels")

d <- d %>%
  mutate(
    Species_label = factor(Species, levels = rev(species_order)),
    Cohort_label = factor(
      Cohort, levels = cohort_order,
      labels = c("PRJEB\n10878", "PRJEB\n27928", "PRJEB\n6070", "PRJNA\n429097")
    )
  )

p <- ggplot(
  d,
  aes(
    x = log2FC, y = Species_label,
    fill = strict_healthy, colour = strict_healthy
  )
) +
  geom_vline(xintercept = 0, colour = "#343434", linewidth = 0.30) +
  geom_vline(
    xintercept = -1, colour = "#8A8A8A", linewidth = 0.22,
    linetype = "22"
  ) +
  geom_col(width = 0.66, linewidth = 0.35) +
  facet_grid(. ~ Cohort_label) +
  scale_fill_manual(
    values = c(`FALSE` = "#ECECEC", `TRUE` = COL_HEALTHY),
    guide = "none"
  ) +
  scale_colour_manual(
    values = c(`FALSE` = "#A8A8A8", `TRUE` = "#202020"),
    guide = "none"
  ) +
  scale_x_continuous(
    breaks = c(-4, -2, 0),
    limits = c(-4.5, 1.0),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(labels = function(x) parse(text = unname(species_plotmath[x]))) +
  labs(
    x = expression(log[2] * "FC (Cancer/Healthy)"), y = NULL
  ) +
  theme_minimal(base_family = "Arial", base_size = 6) +
  theme(
    strip.text = element_text(
      family = "Arial", face = "bold", size = 5.3,
      colour = "#111111", lineheight = 0.90,
      margin = margin(b = 1.0)
    ),
    strip.background = element_blank(),
    axis.text.x = element_text(family = "Arial", size = 4.7, colour = "#222222"),
    axis.text.y = element_text(
      family = "Arial", size = 5.35, colour = "#111111",
      margin = margin(r = 1.5)
    ),
    axis.title.x = element_text(family = "Arial", size = 5.2, margin = margin(t = 1.0)),
    panel.spacing.x = unit(2.2, "pt"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(1.5, 2.0, 1.0, 1.0)
  )

standardize_png <- function(path) {
  img <- magick::image_read(path)
  img <- magick::image_background(img, "white", flatten = TRUE)
  img <- magick::image_convert(img, format = "png", colorspace = "sRGB", depth = 8)
  img <- magick::image_strip(img)
  magick::image_write(img, path = path, format = "png", depth = 8)
  before <- png::readPNG(path, native = FALSE, info = TRUE)

  tmp <- tempfile(fileext = ".png")
  output <- system2(
    "/usr/bin/sips",
    args = c(
      "-s", "dpiWidth", as.character(DPI),
      "-s", "dpiHeight", as.character(DPI),
      "-e", shQuote(SRGB_PROFILE), shQuote(path), "--out", shQuote(tmp)
    ),
    stdout = TRUE, stderr = TRUE
  )
  if (!file.exists(tmp) || file.info(tmp)$size <= 5000) {
    stop("sips failed: ", paste(output, collapse = "\n"))
  }
  if (!file.copy(tmp, path, overwrite = TRUE)) stop("Could not install standardized PNG")
  unlink(tmp)
  after <- png::readPNG(path, native = FALSE, info = TRUE)
  if (!identical(dim(before), dim(after)) ||
      !identical(as.numeric(before), as.numeric(after))) {
    stop("PNG metadata insertion changed decoded pixels")
  }
  metadata <- system2(
    "/usr/bin/sips",
    args = c("-g", "profile", "-g", "dpiWidth", "-g", "dpiHeight", shQuote(path)),
    stdout = TRUE, stderr = TRUE
  )
  if (!any(grepl("profile: sRGB IEC61966-2.1", metadata, fixed = TRUE)) ||
      !any(grepl("dpiWidth: 600.000", metadata, fixed = TRUE)) ||
      !any(grepl("dpiHeight: 600.000", metadata, fixed = TRUE))) {
    stop("PNG metadata validation failed: ", paste(metadata, collapse = "\n"))
  }
}

png_file <- file.path(OUTDIR, "SupFig3e_healthy_enriched_3of4_bar_alternative.png")
grDevices::png(
  png_file,
  width = as.integer(round(WIDTH_IN * DPI)),
  height = as.integer(round(HEIGHT_IN * DPI)),
  units = "px", res = DPI, type = "cairo-png", bg = "white"
)
print(p)
grDevices::dev.off()
standardize_png(png_file)

ggsave(
  file.path(OUTDIR, "SupFig3e_healthy_enriched_3of4_bar_alternative.pdf"),
  p, width = WIDTH_IN, height = HEIGHT_IN, device = cairo_pdf, bg = "white"
)
readr::write_csv(d, file.path(OUTDIR, "bar_plot_source_data.csv"))
writeLines(
  c(
    "Supplementary Figure 3e bar-chart display alternative (2026-08-27)",
    "No statistic was recomputed; all 32 rows came from the verified script-108 output.",
    paste0("Input SHA-256: ", input_sha),
    "Teal bars with black outlines: BH q < 0.05 and log2FC(Cancer/Healthy) < -1.",
    "Grey bars: the strict criterion was not met.",
    paste0("PNG: ", round(WIDTH_IN * DPI), " x ", round(HEIGHT_IN * DPI),
           " px; 600 dpi; 8-bit RGB; embedded sRGB profile.")
  ),
  file.path(OUTDIR, "VALIDATION_REPORT.txt"), useBytes = TRUE
)
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

cat("Bar-chart alternative complete.\n", png_file, "\n", sep = "")
