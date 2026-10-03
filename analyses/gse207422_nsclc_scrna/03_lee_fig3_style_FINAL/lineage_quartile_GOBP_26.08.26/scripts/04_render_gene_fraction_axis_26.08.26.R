#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

sha256 <- function(path) {
  sub(" .*", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
}

plotdata_file <- file.path(analysis_dir, "results", "tables", "plotdata_top_GO_BP_terms.csv")
expected_plotdata_sha <- "3cb9e159efba6856ccca8823f46a4b0995544fdb41599bafe149133a8ffb027f"
if (!file.exists(plotdata_file) || sha256(plotdata_file) != expected_plotdata_sha) {
  stop("Frozen GO BP plot-data identity mismatch")
}

output_dir <- file.path(analysis_dir, "results", "figures_publication", "DE_N_AXIS_26.08.26")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
plot_data <- fread(plotdata_file)
required <- c("lineage", "direction", "Term", "N", "DE", "q", "gene_fraction")
if (!all(required %in% names(plot_data)) || nrow(plot_data) != 40L) {
  stop("Unexpected frozen GO BP plot-data schema or dimensions")
}
if (any(abs(plot_data$gene_fraction - plot_data$DE / plot_data$N) > 1e-15)) {
  stop("gene_fraction is not exactly DE/N")
}
if (!setequal(plot_data$lineage, c("Epithelial", "CAF")) ||
    !setequal(plot_data$direction, c("Top 25%", "Bottom 25%"))) {
  stop("Unexpected lineage/direction values")
}

neglog_range <- range(-log10(pmax(plot_data$q, .Machine$double.xmin)), finite = TRUE)
size_range <- range(plot_data$DE, finite = TRUE)
COL_LOW_Q <- "#D5E5F2"
COL_HIGH_Q <- "#C43C3C"
COL_TEXT <- "#202020"

theme_nc <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = COL_TEXT),
      axis.text = element_text(size = 6, colour = COL_TEXT),
      axis.text.y = element_text(size = 5.8),
      axis.title = element_text(size = 7, colour = COL_TEXT),
      axis.line = element_line(linewidth = 0.35, colour = COL_TEXT),
      axis.ticks = element_line(linewidth = 0.35, colour = COL_TEXT),
      plot.title = element_text(size = 7.5, face = "bold", hjust = 0),
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6),
      plot.margin = margin(4, 5, 4, 4, "pt")
    )
}

make_panel <- function(lineage_name, direction_name) {
  z <- copy(plot_data[lineage == lineage_name & direction == direction_name])
  if (nrow(z) != 10L) stop("Expected ten frozen terms for ", lineage_name, " / ", direction_name)
  z[, Term_plot := factor(Term, levels = rev(Term))]
  direction_short <- if (identical(direction_name, "Top 25%")) "Q4" else "Q1"
  ggplot(z, aes(gene_fraction, Term_plot)) +
    geom_point(
      aes(size = DE, colour = -log10(pmax(q, .Machine$double.xmin))),
      alpha = 0.95
    ) +
    scale_x_continuous(
      limits = c(0, 0.80), breaks = seq(0, 0.8, by = 0.2),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_colour_gradient(
      low = COL_LOW_Q, high = COL_HIGH_Q, limits = neglog_range,
      name = expression(-log[10]("BH q"))
    ) +
    scale_size_continuous(limits = size_range, range = c(1.8, 5.4), name = "DE genes") +
    labs(
      title = paste0(lineage_name, " — ", direction_short),
      x = "DE genes / represented\nterm genes", y = NULL
    ) +
    theme_nc()
}

combined <- (
  make_panel("Epithelial", "Top 25%") |
    make_panel("Epithelial", "Bottom 25%")
) / (
  make_panel("CAF", "Top 25%") |
    make_panel("CAF", "Bottom 25%")
) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Lineage-restricted GO Biological Process enrichment",
    subtitle = "Top and bottom quartiles of the production/degradation ratio within each lineage",
    theme = theme(
      plot.title = element_text(family = "Arial", size = 8.5, face = "bold", colour = COL_TEXT),
      plot.subtitle = element_text(family = "Arial", size = 6.5, colour = "#555555")
    )
  ) & theme(legend.position = "bottom")

stem <- "Fig_Epithelial_CAF_ratio_Q4_vs_Q1_GO_BP_DE_over_N"
png_file <- file.path(output_dir, paste0(stem, ".png"))
pdf_file <- file.path(output_dir, paste0(stem, ".pdf"))
ggsave(png_file, combined, width = 7.20, height = 8.40, units = "in", dpi = 600,
       device = "png", bg = "white", limitsize = FALSE)
ggsave(pdf_file, combined, width = 7.20, height = 8.40, units = "in",
       device = cairo_pdf, bg = "white", limitsize = FALSE)
if (!all(file.exists(c(png_file, pdf_file)))) stop("Figure export failed")

srgb_profile <- "/System/Library/ColorSync/Profiles/sRGB Profile.icc"
if (nzchar(Sys.which("sips")) && file.exists(srgb_profile)) {
  tmp_png <- paste0(png_file, ".srgb.png")
  status <- system2("sips", c("-m", shQuote(srgb_profile), shQuote(png_file),
                              "--out", shQuote(tmp_png)), stdout = TRUE, stderr = TRUE)
  if (!file.exists(tmp_png)) stop("sRGB conversion failed: ", paste(status, collapse = "\n"))
  if (!file.rename(tmp_png, png_file)) stop("Could not install sRGB PNG")
}
invisible(system2("xattr", c("-c", shQuote(png_file)), stdout = TRUE, stderr = TRUE))
invisible(system2("chflags", c("nohidden", shQuote(png_file)), stdout = TRUE, stderr = TRUE))

manifest <- data.table(
  file = c(png_file, pdf_file),
  type = c("PNG", "PDF"),
  sha256 = vapply(c(png_file, pdf_file), sha256, character(1)),
  width_in = 7.20, height_in = 8.40,
  x_axis = "DE genes / represented term genes (DE/N)",
  source_plotdata_sha256 = expected_plotdata_sha
)
fwrite(manifest, file.path(output_dir, "output_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
message("Gene-fraction-axis GO BP figure complete")
print(manifest)
