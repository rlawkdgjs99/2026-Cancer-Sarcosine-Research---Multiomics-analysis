#!/usr/bin/env Rscript

# Aesthetic-only rerender of the validated Lee-Figure-3-style compact UMAP.
# No scores, cell labels, UMAP coordinates, thresholds, or statistics are
# recomputed. The original publication render is intentionally preserved.

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
.libPaths(c(file.path(analysis_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

set.seed(260826)

score_file <- file.path(analysis_dir, "intermediate", "01_cell_paper_style_scores.rds")
cap_file <- file.path(analysis_dir, "results", "tables", "plotdata_01_UMAP_display_caps.csv")
output_dir <- file.path(
  analysis_dir, "results", "figures_publication", "UMAP_DARKER_26.08.26"
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(score_file) || !file.exists(cap_file)) {
  stop("Required validated score or display-cap input is missing")
}

scores <- as.data.table(readRDS(score_file))
caps <- fread(cap_file)
required <- c(
  "cell_id", "final_lineage", "analysis_eligible", "umap_1", "umap_2",
  "production_module_shifted", "degradation_module_shifted",
  "ratio_defined", "production_degradation_ratio", "ratio_group"
)
if (!all(required %in% names(scores))) stop("Validated score table lacks required columns")
if (nrow(scores) != 92053L || anyDuplicated(scores$cell_id)) {
  stop("Unexpected validated score-table dimensions or duplicate cell IDs")
}
if (scores[, sum(analysis_eligible)] != 90512L) {
  stop("Unexpected number of analysis-eligible cells")
}

COL_PROD <- "#C47B3B"
COL_DEG <- "#2E5F8A"
COL_RATIO <- "#6E5A9A"
COL_HIGH <- "#C43C3C"
COL_LOW <- "#202020"

lineage_order <- c(
  "Epithelial", "CAF", "B cell", "Plasma cell", "CD4 T cell",
  "CD8 T cell", "Cycling T cell", "NK cell", "Mast cell",
  "Neutrophil", "Monocyte", "Macrophage", "Conventional DC", "pDC"
)
lineage_cols <- setNames(hcl.colors(length(lineage_order), palette = "Dark 3"), lineage_order)

plot_cells <- copy(scores[analysis_eligible == TRUE])
plot_cells[, final_lineage := factor(final_lineage, levels = lineage_order)]
if (anyNA(plot_cells$final_lineage)) stop("Unexpected analysis-eligible lineage")

# Larger, fully opaque points prevent sub-pixel dilution when the 7.2-inch
# composite is shown or printed at its intended size.
p_lineage <- ggplot(plot_cells[sample(.N)], aes(umap_1, umap_2, colour = final_lineage)) +
  geom_point(size = 0.16, alpha = 1.00, stroke = 0) +
  scale_colour_manual(values = lineage_cols, drop = FALSE) +
  coord_equal() +
  labs(title = "Cell lineages", colour = NULL) +
  theme_void(base_family = "Arial", base_size = 7) +
  theme(
    plot.title = element_text(size = 7.5, face = "bold"),
    legend.text = element_text(size = 5.5),
    legend.key.height = unit(2.8, "mm"),
    legend.key.width = unit(3.4, "mm"),
    legend.position = "right",
    legend.box = "vertical",
    plot.margin = margin(3, 3, 3, 3, "pt")
  ) +
  guides(colour = guide_legend(
    override.aes = list(size = 2.1, alpha = 1, stroke = 0), ncol = 1
  ))

ratio_median <- median(
  scores[ratio_defined == TRUE]$production_degradation_ratio,
  na.rm = TRUE
)
scores[, ratio_group_plot := factor(
  ratio_group,
  levels = c("Low", "High", "At median", "Undefined denominator")
)]
highlow_cols <- c(
  "Low" = COL_LOW,
  "High" = COL_HIGH,
  "At median" = "#8C8C8C",
  "Undefined denominator" = "#BDBDBD"
)
p_highlow <- ggplot(scores[sample(.N)], aes(umap_1, umap_2, colour = ratio_group_plot)) +
  geom_point(size = 0.16, alpha = 1.00, stroke = 0) +
  scale_colour_manual(values = highlow_cols, drop = TRUE) +
  coord_equal() +
  labs(
    title = "Production/degradation ratio",
    subtitle = paste0("Finite-score median = ", format(ratio_median, digits = 4)),
    colour = NULL, x = NULL, y = NULL
  ) +
  theme_void(base_family = "Arial", base_size = 7) +
  theme(
    plot.title = element_text(size = 7.5, face = "bold"),
    plot.subtitle = element_text(size = 6, colour = "#555555"),
    legend.position = "right",
    legend.text = element_text(size = 6),
    legend.key.height = unit(3.0, "mm"),
    plot.margin = margin(3, 3, 3, 3, "pt")
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2.0, alpha = 1)))

umap_long <- rbindlist(list(
  scores[, .(
    cell_id, umap_1, umap_2, axis = "Production",
    value = production_module_shifted
  )],
  scores[, .(
    cell_id, umap_1, umap_2, axis = "Degradation",
    value = degradation_module_shifted
  )],
  scores[ratio_defined == TRUE, .(
    cell_id, umap_1, umap_2, axis = "Production/degradation ratio",
    value = production_degradation_ratio
  )]
))
umap_long[, axis := factor(
  axis,
  levels = c("Production", "Degradation", "Production/degradation ratio")
)]

expected_caps <- umap_long[, .(
  display_cap_99_5 = as.numeric(quantile(value, 0.995, na.rm = TRUE)),
  maximum = max(value, na.rm = TRUE),
  cells = .N
), by = axis]
caps[, axis := factor(axis, levels = levels(umap_long$axis))]
setorder(expected_caps, axis)
setorder(caps, axis)
if (!isTRUE(all.equal(expected_caps, caps, tolerance = 1e-12, check.attributes = FALSE))) {
  stop("Recomputed UMAP display caps differ from the frozen plot-data table")
}
umap_long <- merge(
  umap_long,
  caps[, .(axis, display_cap_99_5)],
  by = "axis", all.x = TRUE, sort = FALSE
)
umap_long[, value_display := pmin(value, display_cap_99_5)]

axis_gradients <- c(
  "Production" = COL_PROD,
  "Degradation" = COL_DEG,
  "Production/degradation ratio" = COL_RATIO
)
make_score_umap <- function(axis_name) {
  d <- copy(umap_long[axis == axis_name])
  # Draw low-score cells first and high-score cells last so saturated signal is
  # not hidden underneath the darker background cloud.
  setorder(d, value_display, cell_id)
  high_col <- axis_gradients[[axis_name]]
  legend_breaks <- pretty(c(0, max(d$value_display, na.rm = TRUE)), n = 3)
  legend_breaks <- legend_breaks[
    legend_breaks >= 0 & legend_breaks <= max(d$value_display, na.rm = TRUE)
  ]
  ggplot(d, aes(umap_1, umap_2, colour = value_display)) +
    geom_point(size = 0.15, alpha = 1.00, stroke = 0) +
    scale_colour_gradient(
      low = "#AFAFAF", high = high_col, name = "Score",
      breaks = legend_breaks, labels = label_number(accuracy = 0.01)
    ) +
    coord_equal() +
    labs(title = axis_name, x = NULL, y = NULL) +
    theme_void(base_family = "Arial", base_size = 7) +
    theme(
      plot.title = element_text(size = 7.5, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6),
      legend.key.width = unit(24, "mm"),
      legend.key.height = unit(1.8, "mm"),
      plot.margin = margin(3, 3, 3, 3, "pt")
    ) +
    guides(colour = guide_colourbar(
      title.position = "top",
      barwidth = unit(24, "mm"), barheight = unit(1.8, "mm")
    ))
}

p_prod <- make_score_umap("Production")
p_deg <- make_score_umap("Degradation")
p_ratio <- make_score_umap("Production/degradation ratio")
p_continuous <- p_prod + p_deg + p_ratio + plot_layout(ncol = 3)

p_compact <- (p_lineage | p_highlow) / p_continuous +
  plot_annotation(
    title = "Lee et al. Figure 3-style sarcosine functional-score analysis",
    subtitle = "Production: GNMT+DMGDH; Degradation: SARDH+PIPOX"
  ) &
  theme(
    plot.title = element_text(family = "Arial", size = 8, face = "bold"),
    plot.subtitle = element_text(family = "Arial", size = 6.5)
  )

png_file <- file.path(output_dir, "Fig_scRNA_Lee2024_Fig3_style_compact_DARKER.png")
pdf_file <- file.path(output_dir, "Fig_scRNA_Lee2024_Fig3_style_compact_DARKER.pdf")
ggsave(
  png_file, p_compact, width = 7.20, height = 6.20,
  units = "in", dpi = 600, device = "png", bg = "white"
)
ggsave(
  pdf_file, p_compact, width = 7.20, height = 6.20,
  units = "in", device = cairo_pdf, bg = "white"
)
if (!all(file.exists(c(png_file, pdf_file)))) stop("Darker compact UMAP write failed")

audit <- data.table(
  item = c(
    "validated_score_cells", "analysis_eligible_lineage_cells",
    "finite_ratio_cells", "lineage_point_size_mm", "lineage_point_alpha",
    "continuous_point_size_mm", "continuous_point_alpha",
    "continuous_low_colour", "score_caps_changed"
  ),
  value = c(
    nrow(scores), nrow(plot_cells), scores[, sum(ratio_defined)],
    0.16, 1.00, 0.15, 1.00, "#AFAFAF", FALSE
  )
)
fwrite(audit, file.path(output_dir, "DARKER_RENDER_AUDIT.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

# Author-requested second candidate: remove only the bottom continuous Ratio
# map, widen Production and Degradation, and use more saturated sequential
# palettes.  The top binary Ratio panel remains unchanged.
output_dir_2axis <- file.path(
  analysis_dir, "results", "figures_publication", "UMAP_DARKER_2AXIS_26.08.26"
)
dir.create(output_dir_2axis, recursive = TRUE, showWarnings = FALSE)

strong_palettes <- list(
  "Production" = c("#D1B49D", "#C47B3B", "#7A2D0B"),
  "Degradation" = c("#AFC2D2", "#2E5F8A", "#0B365A")
)
make_stronger_score_umap <- function(axis_name) {
  if (!axis_name %in% names(strong_palettes)) stop("Unsupported two-axis UMAP")
  d <- copy(umap_long[axis == axis_name])
  setorder(d, value_display, cell_id)
  legend_breaks <- pretty(c(0, max(d$value_display, na.rm = TRUE)), n = 3)
  legend_breaks <- legend_breaks[
    legend_breaks >= 0 & legend_breaks <= max(d$value_display, na.rm = TRUE)
  ]
  ggplot(d, aes(umap_1, umap_2, colour = value_display)) +
    geom_point(size = 0.17, alpha = 1.00, stroke = 0) +
    scale_colour_gradientn(
      colours = strong_palettes[[axis_name]],
      values = c(0, 0.25, 1),
      name = "Score",
      breaks = legend_breaks,
      labels = label_number(accuracy = 0.01)
    ) +
    coord_equal() +
    labs(title = axis_name, x = NULL, y = NULL) +
    theme_void(base_family = "Arial", base_size = 7) +
    theme(
      plot.title = element_text(size = 7.5, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6),
      legend.key.width = unit(31, "mm"),
      legend.key.height = unit(1.8, "mm"),
      plot.margin = margin(3, 3, 3, 3, "pt")
    ) +
    guides(colour = guide_colourbar(
      title.position = "top",
      barwidth = unit(31, "mm"), barheight = unit(1.8, "mm")
    ))
}

p_prod_stronger <- make_stronger_score_umap("Production")
p_deg_stronger <- make_stronger_score_umap("Degradation")
p_continuous_2axis <- p_prod_stronger + p_deg_stronger + plot_layout(ncol = 2)
p_compact_2axis <- (p_lineage | p_highlow) / p_continuous_2axis +
  plot_annotation(
    title = "Lee et al. Figure 3-style sarcosine functional-score analysis",
    subtitle = "Production: GNMT+DMGDH; Degradation: SARDH+PIPOX"
  ) &
  theme(
    plot.title = element_text(family = "Arial", size = 8, face = "bold"),
    plot.subtitle = element_text(family = "Arial", size = 6.5)
  )

png_file_2axis <- file.path(
  output_dir_2axis,
  "Fig_scRNA_Lee2024_Fig3_style_compact_DARKER_2AXIS.png"
)
pdf_file_2axis <- file.path(
  output_dir_2axis,
  "Fig_scRNA_Lee2024_Fig3_style_compact_DARKER_2AXIS.pdf"
)
ggsave(
  png_file_2axis, p_compact_2axis, width = 7.20, height = 6.20,
  units = "in", dpi = 600, device = "png", bg = "white"
)
ggsave(
  pdf_file_2axis, p_compact_2axis, width = 7.20, height = 6.20,
  units = "in", device = cairo_pdf, bg = "white"
)
if (!all(file.exists(c(png_file_2axis, pdf_file_2axis)))) {
  stop("Darker two-axis compact UMAP write failed")
}

audit_2axis <- data.table(
  item = c(
    "validated_score_cells", "analysis_eligible_lineage_cells",
    "finite_ratio_cells", "bottom_axes", "bottom_ratio_removed",
    "continuous_point_size_mm", "continuous_point_alpha",
    "production_palette", "degradation_palette", "gradient_positions",
    "score_values_changed", "score_caps_changed"
  ),
  value = c(
    nrow(scores), nrow(plot_cells), scores[, sum(ratio_defined)],
    "Production|Degradation", TRUE, 0.17, 1.00,
    paste(strong_palettes[["Production"]], collapse = "|"),
    paste(strong_palettes[["Degradation"]], collapse = "|"),
    "0|0.25|1", FALSE, FALSE
  )
)
fwrite(audit_2axis, file.path(output_dir_2axis, "DARKER_2AXIS_RENDER_AUDIT.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir_2axis, "sessionInfo.txt"))

message("Darker compact UMAP render completed: ", png_file)
message("Darker two-axis compact UMAP render completed: ", png_file_2axis)
