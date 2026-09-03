#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
analysis_root <- normalizePath(file.path(analysis_dir, ".."), mustWork = TRUE)
data_root <- normalizePath(file.path(analysis_root, ".."), mustWork = TRUE)
prior_dir <- file.path(analysis_root, "02_lineage_reannotation")
.libPaths(c(file.path(analysis_dir, "R_libs"), file.path(prior_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

options(stringsAsFactors = FALSE, warn = 1, future.globals.maxSize = 32 * 1024^3)
set.seed(260825)

table_dir <- file.path(analysis_dir, "results", "tables")
intermediate_dir <- file.path(analysis_dir, "intermediate")
pub_dir <- file.path(analysis_dir, "results", "figures_publication")
diag_dir <- file.path(analysis_dir, "results", "figures_diagnostic")
log_dir <- file.path(analysis_dir, "logs")
for (d in c(table_dir, intermediate_dir, pub_dir, diag_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

plan_sha <- sub(" .*", "", system2(
  "shasum", c("-a", "256", file.path(analysis_dir, "ANALYSIS_PLAN_FROZEN.md")), stdout = TRUE
))
if (!identical(plan_sha, "6d057689efc65f0b262d70ce12fa471aabac51028de70bfdb62740fa2c550c2f")) {
  stop("Frozen plan hash mismatch")
}

counts_file <- file.path(prior_dir, "intermediate", "01_full_counts_mt20_qc.rds")
lineage_file <- file.path(prior_dir, "results", "tables", "09_final_cell_lineages_FROZEN.csv")
if (!all(file.exists(c(counts_file, lineage_file)))) stop("Required verified input is missing")
sha256 <- function(path) sub(" .*", "", system2("shasum", c("-a", "256", path), stdout = TRUE))
expected_input_sha <- c(
  counts = "7552ba882c306fa68380390a8172e0fc96af2c7be7cb3f5778345204e9fd61b2",
  lineages = "97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1"
)
observed_input_sha <- c(counts = sha256(counts_file), lineages = sha256(lineage_file))
if (!identical(observed_input_sha, expected_input_sha)) stop("Frozen input SHA-256 mismatch")

COL_PROD <- "#C47B3B"
COL_DEG <- "#2E5F8A"
COL_RATIO <- "#6E5A9A"
COL_TEXT <- "#202020"
COL_GREY <- "#D9D9D9"
COL_HIGH <- "#C43C3C"
COL_LOW <- "#303030"

lineage_order <- c(
  "Epithelial", "CAF", "B cell", "Plasma cell", "CD4 T cell",
  "CD8 T cell", "Cycling T cell", "NK cell", "Mast cell",
  "Neutrophil", "Monocyte", "Macrophage", "Conventional DC", "pDC"
)
lineage_cols <- setNames(hcl.colors(length(lineage_order), palette = "Dark 3"), lineage_order)

theme_nc <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = COL_TEXT),
      axis.text = element_text(size = 6, colour = COL_TEXT),
      axis.title = element_text(size = 7, colour = COL_TEXT),
      axis.line = element_line(linewidth = 0.35, colour = COL_TEXT),
      axis.ticks = element_line(linewidth = 0.35, colour = COL_TEXT),
      axis.ticks.length = unit(1.3, "mm"),
      plot.title = element_text(size = 7.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6, hjust = 0, colour = "#555555"),
      strip.background = element_blank(),
      strip.text = element_text(size = 7, face = "bold", colour = COL_TEXT),
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6),
      plot.margin = margin(3, 3, 3, 3, unit = "pt")
    )
}

save_plot <- function(stem, plot, width, height, directory = pub_dir) {
  png_file <- file.path(directory, paste0(stem, ".png"))
  pdf_file <- file.path(directory, paste0(stem, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, units = "in", dpi = 600,
         device = "png", bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, units = "in",
         device = cairo_pdf, bg = "white")
  if (!file.exists(png_file) || !file.exists(pdf_file)) stop("Figure write failed: ", stem)
}

message("Loading verified sparse counts")
counts <- readRDS(counts_file)
lin <- fread(lineage_file)
if (!inherits(counts, "dgCMatrix")) stop("Verified count object is not a dgCMatrix")
if (anyDuplicated(colnames(counts)) || anyDuplicated(lin$cell_id)) {
  stop("Duplicate cell IDs in count matrix or frozen lineage table")
}
count_only <- setdiff(colnames(counts), lin$cell_id)
lineage_only <- setdiff(lin$cell_id, colnames(counts))
if (length(count_only) || length(lineage_only)) {
  stop(
    "Count/lineage cell-ID sets differ: count-only=", length(count_only),
    "; lineage-only=", length(lineage_only)
  )
}
input_alignment_audit <- data.table(
  metric = c(
    "count_cells", "lineage_rows", "count_unique_cell_ids",
    "lineage_unique_cell_ids", "count_only_ids", "lineage_only_ids",
    "positions_identical_before_reorder", "positions_identical_after_reorder"
  ),
  value = c(
    ncol(counts), nrow(lin), uniqueN(colnames(counts)), uniqueN(lin$cell_id),
    length(count_only), length(lineage_only),
    sum(colnames(counts) == lin$cell_id), ncol(counts)
  )
)
# The frozen lineage CSV is grouped for audit readability, whereas the matrix
# retains its original cell order.  Reorder only after proving exact set equality.
lin <- lin[match(colnames(counts), cell_id)]
if (!identical(colnames(counts), lin$cell_id)) stop("Explicit lineage reordering failed")
fwrite(input_alignment_audit, file.path(table_dir, "00_input_cell_alignment_audit.csv"))
if (nrow(counts) != 24292L || ncol(counts) != 92053L) stop("Unexpected count dimensions")
targets <- c("GNMT", "DMGDH", "SARDH", "PIPOX")
if (!all(targets %in% rownames(counts)) || anyDuplicated(rownames(counts))) {
  stop("Target-gene presence or gene uniqueness check failed")
}

md <- as.data.frame(lin[, .(
  cell_id, Sample, Patient, Resource, Pathologic.Response,
  final_lineage, analysis_eligible, umap_1, umap_2
)])
rownames(md) <- md$cell_id
md$cell_id <- NULL

message("Building and LogNormalizing the full RNA assay")
obj <- CreateSeuratObject(counts = counts, meta.data = md, min.cells = 0, min.features = 0)
rm(counts, md)
invisible(gc())
obj <- NormalizeData(
  obj, assay = "RNA", normalization.method = "LogNormalize",
  scale.factor = 10000, margin = 1, verbose = TRUE
)
DefaultAssay(obj) <- "RNA"

message("Running the paper-style Seurat AddModuleScore call")
obj <- AddModuleScore(
  obj,
  features = list(
    Production = c("GNMT", "DMGDH"),
    Degradation = c("SARDH", "PIPOX")
  ),
  pool = rownames(obj), nbin = 24, ctrl = 100, k = FALSE,
  assay = "RNA", name = "SarcosineModule", seed = 260825,
  search = FALSE, slot = "data"
)
score_cols <- grep("^SarcosineModule[12]$", colnames(obj[[]]), value = TRUE)
if (!identical(score_cols, c("SarcosineModule1", "SarcosineModule2"))) {
  stop("Unexpected AddModuleScore output columns: ", paste(score_cols, collapse = ", "))
}

norm_targets <- LayerData(obj, assay = "RNA", layer = "data")[targets, , drop = FALSE]
score <- as.data.table(obj[[]], keep.rownames = "cell_id")
if (!identical(score$cell_id, Cells(obj))) stop("Seurat metadata/cell order mismatch")
score[, `:=`(
  production_module_raw = SarcosineModule1,
  degradation_module_raw = SarcosineModule2,
  GNMT_log1pCP10k = as.numeric(norm_targets["GNMT", ]),
  DMGDH_log1pCP10k = as.numeric(norm_targets["DMGDH", ]),
  SARDH_log1pCP10k = as.numeric(norm_targets["SARDH", ]),
  PIPOX_log1pCP10k = as.numeric(norm_targets["PIPOX", ])
)]
score[, `:=`(
  production_direct_mean = (GNMT_log1pCP10k + DMGDH_log1pCP10k) / 2,
  degradation_direct_mean = (SARDH_log1pCP10k + PIPOX_log1pCP10k) / 2
)]

prod_min <- min(score$production_module_raw)
deg_min <- min(score$degradation_module_raw)
score[, production_module_shifted := production_module_raw - prod_min]
score[, degradation_module_shifted := degradation_module_raw - deg_min]
if (min(score$production_module_shifted) != 0 || min(score$degradation_module_shifted) != 0) {
  stop("Minimum-to-zero transformation failed")
}
score[, ratio_defined := degradation_module_shifted > 0]
score[, production_degradation_ratio := fifelse(
  ratio_defined,
  production_module_shifted / degradation_module_shifted,
  NA_real_
)]
if (score[ratio_defined == TRUE, any(!is.finite(production_degradation_ratio))]) {
  stop("Finite denominator produced non-finite ratio")
}
ratio_median <- median(score$production_degradation_ratio, na.rm = TRUE)
score[, ratio_group := fifelse(
  !ratio_defined, "Undefined denominator",
  fifelse(production_degradation_ratio > ratio_median, "High",
          fifelse(production_degradation_ratio < ratio_median, "Low", "At median"))
)]

ratio_audit <- data.table(
  metric = c(
    "production_raw_min", "degradation_raw_min", "production_shifted_min",
    "degradation_shifted_min", "undefined_denominator_cells", "finite_ratio_cells",
    "finite_ratio_median", "high_cells", "low_cells", "at_median_cells"
  ),
  value = c(
    prod_min, deg_min, min(score$production_module_shifted),
    min(score$degradation_module_shifted), sum(!score$ratio_defined),
    sum(score$ratio_defined), ratio_median,
    sum(score$ratio_group == "High"), sum(score$ratio_group == "Low"),
    sum(score$ratio_group == "At median")
  )
)
fwrite(ratio_audit, file.path(table_dir, "01_ratio_transformation_audit.csv"))
fwrite(score[ratio_defined == FALSE, .(cell_id, Patient, Sample, final_lineage, degradation_module_raw)],
       file.path(table_dir, "01_undefined_ratio_cells.csv"))

score_summary <- score[, .(
  cells = .N,
  production_raw_min = min(production_module_raw),
  production_raw_median = median(production_module_raw),
  production_raw_max = max(production_module_raw),
  degradation_raw_min = min(degradation_module_raw),
  degradation_raw_median = median(degradation_module_raw),
  degradation_raw_max = max(degradation_module_raw),
  production_shifted_median = median(production_module_shifted),
  degradation_shifted_median = median(degradation_module_shifted),
  finite_ratio_cells = sum(ratio_defined),
  ratio_median = median(production_degradation_ratio, na.rm = TRUE),
  ratio_q1 = quantile(production_degradation_ratio, 0.25, na.rm = TRUE),
  ratio_q3 = quantile(production_degradation_ratio, 0.75, na.rm = TRUE),
  ratio_max = max(production_degradation_ratio, na.rm = TRUE)
), by = .(final_lineage, analysis_eligible)]
fwrite(score_summary, file.path(table_dir, "01_score_summary_by_lineage.csv"))

concordance <- data.table(
  module = c("Production", "Degradation"),
  genes = c("GNMT;DMGDH", "SARDH;PIPOX"),
  spearman_module_vs_direct = c(
    cor(score$production_module_raw, score$production_direct_mean, method = "spearman"),
    cor(score$degradation_module_raw, score$degradation_direct_mean, method = "spearman")
  ),
  pearson_module_vs_direct = c(
    cor(score$production_module_raw, score$production_direct_mean, method = "pearson"),
    cor(score$degradation_module_raw, score$degradation_direct_mean, method = "pearson")
  )
)
fwrite(concordance, file.path(table_dir, "01_module_direct_expression_concordance.csv"))

fwrite(score, file.path(table_dir, "01_cell_paper_style_scores.csv.gz"), compress = "gzip")
saveRDS(score, file.path(intermediate_dir, "01_cell_paper_style_scores.rds"), compress = TRUE)

lineage_group <- score[analysis_eligible == TRUE, .(
  cells = .N,
  patients = uniqueN(Patient)
), by = .(final_lineage, ratio_group)]
lineage_total <- score[analysis_eligible == TRUE, .(lineage_cells = .N), by = final_lineage]
lineage_group <- merge(lineage_group, lineage_total, by = "final_lineage", all.x = TRUE)
lineage_group[, cell_fraction_within_lineage := cells / lineage_cells]
fwrite(lineage_group, file.path(table_dir, "01_ratio_group_by_lineage.csv"))
fwrite(score[, .(cells = .N), by = .(Patient, ratio_group)][order(Patient, ratio_group)],
       file.path(table_dir, "01_ratio_group_by_patient.csv"))

# -----------------------------------------------------------------------------
# Figure 3C-style lineage UMAP
# -----------------------------------------------------------------------------
plot_cells <- score[analysis_eligible == TRUE]
plot_cells[, final_lineage := factor(final_lineage, levels = lineage_order)]
if (anyNA(plot_cells$final_lineage)) stop("Unexpected analysis-eligible lineage")

p_lineage <- ggplot(plot_cells, aes(umap_1, umap_2, colour = final_lineage)) +
  geom_point(size = 0.09, alpha = 0.88, stroke = 0) +
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
save_plot("Fig3_style_C_lineage_UMAP", p_lineage, 4.25, 3.75)

# -----------------------------------------------------------------------------
# Figure 3D-style violin distributions for all three author-requested axes
# -----------------------------------------------------------------------------
violin <- rbindlist(list(
  plot_cells[, .(cell_id, final_lineage, axis = "Production", value = production_module_shifted)],
  plot_cells[, .(cell_id, final_lineage, axis = "Degradation", value = degradation_module_shifted)],
  plot_cells[ratio_defined == TRUE, .(
    cell_id, final_lineage, axis = "Production/degradation ratio",
    value = production_degradation_ratio
  )]
))
violin[, axis := factor(axis, levels = c("Production", "Degradation", "Production/degradation ratio"))]
violin_caps <- violin[, .(
  display_cap_99_5 = as.numeric(quantile(value, 0.995, na.rm = TRUE)),
  maximum = max(value, na.rm = TRUE), cells = .N
), by = axis]
violin <- merge(violin, violin_caps, by = "axis", all.x = TRUE, sort = FALSE)
violin[, value_display := pmin(value, display_cap_99_5)]
fwrite(violin_caps, file.path(table_dir, "plotdata_01_violin_display_caps.csv"))
fwrite(violin, file.path(table_dir, "plotdata_01_violin_scores.csv.gz"), compress = "gzip")

p_violin <- ggplot(violin, aes(final_lineage, value_display, fill = final_lineage)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.20, colour = COL_TEXT) +
  facet_wrap(~axis, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = lineage_cols, drop = FALSE) +
  labs(
    title = "Sarcosine functional module scores across cell lineages",
    subtitle = "Seurat AddModuleScore; upper 0.5% clipped within each axis for display",
    x = NULL, y = "Shifted module score / ratio"
  ) +
  theme_nc() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.4),
    strip.text = element_text(size = 6.8, face = "bold"),
    panel.spacing.y = unit(2.5, "mm")
  )
save_plot("Fig3_style_D_three_axis_violins", p_violin, 7.20, 6.10)

# Figure 3D in Lee et al. presents the ratio across annotated cell types.  A
# dedicated ratio panel is more legible here than the three-axis diagnostic
# because Production is zero-inflated and has a very long positive tail.
ratio_violin <- copy(plot_cells[ratio_defined == TRUE])
ratio_display_cap <- as.numeric(quantile(ratio_violin$production_degradation_ratio, 0.99))
ratio_violin[, ratio_display := pmin(production_degradation_ratio, ratio_display_cap)]
fwrite(
  data.table(
    display_rule = "Upper 1% winsorized for display only",
    quantile = 0.99,
    display_cap = ratio_display_cap,
    raw_maximum = max(ratio_violin$production_degradation_ratio),
    cells = nrow(ratio_violin)
  ),
  file.path(table_dir, "plotdata_01_ratio_violin_display_cap.csv")
)
p_ratio_violin <- ggplot(ratio_violin, aes(final_lineage, ratio_display, fill = final_lineage)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, colour = COL_TEXT) +
  geom_boxplot(
    width = 0.08, outlier.shape = NA, fill = "white", alpha = 0.80,
    linewidth = 0.28, colour = COL_TEXT
  ) +
  scale_fill_manual(values = lineage_cols, drop = FALSE) +
  labs(
    title = "Production/degradation ratio across cell lineages",
    subtitle = "Cell-level shifted AddModuleScore ratio; upper 1% winsorized for display only",
    x = NULL, y = "Production/degradation ratio"
  ) +
  theme_nc() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.5)
  )
save_plot("Fig3_style_D_ratio_violin", p_ratio_violin, 7.20, 3.25)

# -----------------------------------------------------------------------------
# Figure 3E-style continuous and high/low UMAPs
# -----------------------------------------------------------------------------
umap_long <- rbindlist(list(
  score[, .(cell_id, umap_1, umap_2, axis = "Production", value = production_module_shifted)],
  score[, .(cell_id, umap_1, umap_2, axis = "Degradation", value = degradation_module_shifted)],
  score[ratio_defined == TRUE, .(
    cell_id, umap_1, umap_2, axis = "Production/degradation ratio",
    value = production_degradation_ratio
  )]
))
umap_long[, axis := factor(axis, levels = c("Production", "Degradation", "Production/degradation ratio"))]
umap_caps <- umap_long[, .(
  display_cap_99_5 = as.numeric(quantile(value, 0.995, na.rm = TRUE)),
  maximum = max(value, na.rm = TRUE), cells = .N
), by = axis]
umap_long <- merge(umap_long, umap_caps, by = "axis", all.x = TRUE, sort = FALSE)
umap_long[, value_display := pmin(value, display_cap_99_5)]
fwrite(umap_caps, file.path(table_dir, "plotdata_01_UMAP_display_caps.csv"))

axis_gradients <- c(
  "Production" = COL_PROD,
  "Degradation" = COL_DEG,
  "Production/degradation ratio" = COL_RATIO
)
make_score_umap <- function(axis_name) {
  d <- umap_long[axis == axis_name][sample(.N)]
  high_col <- axis_gradients[[axis_name]]
  legend_breaks <- pretty(c(0, max(d$value_display, na.rm = TRUE)), n = 3)
  legend_breaks <- legend_breaks[legend_breaks >= 0 & legend_breaks <= max(d$value_display, na.rm = TRUE)]
  ggplot(d, aes(umap_1, umap_2, colour = value_display)) +
    geom_point(size = 0.08, alpha = 0.90, stroke = 0) +
    scale_colour_gradient(
      low = "#D0D0D0", high = high_col, name = "Score",
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
    guides(colour = guide_colourbar(title.position = "top", barwidth = unit(24, "mm"),
                                    barheight = unit(1.8, "mm")))
}
p_prod <- make_score_umap("Production")
p_deg <- make_score_umap("Degradation")
p_ratio <- make_score_umap("Production/degradation ratio")
p_continuous <- p_prod + p_deg + p_ratio + plot_layout(ncol = 3)
save_plot("Fig3_style_E_continuous_UMAPs", p_continuous, 7.20, 2.65)

score[, ratio_group_plot := factor(
  ratio_group,
  levels = c("Low", "High", "At median", "Undefined denominator")
)]
highlow_cols <- c(
  "Low" = COL_LOW, "High" = COL_HIGH,
  "At median" = "#AFAFAF", "Undefined denominator" = "#E5E5E5"
)
p_highlow <- ggplot(score[sample(.N)], aes(umap_1, umap_2, colour = ratio_group_plot)) +
  geom_point(size = 0.085, alpha = 0.86, stroke = 0) +
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
save_plot("Fig3_style_E_ratio_high_low_UMAP", p_highlow, 4.25, 3.75)

# -----------------------------------------------------------------------------
# Figure 3F-style high-versus-low DEG and GO BP enrichment
# -----------------------------------------------------------------------------
obj$ratio_group_fig3 <- score$ratio_group[match(Cells(obj), score$cell_id)]
if (anyNA(obj$ratio_group_fig3)) stop("Ratio group failed to align to Seurat object")
Idents(obj) <- "ratio_group_fig3"

deg_file <- file.path(table_dir, "02_high_vs_low_cell_FindMarkers.csv")
verified_deg_sha <- "af5b154e2d16434be435204ab511449a1713206773fc4d62087a8325e669dc09"
reuse_verified_deg <- identical(Sys.getenv("REUSE_VERIFIED_DEG", unset = "0"), "1")
if (reuse_verified_deg) {
  if (!file.exists(deg_file) || !identical(sha256(deg_file), verified_deg_sha)) {
    stop("Requested DEG reuse, but the twice-reproduced DEG SHA-256 is absent or mismatched")
  }
  message("Reusing twice-reproduced High-versus-Low DEG table after SHA-256 verification")
  deg <- fread(deg_file)
} else {
  message("Running paper-style cell-level High-versus-Low FindMarkers")
  deg <- FindMarkers(
    obj, ident.1 = "High", ident.2 = "Low", assay = "RNA", slot = "data",
    test.use = "wilcox", logfc.threshold = 0.1, min.pct = 0.01,
    min.diff.pct = -Inf, only.pos = FALSE, max.cells.per.ident = Inf,
    random.seed = 260825, densify = FALSE, verbose = TRUE
  )
  deg <- as.data.table(deg, keep.rownames = "gene")
  if (!"avg_log2FC" %in% names(deg)) stop("FindMarkers output lacks avg_log2FC")
  deg[, direction := fifelse(
    p_val_adj < 0.05 & avg_log2FC >= 0.25, "High",
    fifelse(p_val_adj < 0.05 & avg_log2FC <= -0.25, "Low", "Not selected")
  )]
  fwrite(deg, deg_file)
}
if (!all(c("gene", "avg_log2FC", "p_val_adj", "direction") %in% names(deg))) {
  stop("High-versus-Low DEG table schema is incomplete")
}

entrez <- AnnotationDbi::mapIds(
  org.Hs.eg.db, keys = deg$gene, keytype = "SYMBOL", column = "ENTREZID",
  multiVals = "first"
)
deg[, ENTREZID := unname(entrez[gene])]
universe <- unique(na.omit(deg$ENTREZID))
high_genes <- unique(na.omit(deg[direction == "High"]$ENTREZID))
low_genes <- unique(na.omit(deg[direction == "Low"]$ENTREZID))

deg_audit <- data.table(
  metric = c(
    "tested_genes", "mapped_universe_entrez", "selected_high_symbols",
    "selected_low_symbols", "selected_high_entrez", "selected_low_entrez"
  ),
  value = c(
    nrow(deg), length(universe), sum(deg$direction == "High"),
    sum(deg$direction == "Low"), length(high_genes), length(low_genes)
  )
)
fwrite(deg_audit, file.path(table_dir, "02_DEG_enrichment_input_audit.csv"))

if (length(high_genes) > 0L && length(low_genes) > 0L) {
  go <- as.data.table(
    limma::goana(list(High = high_genes, Low = low_genes), universe = universe, species = "Hs"),
    keep.rownames = "GO_ID"
  )
  required_go <- c("GO_ID", "Term", "Ont", "N", "High", "P.High", "Low", "P.Low")
  if (!all(required_go %in% names(go))) stop("Unexpected goana output schema")
  go_bp <- go[Ont == "BP"]
  go_bp[, `:=`(
    BH_High = p.adjust(P.High, method = "BH"),
    BH_Low = p.adjust(P.Low, method = "BH")
  )]
  fwrite(go_bp, file.path(table_dir, "03_high_low_GO_BP_overrepresentation.csv"))

  go_long <- rbindlist(list(
    go_bp[, .(GO_ID, Term, direction = "High-ratio cells", N, DE = High, p = P.High, q = BH_High)],
    go_bp[, .(GO_ID, Term, direction = "Low-ratio cells", N, DE = Low, p = P.Low, q = BH_Low)]
  ))
  go_long <- go_long[N >= 10 & N <= 500 & DE > 0]
  go_top <- go_long[order(q, p), head(.SD, 10L), by = direction]
  go_top[, Term_plot := factor(Term, levels = rev(unique(Term[order(q, p)])))]
  go_top[, gene_fraction := DE / N]
  fwrite(go_top, file.path(table_dir, "plotdata_03_top_GO_BP_terms.csv"))

  if (nrow(go_top) > 0L) {
    p_go <- ggplot(go_top, aes(gene_fraction, Term_plot)) +
      geom_point(aes(size = DE, colour = -log10(pmax(q, .Machine$double.xmin)))) +
      facet_wrap(~direction, scales = "free_y", ncol = 2) +
      scale_colour_gradient(low = "#E9C3A2", high = COL_HIGH, name = expression(-log[10](q))) +
      scale_size_continuous(name = "Genes", range = c(1.8, 5.2)) +
      labs(
        title = "GO Biological Process enrichment",
        subtitle = "Exploratory cell-level High versus Low ratio comparison",
        x = "DE genes / represented term genes", y = NULL
      ) +
      theme_nc() +
      theme(
        axis.text.y = element_text(size = 5.7),
        legend.position = "bottom",
        panel.spacing.x = unit(7, "mm")
      )
    save_plot("Fig3_style_F_high_low_GO_BP", p_go, 7.20, 4.50)
  }
} else {
  fwrite(data.table(note = "No genes were available in both directions for GO BP ORA"),
         file.path(table_dir, "03_high_low_GO_BP_overrepresentation_NOT_RUN.csv"))
}

# A compact Figure 3C-E analogue; GO remains a separate panel because term labels
# become unreadable when forced into the same footprint.
p_compact <- (p_lineage | p_highlow) / p_continuous +
  plot_annotation(
    title = "Lee et al. Figure 3-style sarcosine functional-score analysis",
    subtitle = "Production: GNMT+DMGDH; Degradation: SARDH+PIPOX"
  ) &
  theme(
    plot.title = element_text(family = "Arial", size = 8, face = "bold"),
    plot.subtitle = element_text(family = "Arial", size = 6.5)
  )
save_plot("Fig_scRNA_Lee2024_Fig3_style_compact", p_compact, 7.20, 6.20)

figure_specs <- data.table(
  stem = c(
    "Fig3_style_C_lineage_UMAP", "Fig3_style_D_three_axis_violins", "Fig3_style_D_ratio_violin",
    "Fig3_style_E_continuous_UMAPs", "Fig3_style_E_ratio_high_low_UMAP",
    "Fig3_style_F_high_low_GO_BP", "Fig_scRNA_Lee2024_Fig3_style_compact"
  ),
  width_in = c(4.25, 7.20, 7.20, 7.20, 4.25, 7.20, 7.20),
  height_in = c(3.75, 6.10, 3.25, 2.65, 3.75, 4.50, 6.20),
  dpi = 600L,
  intended_role = c(
    "Figure 3C-style lineage UMAP", "Figure 3D-style lineage violin extension",
    "Figure 3D-style dedicated ratio violin",
    "Figure 3E-style continuous score UMAPs", "Figure 3E-style median High/Low UMAP",
    "Figure 3F-style exploratory GO BP", "compact Figure 3C-E analogue"
  )
)
figure_specs[, png_expected_width_px := as.integer(round(width_in * dpi))]
figure_specs[, png_expected_height_px := as.integer(round(height_in * dpi))]
fwrite(figure_specs, file.path(table_dir, "04_figure_specifications.csv"))

capture.output(sessionInfo(), file = file.path(log_dir, "01_sessionInfo.txt"))

cat("LEE FIGURE 3-STYLE ANALYSIS PASS\n")
cat("Cells scored:", nrow(score), "\n")
cat("Finite ratios:", sum(score$ratio_defined), "\n")
cat("Undefined shifted denominator:", sum(!score$ratio_defined), "\n")
cat("Median ties:", sum(score$ratio_group == "At median"), "\n")
cat("High/Low cells:", sum(score$ratio_group == "High"), "/", sum(score$ratio_group == "Low"), "\n")
cat("DE genes High/Low:", sum(deg$direction == "High"), "/", sum(deg$direction == "Low"), "\n")
