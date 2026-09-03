#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
.libPaths(c(file.path(analysis_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

options(stringsAsFactors = FALSE, warn = 1)

table_dir <- file.path(analysis_dir, "results", "tables")
figure_dir <- file.path(analysis_dir, "results", "figures_diagnostic")
intermediate_dir <- file.path(analysis_dir, "intermediate")
log_dir <- file.path(analysis_dir, "logs")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

counts_file <- file.path(intermediate_dir, "01_full_counts_mt20_qc.rds")
cluster_file <- file.path(table_dir, "05_tnk_clusters_and_umap.csv")
if (!all(file.exists(c(counts_file, cluster_file)))) stop("Required T/NK marker-audit inputs are missing")

counts <- readRDS(counts_file)
meta <- fread(cluster_file)
if (!inherits(counts, "dgCMatrix")) stop("Counts are not a dgCMatrix")
if (nrow(meta) != 38730L || anyDuplicated(meta$cell_id)) stop("Invalid T/NK cluster table")
idx <- match(meta$cell_id, colnames(counts))
if (anyNA(idx)) stop("T/NK cluster cells absent from full count matrix")
counts <- counts[, idx, drop = FALSE]
if (!identical(colnames(counts), meta$cell_id)) stop("T/NK counts/cluster ordering failed")
if (!identical(as.numeric(Matrix::colSums(counts)), as.numeric(meta$library_size))) {
  stop("T/NK library sizes disagree with counts")
}

# Marker genes are copied from the source team's T/NK script where available;
# the contamination panel is used only to identify incompatible major lineages.
marker_panels <- list(
  `T lineage` = c("CD3D", "CD3E", "TRAC"),
  `CD4/naive-memory` = c("CD4", "IL7R", "CCR7", "TCF7", "SELL", "LEF1"),
  `CD8/cytotoxic` = c("CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "GNLY", "IFNG", "PRF1", "NKG7"),
  `NK` = c("FCGR3A", "KLRD1", "FGFBP2", "XCL1", "XCL2", "TYROBP"),
  `Treg` = c("FOXP3", "IL2RA", "IKZF2", "CTLA4", "TNFRSF18"),
  `Exhaustion/TRM` = c("LAG3", "TIGIT", "PDCD1", "HAVCR2", "LAYN", "ENTPD1", "ZNF683", "ITGAE", "CXCL13"),
  `Cycling` = c("MKI67", "STMN1", "TOP2A"),
  `Incompatible lineage audit` = c(
    "LYZ", "FCER1G", "LST1", "C1QC", "CSF3R", "FCGR3B", "MS4A1", "CD79A",
    "MZB1", "JCHAIN", "EPCAM", "KRT19", "COL1A1", "DCN", "TPSB2", "KIT",
    "LILRA4", "GPR183", "TCF4"
  )
)
marker_manifest <- rbindlist(lapply(names(marker_panels), function(panel) {
  data.table(marker_category = panel, gene = marker_panels[[panel]])
}))
marker_manifest[, present := gene %in% rownames(counts)]
fwrite(marker_manifest, file.path(table_dir, "06_tnk_marker_manifest.csv"))
if (!all(marker_manifest$present)) {
  stop("Missing T/NK audit markers: ", paste(marker_manifest[!present]$gene, collapse = ", "))
}

cluster_col <- "tnk_source_snn_res.0.4"
cluster_ids <- sort(unique(as.character(meta[[cluster_col]])))
membership <- sparseMatrix(
  i = seq_len(nrow(meta)), j = match(as.character(meta[[cluster_col]]), cluster_ids), x = 1,
  dims = c(nrow(meta), length(cluster_ids)), dimnames = list(NULL, cluster_ids)
)
n_cells <- as.numeric(Matrix::colSums(membership))

marker_genes <- unique(marker_manifest$gene)
marker_counts <- counts[marker_genes, , drop = FALSE]
marker_norm <- marker_counts
marker_norm@x <- log1p(10000 * marker_norm@x / rep(meta$library_size, diff(marker_norm@p)))
sum_expr <- as.matrix(marker_norm %*% membership)
detected <- marker_counts
detected@x[] <- 1
n_detected <- as.matrix(detected %*% membership)

marker_summary <- rbindlist(lapply(seq_along(cluster_ids), function(j) {
  data.table(
    cluster = cluster_ids[j], cells = n_cells[j], gene = marker_genes,
    avg_log_normalized_expression = sum_expr[, j] / n_cells[j],
    percent_expressing = 100 * n_detected[, j] / n_cells[j]
  )
}))
marker_summary <- merge(marker_summary, marker_manifest[, .(marker_category, gene)], by = "gene", all.x = TRUE)
setcolorder(marker_summary, c("cluster", "cells", "marker_category", "gene",
                             "avg_log_normalized_expression", "percent_expressing"))
fwrite(marker_summary, file.path(table_dir, "06_tnk_cluster_marker_profiles.csv"))
rm(marker_counts, marker_norm, detected, sum_expr, n_detected)
invisible(gc())

# Full-gene descriptive ranking. This is not replicate-level differential
# expression and is used only to support manual cluster annotation.
message("Computing descriptive T/NK top-cluster genes")
norm <- counts
norm@x <- log1p(10000 * norm@x / rep(meta$library_size, diff(norm@p)))
overall_sum <- as.numeric(Matrix::rowSums(norm))
sum_in <- as.matrix(norm %*% membership)
top_descriptive <- rbindlist(lapply(seq_along(cluster_ids), function(j) {
  avg_in <- sum_in[, j] / n_cells[j]
  avg_out <- (overall_sum - sum_in[, j]) / (ncol(norm) - n_cells[j])
  delta <- avg_in - avg_out
  keep <- head(order(delta, decreasing = TRUE, na.last = NA), 75L)
  data.table(
    cluster = cluster_ids[j], cells = n_cells[j], rank = seq_along(keep),
    gene = rownames(norm)[keep], avg_log_normalized_in_cluster = avg_in[keep],
    avg_log_normalized_outside_cluster = avg_out[keep], descriptive_delta = delta[keep]
  )
}))
fwrite(top_descriptive, file.path(table_dir, "06_tnk_cluster_top75_descriptive_genes.csv"))
rm(norm, overall_sum, sum_in)
invisible(gc())

sample_composition <- meta[, .N, by = c(cluster_col, "Sample", "Patient", "clinical_group")]
sample_composition[, cluster_total := sum(N), by = cluster_col]
sample_composition[, fraction_of_cluster := N / cluster_total]
fwrite(sample_composition[order(get(cluster_col), -N)],
       file.path(table_dir, "06_tnk_cluster_sample_composition.csv"))
cluster_audit <- sample_composition[, .(
  cells = unique(cluster_total),
  samples_present = sum(N > 0),
  samples_with_at_least_10_cells = sum(N >= 10),
  max_single_sample_fraction = max(fraction_of_cluster),
  dominant_sample = Sample[which.max(fraction_of_cluster)]
), by = cluster_col]
fwrite(cluster_audit[order(get(cluster_col))], file.path(table_dir, "06_tnk_cluster_technical_audit.csv"))

# Diagnostic dot plot at a fixed physical size. It is an audit graphic, not a
# manuscript-ready inferential panel.
plot_data <- copy(marker_summary)
plot_data[, scaled_average := {
  z <- as.numeric(scale(avg_log_normalized_expression))
  if (all(!is.finite(z))) rep(0, .N) else pmax(-2, pmin(2, z))
}, by = gene]
plot_data[, cluster := factor(cluster, levels = cluster_ids)]
plot_data[, marker_category := factor(marker_category, levels = names(marker_panels))]
plot_data[, gene := factor(gene, levels = rev(unique(marker_manifest$gene)))]

p <- ggplot(plot_data, aes(x = cluster, y = gene)) +
  geom_point(aes(size = percent_expressing, fill = scaled_average), shape = 21,
             colour = "#202020", stroke = 0.22) +
  facet_grid(marker_category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_size_continuous(name = "% expressing", range = c(0.2, 4.2), limits = c(0, 100)) +
  scale_fill_gradient2(name = "Mean expression\n(gene-scaled)", low = "#2E5F8A",
                       mid = "#F4F4F4", high = "#C47B3B", midpoint = 0,
                       limits = c(-2, 2)) +
  labs(x = "T/NK cluster (resolution 0.4)", y = NULL) +
  theme_classic(base_family = "Arial", base_size = 7) +
  theme(
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5),
    axis.text = element_text(size = 6, colour = "black"),
    axis.title.x = element_text(size = 7),
    strip.background = element_blank(),
    strip.text.y.left = element_text(size = 6.5, angle = 0, hjust = 1),
    panel.spacing.y = unit(2, "mm"),
    legend.position = "right"
  )
ggsave(file.path(figure_dir, "06_tnk_marker_dotplot.png"), p,
       width = 7.2, height = 9.0, units = "in", dpi = 300, bg = "white")

capture.output(sessionInfo(), file = file.path(log_dir, "06_sessionInfo.txt"))
cat("T/NK marker audit PASS\n")
cat("Clusters:", length(cluster_ids), "\n")
cat("Cells:", nrow(meta), "\n")
cat("All", nrow(marker_manifest), "audit-marker rows present\n")
