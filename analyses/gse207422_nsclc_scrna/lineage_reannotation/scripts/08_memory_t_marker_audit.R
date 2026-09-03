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

table_dir <- file.path(analysis_dir, "results", "tables")
figure_dir <- file.path(analysis_dir, "results", "figures_diagnostic")
intermediate_dir <- file.path(analysis_dir, "intermediate")
log_dir <- file.path(analysis_dir, "logs")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

counts <- readRDS(file.path(intermediate_dir, "01_full_counts_mt20_qc.rds"))
meta <- fread(file.path(table_dir, "07_memory_t_clusters_and_umap.csv"))
if (!inherits(counts, "dgCMatrix") || nrow(meta) != 7993L || anyDuplicated(meta$cell_id)) {
  stop("Invalid memory-T audit inputs")
}
idx <- match(meta$cell_id, colnames(counts))
if (anyNA(idx)) stop("Memory-T cells absent from full counts")
counts <- counts[, idx, drop = FALSE]
if (!identical(colnames(counts), meta$cell_id) ||
    !identical(as.numeric(Matrix::colSums(counts)), as.numeric(meta$library_size))) {
  stop("Memory-T count/metadata alignment failed")
}

panels <- list(
  `T lineage` = c("CD3D", "CD3E", "TRAC"),
  `CD4 memory` = c("CD4", "IL7R", "CCR7", "TCF7", "MAL", "LTB", "GPR183"),
  `CD8 memory` = c("CD8A", "CD8B", "GZMK", "CCL5", "EOMES", "CRTAM"),
  `Treg/activation` = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "PDCD1", "CXCL13"),
  `Cytotoxic/NK audit` = c("GZMB", "GNLY", "NKG7", "PRF1", "KLRD1", "FCGR3A", "FGFBP2"),
  `Cycling/other lineage audit` = c("MKI67", "STMN1", "LYZ", "FCER1G", "MS4A1", "LILRA4")
)
manifest <- rbindlist(lapply(names(panels), function(panel) {
  data.table(marker_category = panel, gene = panels[[panel]])
}))
manifest[, present := gene %in% rownames(counts)]
fwrite(manifest, file.path(table_dir, "08_memory_t_marker_manifest.csv"))
if (!all(manifest$present)) stop("Missing memory-T audit marker(s)")

cluster_col <- "memory_t_snn_res.0.4"
cluster_ids <- sort(unique(as.character(meta[[cluster_col]])))
membership <- sparseMatrix(
  i = seq_len(nrow(meta)), j = match(as.character(meta[[cluster_col]]), cluster_ids), x = 1,
  dims = c(nrow(meta), length(cluster_ids)), dimnames = list(NULL, cluster_ids)
)
n_cells <- as.numeric(Matrix::colSums(membership))
genes <- unique(manifest$gene)
mc <- counts[genes, , drop = FALSE]
mn <- mc
mn@x <- log1p(10000 * mn@x / rep(meta$library_size, diff(mn@p)))
sum_expr <- as.matrix(mn %*% membership)
det <- mc
det@x[] <- 1
n_det <- as.matrix(det %*% membership)
profiles <- rbindlist(lapply(seq_along(cluster_ids), function(j) {
  data.table(
    cluster = cluster_ids[j], cells = n_cells[j], gene = genes,
    avg_log_normalized_expression = sum_expr[, j] / n_cells[j],
    percent_expressing = 100 * n_det[, j] / n_cells[j]
  )
}))
profiles <- merge(profiles, manifest[, .(marker_category, gene)], by = "gene", all.x = TRUE)
fwrite(profiles, file.path(table_dir, "08_memory_t_cluster_marker_profiles.csv"))
rm(mc, mn, det, sum_expr, n_det)

norm <- counts
norm@x <- log1p(10000 * norm@x / rep(meta$library_size, diff(norm@p)))
overall <- as.numeric(Matrix::rowSums(norm))
sum_in <- as.matrix(norm %*% membership)
top <- rbindlist(lapply(seq_along(cluster_ids), function(j) {
  avg_in <- sum_in[, j] / n_cells[j]
  avg_out <- (overall - sum_in[, j]) / (ncol(norm) - n_cells[j])
  delta <- avg_in - avg_out
  keep <- head(order(delta, decreasing = TRUE), 75L)
  data.table(
    cluster = cluster_ids[j], cells = n_cells[j], rank = seq_along(keep),
    gene = rownames(norm)[keep], avg_in = avg_in[keep], avg_out = avg_out[keep],
    descriptive_delta = delta[keep]
  )
}))
fwrite(top, file.path(table_dir, "08_memory_t_cluster_top75_descriptive_genes.csv"))

sample_comp <- meta[, .N, by = c(cluster_col, "Sample", "clinical_group")]
sample_comp[, cluster_total := sum(N), by = cluster_col]
sample_comp[, fraction_of_cluster := N / cluster_total]
fwrite(sample_comp[order(get(cluster_col), -N)],
       file.path(table_dir, "08_memory_t_cluster_sample_composition.csv"))

plot_data <- copy(profiles)
plot_data[, scaled_average := {
  z <- as.numeric(scale(avg_log_normalized_expression))
  if (all(!is.finite(z))) rep(0, .N) else pmax(-2, pmin(2, z))
}, by = gene]
plot_data[, cluster := factor(cluster, levels = cluster_ids)]
plot_data[, marker_category := factor(marker_category, levels = names(panels))]
plot_data[, gene := factor(gene, levels = rev(unique(manifest$gene)))]
p <- ggplot(plot_data, aes(cluster, gene)) +
  geom_point(aes(size = percent_expressing, fill = scaled_average), shape = 21,
             colour = "#202020", stroke = 0.22) +
  facet_grid(marker_category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_size_continuous(name = "% expressing", range = c(0.2, 4.5), limits = c(0, 100)) +
  scale_fill_gradient2(name = "Mean expression\n(gene-scaled)", low = "#2E5F8A",
                       mid = "#F4F4F4", high = "#C47B3B", limits = c(-2, 2)) +
  labs(x = "Memory-T subcluster (resolution 0.4)", y = NULL) +
  theme_classic(base_family = "Arial", base_size = 7) +
  theme(axis.line = element_line(linewidth = 0.5), axis.ticks = element_line(linewidth = 0.5),
        axis.text = element_text(size = 6, colour = "black"),
        strip.background = element_blank(),
        strip.text.y.left = element_text(size = 6.5, angle = 0, hjust = 1),
        panel.spacing.y = unit(2, "mm"))
ggsave(file.path(figure_dir, "08_memory_t_marker_dotplot.png"), p,
       width = 6.2, height = 7.2, units = "in", dpi = 300, bg = "white")

capture.output(sessionInfo(), file = file.path(log_dir, "08_sessionInfo.txt"))
cat("Memory-T marker audit PASS\n")
cat("Clusters:", length(cluster_ids), "\n")
cat("Cells:", nrow(meta), "\n")
