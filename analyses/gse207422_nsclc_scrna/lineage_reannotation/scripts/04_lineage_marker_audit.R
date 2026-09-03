#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
.libPaths(c(file.path(analysis_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

options(stringsAsFactors = FALSE, warn = 1)

table_dir <- file.path(analysis_dir, "results", "tables")
intermediate_dir <- file.path(analysis_dir, "intermediate")
log_dir <- file.path(analysis_dir, "logs")

counts_file <- file.path(intermediate_dir, "01_full_counts_mt20_qc.rds")
cluster_file <- file.path(table_dir, "03_integrated_clusters_and_umap.csv")
if (!all(file.exists(c(counts_file, cluster_file)))) stop("Required lineage-audit inputs are missing")

counts <- readRDS(counts_file)
clusters <- fread(cluster_file)
if (!inherits(counts, "dgCMatrix")) stop("Counts are not a dgCMatrix")
if (ncol(counts) != 92053L || nrow(counts) != 24292L) stop("Unexpected count-matrix dimensions")
if (nrow(clusters) != ncol(counts) || anyDuplicated(clusters$cell_id)) stop("Invalid cluster table")
ord <- match(colnames(counts), clusters$cell_id)
if (anyNA(ord)) stop("Counts/cluster cell-ID mismatch")
clusters <- clusters[ord]
if (!identical(clusters$cell_id, colnames(counts))) stop("Counts/cluster ordering failed")
if (!identical(as.numeric(clusters$library_size), as.numeric(Matrix::colSums(counts)))) stop("Library-size alignment failed")

# Marker panels copied verbatim from the authors' downloaded analysis scripts.
authors_all_cell <- c(
  "PTPRC", "CD3E", "LYZ", "CD79A", "MS4A1", "IGHG1", "CSF3R", "FGFBP2",
  "KIT", "LILRA4", "VWF", "COL1A1", "EPCAM", "MKI67"
)
authors_tnk_short <- c(
  "CD8A", "CD4", "FCGR3A", "CTLA4", "PDCD1", "TIGIT", "HAVCR2", "GZMA",
  "GZMB", "GZMK", "NKG7", "SELL", "TCF7", "FOXP3", "ZNF683", "ITGAE"
)
authors_tnk_heatmap <- c(
  "CD8A", "CD8B", "CD4", "FCGR3A", "KLRD1", "FOXP3", "IL2RA", "IKZF2",
  "TCF7", "SELL", "LEF1", "CCR7", "LAG3", "TIGIT", "PDCD1", "HAVCR2",
  "CTLA4", "LAYN", "ENTPD1", "GZMA", "GZMB", "GZMK", "GNLY", "IFNG",
  "PRF1", "NKG7", "CD28", "ICOS", "CD40LG", "TNFRSF4", "TNFRSF9",
  "TNFRSF18", "ZNF683", "ITGAE", "CCL3", "CCL5", "CXCL13", "IL21",
  "EOMES", "MAF", "TOX2", "ID2", "TBX21", "HOPX", "MKI67", "STMN1"
)
authors_myeloid <- c(
  "VCAN", "FCN1", "S100A8", "S100A9", "FABP4", "MCEMP1", "MARCO", "C1QA",
  "C1QB", "GPNMB", "APOE", "SPP1", "SELENOP", "MRC1", "TGFB1", "CD163",
  "CCL18", "MSR1", "VEGFA", "TNF", "CXCL9", "CXCL10", "CXCL11", "HLA-DRA",
  "HLA-DQA1", "HLA-DPA1", "CD74", "MKI67", "TOP2A"
)
authors_dc <- c(
  "XCR1", "CLEC9A", "CADM1", "CLNK", "CD1C", "CD1E", "FCER1A", "CLEC10A",
  "LAMP3", "FSCN1", "CCR7", "LAD1", "CCL17", "CCL19", "CCL22", "CX3CR1",
  "CD86", "CD83", "CD80", "CD40", "ICOSLG"
)
authors_b <- c("MS4A1", "CD27", "GPR183", "IGHD", "RGS13", "IGHM")
authors_neutrophil <- c(
  "FCGR3B", "SELL", "CXCR4", "CXCR2", "CXCR1", "CXCL8", "IL1B", "CCL3",
  "CCL4", "CCL4L2", "S100A12", "S100A9", "S100A8", "CYBB", "ELANE", "MMP9",
  "PADI4", "HMGB1", "TNF", "CXCL9", "CXCL10", "OSM", "PTGS2", "ARG1",
  "TGFB1", "VEGFA", "PROK2", "IFIT1", "IFIT2", "IFIT3", "RSAD2", "CD274", "IDO1"
)

marker_manifest <- rbindlist(list(
  data.table(panel = "authors_all_cell", gene = authors_all_cell),
  data.table(panel = "authors_T_NK_short", gene = authors_tnk_short),
  data.table(panel = "authors_T_NK_heatmap", gene = authors_tnk_heatmap),
  data.table(panel = "authors_myeloid_macrophage", gene = authors_myeloid),
  data.table(panel = "authors_DC", gene = authors_dc),
  data.table(panel = "authors_B", gene = authors_b),
  data.table(panel = "authors_neutrophil", gene = authors_neutrophil)
))
marker_manifest[, present := gene %in% rownames(counts)]
fwrite(marker_manifest, file.path(table_dir, "04_authors_marker_manifest.csv"))
if (!all(marker_manifest[panel == "authors_all_cell"]$present)) stop("An authors' all-cell marker is absent")

marker_genes <- unique(marker_manifest[present == TRUE]$gene)
marker_counts <- counts[marker_genes, , drop = FALSE]
marker_norm <- marker_counts
marker_norm@x <- log1p(10000 * marker_norm@x / rep(clusters$library_size, diff(marker_norm@p)))

make_membership <- function(cluster_vector) {
  lev <- sort(unique(as.character(cluster_vector)))
  sparseMatrix(
    i = seq_along(cluster_vector), j = match(as.character(cluster_vector), lev), x = 1,
    dims = c(length(cluster_vector), length(lev)),
    dimnames = list(NULL, lev)
  )
}

marker_summary_one <- function(cluster_col, resolution_name) {
  membership <- make_membership(clusters[[cluster_col]])
  n_cells <- as.numeric(Matrix::colSums(membership))
  sum_expr <- as.matrix(marker_norm %*% membership)
  detected <- marker_counts
  detected@x[] <- 1
  n_detected <- as.matrix(detected %*% membership)
  rm(detected)
  out <- rbindlist(lapply(seq_len(ncol(membership)), function(j) {
    data.table(
      resolution = resolution_name,
      cluster = colnames(membership)[j],
      cells = n_cells[j],
      gene = rownames(marker_norm),
      avg_log_normalized_expression = sum_expr[, j] / n_cells[j],
      percent_expressing = 100 * n_detected[, j] / n_cells[j]
    )
  }))
  out
}

message("Computing authors' marker profiles")
marker_summary <- rbindlist(list(
  marker_summary_one("drugres_snn_res.0.1", "primary_res_0.1"),
  marker_summary_one("source_snn_res.0.6", "validation_res_0.6")
))
fwrite(marker_summary, file.path(table_dir, "04_cluster_authors_marker_profiles.csv"))
rm(marker_counts, marker_norm)
invisible(gc())

# Descriptive, non-inferential top-cluster genes from mean log-normalized
# expression. These are used only to audit annotation and are not treated as
# biological-replicate differential-expression results.
message("Building full-gene log-normalized matrix for descriptive marker ranking")
norm <- counts
norm@x <- log1p(10000 * norm@x / rep(clusters$library_size, diff(norm@p)))
overall_sum <- as.numeric(Matrix::rowSums(norm))

top_descriptive_one <- function(cluster_col, resolution_name, top_n = 50L) {
  membership <- make_membership(clusters[[cluster_col]])
  n_cells <- as.numeric(Matrix::colSums(membership))
  sum_in <- as.matrix(norm %*% membership)
  out <- rbindlist(lapply(seq_len(ncol(membership)), function(j) {
    avg_in <- sum_in[, j] / n_cells[j]
    avg_out <- (overall_sum - sum_in[, j]) / (ncol(norm) - n_cells[j])
    delta <- avg_in - avg_out
    keep <- head(order(delta, decreasing = TRUE, na.last = NA), top_n)
    data.table(
      resolution = resolution_name,
      cluster = colnames(membership)[j],
      cells = n_cells[j],
      rank = seq_along(keep),
      gene = rownames(norm)[keep],
      avg_log_normalized_in_cluster = avg_in[keep],
      avg_log_normalized_outside_cluster = avg_out[keep],
      descriptive_delta = delta[keep]
    )
  }))
  out
}

top_descriptive <- rbindlist(list(
  top_descriptive_one("drugres_snn_res.0.1", "primary_res_0.1"),
  top_descriptive_one("source_snn_res.0.6", "validation_res_0.6")
))
fwrite(top_descriptive, file.path(table_dir, "04_cluster_top50_descriptive_genes.csv"))
rm(norm)
invisible(gc())

crosswalk <- clusters[, .N, by = .(
  primary_cluster = as.character(`drugres_snn_res.0.1`),
  validation_cluster = as.character(`source_snn_res.0.6`)
)][order(primary_cluster, -N)]
crosswalk[, percent_of_primary_cluster := 100 * N / sum(N), by = primary_cluster]
crosswalk[, percent_of_validation_cluster := 100 * N / sum(N), by = validation_cluster]
fwrite(crosswalk, file.path(table_dir, "04_primary_validation_cluster_crosswalk.csv"))

source_map <- data.table(
  source_cluster = as.character(0:25),
  source_label = c(
    "T cell", "B cell", "T cell", "Neutrophil", "Epithelium", "T cell", "T cell",
    "Myeloid cell", "Myeloid cell", "Myeloid cell", "Neutrophil", "Plasma cell",
    "T cell", "Epithelium", "T cell", "Cycling immune cell", "Epithelium",
    "Epithelium", "T cell", "NK cell", "Fibroblast/Endothelium", "Mast cell",
    "Epithelium", "pDC", "Plasma cell", "Myeloid cell"
  )
)
fwrite(source_map, file.path(table_dir, "04_source_study_exact_26_cluster_map_reference.csv"))

expected_2024_labels <- data.table(label = c(
  "B cell", "CAF", "CD4 T", "CD8 T", "Epithelial", "Mast", "Macrophage",
  "Neutrophil", "NK", "pDC", "Plasma"
))
fwrite(expected_2024_labels, file.path(table_dir, "04_drug_resistance_2024_reported_lineages_reference.csv"))

capture.output(sessionInfo(), file = file.path(log_dir, "04_sessionInfo.txt"))
cat("Lineage marker audit PASS\n")
cat("Marker genes present:", length(marker_genes), "of", length(unique(marker_manifest$gene)), "unique authors' markers\n")
cat("Primary clusters:", uniqueN(clusters$`drugres_snn_res.0.1`), "\n")
cat("Validation clusters:", uniqueN(clusters$`source_snn_res.0.6`), "\n")
