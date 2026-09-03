#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
.libPaths(c(file.path(analysis_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(future)
})

options(stringsAsFactors = FALSE, warn = 1, future.globals.maxSize = 15 * 1024^3)
future::plan("multicore", workers = 3)
set.seed(260824)

table_dir <- file.path(analysis_dir, "results", "tables")
intermediate_dir <- file.path(analysis_dir, "intermediate")
log_dir <- file.path(analysis_dir, "logs")

counts_file <- file.path(intermediate_dir, "01_full_counts_mt20_qc.rds")
meta_file <- file.path(table_dir, "01_cell_metadata_mt20_qc.csv")
hvg_file <- file.path(intermediate_dir, "02_sample_hvg_rankings.rds")
feature_file <- file.path(table_dir, "02_integration_features_3000.csv")
if (!all(file.exists(c(counts_file, meta_file, hvg_file, feature_file)))) stop("Required integration inputs are missing")

counts <- readRDS(counts_file)
meta <- fread(meta_file)
hvg_list <- readRDS(hvg_file)
integration_features <- fread(feature_file)$gene
if (!inherits(counts, "dgCMatrix")) stop("Counts are not a dgCMatrix")
if (!identical(colnames(counts), meta$cell_id)) stop("Counts/metadata alignment failed")
if (length(integration_features) != 3000L || anyDuplicated(integration_features)) stop("Invalid integration feature set")
if (!all(integration_features %in% rownames(counts))) stop("Integration features absent from counts")

log_normalize_selected <- function(counts_selected, full_library_size, scale_factor = 10000) {
  if (!inherits(counts_selected, "dgCMatrix")) counts_selected <- as(counts_selected, "dgCMatrix")
  if (length(full_library_size) != ncol(counts_selected) || any(full_library_size <= 0)) stop("Invalid full library sizes")
  out <- counts_selected
  denom <- rep(full_library_size, diff(out@p))
  if (length(denom) != length(out@x)) stop("Sparse normalization index mismatch")
  out@x <- log1p(scale_factor * out@x / denom)
  out
}

# Verify that memory-efficient selected-feature normalization exactly reproduces
# Seurat LogNormalize applied to the full gene universe.
test_cells <- meta$cell_id[seq_len(50L)]
test_full <- CreateSeuratObject(counts[, test_cells, drop = FALSE], min.cells = 0, min.features = 0)
test_full <- NormalizeData(test_full, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
test_manual <- log_normalize_selected(
  counts[integration_features, test_cells, drop = FALSE],
  meta$library_size[match(test_cells, meta$cell_id)]
)
test_reference <- LayerData(test_full, assay = "RNA", layer = "data")[integration_features, , drop = FALSE]
normalization_max_abs_diff <- max(abs(test_reference - test_manual))
if (!is.finite(normalization_max_abs_diff) || normalization_max_abs_diff > 1e-12) {
  stop("Selected-feature normalization did not reproduce full-universe Seurat LogNormalize")
}
rm(test_full, test_manual, test_reference)
invisible(gc())

sample_ids <- sort(unique(meta$Sample))
if (!identical(sample_ids, sort(names(hvg_list)))) stop("Sample/HVG list mismatch")
object_list <- setNames(vector("list", length(sample_ids)), sample_ids)
object_audit <- vector("list", length(sample_ids))

for (idx in seq_along(sample_ids)) {
  sid <- sample_ids[idx]
  cell_idx <- which(meta$Sample == sid)
  cells <- meta$cell_id[cell_idx]
  md <- as.data.frame(meta[cell_idx])
  rownames(md) <- md$cell_id
  message("Preparing exact normalized CCA input: ", sid, " (", idx, "/", length(sample_ids), ")")

  counts_selected <- counts[integration_features, cells, drop = FALSE]
  data_selected <- log_normalize_selected(counts_selected, meta$library_size[cell_idx])
  obj <- CreateSeuratObject(
    counts = counts_selected,
    meta.data = md,
    project = sid,
    min.cells = 0,
    min.features = 0
  )
  obj <- SetAssayData(obj, assay = "RNA", layer = "data", new.data = data_selected)
  VariableFeatures(obj) <- intersect(hvg_list[[sid]], integration_features)
  if (!identical(rownames(LayerData(obj, assay = "RNA", layer = "data")), integration_features)) {
    stop("Feature order changed for ", sid)
  }
  object_list[[sid]] <- obj
  object_audit[[idx]] <- data.table(
    Sample = sid,
    cells = ncol(obj),
    integration_features = nrow(obj),
    original_HVGs_in_integration_set = length(VariableFeatures(obj)),
    full_library_sum = sum(meta$library_size[cell_idx]),
    selected_feature_UMI_sum = sum(counts_selected)
  )
  rm(obj, md, counts_selected, data_selected)
  invisible(gc())
}
fwrite(rbindlist(object_audit), file.path(table_dir, "03_cca_input_object_audit.csv"))
fwrite(data.table(metric = "normalization_max_abs_diff", value = normalization_max_abs_diff),
       file.path(table_dir, "03_selected_feature_normalization_unit_check.csv"))

rm(counts)
invisible(gc())

anchor_file <- file.path(intermediate_dir, "03_cca_anchors_checkpoint.rds")
if (file.exists(anchor_file)) {
  message("Loading verified CCA anchor checkpoint")
  anchors <- readRDS(anchor_file)
  if (!inherits(anchors, "AnchorSet")) stop("CCA anchor checkpoint is not an AnchorSet")
} else {
  message("Finding CCA integration anchors across 15 samples")
  anchors <- FindIntegrationAnchors(
    object.list = object_list,
    anchor.features = integration_features,
    normalization.method = "LogNormalize",
    reduction = "cca",
    dims = 1:20,
    l2.norm = TRUE,
    k.anchor = 5,
    k.filter = 200,
    k.score = 30,
    max.features = 200,
    nn.method = "annoy",
    n.trees = 50,
    eps = 0,
    verbose = TRUE
  )
  saveRDS(anchors, anchor_file, compress = FALSE)
  if (!file.exists(anchor_file) || file.info(anchor_file)$size <= 0) stop("CCA anchor checkpoint was not written")
}

rm(object_list)
invisible(gc())
future::plan("sequential")

message("Integrating data")
integrated <- IntegrateData(
  anchorset = anchors,
  new.assay.name = "integrated",
  normalization.method = "LogNormalize",
  dims = 1:20,
  k.weight = 100,
  sd.weight = 1,
  preserve.order = FALSE,
  eps = 0,
  verbose = TRUE
)
rm(anchors)
invisible(gc())

DefaultAssay(integrated) <- "integrated"
message("Scaling, PCA, neighbor graph, clustering, and UMAP")
integrated <- ScaleData(integrated, features = integration_features, verbose = TRUE)
integrated <- RunPCA(integrated, features = integration_features, npcs = 50, seed.use = 260824, verbose = FALSE)
integrated <- FindNeighbors(integrated, dims = 1:20, nn.method = "annoy", n.trees = 50, verbose = TRUE)
integrated <- FindClusters(
  integrated, graph.name = "integrated_snn", resolution = 0.1, algorithm = 1,
  n.start = 10, n.iter = 10, random.seed = 260824,
  cluster.name = "drugres_snn_res.0.1", verbose = TRUE
)
integrated <- FindClusters(
  integrated, graph.name = "integrated_snn", resolution = 0.6, algorithm = 1,
  n.start = 10, n.iter = 10, random.seed = 260824,
  cluster.name = "source_snn_res.0.6", verbose = TRUE
)
integrated <- RunUMAP(
  integrated, dims = 1:20, reduction = "pca", umap.method = "uwot",
  n.neighbors = 30L, min.dist = 0.3, metric = "cosine",
  seed.use = 260824, verbose = TRUE
)

if (!setequal(Cells(integrated), meta$cell_id) || ncol(integrated) != 92053L) stop("Integrated cell set changed")
if (!all(c("drugres_snn_res.0.1", "source_snn_res.0.6") %in% colnames(integrated[[]]))) {
  stop("Expected cluster columns are absent")
}

umap <- as.data.table(Embeddings(integrated, reduction = "umap"), keep.rownames = "cell_id")
cluster_meta <- as.data.table(integrated[[]], keep.rownames = "seurat_cell_id")
if (!"cell_id" %in% names(cluster_meta)) stop("Original cell_id metadata column is absent")
if (!identical(cluster_meta$seurat_cell_id, cluster_meta$cell_id)) stop("Seurat row names disagree with cell_id metadata")
cluster_meta[, cell_id := NULL]
setnames(cluster_meta, "seurat_cell_id", "cell_id")
out <- merge(cluster_meta, umap, by = "cell_id", sort = FALSE)
out <- out[match(Cells(integrated), cell_id)]
if (!identical(out$cell_id, Cells(integrated))) stop("Cluster/UMAP export order failed")
fwrite(out, file.path(table_dir, "03_integrated_clusters_and_umap.csv"))

for (cluster_col in c("drugres_snn_res.0.1", "source_snn_res.0.6")) {
  fwrite(out[, .N, by = cluster_col][order(get(cluster_col))],
         file.path(table_dir, paste0("03_cluster_sizes_", gsub("[^A-Za-z0-9]+", "_", cluster_col), ".csv")))
  fwrite(out[, .N, by = c(cluster_col, "Sample")][order(get(cluster_col), Sample)],
         file.path(table_dir, paste0("03_cluster_by_sample_", gsub("[^A-Za-z0-9]+", "_", cluster_col), ".csv")))
}

saveRDS(integrated, file.path(intermediate_dir, "03_integrated_cca_clustered.rds"), compress = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "03_sessionInfo.txt"))

cat("CCA integration and clustering PASS\n")
cat("Cells:", ncol(integrated), "\n")
cat("Primary resolution 0.1 clusters:", length(unique(integrated$drugres_snn_res.0.1)), "\n")
cat("Validation resolution 0.6 clusters:", length(unique(integrated$source_snn_res.0.6)), "\n")
cat("Selected-feature normalization max absolute difference:", format(normalization_max_abs_diff, scientific = TRUE), "\n")
