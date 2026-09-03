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

options(stringsAsFactors = FALSE, warn = 1, future.globals.maxSize = 8 * 1024^3)
set.seed(260825)

table_dir <- file.path(analysis_dir, "results", "tables")
intermediate_dir <- file.path(analysis_dir, "intermediate")
log_dir <- file.path(analysis_dir, "logs")
counts_file <- file.path(intermediate_dir, "01_full_counts_mt20_qc.rds")
tnk_file <- file.path(table_dir, "05_tnk_clusters_and_umap.csv")
if (!all(file.exists(c(counts_file, tnk_file)))) stop("Required memory-T inputs are missing")

counts <- readRDS(counts_file)
meta <- fread(tnk_file)
if (!inherits(counts, "dgCMatrix")) stop("Counts are not a dgCMatrix")
if (nrow(meta) != 38730L || anyDuplicated(meta$cell_id)) stop("Invalid T/NK metadata")

# Cluster 1 is the reconstructed counterpart of the authors' mixed T_IL7R
# cluster: IL7R/TCF7-high, with both CD4- and CD8-marker-positive cells. The
# private Tmem.rds creation code is unavailable, so this is explicitly a
# reconstruction using the source paper's reported T-cell parameters.
meta <- meta[as.character(`tnk_source_snn_res.0.4`) == "1"]
if (nrow(meta) != 7993L) stop("Unexpected mixed memory-T input size: ", nrow(meta))
idx <- match(meta$cell_id, colnames(counts))
if (anyNA(idx)) stop("Memory-T cells absent from count matrix")
counts <- counts[, idx, drop = FALSE]
if (!identical(colnames(counts), meta$cell_id)) stop("Memory-T count/metadata ordering failed")
if (!identical(as.numeric(Matrix::colSums(counts)), as.numeric(meta$library_size))) {
  stop("Memory-T library sizes disagree with counts")
}

sample_ids <- sort(unique(meta$Sample))
if (length(sample_ids) != 15L) stop("Expected 15 samples in mixed memory-T cluster")
hvg_list <- setNames(vector("list", length(sample_ids)), sample_ids)
hvg_audit <- vector("list", length(sample_ids))

message("Selecting mixed memory-T variable features")
for (i in seq_along(sample_ids)) {
  sid <- sample_ids[i]
  cell_idx <- which(meta$Sample == sid)
  cells <- meta$cell_id[cell_idx]
  md <- as.data.frame(meta[cell_idx])
  rownames(md) <- md$cell_id
  message("Memory-T VST: ", sid, " (", length(cells), " cells; ", i, "/15)")
  obj <- CreateSeuratObject(counts[, cells, drop = FALSE], meta.data = md,
                            project = sid, min.cells = 0, min.features = 0)
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  vf <- VariableFeatures(obj)
  if (length(vf) != 3000L || anyDuplicated(vf)) stop("Invalid memory-T HVGs for ", sid)
  hvg_list[[sid]] <- vf
  hvg_audit[[i]] <- data.table(
    Sample = sid, cells = length(cells), variable_features = length(vf),
    full_library_sum = sum(meta$library_size[cell_idx]), object_library_sum = sum(obj$nCount_RNA)
  )
  if (hvg_audit[[i]]$full_library_sum != hvg_audit[[i]]$object_library_sum) {
    stop("Memory-T library mismatch for ", sid)
  }
  rm(obj, md)
  invisible(gc())
}
fwrite(rbindlist(hvg_audit), file.path(table_dir, "07_memory_t_sample_feature_selection_audit.csv"))
saveRDS(hvg_list, file.path(intermediate_dir, "07_memory_t_sample_hvg_rankings.rds"), compress = "xz")

var_features <- sort(table(unname(unlist(hvg_list))), decreasing = TRUE)
tie_val <- var_features[min(3000L, length(var_features))]
integration_features <- names(var_features[var_features > tie_val])
median_rank <- function(g) {
  ranks <- unlist(lapply(hvg_list, function(vf) if (g %in% vf) match(g, vf) else NULL), use.names = FALSE)
  median(ranks)
}
if (length(integration_features)) {
  ranks <- vapply(integration_features, median_rank, numeric(1))
  integration_features <- names(sort(ranks))
}
tie_features <- names(var_features[var_features == tie_val])
tie_ranks <- vapply(tie_features, median_rank, numeric(1))
integration_features <- c(
  integration_features,
  names(head(sort(tie_ranks), 3000L - length(integration_features)))
)
if (length(integration_features) != 3000L || anyDuplicated(integration_features)) {
  stop("Memory-T integration-feature ranking failed")
}
fwrite(data.table(
  rank = seq_along(integration_features), gene = integration_features,
  samples_selected = as.integer(var_features[integration_features]),
  median_within_sample_rank = vapply(integration_features, median_rank, numeric(1))
), file.path(table_dir, "07_memory_t_integration_features_3000.csv"))

log_normalize_selected <- function(x, full_library_size) {
  if (!inherits(x, "dgCMatrix")) x <- as(x, "dgCMatrix")
  out <- x
  denom <- rep(full_library_size, diff(out@p))
  if (length(denom) != length(out@x) || any(full_library_size <= 0)) {
    stop("Invalid memory-T normalization denominator")
  }
  out@x <- log1p(10000 * out@x / denom)
  out
}

test_cells <- meta$cell_id[seq_len(50L)]
test_full <- CreateSeuratObject(counts[, test_cells, drop = FALSE], min.cells = 0, min.features = 0)
test_full <- NormalizeData(test_full, verbose = FALSE)
test_manual <- log_normalize_selected(
  counts[integration_features, test_cells, drop = FALSE],
  meta$library_size[match(test_cells, meta$cell_id)]
)
test_reference <- LayerData(test_full, assay = "RNA", layer = "data")[integration_features, , drop = FALSE]
normalization_max_abs_diff <- max(abs(test_reference - test_manual))
if (!is.finite(normalization_max_abs_diff) || normalization_max_abs_diff > 1e-12) {
  stop("Memory-T selected-feature normalization unit check failed")
}
rm(test_full, test_manual, test_reference)
invisible(gc())

object_list <- setNames(vector("list", length(sample_ids)), sample_ids)
for (i in seq_along(sample_ids)) {
  sid <- sample_ids[i]
  cell_idx <- which(meta$Sample == sid)
  cells <- meta$cell_id[cell_idx]
  md <- as.data.frame(meta[cell_idx])
  rownames(md) <- md$cell_id
  count_selected <- counts[integration_features, cells, drop = FALSE]
  data_selected <- log_normalize_selected(count_selected, meta$library_size[cell_idx])
  obj <- CreateSeuratObject(count_selected, meta.data = md, project = sid,
                            min.cells = 0, min.features = 0)
  obj <- SetAssayData(obj, assay = "RNA", layer = "data", new.data = data_selected)
  VariableFeatures(obj) <- intersect(hvg_list[[sid]], integration_features)
  object_list[[sid]] <- obj
  rm(obj, md, count_selected, data_selected)
  invisible(gc())
}
rm(counts, hvg_list)
invisible(gc())

future::plan("sequential")
anchor_file <- file.path(intermediate_dir, "07_memory_t_cca_anchors_checkpoint.rds")
if (file.exists(anchor_file)) {
  anchors <- readRDS(anchor_file)
  if (!inherits(anchors, "AnchorSet")) stop("Invalid memory-T anchor checkpoint")
} else {
  message("Finding mixed memory-T CCA anchors")
  anchors <- FindIntegrationAnchors(
    object.list = object_list, anchor.features = integration_features,
    normalization.method = "LogNormalize", reduction = "cca", dims = 1:20,
    l2.norm = TRUE, k.anchor = 5, k.filter = 200, k.score = 30,
    max.features = 200, nn.method = "annoy", n.trees = 50, eps = 0,
    verbose = TRUE
  )
  saveRDS(anchors, anchor_file, compress = FALSE)
}
rm(object_list)
invisible(gc())

message("Integrating and clustering mixed memory-T cells")
integrated <- IntegrateData(
  anchorset = anchors, new.assay.name = "integrated",
  normalization.method = "LogNormalize", dims = 1:20,
  k.weight = 100, sd.weight = 1, preserve.order = FALSE, eps = 0,
  verbose = TRUE
)
rm(anchors)
invisible(gc())
DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, features = integration_features, verbose = FALSE)
integrated <- RunPCA(integrated, features = integration_features, npcs = 50,
                     seed.use = 260825, verbose = FALSE)
integrated <- FindNeighbors(integrated, dims = 1:15, nn.method = "annoy",
                            n.trees = 50, verbose = FALSE)
integrated <- FindClusters(
  integrated, graph.name = "integrated_snn", resolution = 0.4,
  algorithm = 1, n.start = 10, n.iter = 10, random.seed = 260825,
  cluster.name = "memory_t_snn_res.0.4", verbose = TRUE
)
integrated <- RunUMAP(
  integrated, dims = 1:15, reduction = "pca", umap.method = "uwot",
  n.neighbors = 30L, min.dist = 0.3, metric = "cosine",
  seed.use = 260825, verbose = FALSE
)
if (!setequal(Cells(integrated), meta$cell_id) || ncol(integrated) != 7993L) {
  stop("Memory-T integrated cell set changed")
}

umap <- as.data.table(Embeddings(integrated, "umap"), keep.rownames = "cell_id")
cluster_meta <- as.data.table(integrated[[]], keep.rownames = "seurat_cell_id")
if (!identical(cluster_meta$seurat_cell_id, cluster_meta$cell_id)) {
  stop("Memory-T Seurat rownames/cell_id mismatch")
}
cluster_meta[, cell_id := NULL]
setnames(cluster_meta, "seurat_cell_id", "cell_id")
out <- merge(cluster_meta, umap, by = "cell_id", sort = FALSE)
out <- out[match(Cells(integrated), cell_id)]
if (!identical(out$cell_id, Cells(integrated))) stop("Memory-T export ordering failed")
fwrite(out, file.path(table_dir, "07_memory_t_clusters_and_umap.csv"))
fwrite(out[, .N, by = `memory_t_snn_res.0.4`][order(`memory_t_snn_res.0.4`)],
       file.path(table_dir, "07_memory_t_cluster_sizes.csv"))
fwrite(out[, .N, by = .(`memory_t_snn_res.0.4`, Sample, clinical_group)][order(`memory_t_snn_res.0.4`, Sample)],
       file.path(table_dir, "07_memory_t_cluster_by_sample_and_group.csv"))
saveRDS(integrated, file.path(intermediate_dir, "07_memory_t_integrated_cca_clustered.rds"), compress = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "07_sessionInfo.txt"))

cat("Mixed memory-T source-method reconstruction PASS\n")
cat("Cells:", ncol(integrated), "\n")
cat("Clusters at resolution 0.4:", length(unique(out$`memory_t_snn_res.0.4`)), "\n")
cat("Normalization max absolute difference:", format(normalization_max_abs_diff, scientific = TRUE), "\n")
