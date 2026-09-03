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
set.seed(260825)

table_dir <- file.path(analysis_dir, "results", "tables")
intermediate_dir <- file.path(analysis_dir, "intermediate")
log_dir <- file.path(analysis_dir, "logs")

counts_file <- file.path(intermediate_dir, "01_full_counts_mt20_qc.rds")
cluster_file <- file.path(table_dir, "03_integrated_clusters_and_umap.csv")
if (!all(file.exists(c(counts_file, cluster_file)))) stop("Required T/NK inputs are missing")

counts <- readRDS(counts_file)
meta <- fread(cluster_file)
if (!inherits(counts, "dgCMatrix")) stop("Counts are not a dgCMatrix")
if (nrow(counts) != 24292L || ncol(counts) != 92053L) stop("Unexpected full matrix dimensions")
if (nrow(meta) != ncol(counts) || anyDuplicated(meta$cell_id)) stop("Invalid all-cell metadata")
ord <- match(colnames(counts), meta$cell_id)
if (anyNA(ord)) stop("Counts/all-cell metadata mismatch")
meta <- meta[ord]
if (!identical(meta$cell_id, colnames(counts))) stop("All-cell metadata ordering failed")
if (!identical(as.numeric(meta$library_size), as.numeric(Matrix::colSums(counts)))) {
  stop("Full library sizes disagree with counts")
}

# These broad T/NK compartments were frozen from the authors-marker audit before
# running this script. Cluster 1 is deliberately retained as mixed cytotoxic
# CD8/NK so the source-method T/NK re-clustering can resolve it.
selection_map <- data.table(
  source_snn_res_0_6 = c("0", "1", "2", "3", "10", "14"),
  decision = "include",
  evidence = c(
    "CD3E/CD8A/GZMK/NKG7 cytotoxic lymphocytes",
    "mixed cytotoxic CD8/NK; retain for T/NK re-clustering",
    "CD3E/CD4/IL7R T cells",
    "activated CD4/Treg-like T cells",
    "CD3E/CD8A/GZMB/CXCL13 cytotoxic-exhausted T cells",
    "MKI67/TOP2A cycling T cells"
  )
)
fwrite(selection_map, file.path(table_dir, "05_tnk_input_cluster_decision.csv"))

keep <- as.character(meta$source_snn_res.0.6) %in% selection_map$source_snn_res_0_6
if (sum(keep) != 38730L) stop("Unexpected T/NK input cell count: ", sum(keep))
meta_tnk <- copy(meta[keep])
counts_tnk <- counts[, keep, drop = FALSE]
rm(counts, meta)
invisible(gc())
if (!identical(colnames(counts_tnk), meta_tnk$cell_id)) stop("T/NK counts/metadata ordering failed")

meta_tnk[, clinical_group := fifelse(
  Resource == "Pre-treatment biopsy", "TN",
  fifelse(Pathologic.Response %in% c("MPR", "pCR"), "MPR",
          fifelse(Pathologic.Response == "NMPR", "NMPR", NA_character_))
)]
if (anyNA(meta_tnk$clinical_group)) stop("Unresolved clinical group")
patient_groups <- unique(meta_tnk[, .(Sample, Patient, Resource, Pathologic.Response, clinical_group)])
if (nrow(patient_groups) != 15L || anyDuplicated(patient_groups$Patient)) stop("Expected 15 unique patients")
expected_groups <- c(MPR = 4L, NMPR = 8L, TN = 3L)
observed_groups <- table(patient_groups$clinical_group)
if (!identical(as.integer(observed_groups[names(expected_groups)]), as.integer(expected_groups))) {
  stop("Clinical group reconstruction does not reproduce TN/MPR/NMPR = 3/4/8")
}
fwrite(patient_groups[order(factor(clinical_group, levels = c("TN", "MPR", "NMPR")), Patient)],
       file.path(table_dir, "05_patient_clinical_groups.csv"))

sample_ids <- sort(unique(meta_tnk$Sample))
if (length(sample_ids) != 15L) stop("Expected all 15 samples in T/NK compartment")

message("Selecting 3,000 T/NK variable features per sample")
hvg_list <- setNames(vector("list", length(sample_ids)), sample_ids)
hvg_audit <- vector("list", length(sample_ids))
for (idx in seq_along(sample_ids)) {
  sid <- sample_ids[idx]
  cell_idx <- which(meta_tnk$Sample == sid)
  cells <- meta_tnk$cell_id[cell_idx]
  md <- as.data.frame(meta_tnk[cell_idx])
  rownames(md) <- md$cell_id
  message("T/NK VST: ", sid, " (", length(cells), " cells; ", idx, "/15)")
  obj <- CreateSeuratObject(
    counts = counts_tnk[, cells, drop = FALSE], meta.data = md,
    project = sid, min.cells = 0, min.features = 0
  )
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  vf <- VariableFeatures(obj)
  if (length(vf) != 3000L || anyDuplicated(vf) || !all(vf %in% rownames(counts_tnk))) {
    stop("Invalid T/NK variable features for ", sid)
  }
  hvg_list[[sid]] <- vf
  hvg_audit[[idx]] <- data.table(
    Sample = sid, cells = length(cells), genes = nrow(obj), variable_features = length(vf),
    full_library_sum = sum(meta_tnk$library_size[cell_idx]), object_library_sum = sum(obj$nCount_RNA)
  )
  if (hvg_audit[[idx]]$full_library_sum != hvg_audit[[idx]]$object_library_sum) {
    stop("T/NK library-size mismatch for ", sid)
  }
  rm(obj, md)
  invisible(gc())
}
fwrite(rbindlist(hvg_audit), file.path(table_dir, "05_tnk_sample_feature_selection_audit.csv"))
saveRDS(hvg_list, file.path(intermediate_dir, "05_tnk_sample_hvg_rankings.rds"), compress = "xz")

# Exact reimplementation of Seurat::SelectIntegrationFeatures ranking.
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
  stop("T/NK integration-feature ranking failed")
}
feature_support <- data.table(
  rank = seq_along(integration_features), gene = integration_features,
  samples_selected = as.integer(var_features[integration_features]),
  median_within_sample_rank = vapply(integration_features, median_rank, numeric(1))
)
fwrite(feature_support, file.path(table_dir, "05_tnk_integration_features_3000.csv"))

log_normalize_selected <- function(counts_selected, full_library_size, scale_factor = 10000) {
  if (!inherits(counts_selected, "dgCMatrix")) counts_selected <- as(counts_selected, "dgCMatrix")
  if (length(full_library_size) != ncol(counts_selected) || any(full_library_size <= 0)) {
    stop("Invalid T/NK full library sizes")
  }
  out <- counts_selected
  denom <- rep(full_library_size, diff(out@p))
  if (length(denom) != length(out@x)) stop("Sparse normalization index mismatch")
  out@x <- log1p(scale_factor * out@x / denom)
  out
}

# Unit check selected-feature normalization against full-universe Seurat.
test_cells <- meta_tnk$cell_id[seq_len(min(50L, nrow(meta_tnk)))]
test_full <- CreateSeuratObject(counts_tnk[, test_cells, drop = FALSE], min.cells = 0, min.features = 0)
test_full <- NormalizeData(test_full, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
test_manual <- log_normalize_selected(
  counts_tnk[integration_features, test_cells, drop = FALSE],
  meta_tnk$library_size[match(test_cells, meta_tnk$cell_id)]
)
test_reference <- LayerData(test_full, assay = "RNA", layer = "data")[integration_features, , drop = FALSE]
normalization_max_abs_diff <- max(abs(test_reference - test_manual))
if (!is.finite(normalization_max_abs_diff) || normalization_max_abs_diff > 1e-12) {
  stop("T/NK selected-feature normalization unit check failed")
}
rm(test_full, test_manual, test_reference)
invisible(gc())

object_list <- setNames(vector("list", length(sample_ids)), sample_ids)
object_audit <- vector("list", length(sample_ids))
for (idx in seq_along(sample_ids)) {
  sid <- sample_ids[idx]
  cell_idx <- which(meta_tnk$Sample == sid)
  cells <- meta_tnk$cell_id[cell_idx]
  md <- as.data.frame(meta_tnk[cell_idx])
  rownames(md) <- md$cell_id
  counts_selected <- counts_tnk[integration_features, cells, drop = FALSE]
  data_selected <- log_normalize_selected(counts_selected, meta_tnk$library_size[cell_idx])
  obj <- CreateSeuratObject(
    counts = counts_selected, meta.data = md, project = sid,
    min.cells = 0, min.features = 0
  )
  obj <- SetAssayData(obj, assay = "RNA", layer = "data", new.data = data_selected)
  VariableFeatures(obj) <- intersect(hvg_list[[sid]], integration_features)
  object_list[[sid]] <- obj
  object_audit[[idx]] <- data.table(
    Sample = sid, cells = ncol(obj), integration_features = nrow(obj),
    original_HVGs_in_integration_set = length(VariableFeatures(obj)),
    full_library_sum = sum(meta_tnk$library_size[cell_idx]),
    selected_feature_UMI_sum = sum(counts_selected)
  )
  rm(obj, md, counts_selected, data_selected)
  invisible(gc())
}
fwrite(rbindlist(object_audit), file.path(table_dir, "05_tnk_cca_input_object_audit.csv"))
fwrite(data.table(metric = "normalization_max_abs_diff", value = normalization_max_abs_diff),
       file.path(table_dir, "05_tnk_selected_feature_normalization_unit_check.csv"))

rm(counts_tnk, hvg_list)
invisible(gc())

# The smaller T/NK compartment is run sequentially to avoid exporting large
# Seurat globals while another independent analysis may share this workstation.
future::plan("sequential")
anchor_file <- file.path(intermediate_dir, "05_tnk_cca_anchors_checkpoint.rds")
if (file.exists(anchor_file)) {
  message("Loading T/NK CCA anchor checkpoint")
  anchors <- readRDS(anchor_file)
  if (!inherits(anchors, "AnchorSet")) stop("T/NK anchor checkpoint is invalid")
} else {
  message("Finding source-method T/NK CCA anchors")
  anchors <- FindIntegrationAnchors(
    object.list = object_list, anchor.features = integration_features,
    normalization.method = "LogNormalize", reduction = "cca", dims = 1:20,
    l2.norm = TRUE, k.anchor = 5, k.filter = 200, k.score = 30,
    max.features = 200, nn.method = "annoy", n.trees = 50, eps = 0,
    verbose = TRUE
  )
  saveRDS(anchors, anchor_file, compress = FALSE)
  if (!file.exists(anchor_file) || file.info(anchor_file)$size <= 0) {
    stop("T/NK anchor checkpoint was not written")
  }
}
rm(object_list)
invisible(gc())

message("Integrating T/NK compartment")
integrated <- IntegrateData(
  anchorset = anchors, new.assay.name = "integrated",
  normalization.method = "LogNormalize", dims = 1:20,
  k.weight = 100, sd.weight = 1, preserve.order = FALSE, eps = 0,
  verbose = TRUE
)
rm(anchors)
invisible(gc())

DefaultAssay(integrated) <- "integrated"
message("Scaling, PCA, source-method graph/clustering, and UMAP")
integrated <- ScaleData(integrated, features = integration_features, verbose = TRUE)
integrated <- RunPCA(integrated, features = integration_features, npcs = 50,
                     seed.use = 260825, verbose = FALSE)
integrated <- FindNeighbors(integrated, dims = 1:15, nn.method = "annoy",
                            n.trees = 50, verbose = TRUE)
integrated <- FindClusters(
  integrated, graph.name = "integrated_snn", resolution = 0.4,
  algorithm = 1, n.start = 10, n.iter = 10, random.seed = 260825,
  cluster.name = "tnk_source_snn_res.0.4", verbose = TRUE
)
integrated <- RunUMAP(
  integrated, dims = 1:15, reduction = "pca", umap.method = "uwot",
  n.neighbors = 30L, min.dist = 0.3, metric = "cosine",
  seed.use = 260825, verbose = TRUE
)

if (!setequal(Cells(integrated), meta_tnk$cell_id) || ncol(integrated) != 38730L) {
  stop("T/NK integrated cell set changed")
}
if (!"tnk_source_snn_res.0.4" %in% colnames(integrated[[]])) stop("T/NK cluster column absent")

umap <- as.data.table(Embeddings(integrated, reduction = "umap"), keep.rownames = "cell_id")
cluster_meta <- as.data.table(integrated[[]], keep.rownames = "seurat_cell_id")
if (!"cell_id" %in% names(cluster_meta)) stop("Original T/NK cell_id metadata column absent")
if (!identical(cluster_meta$seurat_cell_id, cluster_meta$cell_id)) stop("T/NK Seurat rownames/cell_id mismatch")
cluster_meta[, cell_id := NULL]
setnames(cluster_meta, "seurat_cell_id", "cell_id")
out <- merge(cluster_meta, umap, by = "cell_id", sort = FALSE)
out <- out[match(Cells(integrated), cell_id)]
if (!identical(out$cell_id, Cells(integrated))) stop("T/NK cluster/UMAP export order failed")

fwrite(out, file.path(table_dir, "05_tnk_clusters_and_umap.csv"))
fwrite(out[, .N, by = `tnk_source_snn_res.0.4`][order(`tnk_source_snn_res.0.4`)],
       file.path(table_dir, "05_tnk_cluster_sizes.csv"))
fwrite(out[, .N, by = .(`tnk_source_snn_res.0.4`, Sample, clinical_group)][order(`tnk_source_snn_res.0.4`, Sample)],
       file.path(table_dir, "05_tnk_cluster_by_sample_and_group.csv"))
saveRDS(integrated, file.path(intermediate_dir, "05_tnk_integrated_cca_clustered.rds"), compress = FALSE)
capture.output(sessionInfo(), file = file.path(log_dir, "05_sessionInfo.txt"))

cat("Source-method T/NK re-clustering PASS\n")
cat("Cells:", ncol(integrated), "\n")
cat("Samples:", length(unique(out$Sample)), "\n")
cat("Clinical groups (patients): TN/MPR/NMPR = 3/4/8\n")
cat("Resolution 0.4 clusters:", length(unique(out$`tnk_source_snn_res.0.4`)), "\n")
cat("Selected-feature normalization max absolute difference:",
    format(normalization_max_abs_diff, scientific = TRUE), "\n")
