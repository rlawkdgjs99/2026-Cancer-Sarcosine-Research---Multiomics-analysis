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
})

options(stringsAsFactors = FALSE, warn = 1)
set.seed(260824)

table_dir <- file.path(analysis_dir, "results", "tables")
intermediate_dir <- file.path(analysis_dir, "intermediate")
log_dir <- file.path(analysis_dir, "logs")

counts_file <- file.path(intermediate_dir, "01_full_counts_mt20_qc.rds")
meta_file <- file.path(table_dir, "01_cell_metadata_mt20_qc.csv")
if (!all(file.exists(c(counts_file, meta_file)))) stop("QC-passed inputs are missing")

counts <- readRDS(counts_file)
meta <- fread(meta_file)
if (!inherits(counts, "dgCMatrix")) stop("Counts are not a dgCMatrix")
if (!identical(colnames(counts), meta$cell_id)) stop("Counts/metadata alignment failed")
if (ncol(counts) != 92053L || nrow(counts) != 24292L) stop("Unexpected QC-passed dimensions")

sample_ids <- sort(unique(meta$Sample))
if (length(sample_ids) != 15L) stop("Expected 15 samples")
hvg_list <- setNames(vector("list", length(sample_ids)), sample_ids)
sample_audit <- vector("list", length(sample_ids))

for (idx in seq_along(sample_ids)) {
  sid <- sample_ids[idx]
  cell_idx <- which(meta$Sample == sid)
  cells <- meta$cell_id[cell_idx]
  md <- as.data.frame(meta[cell_idx])
  rownames(md) <- md$cell_id
  message("VST feature selection: ", sid, " (", length(cells), " cells; ", idx, "/", length(sample_ids), ")")

  obj <- CreateSeuratObject(
    counts = counts[, cells, drop = FALSE],
    meta.data = md,
    project = sid,
    min.cells = 0,
    min.features = 0
  )
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  vf <- VariableFeatures(obj)
  if (length(vf) != 3000L || anyDuplicated(vf) || !all(vf %in% rownames(counts))) {
    stop("Invalid variable-feature result for ", sid)
  }
  hvg_list[[sid]] <- vf
  sample_audit[[idx]] <- data.table(
    Sample = sid,
    cells = length(cells),
    genes = nrow(obj),
    variable_features = length(vf),
    full_library_sum = sum(meta$library_size[cell_idx]),
    object_library_sum = sum(obj$nCount_RNA)
  )
  if (sample_audit[[idx]]$full_library_sum != sample_audit[[idx]]$object_library_sum) {
    stop("Library-size mismatch in ", sid)
  }
  rm(obj, md)
  invisible(gc())
}

hvg_table <- rbindlist(lapply(names(hvg_list), function(sid) {
  data.table(Sample = sid, rank = seq_along(hvg_list[[sid]]), gene = hvg_list[[sid]])
}))
fwrite(hvg_table, file.path(table_dir, "02_sample_hvg_rankings.csv"))
fwrite(rbindlist(sample_audit), file.path(table_dir, "02_sample_feature_selection_audit.csv"))

# Exact reimplementation of Seurat::SelectIntegrationFeatures ranking logic.
var_features <- sort(table(unname(unlist(hvg_list))), decreasing = TRUE)
tie_val <- var_features[min(3000L, length(var_features))]
features <- names(var_features[var_features > tie_val])
median_rank <- function(g) {
  ranks <- unlist(lapply(hvg_list, function(vf) if (g %in% vf) match(g, vf) else NULL), use.names = FALSE)
  median(ranks)
}
if (length(features)) {
  ranks <- vapply(features, median_rank, numeric(1))
  features <- names(sort(ranks))
}
features_tie <- names(var_features[var_features == tie_val])
tie_ranks <- vapply(features_tie, median_rank, numeric(1))
integration_features <- c(features, names(head(sort(tie_ranks), 3000L - length(features))))
if (length(integration_features) != 3000L || anyDuplicated(integration_features)) {
  stop("Integration-feature ranking failed")
}

# Unit check against Seurat's installed implementation using tiny one-cell objects
# with the same gene universe and the observed ranked VariableFeatures.
dummy_list <- lapply(seq_along(sample_ids), function(idx) {
  dummy_counts <- sparseMatrix(
    i = seq_len(nrow(counts)), j = rep.int(1L, nrow(counts)), x = rep.int(1, nrow(counts)),
    dims = c(nrow(counts), 1L), dimnames = list(rownames(counts), paste0("dummy_", idx))
  )
  obj <- CreateSeuratObject(dummy_counts, min.cells = 0, min.features = 0)
  VariableFeatures(obj) <- hvg_list[[idx]]
  obj
})
seurat_features <- SelectIntegrationFeatures(dummy_list, nfeatures = 3000, verbose = FALSE)
if (!identical(integration_features, seurat_features)) stop("Manual integration-feature ranking disagrees with Seurat")
rm(dummy_list, seurat_features)

feature_support <- data.table(
  rank = seq_along(integration_features),
  gene = integration_features,
  samples_selected = as.integer(var_features[integration_features]),
  median_within_sample_rank = vapply(integration_features, median_rank, numeric(1))
)
fwrite(feature_support, file.path(table_dir, "02_integration_features_3000.csv"))
saveRDS(hvg_list, file.path(intermediate_dir, "02_sample_hvg_rankings.rds"), compress = "xz")

authors_all_cell_markers <- c(
  "PTPRC", "CD3E", "LYZ", "CD79A", "MS4A1", "IGHG1", "CSF3R",
  "FGFBP2", "KIT", "LILRA4", "VWF", "COL1A1", "EPCAM", "MKI67"
)
authors_tnk_markers <- c(
  "CD8A", "CD4", "FCGR3A", "CTLA4", "PDCD1", "TIGIT", "HAVCR2", "GZMA",
  "GZMB", "GZMK", "NKG7", "SELL", "TCF7", "FOXP3", "ZNF683", "ITGAE"
)
target_audit <- data.table(
  gene = unique(c(authors_all_cell_markers, authors_tnk_markers, "SARDH", "PIPOX")),
  role = c(
    rep("authors_all_cell_marker", length(authors_all_cell_markers)),
    rep("authors_T_NK_marker", length(authors_tnk_markers)),
    "sarcosine_degradation_target", "sarcosine_degradation_target"
  )
)
target_audit[, present_in_deposited_matrix := gene %in% rownames(counts)]
target_audit[, selected_as_integration_feature := gene %in% integration_features]
fwrite(target_audit, file.path(table_dir, "02_marker_and_target_gene_audit.csv"))
if (!all(target_audit$present_in_deposited_matrix)) stop("One or more authors' marker/target genes are absent")

capture.output(sessionInfo(), file = file.path(log_dir, "02_sessionInfo.txt"))
cat("Integration-feature selection PASS\n")
cat("Samples:", length(sample_ids), "\n")
cat("Features:", length(integration_features), "\n")
cat("Exact SelectIntegrationFeatures unit check: PASS\n")
