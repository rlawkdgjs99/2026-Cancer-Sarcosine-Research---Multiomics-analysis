#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
.libPaths(c(file.path(analysis_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages(library(data.table))

table_dir <- file.path(analysis_dir, "results", "tables")
log_dir <- file.path(analysis_dir, "logs")
all_cells <- fread(file.path(table_dir, "03_integrated_clusters_and_umap.csv"))
tnk <- fread(file.path(table_dir, "05_tnk_clusters_and_umap.csv"))
memory_t <- fread(file.path(table_dir, "07_memory_t_clusters_and_umap.csv"))
if (nrow(all_cells) != 92053L || nrow(tnk) != 38730L || nrow(memory_t) != 7993L) {
  stop("Unexpected lineage-mapping input size")
}
if (anyDuplicated(all_cells$cell_id) || anyDuplicated(tnk$cell_id) || anyDuplicated(memory_t$cell_id)) {
  stop("Duplicated cell IDs in lineage-mapping inputs")
}

all_map <- data.table(
  source_snn_res_0_6 = as.character(0:22),
  broad_decision = c(
    "T/NK_recluster", "T/NK_recluster", "T/NK_recluster", "T/NK_recluster",
    "B cell", "Epithelial", "Neutrophil", "Macrophage", "Macrophage", "Monocyte",
    "T/NK_recluster", "Plasma cell", "Neutrophil", "Epithelial", "T/NK_recluster",
    "Epithelial", "CAF", "Conventional DC", "Mast cell", "Macrophage",
    "Epithelial", "Epithelial", "pDC"
  ),
  evidence = c(
    rep("resolved by source-method T/NK re-clustering", 4),
    "MS4A1/CD79A", "EPCAM/KRT epithelial program", "CSF3R/FCGR3B/S100A8-A9",
    "APOC1/GPNMB/SPP1 macrophage", "C1QC macrophage", "VCAN/IL1B/CD14 monocyte",
    "resolved by source-method T/NK re-clustering", "JCHAIN/IGHG1/XBP1/MZB1",
    "CSF3R/FCGR3B/S100A8-A9", "EPCAM/KRT epithelial program",
    "resolved by source-method T/NK re-clustering", "EPCAM/KRT epithelial program",
    "COL1A1/COL3A1/DCN/LUM", "CD1C/FCER1A/LAMP3/FSCN1",
    "TPSB2/CPA3/TPSAB1/KIT", "C1Q/SPP1/MKI67 cycling macrophage",
    "SFTPA/B alveolar epithelial", "ciliated epithelial program",
    "LILRA4/GZMB/IRF8/GPR183/TCF4"
  )
)

tnk_map <- data.table(
  tnk_source_snn_res_0_4 = as.character(0:8),
  tnk_decision = c(
    "CD8 T cell", "Memory_T_recluster", "CD8 T cell", "CD4 T cell", "CD4 T cell",
    "NK cell", "Cycling T cell", "NK/gamma-delta T unresolved", "Excluded residual doublet"
  ),
  evidence = c(
    "CD3E/CD8A/GZMK", "IL7R/TCF7 with mixed CD4 and CD8; resolved below",
    "CD3E/CD8A/GZMB/HAVCR2/CXCL13", "CD3E/CD4/MAF/CXCL13",
    "CD3E/CD4/FOXP3/IL2RA/CTLA4", "FCGR3A/FGFBP2/KLRD1/GNLY with low CD3D",
    "MKI67/STMN1/TOP2A T-lineage", "TRDC/TRGC with NK program; no defensible binary label",
    "CD3/TCR plus LYZ/FCER1G/CSF3R incompatible programs"
  ),
  analysis_eligible = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
)

memory_map <- data.table(
  memory_t_snn_res_0_4 = as.character(0:3),
  memory_t_decision = c("CD4 T cell", "CD4 T cell", "CD8 T cell", "CD4 T cell"),
  evidence = c(
    "CD4/IL7R activated memory", "CD4/IL7R/CCR7/TCF7 naive-memory",
    "CD8A/CD8B/CCL5/NKG7/GZMK memory", "15-cell IL7R stress cluster; CD4>T8 marker support"
  )
)

fwrite(all_map, file.path(table_dir, "09_all_cell_cluster_mapping_FROZEN.csv"))
fwrite(tnk_map, file.path(table_dir, "09_tnk_cluster_mapping_FROZEN.csv"))
fwrite(memory_map, file.path(table_dir, "09_memory_t_cluster_mapping_FROZEN.csv"))

out <- copy(all_cells)
out[, source_cluster := as.character(source_snn_res.0.6)]
out <- merge(out, all_map, by.x = "source_cluster", by.y = "source_snn_res_0_6", all.x = TRUE, sort = FALSE)
if (anyNA(out$broad_decision)) stop("Unmapped all-cell cluster")
out[, `:=`(
  final_lineage = broad_decision,
  annotation_level = "all-cell resolution 0.6",
  annotation_evidence = evidence,
  analysis_eligible = TRUE
)]
out[, evidence := NULL]

tnk_small <- tnk[, .(
  cell_id,
  tnk_cluster = as.character(`tnk_source_snn_res.0.4`)
)]
tnk_small <- merge(tnk_small, tnk_map, by.x = "tnk_cluster", by.y = "tnk_source_snn_res_0_4",
                   all.x = TRUE, sort = FALSE)
if (anyNA(tnk_small$tnk_decision)) stop("Unmapped T/NK cluster")
setkey(out, cell_id)
setkey(tnk_small, cell_id)
tnk_idx <- which(out$broad_decision == "T/NK_recluster")
if (length(tnk_idx) != 38730L || anyNA(tnk_small[out[tnk_idx], on = "cell_id"]$tnk_decision)) {
  stop("T/NK/all-cell mapping mismatch")
}
tnk_join <- tnk_small[out[tnk_idx], on = "cell_id"]
out[tnk_idx, `:=`(
  tnk_cluster = tnk_join$tnk_cluster,
  final_lineage = tnk_join$tnk_decision,
  annotation_level = "T/NK resolution 0.4",
  annotation_evidence = tnk_join$evidence,
  analysis_eligible = tnk_join$analysis_eligible
)]

memory_small <- memory_t[, .(
  cell_id,
  memory_t_cluster = as.character(`memory_t_snn_res.0.4`)
)]
memory_small <- merge(memory_small, memory_map,
                      by.x = "memory_t_cluster", by.y = "memory_t_snn_res_0_4",
                      all.x = TRUE, sort = FALSE)
if (anyNA(memory_small$memory_t_decision)) stop("Unmapped memory-T cluster")
setkey(memory_small, cell_id)
memory_idx <- which(out$final_lineage == "Memory_T_recluster")
if (length(memory_idx) != 7993L) stop("Unexpected memory-T mapping size")
memory_join <- memory_small[out[memory_idx], on = "cell_id"]
if (anyNA(memory_join$memory_t_decision)) stop("Memory-T/all-cell mapping mismatch")
out[memory_idx, `:=`(
  memory_t_cluster = memory_join$memory_t_cluster,
  final_lineage = memory_join$memory_t_decision,
  annotation_level = "memory-T resolution 0.4 reconstruction",
  annotation_evidence = memory_join$evidence,
  analysis_eligible = TRUE
)]

out[, broad_decision := NULL]
if (nrow(out) != 92053L || anyDuplicated(out$cell_id) || anyNA(out$final_lineage) ||
    any(out$final_lineage == "T/NK_recluster") || any(out$final_lineage == "Memory_T_recluster")) {
  stop("Invalid final cell-lineage table")
}
if (sum(out$final_lineage == "Excluded residual doublet") != 209L) stop("Unexpected residual-doublet count")
if (sum(out$final_lineage == "NK/gamma-delta T unresolved") != 1332L) stop("Unexpected unresolved T/NK count")

out[, original_all_cell_order := match(cell_id, all_cells$cell_id)]
if (anyNA(out$original_all_cell_order) || anyDuplicated(out$original_all_cell_order)) {
  stop("Could not reconstruct original all-cell order")
}
setorder(out, original_all_cell_order)
out[, original_all_cell_order := NULL]
if (!identical(out$cell_id, all_cells$cell_id)) stop("Final lineage order differs from all-cell order")
fwrite(out, file.path(table_dir, "09_final_cell_lineages_FROZEN.csv"))

lineage_counts <- out[, .N, by = .(final_lineage, analysis_eligible)][order(-N)]
sample_counts <- out[, .N, by = .(final_lineage, analysis_eligible, Sample, Patient,
                                  Resource, Pathologic.Response)][order(final_lineage, Sample)]
sample_counts[, clinical_group := fifelse(
  Resource == "Pre-treatment biopsy", "TN",
  fifelse(Pathologic.Response %in% c("MPR", "pCR"), "MPR", "NMPR")
)]
fwrite(lineage_counts, file.path(table_dir, "09_final_lineage_cell_counts.csv"))
fwrite(sample_counts, file.path(table_dir, "09_final_lineage_by_patient_counts.csv"))
capture.output(sessionInfo(), file = file.path(log_dir, "09_sessionInfo.txt"))

cat("Final lineage mapping FROZEN and PASS\n")
cat("Cells:", nrow(out), "\n")
cat("Final labels:", uniqueN(out$final_lineage), "\n")
cat("Analysis-eligible cells:", sum(out$analysis_eligible), "\n")
cat("Excluded residual doublets:", sum(out$final_lineage == "Excluded residual doublet"), "\n")
cat("Unresolved NK/gamma-delta T:", sum(out$final_lineage == "NK/gamma-delta T unresolved"), "\n")
