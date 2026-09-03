#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
analysis_root <- normalizePath(file.path(analysis_dir, ".."), mustWork = TRUE)
data_root <- normalizePath(file.path(analysis_root, ".."), mustWork = TRUE)
prior_dir <- file.path(analysis_root, "02_lineage_reannotation")
.libPaths(c(file.path(prior_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

options(stringsAsFactors = FALSE, warn = 1)

table_dir <- file.path(analysis_dir, "results", "tables")
log_dir <- file.path(analysis_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  raw_count_matrix = file.path(data_root, "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz"),
  raw_metadata = file.path(data_root, "GSE207422_NSCLC_scRNAseq_metadata.xlsx"),
  sparse_counts_mt20 = file.path(prior_dir, "intermediate", "01_full_counts_mt20_qc.rds"),
  integrated_umap = file.path(prior_dir, "results", "tables", "03_integrated_clusters_and_umap.csv"),
  frozen_lineages = file.path(prior_dir, "results", "tables", "09_final_cell_lineages_FROZEN.csv"),
  gene_audit = file.path(prior_dir, "results", "tables", "01_gene_qc_full_matrix.csv"),
  lee_main_pdf = file.path(data_root, "..", "1-s2.0-S1368764624001171-main.pdf"),
  lee_supp_methods = file.path(data_root, "..", "Supp_Material.docx"),
  frozen_plan = file.path(analysis_dir, "ANALYSIS_PLAN_FROZEN.md")
)
if (!all(file.exists(unlist(paths)))) {
  stop("Missing input(s): ", paste(names(paths)[!file.exists(unlist(paths))], collapse = ", "))
}

expected <- c(
  raw_count_matrix = "aba15960fc7ee6a2443511bce5177e4d71b964131b6e98597e7a85d0a213ba36",
  raw_metadata = "d098a750c7ebc595994c929b666177f25c4da6fe3d3e38bb795269a1fb21053e",
  sparse_counts_mt20 = "7552ba882c306fa68380390a8172e0fc96af2c7be7cb3f5778345204e9fd61b2",
  integrated_umap = "6bee4b491ba9f621a1df0231b91ec3dbf467e404b2ac4b0143d2d59787d5f569",
  frozen_lineages = "97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1",
  lee_main_pdf = "7a7e11c364d850b54ae6f083887757fd7def0ad3e17f37556d66e8e511fe475f",
  lee_supp_methods = "67b301ce1712a4ce407b05d1fee87f627961125f37e273382a4486389d1fe2a6",
  frozen_plan = "6d057689efc65f0b262d70ce12fa471aabac51028de70bfdb62740fa2c550c2f"
)

manifest <- rbindlist(lapply(names(paths), function(role) {
  p <- paths[[role]]
  data.table(
    role = role,
    path = normalizePath(p),
    bytes = file.info(p)$size,
    sha256 = digest(p, algo = "sha256", file = TRUE),
    expected_sha256 = expected[role]
  )
}))
manifest[, checksum_status := fifelse(
  is.na(expected_sha256), "RECORDED_NOT_PREDECLARED",
  fifelse(sha256 == expected_sha256, "MATCH", "MISMATCH")
)]
if (manifest[!is.na(expected_sha256), any(checksum_status != "MATCH")]) {
  stop("Immutable input or plan checksum mismatch")
}
fwrite(manifest, file.path(table_dir, "00_input_manifest.csv"))

lin <- fread(paths$frozen_lineages)
required <- c(
  "cell_id", "Sample", "Patient", "Resource", "Pathologic.Response",
  "umap_1", "umap_2", "final_lineage", "analysis_eligible"
)
if (!all(required %in% names(lin))) stop("Frozen lineage table lacks required columns")
if (nrow(lin) != 92053L || uniqueN(lin$cell_id) != 92053L || anyNA(lin$cell_id)) {
  stop("Frozen cell table failed its 92,053-cell contract")
}
if (nrow(lin[analysis_eligible == TRUE]) != 90512L ||
    uniqueN(lin[analysis_eligible == TRUE]$final_lineage) != 14L) {
  stop("Expected 90,512 cells in 14 marker-audited lineages")
}
if (uniqueN(lin$Patient) != 15L || uniqueN(lin$Sample) != 15L) {
  stop("Expected 15 patients and samples")
}

umap <- fread(paths$integrated_umap)
if (nrow(umap) != 92053L || uniqueN(umap$cell_id) != 92053L) {
  stop("Integrated UMAP table failed its cell contract")
}
if (!identical(lin$cell_id, umap$cell_id)) stop("Frozen lineage/UMAP cell order mismatch")
if (max(abs(lin$umap_1 - umap$umap_1)) > 1e-12 ||
    max(abs(lin$umap_2 - umap$umap_2)) > 1e-12) {
  stop("Frozen lineage/UMAP coordinates differ")
}

genes <- fread(paths$gene_audit)
targets <- c("GNMT", "DMGDH", "SARDH", "PIPOX")
target_audit <- genes[gene %in% targets]
if (nrow(target_audit) != 4L || !setequal(target_audit$gene, targets) || anyDuplicated(target_audit$gene)) {
  stop("The four prespecified genes are not present exactly once")
}
target_audit[, module := fifelse(gene %in% c("GNMT", "DMGDH"), "Production", "Degradation")]
setorder(target_audit, module, gene)
fwrite(target_audit, file.path(table_dir, "00_target_gene_audit.csv"))

packages <- c(
  "Seurat", "SeuratObject", "data.table", "digest", "Matrix", "ggplot2",
  "patchwork", "scales", "limma", "AnnotationDbi", "org.Hs.eg.db"
)
pkg <- data.table(
  package = packages,
  available = vapply(packages, requireNamespace, logical(1), quietly = TRUE)
)
pkg[, version := vapply(package, function(x) {
  if (!requireNamespace(x, quietly = TRUE)) return(NA_character_)
  as.character(packageVersion(x))
}, character(1))]
if (pkg[, any(!available)]) stop("Required package unavailable: ", paste(pkg[!available]$package, collapse = ", "))
fwrite(pkg, file.path(table_dir, "00_package_versions.csv"))

summary <- data.table(
  metric = c(
    "cells_mt_le_20", "published_Fig3_cells", "difference_from_published",
    "analysis_eligible_cells", "marker_audited_lineages", "patients",
    "production_genes", "degradation_genes"
  ),
  value = c(
    "92053", "92031", "22", "90512", "14", "15",
    "GNMT;DMGDH", "SARDH;PIPOX"
  )
)
fwrite(summary, file.path(table_dir, "00_preflight_summary.csv"))

method_contract <- data.table(
  step = c(
    "normalization", "module scoring", "minimum shift", "ratio",
    "undefined denominator", "median high/low", "DEG", "enrichment"
  ),
  frozen_implementation = c(
    "Seurat LogNormalize scale.factor=10000",
    "one AddModuleScore call; nbin=24; ctrl=100; seed=260825",
    "score minus global minimum, separately per module",
    "shifted Production divided by shifted Degradation",
    "exclude D0==0 only from ratio-specific analyses; no epsilon",
    "High>finite median; Low<finite median; exact ties excluded",
    "FindMarkers Wilcoxon; logfc.threshold=0.1; min.pct=0.01; no downsampling",
    "GO BP over-representation via limma::goana; BH by direction"
  )
)
fwrite(method_contract, file.path(table_dir, "00_method_contract.csv"))

capture.output(sessionInfo(), file = file.path(log_dir, "00_sessionInfo.txt"))
cat("LEE FIGURE 3-STYLE PRE-FLIGHT PASS\n")
cat("Plan SHA-256:", expected[["frozen_plan"]], "\n")
cat("Cells: 92,053 observed; 92,031 published\n")
cat("Marker-audited display cells: 90,512 in 14 lineages\n")
cat("Modules: GNMT+DMGDH versus SARDH+PIPOX\n")
