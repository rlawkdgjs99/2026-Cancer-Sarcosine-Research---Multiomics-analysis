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

table_dir <- file.path(analysis_dir, "results", "tables")
intermediate_dir <- file.path(analysis_dir, "intermediate")
log_dir <- file.path(analysis_dir, "logs")

checks <- list()
add_check <- function(name, pass, detail) {
  checks[[length(checks) + 1L]] <<- data.table(check = name, pass = isTRUE(pass), detail = as.character(detail))
  if (!isTRUE(pass)) stop("Verification failed: ", name, " — ", detail)
}

sha256 <- function(path) {
  z <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  if (length(z) != 1L) stop("Could not hash ", path)
  strsplit(z, "[[:space:]]+")[[1]][1]
}

input_dir <- normalizePath(file.path(analysis_dir, "..", ".."), mustWork = TRUE)
umi <- file.path(input_dir, "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz")
metadata <- file.path(input_dir, "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
add_check("UMI SHA-256 unchanged",
          sha256(umi) == "aba15960fc7ee6a2443511bce5177e4d71b964131b6e98597e7a85d0a213ba36",
          sha256(umi))
add_check("metadata SHA-256 unchanged",
          sha256(metadata) == "d098a750c7ebc595994c929b666177f25c4da6fe3d3e38bb795269a1fb21053e",
          sha256(metadata))

counts <- readRDS(file.path(intermediate_dir, "01_full_counts_mt20_qc.rds"))
lineages <- fread(file.path(table_dir, "09_final_cell_lineages_FROZEN.csv"))
values <- fread(file.path(table_dir, "10_pseudobulk_target_counts_logCPM.csv"))
results <- fread(file.path(table_dir, "10_edger_target_results.csv"))
loo <- fread(file.path(table_dir, "10_primary_leave_one_patient_out.csv"))
stability <- fread(file.path(table_dir, "10_primary_stability_summary.csv"))
pb <- readRDS(file.path(intermediate_dir, "10_lineage_pseudobulk_counts.rds"))

add_check("count matrix dimensions", nrow(counts) == 24292L && ncol(counts) == 92053L,
          paste(dim(counts), collapse = " x "))
add_check("final lineage table dimensions", nrow(lineages) == 92053L && !anyDuplicated(lineages$cell_id),
          paste("rows", nrow(lineages)))
ord <- match(colnames(counts), lineages$cell_id)
add_check("count/lineage cell IDs align", !anyNA(ord), paste("unmatched", sum(is.na(ord))))
lineages <- lineages[ord]
add_check("count/lineage order exact", identical(colnames(counts), lineages$cell_id), "ordered")
add_check("eligible-cell count", sum(lineages$analysis_eligible) == 90512L,
          sum(lineages$analysis_eligible))
add_check("residual-doublet count", sum(lineages$final_lineage == "Excluded residual doublet") == 209L,
          sum(lineages$final_lineage == "Excluded residual doublet"))
add_check("unresolved T/NK count", sum(lineages$final_lineage == "NK/gamma-delta T unresolved") == 1332L,
          sum(lineages$final_lineage == "NK/gamma-delta T unresolved"))

eligible_lineages <- sort(unique(lineages[analysis_eligible == TRUE]$final_lineage))
add_check("pseudobulk lineage names", identical(sort(names(pb)), eligible_lineages),
          paste(names(pb), collapse = ", "))
for (lineage in eligible_lineages) {
  x <- pb[[lineage]]
  add_check(paste0("pseudobulk structure: ", lineage),
            nrow(x) == nrow(counts) && ncol(x) == 15L && identical(rownames(x), rownames(counts)),
            paste(dim(x), collapse = " x "))
}

# Independent raw target aggregation from the cell matrix, compared with every
# stored patient x lineage target count.
target_genes <- c("DMGDH", "ETFB", "GNMT", "PIPOX", "SARDH", "SHMT1", "SHMT2")
raw_checks <- list()
k <- 0L
for (lineage in eligible_lineages) {
  idx <- which(lineages$analysis_eligible & lineages$final_lineage == lineage)
  samples <- sort(unique(lineages$Sample))
  membership <- sparseMatrix(
    i = seq_along(idx), j = match(lineages$Sample[idx], samples), x = 1,
    dims = c(length(idx), length(samples)), dimnames = list(NULL, samples)
  )
  direct <- as.matrix(counts[target_genes, idx, drop = FALSE] %*% membership)
  direct <- rbind(direct, SARDH_PIPOX_SUM = direct["SARDH", ] + direct["PIPOX", ])
  expected <- values[final_lineage == lineage & target %in% rownames(direct),
                     .(Sample, target, raw_count)]
  observed <- rbindlist(lapply(rownames(direct), function(target) {
    data.table(Sample = colnames(direct), target = target, independently_aggregated = as.numeric(direct[target, ]))
  }))
  z <- merge(expected, observed, by = c("Sample", "target"), all = TRUE)
  k <- k + 1L
  raw_checks[[k]] <- z[, .(
    final_lineage = lineage, max_absolute_difference = max(abs(raw_count - independently_aggregated)),
    compared_rows = .N
  )]
}
raw_audit <- rbindlist(raw_checks)
fwrite(raw_audit, file.path(table_dir, "11_independent_raw_count_audit.csv"))
add_check("independent target-count aggregation", all(raw_audit$max_absolute_difference == 0),
          paste("max difference", max(raw_audit$max_absolute_difference)))

# Check multiplicity calculations exactly within each prespecified family and
# threshold.
for (threshold in c(10L, 50L)) {
  for (family in unique(results$family)) {
    idx <- which(results$cell_threshold == threshold & results$family == family & results$status == "TESTED")
    expected_q <- p.adjust(results$PValue[idx], method = "BH")
    delta <- if (length(idx)) max(abs(expected_q - results$FDR_BH[idx])) else 0
    add_check(paste("BH adjustment", threshold, family), delta < 1e-14,
              paste("max difference", format(delta, scientific = TRUE)))
  }
}

tested <- results[family == "primary_direct_degradation" & cell_threshold == 10L & status == "TESTED"]
add_check("primary threshold-10 tested row count", nrow(tested) == 11L, nrow(tested))
add_check("no primary FDR signal", all(tested$FDR_BH >= 0.05),
          paste("minimum q", signif(min(tested$FDR_BH), 6)))
add_check("tested confidence intervals contain estimates",
          all(tested$CI95_low <= tested$logFC_MPR_minus_NMPR &
              tested$CI95_high >= tested$logFC_MPR_minus_NMPR),
          "all tested rows")

# Recalculate LOO direction counts from fold-level rows.
loo_calc <- loo[status == "TESTED", .(
  LOO_folds_check = .N,
  positive_folds_check = sum(logFC_MPR_minus_NMPR > 0),
  negative_folds_check = sum(logFC_MPR_minus_NMPR < 0),
  min_logFC_check = min(logFC_MPR_minus_NMPR),
  max_logFC_check = max(logFC_MPR_minus_NMPR)
), by = .(final_lineage, target)]
stab <- merge(
  stability[status == "TESTED", .(final_lineage, target, LOO_folds, positive_folds,
                                   negative_folds, min_logFC, max_logFC)],
  loo_calc, by = c("final_lineage", "target"), all.x = TRUE
)
loo_delta <- max(c(
  abs(stab$LOO_folds - stab$LOO_folds_check),
  abs(stab$positive_folds - stab$positive_folds_check),
  abs(stab$negative_folds - stab$negative_folds_check),
  abs(stab$min_logFC - stab$min_logFC_check),
  abs(stab$max_logFC - stab$max_logFC_check)
), na.rm = TRUE)
add_check("LOO stability summaries", loo_delta < 1e-14,
          paste("max difference", format(loo_delta, scientific = TRUE)))

availability <- fread(file.path(table_dir, "10_Lee_9gene_signature_availability.csv"))
add_check("BHMT/BHMT2 absence explicit",
          identical(availability[present_in_deposited_matrix == FALSE]$target, c("BHMT", "BHMT2")),
          paste(availability[present_in_deposited_matrix == FALSE]$target, collapse = ", "))
missing_rows <- results[target %in% c("BHMT", "BHMT2")]
add_check("missing signature genes not modeled",
          nrow(missing_rows) == 56L && all(missing_rows$status == "NOT_AVAILABLE_IN_DEPOSITED_MATRIX"),
          paste("rows", nrow(missing_rows)))

check_table <- rbindlist(checks)
fwrite(check_table, file.path(table_dir, "11_verification_checks.csv"))
capture.output(sessionInfo(), file = file.path(log_dir, "11_sessionInfo.txt"))
cat("Independent pseudobulk verification PASS\n")
cat("Checks:", nrow(check_table), "\n")
cat("Minimum primary q:", signif(min(tested$FDR_BH), 6), "\n")
