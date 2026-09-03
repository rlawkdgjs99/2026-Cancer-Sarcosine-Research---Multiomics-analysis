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
counts <- readRDS(file.path(analysis_dir, "intermediate", "01_full_counts_mt20_qc.rds"))
lineages <- fread(file.path(table_dir, "09_final_cell_lineages_FROZEN.csv"))
reported <- fread(file.path(table_dir, "12_patient_lineage_target_celllevel_summary.csv"))

ord <- match(colnames(counts), lineages$cell_id)
if (anyNA(ord)) stop("Input alignment failed")
lineages <- lineages[ord]
eligible <- which(lineages$analysis_eligible == TRUE)

# Independent sparse-matrix aggregation, distinct from the data.table grouping
# implementation used by script 12.
key <- paste(lineages$Patient[eligible], lineages$Sample[eligible],
             lineages$final_lineage[eligible], sep = "\037")
key_levels <- sort(unique(key))
membership <- sparseMatrix(
  i = seq_along(eligible),
  j = match(key, key_levels),
  x = 1,
  dims = c(length(eligible), length(key_levels))
)
parts <- tstrsplit(key_levels, "\037", fixed = TRUE)
group_cells <- as.numeric(Matrix::colSums(membership))

independent_rows <- lapply(c("SARDH", "PIPOX"), function(target) {
  raw <- as.numeric(counts[target, eligible])
  normalized <- log1p(raw / as.numeric(lineages$library_size[eligible]) * 10000)
  data.table(
    Patient = parts[[1]],
    Sample = parts[[2]],
    final_lineage = parts[[3]],
    target = target,
    cells = as.integer(group_cells),
    detected_cells = as.integer(as.numeric(crossprod(membership, as.numeric(raw > 0)))),
    detection_pct = as.numeric(crossprod(membership, as.numeric(raw > 0))) / group_cells * 100,
    mean_log1p_CP10k = as.numeric(crossprod(membership, normalized)) / group_cells,
    eligible_ge10_cells = group_cells >= 10L
  )
})
independent <- rbindlist(independent_rows)

setkey(independent, Patient, Sample, final_lineage, target)
setkey(reported, Patient, Sample, final_lineage, target)
if (!identical(independent[, .(Patient, Sample, final_lineage, target)],
               reported[, .(Patient, Sample, final_lineage, target)])) {
  stop("Reported and independently aggregated keys differ")
}

checks <- data.table(
  check = c(
    "reported row count",
    "independent cell counts",
    "independent detected-cell counts",
    "independent detection percentages",
    "independent mean log1p(CP10k)",
    "independent >=10-cell flags"
  ),
  pass = c(
    nrow(reported) == nrow(independent),
    identical(reported$cells, independent$cells),
    identical(reported$detected_cells, independent$detected_cells),
    max(abs(reported$detection_pct - independent$detection_pct)) < 1e-12,
    max(abs(reported$mean_log1p_CP10k - independent$mean_log1p_CP10k)) < 1e-12,
    identical(reported$eligible_ge10_cells, independent$eligible_ge10_cells)
  ),
  detail = c(
    paste(nrow(reported), "rows"),
    paste("max difference", max(abs(reported$cells - independent$cells))),
    paste("max difference", max(abs(reported$detected_cells - independent$detected_cells))),
    format(max(abs(reported$detection_pct - independent$detection_pct)), scientific = TRUE),
    format(max(abs(reported$mean_log1p_CP10k - independent$mean_log1p_CP10k)), scientific = TRUE),
    paste(sum(reported$eligible_ge10_cells), "eligible rows")
  )
)

fwrite(checks, file.path(table_dir, "13_lineage_localization_verification_checks.csv"))
if (!all(checks$pass)) stop("At least one localization verification check failed")
message("PASS: all ", nrow(checks), " independent localization checks passed")
