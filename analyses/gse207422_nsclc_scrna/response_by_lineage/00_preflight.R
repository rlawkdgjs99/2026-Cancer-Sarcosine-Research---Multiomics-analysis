#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
parent_dir <- normalizePath(file.path(analysis_dir, ".."), mustWork = TRUE)
analysis_root <- normalizePath(file.path(parent_dir, ".."), mustWork = TRUE)
prior_dir <- file.path(analysis_root, "02_lineage_reannotation")
.libPaths(c(file.path(parent_dir, "R_libs"), file.path(prior_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages(library(data.table))

table_dir <- file.path(analysis_dir, "results", "tables")
log_dir <- file.path(analysis_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

score_file <- file.path(parent_dir, "results", "tables", "01_cell_paper_style_scores.csv.gz")
lineage_file <- file.path(prior_dir, "results", "tables", "09_final_cell_lineages_FROZEN.csv")
if (!all(file.exists(c(score_file, lineage_file)))) stop("Required verified input is missing")

sha256 <- function(path) sub(" .*", "", system2("shasum", c("-a", "256", path), stdout = TRUE))
expected_sha <- c(
  scores = "17724b2b1fcbad436e6ce4cee91a34bda851361aac7b14c26a3557e5d006fb61",
  lineages = "97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1"
)
observed_sha <- c(scores = sha256(score_file), lineages = sha256(lineage_file))
if (!identical(observed_sha, expected_sha)) stop("Verified input SHA-256 mismatch")

lineage_order <- c(
  "Epithelial", "CAF", "B cell", "Plasma cell", "CD4 T cell",
  "CD8 T cell", "Cycling T cell", "NK cell", "Mast cell",
  "Neutrophil", "Monocyte", "Macrophage", "Conventional DC", "pDC"
)

# Structure-only read: no score value is loaded before the analysis plan is frozen.
score_md <- fread(
  score_file,
  select = c(
    "cell_id", "Patient", "Sample", "Resource", "Pathologic.Response",
    "final_lineage", "analysis_eligible"
  )
)
lin_md <- fread(
  lineage_file,
  select = c(
    "cell_id", "Patient", "Sample", "Resource", "Pathologic.Response",
    "final_lineage", "analysis_eligible"
  )
)

if (nrow(score_md) != 92053L || nrow(lin_md) != 92053L) stop("Unexpected input rows")
if (anyDuplicated(score_md$cell_id) || anyDuplicated(lin_md$cell_id)) stop("Duplicate cell IDs")
if (!setequal(score_md$cell_id, lin_md$cell_id)) stop("Score/lineage cell-ID sets differ")
setkey(score_md, cell_id)
setkey(lin_md, cell_id)
metadata_cols <- setdiff(names(score_md), "cell_id")
if (!identical(score_md[, ..metadata_cols], lin_md[, ..metadata_cols])) {
  stop("Score-table metadata differ from frozen-lineage metadata")
}

eligible <- lin_md[
  analysis_eligible == TRUE &
    Resource == "Post-treatment surgery" &
    Pathologic.Response %chin% c("MPR", "pCR", "NMPR")
]
eligible[, response_group := fifelse(Pathologic.Response %chin% c("MPR", "pCR"), "MPR/pCR", "NMPR")]
eligible[, response_group := factor(response_group, levels = c("MPR/pCR", "NMPR"))]
eligible[, final_lineage := factor(final_lineage, levels = lineage_order)]

patient_meta <- unique(eligible[, .(Patient, Sample, Resource, Pathologic.Response, response_group)])
if (nrow(patient_meta) != 12L || uniqueN(patient_meta$Patient) != 12L) {
  stop("Expected 12 post-treatment evaluable patients")
}
group_n <- patient_meta[, .N, by = response_group][order(response_group)]
if (!identical(group_n$N, c(4L, 8L))) stop("Expected MPR/pCR=4 and NMPR=8 patients")
if (any(patient_meta$Resource != "Post-treatment surgery")) stop("Primary cohort includes non-surgical sample")
if (!setequal(unique(eligible$final_lineage), lineage_order)) stop("Primary cohort lacks an eligible lineage entirely")

lineage_counts <- eligible[, .(n_cells = .N), by = .(Patient, Sample, Pathologic.Response, response_group, final_lineage)]
full_grid <- CJ(
  Patient = patient_meta$Patient,
  final_lineage = factor(lineage_order, levels = lineage_order),
  unique = TRUE
)
lineage_counts <- merge(
  full_grid,
  patient_meta[, .(Patient, Sample, Pathologic.Response, response_group)],
  by = "Patient",
  all.x = TRUE
)
lineage_counts <- merge(
  lineage_counts,
  eligible[, .(n_cells = .N), by = .(Patient, final_lineage)],
  by = c("Patient", "final_lineage"),
  all.x = TRUE
)
lineage_counts[is.na(n_cells), n_cells := 0L]
setorder(lineage_counts, final_lineage, response_group, Patient)

thresholds <- c(1L, 10L, 50L, 100L)
threshold_summary <- rbindlist(lapply(thresholds, function(min_cells) {
  lineage_counts[, .(
    patients_eligible = sum(n_cells >= min_cells),
    cells_in_eligible_patients = sum(n_cells[n_cells >= min_cells]),
    min_cells_among_eligible = if (any(n_cells >= min_cells)) as.numeric(min(n_cells[n_cells >= min_cells])) else NA_real_,
    median_cells_among_eligible = if (any(n_cells >= min_cells)) as.numeric(median(n_cells[n_cells >= min_cells])) else NA_real_,
    max_cells_among_eligible = if (any(n_cells >= min_cells)) as.numeric(max(n_cells[n_cells >= min_cells])) else NA_real_
  ), by = .(final_lineage, response_group)][, min_cells := min_cells]
}), use.names = TRUE)
setcolorder(threshold_summary, c(
  "min_cells", "final_lineage", "response_group", "patients_eligible",
  "cells_in_eligible_patients", "min_cells_among_eligible",
  "median_cells_among_eligible", "max_cells_among_eligible"
))
setorder(threshold_summary, min_cells, final_lineage, response_group)

input_manifest <- data.table(
  role = c("verified_cell_score_table", "frozen_lineage_metadata"),
  path = c(normalizePath(score_file), normalizePath(lineage_file)),
  bytes = file.info(c(score_file, lineage_file))$size,
  sha256 = observed_sha,
  expected_sha256 = expected_sha,
  checksum_status = ifelse(observed_sha == expected_sha, "MATCH", "MISMATCH")
)

fwrite(input_manifest, file.path(table_dir, "00_input_manifest.csv"))
fwrite(patient_meta[order(response_group, Patient)], file.path(table_dir, "00_primary_patient_metadata.csv"))
fwrite(lineage_counts, file.path(table_dir, "00_patient_lineage_cell_counts.csv"))
fwrite(threshold_summary, file.path(table_dir, "00_threshold_eligibility_summary.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "00_preflight_sessionInfo.txt"))

message("Structure-only preflight completed without loading score values")
