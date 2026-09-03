#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
input_dir <- normalizePath(file.path(analysis_dir, "..", ".."), mustWork = TRUE)
.libPaths(c(file.path(analysis_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages(library(data.table))

reference_dir <- file.path(input_dir, "해당 연구팀 원본코드 from Github")
if (!dir.exists(reference_dir)) stop("Downloaded reference-code directory is missing")

files <- sort(list.files(reference_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE))
files <- files[file.info(files)$isdir %in% FALSE]
if (!length(files)) stop("No downloaded reference-code files found")

sha256_file <- function(path) {
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE, stderr = TRUE)
  if (!length(out) || !grepl("^[0-9a-f]{64}  ", out[1])) stop("SHA-256 failed for ", path)
  sub("  .*", "", out[1])
}

manifest <- data.table(
  relative_path = substring(files, nchar(reference_dir) + 2L),
  bytes = as.numeric(file.info(files)$size),
  sha256 = vapply(files, sha256_file, character(1))
)

expected <- c(
  "01.2_QC_each_matrix.R" = "761a6d2cdf75e57ef9c6639cc6590afd0620e71545d8f28eb895400e730f20ea",
  "01.3_Scrublet.py" = "047863d4e2be208c472843c29d3bd54213541dcf1ececb972b97ed242241ba9e",
  "02_All_cell_clustering.R" = "7809ab6d55a57ba839317dc2d26d9c28f864179536561053615481ff8cef17e2",
  "05_T_analysis_new.R" = "de584ca78170528ec0aa4b96983222791f6613fab42cc1fff55df0d35c5367d5"
)
observed <- setNames(manifest$sha256, manifest$relative_path)
if (!all(names(expected) %in% names(observed))) stop("One or more key reference scripts are absent")
if (!identical(unname(observed[names(expected)]), unname(expected))) stop("Key reference-script checksum mismatch")

table_dir <- file.path(analysis_dir, "results", "tables")
log_dir <- file.path(analysis_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
fwrite(manifest, file.path(table_dir, "00b_downloaded_reference_code_manifest.csv"))

audit <- data.table(
  item = c(
    "reference_snapshot_commit", "reference_file_count", "key_script_checksums",
    "source_per_sample_qc", "source_doublet_call", "source_all_cell_extra_qc",
    "source_all_cell_integration", "source_all_cell_resolution", "source_lineage_annotation"
  ),
  finding = c(
    "e4c837d72b9726f5dde6a7c1b42f6cdc22b981dd",
    as.character(nrow(manifest)),
    "PASS: 01.2, 01.3 Scrublet, 02, and 05 match the pinned snapshot",
    "nFeature_RNA >= 500; percent.mito <= 0.2; percent.ribo <= 0.5; ACTB+GAPDH+MALAT1 UMI >= 1, applied per sample",
    "Scrublet expected_doublet_rate=0.025; manually called threshold=0.22, applied per sample",
    "NONE after loading count_all_scrublet.txt and before all-cell clustering",
    "per-sample NormalizeData and 3000 vst features; CCA anchors/integration dims 1:20",
    "Louvain resolution 0.6 in the source-study all-cell script",
    "manual cluster annotation from marker-expression plots; exact marker panel and cluster map are present in 02_All_cell_clustering.R"
  )
)
fwrite(audit, file.path(table_dir, "00b_downloaded_reference_code_audit.csv"))
capture.output(sessionInfo(), file = file.path(log_dir, "00b_sessionInfo.txt"))

cat("Downloaded reference-code audit PASS\n")
cat("Files:", nrow(manifest), "\n")
cat("Key script checksums: PASS\n")
