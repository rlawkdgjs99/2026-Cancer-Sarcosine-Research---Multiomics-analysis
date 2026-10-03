#!/usr/bin/env Rscript

# Execute the complete numerical/figure workflow twice and require byte-identical
# machine-readable tables, verification tables, and PNG figures.

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
script_dir <- file.path(analysis_dir, "scripts")
table_dir <- file.path(analysis_dir, "results", "tables")
pub_dir <- file.path(analysis_dir, "results", "figures_publication")
log_dir <- file.path(analysis_dir, "logs")

suppressPackageStartupMessages(library(data.table))

rscript <- file.path(R.home("bin"), "Rscript")
main_script <- file.path(script_dir, "01_run_response_analysis.R")
verify_script <- file.path(script_dir, "02_verify_response_analysis.R")

sha256 <- function(path) sub(" .*", "", system2("shasum", c("-a", "256", path), stdout = TRUE))

run_workflow <- function(label) {
  main_log <- system2(rscript, main_script, stdout = TRUE, stderr = TRUE)
  main_status <- attr(main_log, "status")
  if (is.null(main_status)) main_status <- 0L
  if (main_status != 0L) stop(label, " main run failed: ", paste(main_log, collapse = "\n"))
  verify_log <- system2(rscript, verify_script, stdout = TRUE, stderr = TRUE)
  verify_status <- attr(verify_log, "status")
  if (is.null(verify_status)) verify_status <- 0L
  if (verify_status != 0L) stop(label, " verification failed: ", paste(verify_log, collapse = "\n"))
  writeLines(c(main_log, verify_log), file.path(log_dir, paste0("03_", label, "_workflow.log")))
}

run_workflow("run1")

table_targets <- list.files(table_dir, pattern = "^(01|02|03|04|05|06|07|08)_.*\\.csv$", full.names = TRUE)
png_targets <- list.files(pub_dir, pattern = "\\.png$", full.names = TRUE)
targets <- sort(c(table_targets, png_targets))
if (length(table_targets) < 9L || length(png_targets) != 6L) stop("Unexpected reproducibility target set")
run1 <- data.table(file = targets, sha256_run1 = vapply(targets, sha256, character(1)))

run_workflow("run2")
run2 <- data.table(file = targets, sha256_run2 = vapply(targets, sha256, character(1)))
comparison <- merge(run1, run2, by = "file", all = TRUE)
comparison[, file := sub(paste0("^", analysis_dir, "/"), "", file)]
comparison[, byte_identical := sha256_run1 == sha256_run2]
if (anyNA(comparison$byte_identical) || !all(comparison$byte_identical)) {
  stop("Double-run reproducibility failure")
}
fwrite(comparison, file.path(table_dir, "09_double_run_reproducibility.csv"))
writeLines(
  c(
    paste0("Double-run reproducibility PASS: ", nrow(comparison), "/", nrow(comparison), " targets byte-identical"),
    capture.output(sessionInfo())
  ),
  file.path(log_dir, "03_double_run_reproducibility_sessionInfo.txt")
)
message("Double-run reproducibility PASS: ", nrow(comparison), "/", nrow(comparison))

