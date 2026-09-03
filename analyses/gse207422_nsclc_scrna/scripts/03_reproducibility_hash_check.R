#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
mode <- commandArgs(trailingOnly = TRUE)
if (length(mode) != 1L || !mode %in% c("snapshot", "compare")) {
  stop("Usage: 03_repro_hash_check.R snapshot|compare")
}

log_dir <- file.path(analysis_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
before_file <- file.path(log_dir, "03_repro_hashes_before.csv")
after_file <- file.path(log_dir, "03_repro_hashes_after.csv")
comparison_file <- file.path(log_dir, "03_repro_hash_comparison.csv")

targets <- c(
  list.files(file.path(analysis_dir, "results", "tables"), pattern = "\\.csv(\\.gz)?$", full.names = TRUE),
  list.files(file.path(analysis_dir, "results", "figures_publication"), pattern = "\\.png$", full.names = TRUE)
)
targets <- sort(normalizePath(targets, mustWork = TRUE))
sha256 <- function(path) sub(" .*", "", system2("shasum", c("-a", "256", path), stdout = TRUE))
manifest <- data.table(
  relative_path = sub(paste0("^", analysis_dir, "/"), "", targets),
  bytes = as.numeric(file.info(targets)$size),
  sha256 = vapply(targets, sha256, character(1))
)

if (mode == "snapshot") {
  fwrite(manifest, before_file)
  cat("REPRO SNAPSHOT WRITTEN:", nrow(manifest), "files\n")
} else {
  if (!file.exists(before_file)) stop("Before manifest is missing")
  before <- fread(before_file)
  fwrite(manifest, after_file)
  comparison <- merge(before, manifest, by = "relative_path", all = TRUE,
                      suffixes = c("_before", "_after"))
  comparison[, identical := !is.na(sha256_before) & !is.na(sha256_after) &
               bytes_before == bytes_after & sha256_before == sha256_after]
  fwrite(comparison, comparison_file)
  if (any(!comparison$identical)) {
    print(comparison[identical == FALSE])
    stop("Machine-readable or PNG reproducibility mismatch")
  }
  cat("REPRO HASH COMPARISON PASS:", nrow(comparison), "files byte-identical\n")
}
