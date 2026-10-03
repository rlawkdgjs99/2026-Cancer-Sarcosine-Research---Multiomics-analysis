#!/usr/bin/env Rscript

# Validate the CRC target-KO cache used by the cross-disease concordance analysis
# against the original HGMT KO_relative_abundance.tsv JSON inputs.
# This script does not alter either input.

set.seed(42)
options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})
if (!requireNamespace("digest", quietly = TRUE)) stop("Package digest is required.")

script_file <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
script_dir <- dirname(script_file)
nsc_root <- normalizePath(file.path(script_dir, "..", ".."), mustWork = TRUE)
collection_root <- dirname(nsc_root)
crc_root <- file.path(collection_root, "HGMT_CRC_WGS-Healthy_vs_Cancer",
                      "Healthy_vs_Cancer_4_CRC_cohorts_integrated")
out_dir <- file.path(nsc_root, "pooled_analysis", "results",
                     "Fig2_cross_disease_species_degradation_concordance_CANDIDATE_26.08.23")
cache_file <- file.path(crc_root, "results_integrated", "sarcosine",
                        "sarcosine_KO_per_sample_pooled.csv")

cohort_dirs <- c(
  PRJEB10878 = "PRJEB10878_CRC",
  PRJEB27928 = "PRJEB27928_CRC",
  PRJEB6070 = "PRJEB6070_CRC_AdenomatousPolyps",
  PRJNA429097 = "PRJNA429097_CRC"
)
target_kos <- c("K00301", "K00302", "K00303", "K00305", "K00306")

cache <- fread(cache_file)
stopifnot(!anyDuplicated(cache$Run.ID), all(target_kos %in% names(cache)))
checks <- list()
raw_files <- character()

for (cohort in names(cohort_dirs)) {
  cat("Reading raw KO JSON:", cohort, "\n")
  raw_file <- file.path(crc_root, cohort_dirs[[cohort]], "KO_relative_abundance.tsv")
  raw_files <- c(raw_files, raw_file)
  raw <- as.data.table(jsonlite::fromJSON(raw_file))
  stopifnot(all(c("run_id", "ko", "abundance") %in% names(raw)))
  if (!is.numeric(raw$abundance) || anyNA(raw$abundance) ||
      any(!is.finite(raw$abundance)) || any(raw$abundance < 0)) {
    stop("Invalid raw KO abundance: ", cohort)
  }
  ccache <- cache[Cohort == cohort]
  target <- raw[run_id %in% ccache$Run.ID & ko %in% target_kos,
                .(abundance = sum(abundance)), by = .(run_id, ko)]
  if (anyDuplicated(paste(target$run_id, target$ko))) {
    stop("Duplicate raw run/KO after collapse: ", cohort)
  }
  grid <- CJ(run_id = ccache$Run.ID, ko = target_kos, unique = TRUE)
  target <- target[grid, on = .(run_id, ko)]
  target[is.na(abundance), abundance := 0]
  wide <- dcast(target, run_id ~ ko, value.var = "abundance")
  setorder(wide, run_id)
  setorder(ccache, Run.ID)
  stopifnot(identical(wide$run_id, ccache$Run.ID))
  raw_mat <- as.matrix(wide[, ..target_kos])
  cache_mat <- as.matrix(ccache[, ..target_kos])
  delta <- raw_mat - cache_mat
  checks[[cohort]] <- data.table(
    cohort = cohort, cache_runs = nrow(ccache), target_KOs = length(target_kos),
    compared_cells = length(delta), max_abs_difference = max(abs(delta)),
    cells_abs_diff_gt_1e_15 = sum(abs(delta) > 1e-15),
    cells_abs_diff_gt_1e_12 = sum(abs(delta) > 1e-12)
  )
  rm(raw, target, grid, wide, raw_mat, cache_mat, delta)
  invisible(gc(verbose = FALSE))
}

checks <- rbindlist(checks)
if (any(checks$cells_abs_diff_gt_1e_12 > 0L)) {
  stop("CRC KO cache differs materially from raw HGMT KO JSON input.")
}
fwrite(checks, file.path(out_dir, "CRC_KO_cache_vs_raw_validation.csv"))

hashes <- data.table(
  file = c(cache_file, raw_files, script_file),
  bytes = file.info(c(cache_file, raw_files, script_file))$size,
  sha256 = vapply(c(cache_file, raw_files, script_file), digest::digest,
                  character(1), algo = "sha256", file = TRUE)
)
fwrite(hashes, file.path(out_dir, "CRC_KO_cache_validation_input_sha256.tsv"), sep = "\t")
writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_CRC_KO_cache_validation.txt"))

stable_outputs <- list.files(out_dir, full.names = TRUE)
stable_outputs <- stable_outputs[
  !basename(stable_outputs) %in% c("analysis_log.txt", "output_sha256.tsv")
]
output_hashes <- data.table(
  file = basename(stable_outputs), bytes = file.info(stable_outputs)$size,
  sha256 = vapply(stable_outputs, digest::digest, character(1),
                  algo = "sha256", file = TRUE)
)
fwrite(output_hashes, file.path(out_dir, "output_sha256.tsv"), sep = "\t")
cat("Validation complete. Maximum absolute difference:",
    format(max(checks$max_abs_difference), scientific = TRUE), "\n")
