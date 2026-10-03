#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))
args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
repro_dir <- "/private/tmp/lineage_quartile_GOBP_repro_260826"
if (!dir.exists(repro_dir)) stop("Reproducibility rerun directory is missing")

sha256 <- function(path) {
  sub(" .*", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
}

relative_files <- c(
  file.path("results", "tables", c(
    "00_input_manifest.csv",
    "00_package_versions.csv",
    "01_patient_contribution_by_quartile.csv",
    "01_quartile_thresholds_and_counts.csv",
    "01_selected_cell_quartile_membership.csv.gz",
    "02_DEG_enrichment_input_audit.csv",
    "02_lineage_quartile_FindMarkers_all.csv",
    "03_lineage_quartile_GO_BP_all_terms.csv",
    "03_lineage_quartile_GO_BP_directional_all.csv",
    "plotdata_top_GO_BP_terms.csv"
  )),
  file.path("results", "figures_publication", c(
    "Fig_CAF_ratio_Q4_vs_Q1_GO_BP.png",
    "Fig_Epithelial_ratio_Q4_vs_Q1_GO_BP.png",
    "Fig_Epithelial_CAF_ratio_Q4_vs_Q1_GO_BP_combined.png"
  )),
  file.path("results", "figures_diagnostic", "Diagnostic_ratio_quartile_definition.png")
)

main_files <- file.path(analysis_dir, relative_files)
repro_files <- file.path(repro_dir, relative_files)
if (!all(file.exists(main_files)) || !all(file.exists(repro_files))) {
  stop("One or more comparison files are missing")
}

comparison <- data.table(
  relative_path = relative_files,
  main_bytes = file.info(main_files)$size,
  repro_bytes = file.info(repro_files)$size,
  main_sha256 = vapply(main_files, sha256, character(1)),
  repro_sha256 = vapply(repro_files, sha256, character(1))
)
comparison[, identical := main_bytes == repro_bytes & main_sha256 == repro_sha256]
fwrite(comparison, file.path(analysis_dir, "logs", "03_reproducibility_comparison.csv"))

if (!all(comparison$identical)) {
  stop("Reproducibility check failed for: ",
       paste(comparison[identical == FALSE, relative_path], collapse = "; "))
}
message("Clean rerun reproducibility complete: ", nrow(comparison), "/", nrow(comparison),
        " CSV/PNG files byte-identical")
