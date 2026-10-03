#!/usr/bin/env Rscript

# Export the four CRC WGS cohorts as three analysis-aligned pooled CSV files:
# metadata, species relative-abundance matrix, and KEGG KO relative-abundance
# matrix. The canonical roster is the 1,647 Run IDs in the frozen derived
# functional-score table, which is the complete metadata/species/KO intersection
# used for the pooled functional analysis.

suppressPackageStartupMessages(library(data.table))

options(stringsAsFactors = FALSE, scipen = 999)

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
if (length(script_file) != 1L) stop("Run with Rscript.")

integrated_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
crc_root <- normalizePath(file.path(integrated_root, ".."), mustWork = TRUE)
use_root <- file.path(crc_root, "사용데이터_모음")
cohorts <- c("PRJEB10878", "PRJEB27928", "PRJEB6070", "PRJNA429097")

score_file <- file.path(
  use_root,
  "Derived_Sarcosine_Functional_Scores",
  "CRC_WGS_sarcosine_functional_scores_all_cohorts.csv"
)

metadata_files <- setNames(
  file.path(use_root, cohorts, paste0(cohorts, "_patient_metadata_WGS_Healthy_Cancer.csv")),
  cohorts
)
species_files <- setNames(
  file.path(use_root, cohorts, paste0(cohorts, "_species_relative_abundance_matrix.csv")),
  cohorts
)
ko_files <- setNames(
  file.path(use_root, cohorts, paste0(cohorts, "_KEGG_KO_relative_abundance_matrix.csv")),
  cohorts
)

out_dir <- file.path(use_root, "CRC_WGS_4cohort_pooled_26.09.01")
provenance_dir <- file.path(
  integrated_root,
  "results_integrated",
  "CRC_WGS_4COHORT_POOLED_INPUT_EXPORT_26.09.01"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(provenance_dir, recursive = TRUE, showWarnings = FALSE)

metadata_out <- file.path(out_dir, "CRC_WGS_4cohort_pooled_metadata.csv")
species_out <- file.path(out_dir, "CRC_WGS_4cohort_pooled_species_relative_abundance_matrix.csv")
ko_out <- file.path(out_dir, "CRC_WGS_4cohort_pooled_KEGG_KO_relative_abundance_matrix.csv")

required <- c(score_file, metadata_files, species_files, ko_files)
if (!all(file.exists(required))) {
  stop("Missing input(s): ", paste(required[!file.exists(required)], collapse = "; "))
}

scores <- fread(score_file, check.names = FALSE)
required_score_cols <- c("Run_ID", "Cohort", "Analysis_Group")
if (!all(required_score_cols %in% names(scores))) stop("Canonical score table columns missing.")
scores <- scores[, ..required_score_cols]
scores[, Cohort := factor(Cohort, levels = cohorts)]
setorder(scores, Cohort, Run_ID)
scores[, Cohort := as.character(Cohort)]

stopifnot(
  nrow(scores) == 1647L,
  !anyDuplicated(scores$Run_ID),
  identical(sort(unique(scores$Cohort)), sort(cohorts)),
  all(scores$Analysis_Group %chin% c("Healthy", "Cancer"))
)

metadata_list <- vector("list", length(cohorts))
names(metadata_list) <- cohorts

for (cohort in cohorts) {
  # Read every metadata field as character so literal source entries such as
  # "NA" remain unchanged rather than being recoded as missing values.
  meta <- fread(
    metadata_files[[cohort]],
    check.names = FALSE,
    colClasses = "character",
    na.strings = NULL
  )
  required_meta <- c(
    "Run ID", "Analysis Group", "Cohort",
    "Included in species matrix", "Included in KO matrix"
  )
  if (!all(required_meta %in% names(meta))) stop("Metadata columns missing for ", cohort)
  if (anyDuplicated(meta[["Run ID"]])) stop("Duplicate metadata Run ID for ", cohort)

  ids <- scores[Cohort == cohort, Run_ID]
  if (!all(ids %chin% meta[["Run ID"]])) stop("Canonical Run missing from metadata for ", cohort)
  meta <- meta[match(ids, meta[["Run ID"]])]

  if (!identical(meta[["Run ID"]], ids)) stop("Metadata Run order failed for ", cohort)
  if (!identical(meta[["Cohort"]], rep(cohort, length(ids)))) stop("Metadata Cohort mismatch for ", cohort)
  if (!identical(meta[["Analysis Group"]], scores[Cohort == cohort, Analysis_Group])) {
    stop("Metadata Analysis Group mismatch for ", cohort)
  }
  if (!all(meta[["Included in species matrix"]] == "TRUE")) {
    stop("Canonical metadata includes a Run excluded from species matrix for ", cohort)
  }
  if (!all(meta[["Included in KO matrix"]] == "TRUE")) {
    stop("Canonical metadata includes a Run excluded from KO matrix for ", cohort)
  }
  metadata_list[[cohort]] <- meta
}

pooled_metadata <- rbindlist(metadata_list, use.names = TRUE, fill = FALSE)
if (!identical(pooled_metadata[["Run ID"]], scores$Run_ID)) {
  stop("Pooled metadata does not preserve the canonical Run order.")
}

pool_matrix <- function(paths, matrix_label) {
  pieces <- vector("list", length(cohorts))
  names(pieces) <- cohorts

  for (cohort in cohorts) {
    dat <- fread(paths[[cohort]], check.names = FALSE)
    if (!("Run_ID" %in% names(dat))) stop("Run_ID missing from ", matrix_label, " for ", cohort)
    if (anyDuplicated(dat$Run_ID)) stop("Duplicate Run_ID in ", matrix_label, " for ", cohort)
    feature_cols <- setdiff(names(dat), "Run_ID")
    if (length(feature_cols) == 0L || anyDuplicated(feature_cols)) {
      stop("Invalid feature columns in ", matrix_label, " for ", cohort)
    }
    if (!all(vapply(dat[, ..feature_cols], is.numeric, logical(1)))) {
      stop("Non-numeric feature column in ", matrix_label, " for ", cohort)
    }
    if (anyNA(dat[, ..feature_cols])) stop("Source NA in ", matrix_label, " for ", cohort)
    if (any(vapply(dat[, ..feature_cols], function(x) any(!is.finite(x) | x < 0), logical(1)))) {
      stop("Non-finite or negative value in ", matrix_label, " for ", cohort)
    }

    ids <- scores[Cohort == cohort, Run_ID]
    if (!all(ids %chin% dat$Run_ID)) stop("Canonical Run missing from ", matrix_label, " for ", cohort)
    dat <- dat[match(ids, Run_ID)]
    if (!identical(dat$Run_ID, ids)) stop("Run order failed in ", matrix_label, " for ", cohort)

    dat[, Cohort := cohort]
    dat[, Analysis_Group := scores[Cohort == cohort, Analysis_Group]]
    setcolorder(dat, c("Run_ID", "Cohort", "Analysis_Group", feature_cols))
    pieces[[cohort]] <- dat
  }

  pooled <- rbindlist(pieces, use.names = TRUE, fill = TRUE)
  feature_cols <- setdiff(names(pooled), c("Run_ID", "Cohort", "Analysis_Group"))
  setnafill(pooled, type = "const", fill = 0, cols = feature_cols)
  setcolorder(pooled, c("Run_ID", "Cohort", "Analysis_Group", feature_cols))

  if (!identical(pooled$Run_ID, scores$Run_ID)) stop("Pooled Run order failed for ", matrix_label)
  if (!identical(pooled$Cohort, scores$Cohort)) stop("Pooled Cohort mismatch for ", matrix_label)
  if (!identical(pooled$Analysis_Group, scores$Analysis_Group)) {
    stop("Pooled Analysis_Group mismatch for ", matrix_label)
  }
  if (anyDuplicated(pooled$Run_ID)) stop("Duplicated pooled Run IDs for ", matrix_label)
  if (anyNA(pooled)) stop("NA remained after zero-fill for ", matrix_label)

  pooled
}

pooled_species <- pool_matrix(species_files, "species matrix")
pooled_ko <- pool_matrix(ko_files, "KEGG KO matrix")

atomic_fwrite <- function(x, final_path) {
  partial <- paste0(final_path, ".partial")
  if (file.exists(partial)) unlink(partial)
  fwrite(x, partial, na = "NA", quote = "auto")
  if (!file.exists(partial) || file.info(partial)$size <= 0L) stop("CSV write failed: ", final_path)
  if (file.exists(final_path)) unlink(final_path)
  if (!file.rename(partial, final_path)) stop("Could not install CSV: ", final_path)
}

atomic_fwrite(pooled_metadata, metadata_out)
atomic_fwrite(pooled_species, species_out)
atomic_fwrite(pooled_ko, ko_out)

source_manifest <- data.table(
  input_type = c(
    "canonical_roster", rep("metadata", 4L), rep("species", 4L), rep("KEGG_KO", 4L)
  ),
  cohort = c("ALL", cohorts, cohorts, cohorts),
  file = required,
  md5 = unname(tools::md5sum(required))
)
fwrite(source_manifest, file.path(provenance_dir, "source_input_md5.csv"), quote = "auto")

output_manifest <- data.table(
  file = c(metadata_out, species_out, ko_out),
  rows = c(nrow(pooled_metadata), nrow(pooled_species), nrow(pooled_ko)),
  columns = c(ncol(pooled_metadata), ncol(pooled_species), ncol(pooled_ko)),
  md5 = unname(tools::md5sum(c(metadata_out, species_out, ko_out)))
)
fwrite(output_manifest, file.path(provenance_dir, "output_manifest.csv"), quote = "auto")

excluded_prjeb6070 <- setdiff(
  fread(metadata_files[["PRJEB6070"]], select = "Run ID", colClasses = "character")[["Run ID"]],
  scores[Cohort == "PRJEB6070", Run_ID]
)

validation <- c(
  "CRC WGS four-cohort pooled input export",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Canonical roster: 1,647 Run IDs present in metadata, species and KEGG KO inputs",
  paste0("Cohorts: ", paste(cohorts, collapse = ", ")),
  paste0("Pooled metadata: ", nrow(pooled_metadata), " rows x ", ncol(pooled_metadata), " columns"),
  paste0("Pooled species: ", nrow(pooled_species), " rows x ", ncol(pooled_species), " columns"),
  paste0("Pooled KEGG KO: ", nrow(pooled_ko), " rows x ", ncol(pooled_ko), " columns"),
  "Cohort-absent species/KO features were filled with numeric zero; observed source values were not transformed.",
  paste0("PRJEB6070 metadata Runs outside canonical roster: ", paste(excluded_prjeb6070, collapse = ", "))
)
writeLines(validation, file.path(provenance_dir, "validation_report.txt"))
capture.output(sessionInfo(), file = file.path(provenance_dir, "sessionInfo.txt"))

readme <- c(
  "CRC WGS four-cohort pooled CSV export",
  "",
  "Files in the use-data delivery folder",
  paste0("- ", basename(metadata_out), ": pooled metadata for the 1,647 canonical Runs"),
  paste0("- ", basename(species_out), ": union species relative-abundance matrix"),
  paste0("- ", basename(ko_out), ": union KEGG KO relative-abundance matrix"),
  "",
  "Rules",
  "- Run order is identical across all three CSV files.",
  "- Cohort and Analysis_Group are explicit in the two matrix files.",
  "- A feature absent from an entire cohort-specific source matrix is zero-filled for that cohort.",
  "- No prevalence filtering, normalization, transformation or statistical test was applied during export.",
  "- The canonical roster excludes four PRJEB6070 metadata Runs without complete metadata/species/KO overlap.",
  "",
  paste0("Delivery folder: ", out_dir)
)
writeLines(readme, file.path(provenance_dir, "README.txt"))

cat("Export complete\n")
print(output_manifest)

