#!/usr/bin/env Rscript

# Independent source-fidelity verification for script 117 pooled CSV exports.

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

metadata_out <- file.path(out_dir, "CRC_WGS_4cohort_pooled_metadata.csv")
species_out <- file.path(out_dir, "CRC_WGS_4cohort_pooled_species_relative_abundance_matrix.csv")
ko_out <- file.path(out_dir, "CRC_WGS_4cohort_pooled_KEGG_KO_relative_abundance_matrix.csv")
manifest_file <- file.path(provenance_dir, "output_manifest.csv")

checks <- character()
failures <- character()
record <- function(ok, label) {
  if (isTRUE(ok)) checks <<- c(checks, paste("PASS", label, sep = "\t"))
  else failures <<- c(failures, paste("FAIL", label, sep = "\t"))
}

required <- c(
  score_file, metadata_files, species_files, ko_files,
  metadata_out, species_out, ko_out, manifest_file
)
record(all(file.exists(required)), "all source and pooled export files exist")
if (!all(file.exists(required))) {
  writeLines(c(checks, failures, "TOTAL_FAIL\t1"), file.path(provenance_dir, "independent_verification_report.txt"))
  stop("Required file missing.")
}

scores <- fread(score_file, select = c("Run_ID", "Cohort", "Analysis_Group"))
scores[, Cohort := factor(Cohort, levels = cohorts)]
setorder(scores, Cohort, Run_ID)
scores[, Cohort := as.character(Cohort)]
record(nrow(scores) == 1647L, "canonical roster contains 1,647 Run IDs")
record(!anyDuplicated(scores$Run_ID), "canonical Run IDs are globally unique")

metadata <- fread(metadata_out, check.names = FALSE, colClasses = "character", na.strings = NULL)
species <- fread(species_out, check.names = FALSE)
ko <- fread(ko_out, check.names = FALSE)

record(nrow(metadata) == 1647L && ncol(metadata) == 27L, "pooled metadata dimensions are 1,647 x 27")
record(nrow(species) == 1647L, "pooled species matrix contains 1,647 rows")
record(nrow(ko) == 1647L, "pooled KEGG KO matrix contains 1,647 rows")
record(identical(metadata[["Run ID"]], scores$Run_ID), "metadata Run order matches canonical roster")
record(identical(species$Run_ID, scores$Run_ID), "species Run order matches canonical roster")
record(identical(ko$Run_ID, scores$Run_ID), "KEGG KO Run order matches canonical roster")
record(
  identical(species$Cohort, scores$Cohort) && identical(ko$Cohort, scores$Cohort),
  "matrix Cohort columns match canonical roster"
)
record(
  identical(species$Analysis_Group, scores$Analysis_Group) &&
    identical(ko$Analysis_Group, scores$Analysis_Group),
  "matrix Analysis_Group columns match canonical roster"
)

species_feature_cols <- setdiff(names(species), c("Run_ID", "Cohort", "Analysis_Group"))
ko_feature_cols <- setdiff(names(ko), c("Run_ID", "Cohort", "Analysis_Group"))

source_species_union <- character()
source_ko_union <- character()
for (cohort in cohorts) {
  source_species_union <- union(
    source_species_union,
    setdiff(names(fread(species_files[[cohort]], nrows = 0L, check.names = FALSE)), "Run_ID")
  )
  source_ko_union <- union(
    source_ko_union,
    setdiff(names(fread(ko_files[[cohort]], nrows = 0L, check.names = FALSE)), "Run_ID")
  )
}
record(identical(species_feature_cols, source_species_union), "species columns equal the exact first-seen source union")
record(identical(ko_feature_cols, source_ko_union), "KEGG KO columns equal the exact first-seen source union")
record(length(species_feature_cols) == 2628L, "pooled species matrix contains 2,628 union species")
record(length(ko_feature_cols) == 13098L, "pooled KEGG KO matrix contains 13,098 union KOs")

metadata_exact <- TRUE
species_values_exact <- TRUE
species_zero_fill_exact <- TRUE
ko_values_exact <- TRUE
ko_zero_fill_exact <- TRUE

for (cohort in cohorts) {
  ids <- scores[Cohort == cohort, Run_ID]
  pooled_idx <- match(ids, scores$Run_ID)

  source_meta <- fread(
    metadata_files[[cohort]],
    check.names = FALSE,
    colClasses = "character",
    na.strings = NULL
  )
  source_meta <- source_meta[match(ids, source_meta[["Run ID"]])]
  pooled_meta <- metadata[pooled_idx]
  metadata_exact <- metadata_exact && identical(source_meta, pooled_meta)

  source_species <- fread(species_files[[cohort]], check.names = FALSE)
  source_species <- source_species[match(ids, Run_ID)]
  source_species_cols <- setdiff(names(source_species), "Run_ID")
  species_values_exact <- species_values_exact && isTRUE(all.equal(
    as.matrix(species[pooled_idx, ..source_species_cols]),
    as.matrix(source_species[, ..source_species_cols]),
    tolerance = 1e-15,
    check.attributes = FALSE
  ))
  absent_species <- setdiff(species_feature_cols, source_species_cols)
  if (length(absent_species) > 0L) {
    species_zero_fill_exact <- species_zero_fill_exact && all(vapply(
      species[pooled_idx, ..absent_species], function(x) all(x == 0), logical(1)
    ))
  }

  source_ko <- fread(ko_files[[cohort]], check.names = FALSE)
  source_ko <- source_ko[match(ids, Run_ID)]
  source_ko_cols <- setdiff(names(source_ko), "Run_ID")
  ko_values_exact <- ko_values_exact && isTRUE(all.equal(
    as.matrix(ko[pooled_idx, ..source_ko_cols]),
    as.matrix(source_ko[, ..source_ko_cols]),
    tolerance = 1e-15,
    check.attributes = FALSE
  ))
  absent_ko <- setdiff(ko_feature_cols, source_ko_cols)
  if (length(absent_ko) > 0L) {
    ko_zero_fill_exact <- ko_zero_fill_exact && all(vapply(
      ko[pooled_idx, ..absent_ko], function(x) all(x == 0), logical(1)
    ))
  }
}

record(metadata_exact, "every pooled metadata cell matches its cohort source row")
record(species_values_exact, "every observed species value matches its cohort source")
record(species_zero_fill_exact, "cohort-absent species features are exactly zero-filled")
record(ko_values_exact, "every observed KEGG KO value matches its cohort source")
record(ko_zero_fill_exact, "cohort-absent KEGG KO features are exactly zero-filled")

record(!anyNA(species) && !anyNA(ko), "pooled matrix exports contain no NA cells")
record(
  all(vapply(species[, ..species_feature_cols], function(x) all(is.finite(x) & x >= 0), logical(1))),
  "all pooled species values are finite and non-negative"
)
record(
  all(vapply(ko[, ..ko_feature_cols], function(x) all(is.finite(x) & x >= 0), logical(1))),
  "all pooled KEGG KO values are finite and non-negative"
)

expected_counts <- data.table(
  Cohort = rep(cohorts, each = 2L),
  Analysis_Group = rep(c("Cancer", "Healthy"), times = 4L),
  N = c(74L, 54L, 140L, 120L, 590L, 476L, 98L, 95L)
)
observed_counts <- scores[, .N, by = .(Cohort, Analysis_Group)]
setorder(expected_counts, Cohort, Analysis_Group)
setorder(observed_counts, Cohort, Analysis_Group)
record(identical(expected_counts, observed_counts), "cohort-by-group counts match the frozen 1,647-Run roster")

prjeb6070_meta <- fread(metadata_files[["PRJEB6070"]], select = "Run ID", colClasses = "character")
excluded <- sort(setdiff(prjeb6070_meta[["Run ID"]], scores[Cohort == "PRJEB6070", Run_ID]))
record(
  identical(excluded, sort(c("ERR479028", "ERR479423", "ERR480520", "ERR480521"))),
  "the four excluded PRJEB6070 metadata Runs are explicitly accounted for"
)

manifest <- fread(manifest_file)
output_files <- c(metadata_out, species_out, ko_out)
manifest_match <- identical(
  manifest$md5[match(output_files, manifest$file)],
  unname(tools::md5sum(output_files))
)
record(manifest_match, "output manifest MD5 values match the installed CSV files")

report <- c(
  "Independent verification: CRC WGS four-cohort pooled CSV exports",
  paste0("Verified: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  checks,
  failures,
  paste0("TOTAL_PASS\t", length(checks)),
  paste0("TOTAL_FAIL\t", length(failures))
)
writeLines(report, file.path(provenance_dir, "independent_verification_report.txt"))
cat(paste(report, collapse = "\n"), "\n")

if (length(failures) > 0L) stop("Independent verification failed.")
