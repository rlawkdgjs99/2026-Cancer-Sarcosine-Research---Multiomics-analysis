#!/usr/bin/env Rscript

# Independent verification of the frozen patient-level response analysis.
# This script recomputes every numerical result directly from the locked inputs
# and only then compares it with the files written by 01_run_response_analysis.R.

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
parent_dir <- normalizePath(file.path(analysis_dir, ".."), mustWork = TRUE)
analysis_root <- normalizePath(file.path(parent_dir, ".."), mustWork = TRUE)
prior_dir <- file.path(analysis_root, "02_lineage_reannotation")
.libPaths(c(file.path(parent_dir, "R_libs"), file.path(prior_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(png)
})

options(stringsAsFactors = FALSE, warn = 1)
MASTER_SEED <- 260825L

table_dir <- file.path(analysis_dir, "results", "tables")
pub_dir <- file.path(analysis_dir, "results", "figures_publication")
log_dir <- file.path(analysis_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) sub(" .*", "", system2("shasum", c("-a", "256", path), stdout = TRUE))
plan_file <- file.path(analysis_dir, "ANALYSIS_PLAN_FROZEN.md")
score_file <- file.path(parent_dir, "results", "tables", "01_cell_paper_style_scores.csv.gz")
lineage_file <- file.path(prior_dir, "results", "tables", "09_final_cell_lineages_FROZEN.csv")

expected_hashes <- c(
  plan = "a20203c72f28c344530f1a729885e37c51c32cdb00dd694157ac36c08b976bfe",
  scores = "17724b2b1fcbad436e6ce4cee91a34bda851361aac7b14c26a3557e5d006fb61",
  lineages = "97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1"
)
observed_hashes <- c(plan = sha256(plan_file), scores = sha256(score_file), lineages = sha256(lineage_file))
if (!identical(observed_hashes, expected_hashes)) stop("Locked input or plan hash mismatch")

lineage_order <- c(
  "Epithelial", "CAF", "B cell", "Plasma cell", "CD4 T cell",
  "CD8 T cell", "Cycling T cell", "NK cell", "Mast cell",
  "Neutrophil", "Monocyte", "Macrophage", "Conventional DC", "pDC"
)
score_order <- c("Production", "Degradation", "Production/degradation ratio")

q25 <- function(x) unname(quantile(x, 0.25, type = 7, na.rm = TRUE))
q75 <- function(x) unname(quantile(x, 0.75, type = 7, na.rm = TRUE))

# Alternative implementations to those in the production script.
independent_exact_p <- function(a, b) {
  pooled <- c(a, b)
  n_a <- length(a)
  rank_values <- rank(pooled, ties.method = "average")
  observed_deviation <- abs(sum(rank_values[seq_len(n_a)]) - n_a * (length(pooled) + 1) / 2)
  allocations <- combn(seq_along(pooled), n_a, simplify = FALSE)
  deviations <- vapply(
    allocations,
    function(index) abs(sum(rank_values[index]) - n_a * (length(pooled) + 1) / 2),
    numeric(1)
  )
  sum(deviations + 1e-12 >= observed_deviation) / length(deviations)
}

independent_hl <- function(a, b) median(rep(a, each = length(b)) - rep(b, times = length(a)))

independent_boot_ci <- function(a, b, seed, B = 10000L) {
  set.seed(seed)
  values <- numeric(B)
  for (i in seq_len(B)) {
    aa <- a[sample.int(length(a), length(a), replace = TRUE)]
    bb <- b[sample.int(length(b), length(b), replace = TRUE)]
    values[i] <- independent_hl(aa, bb)
  }
  unname(quantile(values, c(0.025, 0.975), type = 7))
}

score <- fread(score_file)
lineage <- fread(
  lineage_file,
  select = c("cell_id", "Patient", "Sample", "Resource", "Pathologic.Response", "final_lineage", "analysis_eligible")
)
if (nrow(score) != 92053L || nrow(lineage) != 92053L) stop("Unexpected input row count")
if (anyDuplicated(score$cell_id) || anyDuplicated(lineage$cell_id)) stop("Duplicate cell IDs")
setorder(score, cell_id)
setorder(lineage, cell_id)
if (!identical(score$cell_id, lineage$cell_id)) stop("Cell-ID mismatch")
for (nm in c("Patient", "Sample", "Resource", "Pathologic.Response", "final_lineage", "analysis_eligible")) {
  if (!identical(score[[nm]], lineage[[nm]])) stop("Metadata mismatch in ", nm)
}

cell <- score[
  analysis_eligible == TRUE &
    Resource == "Post-treatment surgery" &
    Pathologic.Response %chin% c("MPR", "pCR", "NMPR")
]
cell[, response_group := ifelse(Pathologic.Response %chin% c("MPR", "pCR"), "MPR/pCR", "NMPR")]
expected_patients <- c("P03", "P06", "P11", "P14", "P02", "P04", "P07", "P09", "P10", "P12", "P13", "P15")
if (!setequal(unique(cell$Patient), expected_patients)) stop("Primary patient set mismatch")

# Reconstruct one row per patient x lineage without reading any generated result.
independent_summary <- cell[, {
  ratio <- production_degradation_ratio[ratio_defined & is.finite(production_degradation_ratio)]
  prod <- production_module_shifted
  deg <- degradation_module_shifted
  list(
    n_cells = .N,
    n_finite_ratio_cells = length(ratio),
    production_mean = sum(prod) / length(prod),
    production_median = stats::median(prod),
    production_q1 = q25(prod),
    production_q3 = q75(prod),
    degradation_mean = sum(deg) / length(deg),
    degradation_median = stats::median(deg),
    degradation_q1 = q25(deg),
    degradation_q3 = q75(deg),
    ratio_median = if (length(ratio)) stats::median(ratio) else NA_real_,
    ratio_mean = if (length(ratio)) sum(ratio) / length(ratio) else NA_real_,
    ratio_q1 = if (length(ratio)) q25(ratio) else NA_real_,
    ratio_q3 = if (length(ratio)) q75(ratio) else NA_real_,
    ratio_of_means = (sum(prod) / length(prod)) / (sum(deg) / length(deg))
  )
}, by = .(Patient, Sample, Pathologic.Response, response_group, final_lineage)]
setorder(independent_summary, final_lineage, response_group, Patient)

make_long <- function(summary, threshold, aggregation = "primary") {
  d <- summary[n_cells >= threshold]
  if (aggregation == "primary") {
    values <- list(d$production_mean, d$degradation_mean, d$ratio_median)
  } else if (aggregation == "cell_median") {
    values <- list(d$production_median, d$degradation_median)
  } else if (aggregation == "ratio_of_means") {
    values <- list(d$ratio_of_means)
  } else stop("Unknown aggregation")
  labels <- switch(
    aggregation,
    primary = score_order,
    cell_median = score_order[1:2],
    ratio_of_means = score_order[3]
  )
  rbindlist(lapply(seq_along(values), function(i) {
    d[, .(Patient, response_group, final_lineage, n_cells, score = labels[i], value = values[[i]])]
  }))
}

analyze <- function(long, bootstrap = FALSE, seed_offset = 0L) {
  out <- long[, {
    a <- value[response_group == "MPR/pCR" & is.finite(value)]
    b <- value[response_group == "NMPR" & is.finite(value)]
    is_tested <- length(a) >= 3L && length(b) >= 3L
    ci <- c(NA_real_, NA_real_)
    if (is_tested && bootstrap) {
      ci <- independent_boot_ci(
        a, b,
        MASTER_SEED + seed_offset + 100L * match(score[1], score_order) + match(final_lineage[1], lineage_order)
      )
    }
    list(
      n_MPR_pCR = length(a),
      n_NMPR = length(b),
      MPR_pCR_median = if (length(a)) median(a) else NA_real_,
      MPR_pCR_q1 = if (length(a)) q25(a) else NA_real_,
      MPR_pCR_q3 = if (length(a)) q75(a) else NA_real_,
      NMPR_median = if (length(b)) median(b) else NA_real_,
      NMPR_q1 = if (length(b)) q25(b) else NA_real_,
      NMPR_q3 = if (length(b)) q75(b) else NA_real_,
      HL_shift_MPR_minus_NMPR = if (is_tested) independent_hl(a, b) else NA_real_,
      HL_bootstrap_CI_low = ci[1],
      HL_bootstrap_CI_high = ci[2],
      exact_P = if (is_tested) independent_exact_p(a, b) else NA_real_,
      tested = is_tested
    )
  }, by = .(final_lineage, score)]
  out[tested == TRUE, BH_q := p.adjust(exact_P, method = "BH")]
  out[tested == FALSE, BH_q := NA_real_]
  out[, family_tests := sum(tested)]
  setorder(out, score, final_lineage)
  out
}

independent_primary_long <- make_long(independent_summary, 10L, "primary")
independent_primary <- analyze(independent_primary_long, bootstrap = TRUE)
if (sum(independent_primary$tested) != 39L) stop("Independent computation did not yield 39 tests")

independent_threshold <- rbindlist(lapply(c(1L, 10L, 50L, 100L), function(threshold) {
  z <- analyze(make_long(independent_summary, threshold, "primary"), bootstrap = FALSE)
  z[, min_cells := threshold]
  z
}))
setorder(independent_threshold, min_cells, score, final_lineage)

median_result <- analyze(make_long(independent_summary, 10L, "cell_median"), bootstrap = FALSE)
median_result[, aggregation := "Cell median"]
ratio_result <- analyze(make_long(independent_summary, 10L, "ratio_of_means"), bootstrap = FALSE)
ratio_result[, aggregation := "Ratio of cell-score means"]
independent_aggregation <- rbindlist(list(median_result, ratio_result), fill = TRUE)
independent_aggregation[, family_tests := sum(tested), by = aggregation]
setorder(independent_aggregation, aggregation, score, final_lineage)

independent_loo <- independent_primary_long[, {
  rbindlist(lapply(unique(Patient), function(drop_patient) {
    d <- .SD[Patient != drop_patient]
    a <- d[response_group == "MPR/pCR", value]
    b <- d[response_group == "NMPR", value]
    data.table(
      dropped_patient = drop_patient,
      dropped_group = .SD[Patient == drop_patient, response_group][1],
      n_MPR_pCR = length(a),
      n_NMPR = length(b),
      HL_shift_MPR_minus_NMPR = if (length(a) >= 2L && length(b) >= 2L) independent_hl(a, b) else NA_real_
    )
  }))
}, by = .(final_lineage, score)]
setorder(independent_loo, score, final_lineage, dropped_patient)

assert_same <- function(independent, saved, keys, numeric_columns, label, tolerance = 5e-13) {
  a <- copy(independent)
  b <- copy(saved)
  for (key in keys) {
    a[, (key) := as.character(get(key))]
    b[, (key) := as.character(get(key))]
  }
  setorderv(a, keys)
  setorderv(b, keys)
  if (nrow(a) != nrow(b) || !identical(a[, ..keys], b[, ..keys])) stop(label, ": key mismatch")
  errors <- vapply(numeric_columns, function(column) {
    x <- a[[column]]
    y <- b[[column]]
    if (!identical(is.na(x), is.na(y))) return(Inf)
    finite <- is.finite(x) & is.finite(y)
    if (!any(finite)) return(0)
    max(abs(x[finite] - y[finite]))
  }, numeric(1))
  if (any(!is.finite(errors)) || any(errors > tolerance)) {
    stop(label, ": numeric mismatch; ", paste(names(errors), format(errors, scientific = TRUE), collapse = ", "))
  }
  max(errors)
}

saved_summary <- fread(file.path(table_dir, "01_patient_lineage_score_summaries.csv"))
saved_primary <- fread(file.path(table_dir, "02_primary_exact_results.csv"))
saved_threshold <- fread(file.path(table_dir, "03_cell_threshold_sensitivity.csv"))
saved_aggregation <- fread(file.path(table_dir, "04_aggregation_sensitivity.csv"))
saved_loo <- fread(file.path(table_dir, "05_leave_one_patient_out_HL.csv"))

summary_error <- assert_same(
  independent_summary, saved_summary,
  c("Patient", "Sample", "Pathologic.Response", "response_group", "final_lineage"),
  setdiff(names(independent_summary), c("Patient", "Sample", "Pathologic.Response", "response_group", "final_lineage")),
  "patient summaries"
)
result_columns <- c(
  "n_MPR_pCR", "n_NMPR", "MPR_pCR_median", "MPR_pCR_q1", "MPR_pCR_q3",
  "NMPR_median", "NMPR_q1", "NMPR_q3", "HL_shift_MPR_minus_NMPR",
  "HL_bootstrap_CI_low", "HL_bootstrap_CI_high", "exact_P", "BH_q", "family_tests"
)
primary_error <- assert_same(independent_primary, saved_primary, c("final_lineage", "score"), result_columns, "primary results")
threshold_error <- assert_same(
  independent_threshold, saved_threshold,
  c("min_cells", "final_lineage", "score"),
  setdiff(result_columns, c("HL_bootstrap_CI_low", "HL_bootstrap_CI_high")),
  "threshold sensitivity"
)
aggregation_error <- assert_same(
  independent_aggregation, saved_aggregation,
  c("aggregation", "final_lineage", "score"),
  setdiff(result_columns, c("HL_bootstrap_CI_low", "HL_bootstrap_CI_high")),
  "aggregation sensitivity"
)
loo_error <- assert_same(
  independent_loo, saved_loo,
  c("final_lineage", "score", "dropped_patient", "dropped_group"),
  c("n_MPR_pCR", "n_NMPR", "HL_shift_MPR_minus_NMPR"),
  "leave-one-patient-out"
)

figure_spec <- fread(file.path(table_dir, "07_figure_specifications.csv"))
dimension_checks <- figure_spec[, {
  path <- file.path(pub_dir, paste0(figure, ".png"))
  info <- attr(readPNG(path, native = TRUE, info = TRUE), "info")
  list(
    observed_width_px = as.integer(info$dim[1]),
    observed_height_px = as.integer(info$dim[2]),
    dimensions_match = as.integer(info$dim[1]) == expected_width_px && as.integer(info$dim[2]) == expected_height_px
  )
}, by = .(figure, expected_width_px, expected_height_px)]
if (!all(dimension_checks$dimensions_match)) stop("PNG dimension mismatch")

verification <- rbindlist(list(
  data.table(check = "Locked input and plan SHA-256", status = "PASS", detail = paste(names(observed_hashes), observed_hashes, collapse = "; ")),
  data.table(check = "Patient-level summary independent recomputation", status = "PASS", detail = sprintf("maximum absolute error %.3g", summary_error)),
  data.table(check = "Primary 39-test independent recomputation", status = "PASS", detail = sprintf("maximum absolute error %.3g", primary_error)),
  data.table(check = "Cell-count threshold sensitivity independent recomputation", status = "PASS", detail = sprintf("maximum absolute error %.3g", threshold_error)),
  data.table(check = "Alternative aggregation independent recomputation", status = "PASS", detail = sprintf("maximum absolute error %.3g", aggregation_error)),
  data.table(check = "Leave-one-patient-out independent recomputation", status = "PASS", detail = sprintf("maximum absolute error %.3g", loo_error)),
  data.table(check = "Primary BH-FDR discoveries", status = if (any(independent_primary$BH_q < 0.05, na.rm = TRUE)) "FAIL" else "PASS", detail = paste0("q<0.05: ", sum(independent_primary$BH_q < 0.05, na.rm = TRUE))),
  data.table(check = "PNG dimensions", status = "PASS", detail = paste0(nrow(dimension_checks), "/", nrow(dimension_checks), " match specifications"))
))
if (any(verification$status != "PASS")) stop("At least one verification check failed")
fwrite(verification, file.path(table_dir, "08_independent_verification.csv"))
fwrite(dimension_checks, file.path(table_dir, "08_png_dimension_verification.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "02_verification_sessionInfo.txt"))
message("Independent verification PASS")
