#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
.libPaths(c(file.path(analysis_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
})

options(stringsAsFactors = FALSE, warn = 1)
set.seed(260825)

table_dir <- file.path(analysis_dir, "results", "tables")
intermediate_dir <- file.path(analysis_dir, "intermediate")
log_dir <- file.path(analysis_dir, "logs")
counts_file <- file.path(intermediate_dir, "01_full_counts_mt20_qc.rds")
lineage_file <- file.path(table_dir, "09_final_cell_lineages_FROZEN.csv")
if (!all(file.exists(c(counts_file, lineage_file)))) stop("Required pseudobulk inputs are missing")

counts <- readRDS(counts_file)
lineages <- fread(lineage_file)
if (!inherits(counts, "dgCMatrix") || ncol(counts) != 92053L || nrow(lineages) != 92053L) {
  stop("Unexpected pseudobulk input dimensions")
}
ord <- match(colnames(counts), lineages$cell_id)
if (anyNA(ord)) stop("Count/lineage cell-ID mismatch")
lineages <- lineages[ord]
if (!identical(colnames(counts), lineages$cell_id)) stop("Count/lineage ordering failed")
if (!identical(as.numeric(Matrix::colSums(counts)), as.numeric(lineages$library_size))) {
  stop("Lineage-table library sizes disagree with counts")
}

primary_targets <- c("SARDH", "PIPOX", "SARDH_PIPOX_SUM")
secondary_declared_targets <- c("BHMT", "BHMT2", "DMGDH", "ETFB", "GNMT", "PIPOX",
                                "SARDH", "SHMT1", "SHMT2")
secondary_target_audit <- data.table(
  target = secondary_declared_targets,
  present_in_deposited_matrix = secondary_declared_targets %in% rownames(counts)
)
fwrite(secondary_target_audit, file.path(table_dir, "10_Lee_9gene_signature_availability.csv"))
secondary_targets <- secondary_target_audit[present_in_deposited_matrix == TRUE]$target
required_genes <- unique(c(setdiff(primary_targets, "SARDH_PIPOX_SUM"), secondary_targets))
if (!all(setdiff(primary_targets, "SARDH_PIPOX_SUM") %in% rownames(counts))) {
  stop("Missing target genes: ", paste(setdiff(required_genes, rownames(counts)), collapse = ", "))
}

lineages[, clinical_group := fifelse(
  Resource == "Pre-treatment biopsy", "TN",
  fifelse(Pathologic.Response %in% c("MPR", "pCR"), "MPR",
          fifelse(Pathologic.Response == "NMPR", "NMPR", NA_character_))
)]
if (anyNA(lineages$clinical_group)) stop("Unresolved clinical group in lineage table")

patient_info <- unique(lineages[, .(Sample, Patient, Resource, Pathologic.Response, clinical_group)])
if (nrow(patient_info) != 15L || anyDuplicated(patient_info$Sample) || anyDuplicated(patient_info$Patient)) {
  stop("Expected 15 unique patients/samples")
}
setorder(patient_info, Sample)
sample_ids <- patient_info$Sample

eligible_lineages <- sort(unique(lineages[analysis_eligible == TRUE]$final_lineage))
if (length(eligible_lineages) != 14L) stop("Expected 14 analysis-eligible lineages")
pseudobulk <- setNames(vector("list", length(eligible_lineages)), eligible_lineages)
cell_count_rows <- vector("list", length(eligible_lineages))

message("Aggregating full-transcriptome patient x lineage pseudobulks")
for (i in seq_along(eligible_lineages)) {
  lineage <- eligible_lineages[i]
  cell_idx <- which(lineages$analysis_eligible & lineages$final_lineage == lineage)
  membership <- sparseMatrix(
    i = seq_along(cell_idx), j = match(lineages$Sample[cell_idx], sample_ids), x = 1,
    dims = c(length(cell_idx), length(sample_ids)), dimnames = list(NULL, sample_ids)
  )
  pb <- counts[, cell_idx, drop = FALSE] %*% membership
  pb <- as(pb, "dgCMatrix")
  rownames(pb) <- rownames(counts)
  colnames(pb) <- sample_ids
  pseudobulk[[lineage]] <- pb
  cell_count_rows[[i]] <- data.table(
    final_lineage = lineage, Sample = sample_ids,
    cell_count = as.integer(Matrix::colSums(membership))
  )
  message("  ", lineage, ": ", length(cell_idx), " cells")
}
cell_counts <- rbindlist(cell_count_rows)
cell_counts <- merge(cell_counts, patient_info, by = "Sample", all.x = TRUE, sort = FALSE)
if (anyNA(cell_counts$Patient)) stop("Pseudobulk cell-count metadata merge failed")
fwrite(cell_counts, file.path(table_dir, "10_pseudobulk_cell_counts.csv"))
saveRDS(pseudobulk, file.path(intermediate_dir, "10_lineage_pseudobulk_counts.rds"), compress = "xz")

# Export patient-level target counts and TMM logCPM for audit and plotting.
target_value_rows <- vector("list", length(eligible_lineages))
normalization_rows <- vector("list", length(eligible_lineages))
for (i in seq_along(eligible_lineages)) {
  lineage <- eligible_lineages[i]
  pb <- pseudobulk[[lineage]]
  nonzero <- Matrix::colSums(pb) > 0
  if (!all(nonzero)) warning(lineage, " has zero-library pseudobulk(s); logCPM set to NA there")
  base_counts <- as.matrix(pb[required_genes, , drop = FALSE])
  composite <- colSums(base_counts[c("SARDH", "PIPOX"), , drop = FALSE])
  extended_counts <- rbind(base_counts, SARDH_PIPOX_SUM = composite)
  logcpm <- matrix(NA_real_, nrow = nrow(extended_counts), ncol = ncol(extended_counts),
                   dimnames = dimnames(extended_counts))
  norm_factor <- rep(NA_real_, ncol(pb))
  if (any(nonzero)) {
    y <- DGEList(counts = pb[, nonzero, drop = FALSE])
    y <- calcNormFactors(y, method = "TMM")
    y_ext <- DGEList(
      counts = extended_counts[, nonzero, drop = FALSE],
      lib.size = y$samples$lib.size,
      norm.factors = y$samples$norm.factors
    )
    logcpm[, nonzero] <- cpm(y_ext, log = TRUE, prior.count = 0.5)
    norm_factor[nonzero] <- y$samples$norm.factors
  }
  vals <- rbindlist(lapply(rownames(extended_counts), function(target) {
    data.table(
      final_lineage = lineage, Sample = sample_ids, target = target,
      raw_count = as.numeric(extended_counts[target, ]),
      TMM_log2_CPM = as.numeric(logcpm[target, ])
    )
  }))
  target_value_rows[[i]] <- vals
  normalization_rows[[i]] <- data.table(
    final_lineage = lineage, Sample = sample_ids,
    pseudobulk_library_size = as.numeric(Matrix::colSums(pb)),
    TMM_normalization_factor = norm_factor
  )
}
target_values <- rbindlist(target_value_rows)
target_values <- merge(target_values, cell_counts, by = c("final_lineage", "Sample"), all.x = TRUE, sort = FALSE)
fwrite(target_values, file.path(table_dir, "10_pseudobulk_target_counts_logCPM.csv"))
fwrite(rbindlist(normalization_rows), file.path(table_dir, "10_pseudobulk_TMM_factors.csv"))

fit_lineage_model <- function(lineage, threshold, drop_patient = NA_character_) {
  info <- cell_counts[final_lineage == lineage & clinical_group %in% c("MPR", "NMPR") &
                      cell_count >= threshold]
  if (!is.na(drop_patient)) info <- info[Patient != drop_patient]
  n_mpr <- sum(info$clinical_group == "MPR")
  n_nmpr <- sum(info$clinical_group == "NMPR")
  targets_all <- unique(c(primary_targets, secondary_targets))
  empty_rows <- data.table(
    final_lineage = lineage, cell_threshold = threshold, dropped_patient = drop_patient,
    target = targets_all, n_MPR = n_mpr, n_NMPR = n_nmpr,
    status = "INELIGIBLE_PATIENT_COUNTS", logFC_MPR_minus_NMPR = NA_real_,
    CI95_low = NA_real_, CI95_high = NA_real_, logCPM = NA_real_,
    F = NA_real_, df_total = NA_real_, PValue = NA_real_
  )
  if (n_mpr < 3L || n_nmpr < 6L) return(list(results = empty_rows, sample_factors = NULL))

  setorder(info, clinical_group, Patient)
  pb <- pseudobulk[[lineage]][, info$Sample, drop = FALSE]
  base <- as.matrix(pb)
  composite <- colSums(base[c("SARDH", "PIPOX"), , drop = FALSE])
  extended <- rbind(base, SARDH_PIPOX_SUM = composite)
  group <- factor(info$clinical_group, levels = c("NMPR", "MPR"))
  design <- model.matrix(~ group)
  y_all <- DGEList(counts = extended, lib.size = colSums(base), group = group)
  keep_expr <- filterByExpr(y_all, design = design, group = group)
  if (sum(keep_expr) < 100L) stop("Too few expressed genes for ", lineage, " threshold ", threshold)
  y <- y_all[keep_expr, , keep.lib.sizes = TRUE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = 2)
  tab <- qlf$table
  df_total <- qlf$df.total
  if (length(df_total) == 1L) df_total <- rep(df_total, nrow(tab))
  names(df_total) <- rownames(tab)

  rows <- rbindlist(lapply(targets_all, function(target) {
    if (!target %in% rownames(tab)) {
      return(data.table(
        final_lineage = lineage, cell_threshold = threshold, dropped_patient = drop_patient,
        target = target, n_MPR = n_mpr, n_NMPR = n_nmpr,
        status = "UNTESTABLE_LOW_EXPRESSION", logFC_MPR_minus_NMPR = NA_real_,
        CI95_low = NA_real_, CI95_high = NA_real_, logCPM = NA_real_,
        F = NA_real_, df_total = NA_real_, PValue = NA_real_
      ))
    }
    z <- tab[target, ]
    df <- as.numeric(df_total[target])
    se_from_f <- if (is.finite(z$F) && z$F > 0) abs(z$logFC) / sqrt(z$F) else NA_real_
    crit <- if (is.finite(df) && df > 0) qt(0.975, df = df) else NA_real_
    data.table(
      final_lineage = lineage, cell_threshold = threshold, dropped_patient = drop_patient,
      target = target, n_MPR = n_mpr, n_NMPR = n_nmpr, status = "TESTED",
      logFC_MPR_minus_NMPR = z$logFC,
      CI95_low = z$logFC - crit * se_from_f,
      CI95_high = z$logFC + crit * se_from_f,
      logCPM = z$logCPM, F = z$F, df_total = df, PValue = z$PValue
    )
  }))
  sample_factors <- data.table(
    final_lineage = lineage, cell_threshold = threshold, dropped_patient = drop_patient,
    Sample = info$Sample, Patient = info$Patient, clinical_group = info$clinical_group,
    cell_count = info$cell_count, library_size = y$samples$lib.size,
    TMM_normalization_factor = y$samples$norm.factors
  )
  list(results = rows, sample_factors = sample_factors)
}

message("Fitting prespecified patient-level edgeR QL models")
model_results <- list()
model_factors <- list()
k <- 0L
for (threshold in c(10L, 50L)) {
  for (lineage in eligible_lineages) {
    k <- k + 1L
    ans <- fit_lineage_model(lineage, threshold)
    model_results[[k]] <- ans$results
    model_factors[[k]] <- ans$sample_factors
  }
}
base_results <- rbindlist(model_results, fill = TRUE)
# SARDH and PIPOX belong to both the primary direct-degradation family and the
# secondary published nine-gene context. Duplicate those rows so multiplicity
# adjustment is correct within each explicitly labelled family.
all_results <- rbindlist(list(
  copy(base_results[target %in% primary_targets])[
    , family := "primary_direct_degradation"
  ],
  copy(base_results[target %in% secondary_targets])[
    , family := "secondary_Lee_signature_available_7_of_9"
  ]
), use.names = TRUE)

# Preserve explicit rows for the two genes absent from the deposited matrix.
missing_secondary <- secondary_target_audit[present_in_deposited_matrix == FALSE]$target
if (length(missing_secondary)) {
  missing_rows <- CJ(
    final_lineage = eligible_lineages,
    cell_threshold = c(10L, 50L),
    target = missing_secondary,
    unique = TRUE
  )
  missing_rows[, `:=`(
    dropped_patient = NA_character_, n_MPR = NA_integer_, n_NMPR = NA_integer_,
    status = "NOT_AVAILABLE_IN_DEPOSITED_MATRIX",
    logFC_MPR_minus_NMPR = NA_real_, CI95_low = NA_real_, CI95_high = NA_real_,
    logCPM = NA_real_, F = NA_real_, df_total = NA_real_, PValue = NA_real_,
    family = "secondary_Lee_signature_available_7_of_9"
  )]
  all_results <- rbindlist(list(all_results, missing_rows), use.names = TRUE, fill = TRUE)
}
all_results[, FDR_BH := NA_real_]
all_results[status == "TESTED" & family == "primary_direct_degradation",
            FDR_BH := p.adjust(PValue, method = "BH"), by = cell_threshold]
all_results[status == "TESTED" & family == "secondary_Lee_signature_available_7_of_9",
            FDR_BH := p.adjust(PValue, method = "BH"), by = cell_threshold]
fwrite(all_results[order(cell_threshold, family, final_lineage, target)],
       file.path(table_dir, "10_edger_target_results.csv"))
fwrite(rbindlist(model_factors, fill = TRUE), file.path(table_dir, "10_edger_model_sample_factors.csv"))

# Leave-one-patient-out direction stability for the primary threshold and family.
loo_rows <- list()
k <- 0L
for (lineage in eligible_lineages) {
  base_result <- all_results[cell_threshold == 10L & final_lineage == lineage &
                             family == "primary_direct_degradation" & status == "TESTED"]
  if (!nrow(base_result)) next
  eligible_patients <- cell_counts[
    final_lineage == lineage & clinical_group %in% c("MPR", "NMPR") & cell_count >= 10L,
    unique(Patient)
  ]
  for (patient in eligible_patients) {
    ans <- fit_lineage_model(lineage, 10L, drop_patient = patient)$results
    k <- k + 1L
    loo_rows[[k]] <- ans[target %in% primary_targets]
  }
}
loo <- rbindlist(loo_rows, fill = TRUE)
fwrite(loo[order(final_lineage, target, dropped_patient)],
       file.path(table_dir, "10_primary_leave_one_patient_out.csv"))

stability <- loo[status == "TESTED", .(
  LOO_folds = .N,
  positive_folds = sum(logFC_MPR_minus_NMPR > 0),
  negative_folds = sum(logFC_MPR_minus_NMPR < 0),
  min_logFC = min(logFC_MPR_minus_NMPR),
  max_logFC = max(logFC_MPR_minus_NMPR)
), by = .(final_lineage, target)]
primary10 <- all_results[cell_threshold == 10L & family == "primary_direct_degradation"]
primary50 <- all_results[cell_threshold == 50L & family == "primary_direct_degradation",
                         .(final_lineage, target, status_50 = status, logFC_50 = logFC_MPR_minus_NMPR,
                           PValue_50 = PValue, FDR_BH_50 = FDR_BH, n_MPR_50 = n_MPR, n_NMPR_50 = n_NMPR)]
summary <- merge(primary10, primary50, by = c("final_lineage", "target"), all.x = TRUE)
summary <- merge(summary, stability, by = c("final_lineage", "target"), all.x = TRUE)
fwrite(summary[order(final_lineage, target)],
       file.path(table_dir, "10_primary_stability_summary.csv"))

capture.output(sessionInfo(), file = file.path(log_dir, "10_sessionInfo.txt"))
cat("Sarcosine lineage pseudobulk analysis PASS\n")
cat("Eligible lineages aggregated:", length(eligible_lineages), "\n")
cat("Primary threshold-10 tested rows:",
    nrow(all_results[cell_threshold == 10L & family == "primary_direct_degradation" & status == "TESTED"]), "\n")
cat("Primary threshold-50 tested rows:",
    nrow(all_results[cell_threshold == 50L & family == "primary_direct_degradation" & status == "TESTED"]), "\n")
