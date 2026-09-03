#!/usr/bin/env Rscript

# Matched ccRCC metabolomics/RNA-seq: final exploratory prioritization of
# immune programs associated with lower measured tumour Sarcosine and the
# IFNG-response program.
#
# Scientific boundary
# -------------------
# - Sarcosine is a normalized untargeted GC-MS tissue signal, not an absolute
#   concentration, cell-specific measurement, or metabolic-flux estimate.
# - Candidate programs were selected after prior ccRCC GSEA and are therefore
#   explicitly post-hoc/exploratory. This script ranks internal consistency;
#   it does not produce confirmatory causal evidence.
# - Continuous log2 Sarcosine is the primary exposure. Median Low/High and
#   Q1/Q4 are reader-facing/sensitivity contrasts.
# - Cross-sectional bulk RNA cannot establish temporal direction, same-cell
#   signalling, protein activation, or cytokine secretion.

options(stringsAsFactors = FALSE, width = 180)
set.seed(260902)

analysis_root <- normalizePath(getwd())
# Unicode normalization differs between shells/filesystems; verify the project
# by required subdirectories instead of comparing the rendered basename.
if (!all(dir.exists(file.path(analysis_root, c("RNAseq_Data", "results", "analysis"))))) {
  stop("Run from the matched ccRCC analysis root: ", analysis_root)
}

project_root <- normalizePath(file.path(analysis_root, ".."))
default_gsea_lib <- file.path(
  project_root, "2024_Drug_Res_Updates_NSCLC", "RNA-seq공공데이터_GSE207422",
  "analysis_sarcosine_FINAL_26.08.25", "03_lee_fig3_style_FINAL",
  "hallmark_GSEA_Q4_vs_Q1_26.08.26", "R_libs"
)
gsea_lib <- Sys.getenv("SARCO_GSEA_R_LIB", unset = default_gsea_lib)
if (!dir.exists(gsea_lib)) stop("GSEA R library not found: ", gsea_lib)
.libPaths(unique(c(normalizePath(gsea_lib), .libPaths())))

required_packages <- c(
  "data.table", "limma", "fgsea", "msigdbr", "BiocParallel",
  "ggplot2", "patchwork", "glmnet", "pheatmap"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
  library(glmnet)
  library(pheatmap)
})

expression_path <- file.path(analysis_root, "RNAseq_Data", "bulkRNA_matrix_TPM.csv")
metadata_path <- file.path(
  analysis_root, "results", "00_input_audit_26.09.02", "tumor_sarcosine_group_map.csv"
)
quartile_path <- file.path(
  analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_Q4_vs_Q1_Extreme_26.09.02",
  "tables", "Q1_Q4_all_tumour_group_map.csv"
)
deconv_path <- file.path(
  analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_HighLow_Deconvolution_26.09.02",
  "tables", "core_celltype_values_long.csv"
)
regulon_path <- file.path(
  analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_HighLow_Immune_TF_Integrated_26.09.02",
  "tables", "DoRothEA_AB_signed_regulon_tested.csv"
)
genomewide_median_path <- file.path(
  analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_HighLow_Transcriptome_GSEA_26.09.02",
  "tables", "GSEA_Hallmark_Reactome_GOBP_adjusted_High_vs_Low.csv"
)
genomewide_q4_path <- file.path(
  analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_Q4_vs_Q1_CD28_IFNG_Followup_26.09.02",
  "tables", "03_GSEA_all_collections_all_adjustments.csv"
)
genomewide_q4_exclude_path <- file.path(
  analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_Q4_vs_Q1_Extreme_26.09.02",
  "tables", "GSEA_Hallmark_Reactome_GOBP_adjusted_Q4_vs_Q1_exclude5_Tukey_low.csv"
)
required_inputs <- c(
  expression_path, metadata_path, quartile_path, deconv_path, regulon_path,
  genomewide_median_path, genomewide_q4_path, genomewide_q4_exclude_path
)
for (path in required_inputs) if (!file.exists(path)) stop("Missing input: ", path)

output_root <- file.path(
  analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_Immune_Pathway_Final_Narrowing_26.09.02"
)
table_dir <- file.path(output_root, "tables")
figure_dir <- file.path(output_root, "figures")
log_dir <- file.path(output_root, "logs")
for (d in c(table_dir, figure_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) {
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  strsplit(out, "[[:space:]]+")[[1]][1]
}

zscore <- function(x) {
  out <- as.numeric(scale(as.numeric(x)))
  if (anyNA(out) || !all(is.finite(out))) stop("Non-finite z score")
  out
}

score_gene_set <- function(genes, label, log_expression, minimum_genes = 10L) {
  mapped <- sort(intersect(unique(genes), rownames(log_expression)))
  mapped <- mapped[apply(log_expression[mapped, , drop = FALSE], 1L, var) > 0]
  if (length(mapped) < minimum_genes) {
    stop(label, " has only ", length(mapped), " usable genes; minimum is ", minimum_genes)
  }
  gene_z <- t(scale(t(log_expression[mapped, , drop = FALSE])))
  gene_z[!is.finite(gene_z)] <- 0
  score <- zscore(colMeans(gene_z))
  names(score) <- colnames(log_expression)
  list(score = score, genes = mapped)
}

extract_lm_term <- function(fit, term, outcome, predictor, model, route = NA_character_, edge = NA_character_) {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) stop("Term absent from fitted model: ", term)
  estimate <- unname(sm[term, "Estimate"])
  se <- unname(sm[term, "Std. Error"])
  df <- df.residual(fit)
  crit <- qt(0.975, df)
  data.table(
    route = route, edge = edge, outcome = outcome, predictor = predictor, model = model,
    standardized_beta = estimate, standard_error = se,
    CI_low = estimate - crit * se, CI_high = estimate + crit * se,
    t_value = unname(sm[term, "t value"]), residual_df = df,
    p_value = unname(sm[term, "Pr(>|t|)"]),
    partial_R2 = unname(sm[term, "t value"])^2 /
      (unname(sm[term, "t value"])^2 + df)
  )
}

ulm_from_matrix <- function(mat, net, minsize = 10L, center_rows = TRUE) {
  stopifnot(is.matrix(mat), all(c("TF", "target", "mor") %in% names(net)))
  shared_targets <- sort(intersect(rownames(mat), unique(net$target)))
  net <- net[target %in% shared_targets]
  keep_tf <- net[, uniqueN(target), by = TF][V1 >= minsize, TF]
  net <- net[TF %in% keep_tf]
  shared_targets <- sort(unique(net$target))
  tfs <- sort(unique(net$TF))
  mor_mat <- matrix(0, nrow = length(shared_targets), ncol = length(tfs),
                    dimnames = list(shared_targets, tfs))
  mor_mat[cbind(match(net$target, shared_targets), match(net$TF, tfs))] <- net$mor
  y <- mat[shared_targets, , drop = FALSE]
  if (center_rows) y <- sweep(y, 1L, rowMeans(y), FUN = "-")
  r <- cor(mor_mat, y)
  r <- pmax(pmin(r, 1 - 1e-12), -1 + 1e-12)
  df <- nrow(mor_mat) - 2L
  score <- r * sqrt(df / ((1 - r + 1e-20) * (1 + r + 1e-20)))
  if (anyNA(score) || !all(is.finite(score))) stop("Non-finite ULM score")
  score
}

prepare_covariates <- function(d) {
  d <- copy(d)
  d[, `:=`(
    batch = factor(batch),
    sex = relevel(factor(sex), ref = "female"),
    age = factor(age, levels = c("40-60", "<40", ">60")),
    grade = factor(grade, levels = c("1", "2", "3", "4"))
  )]
  d
}

fit_rank <- function(d, expression, extra_terms = character(), label) {
  d <- droplevels(copy(d))
  rhs <- c("batch", "sex", "age", "grade", extra_terms, "low_exposure_z")
  design <- model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), data = d)
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient design: ", label)
  fit <- eBayes(lmFit(expression[, d$sample_id, drop = FALSE], design), trend = TRUE, robust = TRUE)
  tt <- as.data.table(
    topTable(fit, coef = "low_exposure_z", number = Inf, sort.by = "none", adjust.method = "BH"),
    keep.rownames = "gene_symbol"
  )
  setnames(tt, c("logFC", "t", "P.Value", "adj.P.Val"),
           c("effect_per_SD", "moderated_t", "p_value", "BH_q"))
  stopifnot(nrow(tt) == nrow(expression), all(is.finite(tt$moderated_t)))
  list(table = tt, rank = setNames(tt$moderated_t, tt$gene_symbol), design = design)
}

run_focused_gsea <- function(stats, pathways, label) {
  ans <- as.data.table(fgseaMultilevel(
    pathways = pathways, stats = stats, minSize = 10L, maxSize = 500L,
    eps = 0, nPermSimple = 100000L, nproc = 1L,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  ))
  ans[, `:=`(
    model = label,
    direction = fifelse(is.finite(NES) & NES > 0, "Sarcosine-Low",
                        fifelse(is.finite(NES), "Sarcosine-High", "Not estimable"))
  )]
  ans
}

fit_score_effects <- function(d, score_matrix, extra_terms = character(), label) {
  d <- droplevels(copy(d))
  rhs <- c("batch", "sex", "age", "grade", extra_terms, "low_exposure_z")
  design <- model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), data = d)
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient score design: ", label)
  fit <- eBayes(lmFit(score_matrix[, d$sample_id, drop = FALSE], design), robust = TRUE)
  term <- "low_exposure_z"
  se <- fit$stdev.unscaled[, term] * sqrt(fit$s2.post)
  df <- fit$df.total
  out <- data.table(
    route = rownames(score_matrix), model = label,
    standardized_beta = fit$coefficients[, term], standard_error = se,
    moderated_t = fit$t[, term], p_value = fit$p.value[, term], residual_df = df
  )
  out[, `:=`(
    CI_low = standardized_beta - qt(0.975, residual_df) * standard_error,
    CI_high = standardized_beta + qt(0.975, residual_df) * standard_error,
    BH_q = p.adjust(p_value, method = "BH")
  )]
  out[]
}

fit_joint_competition <- function(d, outcome, route_catalog, extra_terms = character(), label) {
  d <- droplevels(copy(d))
  rhs <- c("low_exposure_z", route_catalog$column, "batch", "sex", "age", "grade", extra_terms)
  form <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
  design <- model.matrix(form, data = d)
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient joint model: ", label)
  fit <- lm(form, data = d)
  out <- rbindlist(lapply(seq_len(nrow(route_catalog)), function(i) {
    extract_lm_term(
      fit, route_catalog$column[i], outcome, route_catalog$route[i], label,
      route_catalog$route[i], "Program -> IFNG response | all programs"
    )
  }))
  out[, BH_q := p.adjust(p_value, method = "BH")]
  out
}

joint_bootstrap <- function(d, outcome, route_catalog, extra_terms = character(),
                            label, B = 5000L, seed = 260902L) {
  d <- droplevels(copy(d))
  rhs <- c("low_exposure_z", route_catalog$column, "batch", "sex", "age", "grade", extra_terms)
  x <- model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), data = d)
  y <- d[[outcome]]
  idx_terms <- match(route_catalog$column, colnames(x))
  observed <- coef(lm.fit(x, y))[idx_terms]
  boot <- matrix(NA_real_, nrow = B, ncol = nrow(route_catalog))
  set.seed(seed)
  for (b in seq_len(B)) {
    idx <- sample.int(nrow(d), nrow(d), replace = TRUE)
    cf <- tryCatch(coef(lm.fit(x[idx, , drop = FALSE], y[idx])),
                   error = function(e) rep(NA_real_, ncol(x)))
    boot[b, ] <- cf[idx_terms]
  }
  rbindlist(lapply(seq_len(nrow(route_catalog)), function(i) {
    vals <- boot[, i]
    vals <- vals[is.finite(vals)]
    if (length(vals) < 0.95 * B) stop("Too many invalid joint-bootstrap replicates")
    ci <- quantile(vals, c(0.025, 0.975), names = FALSE)
    data.table(
      route = route_catalog$route[i], model = label, observed_beta = observed[i],
      bootstrap_CI_low = ci[1], bootstrap_CI_high = ci[2],
      positive_fraction = mean(vals > 0), negative_fraction = mean(vals < 0),
      CI_excludes_zero = ci[1] > 0 | ci[2] < 0,
      bootstrap_replicates_requested = B, bootstrap_replicates_valid = length(vals)
    )
  }))
}

serial_bootstrap <- function(d, mediator, outcome, route, label,
                             extra_terms = character(), B = 5000L, seed = 260902L,
                             reverse = FALSE) {
  d <- droplevels(copy(d))
  covars <- c("batch", "sex", "age", "grade", extra_terms)
  m <- if (!reverse) mediator else outcome
  y <- if (!reverse) outcome else mediator
  x1 <- model.matrix(as.formula(paste("~ low_exposure_z +", paste(covars, collapse = " + "))), data = d)
  x2 <- model.matrix(as.formula(paste("~ low_exposure_z +", m, "+", paste(covars, collapse = " + "))), data = d)
  y1 <- d[[m]]; y2 <- d[[y]]
  ia <- match("low_exposure_z", colnames(x1)); ib <- match(m, colnames(x2))
  one_product <- function(idx) {
    coef(lm.fit(x1[idx, , drop = FALSE], y1[idx]))[ia] *
      coef(lm.fit(x2[idx, , drop = FALSE], y2[idx]))[ib]
  }
  observed <- one_product(seq_len(nrow(d)))
  set.seed(seed)
  vals <- replicate(B, one_product(sample.int(nrow(d), nrow(d), replace = TRUE)))
  vals <- vals[is.finite(vals)]
  if (length(vals) < 0.95 * B) stop("Too many invalid serial-bootstrap replicates")
  ci <- quantile(vals, c(0.025, 0.975), names = FALSE)
  data.table(
    route = route, model = label, order = ifelse(reverse, "reverse", "forward"),
    serial_product = observed, bootstrap_CI_low = ci[1], bootstrap_CI_high = ci[2],
    positive_fraction = mean(vals > 0), negative_fraction = mean(vals < 0),
    CI_excludes_zero = ci[1] > 0 | ci[2] < 0,
    bootstrap_replicates_requested = B, bootstrap_replicates_valid = length(vals)
  )
}

make_stratified_folds <- function(group, k = 5L) {
  fold <- integer(length(group))
  for (lev in unique(group)) {
    idx <- sample(which(group == lev))
    fold[idx] <- rep(seq_len(k), length.out = length(idx))
  }
  fold
}

repeated_cv <- function(d, outcome, route_catalog, repeats = 100L, k = 5L, seed = 261102L) {
  d <- droplevels(copy(d))
  base_rhs <- c("low_exposure_z", "batch", "sex", "age", "grade")
  model_catalog <- rbindlist(list(
    data.table(model = "Baseline", route = "Baseline", rhs = paste(base_rhs, collapse = " + ")),
    route_catalog[, .(model = route, route = route,
                      rhs = vapply(column, function(z) paste(c(base_rhs, z), collapse = " + "), character(1)))],
    data.table(model = "All programs", route = "All programs",
               rhs = paste(c(base_rhs, route_catalog$column), collapse = " + "))
  ))
  x_list <- setNames(lapply(model_catalog$rhs, function(rhs) {
    model.matrix(as.formula(paste("~", rhs)), data = d)[, -1, drop = FALSE]
  }), model_catalog$model)
  y <- d[[outcome]]
  set.seed(seed)
  out <- vector("list", repeats)
  for (r in seq_len(repeats)) {
    fold <- make_stratified_folds(d$sarcosine_group, k)
    pred <- matrix(NA_real_, nrow(d), nrow(model_catalog), dimnames = list(NULL, model_catalog$model))
    for (f in seq_len(k)) {
      test <- which(fold == f); train <- setdiff(seq_len(nrow(d)), test)
      for (m in seq_len(nrow(model_catalog))) {
        x <- x_list[[model_catalog$model[m]]]
        inner_fold <- make_stratified_folds(d$sarcosine_group[train], 5L)
        fit <- tryCatch(
          cv.glmnet(
            x[train, , drop = FALSE], y[train], family = "gaussian", alpha = 0,
            nfolds = 5L, foldid = inner_fold, type.measure = "mse", standardize = TRUE
          ), error = function(e) NULL
        )
        if (is.null(fit)) stop("Nested ridge CV failed in repeat ", r, ", fold ", f)
        pred[test, m] <- as.numeric(predict(fit, newx = x[test, , drop = FALSE], s = "lambda.min"))
      }
    }
    sst <- sum((y - mean(y))^2)
    perf <- rbindlist(lapply(seq_len(nrow(model_catalog)), function(m) {
      err <- y - pred[, m]
      data.table(
        repeat_id = r, model = model_catalog$model[m], route = model_catalog$route[m],
        RMSE = sqrt(mean(err^2)), R2 = 1 - sum(err^2) / sst
      )
    }))
    perf[, `:=`(
      delta_R2_vs_baseline = R2 - R2[model == "Baseline"],
      delta_RMSE_vs_baseline = RMSE - RMSE[model == "Baseline"]
    )]
    out[[r]] <- perf
  }
  rbindlist(out)
}

elastic_net_stability <- function(d, outcome, route_catalog, repetitions = 1000L,
                                  fraction = 0.80, seed = 261202L) {
  d <- droplevels(copy(d))
  rhs <- c("low_exposure_z", route_catalog$column, "batch", "sex", "age", "grade")
  x <- model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), data = d)[, -1, drop = FALSE]
  y <- d[[outcome]]
  penalty <- ifelse(colnames(x) %in% route_catalog$column, 1, 0)
  coef_1se <- matrix(NA_real_, repetitions, nrow(route_catalog))
  coef_min <- matrix(NA_real_, repetitions, nrow(route_catalog))
  set.seed(seed)
  for (b in seq_len(repetitions)) {
    idx <- sort(unlist(lapply(levels(d$sarcosine_group), function(lev) {
      candidates <- which(d$sarcosine_group == lev)
      sample(candidates, max(10L, floor(length(candidates) * fraction)), replace = FALSE)
    })))
    foldid <- make_stratified_folds(d$sarcosine_group[idx], 5L)
    fit <- tryCatch(
      cv.glmnet(
        x[idx, , drop = FALSE], y[idx], family = "gaussian", alpha = 0.5,
        nfolds = 5L, foldid = foldid, type.measure = "mse", standardize = TRUE,
        penalty.factor = penalty
      ), error = function(e) NULL
    )
    if (!is.null(fit)) {
      c1 <- as.matrix(coef(fit, s = "lambda.1se"))[, 1]
      cm <- as.matrix(coef(fit, s = "lambda.min"))[, 1]
      coef_1se[b, ] <- c1[route_catalog$column]
      coef_min[b, ] <- cm[route_catalog$column]
    }
  }
  valid <- rowSums(is.finite(coef_1se)) == nrow(route_catalog)
  if (sum(valid) < 0.95 * repetitions) stop("Too many invalid elastic-net repetitions")
  rbindlist(lapply(seq_len(nrow(route_catalog)), function(i) {
    x1 <- coef_1se[valid, i]; xm <- coef_min[valid, i]
    data.table(
      route = route_catalog$route[i],
      lambda_1se_selection_frequency = mean(abs(x1) > 1e-10),
      lambda_1se_positive_frequency = mean(x1 > 1e-10),
      lambda_1se_negative_frequency = mean(x1 < -1e-10),
      lambda_min_selection_frequency = mean(abs(xm) > 1e-10),
      lambda_min_positive_frequency = mean(xm > 1e-10),
      repetitions_requested = repetitions, repetitions_valid = sum(valid)
    )
  }))
}

vif_one <- function(term, others, d) {
  fit <- lm(as.formula(paste(term, "~", paste(others, collapse = " + "))), data = d)
  1 / (1 - summary(fit)$r.squared)
}

# -------------------------------------------------------------------------
# Data contract and frozen sample definitions
# -------------------------------------------------------------------------
meta <- fread(metadata_path, check.names = FALSE)
quartile <- fread(quartile_path, check.names = FALSE)
stopifnot(
  nrow(meta) == 100L, uniqueN(meta$sample_id) == 100L,
  nrow(quartile) == 100L, uniqueN(quartile$sample_id) == 100L,
  identical(meta$sample_id, quartile$sample_id),
  max(abs(meta$sarcosine_normalized_intensity - quartile$sarcosine_normalized_intensity)) < 1e-12
)
meta <- prepare_covariates(meta)
meta[, sarcosine_group := factor(sarcosine_group, levels = c("Low", "High"))]
sample_ids <- meta$sample_id

expr_dt <- fread(expression_path, check.names = FALSE)
gene_symbols <- trimws(expr_dt[[1]])
if (anyNA(gene_symbols) || any(!nzchar(gene_symbols)) || anyDuplicated(gene_symbols)) stop("Invalid gene symbols")
if (!all(sample_ids %in% names(expr_dt))) stop("Expression matrix lacks matched tumour samples")
tpm <- as.matrix(expr_dt[, ..sample_ids])
storage.mode(tpm) <- "double"
rownames(tpm) <- gene_symbols; colnames(tpm) <- sample_ids
if (anyNA(tpm) || any(!is.finite(tpm)) || any(tpm < 0)) stop("Invalid TPM matrix")
keep <- rowSums(tpm >= 1) >= 10L & apply(tpm, 1L, var) > 0
log_expression <- log2(tpm[keep, , drop = FALSE] + 1)
if (nrow(log_expression) != 15119L) stop("Unexpected expression universe: ", nrow(log_expression))

deconv <- fread(deconv_path, check.names = FALSE)
context <- dcast(
  deconv[method == "MCPcounter" & cell_type == "Pan T cells" |
           method == "ESTIMATE" & cell_type %in% c("Immune score", "Stromal score"),
         .(sample_id, context_name = fifelse(method == "MCPcounter", "Pan_T",
                                             fifelse(cell_type == "Immune score", "Immune", "Stromal")),
           value = analysis_value)],
  sample_id ~ context_name, value.var = "value"
)
if (nrow(context) != 100L || anyNA(context)) stop("Incomplete cell-context scores")
context[, `:=`(Pan_T = zscore(Pan_T), Immune = zscore(Immune), Stromal = zscore(Stromal))]
meta <- merge(meta, context, by = "sample_id", all.x = TRUE, sort = FALSE)
meta <- meta[match(sample_ids, sample_id)]

q1 <- quantile(meta$sarcosine_normalized_intensity, 0.25)
q3 <- quantile(meta$sarcosine_normalized_intensity, 0.75)
tukey_low_cut <- q1 - 1.5 * IQR(meta$sarcosine_normalized_intensity)
outlier_ids <- sort(meta[sarcosine_normalized_intensity < tukey_low_cut, sample_id])
expected_outliers <- sort(c("H46_T", "N39_T", "R25_T", "R82_T", "Z16_T"))
if (!identical(outlier_ids, expected_outliers)) stop("Unexpected Tukey-low sample set")
win_lo <- quantile(meta$log2_sarcosine_intensity, 0.05)
win_hi <- quantile(meta$log2_sarcosine_intensity, 0.95)

# Positive exposure values always mean lower measured Sarcosine.
model_data <- list()
d <- copy(meta); d[, low_exposure_z := zscore(as.integer(sarcosine_group == "Low"))]
model_data[["Median Low vs High"]] <- list(data = d, extras = character(), display = TRUE)
d <- merge(meta, quartile[, .(sample_id, extreme_group)], by = "sample_id", all.x = TRUE, sort = FALSE)
d <- d[match(sample_ids, sample_id)]
d <- d[extreme_group != "Middle 50% excluded"]
d[, sarcosine_group := factor(fifelse(extreme_group == "Sarcosine-Low (Q1)", "Low", "High"), levels = c("Low", "High"))]
d[, low_exposure_z := zscore(as.integer(sarcosine_group == "Low"))]
model_data[["Q1 vs Q4"]] <- list(data = d, extras = character(), display = TRUE)
model_data[["Q1 vs Q4, pan-T sensitivity"]] <- list(data = copy(d), extras = "Pan_T", display = FALSE)
model_data[["Q1 vs Q4, broad TME sensitivity"]] <- list(data = copy(d), extras = c("Immune", "Stromal"), display = FALSE)
d_no <- d[!sample_id %in% outlier_ids]
d_no[, low_exposure_z := zscore(as.integer(sarcosine_group == "Low"))]
model_data[["Q1 vs Q4, exclude 5"]] <- list(data = d_no, extras = character(), display = TRUE)
d <- copy(meta); d[, low_exposure_z := zscore(-log2_sarcosine_intensity)]
model_data[["Continuous lower Sarcosine"]] <- list(data = d, extras = character(), display = FALSE)
d <- copy(meta[!sample_id %in% outlier_ids]); d[, low_exposure_z := zscore(-log2_sarcosine_intensity)]
model_data[["Continuous lower Sarcosine, exclude 5"]] <- list(data = d, extras = character(), display = TRUE)
d <- copy(meta)
d[, log2_winsorized := pmin(pmax(log2_sarcosine_intensity, win_lo), win_hi)]
d[, low_exposure_z := zscore(-log2_winsorized)]
model_data[["Continuous lower Sarcosine, winsorized"]] <- list(data = d, extras = character(), display = TRUE)

model_manifest <- rbindlist(lapply(names(model_data), function(label) {
  z <- model_data[[label]]
  data.table(
    model = label, n = nrow(z$data),
    Low_n = sum(z$data$sarcosine_group == "Low"), High_n = sum(z$data$sarcosine_group == "High"),
    exposure = ifelse(grepl("Continuous", label), "z(-log2 Sarcosine)", "z(Low-group indicator)"),
    additional_covariates = ifelse(length(z$extras), paste(z$extras, collapse = ";"), "none"),
    reader_facing = z$display
  )
}))

# -------------------------------------------------------------------------
# Candidate programs and strictly gene-disjoint scores
# -------------------------------------------------------------------------
candidate_catalog <- data.table(
  route_order = 1:6,
  route = c(
    "TNF/NF-kappaB", "Type-I-IFN/JAK-STAT", "IL6/JAK-STAT3",
    "T-cell activation/proliferation", "Leukocyte cytotoxicity", "IL12-family/STAT4"
  ),
  anchor_pathway = c(
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING", "GOBP_T_CELL_PROLIFERATION",
    "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY", "REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING"
  ),
  column = paste0("program_", 1:6),
  role = c(
    "Inflammatory signaling candidate", "Parallel interferon candidate", "Cytokine signaling candidate",
    "T-cell-state candidate", "Effector-state candidate", "Parallel cytokine candidate"
  )
)

msig <- as.data.table(msigdbr(species = "Homo sapiens", db_species = "HS"))
msigdb_version <- unique(msig$db_version)
if (length(msigdb_version) != 1L) stop("Non-unique MSigDB version")
if (!all(candidate_catalog$anchor_pathway %in% msig$gs_name)) stop("Missing candidate anchor pathway")

response_ids <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "GOBP_RESPONSE_TO_TYPE_II_INTERFERON"
)
if (!all(response_ids %in% msig$gs_name)) stop("Missing IFNG-response source pathway")
response_core <- msig[gs_name %in% response_ids, .N, by = gene_symbol][N >= 2L, gene_symbol]
response_core <- setdiff(intersect(response_core, rownames(log_expression)), "IFNG")
if (length(response_core) < 50L) stop("Unexpectedly small IFNG-response core")

raw_candidates <- setNames(lapply(candidate_catalog$anchor_pathway, function(pid) {
  intersect(unique(msig[gs_name == pid, gene_symbol]), rownames(log_expression))
}), candidate_catalog$route)
without_outcome <- lapply(raw_candidates, setdiff, y = response_core)
membership_count <- table(unlist(without_outcome, use.names = FALSE))
shared_candidate_genes <- names(membership_count[membership_count > 1L])
strict_candidates <- lapply(without_outcome, setdiff, y = shared_candidate_genes)
if (any(lengths(strict_candidates) < 10L)) stop("A strict candidate program has fewer than 10 genes")

program_objects <- lapply(names(strict_candidates), function(route) {
  score_gene_set(strict_candidates[[route]], route, log_expression, 10L)
})
names(program_objects) <- names(strict_candidates)
ifng_outcome <- score_gene_set(response_core, "IFNG-response core", log_expression, 30L)

analysis_data <- copy(meta)
for (i in seq_len(nrow(candidate_catalog))) {
  analysis_data[, (candidate_catalog$column[i]) := program_objects[[candidate_catalog$route[i]]]$score]
}
analysis_data[, `:=`(
  IFNG_response_z = ifng_outcome$score,
  IFNG_expression_z = zscore(log_expression["IFNG", sample_ids])
)]
if (anyNA(analysis_data)) stop("Missing sample-level analysis value")

program_columns <- candidate_catalog$column
program_matrix <- t(as.matrix(analysis_data[, ..program_columns]))
rownames(program_matrix) <- candidate_catalog$route
colnames(program_matrix) <- sample_ids

score_overlap <- rbindlist(lapply(names(strict_candidates), function(a) {
  rbindlist(lapply(names(strict_candidates), function(b) {
    data.table(score_1 = a, score_2 = b,
               overlap_n = length(intersect(strict_candidates[[a]], strict_candidates[[b]])))
  }))
}))
if (score_overlap[score_1 != score_2, max(overlap_n)] != 0L) stop("Candidate scores overlap")
if (max(vapply(strict_candidates, function(g) length(intersect(g, response_core)), integer(1))) != 0L) {
  stop("Candidate/outcome overlap remains")
}

score_manifest <- rbindlist(lapply(seq_len(nrow(candidate_catalog)), function(i) {
  route <- candidate_catalog$route[i]
  data.table(
    route = route, anchor_pathway = candidate_catalog$anchor_pathway[i],
    raw_expressed_genes_n = length(raw_candidates[[route]]),
    removed_IFNG_core_n = length(intersect(raw_candidates[[route]], response_core)),
    removed_candidate_overlap_n = length(intersect(without_outcome[[route]], shared_candidate_genes)),
    scored_genes_n = length(program_objects[[route]]$genes),
    scored_genes = paste(program_objects[[route]]$genes, collapse = ";")
  )
}))

# Focused GSEA uses the full, unmodified anchor sets; BH q is across six
# predeclared post-hoc candidates and is not substituted for genome-wide q.
focus_pathways <- setNames(lapply(candidate_catalog$anchor_pathway, function(pid) {
  unique(msig[gs_name == pid, gene_symbol])
}), candidate_catalog$route)
rank_fits <- lapply(names(model_data), function(label) {
  z <- model_data[[label]]
  fit_rank(z$data, log_expression, z$extras, label)
})
names(rank_fits) <- names(model_data)
focused_gsea <- rbindlist(lapply(names(rank_fits), function(label) {
  run_focused_gsea(rank_fits[[label]]$rank, focus_pathways, label)
}), use.names = TRUE)
focused_gsea[, route := pathway]
focused_gsea[, pathway := candidate_catalog$anchor_pathway[match(route, candidate_catalog$route)]]
focused_gsea_out <- copy(focused_gsea)
focused_gsea_out[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]

# Sample-score effects for the same models.
score_effects <- rbindlist(lapply(names(model_data), function(label) {
  z <- model_data[[label]]
  fit_score_effects(z$data, program_matrix, z$extras, label)
}))

# Genome-wide-q anchors from already validated full-collection analyses.
gw_median <- fread(genomewide_median_path, check.names = FALSE)
gw_q4 <- fread(genomewide_q4_path, check.names = FALSE)
gw_q4_ex <- fread(genomewide_q4_exclude_path, check.names = FALSE)
extract_gw <- function(x, model_name, source_model = NULL) {
  if (!is.null(source_model)) x <- x[model == source_model]
  out <- x[pathway %in% candidate_catalog$anchor_pathway,
           .(pathway, NES, genomewide_BH_q = padj, direction)]
  if (nrow(out) != nrow(candidate_catalog)) stop("Incomplete genome-wide anchor extraction: ", model_name)
  out[, `:=`(route = candidate_catalog$route[match(pathway, candidate_catalog$anchor_pathway)], model = model_name)]
  out
}
genomewide_anchors <- rbindlist(list(
  extract_gw(gw_median, "Median Low vs High"),
  extract_gw(gw_q4, "Q1 vs Q4", "Clinical covariates"),
  extract_gw(gw_q4, "Q1 vs Q4, pan-T sensitivity", "+ Pan-T estimate"),
  extract_gw(gw_q4, "Q1 vs Q4, broad TME sensitivity", "+ Immune/stromal scores"),
  extract_gw(gw_q4_ex, "Q1 vs Q4, exclude 5")
), use.names = TRUE)
genomewide_anchors[, low_direction := NES > 0 | grepl("Low|Q1", direction)]
# Original outputs encode High-minus-Low NES, so normalize displayed sign.
genomewide_anchors[grepl("Low|Q1", direction), normalized_low_NES := abs(NES)]
genomewide_anchors[!grepl("Low|Q1", direction), normalized_low_NES := -abs(NES)]

# -------------------------------------------------------------------------
# Competing programs against a strictly non-overlapping IFNG-response outcome
# -------------------------------------------------------------------------
for (label in names(model_data)) {
  z <- model_data[[label]]$data
  cols <- c("sample_id", "low_exposure_z", "sarcosine_group")
  if ("extreme_group" %in% names(z)) cols <- c(cols, "extreme_group")
  merged <- merge(z, analysis_data[, c("sample_id", candidate_catalog$column,
                                       "IFNG_response_z", "IFNG_expression_z"), with = FALSE],
                  by = "sample_id", all.x = TRUE, sort = FALSE)
  model_data[[label]]$data <- merged[match(z$sample_id, sample_id)]
}

primary_label <- "Continuous lower Sarcosine, exclude 5"
winsor_label <- "Continuous lower Sarcosine, winsorized"
primary_data <- model_data[[primary_label]]$data

joint_results <- rbindlist(lapply(c(primary_label, winsor_label, "Median Low vs High", "Q1 vs Q4"), function(label) {
  z <- model_data[[label]]
  fit_joint_competition(z$data, "IFNG_response_z", candidate_catalog, z$extras, label)
}))
joint_ifng_expression <- rbindlist(lapply(c(primary_label, winsor_label), function(label) {
  z <- model_data[[label]]
  fit_joint_competition(z$data, "IFNG_expression_z", candidate_catalog, z$extras, label)
}))
joint_context <- rbindlist(list(
  fit_joint_competition(primary_data, "IFNG_response_z", candidate_catalog, "Pan_T", "Primary + pan-T sensitivity"),
  fit_joint_competition(primary_data, "IFNG_response_z", candidate_catalog, c("Immune", "Stromal"), "Primary + broad TME sensitivity")
))

bootstrap_B <- 5000L
joint_boot <- rbindlist(list(
  joint_bootstrap(primary_data, "IFNG_response_z", candidate_catalog, character(), primary_label, bootstrap_B, 261300L),
  joint_bootstrap(model_data[[winsor_label]]$data, "IFNG_response_z", candidate_catalog, character(), winsor_label, bootstrap_B, 261301L)
))

serial_boot <- rbindlist(lapply(seq_len(nrow(candidate_catalog)), function(i) {
  route <- candidate_catalog$route[i]; mediator <- candidate_catalog$column[i]
  rbindlist(list(
    serial_bootstrap(primary_data, mediator, "IFNG_response_z", route, primary_label,
                     B = bootstrap_B, seed = 261400L + i, reverse = FALSE),
    serial_bootstrap(primary_data, mediator, "IFNG_response_z", route, primary_label,
                     B = bootstrap_B, seed = 261500L + i, reverse = TRUE),
    serial_bootstrap(model_data[[winsor_label]]$data, mediator, "IFNG_response_z", route, winsor_label,
                     B = bootstrap_B, seed = 261600L + i, reverse = FALSE)
  ))
}))

joint_rhs <- c("low_exposure_z", candidate_catalog$column, "batch", "sex", "age", "grade")
vif_audit <- rbindlist(lapply(seq_len(nrow(candidate_catalog)), function(i) {
  term <- candidate_catalog$column[i]
  data.table(route = candidate_catalog$route[i], predictor = term,
             VIF = vif_one(term, setdiff(joint_rhs, term), primary_data))
}))

cv_repeats <- 100L
cv_results <- repeated_cv(primary_data, "IFNG_response_z", candidate_catalog, cv_repeats, 5L, 261700L)
cv_summary <- cv_results[model != "Baseline", .(
  median_delta_R2 = median(delta_R2_vs_baseline),
  mean_delta_R2 = mean(delta_R2_vs_baseline),
  delta_R2_CI_low = unname(quantile(delta_R2_vs_baseline, 0.025)),
  delta_R2_CI_high = unname(quantile(delta_R2_vs_baseline, 0.975)),
  median_delta_RMSE = median(delta_RMSE_vs_baseline),
  delta_RMSE_CI_low = unname(quantile(delta_RMSE_vs_baseline, 0.025)),
  delta_RMSE_CI_high = unname(quantile(delta_RMSE_vs_baseline, 0.975))
), by = .(model, route)]

stability <- elastic_net_stability(primary_data, "IFNG_response_z", candidate_catalog,
                                   repetitions = 1000L, fraction = 0.80, seed = 261800L)

# Strict secondary TF diagnostics. Only TFs retaining >=10 A/B targets after
# removing all candidate and IFNG-outcome genes are estimable.
regulon <- fread(regulon_path, check.names = FALSE)
focal_tfs <- c("NFKB1", "RELA", "STAT1", "STAT2", "STAT3")
tf_excluded <- unique(c(unlist(raw_candidates, use.names = FALSE), response_core, "IFNG"))
strict_net <- regulon[
  TF %in% focal_tfs & expressed == TRUE & confidence %in% c("A", "B") & !target %in% tf_excluded
]
tf_coverage <- data.table(TF = focal_tfs)
tf_coverage <- merge(
  tf_coverage,
  strict_net[, .(retained_targets_n = uniqueN(target),
                 retained_targets = paste(sort(unique(target)), collapse = ";")), by = TF],
  by = "TF", all.x = TRUE
)
tf_coverage[is.na(retained_targets_n), `:=`(retained_targets_n = 0L, retained_targets = "")]
estimable_tfs <- tf_coverage[retained_targets_n >= 10L, TF]
tf_results <- data.table()
if (length(estimable_tfs)) {
  tf_activity <- ulm_from_matrix(log_expression, strict_net[TF %in% estimable_tfs], minsize = 10L)
  tf_activity <- tf_activity[estimable_tfs, sample_ids, drop = FALSE]
  tf_activity_z <- t(apply(tf_activity, 1L, zscore))
  rownames(tf_activity_z) <- estimable_tfs
  colnames(tf_activity_z) <- sample_ids
  tf_effect <- fit_score_effects(primary_data, tf_activity_z, character(), primary_label)
  setnames(tf_effect, "route", "TF")
  tf_map <- data.table(
    TF = c("NFKB1", "RELA", "STAT1", "STAT2", "STAT3"),
    candidate_axis = c("TNF/NF-kappaB", "TNF/NF-kappaB", "Type-II-IFN response context",
                       "Type-I-IFN/JAK-STAT", "IL6/JAK-STAT3")
  )
  tf_results <- merge(tf_effect, tf_map, by = "TF", all.x = TRUE)
}

# -------------------------------------------------------------------------
# Predeclared internal-consistency ranking
# -------------------------------------------------------------------------
get_gw <- function(label, name) {
  genomewide_anchors[model == label, .(
    route,
    value = normalized_low_NES,
    q = genomewide_BH_q
  )][, c("value", "q") := .(value, q)][]
}

ranking <- copy(candidate_catalog[, .(route, route_order, anchor_pathway, role)])
merge_metric <- function(base, dt, value_name, q_name) {
  x <- copy(dt)
  setnames(x, c("value", "q"), c(value_name, q_name))
  merge(base, x, by = "route", all.x = TRUE)
}
ranking <- merge_metric(ranking, get_gw("Median Low vs High"), "median_NES_low", "median_genomewide_q")
ranking <- merge_metric(ranking, get_gw("Q1 vs Q4"), "q4q1_NES_low", "q4q1_genomewide_q")
ranking <- merge_metric(ranking, get_gw("Q1 vs Q4, pan-T sensitivity"), "panT_NES_low", "panT_genomewide_q")
ranking <- merge_metric(ranking, get_gw("Q1 vs Q4, broad TME sensitivity"), "TME_NES_low", "TME_genomewide_q")
ranking <- merge_metric(ranking, get_gw("Q1 vs Q4, exclude 5"), "q4q1_exclude5_NES_low", "q4q1_exclude5_genomewide_q")

fg_ex <- focused_gsea[model == primary_label, .(route, continuous_exclude5_NES = NES,
                                                continuous_exclude5_focused_q = padj)]
fg_win <- focused_gsea[model == winsor_label, .(route, continuous_winsor_NES = NES,
                                                continuous_winsor_focused_q = padj)]
se_ex <- score_effects[model == primary_label, .(route, score_beta = standardized_beta, score_BH_q = BH_q)]
jr <- joint_results[model == primary_label, .(route, joint_beta = standardized_beta,
                                              joint_BH_q = BH_q, joint_partial_R2 = partial_R2)]
jb <- joint_boot[model == primary_label, .(route, joint_boot_CI_low = bootstrap_CI_low,
                                           joint_boot_CI_high = bootstrap_CI_high,
                                           joint_boot_positive_fraction = positive_fraction,
                                           joint_boot_CI_excludes_zero = CI_excludes_zero)]
sf <- serial_boot[model == primary_label & order == "forward", .(
  route, forward_product = serial_product, forward_CI_low = bootstrap_CI_low,
  forward_CI_high = bootstrap_CI_high, forward_CI_excludes_zero = CI_excludes_zero
)]
sr <- serial_boot[model == primary_label & order == "reverse", .(
  route, reverse_product = serial_product, reverse_CI_low = bootstrap_CI_low,
  reverse_CI_high = bootstrap_CI_high, reverse_CI_excludes_zero = CI_excludes_zero
)]
ranking <- Reduce(function(x, y) merge(x, y, by = "route", all.x = TRUE),
                  list(ranking, fg_ex, fg_win, se_ex, jr, jb, sf, sr, stability,
                       cv_summary[route %in% candidate_catalog$route], vif_audit))

ranking[, `:=`(
  criterion_median_GSEA = median_NES_low > 0 & median_genomewide_q < 0.05,
  criterion_Q4Q1_GSEA = q4q1_NES_low > 0 & q4q1_genomewide_q < 0.05,
  criterion_Q4Q1_exclude5_GSEA = q4q1_exclude5_NES_low > 0 & q4q1_exclude5_genomewide_q < 0.05,
  criterion_panT_GSEA = panT_NES_low > 0 & panT_genomewide_q < 0.05,
  criterion_broad_TME_GSEA = TME_NES_low > 0 & TME_genomewide_q < 0.05,
  criterion_continuous_exclude5_GSEA = continuous_exclude5_NES > 0 & continuous_exclude5_focused_q < 0.05,
  criterion_continuous_winsor_GSEA = continuous_winsor_NES > 0 & continuous_winsor_focused_q < 0.05,
  criterion_sample_score = score_beta > 0 & score_BH_q < 0.05,
  criterion_joint_IFNG = joint_beta > 0 & joint_BH_q < 0.05,
  criterion_joint_bootstrap = joint_boot_CI_low > 0 & joint_boot_positive_fraction >= 0.80,
  criterion_forward_serial = forward_product > 0 & forward_CI_excludes_zero,
  criterion_elastic_stability = lambda_1se_positive_frequency >= 0.70,
  criterion_CV = delta_R2_CI_low > 0,
  reverse_order_also_supported = reverse_CI_excludes_zero
)]
criterion_cols <- grep("^criterion_", names(ranking), value = TRUE)
exposure_criterion_cols <- c(
  "criterion_median_GSEA", "criterion_Q4Q1_GSEA", "criterion_Q4Q1_exclude5_GSEA",
  "criterion_panT_GSEA", "criterion_broad_TME_GSEA",
  "criterion_continuous_exclude5_GSEA", "criterion_continuous_winsor_GSEA",
  "criterion_sample_score"
)
ifng_criterion_cols <- c(
  "criterion_joint_IFNG", "criterion_joint_bootstrap", "criterion_forward_serial",
  "criterion_elastic_stability", "criterion_CV"
)
ranking[, criteria_passed := rowSums(.SD, na.rm = TRUE), .SDcols = criterion_cols]
ranking[, criteria_total := length(criterion_cols)]
ranking[, exposure_criteria_passed := rowSums(.SD, na.rm = TRUE), .SDcols = exposure_criterion_cols]
ranking[, exposure_criteria_total := length(exposure_criterion_cols)]
ranking[, IFNG_link_criteria_passed := rowSums(.SD, na.rm = TRUE), .SDcols = ifng_criterion_cols]
ranking[, IFNG_link_criteria_total := length(ifng_criterion_cols)]
ranking[, complete_association_route := criteria_passed == criteria_total]
# The primary ranking answers which program is most consistently associated
# with measured Sarcosine. IFNG-link evidence is a secondary tie-breaker; it
# must not overwhelm the exposure-to-program evidence.
setorder(ranking, -exposure_criteria_passed, -IFNG_link_criteria_passed,
         -criteria_passed, joint_BH_q, route_order)
ranking[, evidence_rank := seq_len(.N)]
leading_route <- ranking[1, route]
lead <- ranking[1]
leading_ifng_route <- ranking[order(-IFNG_link_criteria_passed, -exposure_criteria_passed,
                                    -criteria_passed, joint_BH_q, route_order), route][1]
lead_ifng <- ranking[route == leading_ifng_route]

criteria_long <- melt(
  ranking[, c("route", criterion_cols), with = FALSE], id.vars = "route",
  variable.name = "criterion", value.name = "passed"
)
criterion_labels <- c(
  criterion_median_GSEA = "Median GSEA",
  criterion_Q4Q1_GSEA = "Q1/Q4 GSEA",
  criterion_Q4Q1_exclude5_GSEA = "Q1/Q4 exclude-5",
  criterion_panT_GSEA = "After pan-T",
  criterion_broad_TME_GSEA = "Broad TME sensitivity",
  criterion_continuous_exclude5_GSEA = "Continuous exclude-5",
  criterion_continuous_winsor_GSEA = "Continuous winsorized",
  criterion_sample_score = "Sample-score association",
  criterion_joint_IFNG = "Joint IFNG competition",
  criterion_joint_bootstrap = "Joint bootstrap",
  criterion_forward_serial = "Forward serial CI",
  criterion_elastic_stability = "Elastic-net stability",
  criterion_CV = "Cross-validated gain"
)
criteria_long[, criterion_label := factor(criterion_labels[criterion], levels = rev(unname(criterion_labels)))]
criteria_long[, route := factor(route, levels = rev(ranking$route))]

# -------------------------------------------------------------------------
# Reproducibility outputs
# -------------------------------------------------------------------------
input_manifest <- data.table(
  input_role = c(
    "TPM expression", "matched metadata", "quartile map", "cell-context estimates",
    "signed DoRothEA A/B regulon", "genome-wide median GSEA", "genome-wide Q4/Q1 GSEA",
    "genome-wide Q4/Q1 exclude-five GSEA"
  ),
  path = required_inputs,
  sha256 = vapply(required_inputs, sha256_file, character(1))
)
method_contract <- data.table(
  item = c(
    "Population", "Exposure", "Primary inference", "Reader contrasts", "Covariates",
    "Candidate status", "GSEA multiplicity", "Program scores", "Outcome",
    "Joint competition", "Resampling", "TF activity", "Interpretation"
  ),
  specification = c(
    "100 matched treatment-naive ccRCC tumours with GC-MS and bulk RNA-seq",
    "Normalized untargeted tissue Sarcosine; not absolute concentration or flux",
    "Continuous z(-log2 Sarcosine) after excluding five prespecified Tukey-low influential samples; winsorized sensitivity",
    "Median Low/High 50/50 and exploratory Q1/Q4 25/25",
    "RNA batch, sex, age category and tumour grade; pan-T and broad TME sensitivities retained in tables",
    "Six pathways selected after prior ccRCC GSEA; exploratory/post-hoc ranking",
    "Genome-wide BH q imported from validated full-collection runs; new focused BH q across six candidates is labelled separately",
    "Mean row-z expression; candidate programs mutually disjoint and disjoint from the IFNG-response outcome",
    "Cross-collection IFNG-response core present in at least two of Hallmark/Reactome/GO:BP; IFNG gene excluded",
    "All six programs entered together with exposure and clinical covariates",
    "5,000 bootstrap replicates; 100 repeated nested five-fold ridge CV; 1,000 elastic-net subsamples; fixed seeds",
    "Signed DoRothEA A/B ULM only when >=10 targets remain after removing all candidate/outcome genes; secondary diagnostic",
    "Internal consistency only; no temporal, causal, same-cell, protein or cytokine-secretion inference"
  )
)
package_versions <- data.table(
  package = c("R", required_packages),
  version = c(R.version.string, vapply(required_packages, function(x) as.character(packageVersion(x)), character(1)))
)
sample_scores <- analysis_data[, c(
  "sample_id", "sarcosine_normalized_intensity", "log2_sarcosine_intensity",
  "sarcosine_group", candidate_catalog$column, "IFNG_response_z", "IFNG_expression_z",
  "Pan_T", "Immune", "Stromal"
), with = FALSE]

validation <- data.table(
  check = c(
    "100 matched tumours", "median groups 50/50", "quartile groups 25/25",
    "five prespecified Tukey-low samples", "15119-gene universe",
    "six candidate pathways", "strict program scores >=10 genes",
    "zero candidate/candidate gene overlap", "zero candidate/outcome gene overlap",
    "all model designs full rank", "all focused GSEA rows returned",
    ">=95% bootstrap replicates valid", "100 CV repeats complete",
    ">=95% elastic-net repetitions valid", "finite VIF values", "unique MSigDB version"
  ),
  passed = c(
    nrow(meta) == 100L && uniqueN(meta$sample_id) == 100L,
    sum(meta$sarcosine_group == "Low") == 50L && sum(meta$sarcosine_group == "High") == 50L,
    sum(model_data[["Q1 vs Q4"]]$data$sarcosine_group == "Low") == 25L &&
      sum(model_data[["Q1 vs Q4"]]$data$sarcosine_group == "High") == 25L,
    identical(outlier_ids, expected_outliers), nrow(log_expression) == 15119L,
    nrow(candidate_catalog) == 6L, all(score_manifest$scored_genes_n >= 10L),
    score_overlap[score_1 != score_2, max(overlap_n)] == 0L,
    max(vapply(strict_candidates, function(g) length(intersect(g, response_core)), integer(1))) == 0L,
    all(vapply(rank_fits, function(x) qr(x$design)$rank == ncol(x$design), logical(1))),
    nrow(focused_gsea) == length(model_data) * nrow(candidate_catalog),
    all(joint_boot$bootstrap_replicates_valid >= 0.95 * bootstrap_B) &&
      all(serial_boot$bootstrap_replicates_valid >= 0.95 * bootstrap_B),
    uniqueN(cv_results$repeat_id) == cv_repeats,
    all(stability$repetitions_valid >= 950L), all(is.finite(vif_audit$VIF)),
    length(msigdb_version) == 1L
  )
)
if (!all(validation$passed)) stop("Validation failure: ", paste(validation[passed == FALSE, check], collapse = "; "))

fwrite(input_manifest, file.path(table_dir, "00_input_manifest.csv"))
fwrite(method_contract, file.path(table_dir, "00_method_contract.csv"))
fwrite(package_versions, file.path(table_dir, "00_package_versions.csv"))
fwrite(model_manifest, file.path(table_dir, "01_model_manifest.csv"))
fwrite(candidate_catalog, file.path(table_dir, "02_candidate_program_catalog.csv"))
fwrite(score_manifest, file.path(table_dir, "03_strict_score_gene_manifest.csv"))
fwrite(score_overlap, file.path(table_dir, "04_strict_score_overlap_audit.csv"))
fwrite(sample_scores, file.path(table_dir, "05_sample_level_program_scores.csv"))
fwrite(focused_gsea_out, file.path(table_dir, "06_focused_candidate_GSEA_all_models.csv"))
fwrite(genomewide_anchors, file.path(table_dir, "07_genomewide_q_anchor_results.csv"))
fwrite(score_effects, file.path(table_dir, "08_sample_program_effects_all_models.csv"))
fwrite(joint_results, file.path(table_dir, "09_joint_program_competition_IFNG_response.csv"))
fwrite(joint_ifng_expression, file.path(table_dir, "10_joint_program_competition_IFNG_expression.csv"))
fwrite(joint_context, file.path(table_dir, "11_joint_program_cell_context_sensitivities.csv"))
fwrite(joint_boot, file.path(table_dir, "12_joint_program_bootstrap.csv"))
fwrite(serial_boot, file.path(table_dir, "13_program_to_IFNG_serial_bootstrap.csv"))
fwrite(cv_results, file.path(table_dir, "14_repeated_CV_all_repeats.csv"))
fwrite(cv_summary, file.path(table_dir, "15_repeated_CV_summary.csv"))
fwrite(stability, file.path(table_dir, "16_elastic_net_stability.csv"))
fwrite(tf_coverage, file.path(table_dir, "17_strict_TF_coverage.csv"))
fwrite(tf_results, file.path(table_dir, "18_strict_TF_diagnostics.csv"))
fwrite(vif_audit, file.path(table_dir, "19_VIF_audit.csv"))
fwrite(ranking, file.path(table_dir, "20_final_evidence_ranking.csv"))
fwrite(criteria_long, file.path(table_dir, "21_evidence_criteria_long.csv"))
fwrite(validation, file.path(table_dir, "22_validation_checks.csv"))

# -------------------------------------------------------------------------
# Reader-facing figures (no ESTIMATE-labelled panels)
# -------------------------------------------------------------------------
low_color <- "#1F9E93"; high_color <- "#C93C3C"; neutral <- "#9A9A9A"; ink <- "#202020"
theme_reader <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 10.5, color = "#555555"),
    axis.title = element_text(face = "bold"), axis.text.y = element_text(color = ink),
    legend.position = "bottom", strip.text = element_text(face = "bold", size = 11)
  )

display_models <- model_manifest[reader_facing == TRUE, model]
gsea_plot <- focused_gsea[model %in% display_models]
gsea_plot[, route := factor(route, levels = rev(ranking$route))]
gsea_plot[, model := factor(model, levels = display_models)]
gsea_plot[, significant := is.finite(padj) & padj < 0.05]
p_gsea <- ggplot(gsea_plot, aes(NES, route)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#777777") +
  geom_point(aes(color = NES > 0, shape = significant, size = -log10(pmax(padj, .Machine$double.xmin))), stroke = 0.7) +
  facet_wrap(~ model, ncol = 1) +
  scale_color_manual(values = c(`TRUE` = low_color, `FALSE` = high_color), guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "Focused BH q < 0.05") +
  scale_size_continuous(name = expression(-log[10](focused~BH~q)), range = c(2.5, 6.5)) +
  labs(
    title = "a  Candidate-pathway enrichment",
    subtitle = "Positive NES denotes Sarcosine-Low; six-candidate focused sensitivity",
    x = "Normalized enrichment score", y = NULL
  ) + theme_reader

score_plot <- score_effects[model %in% display_models]
score_plot[, route := factor(route, levels = rev(ranking$route))]
score_plot[, model := factor(model, levels = display_models)]
p_score <- ggplot(score_plot, aes(standardized_beta, route)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#777777") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0, color = ink) +
  geom_point(aes(color = standardized_beta > 0, shape = BH_q < 0.05), size = 3.3, stroke = 0.7) +
  facet_wrap(~ model, ncol = 1) +
  scale_color_manual(values = c(`TRUE` = low_color, `FALSE` = high_color), guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "BH q < 0.05") +
  labs(
    title = "b  Gene-disjoint program-score associations",
    subtitle = "Positive effects denote Sarcosine-Low; clinical covariates adjusted",
    x = "Low-associated effect (SD)", y = NULL
  ) + theme_reader

joint_plot <- joint_results[model == primary_label]
joint_plot[, route := factor(route, levels = rev(ranking$route))]
p_joint <- ggplot(joint_plot, aes(standardized_beta, route)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#777777") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0, color = ink) +
  geom_point(aes(color = standardized_beta > 0, shape = BH_q < 0.05), size = 3.8, stroke = 0.8) +
  scale_color_manual(values = c(`TRUE` = low_color, `FALSE` = high_color), guide = "none") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "BH q < 0.05") +
  labs(
    title = "c  Joint competition for IFNG-response association",
    subtitle = "All six non-overlapping programs entered together; five influential low values excluded",
    x = "Adjusted association with IFNG-response score (SD)", y = NULL
  ) + theme_reader

p_criteria <- ggplot(criteria_long, aes(passed, criterion_label)) +
  geom_tile(aes(fill = passed), color = "white", linewidth = 0.5) +
  facet_wrap(~ route, nrow = 1) +
  scale_x_discrete(position = "top") +
  scale_fill_manual(values = c(`TRUE` = low_color, `FALSE` = "#E8E8E8"), guide = "none") +
  labs(
    title = "d  Internal-consistency criteria",
    subtitle = "Green = criterion passed; this is not a causal score",
    x = NULL, y = NULL
  ) + theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_blank(), panel.grid = element_blank(),
    strip.text = element_text(face = "bold", size = 9), plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 10.5, color = "#555555")
  )

summary_figure <- (p_gsea | p_score) / (p_joint | p_criteria) +
  plot_annotation(
    title = "Measured tumour Sarcosine: final immune-program competition in ccRCC",
    subtitle = "Positive estimates denote lower measured Sarcosine. Candidate selection is exploratory; associations do not establish signalling direction.",
    theme = theme(plot.title = element_text(face = "bold", size = 21),
                  plot.subtitle = element_text(size = 12, color = "#555555"))
  )
ggsave(file.path(figure_dir, "Fig_Measured_Sarcosine_Immune_Pathway_Final_Competition.png"),
       summary_figure, width = 20, height = 17, dpi = 320, bg = "white")
ggsave(file.path(figure_dir, "Fig_Measured_Sarcosine_Immune_Pathway_Final_Competition.pdf"),
       summary_figure, width = 20, height = 17, device = cairo_pdf)

# Descriptive heatmap, ordered by measured Sarcosine without clustering columns.
heat_mat <- rbind(program_matrix, `IFNG-response program` = analysis_data$IFNG_response_z)
ord <- order(meta$sarcosine_group, meta$sarcosine_normalized_intensity)
heat_mat <- heat_mat[, ord, drop = FALSE]
ann <- data.frame(`Measured Sarcosine` = meta$sarcosine_group[ord], row.names = sample_ids[ord], check.names = FALSE)
ann_colors <- list(`Measured Sarcosine` = c(Low = low_color, High = high_color))
png(file.path(figure_dir, "Fig_Measured_Sarcosine_Immune_Program_Heatmap.png"),
    width = 5000, height = 2400, res = 320, bg = "white")
pheatmap(
  heat_mat, cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = FALSE,
  annotation_col = ann, annotation_colors = ann_colors,
  color = colorRampPalette(c("#2C7FB8", "white", "#D73027"))(100),
  breaks = seq(-2.5, 2.5, length.out = 101), border_color = NA,
  fontsize = 12, fontsize_row = 12,
  main = "Candidate immune-program activities across matched ccRCC tumours"
)
dev.off()

# Leading-program schematic is deliberately program-level. It does not insert
# an unsupported receptor or TF into the route.
schematic_file <- file.path(figure_dir, "Fig_Leading_Measured_Sarcosine_Immune_Candidate_Schematic.png")
png(schematic_file, width = 5600, height = 2200, res = 320, bg = "white")
par(mar = c(0, 0, 0, 0), xpd = NA)
plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
text(0.03, 0.93, "Leading measured-Sarcosine-associated immune program", adj = 0,
     cex = 2.2, font = 2)
subtitle <- sprintf(
  "%s was most robustly associated with measured Sarcosine (%d/%d exposure criteria).",
  leading_route, lead$exposure_criteria_passed, lead$exposure_criteria_total
)
text(0.03, 0.875, subtitle, adj = 0, cex = 1.25, col = "#444444")
box_node <- function(x, y, label, width) {
  rect(x - width/2, y - 0.075, x + width/2, y + 0.075, lwd = 2.2, border = "black")
  text(x, y, label, cex = 1.55, font = 2)
}
box_node(0.20, 0.64, "Lower measured\ntumour Sarcosine", 0.21)
box_node(0.55, 0.64, leading_route, 0.25)
arrows(0.31, 0.64, 0.415, 0.64, lwd = 2.8, lty = 2, col = low_color, length = 0.12)
text(0.375, 0.69, "strongest exposure-associated program", cex = 0.95, col = "#555555")
box_node(0.38, 0.38, leading_ifng_route, 0.25)
box_node(0.76, 0.38, "IFNG-response\nprogram", 0.20)
arrows(0.51, 0.38, 0.655, 0.38, lwd = 2.8, lty = 2, col = "#777777", length = 0.12)
text(0.585, 0.43, "strongest independent IFNG coupling", cex = 0.95, col = "#555555")
text(0.50, 0.235,
     "No candidate formed a complete lower-Sarcosine -> program -> IFNG association chain.",
     cex = 1.15, font = 2, col = "#555555")
text(0.50, 0.165,
     "The two displayed evidence layers must not be joined as a causal pathway without perturbation and cell-resolved validation.",
     cex = 1.02, col = "#555555")
dev.off()

validation[, figures_exist := TRUE]
if (!all(file.exists(c(
  file.path(figure_dir, "Fig_Measured_Sarcosine_Immune_Pathway_Final_Competition.png"),
  file.path(figure_dir, "Fig_Measured_Sarcosine_Immune_Pathway_Final_Competition.pdf"),
  file.path(figure_dir, "Fig_Measured_Sarcosine_Immune_Program_Heatmap.png"), schematic_file
)))) stop("A reader-facing figure is missing")

report_lines <- c(
  "# Final immune-program narrowing for measured tumour Sarcosine in ccRCC", "",
  "## Bottom line", "",
  sprintf("The program most consistently associated with lower measured Sarcosine was **%s**, passing **%d/%d exposure-association criteria**.",
          leading_route, lead$exposure_criteria_passed, lead$exposure_criteria_total),
  sprintf("The strongest independent correlate of the gene-disjoint IFNG-response outcome was **%s** (%d/%d IFNG-link criteria).",
          leading_ifng_route, lead_ifng$IFNG_link_criteria_passed, lead_ifng$IFNG_link_criteria_total),
  "No candidate passed the forward serial-bootstrap criterion, and no candidate passed every criterion. These two evidence layers therefore cannot be joined into a complete measured-Sarcosine-to-IFNG mechanism.", "",
  "## Design", "",
  "Continuous lower Sarcosine after exclusion of five prespecified Tukey-low influential values was the primary inference model. Winsorized continuous, median 50/50 and Q1/Q4 contrasts were sensitivity analyses. Six post-hoc immune candidates were scored with mutually exclusive genes and were also disjoint from a cross-collection IFNG-response outcome. Genome-wide GSEA q values from previously validated runs were kept distinct from the new six-candidate focused q values.", "",
  "## Final ranking", "",
  "| Rank | Candidate | Exposure criteria | IFNG-link criteria | All criteria | Median genome-wide q | Q1/Q4 genome-wide q | Broad-TME q | Joint IFNG q | Forward serial 95% CI |",
  "|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(ranking)), function(i) {
    x <- ranking[i]
    sprintf(
      "| %d | %s | %d/%d | %d/%d | %d/%d | %.3g | %.3g | %.3g | %.3g | [%.4f, %.4f] |",
      x$evidence_rank, x$route,
      x$exposure_criteria_passed, x$exposure_criteria_total,
      x$IFNG_link_criteria_passed, x$IFNG_link_criteria_total,
      x$criteria_passed, x$criteria_total,
      x$median_genomewide_q, x$q4q1_genomewide_q, x$TME_genomewide_q,
      x$joint_BH_q, x$forward_CI_low, x$forward_CI_high
    )
  }, character(1)), "",
  "## TF boundary", "",
  paste0("Strictly estimable TFs after removing all candidate/outcome genes: ",
         ifelse(length(estimable_tfs), paste(estimable_tfs, collapse = ", "), "none"), "."),
  "These inferred activities are secondary diagnostics and do not measure phosphorylation, nuclear localization or chromatin binding. A pathway is not labelled canonical NF-kappaB solely from Hallmark enrichment unless NFKB1/RELA diagnostics independently support that designation.", "",
  "## Interpretation boundary", "",
  "All results are observational and post-hoc. They can prioritize a transcriptome program associated with lower normalized tissue Sarcosine, but cannot show that Sarcosine changes that program, that the events occur in the same cell, or that IFNG protein is produced. A mechanistic conclusion requires perturbation, protein/phosphoprotein measurements and preferably single-cell or spatial validation.", "",
  "## Reader-facing outputs", "",
  "- `figures/Fig_Measured_Sarcosine_Immune_Pathway_Final_Competition.png`",
  "- `figures/Fig_Measured_Sarcosine_Immune_Program_Heatmap.png`",
  "- `figures/Fig_Leading_Measured_Sarcosine_Immune_Candidate_Schematic.png`"
)
writeLines(report_lines, file.path(output_root, "REPORT.md"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))

cat("Completed final measured-Sarcosine immune-program narrowing.\n")
cat("Output:", output_root, "\n")
cat("Leading Sarcosine-associated program:", leading_route,
    lead$exposure_criteria_passed, "/", lead$exposure_criteria_total, "exposure criteria\n")
cat("Leading IFNG-linked program:", leading_ifng_route,
    lead_ifng$IFNG_link_criteria_passed, "/", lead_ifng$IFNG_link_criteria_total, "IFNG-link criteria\n")
cat("Complete association route:", lead$complete_association_route, "\n")
cat("Validation:", sum(validation$passed), "/", nrow(validation), "passed\n")
