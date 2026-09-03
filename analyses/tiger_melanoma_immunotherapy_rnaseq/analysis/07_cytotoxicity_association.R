#!/usr/bin/env Rscript
# =============================================================================
# Candidate Fig. 2d: host sarcosine-metabolism scores versus tumour cytolytic
# activity (CYT) in TIGER melanoma cohort PRJEB23709.
#
# Prespecified design
#   Biological question: Is host sarcosine-degradation transcriptional capacity
#     positively associated with a cytolytic tumour-immune expression state?
#   Primary unit: one patient / one pretreatment biopsy (PRE; expected n = 73).
#   Primary test: two-sided Spearman correlation, degradation score versus CYT.
#   Specificity control: production score versus CYT; Holm correction across the
#     two score-level tests.
#   Adjusted sensitivity: HC3 linear model controlling response, therapy, age,
#     and sex. Response is not treated as a confounder in the primary analysis;
#     adjustment asks whether the association is more than R/NR group separation.
#   Repeated-sample sensitivity: all 91 biopsies, adjusted for timepoint with a
#     patient random intercept.
#
# Data provenance
#   TIGER processed bulk RNA-seq, FPKM values. The source manuscript records
#   FastQC/Cutadapt/STAR-GRCh38/featureCounts as the upstream workflow.
#   This is a continuous-expression association analysis; count-based DE tools
#   such as DESeq2/edgeR are not appropriate for these FPKM inputs.
# =============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(here)
  library(lmerTest)
  library(lmtest)
  library(sandwich)
})

# Locate the project from this script, then register a portable here() root.
file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) != 1L) stop("Run this analysis with Rscript.")
script_file <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
setwd(project_dir)
here::i_am(file.path("analysis", basename(script_file)))
project_dir <- here::here()

expr_file <- here::here("Melanoma-PRJEB23709_ExpressionData.tsv")
clin_file <- here::here("Melanoma-PRJEB23709_ClinicalData.tsv")
saved_rds <- here::here("results", "analysis_data_all91.rds")
out_dir <- here::here(
  "results", "manuscript_figures", "Fig2d_CYT_candidate_26.08.23"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(expr_file, clin_file, saved_rds)
if (!all(file.exists(required_files))) {
  stop("Missing required input(s): ",
       paste(required_files[!file.exists(required_files)], collapse = ", "))
}

log_file <- file.path(out_dir, "analysis_log.txt")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0L) sink(type = "message")
  while (sink.number(type = "output") > 0L) sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("Candidate Fig. 2d CYT analysis\n")
cat("Run time:", format(Sys.time(), tz = "Asia/Seoul"), "KST\n\n")

# =============================================================================
# 1. Read and verify the off-by-one TIGER TSV layouts.
# =============================================================================
expr_header <- strsplit(readLines(expr_file, n = 1L), "\t")[[1]]
stopifnot(expr_header[1] == "GENE_SYMBOL")
sample_ids <- expr_header[-1]

expr_raw <- data.table::fread(
  expr_file, header = FALSE, sep = "\t", skip = 1L,
  showProgress = FALSE, data.table = FALSE
)
stopifnot(ncol(expr_raw) == length(sample_ids) + 2L)
gene_symbols <- as.character(expr_raw[[2]])
if (anyDuplicated(gene_symbols)) stop("Duplicated gene symbols in expression TSV.")

required_genes <- c("SARDH", "PIPOX", "GNMT", "DMGDH", "GZMA", "PRF1")
gene_rows <- match(required_genes, gene_symbols)
if (anyNA(gene_rows)) {
  stop("Required gene(s) absent: ",
       paste(required_genes[is.na(gene_rows)], collapse = ", "))
}
expr_all_fpkm <- as.matrix(expr_raw[, 3:ncol(expr_raw), drop = FALSE])
storage.mode(expr_all_fpkm) <- "double"
colnames(expr_all_fpkm) <- sample_ids
stopifnot(
  !anyNA(expr_all_fpkm), all(is.finite(expr_all_fpkm)),
  all(expr_all_fpkm >= 0)
)
sample_fpkm_totals <- colSums(expr_all_fpkm)
stopifnot(all(is.finite(sample_fpkm_totals)), all(sample_fpkm_totals > 0))

expr_fpkm <- expr_all_fpkm[gene_rows, , drop = FALSE]
rownames(expr_fpkm) <- required_genes
stopifnot(!anyNA(expr_fpkm), all(is.finite(expr_fpkm)), all(expr_fpkm >= 0))

# Exact within-sample conversion from FPKM to TPM using every feature in the
# deposited TIGER expression matrix. This permits the original TPM-based CYT
# definition while retaining the existing Figure 2 FPKM-based sarcosine scores.
expr_tpm <- sweep(expr_fpkm, 2, sample_fpkm_totals, "/") * 1e6
stopifnot(
  all(is.finite(expr_tpm)), all(expr_tpm >= 0),
  all(abs(sample_fpkm_totals / sample_fpkm_totals * 1e6 - 1e6) < 1e-8)
)

clin_header <- strsplit(readLines(clin_file, n = 1L), "\t")[[1]]
clin_raw <- data.table::fread(
  clin_file, header = FALSE, sep = "\t", skip = 1L,
  na.strings = c("NA", ""), showProgress = FALSE, data.table = FALSE
)
stopifnot(ncol(clin_raw) == length(clin_header) + 1L)
colnames(clin_raw) <- c("row_index", clin_header)

clinical <- data.frame(
  sample_id = as.character(clin_raw[["sample_id"]]),
  patient_name = as.character(clin_raw[["patient_name"]]),
  timepoint = as.character(clin_raw[["Treatment"]]),
  recist = as.character(clin_raw[["response"]]),
  response_code = as.character(clin_raw[["response_NR"]]),
  gender = as.character(clin_raw[["Gender"]]),
  therapy = as.character(clin_raw[["Therapy"]]),
  age = as.numeric(clin_raw[["age_start"]]),
  stringsAsFactors = FALSE
)

stopifnot(
  nrow(clinical) == 91L,
  length(sample_ids) == 91L,
  !anyDuplicated(sample_ids),
  !anyDuplicated(clinical$sample_id),
  setequal(sample_ids, clinical$sample_id),
  !anyNA(clinical),
  all(clinical$timepoint %in% c("PRE", "EDT")),
  all(clinical$response_code %in% c("R", "N")),
  all(clinical$gender %in% c("Female", "Male")),
  all(is.finite(clinical$age))
)

clinical <- clinical[match(sample_ids, clinical$sample_id), , drop = FALSE]
stopifnot(identical(clinical$sample_id, sample_ids))

dat <- clinical
for (gene in c("SARDH", "PIPOX", "GNMT", "DMGDH")) {
  dat[[gene]] <- as.numeric(log2(expr_fpkm[gene, ] + 1))
}
dat$GZMA_FPKM <- as.numeric(expr_fpkm["GZMA", ])
dat$PRF1_FPKM <- as.numeric(expr_fpkm["PRF1", ])
dat$GZMA_TPM <- as.numeric(expr_tpm["GZMA", ])
dat$PRF1_TPM <- as.numeric(expr_tpm["PRF1", ])
if (any(dat$GZMA_TPM <= 0) || any(dat$PRF1_TPM <= 0)) {
  stop("CYT genes include zero/non-positive TPM; no-pseudocount formula invalid.")
}
# Rooney-style log-average of the GZMA and PRF1 geometric mean in TPM.
dat$CYT_score <- 0.5 * (log2(dat$GZMA_TPM) + log2(dat$PRF1_TPM))
# Retained only as a unit-sensitivity check, not as the primary CYT definition.
dat$CYT_score_FPKM_sensitivity <- 0.5 * (
  log2(dat$GZMA_FPKM) + log2(dat$PRF1_FPKM)
)
dat$response_group <- factor(
  ifelse(dat$response_code == "R", "R", "NR"), levels = c("R", "NR")
)
dat$therapy_short <- factor(
  ifelse(
    dat$therapy == "anti-PD-1", "antiPD1",
    ifelse(dat$therapy == "anti-CTLA-4+anti-PD-1", "combo", NA_character_)
  ),
  levels = c("antiPD1", "combo")
)
dat$gender <- factor(dat$gender, levels = c("Female", "Male"))
dat$timepoint <- factor(dat$timepoint, levels = c("PRE", "EDT"))
stopifnot(!anyNA(dat$therapy_short), !anyNA(dat$response_group))

# Confirm agreement with the already validated RDS used for current Figure 2.
saved <- readRDS(saved_rds)
saved <- saved[match(dat$sample_id, saved$sample_id), , drop = FALSE]
stopifnot(
  identical(dat$sample_id, saved$sample_id),
  identical(as.character(dat$response_group), as.character(saved$response_group)),
  identical(as.character(dat$timepoint), as.character(saved$timepoint)),
  isTRUE(all.equal(dat$SARDH, saved$SARDH, tolerance = 1e-12)),
  isTRUE(all.equal(dat$PIPOX, saved$PIPOX, tolerance = 1e-12)),
  isTRUE(all.equal(dat$GNMT, saved$GNMT, tolerance = 1e-12)),
  isTRUE(all.equal(dat$DMGDH, saved$DMGDH, tolerance = 1e-12))
)

# =============================================================================
# 2. PRE-only score construction. PRE means/SDs are also applied to all 91
#    samples for the repeated-sample sensitivity, anchoring the scale to baseline.
# =============================================================================
pre_index <- dat$timepoint == "PRE"
pre <- dat[pre_index, , drop = FALSE]
stopifnot(
  nrow(pre) == 73L,
  length(unique(pre$patient_name)) == 73L,
  !anyDuplicated(pre$patient_name),
  sum(pre$response_group == "R") == 40L,
  sum(pre$response_group == "NR") == 33L
)

z_from_pre <- function(x_all, x_pre) {
  stopifnot(is.numeric(x_all), is.numeric(x_pre), !anyNA(x_all), sd(x_pre) > 0)
  (x_all - mean(x_pre)) / sd(x_pre)
}
for (gene in c("SARDH", "PIPOX", "GNMT", "DMGDH")) {
  dat[[paste0(gene, "_z_PRE")]] <- z_from_pre(dat[[gene]], pre[[gene]])
}
dat$Degradation_score <- rowMeans(dat[, c("SARDH_z_PRE", "PIPOX_z_PRE")])
dat$Production_score <- rowMeans(dat[, c("GNMT_z_PRE", "DMGDH_z_PRE")])
pre <- dat[pre_index, , drop = FALSE]

stopifnot(
  abs(mean(pre$SARDH_z_PRE)) < 1e-12,
  abs(mean(pre$PIPOX_z_PRE)) < 1e-12,
  abs(mean(pre$GNMT_z_PRE)) < 1e-12,
  abs(mean(pre$DMGDH_z_PRE)) < 1e-12,
  !anyNA(pre[, c("Degradation_score", "Production_score", "CYT_score")])
)

cat("Data QC\n")
cat("  Expression matrix:", nrow(expr_raw), "genes x", length(sample_ids), "samples\n")
cat("  Required-gene FPKM range:",
    paste(signif(range(expr_fpkm), 5), collapse = " to "), "\n")
cat("  Total FPKM per sample range:",
    paste(signif(range(sample_fpkm_totals), 6), collapse = " to "), "\n")
cat("  PRE patients:", nrow(pre), "(R=", sum(pre$response_group == "R"),
    ", NR=", sum(pre$response_group == "NR"), ")\n", sep = "")
cat("  All samples:", nrow(dat), "from", length(unique(dat$patient_name)), "patients\n")
cat("  CYT gene TPM zeros: GZMA=", sum(dat$GZMA_TPM == 0),
    ", PRF1=", sum(dat$PRF1_TPM == 0), "\n\n", sep = "")

qc <- data.frame(
  check = c(
    "Expression genes", "Expression samples", "Clinical rows", "PRE samples",
    "PRE unique patients", "All unique patients",
    "Missing values in analysis fields", "Minimum total FPKM",
    "Maximum total FPKM", "GZMA zero TPM", "PRF1 zero TPM"
  ),
  value = c(
    nrow(expr_raw), length(sample_ids), nrow(clinical), nrow(pre),
    length(unique(pre$patient_name)), length(unique(dat$patient_name)),
    sum(is.na(dat[, c(
      "sample_id", "patient_name", "timepoint", "response_group",
      "therapy_short", "gender", "age", "Degradation_score",
      "Production_score", "CYT_score"
    )])),
    min(sample_fpkm_totals), max(sample_fpkm_totals),
    sum(dat$GZMA_TPM == 0), sum(dat$PRF1_TPM == 0)
  ),
  stringsAsFactors = FALSE
)
write.csv(qc, file.path(out_dir, "data_QC.csv"), row.names = FALSE)

# =============================================================================
# 3. Primary score-level Spearman tests plus patient bootstrap CIs.
# =============================================================================
bootstrap_spearman <- function(x, y, replicates = 10000L, seed = 42L) {
  stopifnot(length(x) == length(y), length(x) > 2L, !anyNA(x), !anyNA(y))
  set.seed(seed)
  n <- length(x)
  boot_rho <- replicate(replicates, {
    idx <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  boot_rho <- boot_rho[is.finite(boot_rho)]
  if (length(boot_rho) < 0.99 * replicates) {
    stop("More than 1% of bootstrap correlations were non-finite.")
  }
  c(
    lower = unname(quantile(boot_rho, 0.025, names = FALSE)),
    upper = unname(quantile(boot_rho, 0.975, names = FALSE)),
    finite_replicates = length(boot_rho)
  )
}

score_names <- c("Degradation_score", "Production_score")
score_labels <- c(
  Degradation_score = "Degradation score",
  Production_score = "Production score"
)
primary_rows <- lapply(seq_along(score_names), function(i) {
  score <- score_names[i]
  test <- suppressWarnings(cor.test(
    pre[[score]], pre$CYT_score,
    method = "spearman", exact = FALSE, alternative = "two.sided"
  ))
  ci <- bootstrap_spearman(
    pre[[score]], pre$CYT_score,
    replicates = 10000L, seed = 42L + i - 1L
  )
  data.frame(
    score = score_labels[[score]],
    n = nrow(pre),
    spearman_rho = unname(test$estimate),
    bootstrap_CI_lower = ci[["lower"]],
    bootstrap_CI_upper = ci[["upper"]],
    p_nominal = test$p.value,
    bootstrap_replicates = ci[["finite_replicates"]],
    stringsAsFactors = FALSE
  )
})
primary <- do.call(rbind, primary_rows)
primary$p_Holm <- p.adjust(primary$p_nominal, method = "holm")
write.csv(primary, file.path(out_dir, "primary_score_correlations.csv"), row.names = FALSE)
cat("Primary score correlations\n")
print(primary, row.names = FALSE, digits = 5)
cat("\n")

# =============================================================================
# 4. Covariate-adjusted HC3 models.
# =============================================================================
pre$CYT_std <- as.numeric(scale(pre$CYT_score))
for (score in score_names) {
  pre[[paste0(score, "_std")]] <- as.numeric(scale(pre[[score]]))
}

adjusted_models <- list()
adjusted_fits <- list()
for (score in score_names) {
  predictor <- paste0(score, "_std")
  fit <- lm(
    reformulate(
      c(predictor, "response_group", "therapy_short", "age", "gender"),
      response = "CYT_std"
    ),
    data = pre
  )
  robust <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
  estimate <- robust[predictor, "Estimate"]
  se <- robust[predictor, "Std. Error"]
  critical <- qt(0.975, df = df.residual(fit))
  adjusted_models[[score]] <- data.frame(
    score = score_labels[[score]],
    n = nobs(fit),
    standardized_beta = estimate,
    HC3_SE = se,
    HC3_CI_lower = estimate - critical * se,
    HC3_CI_upper = estimate + critical * se,
    HC3_p = robust[predictor, "Pr(>|t|)"],
    classical_p = coef(summary(fit))[predictor, "Pr(>|t|)"],
    adjusted_R_squared = summary(fit)$adj.r.squared,
    stringsAsFactors = FALSE
  )
  adjusted_fits[[score]] <- fit
}
adjusted <- do.call(rbind, adjusted_models)
adjusted$HC3_p_Holm <- p.adjust(adjusted$HC3_p, method = "holm")
write.csv(adjusted, file.path(out_dir, "adjusted_models_HC3.csv"), row.names = FALSE)
cat("Adjusted HC3 models\n")
print(adjusted, row.names = FALSE, digits = 5)
cat("\n")

# =============================================================================
# 5. Response/therapy strata and formal slope-interaction sensitivities.
# =============================================================================
stratified_rows <- list()
for (stratum_variable in c("response_group", "therapy_short")) {
  for (level_value in levels(droplevels(pre[[stratum_variable]]))) {
    keep <- pre[[stratum_variable]] == level_value
    test <- suppressWarnings(cor.test(
      pre$Degradation_score[keep], pre$CYT_score[keep],
      method = "spearman", exact = FALSE, alternative = "two.sided"
    ))
    stratified_rows[[length(stratified_rows) + 1L]] <- data.frame(
      stratum_variable = stratum_variable,
      stratum = level_value,
      n = sum(keep),
      spearman_rho = unname(test$estimate),
      p_nominal = test$p.value,
      stringsAsFactors = FALSE
    )
  }
}
stratified <- do.call(rbind, stratified_rows)
stratified$p_BH_within_4_descriptive_tests <- p.adjust(
  stratified$p_nominal, method = "BH"
)
write.csv(stratified, file.path(out_dir, "stratified_correlations.csv"), row.names = FALSE)

interaction_model <- function(interactor) {
  fit <- lm(
    as.formula(paste(
      "CYT_std ~ Degradation_score_std *", interactor,
      "+ response_group + therapy_short + age + gender"
    )),
    data = pre
  )
  robust <- lmtest::coeftest(fit, vcov. = sandwich::vcovHC(fit, type = "HC3"))
  term <- grep(
    paste0("^Degradation_score_std:", interactor),
    rownames(robust), value = TRUE
  )
  if (length(term) != 1L) stop("Could not identify interaction term for ", interactor)
  data.frame(
    interaction = paste("Degradation score x", interactor),
    term = term,
    estimate = robust[term, "Estimate"],
    HC3_SE = robust[term, "Std. Error"],
    HC3_p = robust[term, "Pr(>|t|)"],
    n = nobs(fit),
    stringsAsFactors = FALSE
  )
}
interactions <- rbind(
  interaction_model("response_group"),
  interaction_model("therapy_short")
)
interactions$HC3_p_Holm <- p.adjust(interactions$HC3_p, method = "holm")
write.csv(interactions, file.path(out_dir, "interaction_sensitivity.csv"), row.names = FALSE)
cat("Stratified correlations\n")
print(stratified, row.names = FALSE, digits = 5)
cat("\nInteraction sensitivity\n")
print(interactions, row.names = FALSE, digits = 5)
cat("\n")

# =============================================================================
# 6. All-91 repeated-sample mixed-model sensitivity.
# =============================================================================
dat$CYT_std_all <- as.numeric(scale(dat$CYT_score))
for (score in score_names) {
  dat[[paste0(score, "_std_all")]] <- as.numeric(scale(dat[[score]]))
}

mixed_rows <- list()
for (score in score_names) {
  predictor <- paste0(score, "_std_all")
  fit <- lmerTest::lmer(
    reformulate(
      c(
        predictor, "response_group", "timepoint", "therapy_short",
        "age", "gender", "(1 | patient_name)"
      ),
      response = "CYT_std_all"
    ),
    data = dat,
    REML = TRUE
  )
  coef_table <- coef(summary(fit))
  ci <- suppressMessages(confint(fit, parm = predictor, method = "Wald"))
  mixed_rows[[score]] <- data.frame(
    score = score_labels[[score]],
    samples = nrow(dat),
    patients = length(unique(dat$patient_name)),
    repeated_rows = sum(duplicated(dat$patient_name)),
    standardized_beta = coef_table[predictor, "Estimate"],
    SE = coef_table[predictor, "Std. Error"],
    df = coef_table[predictor, "df"],
    CI_lower = ci[1, 1],
    CI_upper = ci[1, 2],
    p = coef_table[predictor, "Pr(>|t|)"],
    singular_fit = lme4::isSingular(fit, tol = 1e-4),
    stringsAsFactors = FALSE
  )
}
mixed <- do.call(rbind, mixed_rows)
mixed$p_Holm <- p.adjust(mixed$p, method = "holm")
write.csv(mixed, file.path(out_dir, "all91_mixed_model_sensitivity.csv"), row.names = FALSE)
cat("All-91 mixed-model sensitivity\n")
print(mixed, row.names = FALSE, digits = 5)
cat("\n")

# =============================================================================
# 7. Influence and leave-one-patient-out diagnostics for the primary result.
# =============================================================================
deg_fit <- adjusted_fits[["Degradation_score"]]
cook <- cooks.distance(deg_fit)
cook_threshold <- 4 / nrow(pre)
influence_table <- data.frame(
  sample_id = pre$sample_id,
  patient_name = pre$patient_name,
  cooks_distance = unname(cook),
  above_4_over_n = unname(cook > cook_threshold),
  stringsAsFactors = FALSE
)
influence_table <- influence_table[order(-influence_table$cooks_distance), ]
write.csv(
  influence_table, file.path(out_dir, "influence_Cooks_distance.csv"),
  row.names = FALSE
)

loo <- lapply(seq_len(nrow(pre)), function(i) {
  test <- suppressWarnings(cor.test(
    pre$Degradation_score[-i], pre$CYT_score[-i],
    method = "spearman", exact = FALSE, alternative = "two.sided"
  ))
  data.frame(
    omitted_sample_id = pre$sample_id[i],
    omitted_patient = pre$patient_name[i],
    spearman_rho = unname(test$estimate),
    p_nominal = test$p.value,
    stringsAsFactors = FALSE
  )
})
loo <- do.call(rbind, loo)
write.csv(loo, file.path(out_dir, "leave_one_patient_out.csv"), row.names = FALSE)

influence_summary <- data.frame(
  metric = c(
    "Cook threshold 4/n", "Maximum Cook distance", "Patients above Cook 4/n",
    "LOO minimum rho", "LOO maximum rho", "LOO maximum nominal P",
    "LOO analyses with P < 0.05"
  ),
  value = c(
    cook_threshold, max(cook), sum(cook > cook_threshold),
    min(loo$spearman_rho), max(loo$spearman_rho), max(loo$p_nominal),
    sum(loo$p_nominal < 0.05)
  ),
  stringsAsFactors = FALSE
)
write.csv(influence_summary, file.path(out_dir, "influence_summary.csv"), row.names = FALSE)
cat("Influence summary\n")
print(influence_summary, row.names = FALSE, digits = 5)
cat("\n")

# =============================================================================
# 8. Secondary component checks, BH-adjusted as one exploratory family.
# =============================================================================
secondary_pairs <- list(
  c("SARDH", "CYT_score"), c("PIPOX", "CYT_score"),
  c("GNMT", "CYT_score"), c("DMGDH", "CYT_score"),
  c("Degradation_score", "GZMA_TPM"),
  c("Degradation_score", "PRF1_TPM")
)
secondary <- do.call(rbind, lapply(secondary_pairs, function(pair) {
  test <- suppressWarnings(cor.test(
    pre[[pair[1]]], pre[[pair[2]]],
    method = "spearman", exact = FALSE, alternative = "two.sided"
  ))
  data.frame(
    predictor = pair[1], outcome = pair[2], n = nrow(pre),
    spearman_rho = unname(test$estimate), p_nominal = test$p.value,
    stringsAsFactors = FALSE
  )
}))
secondary$p_BH <- p.adjust(secondary$p_nominal, method = "BH")
write.csv(
  secondary, file.path(out_dir, "secondary_component_correlations.csv"),
  row.names = FALSE
)

# Unit sensitivity: the primary CYT uses TPM, while this table confirms whether
# the conclusion changes if the supplied FPKM units are used directly.
cyt_unit_sensitivity <- do.call(rbind, lapply(
  c("CYT_score", "CYT_score_FPKM_sensitivity"),
  function(outcome) {
    test <- suppressWarnings(cor.test(
      pre$Degradation_score, pre[[outcome]],
      method = "spearman", exact = FALSE, alternative = "two.sided"
    ))
    data.frame(
      CYT_definition = ifelse(
        outcome == "CYT_score", "TPM (primary)", "FPKM (sensitivity)"
      ),
      n = nrow(pre), spearman_rho = unname(test$estimate),
      p_nominal = test$p.value, stringsAsFactors = FALSE
    )
  }
))
write.csv(
  cyt_unit_sensitivity, file.path(out_dir, "CYT_unit_sensitivity.csv"),
  row.names = FALSE
)

# =============================================================================
# 9. Source-data exports.
# =============================================================================
pre_source <- pre[, c(
  "sample_id", "patient_name", "timepoint", "recist", "response_group",
  "therapy_short", "age", "gender", "SARDH", "PIPOX", "GNMT", "DMGDH",
  "GZMA_FPKM", "PRF1_FPKM", "GZMA_TPM", "PRF1_TPM",
  "Degradation_score", "Production_score", "CYT_score",
  "CYT_score_FPKM_sensitivity"
)]
names(pre_source)[names(pre_source) == "SARDH"] <- "SARDH_log2_FPKM_plus1"
names(pre_source)[names(pre_source) == "PIPOX"] <- "PIPOX_log2_FPKM_plus1"
names(pre_source)[names(pre_source) == "GNMT"] <- "GNMT_log2_FPKM_plus1"
names(pre_source)[names(pre_source) == "DMGDH"] <- "DMGDH_log2_FPKM_plus1"
write.csv(
  pre_source, file.path(out_dir, "SourceData_Fig2d_CYT_PRE73.csv"), row.names = FALSE
)

all91_source <- dat[, c(
  "sample_id", "patient_name", "timepoint", "recist", "response_group",
  "therapy_short", "age", "gender", "Degradation_score", "Production_score",
  "CYT_score", "CYT_score_FPKM_sensitivity"
)]
write.csv(
  all91_source, file.path(out_dir, "SourceData_CYT_all91_sensitivity.csv"),
  row.names = FALSE
)

# =============================================================================
# 10. Candidate main panel: raw PRE-only data, no cutpoint.
# =============================================================================
workspace_root <- normalizePath(file.path(project_dir, "..", "..", ".."))
theme_candidates <- Sys.glob(file.path(
  workspace_root, "*Metabolomics*", "_shared", "theme_nc_26.08.18.R"
))
if (length(theme_candidates) != 1L) {
  stop("Expected exactly one shared theme file; found ", length(theme_candidates))
}
source(theme_candidates)

plot_data <- rbind(
  data.frame(pre_source, score_type = "Degradation",
             metabolism_score = pre$Degradation_score),
  data.frame(pre_source, score_type = "Production",
             metabolism_score = pre$Production_score)
)
plot_data$score_type <- factor(
  plot_data$score_type, levels = c("Degradation", "Production")
)
plot_data$response_group <- factor(plot_data$response_group, levels = c("R", "NR"))
plot_data$therapy_short <- factor(
  plot_data$therapy_short, levels = c("antiPD1", "combo")
)

format_p_panel <- function(p) {
  ifelse(p < 0.001, "<0.001", sprintf("%.3f", p))
}
annotation <- merge(
  primary[, c("score", "spearman_rho", "p_Holm")],
  adjusted[, c("score", "HC3_p_Holm")], by = "score", sort = FALSE
)
annotation$score_type <- factor(
  sub(" score$", "", annotation$score), levels = c("Degradation", "Production")
)
annotation$label <- sprintf(
  "rho = %.2f\nHolm P %s\nAdjusted P %s",
  annotation$spearman_rho,
  ifelse(annotation$p_Holm < 0.001, "< 0.001",
         paste("=", format_p_panel(annotation$p_Holm))),
  ifelse(annotation$HC3_p_Holm < 0.001, "< 0.001",
         paste("=", format_p_panel(annotation$HC3_p_Holm)))
)

panel <- ggplot(plot_data, aes(x = metabolism_score, y = CYT_score)) +
  geom_smooth(
    method = "lm", formula = y ~ x, colour = "#383838", fill = "grey80",
    linewidth = 0.45, alpha = 0.35, se = TRUE
  ) +
  geom_point(
    aes(colour = response_group, shape = therapy_short),
    size = 0.9, alpha = 0.72, stroke = 0.18
  ) +
  geom_text(
    data = annotation, aes(x = -Inf, y = Inf, label = label),
    inherit.aes = FALSE, hjust = -0.08, vjust = 1.12,
    family = "Arial", size = mm_text(NC_ANNOT_PT), lineheight = 0.95
  ) +
  facet_wrap(~score_type, nrow = 1) +
  scale_colour_manual(
    values = RESPONSE_COLORS, breaks = c("R", "NR"), labels = c("R", "NR")
  ) +
  scale_shape_manual(
    values = c(antiPD1 = 16, combo = 17), breaks = c("antiPD1", "combo"),
    labels = c("anti-PD-1", "Combination")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.14))) +
  labs(
    x = "Sarcosine-metabolism score (z)",
    y = "Tumour CYT score",
    colour = "Response", shape = "Therapy"
  ) +
  theme_nc() +
  theme(
    legend.position = "top",
    legend.title = element_text(
      family = "Arial", size = NC_LEGEND_PT, face = "bold"
    ),
    legend.spacing.x = grid::unit(4, "pt"),
    legend.box.spacing = grid::unit(0, "pt"),
    plot.margin = margin(1, 3, 1, 1)
  ) +
  guides(
    colour = guide_legend(order = 1, nrow = 1, override.aes = list(size = 1.5)),
    shape = guide_legend(order = 2, nrow = 1, override.aes = list(size = 1.5))
  )

figure_stem <- file.path(out_dir, "Fig2d_host_degradation_CYT_candidate")
ggsave(
  paste0(figure_stem, ".png"), panel,
  width = 4.45, height = 2.35, units = "in", dpi = 600,
  bg = "white", device = ragg::agg_png
)
ggsave(
  paste0(figure_stem, ".svg"), panel,
  width = 4.45, height = 2.35, units = "in",
  bg = "white", device = svglite::svglite
)
# Register the system Arial AFM metrics for R's vector PDF device. cairo_pdf
# fails on this macOS R build with a CID-font error, whereas the Type1 mapping
# below is available and verified before export.
grDevices::pdfFonts(Arial = grDevices::pdfFonts("ArialMT")[[1]])
ggsave(
  paste0(figure_stem, ".pdf"), panel,
  width = 4.45, height = 2.35, units = "in",
  bg = "white", device = grDevices::pdf, family = "Arial",
  useDingbats = FALSE
)
stopifnot(
  file.info(paste0(figure_stem, ".png"))$size > 5000,
  file.info(paste0(figure_stem, ".pdf"))$size > 5000,
  file.info(paste0(figure_stem, ".svg"))$size > 5000
)

# =============================================================================
# 11. Prespecified decision checks and reproducibility records.
# =============================================================================
deg_primary <- primary[primary$score == "Degradation score", ]
prod_primary <- primary[primary$score == "Production score", ]
deg_adjusted <- adjusted[adjusted$score == "Degradation score", ]
deg_mixed <- mixed[mixed$score == "Degradation score", ]
decision <- data.frame(
  criterion = c(
    "Primary degradation rho is positive",
    "Primary degradation Holm P < 0.05",
    "Adjusted degradation HC3 Holm P < 0.05",
    "All-91 degradation mixed-model Holm P < 0.05",
    "All-91 degradation mixed model is not singular",
    "Every leave-one-patient-out rho is positive",
    "Every leave-one-patient-out P < 0.05",
    "Production specificity control Holm P >= 0.05"
  ),
  passed = c(
    deg_primary$spearman_rho > 0,
    deg_primary$p_Holm < 0.05,
    deg_adjusted$HC3_p_Holm < 0.05,
    deg_mixed$p_Holm < 0.05,
    !deg_mixed$singular_fit,
    min(loo$spearman_rho) > 0,
    max(loo$p_nominal) < 0.05,
    prod_primary$p_Holm >= 0.05
  ),
  stringsAsFactors = FALSE
)
write.csv(
  decision, file.path(out_dir, "main_figure_candidate_decision_checks.csv"),
  row.names = FALSE
)

readme <- c(
  "Candidate Figure 2d: host sarcosine-degradation score versus tumour CYT",
  "=========================================================================",
  "",
  "Scope",
  "- TIGER melanoma PRJEB23709 bulk RNA-seq (FPKM).",
  "- Primary unit: one pretreatment biopsy per patient (n=73; R=40, NR=33).",
  "- The full deposited FPKM matrix was converted within each sample to TPM by",
  "  dividing each FPKM by the sample's total FPKM and multiplying by 1e6.",
  "- CYT is the log2 geometric mean of GZMA and PRF1 TPM, matching the original",
  "  Rooney unit definition. Both genes were positive, so no pseudocount was used.",
  "- Degradation and production scores use PRE-anchored per-gene z-scores.",
  "",
  "Inference",
  "- Primary: two-sided Spearman correlation; Holm correction across degradation",
  "  and production score tests.",
  "- Adjusted sensitivity: HC3 linear model controlling response, therapy, age, sex.",
  "- All-91 sensitivity: timepoint-adjusted patient-random-intercept mixed model.",
  "- Bootstrap intervals are percentile intervals from 10,000 patient resamples.",
  "",
  "Interpretive boundary",
  "- This is an association in bulk tumour RNA. It does not establish that sarcosine",
  "  degradation causes cytolytic activity, and it cannot distinguish immune-cell",
  "  abundance from immune-cell activation or exclude tumour-purity confounding.",
  "- The regression line is a visual guide; the displayed unadjusted test is Spearman.",
  "",
  "Main outputs",
  "- Fig2d_host_degradation_CYT_candidate.png/.pdf/.svg",
  "- SourceData_Fig2d_CYT_PRE73.csv",
  "- primary_score_correlations.csv",
  "- adjusted_models_HC3.csv",
  "- all91_mixed_model_sensitivity.csv",
  "- CYT_unit_sensitivity.csv",
  "- main_figure_candidate_decision_checks.csv"
)
writeLines(readme, file.path(out_dir, "README.txt"), useBytes = TRUE)
writeLines(
  capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"),
  useBytes = TRUE
)

input_checksums <- data.frame(
  file = normalizePath(required_files),
  sha256 = vapply(
    required_files, digest::digest, character(1), algo = "sha256", file = TRUE
  ),
  stringsAsFactors = FALSE
)
write.csv(input_checksums, file.path(out_dir, "input_SHA256.csv"), row.names = FALSE)

output_files <- list.files(out_dir, full.names = TRUE)
output_files <- output_files[
  !basename(output_files) %in% c("output_SHA256.csv", "analysis_log.txt") &
    !file.info(output_files)$isdir
]
output_checksums <- data.frame(
  file = basename(output_files),
  sha256 = vapply(
    output_files, digest::digest, character(1), algo = "sha256", file = TRUE
  ),
  stringsAsFactors = FALSE
)
write.csv(output_checksums, file.path(out_dir, "output_SHA256.csv"), row.names = FALSE)

cat("Decision checks\n")
print(decision, row.names = FALSE)
cat("\nOutputs written to:\n", out_dir, "\n", sep = "")
cat("\nSession information\n")
print(sessionInfo())
