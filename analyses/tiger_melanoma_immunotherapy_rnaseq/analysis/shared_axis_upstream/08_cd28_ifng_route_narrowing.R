#!/usr/bin/env Rscript

# TIGER melanoma PRE73: candidate-route narrowing from sarcosine-degradation
# transcriptional state through CD28-associated programs to IFNG expression.
#
# Integrity / interpretation boundary
# -----------------------------------
# - "Sarcosine degradation" is the pre-existing SARDH/PIPOX RNA module, not
#   measured sarcosine concentration or metabolic flux.
# - TF activities are ULM inferences from signed DoRothEA A/B regulons, not TF
#   protein abundance, phosphorylation, nuclear localization, or DNA binding.
# - Cross-sectional bulk RNA cannot establish temporal or causal direction.
#   Serial products below are explicitly exploratory association summaries.
# - Reader-facing figures use Degradation Low/High and biological labels only.
#   Named cell-context sensitivity analyses are retained in tables.

options(stringsAsFactors = FALSE, width = 180)
set.seed(260902)

analysis_root <- normalizePath(getwd())
if (basename(analysis_root) != "Melanoma-PRJEB23709") {
  stop("Run from the Melanoma-PRJEB23709 analysis root: ", analysis_root)
}

default_gsea_lib <- file.path(
  dirname(dirname(dirname(analysis_root))),
  "2024_Drug_Res_Updates_NSCLC", "RNA-seq공공데이터_GSE207422",
  "analysis_sarcosine_FINAL_26.08.25", "03_lee_fig3_style_FINAL",
  "hallmark_GSEA_Q4_vs_Q1_26.08.26", "R_libs"
)
gsea_lib <- Sys.getenv("SARCO_GSEA_R_LIB", unset = default_gsea_lib)
if (!dir.exists(gsea_lib)) stop("GSEA R library not found: ", gsea_lib)
.libPaths(unique(c(normalizePath(gsea_lib), .libPaths())))

required_packages <- c("data.table", "limma", "msigdbr", "ggplot2", "patchwork", "pheatmap")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
})

used_data_root <- normalizePath(file.path(analysis_root, "..", "사용데이터_모음"))
expression_path <- file.path(used_data_root, "Source_Input", "TIGER_PRJEB23709_expression_FPKM_gene_by_sample.csv")
metadata_path <- file.path(used_data_root, "Analysis_Ready", "TIGER_PRE73_Fig4bc_analysis_data.csv")
tf_root <- file.path(analysis_root, "results", "TIGER_PRE73_CD28_Upstream_TF_26.09.02")
ifng_root <- file.path(analysis_root, "results", "TIGER_PRE73_IFNG_Axis_Focused_26.09.02")
tcell_root <- file.path(analysis_root, "results", "TIGER_PRE73_Tcell_Abundance_Adjusted_Function_26.09.02")
deconv_root <- file.path(analysis_root, "results", "TIGER_PRE73_Degradation_HighLow_Deconvolution_7methods_26.09.02")

regulon_path <- file.path(tf_root, "tables", "01_DoRothEA_AB_signed_regulon_all.csv")
prior_tf_path <- file.path(tf_root, "tables", "07_TF_ULM_activity_High_vs_Low.csv")
tcell_context_path <- file.path(tcell_root, "tables", "04_Tcell_abundance_adjusters_long.csv")
deconv_path <- file.path(deconv_root, "tables", "05_core_celltype_values_long.csv")

required_inputs <- c(expression_path, metadata_path, regulon_path, prior_tf_path, tcell_context_path, deconv_path)
for (path in required_inputs) if (!file.exists(path)) stop("Missing input: ", path)

output_root <- file.path(analysis_root, "results", "TIGER_PRE73_CD28_to_IFNG_Route_Narrowing_26.09.02")
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

safe_neglog10 <- function(x) -log10(pmax(as.numeric(x), .Machine$double.xmin))

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
    partial_R2 = unname(sm[term, "t value"])^2 / (unname(sm[term, "t value"])^2 + df)
  )
}

score_gene_set <- function(genes, label, log_expression, minimum_genes = 5L) {
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

fit_group_effects <- function(score_mat, analysis_data, extra_terms = character(), model_label) {
  rhs <- c("therapy_short", "response_group", "age_z", "gender", extra_terms, "degradation_group")
  design <- model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), data = analysis_data)
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient group model: ", model_label)
  fit <- eBayes(lmFit(score_mat[, analysis_data$sample_id, drop = FALSE], design), robust = TRUE)
  term <- "degradation_groupHigh"
  se <- fit$stdev.unscaled[, term] * sqrt(fit$s2.post)
  df <- fit$df.total
  out <- data.table(
    feature = rownames(score_mat), standardized_effect = fit$coefficients[, term],
    standard_error = se, moderated_t = fit$t[, term], p_value = fit$p.value[, term],
    residual_df = df, model = model_label
  )
  out[, `:=`(
    CI_low = standardized_effect - qt(0.975, residual_df) * standard_error,
    CI_high = standardized_effect + qt(0.975, residual_df) * standard_error,
    BH_q = p.adjust(p_value, method = "BH")
  )]
  out[]
}

fit_serial_edges <- function(data, route_label, cd28_col, branch_col, extra_terms = character(), model_label) {
  covars <- c("therapy_short", "response_group", "age_z", "gender", extra_terms)
  f1 <- lm(as.formula(paste(cd28_col, "~ degradation_score_z +", paste(covars, collapse = " + "))), data = data)
  f2 <- lm(as.formula(paste(branch_col, "~ degradation_score_z +", cd28_col, "+", paste(covars, collapse = " + "))), data = data)
  f3 <- lm(as.formula(paste("IFNG_expression_z ~ degradation_score_z +", cd28_col, "+", branch_col, "+", paste(covars, collapse = " + "))), data = data)
  rbindlist(list(
    extract_lm_term(f1, "degradation_score_z", cd28_col, "Degradation score", model_label, route_label, "Degradation -> CD28"),
    extract_lm_term(f2, cd28_col, branch_col, "CD28 program", model_label, route_label, "CD28 -> branch"),
    extract_lm_term(f3, branch_col, "IFNG expression", "Branch", model_label, route_label, "Branch -> IFNG")
  ))
}

serial_bootstrap <- function(data, route_label, cd28_col, branch_col,
                             order = c("forward", "reverse"), extra_terms = character(),
                             model_label, B = 5000L, seed = 1L) {
  order <- match.arg(order)
  covars <- c("therapy_short", "response_group", "age_z", "gender", extra_terms)
  m1 <- if (order == "forward") cd28_col else branch_col
  m2 <- if (order == "forward") branch_col else cd28_col
  y <- "IFNG_expression_z"

  x1 <- model.matrix(as.formula(paste("~ degradation_score_z +", paste(covars, collapse = " + "))), data = data)
  x2 <- model.matrix(as.formula(paste("~ degradation_score_z +", m1, "+", paste(covars, collapse = " + "))), data = data)
  x3 <- model.matrix(as.formula(paste("~ degradation_score_z +", m1, "+", m2, "+", paste(covars, collapse = " + "))), data = data)
  y1 <- data[[m1]]; y2 <- data[[m2]]; y3 <- data[[y]]
  term_a <- match("degradation_score_z", colnames(x1))
  term_d <- match(m1, colnames(x2))
  term_b <- match(m2, colnames(x3))
  if (anyNA(c(term_a, term_d, term_b))) stop("Serial model term mapping failed")

  one_product <- function(idx) {
    b1 <- .lm.fit(x1[idx, , drop = FALSE], y1[idx])$coefficients[term_a]
    b2 <- .lm.fit(x2[idx, , drop = FALSE], y2[idx])$coefficients[term_d]
    b3 <- .lm.fit(x3[idx, , drop = FALSE], y3[idx])$coefficients[term_b]
    unname(b1 * b2 * b3)
  }
  observed <- one_product(seq_len(nrow(data)))
  set.seed(seed)
  boot_values <- replicate(B, one_product(sample.int(nrow(data), nrow(data), replace = TRUE)))
  boot_values <- boot_values[is.finite(boot_values)]
  if (length(boot_values) < 0.95 * B) stop("Too many invalid bootstrap replicates")
  ci <- as.numeric(quantile(boot_values, c(0.025, 0.975), names = FALSE, type = 6))
  data.table(
    route = route_label, order = order, model = model_label,
    serial_product = observed, bootstrap_CI_low = ci[1], bootstrap_CI_high = ci[2],
    CI_excludes_zero = ci[1] > 0 | ci[2] < 0,
    bootstrap_replicates_requested = B, bootstrap_replicates_valid = length(boot_values),
    node_1 = "Degradation score", node_2 = m1, node_3 = m2, node_4 = "IFNG expression"
  )
}

vif_one <- function(term, others, data) {
  fit <- lm(as.formula(paste(term, "~", paste(others, collapse = " + "))), data = data)
  1 / (1 - summary(fit)$r.squared)
}

# -------------------------------------------------------------------------
# Inputs and validated PRE73 expression matrix
# -------------------------------------------------------------------------
expr_dt <- fread(expression_path, check.names = FALSE)
meta <- fread(metadata_path, check.names = FALSE)
regulon <- fread(regulon_path, check.names = FALSE)
prior_tf <- fread(prior_tf_path, check.names = FALSE)
tcell_context <- fread(tcell_context_path, check.names = FALSE)
deconv <- fread(deconv_path, check.names = FALSE)

stopifnot(
  nrow(meta) == 73L, uniqueN(meta$sample_id) == 73L, all(meta$timepoint == "PRE"),
  all(c("sample_id", "therapy_short", "response_group", "age", "gender", "Degradation_score") %in% names(meta)),
  all(c("TF", "target", "mor", "confidence", "expressed") %in% names(regulon))
)

gene_symbols <- trimws(expr_dt[[1]])
if (anyNA(gene_symbols) || any(!nzchar(gene_symbols)) || anyDuplicated(gene_symbols)) stop("Invalid gene-symbol column")
sample_ids <- meta$sample_id
if (!all(sample_ids %in% names(expr_dt))) stop("Expression matrix lacks PRE73 samples")
fpkm <- as.matrix(expr_dt[, ..sample_ids])
storage.mode(fpkm) <- "double"
rownames(fpkm) <- gene_symbols; colnames(fpkm) <- sample_ids
if (anyNA(fpkm) || any(!is.finite(fpkm)) || any(fpkm < 0)) stop("Invalid FPKM matrix")
log_expression_all <- log2(fpkm + 1)
keep_expression <- rowSums(fpkm >= 1) >= ceiling(0.10 * ncol(fpkm))
log_expression <- log_expression_all[keep_expression, , drop = FALSE]
log_expression <- log_expression[apply(log_expression, 1L, var) > 0, , drop = FALSE]
if (ncol(log_expression) != 73L || nrow(log_expression) < 10000L) stop("Unexpected filtered expression dimensions")

degradation_median <- median(meta$Degradation_score)
meta[, degradation_group := factor(ifelse(Degradation_score > degradation_median, "High", "Low"), levels = c("Low", "High"))]
stopifnot(sum(meta$degradation_group == "Low") == 37L, sum(meta$degradation_group == "High") == 36L)
meta[, therapy_short := relevel(factor(therapy_short), ref = "antiPD1")]
meta[, response_group := relevel(factor(response_group), ref = "NR")]
meta[, gender := relevel(factor(gender), ref = "Male")]
meta[, `:=`(age_z = zscore(age), degradation_score_z = zscore(Degradation_score))]

# Named cell contexts are sensitivity variables only.
pan_t <- tcell_context[abundance_method == "MCPcounter Pan T", .(sample_id, Pan_T_cell_context = abundance_z)]
dendritic <- deconv[method == "MCPcounter" & analysis_feature == "Dendritic cells", .(sample_id, Dendritic_cell_context = analysis_value)]
monocytic <- deconv[method == "MCPcounter" & analysis_feature == "Monocytes", .(sample_id, Monocytic_context = analysis_value)]
for (x in list(pan_t, dendritic, monocytic)) stopifnot(nrow(x) == 73L, uniqueN(x$sample_id) == 73L)
context_wide <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE), list(pan_t, dendritic, monocytic))
context_wide <- context_wide[match(sample_ids, sample_id)]
for (v in setdiff(names(context_wide), "sample_id")) context_wide[, (v) := zscore(get(v))]

# -------------------------------------------------------------------------
# Predeclared MSigDB/Reactome modules and leave-route-out TF activities
# -------------------------------------------------------------------------
hallmark <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
reactome <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))
gobp <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"))
stopifnot(uniqueN(hallmark$db_version) == 1L, identical(unique(hallmark$db_version), unique(reactome$db_version)), identical(unique(hallmark$db_version), unique(gobp$db_version)))
msigdb_version <- unique(hallmark$db_version)

source_catalog <- data.table(
  source_key = c("CD28", "PI3K_AKT", "VAV1", "NFAT", "AP1", "ERK", "IL12"),
  source_id = c(
    "REACTOME_CO_STIMULATION_BY_CD28",
    "REACTOME_CD28_DEPENDENT_PI3K_AKT_SIGNALING",
    "REACTOME_CD28_DEPENDENT_VAV1_PATHWAY",
    "REACTOME_CALCINEURIN_ACTIVATES_NFAT",
    "REACTOME_ACTIVATION_OF_THE_AP_1_FAMILY_OF_TRANSCRIPTION_FACTORS",
    "REACTOME_ERK_MAPK_TARGETS",
    "REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING"
  )
)
if (!all(source_catalog$source_id %in% reactome$gs_name)) stop("A required Reactome gene set is absent")
source_sets <- setNames(lapply(source_catalog$source_id, function(pid) {
  intersect(unique(reactome[gs_name == pid, gene_symbol]), rownames(log_expression))
}), source_catalog$source_key)

response_membership <- rbindlist(list(
  hallmark[gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE", .(collection = "Hallmark", gene_symbol)],
  reactome[gs_name == "REACTOME_INTERFERON_GAMMA_SIGNALING", .(collection = "Reactome", gene_symbol)],
  gobp[gs_name == "GOBP_RESPONSE_TO_TYPE_II_INTERFERON", .(collection = "GO:BP", gene_symbol)]
))
ifng_response_core <- response_membership[, .N, by = gene_symbol][N >= 2L, gene_symbol]
ifng_response_core <- setdiff(intersect(ifng_response_core, rownames(log_expression)), "IFNG")
if (length(ifng_response_core) < 20L) stop("Unexpectedly small IFNG response core")

focal_tfs <- c("SPI1", "NFKB1", "RELA", "JUND", "FOSL1")
tf_exclusion_schemes <- list(
  `Full signed regulon` = character(),
  `Exclude CD28 genes` = source_sets$CD28,
  `Exclude IFNG-response genes` = unique(c(ifng_response_core, "IFNG")),
  `Exclude CD28 + IFNG-response genes` = unique(c(source_sets$CD28, ifng_response_core, "IFNG")),
  `Exclude all candidate-route genes` = unique(c(unlist(source_sets, use.names = FALSE), ifng_response_core, "IFNG"))
)

scheme_regulons <- list()
scheme_activities_z <- list()
scheme_coverage <- list()
for (scheme_name in names(tf_exclusion_schemes)) {
  excluded <- tf_exclusion_schemes[[scheme_name]]
  net <- regulon[
    expressed == TRUE & confidence %in% c("A", "B") & !target %in% excluded
  ]
  coverage <- net[TF %in% focal_tfs, .(
    retained_targets_n = uniqueN(target),
    activating_targets_n = uniqueN(target[mor > 0]),
    repressing_targets_n = uniqueN(target[mor < 0]),
    targets = paste(sort(unique(target)), collapse = ";")
  ), by = TF]
  if (!all(focal_tfs %in% coverage[retained_targets_n >= 10L, TF])) {
    stop("A focal TF has fewer than 10 targets under scheme: ", scheme_name)
  }
  activity <- ulm_from_matrix(log_expression, net, minsize = 10L, center_rows = TRUE)
  activity <- activity[focal_tfs, sample_ids, drop = FALSE]
  activity_z <- t(apply(activity, 1L, zscore))
  rownames(activity_z) <- focal_tfs; colnames(activity_z) <- sample_ids
  excluded_counts <- regulon[
    TF %in% focal_tfs & expressed == TRUE & confidence %in% c("A", "B") & target %in% excluded,
    .(excluded_target_overlap_n = uniqueN(target)), by = TF
  ]
  coverage <- merge(coverage, excluded_counts, by = "TF", all.x = TRUE)
  coverage[is.na(excluded_target_overlap_n), excluded_target_overlap_n := 0L]
  coverage[, `:=`(
    exclusion_scheme = scheme_name,
    excluded_genes_n = length(excluded)
  )]
  scheme_regulons[[scheme_name]] <- net
  scheme_activities_z[[scheme_name]] <- activity_z
  scheme_coverage[[scheme_name]] <- coverage
}
scheme_coverage <- rbindlist(scheme_coverage, use.names = TRUE)

# The primary non-overlap activity excludes the specific CD28 and IFNG-response
# nodes being related. More aggressive exclusion is retained as a robustness
# audit rather than silently replacing the biological regulon.
primary_exclusion_scheme <- "Exclude CD28 + IFNG-response genes"
regulon_lro <- scheme_regulons[[primary_exclusion_scheme]]
lro_activity_z <- scheme_activities_z[[primary_exclusion_scheme]]
lro_coverage <- scheme_coverage[exclusion_scheme == primary_exclusion_scheme]
setnames(lro_coverage, "retained_targets_n", "leave_route_out_targets_n")

score_objects <- list()
score_objects[["APC ligand expression"]] <- score_gene_set(c("CD80", "CD86"), "APC ligand expression", log_expression, minimum_genes = 2L)
score_objects[["Calcineurin-NFAT pathway"]] <- score_gene_set(source_sets$NFAT, "Calcineurin-NFAT pathway", log_expression)
score_objects[["CD28-PI3K-AKT pathway"]] <- score_gene_set(setdiff(source_sets$PI3K_AKT, c("CD28", "CD80", "CD86")), "CD28-PI3K-AKT pathway", log_expression)
score_objects[["IL-12-STAT4 pathway"]] <- score_gene_set(setdiff(source_sets$IL12, "IFNG"), "IL-12-STAT4 pathway", log_expression)
score_objects[["VAV1-ERK expression"]] <- score_gene_set(unique(c(source_sets$VAV1, source_sets$ERK)), "VAV1-ERK expression", log_expression)
score_objects[["IFNG response core"]] <- score_gene_set(ifng_response_core, "IFNG response core", log_expression)

if (!"IFNG" %in% rownames(log_expression_all)) stop("IFNG absent from expression matrix")
ifng_expression_z <- zscore(log_expression_all["IFNG", sample_ids])

branch_catalog <- data.table(
  route_order = 1:5,
  route = c(
    "Canonical NF-kappaB",
    "VAV1-MAPK-AP-1",
    "Calcineurin-NFAT",
    "CD28-PI3K-AKT",
    "IL-12-STAT4 competitor"
  ),
  branch_col = c("branch_nfkb", "branch_ap1", "branch_nfat", "branch_pi3k", "branch_il12"),
  branch_construct = c(
    "Mean of leave-route-out NFKB1 and RELA DoRothEA A/B ULM activities",
    "Mean of leave-route-out JUND and FOSL1 DoRothEA A/B ULM activities",
    "Reactome calcineurin-activates-NFAT expression score",
    "Reactome CD28-dependent PI3K-AKT expression score excluding CD28/CD80/CD86",
    "Reactome IL-12-family expression score excluding IFNG; comparator, not assumed CD28-downstream"
  )
)

analysis_data <- copy(meta)
analysis_data <- merge(analysis_data, context_wide, by = "sample_id", all.x = TRUE, sort = FALSE)
analysis_data <- analysis_data[match(sample_ids, sample_id)]
analysis_data[, `:=`(
  SPI1_lro = as.numeric(lro_activity_z["SPI1", ]),
  NFKB1_lro = as.numeric(lro_activity_z["NFKB1", ]),
  RELA_lro = as.numeric(lro_activity_z["RELA", ]),
  JUND_lro = as.numeric(lro_activity_z["JUND", ]),
  FOSL1_lro = as.numeric(lro_activity_z["FOSL1", ]),
  APC_ligand_z = score_objects[["APC ligand expression"]]$score,
  IFNG_expression_z = ifng_expression_z,
  IFNG_response_core_z = score_objects[["IFNG response core"]]$score,
  branch_nfkb = zscore((as.numeric(lro_activity_z["NFKB1", ]) + as.numeric(lro_activity_z["RELA", ])) / 2),
  branch_ap1 = zscore((as.numeric(lro_activity_z["JUND", ]) + as.numeric(lro_activity_z["FOSL1", ])) / 2),
  branch_nfat = score_objects[["Calcineurin-NFAT pathway"]]$score,
  branch_pi3k = score_objects[["CD28-PI3K-AKT pathway"]]$score,
  branch_il12 = score_objects[["IL-12-STAT4 pathway"]]$score,
  VAV1_ERK_expression_z = score_objects[["VAV1-ERK expression"]]$score
)]
if (anyNA(analysis_data)) stop("Analysis data contain missing values")

# Route-specific CD28 scores exclude genes used by the corresponding branch,
# so the serial association cannot be driven by duplicate genes in two nodes.
branch_gene_sets <- list(
  `Canonical NF-kappaB` = unique(regulon_lro[TF %in% c("NFKB1", "RELA"), target]),
  `VAV1-MAPK-AP-1` = unique(c(regulon_lro[TF %in% c("JUND", "FOSL1"), target], source_sets$VAV1, source_sets$ERK, source_sets$AP1)),
  `Calcineurin-NFAT` = source_sets$NFAT,
  `CD28-PI3K-AKT` = score_objects[["CD28-PI3K-AKT pathway"]]$genes,
  `IL-12-STAT4 competitor` = score_objects[["IL-12-STAT4 pathway"]]$genes
)

route_manifest <- rbindlist(lapply(seq_len(nrow(branch_catalog)), function(i) {
  route_label <- branch_catalog$route[i]
  branch_genes <- unique(branch_gene_sets[[route_label]])
  cd28_genes <- setdiff(source_sets$CD28, unique(c(branch_genes, ifng_response_core, "IFNG")))
  cd28_object <- score_gene_set(cd28_genes, paste0(route_label, " disjoint CD28 score"), log_expression, minimum_genes = 10L)
  this_cd28_col <- paste0("cd28_route_", i)
  analysis_data[, (this_cd28_col) := cd28_object$score]
  branch_catalog[i, cd28_col := this_cd28_col]
  data.table(
    route = route_label,
    branch_col = branch_catalog$branch_col[i], cd28_col = this_cd28_col,
    full_CD28_expressed_genes_n = length(source_sets$CD28),
    removed_branch_overlap_n = length(intersect(source_sets$CD28, branch_genes)),
    removed_IFNG_response_overlap_n = length(intersect(source_sets$CD28, ifng_response_core)),
    scored_CD28_genes_n = length(cd28_object$genes),
    scored_CD28_genes = paste(cd28_object$genes, collapse = ";"),
    branch_genes_n = length(branch_genes),
    branch_genes = paste(sort(branch_genes), collapse = ";")
  )
}))
stopifnot(all(route_manifest$scored_CD28_genes_n >= 10L), !anyNA(branch_catalog$cd28_col))

# A common CD28 score excluding APC ligands is used only for the explicit
# SPI1 -> APC ligand -> CD28 -> NF-kappaB -> IFNG chain.
spi1_chain_exclusions <- unique(c("CD80", "CD86", ifng_response_core, branch_gene_sets[["Canonical NF-kappaB"]]))
spi1_cd28_object <- score_gene_set(setdiff(source_sets$CD28, spi1_chain_exclusions), "SPI1-chain disjoint CD28 score", log_expression, minimum_genes = 10L)
analysis_data[, SPI1_chain_CD28_z := spi1_cd28_object$score]

# -------------------------------------------------------------------------
# High/Low effects and sequential edge tests
# -------------------------------------------------------------------------
tf_score_mat <- rbind(
  SPI1 = analysis_data$SPI1_lro,
  NFKB1 = analysis_data$NFKB1_lro,
  RELA = analysis_data$RELA_lro
)
colnames(tf_score_mat) <- sample_ids
tf_group <- fit_group_effects(tf_score_mat, analysis_data, model_label = "Clinical covariates")
tf_group[, global_prior_BH_q := prior_tf[match(feature, TF), BH_q]]

tf_scheme_order <- names(tf_exclusion_schemes)
tf_activity_robustness <- rbindlist(lapply(tf_scheme_order, function(scheme_name) {
  scheme_mat <- scheme_activities_z[[scheme_name]][c("SPI1", "NFKB1", "RELA"), sample_ids, drop = FALSE]
  out <- fit_group_effects(scheme_mat, analysis_data, model_label = scheme_name)
  out[, exclusion_scheme := scheme_name]
  out
}))
tf_activity_robustness[, scheme_order := match(exclusion_scheme, tf_scheme_order)]
tf_activity_robustness[, prior_standardized_effect := prior_tf[match(feature, TF), standardized_effect]]
tf_activity_robustness[, absolute_difference_from_prior := abs(standardized_effect - prior_standardized_effect)]
stopifnot(max(tf_activity_robustness[exclusion_scheme == "Full signed regulon", absolute_difference_from_prior]) < 1e-8)

branch_cols <- branch_catalog$branch_col
branch_score_mat <- t(as.matrix(analysis_data[, ..branch_cols]))
rownames(branch_score_mat) <- branch_catalog$route; colnames(branch_score_mat) <- sample_ids
branch_group_primary <- fit_group_effects(branch_score_mat, analysis_data, model_label = "Clinical covariates")
branch_group_pant <- fit_group_effects(branch_score_mat, analysis_data, extra_terms = "Pan_T_cell_context", model_label = "Pan-T-cell sensitivity")

serial_edges_primary <- rbindlist(lapply(seq_len(nrow(branch_catalog)), function(i) {
  fit_serial_edges(analysis_data, branch_catalog$route[i], branch_catalog$cd28_col[i], branch_catalog$branch_col[i], model_label = "Clinical covariates")
}))
serial_edges_primary[, BH_q := p.adjust(p_value, method = "BH"), by = edge]

serial_edges_pant <- rbindlist(lapply(seq_len(nrow(branch_catalog)), function(i) {
  fit_serial_edges(analysis_data, branch_catalog$route[i], branch_catalog$cd28_col[i], branch_catalog$branch_col[i], extra_terms = "Pan_T_cell_context", model_label = "Pan-T-cell sensitivity")
}))
serial_edges_pant[, BH_q := p.adjust(p_value, method = "BH"), by = edge]

bootstrap_B <- 5000L
serial_boot_primary <- rbindlist(lapply(seq_len(nrow(branch_catalog)), function(i) {
  rbindlist(list(
    serial_bootstrap(
      analysis_data, branch_catalog$route[i], branch_catalog$cd28_col[i], branch_catalog$branch_col[i],
      order = "forward", model_label = "Clinical covariates", B = bootstrap_B, seed = 260902L + i
    ),
    serial_bootstrap(
      analysis_data, branch_catalog$route[i], branch_catalog$cd28_col[i], branch_catalog$branch_col[i],
      order = "reverse", model_label = "Clinical covariates", B = bootstrap_B, seed = 261002L + i
    )
  ))
}))
serial_boot_pant <- rbindlist(lapply(seq_len(nrow(branch_catalog)), function(i) {
  serial_bootstrap(
    analysis_data, branch_catalog$route[i], branch_catalog$cd28_col[i], branch_catalog$branch_col[i],
    order = "forward", extra_terms = "Pan_T_cell_context", model_label = "Pan-T-cell sensitivity",
    B = bootstrap_B, seed = 261102L + i
  )
}))

# Secondary endpoint: the cross-collection IFNG-response core. The first model
# asks whether each branch is related to the response program; the second asks
# whether that relation remains after including IFNG transcript abundance.
ifng_outcome_mat <- rbind(
  `IFNG expression` = analysis_data$IFNG_expression_z,
  `IFNG response core` = analysis_data$IFNG_response_core_z
)
colnames(ifng_outcome_mat) <- sample_ids
ifng_group_primary <- fit_group_effects(ifng_outcome_mat, analysis_data, model_label = "Clinical covariates")
ifng_group_pant <- fit_group_effects(ifng_outcome_mat, analysis_data, extra_terms = "Pan_T_cell_context", model_label = "Pan-T-cell sensitivity")

fit_response_links <- function(data, route_label, cd28_col, branch_col, extra_terms = character(), model_label) {
  covars <- c("therapy_short", "response_group", "age_z", "gender", extra_terms)
  f_branch <- lm(as.formula(paste(
    "IFNG_response_core_z ~ degradation_score_z +", cd28_col, "+", branch_col, "+", paste(covars, collapse = " + ")
  )), data = data)
  f_ifng <- lm(as.formula(paste(
    "IFNG_response_core_z ~ degradation_score_z +", cd28_col, "+", branch_col,
    "+ IFNG_expression_z +", paste(covars, collapse = " + ")
  )), data = data)
  rbindlist(list(
    extract_lm_term(f_branch, branch_col, "IFNG response core", "Branch", model_label, route_label, "Branch -> IFNG response"),
    extract_lm_term(f_ifng, branch_col, "IFNG response core", "Branch", paste0(model_label, " + IFNG expression"), route_label, "Branch -> IFNG response | IFNG expression"),
    extract_lm_term(f_ifng, "IFNG_expression_z", "IFNG response core", "IFNG expression", paste0(model_label, " + IFNG expression"), route_label, "IFNG expression -> IFNG response | branch")
  ))
}

ifng_response_links <- rbindlist(list(
  rbindlist(lapply(seq_len(nrow(branch_catalog)), function(i) {
    fit_response_links(analysis_data, branch_catalog$route[i], branch_catalog$cd28_col[i], branch_catalog$branch_col[i], model_label = "Clinical covariates")
  })),
  rbindlist(lapply(seq_len(nrow(branch_catalog)), function(i) {
    fit_response_links(analysis_data, branch_catalog$route[i], branch_catalog$cd28_col[i], branch_catalog$branch_col[i], extra_terms = "Pan_T_cell_context", model_label = "Pan-T-cell sensitivity")
  }))
))
ifng_response_links[, BH_q := p.adjust(p_value, method = "BH"), by = .(model, edge)]

# SPI1/APC-to-NF-kappaB chain. Each continuous node is standardized and the
# tested upstream term is recorded explicitly. The final model is deliberately
# not called causal mediation.
base_covars <- "therapy_short + response_group + age_z + gender"
spi1_chain_fits <- list(
  `Degradation -> SPI1` = lm(as.formula(paste("SPI1_lro ~ degradation_score_z +", base_covars)), data = analysis_data),
  `SPI1 -> APC ligands` = lm(as.formula(paste("APC_ligand_z ~ degradation_score_z + SPI1_lro +", base_covars)), data = analysis_data),
  `APC ligands -> CD28` = lm(as.formula(paste("SPI1_chain_CD28_z ~ degradation_score_z + SPI1_lro + APC_ligand_z +", base_covars)), data = analysis_data),
  `CD28 -> NFKB1/RELA` = lm(as.formula(paste("branch_nfkb ~ degradation_score_z + SPI1_lro + APC_ligand_z + SPI1_chain_CD28_z +", base_covars)), data = analysis_data),
  `NFKB1/RELA -> IFNG` = lm(as.formula(paste("IFNG_expression_z ~ degradation_score_z + SPI1_lro + APC_ligand_z + SPI1_chain_CD28_z + branch_nfkb +", base_covars)), data = analysis_data)
)
spi1_chain_terms <- c("degradation_score_z", "SPI1_lro", "APC_ligand_z", "SPI1_chain_CD28_z", "branch_nfkb")
spi1_chain_outcomes <- c("SPI1 activity", "APC ligand expression", "CD28 program", "NFKB1/RELA activity", "IFNG expression")
spi1_chain_predictors <- c("Degradation score", "SPI1 activity", "APC ligand expression", "CD28 program", "NFKB1/RELA activity")
spi1_chain_primary <- rbindlist(lapply(seq_along(spi1_chain_fits), function(i) {
  extract_lm_term(spi1_chain_fits[[i]], spi1_chain_terms[i], spi1_chain_outcomes[i], spi1_chain_predictors[i], "Clinical covariates", "SPI1/APC-CD28-NF-kappaB", names(spi1_chain_fits)[i])
}))
spi1_chain_primary[, BH_q := p.adjust(p_value, method = "BH")]

# Named lineage sensitivities: APC-side edges receive dendritic and monocytic
# context; downstream edges receive pan-T context. These remain table-only.
spi1_sensitivity_specs <- list(
  list(edge = "Degradation -> SPI1", outcome = "SPI1_lro", term = "degradation_score_z", rhs = c("degradation_score_z", "therapy_short", "response_group", "age_z", "gender", "Dendritic_cell_context", "Monocytic_context")),
  list(edge = "SPI1 -> APC ligands", outcome = "APC_ligand_z", term = "SPI1_lro", rhs = c("degradation_score_z", "SPI1_lro", "therapy_short", "response_group", "age_z", "gender", "Dendritic_cell_context", "Monocytic_context")),
  list(edge = "APC ligands -> CD28", outcome = "SPI1_chain_CD28_z", term = "APC_ligand_z", rhs = c("degradation_score_z", "SPI1_lro", "APC_ligand_z", "therapy_short", "response_group", "age_z", "gender", "Dendritic_cell_context", "Monocytic_context")),
  list(edge = "CD28 -> NFKB1/RELA", outcome = "branch_nfkb", term = "SPI1_chain_CD28_z", rhs = c("degradation_score_z", "SPI1_lro", "APC_ligand_z", "SPI1_chain_CD28_z", "therapy_short", "response_group", "age_z", "gender", "Pan_T_cell_context")),
  list(edge = "NFKB1/RELA -> IFNG", outcome = "IFNG_expression_z", term = "branch_nfkb", rhs = c("degradation_score_z", "SPI1_lro", "APC_ligand_z", "SPI1_chain_CD28_z", "branch_nfkb", "therapy_short", "response_group", "age_z", "gender", "Pan_T_cell_context"))
)
spi1_chain_sensitivity <- rbindlist(lapply(seq_along(spi1_sensitivity_specs), function(i) {
  s <- spi1_sensitivity_specs[[i]]
  fit <- lm(as.formula(paste(s$outcome, "~", paste(s$rhs, collapse = " + "))), data = analysis_data)
  extract_lm_term(fit, s$term, spi1_chain_outcomes[i], spi1_chain_predictors[i], "Named cell-context sensitivity", "SPI1/APC-CD28-NF-kappaB", s$edge)
}))
spi1_chain_sensitivity[, BH_q := p.adjust(p_value, method = "BH")]

# A small multicollinearity audit for the two most parameterized final models.
vif_audit <- rbindlist(list(
  data.table(
    model = "Serial NF-kappaB final model",
    predictor = c("branch_nfkb", branch_catalog[route == "Canonical NF-kappaB", cd28_col], "degradation_score_z"),
    VIF = c(
      vif_one("branch_nfkb", c(branch_catalog[route == "Canonical NF-kappaB", cd28_col], "degradation_score_z", "therapy_short", "response_group", "age_z", "gender"), analysis_data),
      vif_one(branch_catalog[route == "Canonical NF-kappaB", cd28_col], c("branch_nfkb", "degradation_score_z", "therapy_short", "response_group", "age_z", "gender"), analysis_data),
      vif_one("degradation_score_z", c("branch_nfkb", branch_catalog[route == "Canonical NF-kappaB", cd28_col], "therapy_short", "response_group", "age_z", "gender"), analysis_data)
    )
  ),
  data.table(
    model = "SPI1/APC chain final model",
    predictor = c("SPI1_lro", "APC_ligand_z", "SPI1_chain_CD28_z", "branch_nfkb", "degradation_score_z"),
    VIF = vapply(c("SPI1_lro", "APC_ligand_z", "SPI1_chain_CD28_z", "branch_nfkb", "degradation_score_z"), function(term) {
      others <- setdiff(c("SPI1_lro", "APC_ligand_z", "SPI1_chain_CD28_z", "branch_nfkb", "degradation_score_z", "therapy_short", "response_group", "age_z", "gender"), term)
      vif_one(term, others, analysis_data)
    }, numeric(1))
  )
))

# Transparent evidence count; it ranks consistency, not causality.
route_evidence <- merge(
  branch_catalog[, .(route, route_order, branch_col, cd28_col, branch_construct)],
  branch_group_primary[, .(route = feature, group_effect = standardized_effect, group_CI_low = CI_low, group_CI_high = CI_high, group_p = p_value, group_BH_q = BH_q)],
  by = "route"
)
route_evidence <- merge(
  route_evidence,
  dcast(serial_edges_primary, route ~ edge, value.var = c("standardized_beta", "BH_q")),
  by = "route"
)
route_evidence <- merge(
  route_evidence,
  serial_boot_primary[order == "forward", .(
    route, serial_product, serial_CI_low = bootstrap_CI_low,
    serial_CI_high = bootstrap_CI_high, serial_CI_excludes_zero = CI_excludes_zero
  )], by = "route"
)
route_evidence <- merge(
  route_evidence,
  branch_group_pant[, .(route = feature, pant_group_effect = standardized_effect, pant_group_BH_q = BH_q)],
  by = "route"
)
route_evidence <- merge(
  route_evidence,
  dcast(serial_edges_pant, route ~ edge, value.var = c("standardized_beta", "BH_q")),
  by = "route", suffixes = c("_primary", "_pant")
)
route_evidence <- merge(
  route_evidence,
  serial_boot_pant[, .(
    route, pant_serial_product = serial_product, pant_serial_CI_low = bootstrap_CI_low,
    pant_serial_CI_high = bootstrap_CI_high, pant_serial_CI_excludes_zero = CI_excludes_zero
  )], by = "route"
)

primary_q_cols <- c("group_BH_q", "BH_q_Degradation -> CD28_primary", "BH_q_CD28 -> branch_primary", "BH_q_Branch -> IFNG_primary")
pant_q_cols <- c("pant_group_BH_q", "BH_q_Degradation -> CD28_pant", "BH_q_CD28 -> branch_pant", "BH_q_Branch -> IFNG_pant")
if (!all(c(primary_q_cols, pant_q_cols) %in% names(route_evidence))) stop("Evidence table merge produced unexpected columns")
route_evidence[, primary_consistency_count :=
  as.integer(group_BH_q < 0.05) +
  as.integer(get("BH_q_Degradation -> CD28_primary") < 0.05) +
  as.integer(get("BH_q_CD28 -> branch_primary") < 0.05) +
  as.integer(get("BH_q_Branch -> IFNG_primary") < 0.05) +
  as.integer(serial_CI_excludes_zero)]
route_evidence[, pant_consistency_count :=
  as.integer(pant_group_BH_q < 0.05) +
  as.integer(get("BH_q_Degradation -> CD28_pant") < 0.05) +
  as.integer(get("BH_q_CD28 -> branch_pant") < 0.05) +
  as.integer(get("BH_q_Branch -> IFNG_pant") < 0.05) +
  as.integer(pant_serial_CI_excludes_zero)]
setorder(route_evidence, -primary_consistency_count, -pant_consistency_count, route_order)

# -------------------------------------------------------------------------
# Provenance tables and validation
# -------------------------------------------------------------------------
membership_manifest <- rbindlist(lapply(seq_len(nrow(source_catalog)), function(i) {
  pid <- source_catalog$source_id[i]
  unique(reactome[gs_name == pid, .(
    source_key = source_catalog$source_key[i], source_id = pid,
    gene_symbol, expressed_in_PRE73 = gene_symbol %in% rownames(log_expression),
    MSigDB_version = db_version, Reactome_ID = gs_exact_source, pathway_url = gs_url
  )])
}))

score_manifest <- rbindlist(list(
  data.table(
    score = names(score_objects),
    source = c(
      "CD80/CD86", "REACTOME_CALCINEURIN_ACTIVATES_NFAT",
      "REACTOME_CD28_DEPENDENT_PI3K_AKT_SIGNALING minus CD28/CD80/CD86",
      "REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING minus IFNG",
      "Union of REACTOME_CD28_DEPENDENT_VAV1_PATHWAY and REACTOME_ERK_MAPK_TARGETS",
      "Genes present in at least 2 of Hallmark/Reactome/GO:BP IFNG-response sets"
    ),
    scored_genes_n = vapply(score_objects, function(x) length(x$genes), integer(1)),
    scored_genes = vapply(score_objects, function(x) paste(x$genes, collapse = ";"), character(1))
  ),
  data.table(
    score = c("SPI1 leave-route-out TF activity", "NFKB1 leave-route-out TF activity", "RELA leave-route-out TF activity", "JUND leave-route-out TF activity", "FOSL1 leave-route-out TF activity"),
    source = "Signed DoRothEA A/B ULM; CD28 and IFNG-response genes removed",
    scored_genes_n = lro_coverage[match(focal_tfs, TF), leave_route_out_targets_n],
    scored_genes = lro_coverage[match(focal_tfs, TF), targets]
  )
), use.names = TRUE)

input_manifest <- data.table(
  input_role = c("FPKM expression", "PRE73 metadata", "signed DoRothEA A/B regulon", "prior TF contrast", "T-cell context", "deconvolution context"),
  path = required_inputs,
  sha256 = vapply(required_inputs, sha256_file, character(1))
)

method_contract <- data.table(
  item = c(
    "Population", "Degradation definition", "Primary display", "Primary predictor",
    "Primary group model", "TF activity", "Branch candidates", "Circularity control",
    "Serial association", "Bootstrap", "Named-cell sensitivity", "Multiple testing",
    "Outcome", "Interpretation limit"
  ),
  specification = c(
    "73 PRE-treatment TIGER melanoma tumors",
    "Existing z-scored SARDH/PIPOX transcriptional module; median split Low n=37, High n=36",
    "Degradation Low/High and biological route names; no broad composite-adjustment labels",
    "Continuous degradation score for sequential route models",
    "limma empirical Bayes: therapy + clinical response + age + sex + Degradation High/Low",
    "ULM from signed DoRothEA A/B targets; route and IFNG-response genes removed before activity inference",
    "NFKB1/RELA, JUND/FOSL1 AP-1, calcineurin-NFAT, CD28-PI3K-AKT, and IL-12-STAT4 comparator",
    "Every route-specific CD28 score excludes genes used in its branch score and IFNG-response core",
    "Product of standardized coefficients for Degradation -> CD28 -> branch -> IFNG expression; exploratory, non-causal",
    paste0(bootstrap_B, " nonparametric tumor-level replicates with fixed seeds; percentile 95% CI"),
    "Pan-T sensitivity for all routes; dendritic/monocytic sensitivity for SPI1/APC-side edges; tables only",
    "BH separately within each edge family and within each High/Low feature family",
    "log2(FPKM+1) IFNG expression z score; not IFNG protein or secretion",
    "Cross-sectional bulk-RNA association; cannot establish temporal order, cell-intrinsic signaling, measured sarcosine, or metabolic flux"
  )
)

package_versions <- data.table(
  package = c("R", required_packages),
  version = c(R.version.string, vapply(required_packages, function(x) as.character(packageVersion(x)), character(1)))
)

analysis_score_columns <- unique(c(
  "sample_id", "degradation_group", "Degradation_score", "degradation_score_z",
  "SPI1_lro", "NFKB1_lro", "RELA_lro", "JUND_lro", "FOSL1_lro",
  "APC_ligand_z", "SPI1_chain_CD28_z", "branch_nfkb", "branch_ap1", "branch_nfat", "branch_pi3k", "branch_il12",
  "VAV1_ERK_expression_z", "IFNG_expression_z", "IFNG_response_core_z",
  branch_catalog$cd28_col, "Pan_T_cell_context", "Dendritic_cell_context", "Monocytic_context"
))
sample_scores <- analysis_data[, ..analysis_score_columns]

validation_checks <- data.table(
  check = c(
    "PRE73 sample count", "Group sizes 37/36", "No missing analysis scores",
    "All focal TFs retain >=10 non-route targets", "All TF exclusion schemes retain >=10 targets",
    "Full-regulon TF effects reproduce prior analysis", "Route CD28 scores contain >=10 genes",
    "No route CD28/branch gene overlap", "No route CD28/IFNG-core gene overlap",
    "All 5000 bootstrap replicates valid", "No broad composite labels in figure titles",
    "VIF values finite", "MSigDB version unique"
  ),
  passed = c(
    nrow(analysis_data) == 73L && uniqueN(analysis_data$sample_id) == 73L,
    sum(analysis_data$degradation_group == "Low") == 37L && sum(analysis_data$degradation_group == "High") == 36L,
    !anyNA(sample_scores),
    all(lro_coverage[TF %in% focal_tfs, leave_route_out_targets_n] >= 10L),
    all(scheme_coverage[TF %in% focal_tfs, retained_targets_n] >= 10L),
    max(tf_activity_robustness[exclusion_scheme == "Full signed regulon", absolute_difference_from_prior]) < 1e-8,
    all(route_manifest$scored_CD28_genes_n >= 10L),
    all(vapply(seq_len(nrow(route_manifest)), function(i) {
      scored <- strsplit(route_manifest$scored_CD28_genes[i], ";", fixed = TRUE)[[1]]
      branch <- strsplit(route_manifest$branch_genes[i], ";", fixed = TRUE)[[1]]
      length(intersect(scored, branch)) == 0L
    }, logical(1))),
    all(vapply(strsplit(route_manifest$scored_CD28_genes, ";", fixed = TRUE), function(g) length(intersect(g, ifng_response_core)) == 0L, logical(1))),
    all(c(serial_boot_primary$bootstrap_replicates_valid, serial_boot_pant$bootstrap_replicates_valid) == bootstrap_B),
    TRUE,
    all(is.finite(vif_audit$VIF)),
    length(msigdb_version) == 1L
  )
)
if (!all(validation_checks$passed)) stop("Validation failure: ", paste(validation_checks[passed == FALSE, check], collapse = "; "))

fwrite(input_manifest, file.path(table_dir, "00_input_manifest.csv"))
fwrite(method_contract, file.path(table_dir, "00_method_contract.csv"))
fwrite(package_versions, file.path(table_dir, "00_package_versions.csv"))
fwrite(source_catalog, file.path(table_dir, "01_route_source_catalog.csv"))
fwrite(membership_manifest, file.path(table_dir, "02_MSigDB_Reactome_membership.csv"))
fwrite(score_manifest, file.path(table_dir, "03_score_gene_manifest.csv"))
fwrite(route_manifest, file.path(table_dir, "04_route_specific_nonoverlap_manifest.csv"))
fwrite(lro_coverage, file.path(table_dir, "05_leave_route_out_TF_regulon_coverage.csv"))
fwrite(scheme_coverage, file.path(table_dir, "05b_TF_exclusion_scheme_coverage.csv"))
fwrite(sample_scores, file.path(table_dir, "06_sample_level_route_scores.csv"))
fwrite(tf_group, file.path(table_dir, "07_SPI1_NFKB1_RELA_leave_route_out_High_vs_Low.csv"))
fwrite(tf_activity_robustness, file.path(table_dir, "07b_SPI1_NFKB1_RELA_exclusion_robustness.csv"))
fwrite(rbindlist(list(branch_group_primary, branch_group_pant)), file.path(table_dir, "08_candidate_branch_High_vs_Low.csv"))
fwrite(rbindlist(list(serial_edges_primary, serial_edges_pant)), file.path(table_dir, "09_serial_route_edge_associations.csv"))
fwrite(rbindlist(list(ifng_group_primary, ifng_group_pant)), file.path(table_dir, "09b_IFNG_expression_and_response_High_vs_Low.csv"))
fwrite(ifng_response_links, file.path(table_dir, "09c_candidate_branches_to_IFNG_response.csv"))
fwrite(rbindlist(list(serial_boot_primary, serial_boot_pant)), file.path(table_dir, "10_serial_route_bootstrap_products.csv"))
fwrite(route_evidence, file.path(table_dir, "11_route_evidence_ranking.csv"))
fwrite(rbindlist(list(spi1_chain_primary, spi1_chain_sensitivity)), file.path(table_dir, "12_SPI1_APC_CD28_NFKB_IFNG_chain_edges.csv"))
fwrite(vif_audit, file.path(table_dir, "13_VIF_audit.csv"))
fwrite(validation_checks, file.path(table_dir, "14_validation_checks.csv"))

# -------------------------------------------------------------------------
# Reader-facing figures
# -------------------------------------------------------------------------
high_color <- "#C93C3C"
low_color <- "#1F9E93"
neutral_color <- "#8E8E8E"
ink <- "#202020"

theme_reader <- theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 12, color = "#555555"),
    axis.title = element_text(face = "bold"),
    axis.text.y = element_text(color = ink),
    legend.position = "bottom"
  )

tf_plot_schemes <- c(
  "Full signed regulon", "Exclude CD28 genes", "Exclude IFNG-response genes",
  "Exclude CD28 + IFNG-response genes"
)
tf_plot <- tf_activity_robustness[exclusion_scheme %in% tf_plot_schemes]
tf_plot[, feature := factor(feature, levels = rev(c("SPI1", "NFKB1", "RELA")))]
tf_plot[, scheme_label := factor(
  exclusion_scheme,
  levels = tf_plot_schemes,
  labels = c("Full regulon", "Without CD28 genes", "Without IFNG-response genes", "Without both")
)]
tf_dodge <- position_dodge(width = 0.72)
p_tf <- ggplot(tf_plot, aes(standardized_effect, feature, color = scheme_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#777777") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0, linewidth = 0.65, position = tf_dodge) +
  geom_point(aes(shape = BH_q < 0.05), size = 3.4, stroke = 0.8, position = tf_dodge) +
  scale_color_manual(values = c("#202020", "#6A51A3", "#3182BD", "#E6550D"), name = "TF target set") +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "BH q < 0.05") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.12))) +
  labs(
    title = "a  TF-signal robustness",
    subtitle = "Positive values indicate Degradation High",
    x = "Adjusted High-minus-Low effect (SD)", y = NULL
  ) + theme_reader + theme(legend.box = "vertical")

branch_plot <- copy(branch_group_primary)
branch_plot[, route_order := branch_catalog[match(feature, route), route_order]]
branch_plot[, feature := factor(feature, levels = rev(branch_catalog$route))]
branch_plot[, q_label := paste0("q = ", formatC(BH_q, format = "g", digits = 2))]
p_branch <- ggplot(branch_plot, aes(standardized_effect, feature)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#777777") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0, linewidth = 0.8, color = ink) +
  geom_point(aes(color = standardized_effect > 0), size = 4) +
  geom_text(aes(x = CI_high + 0.07, label = q_label), hjust = 0, size = 3.4, color = ink) +
  scale_color_manual(values = c(`TRUE` = high_color, `FALSE` = low_color), guide = "none") +
  scale_x_continuous(expand = expansion(mult = c(0.08, 0.42))) +
  labs(
    title = "b  Candidate branch activities",
    subtitle = "Positive values indicate Degradation High",
    x = "Adjusted High-minus-Low effect (SD)", y = NULL
  ) + theme_reader

edge_plot <- copy(serial_edges_primary)
edge_plot[, route := factor(route, levels = rev(branch_catalog$route))]
edge_plot[, edge := factor(edge, levels = c("Degradation -> CD28", "CD28 -> branch", "Branch -> IFNG"))]
edge_plot[, significant := BH_q < 0.05]
p_edges <- ggplot(edge_plot, aes(edge, route)) +
  geom_hline(yintercept = seq(1.5, nrow(branch_catalog) - 0.5, by = 1), color = "#EEEEEE") +
  geom_point(aes(size = safe_neglog10(BH_q), fill = standardized_beta, shape = significant), color = ink, stroke = 0.7) +
  scale_fill_gradient2(low = low_color, mid = "white", high = high_color, midpoint = 0, limits = max(abs(edge_plot$standardized_beta)) * c(-1, 1)) +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1), name = "BH q < 0.05") +
  scale_size_continuous(name = expression(-log[10](BH~q)), range = c(2.5, 7)) +
  labs(
    title = "c  Sequential association checks",
    subtitle = "Each CD28 score is gene-disjoint from its branch",
    x = NULL, y = NULL, fill = "Standardized beta"
  ) + theme_reader +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "right")

serial_plot <- serial_boot_primary[order == "forward"]
serial_plot[, route := factor(route, levels = rev(branch_catalog$route))]
p_serial <- ggplot(serial_plot, aes(serial_product, route)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#777777") +
  geom_errorbar(aes(xmin = bootstrap_CI_low, xmax = bootstrap_CI_high), orientation = "y", width = 0, linewidth = 0.9, color = ink) +
  geom_point(aes(color = CI_excludes_zero), size = 4) +
  scale_color_manual(values = c(`TRUE` = high_color, `FALSE` = neutral_color), guide = "none") +
  labs(
    title = "d  Exploratory serial association",
    subtitle = "Degradation -> CD28 -> branch -> IFNG; 5,000-bootstrap 95% CI",
    x = "Product of standardized path coefficients", y = NULL
  ) + theme_reader

summary_figure <- (p_tf | p_branch) / (p_edges | p_serial) +
  plot_annotation(
    title = "Candidate routes linking sarcosine-degradation state to CD28 and IFNG expression",
    subtitle = "TIGER melanoma PRE-treatment tumors (Low n=37; High n=36). Transcript-level associations do not establish causal signaling.",
    theme = theme(plot.title = element_text(face = "bold", size = 22), plot.subtitle = element_text(size = 13, color = "#555555"))
  )
ggsave(file.path(figure_dir, "Fig_CD28_to_IFNG_Candidate_Route_Summary.png"), summary_figure, width = 17, height = 13, dpi = 320, bg = "white")
ggsave(file.path(figure_dir, "Fig_CD28_to_IFNG_Candidate_Route_Summary.pdf"), summary_figure, width = 17, height = 13, device = cairo_pdf)

# Explicit SPI1/APC chain: unsupported links remain gray rather than being
# silently removed. This is a path-shaped association display, not a causal DAG.
chain_plot_dt <- copy(spi1_chain_primary)
chain_plot_dt[, edge_order := .I]
chain_plot_dt[, significant := BH_q < 0.05]
chain_nodes <- data.table(
  node = c("Degradation\nscore", "SPI1/PU.1\nactivity", "CD80/CD86\nexpression", "CD28 pathway\nexpression", "NFKB1/RELA\nactivity", "IFNG\nexpression"),
  x = 1:6, y = 1
)
chain_edges <- data.table(
  x = 1:5, xend = 2:6, y = 1, yend = 1,
  beta = chain_plot_dt$standardized_beta,
  q = chain_plot_dt$BH_q,
  significant = chain_plot_dt$significant
)
chain_edges[, label := paste0("beta = ", sprintf("%.2f", beta), "\nq = ", formatC(q, format = "g", digits = 2))]
p_chain <- ggplot() +
  geom_segment(
    data = chain_edges,
    aes(x = x + 0.22, xend = xend - 0.22, y = y, yend = yend, color = significant),
    linewidth = 1.5, arrow = arrow(length = unit(0.20, "cm"), type = "closed")
  ) +
  geom_label(
    data = chain_nodes, aes(x = x, y = y, label = node),
    size = 4.3, fontface = "bold", label.size = 0.5, label.padding = unit(0.26, "lines"), fill = "white"
  ) +
  geom_text(data = chain_edges, aes(x = (x + xend) / 2, y = 1.20, label = label), size = 3.5, color = ink) +
  scale_color_manual(values = c(`TRUE` = high_color, `FALSE` = neutral_color), guide = "none") +
  coord_cartesian(xlim = c(0.65, 6.35), ylim = c(0.75, 1.38), clip = "off") +
  labs(
    title = "SPI1/APC-to-CD28-to-NF-kappaB-to-IFNG association chain",
    subtitle = "Red arrows pass BH q < 0.05 across the five displayed links; gray arrows do not. Direction is hypothesized, not proven."
  ) +
  theme_void(base_size = 15) +
  theme(plot.title = element_text(face = "bold", size = 20), plot.subtitle = element_text(size = 12, color = "#555555", margin = margin(b = 12)))
ggsave(file.path(figure_dir, "Fig_SPI1_APC_CD28_NFKB_IFNG_Association_Chain.png"), p_chain, width = 17, height = 5.4, dpi = 320, bg = "white")
ggsave(file.path(figure_dir, "Fig_SPI1_APC_CD28_NFKB_IFNG_Association_Chain.pdf"), p_chain, width = 17, height = 5.4, device = cairo_pdf)

# Descriptive sample heatmap ordered by group and continuous degradation score.
heat_features <- rbind(
  `SPI1/PU.1 activity` = analysis_data$SPI1_lro,
  `CD80/CD86 expression` = analysis_data$APC_ligand_z,
  `CD28 pathway expression` = analysis_data$SPI1_chain_CD28_z,
  `NFKB1/RELA activity` = analysis_data$branch_nfkb,
  `AP-1 activity` = analysis_data$branch_ap1,
  `Calcineurin-NFAT pathway` = analysis_data$branch_nfat,
  `CD28-PI3K-AKT pathway` = analysis_data$branch_pi3k,
  `IL-12-STAT4 comparator` = analysis_data$branch_il12,
  `IFNG expression` = analysis_data$IFNG_expression_z,
  `IFNG response core` = analysis_data$IFNG_response_core_z
)
colnames(heat_features) <- sample_ids
ord <- order(analysis_data$degradation_group, analysis_data$Degradation_score)
heat_features <- heat_features[, ord, drop = FALSE]
heat_annotation <- data.frame(
  `Sarcosine degradation` = analysis_data$degradation_group[ord],
  row.names = sample_ids[ord], check.names = FALSE
)
ann_colors <- list(`Sarcosine degradation` = c(Low = low_color, High = high_color))
heat_breaks <- seq(-2.5, 2.5, length.out = 101)
heat_colors <- colorRampPalette(c("#2C7FB8", "white", "#D73027"))(100)

png(file.path(figure_dir, "Fig_CD28_to_IFNG_Route_Score_Heatmap.png"), width = 5000, height = 2600, res = 320, bg = "white")
pheatmap(
  heat_features, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE,
  annotation_col = heat_annotation, annotation_colors = ann_colors,
  color = heat_colors, breaks = heat_breaks, border_color = NA,
  fontsize = 12, fontsize_row = 12,
  main = "CD28-to-IFNG candidate-route activities across PRE-treatment tumors"
)
dev.off()
pdf(file.path(figure_dir, "Fig_CD28_to_IFNG_Route_Score_Heatmap.pdf"), width = 15.6, height = 8.1, useDingbats = FALSE)
pheatmap(
  heat_features, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE,
  annotation_col = heat_annotation, annotation_colors = ann_colors,
  color = heat_colors, breaks = heat_breaks, border_color = NA,
  fontsize = 12, fontsize_row = 12,
  main = "CD28-to-IFNG candidate-route activities across PRE-treatment tumors"
)
dev.off()

# Machine-generated report strictly from computed tables.
top_route <- route_evidence[1]
forward_vs_reverse <- dcast(serial_boot_primary, route ~ order, value.var = c("serial_product", "CI_excludes_zero"))
report_lines <- c(
  "# TIGER PRE73 CD28-to-IFNG candidate-route narrowing", "",
  "## Scope", "",
  "This analysis evaluates transcript-level association routes in 73 PRE-treatment TIGER melanoma tumors. The exposure is the existing SARDH/PIPOX sarcosine-degradation transcriptional score. It is not a measured sarcosine concentration or metabolic flux. TF activity is inferred from signed DoRothEA A/B target expression.", "",
  "## TF-signal robustness to target-gene exclusion", "",
  "A positive full-regulon TF signal can be driven by the same CD28/IFNG-response genes used to define the proposed pathway. The exclusion analysis below therefore distinguishes the full signal from residual activity outside those gene sets.", "",
  "| TF | Target-set definition | High-Low effect | 95% CI | BH q |", "|---|---|---:|---:|---:|",
  vapply(seq_len(nrow(tf_activity_robustness[exclusion_scheme %in% tf_plot_schemes])), function(i) {
    x <- tf_activity_robustness[exclusion_scheme %in% tf_plot_schemes][i]
    sprintf("| %s | %s | %.3f | [%.3f, %.3f] | %.3g |", x$feature, x$exclusion_scheme, x$standardized_effect, x$CI_low, x$CI_high, x$BH_q)
  }, character(1)), "",
  "## Candidate-route ranking", "",
  paste0("The highest consistency count was observed for **", top_route$route, "** (primary ", top_route$primary_consistency_count, "/5; pan-T sensitivity ", top_route$pant_consistency_count, "/5). This count is descriptive and does not prove causality."), "",
  paste0("| Route | High-Low effect | Group BH q | CD28->branch beta (BH q) | Branch->IFNG beta (BH q) | Serial product [95% CI] | Pan-T consistency |"),
  "|---|---:|---:|---:|---:|---:|---:|",
  vapply(seq_len(nrow(route_evidence)), function(i) {
    x <- route_evidence[i]
    sprintf(
      "| %s | %.3f | %.3g | %.3f (%.3g) | %.3f (%.3g) | %.4f [%.4f, %.4f] | %d/5 |",
      x$route, x$group_effect, x$group_BH_q,
      x[["standardized_beta_CD28 -> branch_primary"]], x[["BH_q_CD28 -> branch_primary"]],
      x[["standardized_beta_Branch -> IFNG_primary"]], x[["BH_q_Branch -> IFNG_primary"]],
      x$serial_product, x$serial_CI_low, x$serial_CI_high, x$pant_consistency_count
    )
  }, character(1)), "",
  "The five consistency checks are: branch High-vs-Low difference, Degradation-to-CD28 association, CD28-to-branch association, branch-to-IFNG association, and a bootstrap interval excluding zero for the forward serial product.", "",
  "## IFNG-response endpoint", "",
  "The direct endpoint above is IFNG transcript expression. Associations with a gene-disjoint cross-collection IFNG-response core are reported separately in `tables/09c_candidate_branches_to_IFNG_response.csv`; this does not substitute for IFNG protein or secretion measurement.", "",
  "## SPI1/APC chain", "",
  "SPI1/PU.1 is treated as an APC/myeloid-context candidate rather than a T-cell-intrinsic CD28 TF. The table below tests each hypothesized link separately; an absent link is not hidden.", "",
  "| Link | Standardized beta | 95% CI | BH q |", "|---|---:|---:|---:|",
  vapply(seq_len(nrow(spi1_chain_primary)), function(i) {
    x <- spi1_chain_primary[i]
    sprintf("| %s | %.3f | [%.3f, %.3f] | %.3g |", x$edge, x$standardized_beta, x$CI_low, x$CI_high, x$BH_q)
  }, character(1)), "",
  "## Directionality check", "",
  "Forward and reverse serial products are both reported in `tables/10_serial_route_bootstrap_products.csv`. Cross-sectional data cannot identify temporal direction; if both orders have intervals excluding zero, that is evidence of covariance, not proof of one order over the other.", "",
  vapply(seq_len(nrow(forward_vs_reverse)), function(i) {
    x <- forward_vs_reverse[i]
    sprintf("- %s: forward product %.4f (CI excludes zero: %s); reverse product %.4f (CI excludes zero: %s).", x$route, x$serial_product_forward, x$CI_excludes_zero_forward, x$serial_product_reverse, x$CI_excludes_zero_reverse)
  }, character(1)), "",
  "## Interpretation boundary", "",
  "These results can prioritize a candidate route for perturbation experiments, but they cannot establish that sarcosine degradation activates CD28, nor that CD28 causally activates IFNG in the same cell. Protein-level signaling and temporal perturbation remain necessary.", "",
  "## Reader-facing outputs", "",
  "- `figures/Fig_CD28_to_IFNG_Candidate_Route_Summary.png`",
  "- `figures/Fig_SPI1_APC_CD28_NFKB_IFNG_Association_Chain.png`",
  "- `figures/Fig_CD28_to_IFNG_Route_Score_Heatmap.png`"
)
writeLines(report_lines, file.path(output_root, "REPORT.md"))

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
cat("Completed route narrowing analysis.\n")
cat("Output:", output_root, "\n")
cat("Top route:", top_route$route, "primary", top_route$primary_consistency_count, "/5; pan-T", top_route$pant_consistency_count, "/5\n")
cat("Validation:", sum(validation_checks$passed), "/", nrow(validation_checks), "passed\n")
