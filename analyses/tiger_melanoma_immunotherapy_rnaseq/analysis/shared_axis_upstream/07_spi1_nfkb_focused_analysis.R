# TIGER melanoma PRE73: focused SPI1/PU.1 and canonical NF-kappaB analysis.
#
# Reader-facing figures show only the Degradation Low/High comparison and
# biologically named axes. TF activity is inferred from signed DoRothEA A/B
# targets; it is not protein abundance, phosphorylation, DNA binding, measured
# sarcosine, or metabolic flux.

options(stringsAsFactors = FALSE, width = 180)
set.seed(42)

analysis_root <- normalizePath(getwd())
if (basename(analysis_root) != "Melanoma-PRJEB23709") {
  stop("Run from the Melanoma-PRJEB23709 analysis root: ", analysis_root)
}

default_gsea_lib <- file.path(
  dirname(dirname(dirname(analysis_root))),
  "2024_Drug_Res_Updates_NSCLC", "RNA-seq공공데이터_GSE207422",
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
  library(data.table); library(limma); library(msigdbr)
  library(ggplot2); library(patchwork); library(pheatmap)
})

used_data_root <- normalizePath(file.path(analysis_root, "..", "사용데이터_모음"))
expression_path <- file.path(used_data_root, "Source_Input", "TIGER_PRJEB23709_expression_FPKM_gene_by_sample.csv")
metadata_path <- file.path(used_data_root, "Analysis_Ready", "TIGER_PRE73_Fig4bc_analysis_data.csv")
tf_root <- file.path(analysis_root, "results", "TIGER_PRE73_CD28_Upstream_TF_26.09.02")
tcell_root <- file.path(analysis_root, "results", "TIGER_PRE73_Tcell_Abundance_Adjusted_Function_26.09.02")
ifng_root <- file.path(analysis_root, "results", "TIGER_PRE73_IFNG_Axis_Focused_26.09.02")
deconv_root <- file.path(analysis_root, "results", "TIGER_PRE73_Degradation_HighLow_Deconvolution_7methods_26.09.02")

tf_activity_path <- file.path(tf_root, "tables", "09_sample_level_TF_ULM_activities_long.csv")
tf_primary_path <- file.path(tf_root, "tables", "07_TF_ULM_activity_High_vs_Low.csv")
regulon_path <- file.path(tf_root, "tables", "01_DoRothEA_AB_signed_regulon_all.csv")
direct_edge_path <- file.path(tf_root, "tables", "15_candidate_TF_to_CD28_pathway_edges.csv")
program_membership_path <- file.path(tcell_root, "tables", "02_MSigDB_program_gene_membership.csv")
program_score_path <- file.path(tcell_root, "tables", "03_sample_functional_program_scores_long.csv")
tcell_context_path <- file.path(tcell_root, "tables", "04_Tcell_abundance_adjusters_long.csv")
ifng_score_path <- file.path(ifng_root, "tables", "05_IFNG_axis_sample_level_scores.csv")
deconv_path <- file.path(deconv_root, "tables", "05_core_celltype_values_long.csv")
required_inputs <- c(
  expression_path, metadata_path, tf_activity_path, tf_primary_path, regulon_path,
  direct_edge_path, program_membership_path, program_score_path,
  tcell_context_path, ifng_score_path, deconv_path
)
for (path in required_inputs) if (!file.exists(path)) stop("Missing input: ", path)

output_root <- file.path(analysis_root, "results", "TIGER_PRE73_SPI1_NFKB1_RELA_Focused_26.09.02")
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
  if (anyNA(out) || !all(is.finite(out))) stop("Non-finite standardized vector")
  out
}
extract_lm_term <- function(fit, term, outcome, predictor, model_label) {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) stop("Term not found: ", term)
  df <- df.residual(fit)
  estimate <- unname(sm[term, "Estimate"])
  se <- unname(sm[term, "Std. Error"])
  t_value <- unname(sm[term, "t value"])
  p_value <- unname(sm[term, "Pr(>|t|)"])
  critical <- qt(0.975, df = df)
  data.table(
    outcome = outcome, predictor = predictor, model = model_label,
    standardized_beta = estimate, standard_error = se,
    CI_low = estimate - critical * se, CI_high = estimate + critical * se,
    t_value = t_value, residual_df = df, p_value = p_value,
    partial_R2 = t_value^2 / (t_value^2 + df)
  )
}
fit_group_limma <- function(outcome_mat, meta_table, adjuster = NULL, model_label) {
  if (is.null(adjuster)) {
    design <- model.matrix(
      ~ therapy_short + response_group + age_z + gender + degradation_group,
      data = meta_table
    )
  } else {
    stopifnot(adjuster %in% names(meta_table))
    design <- model.matrix(
      as.formula(paste0(
        "~ therapy_short + response_group + age_z + gender + ", adjuster,
        " + degradation_group"
      )), data = meta_table
    )
  }
  stopifnot(qr(design)$rank == ncol(design))
  fit <- eBayes(lmFit(outcome_mat[, meta_table$sample_id, drop = FALSE], design), robust = TRUE)
  coef_name <- "degradation_groupHigh"
  out <- data.table(
    feature = rownames(outcome_mat),
    standardized_effect = fit$coefficients[, coef_name],
    standard_error = fit$stdev.unscaled[, coef_name] * sqrt(fit$s2.post),
    moderated_t = fit$t[, coef_name], p_value = fit$p.value[, coef_name],
    residual_df = fit$df.total, model = model_label
  )
  out[, `:=`(
    CI_low = standardized_effect - qt(0.975, residual_df) * standard_error,
    CI_high = standardized_effect + qt(0.975, residual_df) * standard_error,
    BH_q_within_model = p.adjust(p_value, method = "BH")
  )]
  out[]
}

# Fixed PRE73 design and expression matrix.
expr_dt <- fread(expression_path, check.names = FALSE)
meta <- fread(metadata_path, check.names = FALSE)
tf_activity <- fread(tf_activity_path, check.names = FALSE)
tf_primary_prior <- fread(tf_primary_path, check.names = FALSE)
regulon <- fread(regulon_path, check.names = FALSE)
direct_edges <- fread(direct_edge_path, check.names = FALSE)
program_membership <- fread(program_membership_path, check.names = FALSE)
program_scores <- fread(program_score_path, check.names = FALSE)
tcell_context <- fread(tcell_context_path, check.names = FALSE)
ifng_scores <- fread(ifng_score_path, check.names = FALSE)
deconv <- fread(deconv_path, check.names = FALSE)

stopifnot(
  nrow(meta) == 73L, uniqueN(meta$sample_id) == 73L, all(meta$timepoint == "PRE"),
  all(c("sample_id", "therapy_short", "response_group", "age", "gender", "Degradation_score") %in% names(meta)),
  all(c("sample_id", "TF", "ULM_activity") %in% names(tf_activity)),
  all(c("TF", "target", "mor", "confidence", "expressed") %in% names(regulon))
)
degradation_median <- median(meta$Degradation_score)
meta[, degradation_group := factor(ifelse(Degradation_score > degradation_median, "High", "Low"), levels = c("Low", "High"))]
stopifnot(sum(meta$degradation_group == "Low") == 37L, sum(meta$degradation_group == "High") == 36L)
meta[, therapy_short := relevel(factor(therapy_short), ref = "antiPD1")]
meta[, response_group := relevel(factor(response_group), ref = "NR")]
meta[, gender := relevel(factor(gender), ref = "Male")]
meta[, `:=`(age_z = zscore(age), degradation_score_z = zscore(Degradation_score))]

gene_symbols <- trimws(expr_dt[[1]])
stopifnot(!anyNA(gene_symbols), all(nzchar(gene_symbols)), !anyDuplicated(gene_symbols))
sample_ids <- meta$sample_id
stopifnot(all(sample_ids %in% names(expr_dt)))
fpkm <- as.matrix(expr_dt[, ..sample_ids])
storage.mode(fpkm) <- "double"
rownames(fpkm) <- gene_symbols; colnames(fpkm) <- sample_ids
stopifnot(!anyNA(fpkm), all(is.finite(fpkm)), all(fpkm >= 0))
log_expression_all <- log2(fpkm + 1)
keep_expression <- rowSums(fpkm >= 1) >= ceiling(0.10 * ncol(fpkm))
log_expression <- log_expression_all[keep_expression, , drop = FALSE]
log_expression <- log_expression[apply(log_expression, 1L, var) > 0, , drop = FALSE]
stopifnot(ncol(log_expression) == 73L, nrow(log_expression) > 10000L)

# Three TF activities; NFKB1 and RELA are combined for downstream models.
focal_tfs <- c("SPI1", "NFKB1", "RELA")
focal_activity <- tf_activity[TF %in% focal_tfs]
stopifnot(nrow(focal_activity) == length(focal_tfs) * nrow(meta))
activity_wide <- dcast(focal_activity, sample_id ~ TF, value.var = "ULM_activity")
activity_wide <- activity_wide[match(sample_ids, sample_id)]
stopifnot(identical(activity_wide$sample_id, sample_ids), !anyNA(activity_wide[, ..focal_tfs]))
for (tf in focal_tfs) activity_wide[, paste0(tf, "_z") := zscore(get(tf))]
activity_wide[, NFKB_composite_z := zscore((NFKB1_z + RELA_z) / 2)]

analysis_data <- merge(meta, activity_wide, by = "sample_id", all.x = TRUE, sort = FALSE)
analysis_data <- analysis_data[match(sample_ids, sample_id)]
stopifnot(identical(analysis_data$sample_id, sample_ids))
activity_mat <- t(as.matrix(activity_wide[, .(SPI1_z, NFKB1_z, RELA_z, NFKB_composite_z)]))
colnames(activity_mat) <- sample_ids
rownames(activity_mat) <- c("SPI1", "NFKB1", "RELA", "NFKB1-RELA composite")
primary_focused <- fit_group_limma(activity_mat, analysis_data, model_label = "Clinical-covariate adjusted")
primary_focused[, BH_q_focused_family := p.adjust(p_value, method = "BH")]

prior_crosscheck <- merge(
  primary_focused[feature %in% focal_tfs, .(TF = feature, recomputed_standardized_effect = standardized_effect)],
  tf_primary_prior[TF %in% focal_tfs, .(
    TF, prior_standardized_effect = standardized_effect,
    prior_standardized_CI_low = standardized_CI_low,
    prior_standardized_CI_high = standardized_CI_high,
    prior_p_value = p_value, prior_global_BH_q = BH_q
  )],
  by = "TF", all = TRUE
)
prior_crosscheck[, absolute_difference := abs(recomputed_standardized_effect - prior_standardized_effect)]
stopifnot(nrow(prior_crosscheck) == 3L, max(prior_crosscheck$absolute_difference) < 1e-8)
primary_focused <- merge(
  primary_focused, prior_crosscheck[, .(
    feature = TF, prior_standardized_CI_low, prior_standardized_CI_high,
    prior_p_value, prior_global_BH_q
  )],
  by = "feature", all.x = TRUE, sort = FALSE
)
primary_focused[, feature_order := match(feature, rownames(activity_mat))]
setorder(primary_focused, feature_order)
primary_focused[, `:=`(
  display_CI_low = fifelse(!is.na(prior_standardized_CI_low), prior_standardized_CI_low, CI_low),
  display_CI_high = fifelse(!is.na(prior_standardized_CI_high), prior_standardized_CI_high, CI_high),
  display_p_value = fifelse(!is.na(prior_p_value), prior_p_value, p_value),
  display_BH_q = fifelse(!is.na(prior_global_BH_q), prior_global_BH_q, BH_q_focused_family)
)]

activity_cor <- cor(activity_wide[, .(SPI1_z, NFKB1_z, RELA_z, NFKB_composite_z)])
activity_cor_long <- as.data.table(as.table(activity_cor))
setnames(activity_cor_long, c("feature_1", "feature_2", "pearson_r"))

# Four named cell contexts; broad composite adjustment labels are not used.
context_cd8 <- tcell_context[abundance_method == "CD8 consensus", .(sample_id, CD8_T_cell_context = abundance_z)]
context_pant <- tcell_context[abundance_method == "MCPcounter Pan T", .(sample_id, Pan_T_cell_context = abundance_z)]
context_dc <- deconv[method == "MCPcounter" & analysis_feature == "Dendritic cells", .(sample_id, Dendritic_cell_context = analysis_value)]
context_mono <- deconv[method == "MCPcounter" & analysis_feature == "Monocytes", .(sample_id, Monocytic_context = analysis_value)]
for (x in list(context_cd8, context_pant, context_dc, context_mono)) stopifnot(nrow(x) == 73L, uniqueN(x$sample_id) == 73L)
context_wide <- Reduce(function(x, y) merge(x, y, by = "sample_id", all = TRUE), list(context_cd8, context_pant, context_dc, context_mono))
context_wide <- context_wide[match(sample_ids, sample_id)]
context_cols <- setdiff(names(context_wide), "sample_id")
for (v in context_cols) context_wide[, (v) := zscore(get(v))]
analysis_data <- merge(analysis_data, context_wide, by = "sample_id", all.x = TRUE, sort = FALSE)
analysis_data <- analysis_data[match(sample_ids, sample_id)]
stopifnot(!anyNA(analysis_data[, ..context_cols]))
context_labels <- c(
  CD8_T_cell_context = "CD8 T-cell context", Pan_T_cell_context = "Pan-T-cell context",
  Dendritic_cell_context = "Dendritic-cell context", Monocytic_context = "Monocytic context"
)

activity_context_results <- rbindlist(lapply(focal_tfs, function(tf) {
  outcome_col <- paste0(tf, "_z")
  rbindlist(lapply(context_cols, function(ctx) {
    fit <- lm(as.formula(paste0(
      outcome_col, " ~ ", ctx,
      " + degradation_score_z + therapy_short + response_group + age_z + gender"
    )), data = analysis_data)
    extract_lm_term(fit, ctx, tf, context_labels[[ctx]], "Clinical + degradation-score adjusted")
  }))
}))
activity_context_results[, BH_q := p.adjust(p_value, method = "BH")]

group_context_sensitivity <- rbindlist(lapply(context_cols, function(ctx) {
  fit_group_limma(activity_mat, analysis_data, adjuster = ctx, model_label = paste0("Adjusted for ", context_labels[[ctx]]))
}))
group_context_sensitivity[, BH_q_global := p.adjust(p_value, method = "BH")]

# Remove the union of the three expressed signed A/B regulons from outcomes.
candidate_regulon <- regulon[TF %in% focal_tfs & expressed == TRUE & confidence %in% c("A", "B")]
candidate_targets <- sort(unique(candidate_regulon$target))
stopifnot(length(candidate_targets) > 100L)
program_catalog <- data.table(
  program_id = c(
    "REACTOME_CO_STIMULATION_BY_CD28", "REACTOME_CD28_DEPENDENT_PI3K_AKT_SIGNALING",
    "GOBP_T_CELL_PROLIFERATION", "GOBP_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE",
    "GOBP_T_CELL_MEDIATED_CYTOTOXICITY"
  ),
  outcome_label = c(
    "CD28 co-stimulation", "CD28-PI3K/AKT signaling", "T-cell proliferation",
    "T-cell differentiation", "T-cell-mediated cytotoxicity"
  ),
  outcome_order = 1:5
)

hallmark_membership <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
reactome_membership <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))
gobp_membership <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"))
stopifnot(
  uniqueN(hallmark_membership$db_version) == 1L,
  identical(unique(hallmark_membership$db_version), unique(reactome_membership$db_version)),
  identical(unique(hallmark_membership$db_version), unique(gobp_membership$db_version))
)
msigdb_version <- unique(hallmark_membership$db_version)
response_membership <- rbindlist(list(
  hallmark_membership[gs_name == "HALLMARK_INTERFERON_GAMMA_RESPONSE", .(collection = "Hallmark", gene_symbol)],
  reactome_membership[gs_name == "REACTOME_INTERFERON_GAMMA_SIGNALING", .(collection = "Reactome", gene_symbol)],
  gobp_membership[gs_name == "GOBP_RESPONSE_TO_TYPE_II_INTERFERON", .(collection = "GO:BP", gene_symbol)]
))
ifng_core_genes <- setdiff(response_membership[, .N, by = gene_symbol][N >= 2L, gene_symbol], "IFNG")

score_gene_set <- function(genes, label, minimum_genes = 5L) {
  mapped <- intersect(unique(genes), rownames(log_expression))
  mapped <- mapped[apply(log_expression[mapped, , drop = FALSE], 1L, var) > 0]
  if (length(mapped) < minimum_genes) stop(label, " has only ", length(mapped), " usable genes")
  gene_z <- t(scale(t(log_expression[mapped, , drop = FALSE])))
  gene_z[!is.finite(gene_z)] <- 0
  score <- zscore(colMeans(gene_z)); names(score) <- colnames(log_expression)
  list(score = score, genes = mapped)
}

outcome_score_objects <- list(); outcome_manifest <- list()
for (i in seq_len(nrow(program_catalog))) {
  pid <- program_catalog$program_id[i]; label <- program_catalog$outcome_label[i]
  full_genes <- unique(program_membership[program_id == pid & expressed_in_PRE73 == TRUE, gene_symbol])
  nonoverlap_genes <- setdiff(full_genes, candidate_targets)
  score_object <- score_gene_set(nonoverlap_genes, label)
  outcome_score_objects[[label]] <- score_object
  outcome_manifest[[label]] <- data.table(
    outcome_label = label, source_id = pid, full_expressed_genes_n = length(full_genes),
    removed_candidate_regulon_targets_n = length(intersect(full_genes, candidate_targets)),
    scored_nonoverlap_genes_n = length(score_object$genes),
    removed_genes = paste(sort(intersect(full_genes, candidate_targets)), collapse = ";"),
    scored_genes = paste(sort(score_object$genes), collapse = ";")
  )
}
ifng_full_genes <- intersect(ifng_core_genes, rownames(log_expression))
ifng_nonoverlap_genes <- setdiff(ifng_full_genes, candidate_targets)
ifng_score_object <- score_gene_set(ifng_nonoverlap_genes, "IFN-gamma response")
outcome_score_objects[["IFN-gamma response"]] <- ifng_score_object
outcome_manifest[["IFN-gamma response"]] <- data.table(
  outcome_label = "IFN-gamma response", source_id = "Cross-collection core (>=2 of Hallmark/Reactome/GO:BP)",
  full_expressed_genes_n = length(ifng_full_genes),
  removed_candidate_regulon_targets_n = length(intersect(ifng_full_genes, candidate_targets)),
  scored_nonoverlap_genes_n = length(ifng_score_object$genes),
  removed_genes = paste(sort(intersect(ifng_full_genes, candidate_targets)), collapse = ";"),
  scored_genes = paste(sort(ifng_score_object$genes), collapse = ";")
)
outcome_manifest <- rbindlist(outcome_manifest, use.names = TRUE)
outcome_manifest[, outcome_order := match(outcome_label, c(program_catalog$outcome_label, "IFN-gamma response"))]
setorder(outcome_manifest, outcome_order)
stopifnot(all(outcome_manifest$scored_nonoverlap_genes_n >= 5L))

outcome_scores <- rbindlist(lapply(names(outcome_score_objects), function(label) {
  data.table(sample_id = names(outcome_score_objects[[label]]$score), outcome_label = label, nonoverlap_score = outcome_score_objects[[label]]$score)
}))
stopifnot(nrow(outcome_scores) == 6L * 73L, !anyNA(outcome_scores$nonoverlap_score))
outcome_wide <- dcast(outcome_scores, sample_id ~ outcome_label, value.var = "nonoverlap_score")
analysis_outcomes <- merge(analysis_data, outcome_wide, by = "sample_id", all.x = TRUE, sort = FALSE)
analysis_outcomes <- analysis_outcomes[match(sample_ids, sample_id)]

predictor_terms <- c(SPI1 = "SPI1_z", `NFKB1-RELA composite` = "NFKB_composite_z")
fit_outcome_model <- function(outcome, extra_terms = character(), model_label) {
  rhs <- c(unname(predictor_terms), "degradation_score_z", "therapy_short", "response_group", "age_z", "gender", extra_terms)
  fit <- lm(as.formula(paste0("`", outcome, "` ~ ", paste(rhs, collapse = " + "))), data = analysis_outcomes)
  rbindlist(lapply(names(predictor_terms), function(pred_label) {
    extract_lm_term(fit, predictor_terms[[pred_label]], outcome, pred_label, model_label)
  }))
}
primary_outcome_associations <- rbindlist(lapply(
  names(outcome_score_objects), fit_outcome_model,
  extra_terms = character(), model_label = "Joint SPI1 and NF-kappaB model"
))
primary_outcome_associations[, BH_q := p.adjust(p_value, method = "BH")]
sensitivity_outcome_associations <- rbindlist(list(
  rbindlist(lapply(
    names(outcome_score_objects), fit_outcome_model,
    extra_terms = c("CD8_T_cell_context", "Dendritic_cell_context"),
    model_label = "CD8 T-cell and dendritic-cell context"
  )),
  rbindlist(lapply(
    names(outcome_score_objects), fit_outcome_model,
    extra_terms = c("Pan_T_cell_context", "Monocytic_context"),
    model_label = "Pan-T-cell and monocytic context"
  ))
))
sensitivity_outcome_associations[, BH_q_within_model := p.adjust(p_value, method = "BH"), by = model]
sensitivity_outcome_associations[, BH_q_global := p.adjust(p_value, method = "BH")]

vif_for_term <- function(term, other_terms, data) {
  fit <- lm(as.formula(paste0(term, " ~ ", paste(other_terms, collapse = " + "))), data = data)
  1 / (1 - summary(fit)$r.squared)
}
vif_table <- data.table(
  predictor = names(predictor_terms),
  VIF = c(
    vif_for_term("SPI1_z", c("NFKB_composite_z", "degradation_score_z", "therapy_short", "response_group", "age_z", "gender"), analysis_outcomes),
    vif_for_term("NFKB_composite_z", c("SPI1_z", "degradation_score_z", "therapy_short", "response_group", "age_z", "gender"), analysis_outcomes)
  )
)

# Descriptive sample heatmap: standardized activities, focal genes and scores.
heat_feature_list <- list(
  "SPI1/PU.1 inferred activity" = activity_wide$SPI1_z,
  "NFKB1 inferred activity" = activity_wide$NFKB1_z,
  "RELA inferred activity" = activity_wide$RELA_z
)
focal_genes <- c("CD28", "CD80", "CD86", "PIK3CG", "VAV1", "TRIB3", "IFNG")
missing_focal_genes <- setdiff(focal_genes, rownames(log_expression_all))
if (length(missing_focal_genes)) stop("Missing focal genes: ", paste(missing_focal_genes, collapse = ", "))
for (gene in focal_genes) heat_feature_list[[gene]] <- zscore(log_expression_all[gene, sample_ids])
full_program_wide <- dcast(
  program_scores[program_id %in% c(
    "REACTOME_CO_STIMULATION_BY_CD28", "GOBP_T_CELL_PROLIFERATION", "GOBP_T_CELL_MEDIATED_CYTOTOXICITY"
  )], sample_id ~ program_id, value.var = "program_score"
)
full_program_wide <- full_program_wide[match(sample_ids, sample_id)]
heat_feature_list[["CD28 co-stimulation score"]] <- full_program_wide$REACTOME_CO_STIMULATION_BY_CD28
heat_feature_list[["T-cell proliferation score"]] <- full_program_wide$GOBP_T_CELL_PROLIFERATION
heat_feature_list[["T-cell cytotoxicity score"]] <- full_program_wide$GOBP_T_CELL_MEDIATED_CYTOTOXICITY
ifng_ordered <- ifng_scores[match(sample_ids, sample_id)]
stopifnot(identical(ifng_ordered$sample_id, sample_ids))
heat_feature_list[["IFN-gamma response score"]] <- zscore(ifng_ordered$IFNG_response_core_z)
heat_mat <- do.call(rbind, heat_feature_list); colnames(heat_mat) <- sample_ids
stopifnot(!anyNA(heat_mat), all(is.finite(heat_mat)))
column_order <- order(analysis_data$degradation_group, analysis_data$Degradation_score)
heat_mat <- heat_mat[, column_order, drop = FALSE]
heat_group <- data.frame(
  `Sarcosine degradation` = analysis_data$degradation_group[column_order],
  row.names = sample_ids[column_order], check.names = FALSE
)

# Provenance tables.
input_manifest <- data.table(
  input_role = c(
    "expression", "metadata", "sample TF activities", "prior TF contrasts",
    "DoRothEA regulon", "curated CD28 edges", "program membership",
    "program scores", "T-cell context", "IFNG scores", "cell-context estimates"
  ), path = required_inputs, sha256 = vapply(required_inputs, sha256_file, character(1))
)
method_contract <- data.table(
  item = c(
    "Population", "Group definition", "TF activity", "NF-kappaB summary",
    "Primary group model", "TF-to-outcome model", "Circularity control",
    "Cell-context sensitivity", "Multiple testing", "Interpretation limit", "Reader-facing labels"
  ),
  specification = c(
    "73 PRE-treatment TIGER melanoma tumors",
    "Median split of pre-existing z-scored SARDH/PIPOX transcriptional degradation module; Low n=37, High n=36",
    "ULM inferred activity from signed DoRothEA A/B regulons",
    "Mean of standardized NFKB1 and RELA activities, then standardized across tumors",
    "limma empirical Bayes; therapy + clinical response + age + sex + Degradation High/Low",
    "Each non-overlap outcome ~ SPI1 activity + NFKB1-RELA composite + continuous degradation score + therapy + clinical response + age + sex",
    "All expressed signed A/B targets of SPI1/NFKB1/RELA removed from each outcome gene set before scoring",
    "Named CD8 T-cell, pan-T-cell, dendritic-cell and monocytic scores only; sensitivity evidence from bulk RNA",
    "BH within each declared family; prior three-TF q values retain genome-wide 83-TF correction",
    "Associational bulk-RNA evidence; not causal, cell-intrinsic, protein activation, measured sarcosine or flux",
    "Degradation Low/High and biologically named axes only; broad composite adjustment labels are absent from figures"
  )
)
package_versions <- data.table(
  package = c("R", required_packages),
  version = c(R.version.string, vapply(required_packages, function(x) as.character(packageVersion(x)), character(1)))
)
direct_focal_edges <- direct_edges[TF %in% focal_tfs]
candidate_regulon_summary <- candidate_regulon[, .(
  expressed_targets_n = uniqueN(target), activating_targets_n = uniqueN(target[mor > 0]),
  repressing_targets_n = uniqueN(target[mor < 0]), confidence_A_targets_n = uniqueN(target[confidence == "A"]),
  confidence_B_targets_n = uniqueN(target[confidence == "B"])
), by = TF]

fwrite(input_manifest, file.path(table_dir, "00_input_manifest.csv"))
fwrite(method_contract, file.path(table_dir, "00_method_contract.csv"))
fwrite(package_versions, file.path(table_dir, "00_package_versions.csv"))
fwrite(primary_focused, file.path(table_dir, "01_focal_TF_activity_High_vs_Low.csv"))
fwrite(prior_crosscheck, file.path(table_dir, "02_prior_TF_effect_reproduction_check.csv"))
fwrite(activity_cor_long, file.path(table_dir, "03_focal_TF_activity_correlations.csv"))
fwrite(activity_wide, file.path(table_dir, "04_sample_focal_TF_activities.csv"))
fwrite(candidate_regulon_summary, file.path(table_dir, "05_focal_TF_regulon_summary.csv"))
fwrite(direct_focal_edges, file.path(table_dir, "06_focal_TF_to_CD28_direct_edges.csv"))
fwrite(activity_context_results, file.path(table_dir, "07_TF_activity_cell_context_associations.csv"))
fwrite(group_context_sensitivity, file.path(table_dir, "08_group_effect_after_named_cell_context_adjustment.csv"))
fwrite(outcome_manifest, file.path(table_dir, "09_nonoverlap_outcome_gene_sets.csv"))
fwrite(outcome_scores, file.path(table_dir, "10_nonoverlap_outcome_scores_long.csv"))
fwrite(primary_outcome_associations, file.path(table_dir, "11_SPI1_NFKB_to_nonoverlap_outcomes_primary.csv"))
fwrite(sensitivity_outcome_associations, file.path(table_dir, "12_SPI1_NFKB_to_nonoverlap_outcomes_cell_context_sensitivity.csv"))
fwrite(vif_table, file.path(table_dir, "13_joint_model_collinearity_VIF.csv"))

# Reader-facing figures.
low_color <- "#2A9D8F"; high_color <- "#CC3D3D"
nfkb_color <- "#D55E00"; spi1_color <- "#247BA0"
box_data <- melt(
  analysis_data[, .(sample_id, degradation_group, SPI1_z, NFKB1_z, RELA_z)],
  id.vars = c("sample_id", "degradation_group"), variable.name = "TF", value.name = "activity_z"
)
box_data[, TF := factor(TF, levels = c("SPI1_z", "NFKB1_z", "RELA_z"), labels = c("SPI1 / PU.1", "NFKB1", "RELA"))]
box_annotations <- primary_focused[feature %in% focal_tfs]
box_annotations[, TF := factor(feature, levels = focal_tfs, labels = c("SPI1 / PU.1", "NFKB1", "RELA"))]
box_annotations[, label := sprintf("effect = %.2f\nBH q = %.2g", standardized_effect, display_BH_q)]
box_annotations[, y := 3.55]
p_box <- ggplot(box_data, aes(degradation_group, activity_z, color = degradation_group)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, color = "#222222", linewidth = 0.55) +
  geom_jitter(width = 0.13, size = 1.35, alpha = 0.78) +
  geom_text(data = box_annotations, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 3.2, lineheight = 0.95) +
  facet_wrap(~ TF, nrow = 1) +
  scale_color_manual(values = c(Low = low_color, High = high_color), guide = "none") +
  coord_cartesian(ylim = c(-3.7, 4.2), clip = "off") +
  labs(
    title = "a  Inferred TF activities in PRE-treatment tumors",
    subtitle = "Signed DoRothEA A/B target-gene programs; effects are Degradation High minus Low",
    x = "Sarcosine-degradation group", y = "Inferred TF activity (z score)"
  ) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_rect(fill = "white", color = "#333333"), strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))

association_plot_data <- copy(primary_outcome_associations)
association_plot_data[, outcome := factor(outcome, levels = rev(c(program_catalog$outcome_label, "IFN-gamma response")))]
association_plot_data[, predictor_label := factor(
  predictor, levels = c("SPI1", "NFKB1-RELA composite"), labels = c("SPI1 / PU.1", "NFKB1-RELA (NF-kappaB)")
)]
association_plot_data[, significant := BH_q < 0.05]
p_assoc <- ggplot(association_plot_data, aes(standardized_beta, outcome, color = predictor_label)) +
  geom_vline(xintercept = 0, color = "#777777", linewidth = 0.45, linetype = 2) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0.16, position = position_dodge(width = 0.48), linewidth = 0.65) +
  geom_point(aes(shape = significant), size = 2.8, position = position_dodge(width = 0.48), stroke = 0.8) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "BH q < 0.05", `FALSE` = "BH q >= 0.05")) +
  scale_color_manual(values = c("SPI1 / PU.1" = spi1_color, "NFKB1-RELA (NF-kappaB)" = nfkb_color)) +
  labs(
    title = "b  Associations with non-overlapping functional programs",
    subtitle = "Both axes entered together; shared regulon targets were removed from every outcome score",
    x = "Standardized association coefficient", y = NULL, color = "Regulatory axis", shape = "Multiplicity"
  ) +
  theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

context_plot_data <- activity_context_results[outcome %in% focal_tfs]
context_plot_data[, TF := factor(outcome, levels = rev(focal_tfs), labels = rev(c("SPI1 / PU.1", "NFKB1", "RELA")))]
context_plot_data[, context := factor(predictor, levels = unname(context_labels))]
context_plot_data[, tile_label := sprintf("%.2f%s", standardized_beta, ifelse(BH_q < 0.05, "*", ""))]
p_context <- ggplot(context_plot_data, aes(context, TF, fill = standardized_beta)) +
  geom_tile(color = "white", linewidth = 1) + geom_text(aes(label = tile_label), size = 3.2) +
  scale_fill_gradient2(low = low_color, mid = "white", high = high_color, midpoint = 0, limits = c(-1.2, 1.2), oob = scales::squish) +
  labs(
    title = "c  Cell-context associations",
    subtitle = "Clinical and continuous degradation-score adjusted; * BH q < 0.05",
    x = NULL, y = NULL, fill = "Standardized\ncoefficient"
  ) +
  theme_minimal(base_size = 10.5) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"), axis.text.x = element_text(angle = 0, hjust = 0.5), legend.position = "right")

summary_figure <- p_box / p_assoc / p_context +
  plot_layout(heights = c(1.00, 0.82, 0.66)) +
  plot_annotation(
    title = "SPI1/PU.1 and canonical NF-kappaB programs in Sarcosine Degradation High versus Low tumors",
    subtitle = "TIGER melanoma, PRE-treatment only (Low n=37; High n=36)",
    theme = theme(plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 11))
  )
ggsave(file.path(figure_dir, "Fig_SPI1_NFKB1_RELA_Focused_Summary.png"), summary_figure, width = 14.5, height = 15.2, dpi = 400, bg = "white")
ggsave(file.path(figure_dir, "Fig_SPI1_NFKB1_RELA_Focused_Summary.pdf"), summary_figure, width = 14.5, height = 15.2, device = cairo_pdf, bg = "white")

pheatmap(
  heat_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE,
  fontsize_row = 9.5, gaps_row = c(3, 10), gaps_col = 37, annotation_col = heat_group,
  annotation_colors = list(`Sarcosine degradation` = c(Low = low_color, High = high_color)),
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks = seq(-2.5, 2.5, length.out = 102), border_color = NA,
  main = "SPI1/NF-kappaB-CD28-IFN-gamma features across PRE-treatment tumors",
  filename = file.path(figure_dir, "Fig_SPI1_NFKB1_RELA_Sample_Heatmap.png"), width = 13.5, height = 7.6
)
pheatmap(
  heat_mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE,
  fontsize_row = 9.5, gaps_row = c(3, 10), gaps_col = 37, annotation_col = heat_group,
  annotation_colors = list(`Sarcosine degradation` = c(Low = low_color, High = high_color)),
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks = seq(-2.5, 2.5, length.out = 102), border_color = NA,
  main = "SPI1/NF-kappaB-CD28-IFN-gamma features across PRE-treatment tumors",
  filename = file.path(figure_dir, "Fig_SPI1_NFKB1_RELA_Sample_Heatmap.pdf"), width = 13.5, height = 7.6
)

# Reports and calibrated interpretation.
sig_primary <- primary_outcome_associations[BH_q < 0.05][order(BH_q)]
context_sig <- activity_context_results[BH_q < 0.05][order(BH_q)]
format_assoc_lines <- function(x) {
  if (!nrow(x)) return("- None at BH q < 0.05 in the declared family.")
  vapply(seq_len(nrow(x)), function(i) sprintf(
    "- %s -> %s: beta %.2f (95%% CI %.2f to %.2f), BH q %.3g.",
    x$predictor[i], x$outcome[i], x$standardized_beta[i], x$CI_low[i], x$CI_high[i], x$BH_q[i]
  ), character(1))
}
report_lines <- c(
  "# Focused SPI1/PU.1-NF-kappaB-CD28-IFN-gamma analysis", "",
  "## Scope", "",
  "This reproducible analysis narrows the earlier CD28-linked TF screen to SPI1, NFKB1 and RELA in 73 PRE-treatment TIGER melanoma tumors. Reader-facing figures show only Degradation Low/High and biologically named regulatory or cell-context axes.", "",
  "The grouping variable is a median split of a SARDH/PIPOX transcriptional degradation module. It is not a direct measurement of sarcosine concentration or metabolic flux.", "",
  "## Primary TF result", "",
  sprintf("- SPI1/PU.1: standardized High-minus-Low effect %.2f (95%% CI %.2f to %.2f), prior genome-wide BH q %.3g.", primary_focused[feature == "SPI1", standardized_effect], primary_focused[feature == "SPI1", display_CI_low], primary_focused[feature == "SPI1", display_CI_high], primary_focused[feature == "SPI1", display_BH_q]),
  sprintf("- NFKB1: standardized effect %.2f (95%% CI %.2f to %.2f), prior genome-wide BH q %.3g.", primary_focused[feature == "NFKB1", standardized_effect], primary_focused[feature == "NFKB1", display_CI_low], primary_focused[feature == "NFKB1", display_CI_high], primary_focused[feature == "NFKB1", display_BH_q]),
  sprintf("- RELA: standardized effect %.2f (95%% CI %.2f to %.2f), prior genome-wide BH q %.3g.", primary_focused[feature == "RELA", standardized_effect], primary_focused[feature == "RELA", display_CI_low], primary_focused[feature == "RELA", display_CI_high], primary_focused[feature == "RELA", display_BH_q]),
  sprintf("- NFKB1-RELA activity correlation: Pearson r = %.3f. They were combined for downstream models.", activity_cor["NFKB1_z", "RELA_z"]), "",
  "## Non-overlapping outcome analysis", "",
  "SPI1/PU.1 and the NFKB1-RELA composite were entered together with the continuous degradation score and clinical covariates. Every expressed signed A/B target of the three TFs was removed from each outcome gene set before scoring, reducing direct target-overlap circularity.", "",
  format_assoc_lines(sig_primary), "",
  sprintf("Joint-model VIFs were %.2f for SPI1/PU.1 and %.2f for the NFKB1-RELA composite; estimates should be read with this collinearity in mind.", vif_table[predictor == "SPI1", VIF], vif_table[predictor == "NFKB1-RELA composite", VIF]), "",
  "The direct CD28-program associations did not pass the declared 12-test BH threshold after target removal and joint modeling. The stronger evidence is therefore for association with downstream T-cell/IFN-gamma functional state, not for unique causal control of CD28 itself.", "",
  "## Curated TF-to-CD28 links", "",
  "The signed DoRothEA A/B prior-knowledge edges connecting these TFs to the Reactome CD28 set were SPI1 -> CD86, PIK3CG and VAV1; NFKB1 -> CD80 and CD86; and RELA -> CD86 with a negative edge to TRIB3. These are literature-curated regulatory priors, not cohort-specific causal edges.", "",
  "## Cell-context findings", "",
  "The following are associations with named bulk-RNA cell-context scores, not measured cell fractions and not proof of a cell of origin:", "",
  if (nrow(context_sig)) vapply(seq_len(nrow(context_sig)), function(i) sprintf(
    "- %s with %s: beta %.2f, BH q %.3g.", context_sig$outcome[i], context_sig$predictor[i], context_sig$standardized_beta[i], context_sig$BH_q[i]
  ), character(1)) else "- None at BH q < 0.05.", "",
  "## High-versus-Low sensitivity to named cell contexts", "",
  sprintf(
    "- After CD8 T-cell-context adjustment, all three TF effects remained significant (within-model BH q <= %.3g).",
    max(group_context_sensitivity[model == "Adjusted for CD8 T-cell context" & feature %in% focal_tfs, BH_q_within_model])
  ),
  sprintf(
    "- After pan-T-cell-context adjustment, none of the three passed BH q < 0.05 (BH q range %.3g to %.3g).",
    min(group_context_sensitivity[model == "Adjusted for Pan-T-cell context" & feature %in% focal_tfs, BH_q_within_model]),
    max(group_context_sensitivity[model == "Adjusted for Pan-T-cell context" & feature %in% focal_tfs, BH_q_within_model])
  ),
  sprintf(
    "- After dendritic-cell-context adjustment, all three remained significant (within-model BH q <= %.3g).",
    max(group_context_sensitivity[model == "Adjusted for Dendritic-cell context" & feature %in% focal_tfs, BH_q_within_model])
  ),
  sprintf(
    "- After monocytic-context adjustment, all three remained significant within the focused family (BH q <= %.3g; global 16-test BH q <= %.3g).",
    max(group_context_sensitivity[model == "Adjusted for Monocytic context" & feature %in% focal_tfs, BH_q_within_model]),
    max(group_context_sensitivity[model == "Adjusted for Monocytic context" & feature %in% focal_tfs, BH_q_global])
  ), "",
  "Because all context estimates and TF activities come from the same bulk RNA profiles, attenuation after pan-T-cell adjustment is compatible with shared T-cell abundance and state; it does not identify a causal mediator.", "",
  "## Biological interpretation", "",
  "NFKB1/p50 and RELA/p65 are most defensibly treated here as canonical NF-kappaB transcriptional readouts compatible with TCR/CD28 signaling, rather than as independent upstream causes of CD28 activation. SPI1/PU.1 is a distinct lineage-context candidate, especially relevant to antigen-presenting/myeloid and dendritic-cell programs. The bulk-tumor data support a coordinated immune-rich state; they do not establish the causal chain sarcosine degradation -> CD28 -> NF-kappaB -> IFNG.", "",
  "## Outputs", "",
  "- `figures/Fig_SPI1_NFKB1_RELA_Focused_Summary.png`: group effects, non-overlap functional associations and named cell contexts.",
  "- `figures/Fig_SPI1_NFKB1_RELA_Sample_Heatmap.png`: descriptive tumor-by-feature heatmap ordered Low then High.",
  "- `tables/`: inputs, method contract, activities, correlations, gene-set manifests, association models and validation.", "",
  "## Limitations", "",
  "All TF, pathway and cell-context quantities are inferred from the same bulk RNA-seq data. The results are observational and may reflect immune-cell abundance, activation state, or both. Direct confirmation requires orthogonal protein/phospho measurements, cell-resolved data or perturbation experiments."
)
writeLines(report_lines, file.path(output_root, "FINAL_SPI1_NFKB1_RELA_FOCUSED_REPORT.md"))

literature_lines <- c(
  "# Biological context for the three focused TFs", "", "## NFKB1 and RELA", "",
  "- CD28-responsive NF-kappaB complexes containing p50/NFKB1 and p65/RELA have been demonstrated in T cells. Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC45946/",
  "- CD28 co-stimulation can accelerate nuclear translocation of p50 and p65, supporting their interpretation as downstream pathway readouts. Source: https://pmc.ncbi.nlm.nih.gov/articles/PMC359332/", "",
  "## SPI1 / PU.1", "",
  "- PU.1 is required for conventional dendritic-cell identity and function, making SPI1 activity a plausible antigen-presenting/myeloid-context signal in bulk tumors. Source: https://www.sciencedirect.com/science/article/pii/S1074761318304928",
  "- Independent developmental evidence supports a role for PU.1 in dendritic-cell development. Source: https://pubmed.ncbi.nlm.nih.gov/20510871/", "",
  "## Inference framework", "",
  "- DoRothEA infers TF activity from signed regulon targets and was benchmarked against perturbation evidence; activity scores are not equivalent to TF expression or protein activation. Source: https://pubmed.ncbi.nlm.nih.gov/31340985/",
  "- CollecTRI provides a complementary literature-curated TF-target resource for future robustness analyses. Source: https://pubmed.ncbi.nlm.nih.gov/37843125/", "",
  "## Calibrated conclusion", "",
  "The current data justify a two-axis hypothesis: an SPI1/PU.1-associated antigen-presenting-cell context and a NFKB1/RELA canonical NF-kappaB activation state. They do not prove that SPI1 regulates CD28 in T cells or that sarcosine degradation directly activates either axis."
)
writeLines(literature_lines, file.path(output_root, "TF_BIOLOGICAL_CONTEXT.md"))

caption_lines <- c(
  "# Reader-facing figure captions", "", "## Fig_SPI1_NFKB1_RELA_Focused_Summary", "",
  "Focused SPI1/PU.1 and canonical NF-kappaB analysis in 73 PRE-treatment TIGER melanoma tumors classified by the median SARDH/PIPOX transcriptional degradation score (Low, n=37; High, n=36). (a) Inferred SPI1, NFKB1 and RELA activities from signed DoRothEA A/B regulons. Effects are adjusted High-minus-Low differences controlling for therapy, clinical response, age and sex; BH q values retain the original genome-wide 83-TF correction. (b) Joint associations of SPI1/PU.1 and a NFKB1-RELA composite with CD28, T-cell and IFN-gamma scores after removing every expressed signed target of the three TFs from each outcome gene set. Models also included continuous degradation score and clinical covariates. Filled points indicate BH q<0.05. (c) Associations between each TF activity and biologically named bulk-RNA cell-context scores after accounting for continuous degradation score and clinical covariates. These are computational estimates, not measured cell fractions.", "",
  "## Fig_SPI1_NFKB1_RELA_Sample_Heatmap", "",
  "Descriptive standardized feature heatmap across PRE-treatment tumors ordered by Sarcosine Degradation Low then High and by the continuous degradation score within each group. Rows show inferred TF activities, selected CD28-linked genes and functional scores. This panel visualizes co-variation and is not an independent test of causality or cell-intrinsic signaling."
)
writeLines(caption_lines, file.path(output_root, "FIGURE_CAPTIONS.md"))

validation <- data.table(
  check = c(
    "PRE tumors n=73", "Low n=37", "High n=36", "TF prior effects exactly reproduced",
    "NFKB1/RELA correlation finite", "All non-overlap scores >=5 genes",
    "No focal regulon target remains in scored outcomes", "Joint model VIF finite",
    "Summary PNG exists", "Heatmap PNG exists", "Reader-facing text excludes prohibited model label"
  ),
  passed = c(
    nrow(meta) == 73L, sum(meta$degradation_group == "Low") == 37L,
    sum(meta$degradation_group == "High") == 36L,
    max(prior_crosscheck$absolute_difference) < 1e-8 &&
      all(prior_crosscheck$prior_standardized_CI_low < prior_crosscheck$prior_standardized_CI_high),
    is.finite(activity_cor["NFKB1_z", "RELA_z"]),
    all(outcome_manifest$scored_nonoverlap_genes_n >= 5L),
    all(vapply(outcome_score_objects, function(x) length(intersect(x$genes, candidate_targets)) == 0L, logical(1))),
    all(is.finite(vif_table$VIF)),
    file.exists(file.path(figure_dir, "Fig_SPI1_NFKB1_RELA_Focused_Summary.png")),
    file.exists(file.path(figure_dir, "Fig_SPI1_NFKB1_RELA_Sample_Heatmap.png")),
    !any(grepl("ESTIMATE", c(caption_lines, report_lines), fixed = TRUE))
  )
)
fwrite(validation, file.path(table_dir, "14_validation_checks.csv"))
if (!all(validation$passed)) stop("Validation failure: ", paste(validation[passed == FALSE, check], collapse = "; "))
sink(file.path(log_dir, "sessionInfo.txt")); print(sessionInfo()); sink()
message("Focused SPI1/NFKB1/RELA analysis completed: ", output_root)
