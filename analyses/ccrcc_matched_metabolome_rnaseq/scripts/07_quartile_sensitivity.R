#!/usr/bin/env Rscript

# Exploratory extreme-group sensitivity: top versus bottom quartile of measured
# tumour GC-MS Sarcosine. This analysis does not replace the frozen median split.

options(stringsAsFactors = FALSE, width = 180)
set.seed(260902)

analysis_root <- normalizePath(getwd())
project_root <- normalizePath(file.path(analysis_root, ".."))
gsea_lib <- file.path(
  project_root, "2024_Drug_Res_Updates_NSCLC", "RNA-seq공공데이터_GSE207422",
  "analysis_sarcosine_FINAL_26.08.25", "03_lee_fig3_style_FINAL",
  "hallmark_GSEA_Q4_vs_Q1_26.08.26", "R_libs"
)
.libPaths(unique(c(normalizePath(gsea_lib), .libPaths())))

required <- c("data.table", "limma", "fgsea", "msigdbr", "ggplot2", "pheatmap", "BiocParallel")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table); library(limma); library(fgsea); library(msigdbr)
  library(ggplot2); library(pheatmap)
})

expression_path <- file.path(analysis_root, "RNAseq_Data", "bulkRNA_matrix_TPM.csv")
group_path <- file.path(analysis_root, "results", "00_input_audit_26.09.02", "tumor_sarcosine_group_map.csv")
deconv_path <- file.path(analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_HighLow_Deconvolution_26.09.02", "tables", "core_celltype_values_long.csv")
cd8_path <- file.path(analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_HighLow_Deconvolution_26.09.02", "tables", "CD8_consensus_per_tumour.csv")
program_path <- file.path(analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_HighLow_Immune_TF_Integrated_26.09.02", "tables", "immune_program_scores_per_tumour.csv")
tf_path <- file.path(analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_HighLow_Immune_TF_Integrated_26.09.02", "tables", "sample_level_TF_activities_long.csv")
tf_meta_path <- file.path(analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_HighLow_Immune_TF_Integrated_26.09.02", "tables", "DoRothEA_TF_activity_all_eligible_High_vs_Low.csv")
inputs <- c(expression_path, group_path, deconv_path, cd8_path, program_path, tf_path, tf_meta_path)
if (!all(file.exists(inputs))) stop("Missing input(s): ", paste(inputs[!file.exists(inputs)], collapse = ", "))

out_root <- file.path(analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_Q4_vs_Q1_Extreme_26.09.02")
table_dir <- file.path(out_root, "tables"); figure_dir <- file.path(out_root, "figures"); log_dir <- file.path(out_root, "logs")
for (d in c(table_dir, figure_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

meta_all <- fread(group_path, check.names = FALSE)
stopifnot(nrow(meta_all) == 100L, uniqueN(meta_all$sample_id) == 100L)
cut_q1 <- as.numeric(quantile(meta_all$sarcosine_normalized_intensity, 0.25, type = 7))
cut_q4 <- as.numeric(quantile(meta_all$sarcosine_normalized_intensity, 0.75, type = 7))
meta_all[, extreme_group := fifelse(
  sarcosine_normalized_intensity <= cut_q1, "Sarcosine-Low (Q1)",
  fifelse(sarcosine_normalized_intensity >= cut_q4, "Sarcosine-High (Q4)", "Middle 50% excluded")
)]
meta <- copy(meta_all[extreme_group != "Middle 50% excluded"])
meta[, extreme_group := factor(extreme_group, levels = c("Sarcosine-Low (Q1)", "Sarcosine-High (Q4)"))]
stopifnot(nrow(meta) == 50L, sum(meta$extreme_group == "Sarcosine-Low (Q1)") == 25L, sum(meta$extreme_group == "Sarcosine-High (Q4)") == 25L)
stopifnot(max(meta[extreme_group == "Sarcosine-Low (Q1)", sarcosine_normalized_intensity]) < cut_q1)
stopifnot(min(meta[extreme_group == "Sarcosine-High (Q4)", sarcosine_normalized_intensity]) > cut_q4)

meta[, batch := factor(batch)]
meta[, sex := relevel(factor(sex), ref = "female")]
meta[, age := factor(age, levels = c("40-60", "<40", ">60"))]
meta[, grade := factor(grade, levels = c("1", "2", "3", "4"))]
meta[, stage := factor(stage, levels = c("stage I", "stage II", "stage III", "stage IV"))]
sample_ids <- meta$sample_id

group_manifest <- rbind(
  data.table(group = "Sarcosine-Low (Q1)", rule = sprintf("intensity <= 25th percentile %.10f", cut_q1), n = 25L),
  data.table(group = "Middle 50% excluded", rule = sprintf("%.10f < intensity < %.10f", cut_q1, cut_q4), n = 50L),
  data.table(group = "Sarcosine-High (Q4)", rule = sprintf("intensity >= 75th percentile %.10f", cut_q4), n = 25L)
)
fwrite(group_manifest, file.path(table_dir, "Q1_Q4_group_definition.csv"))
fwrite(meta_all[, .(sample_id, sarcosine_normalized_intensity, log2_sarcosine_intensity, extreme_group, batch, sex, age, stage, grade)], file.path(table_dir, "Q1_Q4_all_tumour_group_map.csv"))

expr_dt <- fread(expression_path, check.names = FALSE)
gene <- trimws(expr_dt[[1]])
all_sample_ids <- meta_all$sample_id
expr_all <- as.matrix(expr_dt[, ..all_sample_ids]); storage.mode(expr_all) <- "double"
rownames(expr_all) <- gene; colnames(expr_all) <- all_sample_ids
# Preserve the exact 15,119-gene universe used in the median-split analysis,
# then subset samples. This prevents an extreme-group-specific filter from
# creating a different GSEA universe.
keep <- rowSums(expr_all >= 1) >= 10 & apply(expr_all, 1L, var) > 0
expr <- expr_all[keep, sample_ids, drop = FALSE]
stopifnot(nrow(expr) == 15119L, !anyNA(expr), all(is.finite(expr)))

designs <- list(
  adjusted_Q4_vs_Q1 = model.matrix(~ batch + sex + age + grade + extreme_group, data = meta),
  clinical_no_batch_Q4_vs_Q1 = model.matrix(~ sex + age + grade + extreme_group, data = meta),
  unadjusted_Q4_vs_Q1 = model.matrix(~ extreme_group, data = meta),
  stage_instead_of_grade_Q4_vs_Q1 = model.matrix(~ batch + sex + age + stage + extreme_group, data = meta)
)
coef_name <- "extreme_groupSarcosine-High (Q4)"
for (nm in names(designs)) {
  if (qr(designs[[nm]])$rank != ncol(designs[[nm]])) stop("Rank-deficient design: ", nm)
  if (!coef_name %in% colnames(designs[[nm]])) stop("Missing contrast coefficient: ", nm)
}

fit_matrix <- function(mat, design, label) {
  fit <- eBayes(lmFit(mat, design), trend = TRUE, robust = TRUE)
  tt <- as.data.table(topTable(fit, coef = coef_name, number = Inf, sort.by = "none", adjust.method = "BH"), keep.rownames = "feature")
  se <- fit$stdev.unscaled[, coef_name] * sqrt(fit$s2.post)
  crit <- qt(0.975, fit$df.total); if (length(crit) == 1L) crit <- rep(crit, nrow(tt))
  tt[, standard_error := se[feature]]
  tt[, `:=`(CI_low = logFC - crit * standard_error, CI_high = logFC + crit * standard_error, model = label)]
  setnames(tt, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"), c("effect_Q4_minus_Q1", "average", "moderated_t", "p_value", "BH_q", "B_statistic"))
  stopifnot(nrow(tt) == nrow(mat), all(is.finite(tt$moderated_t)))
  list(fit = fit, table = tt)
}

gene_fits <- lapply(names(designs), function(nm) fit_matrix(expr, designs[[nm]], nm)); names(gene_fits) <- names(designs)
for (nm in names(gene_fits)) fwrite(gene_fits[[nm]]$table, file.path(table_dir, paste0("limma_", nm, "_all_genes.csv")))
gene_summary <- rbindlist(lapply(names(gene_fits), function(nm) {
  x <- gene_fits[[nm]]$table
  data.table(model = nm, tested_genes = nrow(x), positive_BH_q_lt_0_05 = sum(x$BH_q < 0.05 & x$effect_Q4_minus_Q1 > 0), negative_BH_q_lt_0_05 = sum(x$BH_q < 0.05 & x$effect_Q4_minus_Q1 < 0))
}))
fwrite(gene_summary, file.path(table_dir, "limma_model_summary.csv"))

# GSEA: adjusted extreme-group contrast, three collections, fixed MSigDB release.
msig <- as.data.table(msigdbr(species = "Homo sapiens", db_species = "HS"))
db_version <- unique(msig$db_version); stopifnot(length(db_version) == 1L)
msig[, collection := fifelse(gs_collection == "H", "Hallmark", fifelse(gs_collection == "C2" & gs_subcollection == "CP:REACTOME", "Reactome", fifelse(gs_collection == "C5" & gs_subcollection == "GO:BP", "GO:BP", NA_character_)))]
msig <- unique(msig[!is.na(collection), .(collection, pathway = gs_name, gene_symbol)])
pathway_lists <- lapply(split(msig, msig$collection), function(x) split(x$gene_symbol, x$pathway))
rank_stats <- setNames(gene_fits$adjusted_Q4_vs_Q1$table$moderated_t, gene_fits$adjusted_Q4_vs_Q1$table$feature)

run_gsea <- function(paths, collection) {
  min_size <- if (collection == "Hallmark") 10L else 15L
  ans <- as.data.table(fgseaMultilevel(
    pathways = paths, stats = rank_stats, minSize = min_size, maxSize = 500L,
    eps = 0, nPermSimple = 100000L, nproc = 1,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  ))
  ans[, `:=`(collection = collection, direction = fifelse(NES > 0, "Sarcosine-High (Q4)", "Sarcosine-Low (Q1)"))]
  ans
}
gsea <- rbindlist(lapply(names(pathway_lists), function(collection) run_gsea(pathway_lists[[collection]], collection)), use.names = TRUE)
gsea_out <- copy(gsea); gsea_out[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
gsea_out[, abs_NES_sort := abs(NES)]; setorder(gsea_out, padj, -abs_NES_sort); gsea_out[, abs_NES_sort := NULL]
fwrite(gsea_out, file.path(table_dir, "GSEA_Hallmark_Reactome_GOBP_adjusted_Q4_vs_Q1.csv"))
fwrite(gsea_out[is.na(NES) | is.na(padj), .(collection, pathway, reason = "fgsea unbalanced statistic; NES/p/q not estimable")], file.path(table_dir, "GSEA_unbalanced_NA_audit.csv"))
gsea_summary <- gsea[, .(
  tested = .N,
  positive_BH_q_lt_0_05 = sum(padj < 0.05 & NES > 0, na.rm = TRUE),
  negative_BH_q_lt_0_05 = sum(padj < 0.05 & NES < 0, na.rm = TRUE),
  NA_pathways = sum(is.na(NES) | is.na(padj))
), by = collection]
fwrite(gsea_summary, file.path(table_dir, "GSEA_collection_summary.csv"))

focus_ids <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "REACTOME_INTERFERON_GAMMA_SIGNALING", "REACTOME_CO_STIMULATION_BY_CD28", "REACTOME_CD28_DEPENDENT_PI3K_AKT_SIGNALING", "REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING",
  "GOBP_T_CELL_PROLIFERATION", "GOBP_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE", "GOBP_T_CELL_MEDIATED_CYTOTOXICITY", "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY",
  "GOBP_NATURAL_KILLER_CELL_ACTIVATION", "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I"
)
focus <- gsea[pathway %in% focus_ids]
stopifnot(uniqueN(focus$pathway) == length(focus_ids))
fwrite(copy(focus)[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")], file.path(table_dir, "GSEA_IFNG_CD28_IL12_focus_Q4_vs_Q1.csv"))

# Reuse already-computed sample-level deconvolution values; only the contrast is new.
deconv <- fread(deconv_path)
deconv <- deconv[sample_id %in% sample_ids]
feature_meta <- unique(deconv[, .(feature_id, method, cell_type, measure_type, compartment)])
deconv_wide <- dcast(deconv, feature_id ~ sample_id, value.var = "analysis_value")
deconv_mat <- as.matrix(deconv_wide[, ..sample_ids]); storage.mode(deconv_mat) <- "double"; rownames(deconv_mat) <- deconv_wide$feature_id
stopifnot(!anyNA(deconv_mat), all(is.finite(deconv_mat)))
deconv_fit <- fit_matrix(deconv_mat, designs$adjusted_Q4_vs_Q1, "Adjusted Q4 vs Q1")$table
deconv_fit <- merge(deconv_fit, feature_meta, by.x = "feature", by.y = "feature_id", all.x = TRUE)
deconv_fit[, BH_q_all_cell_estimates := p.adjust(p_value, method = "BH")]
fwrite(deconv_fit, file.path(table_dir, "deconvolution_adjusted_Q4_vs_Q1.csv"))

cd8 <- fread(cd8_path)[sample_id %in% sample_ids]
cd8 <- cd8[match(sample_ids, sample_id)]
cd8_mat <- rbind(
  `CD8 consensus (six methods)` = cd8$CD8_consensus_z,
  `CD8 consensus (excluding CIBERSORT)` = cd8$CD8_consensus_no_CIBERSORT_z
)
colnames(cd8_mat) <- sample_ids
cd8_fit <- fit_matrix(cd8_mat, designs$adjusted_Q4_vs_Q1, "Adjusted Q4 vs Q1")$table
fwrite(cd8_fit, file.path(table_dir, "CD8_consensus_adjusted_Q4_vs_Q1.csv"))

# Program-score models with abundance/TME sensitivity adjustments.
program <- fread(program_path)[sample_id %in% sample_ids]
program_meta <- unique(program[, .(program_id, program_label, family, collection, display_order)])
program_wide <- dcast(program, program_id ~ sample_id, value.var = "program_score_z")
program_mat <- as.matrix(program_wide[, ..sample_ids]); storage.mode(program_mat) <- "double"; rownames(program_mat) <- program_wide$program_id

pan_t <- deconv[method == "MCPcounter" & cell_type == "Pan T cells", .(sample_id, pan_T_z = analysis_value)]
estimate <- dcast(deconv[method == "ESTIMATE" & cell_type %in% c("Immune score", "Stromal score"), .(sample_id, cell_type, analysis_value)], sample_id ~ cell_type, value.var = "analysis_value")
setnames(estimate, c("Immune score", "Stromal score"), c("immune_score_z", "stromal_score_z"))
model_meta <- Reduce(function(x, y) merge(x, y, by = "sample_id", all.x = TRUE, sort = FALSE), list(
  meta, cd8[, .(sample_id, CD8_consensus_z, CD8_consensus_no_CIBERSORT_z)], pan_t, estimate
))
model_meta <- model_meta[match(sample_ids, sample_id)]
stopifnot(!anyNA(model_meta[, .(CD8_consensus_z, CD8_consensus_no_CIBERSORT_z, pan_T_z, immune_score_z, stromal_score_z)]))
program_formulas <- list(
  `Q4 vs Q1` = ~ batch + sex + age + grade + extreme_group,
  `Adjusted for CD8 consensus excluding CIBERSORT` = ~ batch + sex + age + grade + CD8_consensus_no_CIBERSORT_z + extreme_group,
  `Adjusted for pan-T estimate` = ~ batch + sex + age + grade + pan_T_z + extreme_group,
  `Adjusted for immune/stromal scores` = ~ batch + sex + age + grade + immune_score_z + stromal_score_z + extreme_group
)
program_effects <- rbindlist(lapply(names(program_formulas), function(label) {
  d <- model.matrix(program_formulas[[label]], data = model_meta)
  if (qr(d)$rank != ncol(d)) stop("Rank-deficient program model: ", label)
  x <- fit_matrix(program_mat, d, label)$table
  x[, BH_q := p.adjust(p_value, "BH")]
  x[, adjustment := label]
  x
}))
program_effects <- merge(program_effects, program_meta, by.x = "feature", by.y = "program_id", all.x = TRUE)
fwrite(program_effects, file.path(table_dir, "immune_program_Q4_vs_Q1_abundance_adjusted_effects.csv"))

# Signed DoRothEA TF activities were already computed per sample; fit the new contrast.
tf_long <- fread(tf_path)[sample_id %in% sample_ids]
tf_wide <- dcast(tf_long, TF ~ sample_id, value.var = "ULM_activity")
tf_mat <- as.matrix(tf_wide[, ..sample_ids]); storage.mode(tf_mat) <- "double"; rownames(tf_mat) <- tf_wide$TF
tf_fit <- fit_matrix(tf_mat, designs$adjusted_Q4_vs_Q1, "Adjusted Q4 vs Q1")$table
tf_sd <- apply(tf_mat, 1L, sd)
tf_fit[, `:=`(
  standardized_effect = effect_Q4_minus_Q1 / tf_sd[feature],
  standardized_CI_low = CI_low / tf_sd[feature],
  standardized_CI_high = CI_high / tf_sd[feature],
  BH_q = p.adjust(p_value, "BH")
)]
tf_meta <- unique(fread(tf_meta_path)[, .(TF, expressed_targets_n, CD28_pathway_targets_n, CD28_pathway_targets, direct_CD28_gene_regulator)])
tf_fit <- merge(tf_fit, tf_meta, by.x = "feature", by.y = "TF", all.x = TRUE)
setnames(tf_fit, "feature", "TF")
tf_fit[, abs_effect_sort := abs(standardized_effect)]; setorder(tf_fit, BH_q, -abs_effect_sort); tf_fit[, abs_effect_sort := NULL]
fwrite(tf_fit, file.path(table_dir, "DoRothEA_TF_activity_adjusted_Q4_vs_Q1.csv"))
tf_linked <- tf_fit[CD28_pathway_targets_n > 0 | direct_CD28_gene_regulator == TRUE]
fwrite(tf_linked, file.path(table_dir, "CD28_linked_TF_activity_adjusted_Q4_vs_Q1.csv"))

# Prespecified influence check: all five Tukey-low Sarcosine observations fall
# in Q1. Refit the same Q4-vs-Q1 models after excluding them (Q1 n=20;
# Q4 n=25). This does not redefine the quartile cutpoints or the gene universe.
tukey_low_ids <- c("H46_T", "N39_T", "R25_T", "R82_T", "Z16_T")
stopifnot(all(tukey_low_ids %in% meta[extreme_group == "Sarcosine-Low (Q1)", sample_id]))
meta_no_low <- copy(meta[!sample_id %in% tukey_low_ids])
meta_no_low[, extreme_group := droplevels(extreme_group)]
no_low_ids <- meta_no_low$sample_id
stopifnot(
  nrow(meta_no_low) == 45L,
  sum(meta_no_low$extreme_group == "Sarcosine-Low (Q1)") == 20L,
  sum(meta_no_low$extreme_group == "Sarcosine-High (Q4)") == 25L
)
design_no_low <- model.matrix(~ batch + sex + age + grade + extreme_group, data = meta_no_low)
if (qr(design_no_low)$rank != ncol(design_no_low)) stop("Rank-deficient Tukey-low-exclusion design")

gene_no_low <- fit_matrix(expr[, no_low_ids, drop = FALSE], design_no_low, "Adjusted Q4 vs Q1; exclude five Tukey-low values")$table
fwrite(gene_no_low, file.path(table_dir, "limma_adjusted_Q4_vs_Q1_exclude5_Tukey_low_all_genes.csv"))
gene_summary_no_low <- data.table(
  model = "adjusted_Q4_vs_Q1_exclude5_Tukey_low",
  tested_genes = nrow(gene_no_low),
  positive_BH_q_lt_0_05 = sum(gene_no_low$BH_q < 0.05 & gene_no_low$effect_Q4_minus_Q1 > 0),
  negative_BH_q_lt_0_05 = sum(gene_no_low$BH_q < 0.05 & gene_no_low$effect_Q4_minus_Q1 < 0)
)
fwrite(rbind(gene_summary, gene_summary_no_low), file.path(table_dir, "limma_model_summary_with_outlier_sensitivity.csv"))

main_gene <- gene_fits$adjusted_Q4_vs_Q1$table[, .(
  feature, main_effect = effect_Q4_minus_Q1, main_p = p_value, main_BH_q = BH_q,
  main_significant = BH_q < 0.05
)]
robust_gene <- gene_no_low[, .(
  feature, exclude5_effect = effect_Q4_minus_Q1, exclude5_p = p_value,
  exclude5_BH_q = BH_q, exclude5_significant = BH_q < 0.05
)]
gene_robustness <- merge(main_gene, robust_gene, by = "feature")
gene_robustness <- gene_robustness[main_significant | exclude5_significant]
gene_robustness[, minimum_BH_q := pmin(main_BH_q, exclude5_BH_q)]
setorder(gene_robustness, minimum_BH_q, feature)
fwrite(gene_robustness, file.path(table_dir, "Q4_Q1_FDR_gene_outlier_sensitivity.csv"))

rank_stats_main <- rank_stats
rank_stats <- setNames(gene_no_low$moderated_t, gene_no_low$feature)
gsea_no_low <- rbindlist(lapply(names(pathway_lists), function(collection) run_gsea(pathway_lists[[collection]], collection)), use.names = TRUE)
rank_stats <- rank_stats_main
gsea_no_low_out <- copy(gsea_no_low)
gsea_no_low_out[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
gsea_no_low_out[, abs_NES_sort := abs(NES)]
setorder(gsea_no_low_out, padj, -abs_NES_sort)
gsea_no_low_out[, abs_NES_sort := NULL]
fwrite(gsea_no_low_out, file.path(table_dir, "GSEA_Hallmark_Reactome_GOBP_adjusted_Q4_vs_Q1_exclude5_Tukey_low.csv"))
fwrite(
  gsea_no_low_out[is.na(NES) | is.na(padj), .(collection, pathway, reason = "fgsea unbalanced statistic; NES/p/q not estimable")],
  file.path(table_dir, "GSEA_unbalanced_NA_audit_exclude5_Tukey_low.csv")
)
gsea_no_low_summary <- gsea_no_low[, .(
  tested = .N,
  positive_BH_q_lt_0_05 = sum(padj < 0.05 & NES > 0, na.rm = TRUE),
  negative_BH_q_lt_0_05 = sum(padj < 0.05 & NES < 0, na.rm = TRUE),
  NA_pathways = sum(is.na(NES) | is.na(padj))
), by = collection]
fwrite(gsea_no_low_summary, file.path(table_dir, "GSEA_collection_summary_exclude5_Tukey_low.csv"))
focus_no_low <- gsea_no_low[pathway %in% focus_ids]
stopifnot(uniqueN(focus_no_low$pathway) == length(focus_ids))

focus_comparison <- rbind(
  copy(focus)[, analysis := "All Q1/Q4 (n=25/25)"],
  copy(focus_no_low)[, analysis := "Exclude five Tukey-low values (n=20/25)"]
)
focus_comparison_out <- copy(focus_comparison)
focus_comparison_out[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
fwrite(focus_comparison_out, file.path(table_dir, "GSEA_IFNG_CD28_IL12_focus_outlier_sensitivity.csv"))

deconv_no_low_fit <- fit_matrix(
  deconv_mat[, no_low_ids, drop = FALSE], design_no_low,
  "Adjusted Q4 vs Q1; exclude five Tukey-low values"
)$table
deconv_no_low_fit <- merge(deconv_no_low_fit, feature_meta, by.x = "feature", by.y = "feature_id", all.x = TRUE)
deconv_no_low_fit[, BH_q_all_cell_estimates := p.adjust(p_value, method = "BH")]
fwrite(deconv_no_low_fit, file.path(table_dir, "deconvolution_adjusted_Q4_vs_Q1_exclude5_Tukey_low.csv"))

cd8_no_low_fit <- fit_matrix(
  cd8_mat[, no_low_ids, drop = FALSE], design_no_low,
  "Adjusted Q4 vs Q1; exclude five Tukey-low values"
)$table
fwrite(cd8_no_low_fit, file.path(table_dir, "CD8_consensus_adjusted_Q4_vs_Q1_exclude5_Tukey_low.csv"))

model_meta_no_low <- model_meta[sample_id %in% no_low_ids]
model_meta_no_low <- model_meta_no_low[match(no_low_ids, sample_id)]
program_effects_no_low <- rbindlist(lapply(names(program_formulas), function(label) {
  d <- model.matrix(program_formulas[[label]], data = model_meta_no_low)
  if (qr(d)$rank != ncol(d)) stop("Rank-deficient no-low program model: ", label)
  x <- fit_matrix(program_mat[, no_low_ids, drop = FALSE], d, paste0(label, "; exclude five Tukey-low values"))$table
  x[, BH_q := p.adjust(p_value, "BH")]
  x[, adjustment := label]
  x
}))
program_effects_no_low <- merge(program_effects_no_low, program_meta, by.x = "feature", by.y = "program_id", all.x = TRUE)
fwrite(program_effects_no_low, file.path(table_dir, "immune_program_Q4_vs_Q1_abundance_adjusted_effects_exclude5_Tukey_low.csv"))

tf_no_low_fit <- fit_matrix(
  tf_mat[, no_low_ids, drop = FALSE], design_no_low,
  "Adjusted Q4 vs Q1; exclude five Tukey-low values"
)$table
tf_no_low_sd <- apply(tf_mat[, no_low_ids, drop = FALSE], 1L, sd)
tf_no_low_fit[, `:=`(
  standardized_effect = effect_Q4_minus_Q1 / tf_no_low_sd[feature],
  standardized_CI_low = CI_low / tf_no_low_sd[feature],
  standardized_CI_high = CI_high / tf_no_low_sd[feature],
  BH_q = p.adjust(p_value, "BH")
)]
tf_no_low_fit <- merge(tf_no_low_fit, tf_meta, by.x = "feature", by.y = "TF", all.x = TRUE)
setnames(tf_no_low_fit, "feature", "TF")
tf_no_low_fit[, abs_effect_sort := abs(standardized_effect)]
setorder(tf_no_low_fit, BH_q, -abs_effect_sort)
tf_no_low_fit[, abs_effect_sort := NULL]
fwrite(tf_no_low_fit, file.path(table_dir, "DoRothEA_TF_activity_adjusted_Q4_vs_Q1_exclude5_Tukey_low.csv"))
fwrite(
  tf_no_low_fit[CD28_pathway_targets_n > 0 | direct_CD28_gene_regulator == TRUE],
  file.path(table_dir, "CD28_linked_TF_activity_adjusted_Q4_vs_Q1_exclude5_Tukey_low.csv")
)

# Reader-facing figures.
COL_LOW <- "#1B9E8F"; COL_HIGH <- "#C43C3C"; COL_NS <- "#8F8F8F"
theme_pub <- theme_classic(base_size = 13) + theme(plot.title = element_text(face = "bold", size = 17), plot.subtitle = element_text(colour = "#444444"), strip.text = element_text(face = "bold"), legend.position = "bottom")

focus_plot <- copy(focus)
focus_plot[, label := gsub("_", " ", gsub("^(HALLMARK_|REACTOME_|GOBP_)", "", pathway))]
focus_plot[, plot_NES := fifelse(is.finite(NES), NES, 0)]
focus_plot[, significance := fifelse(!is.finite(padj), "Not estimable", fifelse(padj < 0.05, "BH q<0.05", "BH q>=0.05"))]
focus_plot[, plot_direction := fifelse(!is.finite(NES), "Not estimable", direction)]
focus_plot[, point_size := fifelse(is.finite(padj), -log10(pmax(padj, 1e-300)), 1)]
p_focus <- ggplot(focus_plot, aes(plot_NES, reorder(label, plot_NES), colour = plot_direction, shape = significance)) +
  geom_vline(xintercept = 0, colour = "#777777") + geom_point(aes(size = point_size), stroke = 1.1) +
  scale_colour_manual(values = c("Sarcosine-Low (Q1)" = COL_LOW, "Sarcosine-High (Q4)" = COL_HIGH, "Not estimable" = COL_NS)) +
  scale_shape_manual(values = c("BH q<0.05" = 16, "BH q>=0.05" = 1, "Not estimable" = 4)) +
  scale_size_continuous(range = c(2.5, 7)) +
  labs(title = "Immune-pathway enrichment in extreme tissue-Sarcosine groups", subtitle = "Top quartile versus bottom quartile; n=25 per group; positive NES = Sarcosine-High (Q4)", x = "Normalized enrichment score (Q4 vs Q1)", y = NULL, colour = NULL, shape = NULL, size = expression(-log[10](BH~q))) +
  guides(size = guide_legend(order = 1), shape = guide_legend(order = 2), colour = guide_legend(order = 3)) +
  theme_pub + theme(legend.box = "vertical")
ggsave(file.path(figure_dir, "Fig_Q4_vs_Q1_IFNG_CD28_IL12_Focused_GSEA.png"), p_focus, width = 14, height = 10, dpi = 400)

gsea_top <- gsea[is.finite(NES) & is.finite(padj)]
gsea_top[, rank_key := frank(padj, ties.method = "first"), by = .(collection, NES > 0)]
gsea_top <- gsea_top[rank_key <= 10]
gsea_top[, label := gsub("_", " ", gsub("^(HALLMARK_|REACTOME_|GOBP_)", "", pathway))]
gsea_top[, significance := fifelse(padj < 0.05, "BH q<0.05", "BH q>=0.05")]
p_gsea <- ggplot(gsea_top, aes(NES, reorder(label, NES), colour = direction, shape = significance)) +
  geom_vline(xintercept = 0, colour = "#777777") + geom_point(size = 3) + facet_wrap(~ collection, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c("Sarcosine-Low (Q1)" = COL_LOW, "Sarcosine-High (Q4)" = COL_HIGH)) + scale_shape_manual(values = c("BH q<0.05" = 16, "BH q>=0.05" = 1)) +
  labs(title = "Hallmark, Reactome and GO:BP extreme-group GSEA", subtitle = "Covariate-adjusted Q4 versus Q1; ten lowest-q pathways per direction and collection", x = "Normalized enrichment score (Q4 vs Q1)", y = NULL, colour = NULL, shape = NULL) + theme_pub
ggsave(file.path(figure_dir, "Fig_Q4_vs_Q1_ThreeCollection_GSEA.png"), p_gsea, width = 14, height = 18, dpi = 400)

cd8_plot <- merge(cd8[, .(sample_id, CD8_consensus_z, CD8_consensus_no_CIBERSORT_z)], meta[, .(sample_id, extreme_group)], by = "sample_id")
cd8_plot <- melt(cd8_plot, id.vars = c("sample_id", "extreme_group"), variable.name = "consensus", value.name = "z")
cd8_plot[, consensus := factor(consensus, levels = c("CD8_consensus_z", "CD8_consensus_no_CIBERSORT_z"), labels = c("Six-method consensus", "Excluding CIBERSORT"))]
p_cd8 <- ggplot(cd8_plot, aes(extreme_group, z, colour = extreme_group)) + geom_boxplot(outlier.shape = NA, colour = "black", fill = "white", width = 0.55) + geom_jitter(width = 0.12, alpha = 0.55, size = 1.7) + facet_wrap(~ consensus) + scale_colour_manual(values = c("Sarcosine-Low (Q1)" = COL_LOW, "Sarcosine-High (Q4)" = COL_HIGH)) +
  labs(title = "Cross-method CD8 T-cell estimates in extreme Sarcosine groups", subtitle = sprintf("Adjusted six-method effect %.2f SD (P=%.3g); CIBERSORT-excluded %.2f SD (P=%.3g)", cd8_fit[feature == "CD8 consensus (six methods)", effect_Q4_minus_Q1], cd8_fit[feature == "CD8 consensus (six methods)", p_value], cd8_fit[feature == "CD8 consensus (excluding CIBERSORT)", effect_Q4_minus_Q1], cd8_fit[feature == "CD8 consensus (excluding CIBERSORT)", p_value]), x = NULL, y = "CD8 consensus (z)", colour = NULL) + theme_pub + theme(axis.text.x = element_text(angle = 15, hjust = 1))
ggsave(file.path(figure_dir, "Fig_Q4_vs_Q1_CD8_Consensus.png"), p_cd8, width = 12, height = 7.5, dpi = 400)

program_plot <- program_effects[display_order <= 10]
program_plot[, significance := fifelse(BH_q < 0.05, "BH q<0.05", "BH q>=0.05")]
p_program <- ggplot(program_plot, aes(effect_Q4_minus_Q1, reorder(program_label, effect_Q4_minus_Q1))) + geom_vline(xintercept = 0, colour = "#777777") + geom_errorbar(aes(xmin = CI_low, xmax = CI_high, colour = significance), width = 0, orientation = "y") + geom_point(aes(colour = significance), size = 2.4) + facet_wrap(~ adjustment, ncol = 2) + scale_colour_manual(values = c("BH q<0.05" = COL_HIGH, "BH q>=0.05" = COL_NS)) + labs(title = "Extreme-group immune programs after T-cell/TME adjustment", subtitle = "Positive effect = higher in Sarcosine-High (Q4); n=25 versus 25", x = "Adjusted Q4 - Q1 difference (program-score SD; 95% CI)", y = NULL, colour = NULL) + theme_pub
ggsave(file.path(figure_dir, "Fig_Q4_vs_Q1_Immune_Programs_Tcell_TME_Adjustment.png"), p_program, width = 15, height = 12.5, dpi = 400)

tf_display <- unique(rbind(tf_linked[order(BH_q, -abs(standardized_effect))][1:min(.N, 20L)], tf_fit[TF %in% c("NFKB1", "RELA", "STAT1", "JUN", "FOS", "EGR1", "SP1")]), by = "TF")
tf_display[, significance := fifelse(BH_q < 0.05, "BH q<0.05", "BH q>=0.05")]
fwrite(tf_display, file.path(table_dir, "Figure_Source_CD28_candidate_TFs_Q4_vs_Q1.csv"))
p_tf <- ggplot(tf_display, aes(standardized_effect, reorder(TF, standardized_effect))) + geom_vline(xintercept = 0, colour = "#777777") + geom_errorbar(aes(xmin = standardized_CI_low, xmax = standardized_CI_high, colour = significance), width = 0, orientation = "y") + geom_point(aes(colour = significance, size = pmax(CD28_pathway_targets_n, 1L))) + scale_colour_manual(values = c("BH q<0.05" = COL_HIGH, "BH q>=0.05" = COL_NS)) + scale_size_continuous(range = c(2.5, 6)) + labs(title = "Candidate TF programs linked to the CD28 pathway", subtitle = "Signed DoRothEA A/B activity; adjusted Q4 versus Q1; positive = Sarcosine-High", x = "Adjusted Q4 - Q1 TF-activity difference (SD; 95% CI)", y = NULL, colour = NULL, size = "CD28-pathway\ntargets") + theme_pub
ggsave(file.path(figure_dir, "Fig_Q4_vs_Q1_CD28_Linked_TF_Activity.png"), p_tf, width = 10, height = 10, dpi = 400)

focus_sens_plot <- copy(focus_comparison)
focus_sens_plot[, label := gsub("_", " ", gsub("^(HALLMARK_|REACTOME_|GOBP_)", "", pathway))]
focus_sens_plot[, plot_NES := fifelse(is.finite(NES), NES, 0)]
focus_sens_plot[, significance := fifelse(!is.finite(padj), "Not estimable", fifelse(padj < 0.05, "BH q<0.05", "BH q>=0.05"))]
focus_sens_plot[, plot_direction := fifelse(!is.finite(NES), "Not estimable", direction)]
p_focus_sens <- ggplot(focus_sens_plot, aes(plot_NES, reorder(label, plot_NES), colour = plot_direction, shape = significance)) +
  geom_vline(xintercept = 0, colour = "#777777") + geom_point(size = 3.1, stroke = 1.1) +
  facet_wrap(~ analysis, ncol = 2) +
  scale_colour_manual(values = c("Sarcosine-Low (Q1)" = COL_LOW, "Sarcosine-High (Q4)" = COL_HIGH, "Not estimable" = COL_NS)) +
  scale_shape_manual(values = c("BH q<0.05" = 16, "BH q>=0.05" = 1, "Not estimable" = 4)) +
  labs(
    title = "Focused immune GSEA: influence of five Tukey-low Sarcosine values",
    subtitle = "Fixed original quartile cutpoints; positive NES = Sarcosine-High (Q4)",
    x = "Normalized enrichment score (Q4 vs Q1)", y = NULL, colour = NULL, shape = NULL
  ) + theme_pub
ggsave(file.path(figure_dir, "Fig_Q4_vs_Q1_Focused_GSEA_Outlier_Sensitivity.png"), p_focus_sens, width = 17, height = 10, dpi = 400)

fwrite(data.table(input = inputs, md5 = unname(tools::md5sum(inputs)), size_bytes = file.info(inputs)$size), file.path(log_dir, "input_manifest.csv"))
validation <- data.table(
  check = c("Q1 n=25", "Q4 n=25", "middle n=50", "strict cutpoint separation", ">5000 filtered genes", "four full-rank gene designs", "three collections returned", "all focused pathways returned", "deconvolution finite", "CD8 finite", "four full-rank program designs", "TF activities finite", "22 CD28-linked TFs retained", "all five Tukey-low values are in Q1", "outlier sensitivity n=20/25", "outlier sensitivity full-rank design", "all focused pathways returned after exclusion", "all six figures exist"),
  passed = c(
    sum(meta_all$extreme_group == "Sarcosine-Low (Q1)") == 25L,
    sum(meta_all$extreme_group == "Sarcosine-High (Q4)") == 25L,
    sum(meta_all$extreme_group == "Middle 50% excluded") == 50L,
    max(meta[extreme_group == "Sarcosine-Low (Q1)", sarcosine_normalized_intensity]) < cut_q1 && min(meta[extreme_group == "Sarcosine-High (Q4)", sarcosine_normalized_intensity]) > cut_q4,
    nrow(expr) == 15119L, all(vapply(designs, function(d) qr(d)$rank == ncol(d), logical(1))),
    setequal(unique(gsea$collection), c("Hallmark", "Reactome", "GO:BP")), uniqueN(focus$pathway) == length(focus_ids),
    !anyNA(deconv_mat) && all(is.finite(deconv_mat)), !anyNA(cd8_mat) && all(is.finite(cd8_mat)),
    all(vapply(program_formulas, function(f) { d <- model.matrix(f, data = model_meta); qr(d)$rank == ncol(d) }, logical(1))),
    !anyNA(tf_mat) && all(is.finite(tf_mat)), nrow(tf_linked) == 22L,
    all(tukey_low_ids %in% meta[extreme_group == "Sarcosine-Low (Q1)", sample_id]),
    nrow(meta_no_low) == 45L && sum(meta_no_low$extreme_group == "Sarcosine-Low (Q1)") == 20L && sum(meta_no_low$extreme_group == "Sarcosine-High (Q4)") == 25L,
    qr(design_no_low)$rank == ncol(design_no_low), uniqueN(focus_no_low$pathway) == length(focus_ids),
    all(file.exists(file.path(figure_dir, c("Fig_Q4_vs_Q1_IFNG_CD28_IL12_Focused_GSEA.png", "Fig_Q4_vs_Q1_ThreeCollection_GSEA.png", "Fig_Q4_vs_Q1_CD8_Consensus.png", "Fig_Q4_vs_Q1_Immune_Programs_Tcell_TME_Adjustment.png", "Fig_Q4_vs_Q1_CD28_Linked_TF_Activity.png", "Fig_Q4_vs_Q1_Focused_GSEA_Outlier_Sensitivity.png"))))
  )
)
fwrite(validation, file.path(log_dir, "validation_checks.csv"))
if (!all(validation$passed)) stop("Q4/Q1 validation failed")
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))

cat("Q4 versus Q1 analysis complete: ", out_root, "\n", sep = "")
cat(sprintf("Cutpoints: Q1 <= %.10f; Q4 >= %.10f; n=25/25\n", cut_q1, cut_q4))
print(gene_summary); print(gsea_summary)
print(focus[, .(collection, pathway, NES, pval, padj, direction)])
print(cd8_fit[, .(feature, effect_Q4_minus_Q1, CI_low, CI_high, p_value, BH_q)])
print(program_effects[adjustment == "Q4 vs Q1", .(program_label, effect_Q4_minus_Q1, CI_low, CI_high, p_value, BH_q)])
print(tf_linked[, .(TF, standardized_effect, standardized_CI_low, standardized_CI_high, p_value, BH_q, CD28_pathway_targets_n)])
cat("\nExclude-five Tukey-low sensitivity:\n")
print(gene_summary_no_low)
print(gsea_no_low_summary)
print(focus_no_low[, .(collection, pathway, NES, pval, padj, direction)])
print(cd8_no_low_fit[, .(feature, effect_Q4_minus_Q1, CI_low, CI_high, p_value, BH_q)])
