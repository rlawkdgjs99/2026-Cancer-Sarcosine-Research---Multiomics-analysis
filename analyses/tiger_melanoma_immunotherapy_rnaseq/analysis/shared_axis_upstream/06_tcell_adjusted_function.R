# TIGER melanoma PRE73: distinguish T-cell abundance from functional state.
#
# Scientific question
# -------------------
# Is the Degradation-High versus -Low functional-program difference explained
# only by inferred T-cell abundance, or does an association remain after
# adjusting for independently generated bulk-RNA deconvolution estimates?
#
# Reader-facing figures show the biological comparison (Degradation High/Low).
# Clinical covariates and abundance estimates are model adjustments, not groups.
#
# Limits
# ------
# Both pathway scores and deconvolution estimates derive from the same bulk RNA.
# Therefore abundance-adjusted persistence is sensitivity evidence compatible
# with a cell-state difference, not proof of a T-cell-intrinsic mechanism.

options(stringsAsFactors = FALSE, width = 180)
set.seed(42)

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

required_packages <- c("data.table", "limma", "msigdbr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
})

used_data_root <- normalizePath(file.path(analysis_root, "..", "사용데이터_모음"))
expression_path <- file.path(
  used_data_root, "Source_Input", "TIGER_PRJEB23709_expression_FPKM_gene_by_sample.csv"
)
metadata_path <- file.path(
  used_data_root, "Analysis_Ready", "TIGER_PRE73_Fig4bc_analysis_data.csv"
)
deconv_path <- file.path(
  analysis_root, "results", "TIGER_PRE73_Degradation_HighLow_Deconvolution_7methods_26.09.02",
  "tables", "05_core_celltype_values_long.csv"
)
tf_activity_path <- file.path(
  analysis_root, "results", "TIGER_PRE73_CD28_Upstream_TF_26.09.02",
  "tables", "09_sample_level_TF_ULM_activities_long.csv"
)
tf_primary_path <- file.path(
  analysis_root, "results", "TIGER_PRE73_CD28_Upstream_TF_26.09.02",
  "tables", "07_TF_ULM_activity_High_vs_Low.csv"
)
required_inputs <- c(expression_path, metadata_path, deconv_path, tf_activity_path, tf_primary_path)
for (path in required_inputs) if (!file.exists(path)) stop("Missing input: ", path)

output_root <- file.path(
  analysis_root, "results", "TIGER_PRE73_Tcell_Abundance_Adjusted_Function_26.09.02"
)
table_dir <- file.path(output_root, "tables")
figure_dir <- file.path(output_root, "figures")
log_dir <- file.path(output_root, "logs")
for (d in c(table_dir, figure_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

sha256_file <- function(path) {
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  strsplit(out, "[[:space:]]+")[[1]][1]
}

# -------------------------------------------------------------------------
# Input and design reconstruction
# -------------------------------------------------------------------------
expr_dt <- fread(expression_path, check.names = FALSE)
meta <- fread(metadata_path, check.names = FALSE)
deconv <- fread(deconv_path, check.names = FALSE)
tf_activity <- fread(tf_activity_path, check.names = FALSE)
tf_primary <- fread(tf_primary_path, check.names = FALSE)

stopifnot(
  nrow(meta) == 73L,
  uniqueN(meta$sample_id) == 73L,
  all(meta$timepoint == "PRE"),
  all(c("sample_id", "therapy_short", "response_group", "age", "gender", "Degradation_score") %in% names(meta)),
  all(c("method", "analysis_feature", "sample_id", "analysis_value") %in% names(deconv)),
  all(c("sample_id", "TF", "ULM_activity") %in% names(tf_activity))
)

degradation_median <- median(meta$Degradation_score)
meta[, degradation_group := factor(
  ifelse(Degradation_score > degradation_median, "High", "Low"),
  levels = c("Low", "High")
)]
stopifnot(sum(meta$degradation_group == "Low") == 37L, sum(meta$degradation_group == "High") == 36L)

meta[, therapy_short := relevel(factor(therapy_short), ref = "antiPD1")]
meta[, response_group := relevel(factor(response_group), ref = "NR")]
meta[, gender := relevel(factor(gender), ref = "Male")]
meta[, age_z := as.numeric(scale(age))]

gene_symbols <- trimws(expr_dt[[1]])
stopifnot(!anyNA(gene_symbols), all(nzchar(gene_symbols)), !anyDuplicated(gene_symbols))
stopifnot(all(meta$sample_id %in% names(expr_dt)))
sample_ids <- meta$sample_id
fpkm <- as.matrix(expr_dt[, ..sample_ids])
storage.mode(fpkm) <- "double"
rownames(fpkm) <- gene_symbols
colnames(fpkm) <- sample_ids
stopifnot(ncol(fpkm) == 73L, !anyNA(fpkm), all(is.finite(fpkm)), all(fpkm >= 0))

minimum_samples <- ceiling(0.10 * ncol(fpkm))
keep_expression <- rowSums(fpkm >= 1) >= minimum_samples
log_expression <- log2(fpkm[keep_expression, , drop = FALSE] + 1)
log_expression <- log_expression[apply(log_expression, 1L, var) > 0, , drop = FALSE]

# -------------------------------------------------------------------------
# Predeclared functional programs
# -------------------------------------------------------------------------
program_catalog <- data.table(
  program_id = c(
    "GOBP_T_CELL_PROLIFERATION",
    "GOBP_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE",
    "GOBP_T_CELL_MEDIATED_CYTOTOXICITY",
    "REACTOME_CO_STIMULATION_BY_CD28",
    "REACTOME_CD28_DEPENDENT_PI3K_AKT_SIGNALING",
    "DOROTHEA_NFKB1_ACTIVITY",
    "DOROTHEA_RELA_ACTIVITY"
  ),
  program_label = c(
    "T-cell proliferation",
    "Effector-associated T-cell differentiation",
    "T-cell-mediated cytotoxicity",
    "CD28 co-stimulation",
    "CD28–PI3K/AKT signaling",
    "NFKB1 regulon activity",
    "RELA regulon activity"
  ),
  program_family = c(
    "Proliferation", "Differentiation", "Effector function",
    "CD28 signaling", "CD28 signaling", "NF-κB activity", "NF-κB activity"
  ),
  source = c("GO:BP", "GO:BP", "GO:BP", "Reactome", "Reactome", "DoRothEA A/B ULM", "DoRothEA A/B ULM"),
  program_order = 1:7
)

go_bp <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"))
reactome <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))
stopifnot(uniqueN(c(go_bp$db_version, reactome$db_version)) == 1L)
msigdb_version <- unique(c(go_bp$db_version, reactome$db_version))

membership <- rbindlist(list(
  go_bp[gs_name %in% program_catalog[source == "GO:BP", program_id],
        .(program_id = gs_name, gene_symbol, MSigDB_version = db_version, exact_source = gs_exact_source, url = gs_url)],
  reactome[gs_name %in% program_catalog[source == "Reactome", program_id],
           .(program_id = gs_name, gene_symbol, MSigDB_version = db_version, exact_source = gs_exact_source, url = gs_url)]
), use.names = TRUE)
membership <- unique(membership, by = c("program_id", "gene_symbol"))
stopifnot(setequal(unique(membership$program_id), program_catalog[source %in% c("GO:BP", "Reactome"), program_id]))
membership[, expressed_in_PRE73 := gene_symbol %in% rownames(log_expression)]

score_one <- function(genes) {
  genes <- intersect(unique(genes), rownames(log_expression))
  if (length(genes) < 10L) stop("Fewer than 10 expressed genes in a program")
  z <- t(scale(t(log_expression[genes, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  score <- colMeans(z)
  as.numeric(scale(score))
}

set_score_ids <- program_catalog[source %in% c("GO:BP", "Reactome"), program_id]
set_score_mat <- do.call(rbind, lapply(set_score_ids, function(pid) {
  score_one(membership[program_id == pid, gene_symbol])
}))
rownames(set_score_mat) <- set_score_ids
colnames(set_score_mat) <- sample_ids

tf_subset <- tf_activity[TF %in% c("NFKB1", "RELA")]
stopifnot(nrow(tf_subset) == 2L * nrow(meta), uniqueN(tf_subset$sample_id) == nrow(meta))
tf_wide <- dcast(tf_subset, TF ~ sample_id, value.var = "ULM_activity")
tf_mat <- as.matrix(tf_wide[, -"TF"])
rownames(tf_mat) <- paste0("DOROTHEA_", tf_wide$TF, "_ACTIVITY")
tf_mat <- tf_mat[, sample_ids, drop = FALSE]
tf_mat <- t(scale(t(tf_mat)))

program_mat <- rbind(set_score_mat, tf_mat)
program_mat <- program_mat[program_catalog$program_id, , drop = FALSE]
stopifnot(!anyNA(program_mat), all(is.finite(program_mat)))

program_manifest <- merge(
  program_catalog,
  membership[expressed_in_PRE73 == TRUE, .(expression_filtered_genes_n = uniqueN(gene_symbol)), by = program_id],
  by = "program_id", all.x = TRUE
)
program_manifest[source == "DoRothEA A/B ULM", expression_filtered_genes_n := NA_integer_]
program_manifest[, score_definition := fifelse(
  source %in% c("GO:BP", "Reactome"),
  "Mean of gene-wise z-scores across expression-filtered set members; standardized across 73 tumors",
  "Signed DoRothEA A/B ULM TF-activity score; standardized across 73 tumors"
)]
setorder(program_manifest, program_order)

program_long <- as.data.table(as.table(program_mat))
setnames(program_long, c("program_id", "sample_id", "program_score"))
program_long <- merge(program_long, program_catalog, by = "program_id", all.x = TRUE)
program_long <- merge(
  program_long,
  meta[, .(sample_id, degradation_group, Degradation_score)],
  by = "sample_id", all.x = TRUE
)

# -------------------------------------------------------------------------
# CD8 and pan-T abundance adjusters: one at a time to avoid collinearity.
# -------------------------------------------------------------------------
cd8_long <- deconv[analysis_feature == "CD8 T cells" & method %in% c(
  "CIBERSORT", "EPIC", "MCPcounter", "quanTIseq", "TIMER", "xCell"
), .(sample_id, abundance_method = method, abundance_z = analysis_value)]
pan_long <- deconv[method == "MCPcounter" & analysis_feature == "Pan T cells",
                   .(sample_id, abundance_method = "MCPcounter Pan T", abundance_z = analysis_value)]
stopifnot(nrow(cd8_long) == 6L * nrow(meta), nrow(pan_long) == nrow(meta))

cd8_wide <- dcast(cd8_long, sample_id ~ abundance_method, value.var = "abundance_z")
cd8_matrix <- as.matrix(cd8_wide[, -"sample_id"])
consensus <- apply(cd8_matrix, 1L, median)
consensus <- as.numeric(scale(consensus))
consensus_long <- data.table(
  sample_id = cd8_wide$sample_id,
  abundance_method = "CD8 consensus",
  abundance_z = consensus
)
abundance_long <- rbindlist(list(cd8_long, consensus_long, pan_long), use.names = TRUE)
abundance_long[, abundance_method := factor(
  abundance_method,
  levels = c("CIBERSORT", "EPIC", "MCPcounter", "quanTIseq", "TIMER", "xCell", "CD8 consensus", "MCPcounter Pan T")
)]
stopifnot(nrow(abundance_long) == 8L * nrow(meta), !anyNA(abundance_long$abundance_z))

# -------------------------------------------------------------------------
# Limma models. Effects are High minus Low in outcome-SD units.
# -------------------------------------------------------------------------
fit_programs <- function(outcome_mat, meta, abundance = NULL, model_label) {
  model_meta <- copy(meta)
  if (is.null(abundance)) {
    design <- model.matrix(
      ~ therapy_short + response_group + age_z + gender + degradation_group,
      data = model_meta
    )
  } else {
    abundance <- abundance[match(model_meta$sample_id, sample_id)]
    stopifnot(identical(abundance$sample_id, model_meta$sample_id), !anyNA(abundance$abundance_z))
    model_meta[, abundance_z := abundance$abundance_z]
    design <- model.matrix(
      ~ therapy_short + response_group + age_z + gender + abundance_z + degradation_group,
      data = model_meta
    )
  }
  stopifnot(qr(design)$rank == ncol(design))
  coef_name <- "degradation_groupHigh"
  fit <- eBayes(lmFit(outcome_mat, design), robust = TRUE)
  tt <- as.data.table(topTable(fit, coef = coef_name, number = Inf, sort.by = "none", adjust.method = "BH"), keep.rownames = "program_id")
  se <- fit$stdev.unscaled[, coef_name] * sqrt(fit$s2.post)
  out <- tt[, .(
    program_id, effect_High_minus_Low = logFC, moderated_t = t,
    p_value = P.Value, BH_q_within_adjuster = adj.P.Val
  )]
  out[, standard_error := se[program_id]]
  df_total <- fit$df.total
  if (length(df_total) == 1L) {
    df_for_out <- rep(df_total, nrow(out))
  } else {
    df_for_out <- df_total[match(out$program_id, rownames(fit$coefficients))]
  }
  stopifnot(length(df_for_out) == nrow(out), !anyNA(df_for_out), all(df_for_out > 0))
  ci_multiplier <- qt(0.975, df = df_for_out)
  out[, `:=`(
    CI_low = effect_High_minus_Low - ci_multiplier * standard_error,
    CI_high = effect_High_minus_Low + ci_multiplier * standard_error,
    adjustment = model_label,
    design_rank = qr(design)$rank,
    design_columns_n = ncol(design)
  )]
  out
}

clinical_results <- fit_programs(program_mat, meta, abundance = NULL, model_label = "No abundance adjustment")
adjusted_results <- rbindlist(lapply(levels(abundance_long$abundance_method), function(method_name) {
  fit_programs(
    program_mat, meta,
    abundance = abundance_long[abundance_method == method_name, .(sample_id, abundance_z)],
    model_label = method_name
  )
}))
adjusted_results[, BH_q_global := p.adjust(p_value, method = "BH")]

all_program_results <- rbindlist(list(
  clinical_results[, BH_q_global := p.adjust(p_value, method = "BH")],
  adjusted_results
), use.names = TRUE)
all_program_results <- merge(all_program_results, program_catalog, by = "program_id", all.x = TRUE)
setorder(all_program_results, program_order, adjustment)

# Refit abundance outcomes themselves to answer whether cell abundance differs.
abundance_wide <- dcast(abundance_long, abundance_method ~ sample_id, value.var = "abundance_z")
abundance_mat <- as.matrix(abundance_wide[, -"abundance_method"])
rownames(abundance_mat) <- as.character(abundance_wide$abundance_method)
abundance_mat <- abundance_mat[, sample_ids, drop = FALSE]
abundance_results <- fit_programs(
  abundance_mat, meta, abundance = NULL, model_label = "Cell-abundance estimate"
)
setnames(abundance_results, "program_id", "abundance_method")
abundance_results[, BH_q_across_abundance_estimates := p.adjust(p_value, method = "BH")]

# Cross-method robustness summary, keeping individual CD8 methods separate.
individual_cd8_methods <- c("CIBERSORT", "EPIC", "MCPcounter", "quanTIseq", "TIMER", "xCell")
robustness <- adjusted_results[adjustment %in% individual_cd8_methods, .(
  CD8_methods_n = .N,
  positive_direction_n = sum(effect_High_minus_Low > 0),
  CI_above_zero_n = sum(CI_low > 0),
  nominal_p_below_0.05_n = sum(p_value < 0.05),
  median_adjusted_effect = median(effect_High_minus_Low),
  min_adjusted_effect = min(effect_High_minus_Low),
  max_adjusted_effect = max(effect_High_minus_Low)
), by = program_id]
consensus_result <- adjusted_results[adjustment == "CD8 consensus", .(
  program_id,
  consensus_effect = effect_High_minus_Low,
  consensus_CI_low = CI_low,
  consensus_CI_high = CI_high,
  consensus_p = p_value,
  consensus_BH_q = BH_q_within_adjuster
)]
pan_result <- adjusted_results[adjustment == "MCPcounter Pan T", .(
  program_id,
  PanT_effect = effect_High_minus_Low,
  PanT_CI_low = CI_low,
  PanT_CI_high = CI_high,
  PanT_p = p_value,
  PanT_BH_q = BH_q_within_adjuster
)]
robustness <- Reduce(function(x, y) merge(x, y, by = "program_id", all = TRUE), list(
  robustness, consensus_result, pan_result, program_catalog
))
robustness[, robustness_class := fifelse(
  positive_direction_n >= 5L & CI_above_zero_n >= 4L & consensus_BH_q < 0.05,
  "Persists after CD8-abundance adjustment",
  fifelse(
    positive_direction_n >= 5L,
    "Directionally persistent; statistical support incomplete",
    "Not robust to CD8-abundance adjustment"
  )
)]
setorder(robustness, program_order)

# Quantify outcome/adjuster dependence. A high correlation is expected because
# both quantities are derived from bulk RNA and cautions against interpreting a
# loss of significance after adjustment as definitive biological mediation.
program_abundance_correlations <- merge(
  program_long[, .(sample_id, program_id, program_label, program_score)],
  abundance_long[, .(sample_id, abundance_method, abundance_z)],
  by = "sample_id", allow.cartesian = TRUE
)[, .(
  Pearson_r = cor(program_score, abundance_z, method = "pearson"),
  Spearman_rho = cor(program_score, abundance_z, method = "spearman")
), by = .(program_id, program_label, abundance_method)]
setorder(program_abundance_correlations, program_id, abundance_method)

# Exact cross-check for NFKB1/RELA unadjusted standardized effects.
nfkb_check <- merge(
  clinical_results[program_id %in% c("DOROTHEA_NFKB1_ACTIVITY", "DOROTHEA_RELA_ACTIVITY"),
                   .(TF = sub("DOROTHEA_(.*)_ACTIVITY", "\\1", program_id), rerun_effect = effect_High_minus_Low)],
  tf_primary[TF %in% c("NFKB1", "RELA"), .(TF, prior_effect = standardized_effect)],
  by = "TF", all = TRUE
)
nfkb_check[, absolute_difference := abs(rerun_effect - prior_effect)]

# -------------------------------------------------------------------------
# Tables
# -------------------------------------------------------------------------
fwrite(program_manifest, file.path(table_dir, "01_functional_program_manifest.csv"))
fwrite(membership, file.path(table_dir, "02_MSigDB_program_gene_membership.csv"))
fwrite(program_long, file.path(table_dir, "03_sample_functional_program_scores_long.csv"))
fwrite(abundance_long, file.path(table_dir, "04_Tcell_abundance_adjusters_long.csv"))
fwrite(abundance_results, file.path(table_dir, "05_Tcell_abundance_High_vs_Low.csv"))
fwrite(clinical_results, file.path(table_dir, "06_functional_programs_High_vs_Low_before_abundance_adjustment.csv"))
fwrite(adjusted_results, file.path(table_dir, "07_functional_programs_after_each_Tcell_abundance_adjustment.csv"))
fwrite(all_program_results, file.path(table_dir, "08_all_functional_program_models.csv"))
fwrite(robustness, file.path(table_dir, "09_cross_method_robustness_summary.csv"))
fwrite(nfkb_check, file.path(table_dir, "10_NFKB_TF_activity_rerun_crosscheck.csv"))
fwrite(program_abundance_correlations, file.path(table_dir, "10_program_abundance_correlations.csv"))

# -------------------------------------------------------------------------
# Reader-facing figures: actual High/Low groups, not model names as groups.
# -------------------------------------------------------------------------
low_color <- "#1F9E89"
high_color <- "#CC3D3D"
neutral_color <- "#777777"

plot_program_order <- rev(program_catalog$program_label)
clinical_plot <- merge(clinical_results, program_catalog, by = "program_id", all.x = TRUE)
clinical_plot[, program_label := factor(program_label, levels = plot_program_order)]

p_before <- ggplot(clinical_plot, aes(y = program_label, x = effect_High_minus_Low)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#777777") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0.22, color = "#444444") +
  geom_point(aes(fill = BH_q_within_adjuster < 0.05), shape = 21, size = 4.2, color = "#222222") +
  scale_fill_manual(values = c(`TRUE` = high_color, `FALSE` = "white"), guide = "none") +
  labs(
    title = "Functional programs before T-cell-abundance adjustment",
    subtitle = "High > Low to the right; filled points: BH q < 0.05",
    x = "Adjusted High − Low difference (outcome SD)", y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 16), axis.text.y = element_text(face = "bold"))

focus_adjusted <- merge(
  adjusted_results[adjustment %in% c("CD8 consensus", "MCPcounter Pan T")],
  program_catalog, by = "program_id", all.x = TRUE
)
focus_adjusted[, program_label := factor(program_label, levels = plot_program_order)]
focus_adjusted[, adjustment := factor(
  adjustment,
  levels = c("CD8 consensus", "MCPcounter Pan T"),
  labels = c("Adjusted for CD8 consensus", "Adjusted for pan-T estimate")
)]
p_after <- ggplot(focus_adjusted, aes(y = program_label, x = effect_High_minus_Low, color = adjustment)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#777777") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0.18, position = position_dodge(width = 0.48)) +
  geom_point(aes(shape = BH_q_within_adjuster < 0.05), position = position_dodge(width = 0.48), size = 3.5) +
  scale_color_manual(values = c("#3366AA", "#AA4499"), name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), name = "BH q < 0.05") +
  labs(
    title = "Functional programs after T-cell-abundance adjustment",
    subtitle = "High > Low to the right; filled points: BH q < 0.05",
    x = "Abundance-adjusted High − Low difference (outcome SD)", y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
    legend.position = "bottom"
  )

abundance_labels <- c(
  CIBERSORT = "CIBERSORT CD8", EPIC = "EPIC CD8", MCPcounter = "MCP-counter CD8",
  quanTIseq = "quanTIseq CD8", TIMER = "TIMER CD8", xCell = "xCell CD8",
  `CD8 consensus` = "CD8 consensus", `MCPcounter Pan T` = "MCP-counter pan T"
)
abundance_plot <- copy(abundance_results)
abundance_plot[, abundance_label := factor(
  abundance_labels[abundance_method], levels = rev(unname(abundance_labels))
)]
p_abundance <- ggplot(abundance_plot, aes(y = abundance_label, x = effect_High_minus_Low)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#777777") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0.20, color = "#444444") +
  geom_point(aes(fill = BH_q_across_abundance_estimates < 0.05), shape = 21, size = 4.1, color = "#222222") +
  scale_fill_manual(values = c(`TRUE` = high_color, `FALSE` = "white"), guide = "none") +
  labs(
    title = "CD8/T-cell abundance estimates",
    subtitle = "High > Low to the right; filled points: BH q < 0.05",
    x = "Adjusted High − Low difference (outcome SD)", y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 16), axis.text.y = element_text(face = "bold"))

main_figure <- (p_before | p_after) / p_abundance +
  plot_layout(heights = c(1.2, 0.9), guides = "collect") +
  plot_annotation(
    title = "T-cell abundance and functional programs by sarcosine-degradation state",
    subtitle = "TIGER melanoma PRE-only: Degradation-Low n = 37; Degradation-High n = 36. Clinical covariates are adjusted in all estimates.",
    tag_levels = "a",
    theme = theme(plot.title = element_text(face = "bold", size = 21), plot.subtitle = element_text(size = 12))
  ) & theme(legend.position = "bottom")
ggsave(
  file.path(figure_dir, "Fig_Tcell_Abundance_and_Function_High_vs_Low.png"),
  main_figure, width = 16, height = 12.5, dpi = 400, bg = "white"
)
ggsave(
  file.path(figure_dir, "Fig_Tcell_Abundance_and_Function_High_vs_Low.pdf"),
  main_figure, width = 16, height = 12.5, device = cairo_pdf, bg = "white"
)

# Group-only boxplots provide the most immediate reader-level view.
box_dt <- copy(program_long)
box_dt[, program_label := factor(program_label, levels = program_catalog$program_label)]
box_dt[, degradation_group := factor(degradation_group, levels = c("Low", "High"))]
q_labels <- clinical_plot[, .(
  program_id, program_label,
  label = paste0("q = ", formatC(BH_q_within_adjuster, format = "g", digits = 2))
)]
box_y <- box_dt[, .(y = max(program_score) + 0.12 * diff(range(program_score))), by = .(program_id, program_label)]
q_labels <- merge(q_labels, box_y, by = c("program_id", "program_label"), all.x = TRUE)
p_box <- ggplot(box_dt, aes(x = degradation_group, y = program_score, color = degradation_group)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, color = "#222222", fill = "white") +
  geom_jitter(width = 0.12, height = 0, size = 1.3, alpha = 0.72) +
  geom_text(data = q_labels, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, color = "#333333", size = 3.1) +
  facet_wrap(~ program_label, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c(Low = low_color, High = high_color), guide = "none") +
  labs(
    title = "T-cell functional-program scores in Degradation-High and -Low tumors",
    subtitle = "Points are PRE-treatment tumors; q values are covariate-adjusted High-versus-Low tests",
    x = "Sarcosine-degradation group", y = "Standardized program score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 19),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(face = "bold")
  )
ggsave(
  file.path(figure_dir, "Fig_Tcell_Function_Scores_High_vs_Low.png"),
  p_box, width = 15, height = 8.8, dpi = 400, bg = "white"
)
ggsave(
  file.path(figure_dir, "Fig_Tcell_Function_Scores_High_vs_Low.pdf"),
  p_box, width = 15, height = 8.8, device = cairo_pdf, bg = "white"
)

# Full six-method CD8 sensitivity figure.
sensitivity_plot <- merge(
  adjusted_results[adjustment %in% individual_cd8_methods],
  program_catalog, by = "program_id", all.x = TRUE
)
sensitivity_plot[, program_label := factor(program_label, levels = plot_program_order)]
sensitivity_plot[, adjustment := factor(adjustment, levels = individual_cd8_methods)]
p_sensitivity <- ggplot(sensitivity_plot, aes(y = program_label, x = effect_High_minus_Low)) +
  geom_vline(xintercept = 0, linetype = 2, color = "#777777") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0.16, color = "#555555") +
  geom_point(aes(fill = BH_q_within_adjuster < 0.05), shape = 21, size = 3.2, color = "#222222") +
  facet_wrap(~ adjustment, ncol = 3) +
  scale_fill_manual(values = c(`TRUE` = high_color, `FALSE` = "white"), name = "BH q < 0.05") +
  labs(
    title = "Sensitivity across six CD8 T-cell abundance estimates",
    subtitle = "Every panel compares the same Degradation-High versus -Low tumors after adjusting one CD8 estimate at a time",
    x = "CD8-abundance-adjusted High − Low difference (outcome SD)", y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 19), strip.text = element_text(face = "bold"))
ggsave(
  file.path(figure_dir, "Fig_Tcell_Function_CD8_Adjustment_Sensitivity.png"),
  p_sensitivity, width = 15, height = 10.5, dpi = 400, bg = "white"
)
ggsave(
  file.path(figure_dir, "Fig_Tcell_Function_CD8_Adjustment_Sensitivity.pdf"),
  p_sensitivity, width = 15, height = 10.5, device = cairo_pdf, bg = "white"
)

# -------------------------------------------------------------------------
# Reproducibility records and validation
# -------------------------------------------------------------------------
input_manifest <- data.table(
  role = c("FPKM expression", "PRE73 metadata", "deconvolution estimates", "sample TF activities", "prior TF model cross-check"),
  path = normalizePath(required_inputs),
  sha256 = vapply(required_inputs, sha256_file, character(1)),
  size_bytes = file.info(required_inputs)$size
)
fwrite(input_manifest, file.path(table_dir, "00_input_manifest.csv"))

method_contract <- data.table(
  field = c(
    "analysis_unit", "population", "group_definition", "estimand", "clinical_covariates",
    "abundance_adjustment", "expression_space", "expression_filter", "program_scoring",
    "TF_activity", "multiplicity", "robustness_rule", "figure_grouping", "interpretation_limit"
  ),
  value = c(
    "one PRE-treatment tumor per patient", "TIGER melanoma PRE-only n=73",
    "median split of z-scored SARDH/PIPOX degradation module: High n=36, Low n=37",
    "functional score difference High minus Low in outcome-SD units",
    "therapy + clinical response + age + sex",
    "six CD8 methods one at a time; median-z CD8 consensus and MCPcounter pan-T as focused summaries",
    "log2(FPKM+1)", "FPKM >=1 in at least ceiling(10% of 73) samples; non-zero variance",
    paste0("gene-wise z-score mean from predeclared MSigDB ", msigdb_version, " GO:BP/Reactome sets, standardized by program"),
    "signed DoRothEA A/B ULM scores from analysis 110, standardized across samples",
    "BH within seven programs for each adjuster; global BH across all abundance-adjusted program tests also reported",
    "positive direction in >=5/6 CD8 methods, CI>0 in >=4/6, and consensus BH q<0.05",
    "Sarcosine Degradation Low and High only; adjustments are not displayed as biological groups",
    "bulk-RNA sensitivity analysis; not proof of T-cell-intrinsic state or causal CD28 signaling"
  )
)
fwrite(method_contract, file.path(table_dir, "00_method_contract.csv"))

package_versions <- rbind(
  data.table(package = "R", version = R.version.string),
  data.table(package = required_packages, version = vapply(required_packages, function(x) as.character(packageVersion(x)), character(1)))
)
fwrite(package_versions, file.path(table_dir, "00_package_versions.csv"))

check_values <- list(
  nrow(meta) == 73L,
  uniqueN(meta$sample_id) == 73L,
  paste(sum(meta$degradation_group == "Low"), sum(meta$degradation_group == "High"), sep = "/") == "37/36",
  nrow(log_expression) == 17225L,
  nrow(program_mat) == 7L && ncol(program_mat) == 73L,
  all(is.finite(program_mat)),
  max(abs(rowMeans(program_mat))) < 1e-10,
  max(abs(apply(program_mat, 1L, sd) - 1)) < 1e-10,
  nrow(abundance_long) == 8L * 73L,
  !anyNA(abundance_long$abundance_z),
  nrow(adjusted_results) == 8L * 7L,
  all(adjusted_results$design_rank == adjusted_results$design_columns_n),
  all(adjusted_results$p_value >= 0 & adjusted_results$p_value <= 1),
  all(nfkb_check$absolute_difference < 1e-10),
  !any(grepl("Primary clinical model|ESTIMATE immune/stromal", c(
    p_before$labels$title, p_before$labels$subtitle,
    p_after$labels$title, p_after$labels$subtitle,
    p_abundance$labels$title, p_abundance$labels$subtitle,
    p_box$labels$title, p_box$labels$subtitle
  ), ignore.case = TRUE))
)
validation <- data.table(
  check = c(
    "PRE sample count", "unique sample IDs", "group sizes Low/High", "expression-filtered gene count",
    "program matrix dimensions", "program scores finite", "program means are zero", "program SDs are one",
    "eight abundance adjusters complete", "abundance values complete", "8 x 7 adjusted results",
    "all adjusted designs full rank", "all p values valid", "NFKB1/RELA rerun matches prior analysis",
    "reader figures contain no model-group labels"
  ),
  pass = unlist(check_values)
)
validation[, value := c(
  nrow(meta), uniqueN(meta$sample_id), "37/36", nrow(log_expression),
  paste(dim(program_mat), collapse = " x "), all(is.finite(program_mat)),
  max(abs(rowMeans(program_mat))), max(abs(apply(program_mat, 1L, sd) - 1)),
  nrow(abundance_long), sum(is.na(abundance_long$abundance_z)), nrow(adjusted_results),
  all(adjusted_results$design_rank == adjusted_results$design_columns_n),
  paste(signif(range(adjusted_results$p_value), 6), collapse = " to "),
  max(nfkb_check$absolute_difference), check_values[[15]]
)]
fwrite(validation, file.path(table_dir, "11_validation_checks.csv"))
if (!all(validation$pass)) {
  print(validation)
  stop("One or more validation checks failed")
}

summary_lines <- c(
  "TIGER PRE73 T-cell abundance-adjusted functional analysis completed.",
  "Groups: Degradation-Low n=37; Degradation-High n=36.",
  "Programs: T-cell proliferation, effector-associated differentiation, cytotoxicity, two CD28 sets, NFKB1 and RELA activities.",
  "Abundance adjustment: six CD8 methods separately, CD8 consensus, and MCPcounter pan-T.",
  paste0("Programs meeting the strict persistence rule: ", paste(robustness[robustness_class == "Persists after CD8-abundance adjustment", program_label], collapse = "; ")),
  paste0(
    "Program versus MCPcounter pan-T Pearson-r range: ",
    sprintf("%.3f to %.3f", min(program_abundance_correlations[abundance_method == "MCPcounter Pan T", Pearson_r]), max(program_abundance_correlations[abundance_method == "MCPcounter Pan T", Pearson_r]))
  ),
  "Reader-facing figures use only Degradation High/Low as biological groups.",
  "Interpretation remains association-level because outcomes and abundance estimates derive from bulk RNA."
)
writeLines(summary_lines, file.path(log_dir, "analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
