#!/usr/bin/env Rscript

# Seven-method tumour-microenvironment estimation for matched TJ-RCC tumours.
# Biological comparison: GC-MS Sarcosine High (n=50) versus Low (n=50).
# Positive standardized effects mean higher estimated abundance in Sarcosine-High.

options(stringsAsFactors = FALSE, width = 180)
set.seed(260902)

analysis_root <- normalizePath(getwd())
if (!all(file.exists(file.path(analysis_root, c("RNAseq_Data", "Metabolomics_Data", "analysis"))))) {
  stop("Run from the ccRCC matched multi-omics analysis root: ", analysis_root)
}
project_root <- normalizePath(file.path(analysis_root, ".."))
tiger_root <- file.path(
  project_root, "공공_BulkRNAseq_ICI반응성관련_분석모음",
  "Sarcosine-TIGER상_ICI반응성_RNAseq분석", "Melanoma-PRJEB23709"
)
project_lib <- file.path(tiger_root, "R_libs")
reference_dir <- file.path(tiger_root, "reference_data", "IOBR_v2.2.3_data-v1.0")
if (dir.exists(project_lib)) .libPaths(unique(c(normalizePath(project_lib), .libPaths())))

required_packages <- c("data.table", "limma", "ggplot2", "IOBR", "GSVA", "pracma")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  library(IOBR)
})

required_reference_files <- file.path(reference_dir, c(
  "lm22.rda", "immuneCuratedData.rda", "cancer_type_genes.rda", "TRef.rda",
  "mRNA_cell_default.rda", "quantiseq_data.rda", "xCell.data.rda",
  "common_genes.rda", "SI_geneset.rda"
))
for (path in required_reference_files) if (!file.exists(path)) stop("Missing IOBR reference: ", path)
options(IOBR.cache_dir = reference_dir)

expression_path <- file.path(analysis_root, "RNAseq_Data", "bulkRNA_matrix_TPM.csv")
group_path <- file.path(
  analysis_root, "results", "00_input_audit_26.09.02", "tumor_sarcosine_group_map.csv"
)
output_root <- file.path(
  analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_HighLow_Deconvolution_26.09.02"
)
raw_dir <- file.path(output_root, "raw_deconvolution")
table_dir <- file.path(output_root, "tables")
figure_dir <- file.path(output_root, "figures")
log_dir <- file.path(output_root, "logs")
for (d in c(raw_dir, table_dir, figure_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

cibersort_permutations <- as.integer(Sys.getenv("SARCO_CIBERSORT_PERM", unset = "1000"))
if (!is.finite(cibersort_permutations) || cibersort_permutations < 100L) {
  stop("SARCO_CIBERSORT_PERM must be an integer >=100")
}
force_deconvolution <- identical(Sys.getenv("SARCO_FORCE_DECONV", unset = "0"), "1")

expr_dt <- fread(expression_path, check.names = FALSE)
meta <- fread(group_path, check.names = FALSE)
stopifnot(nrow(meta) == 100L, uniqueN(meta$sample_id) == 100L)
stopifnot(all(meta$sarcosine_group == ifelse(
  meta$sarcosine_normalized_intensity > 134.4850986, "High", "Low"
)))
meta[, sarcosine_group := factor(sarcosine_group, levels = c("Low", "High"))]
stopifnot(sum(meta$sarcosine_group == "Low") == 50L, sum(meta$sarcosine_group == "High") == 50L)

gene_symbols <- trimws(expr_dt[[1]])
stopifnot(!anyNA(gene_symbols), all(nzchar(gene_symbols)), !anyDuplicated(gene_symbols))
sample_ids <- meta$sample_id
stopifnot(all(sample_ids %in% names(expr_dt)))
log_tpm <- as.matrix(expr_dt[, ..sample_ids])
storage.mode(log_tpm) <- "double"
rownames(log_tpm) <- gene_symbols
colnames(log_tpm) <- sample_ids
stopifnot(!anyNA(log_tpm), all(is.finite(log_tpm)), all(log_tpm >= 0))

# The supplied official TPM file has an exact log2(x+1)-like numeric pattern.
# Inverse values are used only for methods requiring non-log expression.
linear_tpm <- pmax(2^log_tpm - 1, 0)
stopifnot(max(abs(log2(linear_tpm + 1) - log_tpm)) < 1e-10)

meta[, batch := factor(batch)]
meta[, sex := relevel(factor(sex), ref = "female")]
meta[, age := factor(age, levels = c("40-60", "<40", ">60"))]
meta[, grade := factor(grade, levels = c("1", "2", "3", "4"))]
meta[, sarcosine_z := as.numeric(scale(log2_sarcosine_intensity))]

design_primary <- model.matrix(~ batch + sex + age + grade + sarcosine_group, data = meta)
design_unadjusted <- model.matrix(~ sarcosine_group, data = meta)
design_continuous <- model.matrix(~ batch + sex + age + grade + sarcosine_z, data = meta)
stopifnot(
  qr(design_primary)$rank == ncol(design_primary),
  qr(design_unadjusted)$rank == ncol(design_unadjusted),
  qr(design_continuous)$rank == ncol(design_continuous)
)

run_or_load <- function(method_name, expression) {
  output_csv <- file.path(raw_dir, paste0("deconvolution_", method_name, ".csv"))
  if (file.exists(output_csv) && !force_deconvolution) {
    message("Reusing existing checkpoint: ", output_csv)
    return(fread(output_csv, check.names = FALSE))
  }
  message("Running ", method_name, " deconvolution")
  result <- switch(
    method_name,
    CIBERSORT = deconvo_cibersort(
      expression, arrays = FALSE, perm = cibersort_permutations,
      absolute = FALSE, parallel = FALSE, seed = 260902
    ),
    TIMER = deconvo_timer(expression, indications = rep("kirc", ncol(expression))),
    quanTIseq = deconvo_quantiseq(expression, tumor = TRUE, arrays = FALSE, scale_mrna = TRUE),
    EPIC = deconvo_epic(expression, tumor = TRUE),
    MCPcounter = deconvo_mcpcounter(expression),
    xCell = deconvo_xcell(expression, arrays = FALSE),
    ESTIMATE = deconvo_estimate(expression, platform = "affymetrix"),
    stop("Unknown method: ", method_name)
  )
  if (is.null(result)) stop(method_name, " returned NULL")
  result <- as.data.table(result, check.names = FALSE)
  if (!"ID" %in% names(result)) stop(method_name, " output lacks ID")
  stopifnot(nrow(result) == 100L, uniqueN(result$ID) == 100L, setequal(result$ID, sample_ids))
  result <- result[match(sample_ids, ID)]
  stopifnot(identical(result$ID, sample_ids))
  numeric_columns <- setdiff(names(result), "ID")
  stopifnot(all(vapply(result[, ..numeric_columns], is.numeric, logical(1))))
  stopifnot(!anyNA(result[, ..numeric_columns]), all(is.finite(as.matrix(result[, ..numeric_columns]))))
  fwrite(result, output_csv)
  result
}

deconvolution <- list(
  CIBERSORT = run_or_load("CIBERSORT", linear_tpm),
  TIMER = run_or_load("TIMER", linear_tpm),
  quanTIseq = run_or_load("quanTIseq", linear_tpm),
  EPIC = run_or_load("EPIC", linear_tpm),
  MCPcounter = run_or_load("MCPcounter", log_tpm),
  xCell = run_or_load("xCell", linear_tpm),
  ESTIMATE = run_or_load("ESTIMATE", log_tpm)
)

add_direct <- function(method, mapping, measure_type, compartment = "Immune") {
  source <- deconvolution[[method]]
  stopifnot(all(unname(mapping) %in% names(source)))
  rbindlist(lapply(names(mapping), function(cell_type) {
    data.table(
      sample_id = source$ID, method = method, cell_type = cell_type,
      raw_value = source[[mapping[[cell_type]]]], measure_type = measure_type,
      compartment = compartment
    )
  }))
}

add_sum <- function(method, cell_type, columns, measure_type = "fraction", compartment = "Immune") {
  source <- deconvolution[[method]]
  stopifnot(all(columns %in% names(source)))
  data.table(
    sample_id = source$ID, method = method, cell_type = cell_type,
    raw_value = rowSums(as.matrix(source[, ..columns])),
    measure_type = measure_type, compartment = compartment
  )
}

core_long <- rbindlist(list(
  add_sum("CIBERSORT", "B cells", c("B_cells_naive_CIBERSORT", "B_cells_memory_CIBERSORT", "Plasma_cells_CIBERSORT")),
  add_sum("CIBERSORT", "CD8 T cells", "T_cells_CD8_CIBERSORT"),
  add_sum("CIBERSORT", "CD4 T cells", c(
    "T_cells_CD4_naive_CIBERSORT", "T_cells_CD4_memory_resting_CIBERSORT",
    "T_cells_CD4_memory_activated_CIBERSORT", "T_cells_follicular_helper_CIBERSORT"
  )),
  add_sum("CIBERSORT", "Tregs", "T_cells_regulatory_(Tregs)_CIBERSORT"),
  add_sum("CIBERSORT", "NK cells", c("NK_cells_resting_CIBERSORT", "NK_cells_activated_CIBERSORT")),
  add_sum("CIBERSORT", "Monocytes", "Monocytes_CIBERSORT"),
  add_sum("CIBERSORT", "Macrophages", c("Macrophages_M0_CIBERSORT", "Macrophages_M1_CIBERSORT", "Macrophages_M2_CIBERSORT")),
  add_sum("CIBERSORT", "Dendritic cells", c("Dendritic_cells_resting_CIBERSORT", "Dendritic_cells_activated_CIBERSORT")),
  add_sum("CIBERSORT", "Neutrophils", "Neutrophils_CIBERSORT"),

  add_direct("TIMER", c(
    "B cells" = "B_cell_TIMER", "CD4 T cells" = "T_cell_CD4_TIMER",
    "CD8 T cells" = "T_cell_CD8_TIMER", "Neutrophils" = "Neutrophil_TIMER",
    "Macrophages" = "Macrophage_TIMER", "Dendritic cells" = "DC_TIMER"
  ), "abundance score"),

  add_direct("quanTIseq", c(
    "B cells" = "B_cells_quantiseq", "CD4 T cells" = "T_cells_CD4_quantiseq",
    "CD8 T cells" = "T_cells_CD8_quantiseq", "Tregs" = "Tregs_quantiseq",
    "NK cells" = "NK_cells_quantiseq", "Monocytes" = "Monocytes_quantiseq",
    "Macrophages M1" = "Macrophages_M1_quantiseq", "Macrophages M2" = "Macrophages_M2_quantiseq",
    "Dendritic cells" = "Dendritic_cells_quantiseq", "Neutrophils" = "Neutrophils_quantiseq"
  ), "fraction"),
  add_direct("quanTIseq", c("Other/non-reference cells" = "Other_quantiseq"), "fraction", "Tumour/non-reference"),

  add_direct("EPIC", c(
    "B cells" = "Bcells_EPIC", "CD4 T cells" = "CD4_Tcells_EPIC",
    "CD8 T cells" = "CD8_Tcells_EPIC", "NK cells" = "NKcells_EPIC",
    "Macrophages" = "Macrophages_EPIC"
  ), "fraction"),
  add_direct("EPIC", c("CAFs" = "CAFs_EPIC", "Endothelial cells" = "Endothelial_EPIC"), "fraction", "Stromal"),
  add_direct("EPIC", c("Other/non-reference cells" = "otherCells_EPIC"), "fraction", "Tumour/non-reference"),

  add_direct("MCPcounter", c(
    "Pan T cells" = "T_cells_MCPcounter", "CD8 T cells" = "CD8_T_cells_MCPcounter",
    "Cytotoxic lymphocytes" = "Cytotoxic_lymphocytes_MCPcounter",
    "B cells" = "B_lineage_MCPcounter", "NK cells" = "NK_cells_MCPcounter",
    "Monocytes" = "Monocytic_lineage_MCPcounter",
    "Dendritic cells" = "Myeloid_dendritic_cells_MCPcounter",
    "Neutrophils" = "Neutrophils_MCPcounter"
  ), "abundance score"),
  add_direct("MCPcounter", c(
    "Endothelial cells" = "Endothelial_cells_MCPcounter", "Fibroblasts" = "Fibroblasts_MCPcounter"
  ), "abundance score", "Stromal"),

  add_direct("xCell", c(
    "B cells" = "B-cells_xCell", "CD4 T cells" = "CD4+_T-cells_xCell",
    "CD8 T cells" = "CD8+_T-cells_xCell", "Tregs" = "Tregs_xCell",
    "NK cells" = "NK_cells_xCell", "Monocytes" = "Monocytes_xCell",
    "Macrophages" = "Macrophages_xCell", "Dendritic cells" = "DC_xCell",
    "Neutrophils" = "Neutrophils_xCell"
  ), "enrichment score"),
  add_direct("xCell", c(
    "Fibroblasts" = "Fibroblasts_xCell", "Endothelial cells" = "Endothelial_cells_xCell"
  ), "enrichment score", "Stromal"),
  add_direct("xCell", c(
    "Immune score" = "ImmuneScore_xCell", "Stromal score" = "StromaScore_xCell",
    "Microenvironment score" = "MicroenvironmentScore_xCell"
  ), "enrichment score", "Summary"),

  add_direct("ESTIMATE", c(
    "Immune score" = "ImmuneScore_estimate", "Stromal score" = "StromalScore_estimate",
    "Microenvironment score" = "ESTIMATEScore_estimate"
  ), "abundance score", "Summary"),
  add_direct("ESTIMATE", c("Tumour purity" = "TumorPurity_estimate"), "fraction", "Tumour/non-reference")
), use.names = TRUE)
stopifnot(uniqueN(core_long$sample_id) == 100L, !anyNA(core_long$raw_value))

core_long[measure_type == "fraction", transformed_value := asin(sqrt(pmin(pmax(raw_value, 0), 1)))]
core_long[measure_type != "fraction", transformed_value := raw_value]
core_long[, value_sd := sd(transformed_value), by = .(method, cell_type)]
core_long <- core_long[is.finite(value_sd) & value_sd > 0]
core_long[, analysis_value := as.numeric(scale(transformed_value)), by = .(method, cell_type)]
core_long[, feature_id := paste(method, cell_type, sep = "||")]

wide <- dcast(core_long, sample_id ~ feature_id, value.var = "analysis_value")
wide <- wide[match(sample_ids, sample_id)]
stopifnot(identical(wide$sample_id, sample_ids))
score_matrix <- t(as.matrix(wide[, -"sample_id"]))
storage.mode(score_matrix) <- "double"

fit_scores <- function(design, coefficient, model_label) {
  stopifnot(coefficient %in% colnames(design))
  fit <- eBayes(lmFit(score_matrix, design), robust = TRUE)
  tt <- as.data.table(topTable(fit, coef = coefficient, number = Inf, sort.by = "none", adjust.method = "BH"), keep.rownames = "feature_id")
  se <- fit$stdev.unscaled[, coefficient] * sqrt(fit$s2.post)
  crit <- qt(0.975, df = fit$df.total)
  if (length(crit) == 1L) crit <- rep(crit, nrow(tt))
  tt[, standard_error := se[feature_id]]
  tt[, `:=`(CI_low = logFC - crit * standard_error, CI_high = logFC + crit * standard_error)]
  setnames(tt, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"),
           c("standardized_effect", "average_z", "moderated_t", "p_value", "global_BH_q", "B_statistic"))
  tt[, c("method", "cell_type") := tstrsplit(feature_id, "\\|\\|", fixed = FALSE)]
  tt[, model := model_label]
  tt[, BH_q_within_method := p.adjust(p_value, method = "BH"), by = method]
  merge(tt, unique(core_long[, .(method, cell_type, compartment, measure_type)]), by = c("method", "cell_type"), all.x = TRUE)
}

stats_primary <- fit_scores(design_primary, "sarcosine_groupHigh", "Adjusted Sarcosine High vs Low")
stats_unadjusted <- fit_scores(design_unadjusted, "sarcosine_groupHigh", "Unadjusted Sarcosine High vs Low")
stats_continuous <- fit_scores(design_continuous, "sarcosine_z", "Per 1 SD higher log2 Sarcosine")
fwrite(core_long, file.path(table_dir, "core_celltype_values_long.csv"))
fwrite(stats_primary, file.path(table_dir, "core_celltype_adjusted_High_vs_Low.csv"))
fwrite(stats_unadjusted, file.path(table_dir, "core_celltype_unadjusted_High_vs_Low.csv"))
fwrite(stats_continuous, file.path(table_dir, "core_celltype_adjusted_continuous_log2_Sarcosine.csv"))

immune_consensus <- stats_primary[compartment == "Immune", .(
  methods_tested = uniqueN(method),
  method_list = paste(sort(unique(method)), collapse = ";"),
  positive_methods = sum(standardized_effect > 0),
  negative_methods = sum(standardized_effect < 0),
  positive_q_lt_0_05 = sum(standardized_effect > 0 & BH_q_within_method < 0.05),
  negative_q_lt_0_05 = sum(standardized_effect < 0 & BH_q_within_method < 0.05),
  median_standardized_effect = median(standardized_effect),
  min_standardized_effect = min(standardized_effect),
  max_standardized_effect = max(standardized_effect),
  direction_agreement_fraction = max(sum(standardized_effect > 0), sum(standardized_effect < 0)) / .N
), by = cell_type]
immune_consensus[, consensus_direction := fifelse(
  direction_agreement_fraction >= 0.75 & median_standardized_effect > 0, "Sarcosine-High",
  fifelse(direction_agreement_fraction >= 0.75 & median_standardized_effect < 0, "Sarcosine-Low", "Mixed")
)]
setorder(immune_consensus, -negative_q_lt_0_05, -positive_q_lt_0_05, -direction_agreement_fraction)
fwrite(immune_consensus, file.path(table_dir, "immune_cell_cross_method_consensus.csv"))

# Six-method CD8 consensus: median of per-method z scores, then re-standardized.
cd8_long <- core_long[cell_type == "CD8 T cells" & method %in% c(
  "CIBERSORT", "EPIC", "MCPcounter", "quanTIseq", "TIMER", "xCell"
), .(sample_id, method, analysis_value)]
stopifnot(nrow(cd8_long) == 600L)
cd8_wide <- dcast(cd8_long, sample_id ~ method, value.var = "analysis_value")
cd8_wide <- cd8_wide[match(sample_ids, sample_id)]
cd8_wide[, CD8_consensus_z := as.numeric(scale(apply(as.matrix(.SD), 1L, median))), .SDcols = setdiff(names(cd8_wide), "sample_id")]
cd8_wide[, CD8_consensus_no_CIBERSORT_z := as.numeric(scale(apply(
  as.matrix(.SD), 1L, median
))), .SDcols = setdiff(names(cd8_wide), c("sample_id", "CIBERSORT", "CD8_consensus_z"))]
cd8_meta <- merge(meta, cd8_wide[, .(sample_id, CD8_consensus_z)], by = "sample_id", sort = FALSE)
cd8_meta <- merge(cd8_meta, cd8_wide[, .(sample_id, CD8_consensus_no_CIBERSORT_z)], by = "sample_id", sort = FALSE)
cd8_meta <- cd8_meta[match(sample_ids, sample_id)]
cd8_matrix <- rbind(
  `CD8 consensus (six methods)` = cd8_meta$CD8_consensus_z,
  `CD8 consensus (excluding CIBERSORT)` = cd8_meta$CD8_consensus_no_CIBERSORT_z
)
colnames(cd8_matrix) <- sample_ids
rownames(cd8_matrix) <- c("CD8 consensus (six methods)", "CD8 consensus (excluding CIBERSORT)")
fit_cd8 <- eBayes(lmFit(cd8_matrix, design_primary), robust = TRUE)
cd8_tt <- as.data.table(topTable(fit_cd8, coef = "sarcosine_groupHigh", number = Inf, sort.by = "none", adjust.method = "BH"), keep.rownames = "feature")
cd8_se <- fit_cd8$stdev.unscaled[, "sarcosine_groupHigh"] * sqrt(fit_cd8$s2.post)
cd8_crit <- qt(0.975, df = fit_cd8$df.total)
cd8_tt[, `:=`(
  standard_error = as.numeric(cd8_se),
  CI_low = logFC - cd8_crit * as.numeric(cd8_se),
  CI_high = logFC + cd8_crit * as.numeric(cd8_se)
)]
setnames(cd8_tt, c("logFC", "t", "P.Value", "adj.P.Val"), c("standardized_effect", "moderated_t", "p_value", "BH_q"))
fwrite(cd8_meta, file.path(table_dir, "CD8_consensus_per_tumour.csv"))
fwrite(cd8_tt, file.path(table_dir, "CD8_consensus_adjusted_High_vs_Low.csv"))

# Deconvolution QC.
cib <- deconvolution$CIBERSORT
cib_cols <- setdiff(grep("_CIBERSORT$", names(cib), value = TRUE), c("P-value_CIBERSORT", "Correlation_CIBERSORT", "RMSE_CIBERSORT"))
q_cols <- grep("_quantiseq$", names(deconvolution$quanTIseq), value = TRUE)
e_cols <- grep("_EPIC$", names(deconvolution$EPIC), value = TRUE)
qc_summary <- data.table(
  metric = c(
    "CIBERSORT permutations", "CIBERSORT samples P<0.05", "CIBERSORT median correlation",
    "CIBERSORT maximum fraction-sum deviation", "quanTIseq maximum fraction-sum deviation",
    "EPIC maximum fraction-sum deviation", "xCell feature count", "all seven methods completed"
  ),
  value = c(
    cibersort_permutations,
    sum(cib[["P-value_CIBERSORT"]] < 0.05),
    median(cib[["Correlation_CIBERSORT"]]),
    max(abs(rowSums(as.matrix(cib[, ..cib_cols])) - 1)),
    max(abs(rowSums(as.matrix(deconvolution$quanTIseq[, ..q_cols])) - 1)),
    max(abs(rowSums(as.matrix(deconvolution$EPIC[, ..e_cols])) - 1)),
    ncol(deconvolution$xCell) - 1L,
    length(deconvolution) == 7L
  )
)
fwrite(qc_summary, file.path(log_dir, "deconvolution_QC_summary.csv"))

# Figures.
COL_LOW <- "#1B9E8F"
COL_HIGH <- "#C43C3C"
COL_NS <- "#8F8F8F"
theme_publication <- theme_classic(base_size = 14) + theme(
  plot.title = element_text(face = "bold", size = 17),
  plot.subtitle = element_text(size = 11, colour = "#444444"),
  axis.title = element_text(face = "bold"),
  strip.text = element_text(face = "bold"),
  legend.position = "bottom"
)

display_cells <- c(
  "CD8 T cells", "Pan T cells", "Cytotoxic lymphocytes", "CD4 T cells", "Tregs",
  "NK cells", "B cells", "Monocytes", "Macrophages", "Dendritic cells",
  "Neutrophils", "Fibroblasts", "Endothelial cells", "Immune score", "Tumour purity"
)
plot_stats <- stats_primary[cell_type %in% display_cells]
plot_stats[, significance := fifelse(BH_q_within_method < 0.05, "BH q<0.05", "BH q>=0.05")]
plot_stats[, display := paste(method, cell_type, sep = " - ")]
setorder(plot_stats, cell_type, standardized_effect)
fwrite(plot_stats, file.path(table_dir, "Figure_Source_Immune_Deconvolution_Effects.csv"))

p_effects <- ggplot(plot_stats, aes(standardized_effect, reorder(display, standardized_effect))) +
  geom_vline(xintercept = 0, colour = "#777777") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high, colour = significance), width = 0, orientation = "y") +
  geom_point(aes(colour = significance), size = 2.3) +
  scale_colour_manual(values = c("BH q<0.05" = COL_HIGH, "BH q>=0.05" = COL_NS)) +
  labs(
    title = "Immune and stromal estimates by measured tissue Sarcosine",
    subtitle = "Positive effect = higher in Sarcosine-High; each method is tested separately",
    x = "Adjusted standardized High - Low difference (95% CI)", y = NULL, colour = NULL
  ) + theme_publication
ggsave(file.path(figure_dir, "Fig_Immune_Cell_Deconvolution_High_vs_Low.png"), p_effects, width = 13, height = 17, dpi = 400)

cd8_plot <- melt(
  cd8_meta,
  id.vars = c("sample_id", "sarcosine_group"),
  measure.vars = c("CD8_consensus_z", "CD8_consensus_no_CIBERSORT_z"),
  variable.name = "consensus", value.name = "CD8_consensus_z"
)
cd8_plot[, consensus := factor(
  consensus,
  levels = c("CD8_consensus_z", "CD8_consensus_no_CIBERSORT_z"),
  labels = c("Six-method consensus", "Sensitivity excluding CIBERSORT")
)]
p_cd8 <- ggplot(cd8_plot, aes(sarcosine_group, CD8_consensus_z, colour = sarcosine_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "black", fill = "white") +
  geom_jitter(width = 0.14, alpha = 0.65, size = 1.8) +
  facet_wrap(~ consensus) +
  scale_colour_manual(values = c(Low = COL_LOW, High = COL_HIGH)) +
  labs(
    title = "Cross-method CD8 T-cell estimate",
    subtitle = sprintf(
      "Adjusted High-Low effect %.2f SD (95%% CI %.2f to %.2f); P=%.3g",
      cd8_tt[feature == "CD8 consensus (six methods)", standardized_effect],
      cd8_tt[feature == "CD8 consensus (six methods)", CI_low],
      cd8_tt[feature == "CD8 consensus (six methods)", CI_high],
      cd8_tt[feature == "CD8 consensus (six methods)", p_value]
    ),
    x = "Tissue Sarcosine group", y = "Six-method CD8 consensus (z score)", colour = NULL
  ) + theme_publication
ggsave(file.path(figure_dir, "Fig_CD8_Consensus_Sarcosine_High_vs_Low.png"), p_cd8, width = 7.5, height = 6.5, dpi = 400)

summary_values <- core_long[cell_type %in% c("Immune score", "Stromal score", "Tumour purity") & method %in% c("ESTIMATE", "xCell")]
summary_values <- merge(summary_values, meta[, .(sample_id, sarcosine_group)], by = "sample_id", all.x = TRUE)
p_tme <- ggplot(summary_values, aes(sarcosine_group, analysis_value, colour = sarcosine_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "black", fill = "white") +
  geom_jitter(width = 0.14, alpha = 0.6, size = 1.3) +
  facet_wrap(method ~ cell_type, scales = "free_y") +
  scale_colour_manual(values = c(Low = COL_LOW, High = COL_HIGH)) +
  labs(
    title = "Tumour microenvironment summary estimates",
    subtitle = "Reader-facing groups use the fixed tumour median split",
    x = "Tissue Sarcosine group", y = "Within-feature standardized estimate", colour = NULL
  ) + theme_publication
ggsave(file.path(figure_dir, "Fig_TME_Summary_Sarcosine_High_vs_Low.png"), p_tme, width = 12, height = 8, dpi = 400)

input_manifest <- data.table(
  input = c(expression_path, group_path, required_reference_files),
  md5 = unname(tools::md5sum(c(expression_path, group_path, required_reference_files))),
  size_bytes = file.info(c(expression_path, group_path, required_reference_files))$size
)
fwrite(input_manifest, file.path(log_dir, "input_manifest.csv"))

validation <- data.table(
  check = c(
    "100 unique tumour samples", "High n=50", "Low n=50", "inverse log transform exact",
    "three designs full rank", "all seven methods returned", "all method IDs aligned",
    "core values finite", "six CD8 methods complete", "CD8 consensus finite",
    "all three figures exist"
  ),
  passed = c(
    nrow(meta) == 100L && uniqueN(meta$sample_id) == 100L,
    sum(meta$sarcosine_group == "High") == 50L,
    sum(meta$sarcosine_group == "Low") == 50L,
    max(abs(log2(linear_tpm + 1) - log_tpm)) < 1e-10,
    all(c(qr(design_primary)$rank == ncol(design_primary), qr(design_unadjusted)$rank == ncol(design_unadjusted), qr(design_continuous)$rank == ncol(design_continuous))),
    length(deconvolution) == 7L,
    all(vapply(deconvolution, function(x) identical(x$ID, sample_ids), logical(1))),
    !anyNA(core_long$analysis_value) && all(is.finite(core_long$analysis_value)),
    nrow(cd8_long) == 600L,
    !anyNA(cd8_meta[, .(CD8_consensus_z, CD8_consensus_no_CIBERSORT_z)]) && all(is.finite(as.matrix(cd8_meta[, .(CD8_consensus_z, CD8_consensus_no_CIBERSORT_z)]))),
    all(file.exists(file.path(figure_dir, c(
      "Fig_Immune_Cell_Deconvolution_High_vs_Low.png",
      "Fig_CD8_Consensus_Sarcosine_High_vs_Low.png",
      "Fig_TME_Summary_Sarcosine_High_vs_Low.png"
    ))))
  )
)
fwrite(validation, file.path(log_dir, "validation_checks.csv"))
if (!all(validation$passed)) stop("One or more validation checks failed")
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))

cat("Deconvolution complete: ", output_root, "\n", sep = "")
print(qc_summary)
print(immune_consensus)
print(cd8_tt)
