#!/usr/bin/env Rscript

# TIGER PRJEB23709 PRE-only tumour-microenvironment deconvolution.
# Primary comparison: sarcosine-degradation transcriptional score High vs Low.
# Positive standardized effects denote higher estimates in Degradation-High.

options(stringsAsFactors = FALSE, width = 180)
set.seed(42)

analysis_root <- normalizePath(getwd())
if (basename(analysis_root) != "Melanoma-PRJEB23709") {
  stop("Run from the Melanoma-PRJEB23709 analysis root: ", analysis_root)
}

project_lib <- file.path(analysis_root, "R_libs")
if (dir.exists(project_lib)) .libPaths(unique(c(normalizePath(project_lib), .libPaths())))

required_packages <- c("data.table", "limma", "ggplot2", "IOBR", "GSVA", "pracma")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(ggplot2)
  library(IOBR)
})

used_data_root <- normalizePath(file.path(analysis_root, "..", "사용데이터_모음"))
expression_path <- file.path(used_data_root, "Source_Input", "TIGER_PRJEB23709_expression_FPKM_gene_by_sample.csv")
pre73_path <- file.path(used_data_root, "Analysis_Ready", "TIGER_PRE73_Fig4bc_analysis_data.csv")
reference_dir <- normalizePath(file.path(analysis_root, "reference_data", "IOBR_v2.2.3_data-v1.0"))

required_reference_files <- file.path(reference_dir, c(
  "lm22.rda", "immuneCuratedData.rda", "cancer_type_genes.rda", "TRef.rda",
  "mRNA_cell_default.rda", "quantiseq_data.rda", "xCell.data.rda",
  "common_genes.rda", "SI_geneset.rda"
))
for (path in c(expression_path, pre73_path, required_reference_files)) {
  if (!file.exists(path)) stop("Missing required input/reference: ", path)
}
options(IOBR.cache_dir = reference_dir)

output_root <- file.path(
  analysis_root, "results", "TIGER_PRE73_Degradation_HighLow_Deconvolution_7methods_26.09.02"
)
raw_dir <- file.path(output_root, "raw_deconvolution")
table_dir <- file.path(output_root, "tables")
figure_dir <- file.path(output_root, "figures")
log_dir <- file.path(output_root, "logs")
for (d in c(raw_dir, table_dir, figure_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

cibersort_permutations <- as.integer(Sys.getenv("SARCO_CIBERSORT_PERM", unset = "1000"))
if (!is.finite(cibersort_permutations) || cibersort_permutations < 10L) {
  stop("SARCO_CIBERSORT_PERM must be an integer >=10")
}
force_deconvolution <- identical(Sys.getenv("SARCO_FORCE_DECONV", unset = "0"), "1")

expr_dt <- fread(expression_path, check.names = FALSE)
meta <- fread(pre73_path, check.names = FALSE)
stopifnot(nrow(meta) == 73L, uniqueN(meta$sample_id) == 73L, uniqueN(meta$patient_name) == 73L)
stopifnot(all(meta$timepoint == "PRE"))
stopifnot(all(c(
  "sample_id", "patient_name", "timepoint", "response_group", "therapy_short",
  "age", "gender", "Degradation_score"
) %in% names(meta)))
stopifnot(!anyNA(meta[, .(sample_id, patient_name, response_group, therapy_short, age, gender, Degradation_score)]))

gene_symbols <- trimws(expr_dt[[1]])
stopifnot(!anyNA(gene_symbols), all(nzchar(gene_symbols)), !anyDuplicated(gene_symbols))
stopifnot(all(meta$sample_id %in% names(expr_dt)))

sample_ids <- meta$sample_id
fpkm <- as.matrix(expr_dt[, ..sample_ids])
storage.mode(fpkm) <- "double"
rownames(fpkm) <- gene_symbols
colnames(fpkm) <- meta$sample_id
stopifnot(identical(colnames(fpkm), meta$sample_id), !anyNA(fpkm), all(is.finite(fpkm)), all(fpkm >= 0))
log_fpkm <- log2(fpkm + 1)

degradation_median <- median(meta$Degradation_score)
meta[, degradation_group := factor(
  ifelse(Degradation_score > degradation_median, "High", "Low"),
  levels = c("Low", "High")
)]
stopifnot(sum(meta$degradation_group == "High") == 36L, sum(meta$degradation_group == "Low") == 37L)

meta[, therapy_short := relevel(factor(therapy_short), ref = "antiPD1")]
meta[, response_group := relevel(factor(response_group), ref = "NR")]
meta[, gender := relevel(factor(gender), ref = "Male")]
meta[, age_z := as.numeric(scale(age))]
meta[, degradation_score_z := as.numeric(scale(Degradation_score))]

design_adjusted <- model.matrix(
  ~ therapy_short + response_group + age_z + gender + degradation_group,
  data = meta
)
design_unadjusted <- model.matrix(~ degradation_group, data = meta)
design_continuous <- model.matrix(
  ~ therapy_short + response_group + age_z + gender + degradation_score_z,
  data = meta
)
stopifnot(qr(design_adjusted)$rank == ncol(design_adjusted))
stopifnot(qr(design_unadjusted)$rank == ncol(design_unadjusted))
stopifnot(qr(design_continuous)$rank == ncol(design_continuous))
adjusted_coef <- match("degradation_groupHigh", colnames(design_adjusted))
unadjusted_coef <- match("degradation_groupHigh", colnames(design_unadjusted))
continuous_coef <- match("degradation_score_z", colnames(design_continuous))
stopifnot(all(is.finite(c(adjusted_coef, unadjusted_coef, continuous_coef))))

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
      absolute = FALSE, parallel = FALSE, seed = 42
    ),
    TIMER = deconvo_timer(expression, indications = rep("skcm", ncol(expression))),
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
  stopifnot(nrow(result) == 73L, uniqueN(result$ID) == 73L, setequal(result$ID, meta$sample_id))
  result <- result[match(meta$sample_id, ID)]
  stopifnot(identical(result$ID, meta$sample_id))
  numeric_columns <- setdiff(names(result), "ID")
  stopifnot(all(vapply(result[, ..numeric_columns], is.numeric, logical(1))))
  stopifnot(!anyNA(result[, ..numeric_columns]), all(is.finite(as.matrix(result[, ..numeric_columns]))))
  fwrite(result, output_csv)
  result
}

deconvolution <- list(
  CIBERSORT = run_or_load("CIBERSORT", fpkm),
  TIMER = run_or_load("TIMER", fpkm),
  quanTIseq = run_or_load("quanTIseq", fpkm),
  EPIC = run_or_load("EPIC", fpkm),
  MCPcounter = run_or_load("MCPcounter", log_fpkm),
  xCell = run_or_load("xCell", fpkm),
  ESTIMATE = run_or_load("ESTIMATE", log_fpkm)
)

fraction_features <- list(
  CIBERSORT = grep("_CIBERSORT$", names(deconvolution$CIBERSORT), value = TRUE),
  quanTIseq = grep("_quantiseq$", names(deconvolution$quanTIseq), value = TRUE),
  EPIC = grep("_EPIC$", names(deconvolution$EPIC), value = TRUE),
  ESTIMATE = "TumorPurity_estimate"
)
fraction_features$CIBERSORT <- setdiff(
  fraction_features$CIBERSORT,
  c("P-value_CIBERSORT", "Correlation_CIBERSORT", "RMSE_CIBERSORT")
)

method_measure <- c(
  CIBERSORT = "relative fraction", TIMER = "abundance score",
  quanTIseq = "relative fraction", EPIC = "relative fraction",
  MCPcounter = "abundance score", xCell = "enrichment score",
  ESTIMATE = "score/purity estimate"
)

deconv_long <- rbindlist(lapply(names(deconvolution), function(method) {
  x <- copy(deconvolution[[method]])
  features <- setdiff(names(x), "ID")
  if (method == "CIBERSORT") {
    features <- setdiff(features, c("P-value_CIBERSORT", "Correlation_CIBERSORT", "RMSE_CIBERSORT"))
  }
  out <- melt(x[, c("ID", features), with = FALSE], id.vars = "ID", variable.name = "feature", value.name = "raw_value")
  setnames(out, "ID", "sample_id")
  out[, `:=`(
    method = method,
    measure_type = fifelse(feature %in% fraction_features[[method]], "fraction", method_measure[[method]])
  )]
  out
}), use.names = TRUE)
stopifnot(uniqueN(deconv_long$sample_id) == 73L, !anyNA(deconv_long$raw_value))

fit_score_long <- function(score_long, feature_column = "feature") {
  x <- copy(score_long)
  setnames(x, feature_column, "analysis_feature")
  x[measure_type == "fraction", transformed_value := asin(sqrt(pmin(pmax(raw_value, 0), 1)))]
  x[measure_type != "fraction", transformed_value := raw_value]
  x[, transformed_sd := sd(transformed_value), by = .(method, analysis_feature)]
  x <- x[is.finite(transformed_sd) & transformed_sd > 0]
  x[, analysis_value := as.numeric(scale(transformed_value)), by = .(method, analysis_feature)]
  x[, feature_id := paste(method, analysis_feature, sep = "||")]

  wide <- dcast(x, sample_id ~ feature_id, value.var = "analysis_value")
  wide <- wide[match(meta$sample_id, sample_id)]
  stopifnot(identical(wide$sample_id, meta$sample_id))
  score_matrix <- t(as.matrix(wide[, -"sample_id"]))
  storage.mode(score_matrix) <- "double"

  fit_one <- function(design, coefficient, label) {
    fit <- eBayes(lmFit(score_matrix, design), robust = TRUE)
    beta <- fit$coefficients[, coefficient]
    standard_error <- fit$stdev.unscaled[, coefficient] * sqrt(fit$s2.post)
    critical <- qt(0.975, df = fit$df.total)
    data.table(
      feature_id = rownames(fit$coefficients), model = label,
      standardized_effect = beta,
      standard_error = standard_error,
      CI_low = beta - critical * standard_error,
      CI_high = beta + critical * standard_error,
      moderated_t = fit$t[, coefficient],
      p_value = fit$p.value[, coefficient]
    )
  }

  adjusted <- fit_one(design_adjusted, adjusted_coef, "Adjusted High vs Low")
  unadjusted <- fit_one(design_unadjusted, unadjusted_coef, "Unadjusted High vs Low")
  continuous <- fit_one(design_continuous, continuous_coef, "Adjusted continuous Degradation score")

  lookup <- unique(x[, .(feature_id, method, analysis_feature, measure_type)])
  adjusted <- merge(adjusted, lookup, by = "feature_id", all.x = TRUE)
  unadjusted <- merge(unadjusted, lookup, by = "feature_id", all.x = TRUE)
  continuous <- merge(continuous, lookup, by = "feature_id", all.x = TRUE)

  raw_summary <- x[, .(
    High_n = sum(meta$degradation_group[match(sample_id, meta$sample_id)] == "High"),
    Low_n = sum(meta$degradation_group[match(sample_id, meta$sample_id)] == "Low"),
    High_mean_raw = mean(raw_value[meta$degradation_group[match(sample_id, meta$sample_id)] == "High"]),
    Low_mean_raw = mean(raw_value[meta$degradation_group[match(sample_id, meta$sample_id)] == "Low"]),
    High_median_raw = median(raw_value[meta$degradation_group[match(sample_id, meta$sample_id)] == "High"]),
    Low_median_raw = median(raw_value[meta$degradation_group[match(sample_id, meta$sample_id)] == "Low"]),
    raw_mean_difference = mean(raw_value[meta$degradation_group[match(sample_id, meta$sample_id)] == "High"]) -
      mean(raw_value[meta$degradation_group[match(sample_id, meta$sample_id)] == "Low"]),
    wilcoxon_p = suppressWarnings(wilcox.test(
      raw_value[meta$degradation_group[match(sample_id, meta$sample_id)] == "High"],
      raw_value[meta$degradation_group[match(sample_id, meta$sample_id)] == "Low"],
      exact = FALSE
    )$p.value)
  ), by = .(method, analysis_feature, measure_type)]

  adjusted <- merge(adjusted, raw_summary, by = c("method", "analysis_feature", "measure_type"), all.x = TRUE)
  adjusted[, BH_q_within_method := p.adjust(p_value, method = "BH"), by = method]
  adjusted[, BH_q_global := p.adjust(p_value, method = "BH")]
  adjusted[, wilcoxon_BH_q_within_method := p.adjust(wilcoxon_p, method = "BH"), by = method]
  adjusted[, direction := fifelse(standardized_effect > 0, "Higher in Degradation-High", "Higher in Degradation-Low")]
  setorder(adjusted, BH_q_global, p_value)

  unadjusted[, BH_q_within_method := p.adjust(p_value, method = "BH"), by = method]
  continuous[, BH_q_within_method := p.adjust(p_value, method = "BH"), by = method]
  setorder(unadjusted, BH_q_within_method, p_value)
  setorder(continuous, BH_q_within_method, p_value)

  list(values = x, adjusted = adjusted, unadjusted = unadjusted, continuous = continuous)
}

all_feature_stats <- fit_score_long(deconv_long, "feature")
fwrite(all_feature_stats$values, file.path(table_dir, "01_all_deconvolution_values_long.csv"))
fwrite(all_feature_stats$adjusted, file.path(table_dir, "02_all_features_adjusted_High_vs_Low.csv"))
fwrite(all_feature_stats$unadjusted, file.path(table_dir, "03_all_features_unadjusted_High_vs_Low.csv"))
fwrite(all_feature_stats$continuous, file.path(table_dir, "04_all_features_adjusted_continuous_score.csv"))

add_direct_core <- function(method, mapping, measure_type, compartment = "Immune") {
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

sum_core <- function(method, cell_type, columns, measure_type = "fraction", compartment = "Immune") {
  source <- deconvolution[[method]]
  stopifnot(all(columns %in% names(source)))
  data.table(
    sample_id = source$ID, method = method, cell_type = cell_type,
    raw_value = rowSums(as.matrix(source[, ..columns])),
    measure_type = measure_type, compartment = compartment
  )
}

core_long <- rbindlist(list(
  sum_core("CIBERSORT", "B cells", c("B_cells_naive_CIBERSORT", "B_cells_memory_CIBERSORT", "Plasma_cells_CIBERSORT")),
  sum_core("CIBERSORT", "CD8 T cells", "T_cells_CD8_CIBERSORT"),
  sum_core("CIBERSORT", "CD4 T cells", c(
    "T_cells_CD4_naive_CIBERSORT", "T_cells_CD4_memory_resting_CIBERSORT",
    "T_cells_CD4_memory_activated_CIBERSORT", "T_cells_follicular_helper_CIBERSORT"
  )),
  sum_core("CIBERSORT", "Tregs", "T_cells_regulatory_(Tregs)_CIBERSORT"),
  sum_core("CIBERSORT", "NK cells", c("NK_cells_resting_CIBERSORT", "NK_cells_activated_CIBERSORT")),
  sum_core("CIBERSORT", "Monocytes", "Monocytes_CIBERSORT"),
  sum_core("CIBERSORT", "Macrophages", c("Macrophages_M0_CIBERSORT", "Macrophages_M1_CIBERSORT", "Macrophages_M2_CIBERSORT")),
  sum_core("CIBERSORT", "Dendritic cells", c("Dendritic_cells_resting_CIBERSORT", "Dendritic_cells_activated_CIBERSORT")),
  sum_core("CIBERSORT", "Neutrophils", "Neutrophils_CIBERSORT"),

  add_direct_core("TIMER", c(
    "B cells" = "B_cell_TIMER", "CD4 T cells" = "T_cell_CD4_TIMER",
    "CD8 T cells" = "T_cell_CD8_TIMER", "Neutrophils" = "Neutrophil_TIMER",
    "Macrophages" = "Macrophage_TIMER", "Dendritic cells" = "DC_TIMER"
  ), "abundance score"),

  add_direct_core("quanTIseq", c(
    "B cells" = "B_cells_quantiseq", "CD4 T cells" = "T_cells_CD4_quantiseq",
    "CD8 T cells" = "T_cells_CD8_quantiseq", "Tregs" = "Tregs_quantiseq",
    "NK cells" = "NK_cells_quantiseq", "Monocytes" = "Monocytes_quantiseq",
    "Macrophages M1" = "Macrophages_M1_quantiseq", "Macrophages M2" = "Macrophages_M2_quantiseq",
    "Dendritic cells" = "Dendritic_cells_quantiseq", "Neutrophils" = "Neutrophils_quantiseq"
  ), "fraction"),
  add_direct_core("quanTIseq", c("Other/non-reference cells" = "Other_quantiseq"), "fraction", "Tumour/non-reference"),

  add_direct_core("EPIC", c(
    "B cells" = "Bcells_EPIC", "CD4 T cells" = "CD4_Tcells_EPIC",
    "CD8 T cells" = "CD8_Tcells_EPIC", "NK cells" = "NKcells_EPIC",
    "Macrophages" = "Macrophages_EPIC"
  ), "fraction"),
  add_direct_core("EPIC", c("CAFs" = "CAFs_EPIC", "Endothelial cells" = "Endothelial_EPIC"), "fraction", "Stromal"),
  add_direct_core("EPIC", c("Other/non-reference cells" = "otherCells_EPIC"), "fraction", "Tumour/non-reference"),

  add_direct_core("MCPcounter", c(
    "Pan T cells" = "T_cells_MCPcounter", "CD8 T cells" = "CD8_T_cells_MCPcounter",
    "Cytotoxic lymphocytes" = "Cytotoxic_lymphocytes_MCPcounter",
    "B cells" = "B_lineage_MCPcounter", "NK cells" = "NK_cells_MCPcounter",
    "Monocytes" = "Monocytic_lineage_MCPcounter",
    "Dendritic cells" = "Myeloid_dendritic_cells_MCPcounter",
    "Neutrophils" = "Neutrophils_MCPcounter"
  ), "abundance score"),
  add_direct_core("MCPcounter", c(
    "Endothelial cells" = "Endothelial_cells_MCPcounter",
    "Fibroblasts" = "Fibroblasts_MCPcounter"
  ), "abundance score", "Stromal"),

  add_direct_core("xCell", c(
    "B cells" = "B-cells_xCell", "CD4 T cells" = "CD4+_T-cells_xCell",
    "CD8 T cells" = "CD8+_T-cells_xCell", "Tregs" = "Tregs_xCell",
    "NK cells" = "NK_cells_xCell", "Monocytes" = "Monocytes_xCell",
    "Macrophages" = "Macrophages_xCell", "Dendritic cells" = "DC_xCell",
    "Neutrophils" = "Neutrophils_xCell"
  ), "enrichment score"),
  add_direct_core("xCell", c(
    "Fibroblasts" = "Fibroblasts_xCell", "Endothelial cells" = "Endothelial_cells_xCell"
  ), "enrichment score", "Stromal"),
  add_direct_core("xCell", c(
    "Immune score" = "ImmuneScore_xCell", "Stromal score" = "StromaScore_xCell",
    "Microenvironment score" = "MicroenvironmentScore_xCell"
  ), "enrichment score", "Summary"),

  add_direct_core("ESTIMATE", c(
    "Immune score" = "ImmuneScore_estimate", "Stromal score" = "StromalScore_estimate",
    "Microenvironment score" = "ESTIMATEScore_estimate"
  ), "abundance score", "Summary"),
  add_direct_core("ESTIMATE", c("Tumour purity" = "TumorPurity_estimate"), "fraction", "Tumour/non-reference")
), use.names = TRUE)

stopifnot(uniqueN(core_long$sample_id) == 73L, !anyNA(core_long$raw_value))
core_stats <- fit_score_long(core_long, "cell_type")
core_values <- merge(
  core_stats$values,
  unique(core_long[, .(method, cell_type, compartment)]),
  by.x = c("method", "analysis_feature"), by.y = c("method", "cell_type"), all.x = TRUE
)
core_adjusted <- merge(
  core_stats$adjusted,
  unique(core_long[, .(method, cell_type, compartment)]),
  by.x = c("method", "analysis_feature"), by.y = c("method", "cell_type"), all.x = TRUE
)
core_unadjusted <- merge(
  core_stats$unadjusted,
  unique(core_long[, .(method, cell_type, compartment)]),
  by.x = c("method", "analysis_feature"), by.y = c("method", "cell_type"), all.x = TRUE
)
core_continuous <- merge(
  core_stats$continuous,
  unique(core_long[, .(method, cell_type, compartment)]),
  by.x = c("method", "analysis_feature"), by.y = c("method", "cell_type"), all.x = TRUE
)

fwrite(core_values, file.path(table_dir, "05_core_celltype_values_long.csv"))
fwrite(core_adjusted, file.path(table_dir, "06_core_celltype_adjusted_High_vs_Low.csv"))
fwrite(core_unadjusted, file.path(table_dir, "07_core_celltype_unadjusted_High_vs_Low.csv"))
fwrite(core_continuous, file.path(table_dir, "08_core_celltype_adjusted_continuous_score.csv"))

immune_consensus <- core_adjusted[compartment == "Immune", .(
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
), by = .(cell_type = analysis_feature)]
immune_consensus[, consensus_direction := fifelse(
  direction_agreement_fraction >= 0.75 & median_standardized_effect > 0,
  "Higher in Degradation-High",
  fifelse(
    direction_agreement_fraction >= 0.75 & median_standardized_effect < 0,
    "Higher in Degradation-Low", "Mixed"
  )
)]
setorder(immune_consensus, -positive_q_lt_0_05, -direction_agreement_fraction, -median_standardized_effect)
fwrite(immune_consensus, file.path(table_dir, "09_immune_cell_consensus_summary.csv"))

cibersort <- deconvolution$CIBERSORT
cibersort_fraction_columns <- fraction_features$CIBERSORT
quantiseq_fraction_columns <- fraction_features$quanTIseq
epic_fraction_columns <- fraction_features$EPIC
cibersort_fraction_sum <- rowSums(as.matrix(cibersort[, ..cibersort_fraction_columns]))
quantiseq_fraction_sum <- rowSums(as.matrix(deconvolution$quanTIseq[, ..quantiseq_fraction_columns]))
epic_fraction_sum <- rowSums(as.matrix(deconvolution$EPIC[, ..epic_fraction_columns]))

qc_summary <- data.table(
  metric = c(
    "CIBERSORT permutations", "CIBERSORT samples P<0.05", "CIBERSORT median correlation",
    "CIBERSORT max fraction-sum deviation", "quanTIseq max fraction-sum deviation",
    "EPIC max fraction-sum deviation", "xCell features", "all methods completed"
  ),
  value = c(
    as.character(cibersort_permutations),
    as.character(sum(cibersort[["P-value_CIBERSORT"]] < 0.05)),
    format(median(cibersort[["Correlation_CIBERSORT"]]), digits = 8),
    format(max(abs(cibersort_fraction_sum - 1)), scientific = TRUE),
    format(max(abs(quantiseq_fraction_sum - 1)), scientific = TRUE),
    format(max(abs(epic_fraction_sum - 1)), scientific = TRUE),
    as.character(ncol(deconvolution$xCell) - 1L),
    as.character(length(deconvolution) == 7L)
  )
)
fwrite(qc_summary, file.path(table_dir, "10_deconvolution_QC_summary.csv"))

method_contract <- data.table(
  method = names(method_measure),
  primary_quantity = unname(method_measure),
  input_scale = c(
    "linear FPKM; RNA-seq QN disabled", "linear FPKM; SKCM indication",
    "linear FPKM; tumour mode; mRNA scaling", "linear FPKM; tumour reference",
    "log2(FPKM+1)", "linear FPKM; RNA-seq mode", "log2(FPKM+1)"
  ),
  interpretation = c(
    "Relative fractions among LM22 immune cell types; QC P/correlation retained separately",
    "Tumour-type-aware immune abundance scores; not compositional fractions",
    "Relative immune fractions plus Other; sum constrained to one",
    "Fractions for reference immune/stromal types plus non-reference OtherCells",
    "Cell-population abundance scores; not fractions",
    "Cell-type enrichment scores after spillover correction; not fractions",
    "Immune/stromal/ESTIMATE scores and algorithmic tumour-purity estimate"
  )
)
fwrite(method_contract, file.path(table_dir, "00_method_contract.csv"))

reference_manifest <- data.table(
  file = basename(required_reference_files),
  path = required_reference_files,
  bytes = file.info(required_reference_files)$size,
  md5 = unname(tools::md5sum(required_reference_files)),
  sha256 = vapply(required_reference_files, function(path) {
    sha_line <- system2("shasum", c("-a", "256", path), stdout = TRUE)
    sub("[[:space:]].*$", "", sha_line[[1]])
  }, character(1)),
  source = "IOBR/IOBR GitHub release data-v1.0"
)
fwrite(reference_manifest, file.path(table_dir, "00_reference_manifest.csv"))

input_manifest <- data.table(
  input = c("expression_FPKM", "PRE73_metadata"),
  path = c(expression_path, pre73_path),
  bytes = file.info(c(expression_path, pre73_path))$size,
  md5 = unname(tools::md5sum(c(expression_path, pre73_path)))
)
fwrite(input_manifest, file.path(table_dir, "00_input_manifest.csv"))

package_versions <- data.table(
  package = c("R", required_packages),
  version = c(as.character(getRversion()), vapply(required_packages, function(x) as.character(packageVersion(x)), character(1)))
)
fwrite(package_versions, file.path(table_dir, "00_package_versions.csv"))

plot_stats <- core_adjusted[compartment == "Immune"]
plot_stats[, significance := fifelse(BH_q_within_method < 0.001, "***",
  fifelse(BH_q_within_method < 0.01, "**", fifelse(BH_q_within_method < 0.05, "*", ""))
)]
cell_order <- immune_consensus$cell_type
method_order <- c("CIBERSORT", "TIMER", "quanTIseq", "EPIC", "MCPcounter", "xCell")
plot_stats[, analysis_feature := factor(analysis_feature, levels = rev(cell_order))]
plot_stats[, method := factor(method, levels = method_order)]

p_heatmap <- ggplot(plot_stats, aes(x = method, y = analysis_feature, fill = standardized_effect)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = significance), size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "#1B9E8F", mid = "white", high = "#C43C3C", midpoint = 0, name = "Adjusted\nstandardized effect") +
  labs(
    title = "Immune-cell deconvolution across six algorithms",
    subtitle = "Degradation High minus Low; adjusted; blank = not estimated by that method; * within-method BH q<0.05",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "right"
  )
ggsave(file.path(figure_dir, "Fig_Deconvolution_core_immune_heatmap.png"), p_heatmap, width = 9.5, height = 6.5, dpi = 600, bg = "white")
ggsave(file.path(figure_dir, "Fig_Deconvolution_core_immune_heatmap.pdf"), p_heatmap, width = 9.5, height = 6.5, device = cairo_pdf)

tme_plot <- core_adjusted[compartment != "Immune"]
tme_plot[, display_label := paste(analysis_feature, method, sep = " — ")]
setorder(tme_plot, standardized_effect)
tme_plot[, display_label := factor(display_label, levels = display_label)]
p_tme <- ggplot(tme_plot, aes(x = standardized_effect, y = display_label, color = direction)) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.45) +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0, linewidth = 0.6) +
  geom_point(aes(shape = BH_q_within_method < 0.05), size = 3) +
  scale_color_manual(values = c(
    "Higher in Degradation-High" = "#C43C3C",
    "Higher in Degradation-Low" = "#1B9E8F"
  )) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "BH q<0.05", `FALSE` = "BH q>=0.05")) +
  labs(
    title = "Stromal and tumour/non-reference estimates",
    subtitle = "Standardized Degradation High-minus-Low effects; estimates are method-specific",
    x = "Adjusted standardized effect", y = NULL, color = NULL, shape = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
ggsave(file.path(figure_dir, "Fig_Deconvolution_TME_summary_forest.png"), p_tme, width = 10, height = 7.5, dpi = 600, bg = "white")
ggsave(file.path(figure_dir, "Fig_Deconvolution_TME_summary_forest.pdf"), p_tme, width = 10, height = 7.5, device = cairo_pdf)

checks <- data.table(
  check = c(
    "PRE73 rows", "unique samples", "all PRE", "High count", "Low count",
    "seven methods completed", "all raw results aligned", "finite long values",
    "CIBERSORT fractions sum to one", "quanTIseq fractions sum to one",
    "EPIC fractions sum to one", "adjusted design full rank", "all adjusted q in range",
    "all core q in range", "figures created", "reference files present"
  ),
  pass = c(
    nrow(meta) == 73L,
    uniqueN(meta$sample_id) == 73L,
    all(meta$timepoint == "PRE"),
    sum(meta$degradation_group == "High") == 36L,
    sum(meta$degradation_group == "Low") == 37L,
    length(deconvolution) == 7L,
    all(vapply(deconvolution, function(x) identical(x$ID, meta$sample_id), logical(1))),
    all(is.finite(deconv_long$raw_value)),
    max(abs(cibersort_fraction_sum - 1)) < 1e-8,
    max(abs(quantiseq_fraction_sum - 1)) < 1e-8,
    max(abs(epic_fraction_sum - 1)) < 1e-6,
    qr(design_adjusted)$rank == ncol(design_adjusted),
    all(all_feature_stats$adjusted$BH_q_within_method >= 0 & all_feature_stats$adjusted$BH_q_within_method <= 1),
    all(core_adjusted$BH_q_within_method >= 0 & core_adjusted$BH_q_within_method <= 1),
    all(file.exists(file.path(figure_dir, c(
      "Fig_Deconvolution_core_immune_heatmap.png", "Fig_Deconvolution_TME_summary_forest.png"
    )))),
    all(file.exists(required_reference_files))
  )
)
fwrite(checks, file.path(table_dir, "11_validation_checks.csv"))
if (!all(checks$pass)) stop("Validation failure; inspect 11_validation_checks.csv")

report <- c(
  "# TIGER PRE73 seven-method transcriptomic deconvolution",
  "",
  paste0("- PRE patients: n=73; Degradation-High n=36; Low n=37; median cutoff=", format(degradation_median, digits = 8), "."),
  paste0("- CIBERSORT permutations: ", cibersort_permutations, "; RNA-seq quantile normalization disabled."),
  "- Adjusted standardized effects include therapy, response, age and sex.",
  "- Fractions, abundance scores and enrichment scores are retained as distinct method-specific quantities.",
  "- Consensus is descriptive direction agreement; algorithms are not independent studies.",
  "- EPIC/quanTIseq Other and ESTIMATE purity are algorithmic tumour/non-reference estimates, not direct cancer-cell measurements.",
  "",
  "## CIBERSORT QC",
  paste(capture.output(print(qc_summary)), collapse = "\n"),
  "",
  "## Immune-cell consensus",
  paste(capture.output(print(immune_consensus)), collapse = "\n"),
  "",
  "## Interpretation boundary",
  "Bulk-RNA deconvolution cannot prove cell counts, cell-specific transcriptional mechanisms or sarcosine flux. Differences may reflect tumour purity, cell-state expression, reference mismatch or correlated clinical biology."
)
writeLines(report, file.path(output_root, "FINAL_DECONVOLUTION_REPORT.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))

cat("Deconvolution analysis complete\n")
cat("Output:", output_root, "\n")
print(qc_summary)
print(immune_consensus)
