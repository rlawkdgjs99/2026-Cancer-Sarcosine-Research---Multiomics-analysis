#!/usr/bin/env Rscript

# TIGER PRJEB23709 PRE-only degradation High vs Low pathway analysis
# Primary unit: 73 unique pretreatment patients.
# Positive enrichment statistics denote higher expression in Degradation-High.

options(stringsAsFactors = FALSE, width = 180)
set.seed(42)

extra_lib <- Sys.getenv("SARCO_GSEA_R_LIB", unset = "")
if (nzchar(extra_lib)) {
  stopifnot(dir.exists(extra_lib))
  .libPaths(unique(c(normalizePath(extra_lib), .libPaths())))
}

required_packages <- c(
  "data.table", "limma", "fgsea", "msigdbr", "ggplot2", "org.Hs.eg.db",
  "AnnotationDbi", "statmod", "BiocParallel"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing required R packages: ", paste(missing_packages, collapse = ", "),
       ". If fgsea/msigdbr are available in a project library, set SARCO_GSEA_R_LIB.")
}

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

analysis_root <- normalizePath(getwd())
if (basename(analysis_root) != "Melanoma-PRJEB23709") {
  stop("Run from the Melanoma-PRJEB23709 analysis root. Current directory: ", analysis_root)
}

used_data_root <- normalizePath(file.path(analysis_root, "..", "사용데이터_모음"))
expression_path <- file.path(
  used_data_root, "Source_Input", "TIGER_PRJEB23709_expression_FPKM_gene_by_sample.csv"
)
pre73_path <- file.path(
  used_data_root, "Analysis_Ready", "TIGER_PRE73_Fig4bc_analysis_data.csv"
)

output_root <- file.path(
  analysis_root, "results", "TIGER_PRE73_Degradation_HighLow_GSEA_GO_BP_26.09.02"
)
table_dir <- file.path(output_root, "tables")
figure_dir <- file.path(output_root, "figures")
log_dir <- file.path(output_root, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(expression_path, pre73_path)) {
  if (!file.exists(path)) stop("Missing input: ", path)
}

expr_dt <- fread(expression_path, check.names = FALSE)
meta <- fread(pre73_path, check.names = FALSE)

stopifnot(nrow(meta) == 73L)
stopifnot(uniqueN(meta$sample_id) == 73L, uniqueN(meta$patient_name) == 73L)
stopifnot(all(c(
  "sample_id", "patient_name", "timepoint", "response_group", "therapy_short",
  "age", "gender", "Degradation_score"
) %in% names(meta)))
stopifnot(all(meta$timepoint == "PRE"))
stopifnot(!anyNA(meta[, .(sample_id, patient_name, response_group, therapy_short, age, gender, Degradation_score)]))

gene_symbols <- trimws(expr_dt[[1]])
stopifnot(!anyNA(gene_symbols), all(nzchar(gene_symbols)), !anyDuplicated(gene_symbols))
sample_columns <- names(expr_dt)[-1]
stopifnot(all(meta$sample_id %in% sample_columns))

sample_ids <- meta$sample_id
fpkm <- as.matrix(expr_dt[, ..sample_ids])
storage.mode(fpkm) <- "double"
rownames(fpkm) <- gene_symbols
colnames(fpkm) <- meta$sample_id
stopifnot(ncol(fpkm) == 73L, !anyNA(fpkm), all(is.finite(fpkm)), all(fpkm >= 0))

# Predeclared low-expression filter: FPKM >= 1 in at least 10% of PRE patients.
minimum_samples <- ceiling(0.10 * ncol(fpkm))
keep_expression <- rowSums(fpkm >= 1) >= minimum_samples
log_expression <- log2(fpkm[keep_expression, , drop = FALSE] + 1)
keep_variance <- apply(log_expression, 1, var) > 0
log_expression <- log_expression[keep_variance, , drop = FALSE]
stopifnot(nrow(log_expression) > 5000L)

degradation_median <- median(meta$Degradation_score)
meta[, degradation_group := factor(
  ifelse(Degradation_score > degradation_median, "High", "Low"),
  levels = c("Low", "High")
)]
stopifnot(sum(meta$degradation_group == "High") == 36L)
stopifnot(sum(meta$degradation_group == "Low") == 37L)

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

adjusted_coef <- which(colnames(design_adjusted) == "degradation_groupHigh")
unadjusted_coef <- which(colnames(design_unadjusted) == "degradation_groupHigh")
continuous_coef <- which(colnames(design_continuous) == "degradation_score_z")
stopifnot(length(adjusted_coef) == 1L, length(unadjusted_coef) == 1L, length(continuous_coef) == 1L)

fit_limma <- function(expression_matrix, design, coefficient) {
  fit <- lmFit(expression_matrix, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  result <- topTable(fit, coef = coefficient, number = Inf, sort.by = "none")
  result <- as.data.table(result, keep.rownames = "gene_symbol")
  setnames(result, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"),
           c("log2FC", "average_log2_FPKM_plus1", "moderated_t", "p_value", "BH_q", "B_statistic"))
  stopifnot(nrow(result) == nrow(expression_matrix), !anyNA(result$moderated_t), all(is.finite(result$moderated_t)))
  list(fit = fit, table = result)
}

adjusted <- fit_limma(log_expression, design_adjusted, adjusted_coef)
unadjusted <- fit_limma(log_expression, design_unadjusted, unadjusted_coef)
continuous <- fit_limma(log_expression, design_continuous, continuous_coef)

# The unadjusted limma coefficient must equal the observed High-Low mean difference.
observed_difference <- rowMeans(log_expression[, meta$degradation_group == "High", drop = FALSE]) -
  rowMeans(log_expression[, meta$degradation_group == "Low", drop = FALSE])
unadjusted_lookup <- setNames(unadjusted$table$log2FC, unadjusted$table$gene_symbol)
stopifnot(max(abs(unadjusted_lookup[names(observed_difference)] - observed_difference)) < 1e-10)

write_model_table <- function(result_dt, filename, estimand) {
  out <- copy(result_dt)
  out[, estimand := estimand]
  setcolorder(out, c("gene_symbol", "estimand", setdiff(names(out), c("gene_symbol", "estimand"))))
  fwrite(out, file.path(table_dir, filename))
}

write_model_table(adjusted$table, "02_limma_adjusted_High_vs_Low_all_genes.csv",
                  "Degradation High minus Low; adjusted for therapy, response, age, and sex")
write_model_table(unadjusted$table, "03_limma_unadjusted_High_vs_Low_all_genes.csv",
                  "Degradation High minus Low; unadjusted")
write_model_table(continuous$table, "04_limma_adjusted_continuous_score_all_genes.csv",
                  "Per 1 SD higher continuous Degradation score; adjusted for therapy, response, age, and sex")

rank_from_table <- function(result_dt) {
  stats <- result_dt$moderated_t
  names(stats) <- result_dt$gene_symbol
  stats <- sort(stats, decreasing = TRUE)
  stopifnot(!anyNA(stats), all(is.finite(stats)), !anyDuplicated(names(stats)))
  stats
}

rank_adjusted <- rank_from_table(adjusted$table)
rank_unadjusted <- rank_from_table(unadjusted$table)
rank_continuous <- rank_from_table(continuous$table)

hallmark_membership <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
gobp_membership <- as.data.table(msigdbr(
  species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"
))
stopifnot(uniqueN(hallmark_membership$db_version) == 1L, uniqueN(gobp_membership$db_version) == 1L)
stopifnot(unique(hallmark_membership$db_version) == unique(gobp_membership$db_version))

hallmark_pathways <- split(hallmark_membership$gene_symbol, hallmark_membership$gs_name)
hallmark_pathways <- lapply(hallmark_pathways, unique)
gobp_pathways <- split(gobp_membership$gene_symbol, gobp_membership$gs_name)
gobp_pathways <- lapply(gobp_pathways, unique)

term_metadata <- rbind(
  unique(hallmark_membership[, .(
    collection = "Hallmark", pathway = gs_name, description = gs_description,
    MSigDB_version = db_version
  )]),
  unique(gobp_membership[, .(
    collection = "GO:BP", pathway = gs_name, description = gs_description,
    MSigDB_version = db_version
  )])
)

merge_metadata <- function(result, collection_label) {
  metadata_subset <- term_metadata[term_metadata$collection == collection_label]
  result <- merge(result, metadata_subset, by = c("pathway", "collection"), all.x = TRUE)
  result[, abs_NES_sort := abs(NES)]
  setorder(result, padj, -abs_NES_sort, pathway)
  result[, abs_NES_sort := NULL]
  result[]
}

run_fgsea <- function(pathways, stats, collection_label, min_size, max_size, model_label) {
  set.seed(42)
  result <- as.data.table(fgseaMultilevel(
    pathways = pathways,
    stats = stats,
    minSize = min_size,
    maxSize = max_size,
    eps = 0,
    scoreType = "std",
    nproc = 1,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  ))
  result[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
  result[, `:=`(
    collection = collection_label,
    model = model_label,
    direction = fifelse(NES > 0, "Higher in Degradation-High", "Higher in Degradation-Low")
  )]
  merge_metadata(result, collection_label)
}

fgsea_hallmark <- run_fgsea(
  hallmark_pathways, rank_adjusted, "Hallmark", 10L, 500L,
  "Adjusted Degradation High vs Low"
)
fgsea_gobp <- run_fgsea(
  gobp_pathways, rank_adjusted, "GO:BP", 15L, 500L,
  "Adjusted Degradation High vs Low"
)

# Sensitivity 1: remove SARDH and PIPOX from the ranked universe to test direct score-gene influence.
rank_adjusted_no_score_genes <- rank_adjusted[!names(rank_adjusted) %in% c("SARDH", "PIPOX")]
fgsea_hallmark_no_score_genes <- run_fgsea(
  hallmark_pathways, rank_adjusted_no_score_genes, "Hallmark", 10L, 500L,
  "Adjusted High vs Low; SARDH and PIPOX excluded"
)
fgsea_gobp_no_score_genes <- run_fgsea(
  gobp_pathways, rank_adjusted_no_score_genes, "GO:BP", 15L, 500L,
  "Adjusted High vs Low; SARDH and PIPOX excluded"
)

run_camera_pr <- function(pathways, stats, collection_label, model_label, min_size, max_size) {
  index <- ids2indices(pathways, names(stats), remove.empty = TRUE)
  index <- index[lengths(index) >= min_size & lengths(index) <= max_size]
  result <- as.data.table(cameraPR(
    statistic = stats,
    index = index,
    use.ranks = FALSE,
    inter.gene.cor = 0.01,
    sort = FALSE,
    directional = TRUE
  ), keep.rownames = "pathway")
  result[, `:=`(collection = collection_label, model = model_label)]
  metadata_subset <- term_metadata[term_metadata$collection == collection_label]
  result <- merge(result, metadata_subset, by = c("pathway", "collection"), all.x = TRUE)
  setorder(result, FDR, PValue, pathway)
  result[]
}

camera_hallmark_adjusted <- run_camera_pr(
  hallmark_pathways, rank_adjusted, "Hallmark", "Adjusted High vs Low", 10L, 500L
)
camera_gobp_adjusted <- run_camera_pr(
  gobp_pathways, rank_adjusted, "GO:BP", "Adjusted High vs Low", 15L, 500L
)
camera_hallmark_unadjusted <- run_camera_pr(
  hallmark_pathways, rank_unadjusted, "Hallmark", "Unadjusted High vs Low", 10L, 500L
)
camera_gobp_unadjusted <- run_camera_pr(
  gobp_pathways, rank_unadjusted, "GO:BP", "Unadjusted High vs Low", 15L, 500L
)
camera_hallmark_continuous <- run_camera_pr(
  hallmark_pathways, rank_continuous, "Hallmark", "Adjusted continuous Degradation score", 10L, 500L
)
camera_gobp_continuous <- run_camera_pr(
  gobp_pathways, rank_continuous, "GO:BP", "Adjusted continuous Degradation score", 15L, 500L
)

immune_pattern <- paste(c(
  "IMMUN", "(^|_)T_CELL($|_)", "(^|_)B_CELL($|_)", "LEUKOCYTE", "LYMPHOCYTE", "INTERFERON",
  "CYTOKINE(S)?($|_|\\b)", "ANTIGEN", "NATURAL_KILLER", "NK_CELL", "MACROPHAGE",
  "MONOCYTE", "NEUTROPHIL", "DENDRITIC_CELL", "INFLAM", "(^|_)COMPLEMENT($|_)",
  "CHEMOKINE", "TOLL_LIKE", "TNF", "INTERLEUKIN", "MHC", "CYTOTOX",
  "PHAGOCYT", "GRANULOCYTE", "MYELOID", "ALLOGRAFT", "(^|_)IL[0-9]+",
  "JAK_STAT", "NF.?KB", "(^|_)TCR($|_)", "(^|_)BCR($|_)"
), collapse = "|")

add_immune_flag <- function(result) {
  result[, immune_related := grepl(
    immune_pattern,
    pathway,
    ignore.case = TRUE
  )]
  result[]
}

fgsea_hallmark <- add_immune_flag(fgsea_hallmark)
fgsea_gobp <- add_immune_flag(fgsea_gobp)
fgsea_hallmark_no_score_genes <- add_immune_flag(fgsea_hallmark_no_score_genes)
fgsea_gobp_no_score_genes <- add_immune_flag(fgsea_gobp_no_score_genes)

for (x in c(
  "camera_hallmark_adjusted", "camera_gobp_adjusted", "camera_hallmark_unadjusted",
  "camera_gobp_unadjusted", "camera_hallmark_continuous", "camera_gobp_continuous"
)) {
  assign(x, add_immune_flag(get(x)))
}

fwrite(fgsea_hallmark, file.path(table_dir, "05_fgsea_Hallmark_adjusted_High_vs_Low.csv"))
fwrite(fgsea_gobp, file.path(table_dir, "06_fgsea_GO_BP_adjusted_High_vs_Low.csv"))
fwrite(fgsea_hallmark[immune_related == TRUE], file.path(table_dir, "07_fgsea_Hallmark_immune_subset.csv"))
fwrite(fgsea_gobp[immune_related == TRUE], file.path(table_dir, "08_fgsea_GO_BP_immune_subset.csv"))
fwrite(fgsea_hallmark_no_score_genes, file.path(table_dir, "09_fgsea_Hallmark_score_genes_excluded.csv"))
fwrite(fgsea_gobp_no_score_genes, file.path(table_dir, "10_fgsea_GO_BP_score_genes_excluded.csv"))
fwrite(camera_hallmark_adjusted, file.path(table_dir, "11_cameraPR_Hallmark_adjusted.csv"))
fwrite(camera_gobp_adjusted, file.path(table_dir, "12_cameraPR_GO_BP_adjusted.csv"))
fwrite(camera_hallmark_unadjusted, file.path(table_dir, "13_cameraPR_Hallmark_unadjusted.csv"))
fwrite(camera_gobp_unadjusted, file.path(table_dir, "14_cameraPR_GO_BP_unadjusted.csv"))
fwrite(camera_hallmark_continuous, file.path(table_dir, "15_cameraPR_Hallmark_continuous_score.csv"))
fwrite(camera_gobp_continuous, file.path(table_dir, "16_cameraPR_GO_BP_continuous_score.csv"))

group_membership <- meta[, .(
  sample_id, patient_name, response_group, therapy_short, age, gender,
  Degradation_score, degradation_group
)]
fwrite(group_membership, file.path(table_dir, "01_PRE73_Degradation_HighLow_membership.csv"))

group_balance <- rbindlist(list(
  meta[, .(variable = "group_size", level = as.character(degradation_group), N = .N), by = degradation_group][, degradation_group := NULL],
  meta[, .N, by = .(degradation_group, response_group)][, .(variable = "response_group", level = paste(degradation_group, response_group, sep = ":"), N)],
  meta[, .N, by = .(degradation_group, therapy_short)][, .(variable = "therapy_short", level = paste(degradation_group, therapy_short, sep = ":"), N)],
  meta[, .N, by = .(degradation_group, gender)][, .(variable = "gender", level = paste(degradation_group, gender, sep = ":"), N)]
), fill = TRUE)
fwrite(group_balance, file.path(table_dir, "00_group_balance_counts.csv"))

method_contract <- data.table(
  field = c(
    "observational_unit", "analysis_population", "group_definition", "group_counts",
    "expression_scale", "expression_filter", "primary_gene_model", "primary_rank_metric",
    "primary_enrichment", "gene_set_database", "gene_set_size_filter", "multiple_testing",
    "positive_direction", "immune_subset_rule", "score_gene_sensitivity", "interpretive_boundary"
  ),
  value = c(
    "One unique pretreatment patient/biopsy",
    "TIGER melanoma PRJEB23709 PRE-only n=73",
    paste0("High if Degradation_score > pooled PRE73 median ", format(degradation_median, digits = 17), "; otherwise Low"),
    "High n=36; Low n=37",
    "log2(FPKM+1)",
    paste0("FPKM >=1 in at least ", minimum_samples, " of 73 samples (>=10%), then non-zero variance"),
    "limma empirical-Bayes model adjusted for therapy, response, standardized age, and sex",
    "Moderated t statistic for High minus Low",
    "fgseaMultilevel preranked GSEA; cameraPR as competitive sensitivity",
    paste0("MSigDB ", unique(hallmark_membership$db_version), " Hallmark and C5 GO:BP via msigdbr ", packageVersion("msigdbr")),
    "Hallmark 10-500; GO:BP 15-500 genes after intersection with ranked universe",
    "BH FDR within each tested collection/model",
    "Positive NES/t statistic = higher in Degradation-High",
    paste0("Predeclared case-insensitive pathway-name regex: ", immune_pattern),
    "Repeat primary GSEA after removing SARDH and PIPOX from ranked universe",
    "Observational bulk-tumour association; not causal, predictive, flux, or cell-type-specific evidence"
  )
)
fwrite(method_contract, file.path(table_dir, "00_method_contract.csv"))

rank_audit <- data.table(
  model = c("adjusted_High_vs_Low", "unadjusted_High_vs_Low", "adjusted_continuous_score"),
  genes = c(length(rank_adjusted), length(rank_unadjusted), length(rank_continuous)),
  exact_tie_values = c(sum(duplicated(rank_adjusted)), sum(duplicated(rank_unadjusted)), sum(duplicated(rank_continuous))),
  finite = c(all(is.finite(rank_adjusted)), all(is.finite(rank_unadjusted)), all(is.finite(rank_continuous)))
)
fwrite(rank_audit, file.path(table_dir, "00_rank_tie_audit.csv"))

package_versions <- data.table(
  package = c("R", required_packages),
  version = c(as.character(getRversion()), vapply(required_packages, function(x) as.character(packageVersion(x)), character(1)))
)
fwrite(package_versions, file.path(table_dir, "00_package_versions.csv"))

input_manifest <- data.table(
  input = c("expression_FPKM", "PRE73_metadata"),
  path = c(expression_path, pre73_path),
  bytes = file.info(c(expression_path, pre73_path))$size,
  md5 = unname(tools::md5sum(c(expression_path, pre73_path)))
)
fwrite(input_manifest, file.path(table_dir, "00_input_manifest.csv"))

clean_label <- function(x, width = 52L) {
  x <- sub("^HALLMARK_", "", x)
  x <- sub("^GOBP_", "", x)
  x <- gsub("_", " ", x)
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

select_plot_terms <- function(result, immune_only = FALSE, n_each_direction = 10L) {
  x <- copy(result)
  if (immune_only) x <- x[immune_related == TRUE]
  x <- x[is.finite(NES) & is.finite(padj) & padj < 0.05]
  positive <- x[NES > 0][order(padj, -NES)][seq_len(min(.N, n_each_direction))]
  negative <- x[NES < 0][order(padj, NES)][seq_len(min(.N, n_each_direction))]
  unique(rbind(positive, negative), by = "pathway")
}

plot_enrichment <- function(plot_data, title, subtitle, filename_stub, wrap_width = 52L) {
  if (!nrow(plot_data)) return(invisible(NULL))
  plot_data <- copy(plot_data)
  plot_data[, label := clean_label(pathway, width = wrap_width)]
  plot_data[, label := factor(label, levels = label[order(NES)])]
  plot_data[, q_for_plot := pmax(padj, 1e-300)]
  p <- ggplot(plot_data, aes(x = NES, y = label, color = direction, size = -log10(q_for_plot))) +
    geom_vline(xintercept = 0, linewidth = 0.45, color = "grey55") +
    geom_point(alpha = 0.9) +
    scale_color_manual(values = c(
      "Higher in Degradation-High" = "#C43C3C",
      "Higher in Degradation-Low" = "#1B9E8F"
    ), labels = c(
      "Higher in Degradation-High" = "Degradation-High",
      "Higher in Degradation-Low" = "Degradation-Low"
    )) +
    scale_size_continuous(range = c(2.5, 7), name = expression(-log[10](BH~q))) +
    labs(x = "Normalized enrichment score (High vs Low)", y = NULL, title = title, subtitle = subtitle, color = NULL) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 10),
      axis.text.y = element_text(size = 9),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.justification = "center",
      legend.text = element_text(size = 9),
      plot.margin = margin(8, 18, 8, 8)
    ) +
    guides(
      size = guide_legend(order = 1, nrow = 1),
      color = guide_legend(order = 2, nrow = 1, byrow = TRUE)
    )
  figure_height <- max(7.5, 3.2 + 0.36 * nrow(plot_data))
  ggsave(file.path(figure_dir, paste0(filename_stub, ".png")), p, width = 11.5, height = figure_height, dpi = 600, bg = "white")
  ggsave(file.path(figure_dir, paste0(filename_stub, ".pdf")), p, width = 11.5, height = figure_height, device = cairo_pdf)
}

plot_enrichment(
  select_plot_terms(fgsea_hallmark, immune_only = FALSE),
  "Hallmark pathways by sarcosine-degradation state",
  "TIGER melanoma PRE-only; covariate-adjusted limma rank; BH q < 0.05; positive NES = Degradation-High",
  "Fig_Hallmark_GSEA_Degradation_High_vs_Low",
  wrap_width = 42L
)
plot_enrichment(
  select_plot_terms(fgsea_gobp, immune_only = TRUE),
  "Immune-related GO Biological Process enrichment",
  "BH q < 0.05 immune terms selected by a predeclared keyword rule after testing all GO:BP sets",
  "Fig_Immune_GO_BP_GSEA_Degradation_High_vs_Low",
  wrap_width = 48L
)
plot_enrichment(
  select_plot_terms(fgsea_gobp, immune_only = FALSE),
  "Top GO Biological Process pathways",
  "TIGER melanoma PRE-only; covariate-adjusted High versus Low comparison; BH q < 0.05",
  "Fig_Top_GO_BP_GSEA_Degradation_High_vs_Low",
  wrap_width = 48L
)

count_summary <- data.table(
  analysis = c(
    "Hallmark fgsea", "GO:BP fgsea", "Hallmark immune subset", "GO:BP immune subset",
    "Hallmark score-gene-excluded fgsea", "GO:BP score-gene-excluded fgsea"
  ),
  tested_sets = c(
    nrow(fgsea_hallmark), nrow(fgsea_gobp), sum(fgsea_hallmark$immune_related), sum(fgsea_gobp$immune_related),
    nrow(fgsea_hallmark_no_score_genes), nrow(fgsea_gobp_no_score_genes)
  ),
  BH_q_lt_0_05 = c(
    sum(fgsea_hallmark$padj < 0.05, na.rm = TRUE), sum(fgsea_gobp$padj < 0.05, na.rm = TRUE),
    sum(fgsea_hallmark$padj < 0.05 & fgsea_hallmark$immune_related, na.rm = TRUE),
    sum(fgsea_gobp$padj < 0.05 & fgsea_gobp$immune_related, na.rm = TRUE),
    sum(fgsea_hallmark_no_score_genes$padj < 0.05, na.rm = TRUE), sum(fgsea_gobp_no_score_genes$padj < 0.05, na.rm = TRUE)
  ),
  BH_q_lt_0_25 = c(
    sum(fgsea_hallmark$padj < 0.25, na.rm = TRUE), sum(fgsea_gobp$padj < 0.25, na.rm = TRUE),
    sum(fgsea_hallmark$padj < 0.25 & fgsea_hallmark$immune_related, na.rm = TRUE),
    sum(fgsea_gobp$padj < 0.25 & fgsea_gobp$immune_related, na.rm = TRUE),
    sum(fgsea_hallmark_no_score_genes$padj < 0.25, na.rm = TRUE), sum(fgsea_gobp_no_score_genes$padj < 0.25, na.rm = TRUE)
  )
)
fwrite(count_summary, file.path(table_dir, "17_enrichment_count_summary.csv"))

top_hallmark <- fgsea_hallmark[order(padj, -abs(NES))][1:min(.N, 10)]
top_immune_go <- fgsea_gobp[immune_related == TRUE][order(padj, -abs(NES))][1:min(.N, 20)]

report_lines <- c(
  "# TIGER PRE73 sarcosine-degradation High vs Low pathway analysis",
  "",
  "## Design",
  paste0("- Independent unit: one pretreatment patient (n=73; High=36, Low=37)."),
  paste0("- Degradation High: score > PRE73 median ", format(degradation_median, digits = 6), "."),
  paste0("- Expression: log2(FPKM+1), filtered to FPKM >=1 in at least ", minimum_samples, " samples; ", nrow(log_expression), " genes retained."),
  "- Primary model: limma High vs Low adjusted for therapy, response, age and sex.",
  paste0("- Gene sets: MSigDB ", unique(hallmark_membership$db_version), " Hallmark and GO:BP; fgseaMultilevel; BH correction within collection."),
  "- Positive NES denotes enrichment in Degradation-High.",
  "",
  "## Multiplicity-controlled result counts",
  paste(capture.output(print(count_summary)), collapse = "\n"),
  "",
  "## Top Hallmark results",
  paste(capture.output(print(top_hallmark[, .(pathway, NES, pval, padj, direction)])), collapse = "\n"),
  "",
  "## Top immune-related GO:BP results",
  paste(capture.output(print(top_immune_go[, .(pathway, NES, pval, padj, direction)])), collapse = "\n"),
  "",
  "## Interpretation boundary",
  "This is an observational bulk-tumour transcriptomic association. It does not establish causality, sarcosine flux, clinical prediction, or a cell-type-specific mechanism. Immune-pathway enrichment may reflect immune-cell abundance, activation state, tumour purity, or correlated biology.",
  "",
  "## Sensitivity analyses",
  "- The primary fgsea analysis was repeated after removing SARDH and PIPOX from the ranked universe.",
  "- cameraPR was run for adjusted High/Low, unadjusted High/Low and adjusted continuous Degradation-score ranks.",
  "- Full tables, not only selected immune terms, are retained under tables/."
)
writeLines(report_lines, file.path(output_root, "FINAL_ANALYSIS_REPORT.md"), useBytes = TRUE)

checks <- data.table(
  check = c(
    "PRE73 rows", "unique samples", "unique patients", "all PRE", "expression sample alignment",
    "no duplicated gene symbols", "no missing expression", "nonnegative FPKM", "High count",
    "Low count", "adjusted design full rank", "unadjusted coefficient equals observed difference",
    "finite adjusted ranks", "Hallmark sets tested", "GO:BP sets tested", "score genes removed in sensitivity",
    "all fgsea q in range", "all camera FDR in range", "output figures created"
  ),
  pass = c(
    nrow(meta) == 73L,
    uniqueN(meta$sample_id) == 73L,
    uniqueN(meta$patient_name) == 73L,
    all(meta$timepoint == "PRE"),
    identical(colnames(fpkm), meta$sample_id),
    !anyDuplicated(gene_symbols),
    !anyNA(fpkm),
    all(fpkm >= 0),
    sum(meta$degradation_group == "High") == 36L,
    sum(meta$degradation_group == "Low") == 37L,
    qr(design_adjusted)$rank == ncol(design_adjusted),
    max(abs(unadjusted_lookup[names(observed_difference)] - observed_difference)) < 1e-10,
    all(is.finite(rank_adjusted)),
    nrow(fgsea_hallmark) == 50L,
    nrow(fgsea_gobp) > 3000L,
    !any(names(rank_adjusted_no_score_genes) %in% c("SARDH", "PIPOX")),
    all(is.na(fgsea_hallmark$padj) | (fgsea_hallmark$padj >= 0 & fgsea_hallmark$padj <= 1)) &&
      all(is.na(fgsea_gobp$padj) | (fgsea_gobp$padj >= 0 & fgsea_gobp$padj <= 1)),
    all(is.na(camera_hallmark_adjusted$FDR) | (camera_hallmark_adjusted$FDR >= 0 & camera_hallmark_adjusted$FDR <= 1)) &&
      all(is.na(camera_gobp_adjusted$FDR) | (camera_gobp_adjusted$FDR >= 0 & camera_gobp_adjusted$FDR <= 1)),
    all(file.exists(file.path(figure_dir, c(
      "Fig_Hallmark_GSEA_Degradation_High_vs_Low.png",
      "Fig_Immune_GO_BP_GSEA_Degradation_High_vs_Low.png",
      "Fig_Top_GO_BP_GSEA_Degradation_High_vs_Low.png"
    ))))
  )
)
fwrite(checks, file.path(table_dir, "18_validation_checks.csv"))
if (!all(checks$pass)) stop("One or more validation checks failed; inspect 18_validation_checks.csv")

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))

cat("Analysis complete\n")
cat("Output:", output_root, "\n")
print(count_summary)
