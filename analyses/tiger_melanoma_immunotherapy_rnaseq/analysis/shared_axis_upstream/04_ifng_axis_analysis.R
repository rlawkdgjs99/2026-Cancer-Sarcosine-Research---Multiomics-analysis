#!/usr/bin/env Rscript

# Focused IFN-gamma-axis analysis for TIGER melanoma PRE-only bulk RNA-seq.
#
# The analysis separates:
#   (1) candidate upstream/production context (TCR/CD28, T/NK cytotoxicity,
#       IL-12/IL-18 context, type-II-interferon production), and
#   (2) the downstream IFNGR-JAK-STAT1-IRF1 response program.
#
# This is an observational bulk-RNA prioritization analysis. It cannot prove
# that one upstream pathway caused IFNG production or identify the producing
# cell without orthogonal perturbation and cell-resolved validation.

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

required_packages <- c("data.table", "limma", "fgsea", "msigdbr", "ggplot2", "BiocParallel")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
})

prior_root <- file.path(
  analysis_root, "results", "TIGER_PRE73_Degradation_HighLow_GSEA_GO_BP_26.09.02"
)
integrated_root <- file.path(
  analysis_root, "results", "TIGER_PRE73_ThreeCollection_Immune_TME_Integrated_26.09.02"
)
deconv_root <- file.path(
  analysis_root, "results", "TIGER_PRE73_Degradation_HighLow_Deconvolution_7methods_26.09.02"
)
output_root <- file.path(
  analysis_root, "results", "TIGER_PRE73_IFNG_Axis_Focused_26.09.02"
)
table_dir <- file.path(output_root, "tables")
figure_dir <- file.path(output_root, "figures")
log_dir <- file.path(output_root, "logs")
for (d in c(table_dir, figure_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

used_data_root <- normalizePath(file.path(analysis_root, "..", "사용데이터_모음"))
expression_path <- file.path(
  used_data_root, "Source_Input", "TIGER_PRJEB23709_expression_FPKM_gene_by_sample.csv"
)
metadata_path <- file.path(
  used_data_root, "Analysis_Ready", "TIGER_PRE73_Fig4bc_analysis_data.csv"
)
primary_gene_path <- file.path(prior_root, "tables", "02_limma_adjusted_High_vs_Low_all_genes.csv")
hallmark_primary_path <- file.path(prior_root, "tables", "05_fgsea_Hallmark_adjusted_High_vs_Low.csv")
gobp_primary_path <- file.path(prior_root, "tables", "06_fgsea_GO_BP_adjusted_High_vs_Low.csv")
reactome_primary_path <- file.path(integrated_root, "tables", "01_fgsea_Reactome_adjusted_High_vs_Low.csv")
tme_gene_path <- file.path(integrated_root, "tables", "08_limma_TME_adjusted_High_vs_Low_all_genes.csv")
hallmark_tme_path <- file.path(integrated_root, "tables", "09_fgsea_Hallmark_TME_adjusted.csv")
reactome_tme_path <- file.path(integrated_root, "tables", "10_fgsea_Reactome_TME_adjusted.csv")
gobp_tme_path <- file.path(integrated_root, "tables", "11_fgsea_GO_BP_TME_adjusted.csv")
estimate_path <- file.path(deconv_root, "raw_deconvolution", "deconvolution_ESTIMATE.csv")
deconv_core_path <- file.path(deconv_root, "tables", "05_core_celltype_values_long.csv")

required_inputs <- c(
  expression_path, metadata_path, primary_gene_path, hallmark_primary_path,
  gobp_primary_path, reactome_primary_path, tme_gene_path, hallmark_tme_path,
  reactome_tme_path, gobp_tme_path, estimate_path, deconv_core_path
)
for (path in required_inputs) if (!file.exists(path)) stop("Missing input: ", path)

primary_genes <- fread(primary_gene_path, check.names = FALSE)
tme_genes <- fread(tme_gene_path, check.names = FALSE)
primary_enrichment <- rbindlist(list(
  fread(hallmark_primary_path, check.names = FALSE)[, collection := "Hallmark"],
  fread(reactome_primary_path, check.names = FALSE)[, collection := "Reactome"],
  fread(gobp_primary_path, check.names = FALSE)[, collection := "GO:BP"]
), fill = TRUE, use.names = TRUE)
tme_enrichment <- rbindlist(list(
  fread(hallmark_tme_path, check.names = FALSE)[, collection := "Hallmark"],
  fread(reactome_tme_path, check.names = FALSE)[, collection := "Reactome"],
  fread(gobp_tme_path, check.names = FALSE)[, collection := "GO:BP"]
), fill = TRUE, use.names = TRUE)

rank_from_table <- function(x) {
  stopifnot(all(c("gene_symbol", "moderated_t") %in% names(x)))
  out <- x$moderated_t
  names(out) <- x$gene_symbol
  out <- sort(out, decreasing = TRUE)
  stopifnot(!anyNA(out), all(is.finite(out)), !anyDuplicated(names(out)))
  out
}
rank_primary <- rank_from_table(primary_genes)
rank_tme <- rank_from_table(tme_genes)

# MSigDB 2026.1.Hs memberships used by the preceding validated analyses.
hallmark_membership <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
reactome_membership <- as.data.table(msigdbr(
  species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"
))
gobp_membership <- as.data.table(msigdbr(
  species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"
))
tft_membership <- as.data.table(msigdbr(
  species = "Homo sapiens", collection = "C3", subcollection = "TFT:TFT_LEGACY"
))
stopifnot(
  uniqueN(hallmark_membership$db_version) == 1L,
  identical(unique(hallmark_membership$db_version), unique(reactome_membership$db_version)),
  identical(unique(hallmark_membership$db_version), unique(gobp_membership$db_version)),
  identical(unique(hallmark_membership$db_version), unique(tft_membership$db_version))
)
msigdb_version <- unique(hallmark_membership$db_version)

membership_all <- rbindlist(list(
  hallmark_membership[, .(collection = "Hallmark", pathway = gs_name, gene_symbol, description = gs_description)],
  reactome_membership[, .(collection = "Reactome", pathway = gs_name, gene_symbol, description = gs_description)],
  gobp_membership[, .(collection = "GO:BP", pathway = gs_name, gene_symbol, description = gs_description)]
), use.names = TRUE)

get_genes <- function(collection_value, pathway_value) {
  unique(membership_all[
    collection == collection_value & pathway == pathway_value,
    gene_symbol
  ])
}

# Predeclared mechanistic map. "Production" entries are candidate upstream
# contexts/endpoints; "Response" entries read out IFN-gamma receptor signaling.
pathway_catalog <- data.table(
  axis = c(
    rep("IFN-gamma production/upstream context", 10L),
    rep("IFN-gamma response/downstream", 6L),
    rep("Inflammatory comparators", 2L)
  ),
  label = c(
    "Type II IFN production", "Positive regulation of type II IFN production",
    "T-cell receptor signaling", "CD28-family T-cell activation",
    "CD28 co-stimulation", "NK-cell activation", "Leukocyte cytotoxicity",
    "IL-12-family signaling", "IL-12 production", "IL-18 signaling",
    "Hallmark IFN-gamma response", "Reactome IFN-gamma signaling",
    "Response to type II IFN", "Cellular response to type II IFN",
    "MHC-I peptide loading", "MHC-II antigen presentation",
    "IL6-JAK-STAT3", "TNFA-NFKB"
  ),
  collection = c(
    "GO:BP", "GO:BP", "GO:BP", "Reactome", "Reactome", "GO:BP", "GO:BP",
    "Reactome", "GO:BP", "Reactome",
    "Hallmark", "Reactome", "GO:BP", "GO:BP", "Reactome", "Reactome",
    "Hallmark", "Hallmark"
  ),
  pathway = c(
    "GOBP_TYPE_II_INTERFERON_PRODUCTION",
    "GOBP_POSITIVE_REGULATION_OF_TYPE_II_INTERFERON_PRODUCTION",
    "GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
    "REACTOME_REGULATION_OF_T_CELL_ACTIVATION_BY_CD28_FAMILY",
    "REACTOME_CO_STIMULATION_BY_CD28",
    "GOBP_NATURAL_KILLER_CELL_ACTIVATION",
    "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY",
    "REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING",
    "GOBP_INTERLEUKIN_12_PRODUCTION",
    "REACTOME_INTERLEUKIN_18_SIGNALING",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "REACTOME_INTERFERON_GAMMA_SIGNALING",
    "GOBP_RESPONSE_TO_TYPE_II_INTERFERON",
    "GOBP_CELLULAR_RESPONSE_TO_TYPE_II_INTERFERON",
    "REACTOME_ANTIGEN_PRESENTATION_FOLDING_ASSEMBLY_AND_PEPTIDE_LOADING_OF_CLASS_I_MHC",
    "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB"
  )
)
stopifnot(nrow(pathway_catalog) == 18L)

extract_pathway_model <- function(catalog, enrichment, model_label) {
  x <- merge(
    catalog,
    enrichment[, .(collection, pathway, NES, pval, padj, size, leadingEdge)],
    by = c("collection", "pathway"), all.x = TRUE
  )
  x[, `:=`(
    model = model_label,
    mapped_gene_count = vapply(seq_len(.N), function(i) {
      sum(get_genes(collection[i], pathway[i]) %in% names(rank_primary))
    }, integer(1)),
    formally_tested = !is.na(NES)
  )]
  x[]
}
pathway_evidence <- rbindlist(list(
  extract_pathway_model(pathway_catalog, primary_enrichment, "Clinical-covariate adjusted"),
  extract_pathway_model(pathway_catalog, tme_enrichment, "Clinical + ESTIMATE immune/stromal adjusted")
), use.names = TRUE)
pathway_evidence[, direction := fifelse(
  is.na(NES), "Not formally tested (mapped size below analysis threshold)",
  fifelse(NES > 0, "Degradation-High", "Degradation-Low")
)]
setorder(pathway_evidence, axis, label, model)
fwrite(pathway_evidence, file.path(table_dir, "01_IFNG_axis_pathway_evidence.csv"))

# Focal gene-level evidence. Missing entries were not retained by the predeclared
# expression filter and are reported as missing rather than treated as zero.
focal_gene_map <- data.table(
  gene_symbol = c(
    "IFNG", "IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1", "IRF1",
    "STAT4", "TBX21", "EOMES", "IL12A", "IL12B", "IL12RB1", "IL12RB2",
    "IL18", "IL18R1", "IL18RAP", "CD3D", "CD3E", "CD8A", "NKG7",
    "GNLY", "PRF1", "GZMB", "CXCL9", "CXCL10", "CXCL11", "HLA-A",
    "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "NLRC5", "CIITA",
    "HLA-DRA", "HLA-DRB1", "CD274"
  ),
  role = c(
    "Ligand", "Receptor", "Receptor", "Response signaling", "Response signaling",
    "Response TF", "Response TF", "Production TF", "Production TF", "Production TF",
    "IL-12 ligand", "IL-12 ligand", "IL-12 receptor", "IL-12 receptor",
    "IL-18 ligand", "IL-18 receptor", "IL-18 receptor", "T-cell marker",
    "T-cell marker", "CD8 T-cell marker", "NK/cytotoxic marker", "NK/cytotoxic marker",
    "Cytotoxic effector", "Cytotoxic effector", "IFN-gamma target chemokine",
    "IFN-gamma target chemokine", "IFN-gamma target chemokine", "MHC-I",
    "MHC-I", "MHC-I", "MHC-I", "Antigen processing", "Antigen processing",
    "MHC-I regulator", "MHC-II regulator", "MHC-II", "MHC-II", "Checkpoint ligand"
  )
)
focal_genes <- merge(
  focal_gene_map,
  primary_genes[, .(
    gene_symbol, primary_log2FC = log2FC, primary_moderated_t = moderated_t,
    primary_p = p_value, primary_BH_q = BH_q,
    average_log2_FPKM_plus1
  )],
  by = "gene_symbol", all.x = TRUE
)
focal_genes <- merge(
  focal_genes,
  tme_genes[, .(
    gene_symbol, TME_adjusted_log2FC = log2FC,
    TME_adjusted_moderated_t = moderated_t, TME_adjusted_p = p_value,
    TME_adjusted_BH_q = BH_q
  )],
  by = "gene_symbol", all.x = TRUE
)
focal_genes[, retained_by_expression_filter := !is.na(primary_moderated_t)]
focal_genes[, role_order := match(role, unique(focal_gene_map$role))]
setorder(focal_genes, role_order, gene_symbol)
focal_genes[, role_order := NULL]
fwrite(focal_genes, file.path(table_dir, "02_IFNG_axis_focal_gene_evidence.csv"))

parse_leading_edge <- function(x) {
  if (is.na(x) || !nzchar(x)) character() else unique(strsplit(x, ";", fixed = TRUE)[[1]])
}
leading_terms <- pathway_catalog[
  pathway %in% c(
    "HALLMARK_INTERFERON_GAMMA_RESPONSE", "REACTOME_INTERFERON_GAMMA_SIGNALING",
    "GOBP_RESPONSE_TO_TYPE_II_INTERFERON", "GOBP_TYPE_II_INTERFERON_PRODUCTION",
    "GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
    "REACTOME_REGULATION_OF_T_CELL_ACTIVATION_BY_CD28_FAMILY",
    "GOBP_NATURAL_KILLER_CELL_ACTIVATION", "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY"
  )
]
leading_term_rows <- merge(
  leading_terms, primary_enrichment[, .(collection, pathway, NES, padj, leadingEdge)],
  by = c("collection", "pathway"), all.x = TRUE
)
leading_membership <- rbindlist(lapply(seq_len(nrow(leading_term_rows)), function(i) {
  genes <- parse_leading_edge(leading_term_rows$leadingEdge[i])
  data.table(
    gene_symbol = genes, axis = leading_term_rows$axis[i],
    pathway_label = leading_term_rows$label[i], pathway = leading_term_rows$pathway[i],
    NES = leading_term_rows$NES[i], BH_q = leading_term_rows$padj[i]
  )
}))
leading_summary <- leading_membership[, .(
  leading_edge_pathways_n = uniqueN(pathway),
  response_pathways_n = uniqueN(pathway[axis == "IFN-gamma response/downstream"]),
  upstream_pathways_n = uniqueN(pathway[axis == "IFN-gamma production/upstream context"]),
  pathway_labels = paste(sort(unique(pathway_label)), collapse = ";")
), by = gene_symbol]
leading_summary <- merge(
  leading_summary,
  primary_genes[, .(gene_symbol, primary_log2FC = log2FC, primary_moderated_t = moderated_t, primary_BH_q = BH_q)],
  by = "gene_symbol", all.x = TRUE
)
leading_summary <- merge(
  leading_summary,
  tme_genes[, .(gene_symbol, TME_adjusted_log2FC = log2FC, TME_adjusted_moderated_t = moderated_t, TME_adjusted_BH_q = BH_q)],
  by = "gene_symbol", all.x = TRUE
)
setorder(leading_summary, -leading_edge_pathways_n, -response_pathways_n, -upstream_pathways_n, primary_BH_q)
fwrite(leading_summary, file.path(table_dir, "03_IFNG_axis_leading_edge_gene_overlap.csv"))

# Sample-level module scores: mean of gene-wise z scores across PRE samples.
expr_dt <- fread(expression_path, check.names = FALSE)
meta <- fread(metadata_path, check.names = FALSE)
estimate <- fread(estimate_path, check.names = FALSE)
stopifnot(nrow(meta) == 73L, uniqueN(meta$sample_id) == 73L, all(meta$timepoint == "PRE"))
stopifnot(setequal(estimate$ID, meta$sample_id))
estimate <- estimate[match(meta$sample_id, ID)]
stopifnot(identical(estimate$ID, meta$sample_id))

gene_symbols <- trimws(expr_dt[[1]])
stopifnot(!anyNA(gene_symbols), all(nzchar(gene_symbols)), !anyDuplicated(gene_symbols))
sample_ids <- meta$sample_id
fpkm <- as.matrix(expr_dt[, ..sample_ids])
storage.mode(fpkm) <- "double"
rownames(fpkm) <- gene_symbols
colnames(fpkm) <- sample_ids
stopifnot(!anyNA(fpkm), all(is.finite(fpkm)), all(fpkm >= 0))
log_expression_all <- log2(fpkm + 1)

score_gene_set <- function(genes, score_name, minimum_genes = 5L) {
  genes <- intersect(unique(genes), rownames(log_expression_all))
  genes <- genes[apply(log_expression_all[genes, , drop = FALSE], 1, var) > 0]
  if (length(genes) < minimum_genes) {
    stop(score_name, " has only ", length(genes), " usable genes; minimum is ", minimum_genes)
  }
  gene_z <- t(scale(t(log_expression_all[genes, , drop = FALSE])))
  stopifnot(!anyNA(gene_z), all(is.finite(gene_z)))
  score <- colMeans(gene_z)
  score_z <- as.numeric(scale(score))
  names(score_z) <- names(score)
  list(score = score_z, genes = genes, n = length(genes))
}

response_pathways <- c(
  Hallmark = "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  Reactome = "REACTOME_INTERFERON_GAMMA_SIGNALING",
  `GO:BP` = "GOBP_RESPONSE_TO_TYPE_II_INTERFERON"
)
response_membership <- rbindlist(lapply(names(response_pathways), function(collection) {
  pathway <- response_pathways[[collection]]
  data.table(collection = collection, pathway = pathway, gene_symbol = get_genes(collection, pathway))
}))
response_core_genes <- response_membership[, .N, by = gene_symbol][N >= 2L, gene_symbol]
response_union_genes <- unique(response_membership$gene_symbol)
response_core_genes <- setdiff(response_core_genes, "IFNG")
response_score <- score_gene_set(response_core_genes, "Cross-collection IFN-gamma response core")

upstream_catalog <- pathway_catalog[
  axis == "IFN-gamma production/upstream context",
  .(module_label = label, collection, pathway)
]
upstream_catalog[, full_genes := lapply(seq_len(.N), function(i) get_genes(collection[i], pathway[i]))]
upstream_catalog[, mapped_genes_n := lengths(lapply(full_genes, intersect, y = rownames(log_expression_all)))]
upstream_catalog[, formal_GSEA_eligible := mapped_genes_n >= 15L]

upstream_score_objects <- list()
upstream_score_manifest <- list()
for (i in seq_len(nrow(upstream_catalog))) {
  label <- upstream_catalog$module_label[i]
  genes_full <- upstream_catalog$full_genes[[i]]
  genes_for_ifng <- setdiff(genes_full, "IFNG")
  genes_for_response <- setdiff(genes_full, response_union_genes)
  if (length(intersect(genes_for_ifng, rownames(log_expression_all))) >= 5L) {
    upstream_score_objects[[paste(label, "IFNG", sep = "||")]] <- score_gene_set(
      genes_for_ifng, paste0(label, " predictor for IFNG")
    )
  }
  if (length(intersect(genes_for_response, rownames(log_expression_all))) >= 5L) {
    upstream_score_objects[[paste(label, "Response", sep = "||")]] <- score_gene_set(
      genes_for_response, paste0(label, " non-overlap predictor for IFN-gamma response")
    )
  }
  upstream_score_manifest[[length(upstream_score_manifest) + 1L]] <- data.table(
    module_label = label,
    collection = upstream_catalog$collection[i],
    pathway = upstream_catalog$pathway[i],
    mapped_genes_n = upstream_catalog$mapped_genes_n[i],
    formal_GSEA_eligible = upstream_catalog$formal_GSEA_eligible[i],
    IFNG_predictor_genes_n = length(intersect(genes_for_ifng, rownames(log_expression_all))),
    response_nonoverlap_predictor_genes_n = length(intersect(genes_for_response, rownames(log_expression_all))),
    response_overlap_removed_n = length(intersect(genes_full, response_union_genes))
  )
}
upstream_score_manifest <- rbindlist(upstream_score_manifest)
fwrite(upstream_score_manifest, file.path(table_dir, "04_module_gene_set_manifest.csv"))

meta[, degradation_group := factor(
  ifelse(Degradation_score > median(Degradation_score), "High", "Low"),
  levels = c("Low", "High")
)]
meta[, therapy_short := relevel(factor(therapy_short), ref = "antiPD1")]
meta[, response_group := relevel(factor(response_group), ref = "NR")]
meta[, gender := relevel(factor(gender), ref = "Male")]
meta[, age_z := as.numeric(scale(age))]
meta[, immune_score_z := as.numeric(scale(estimate$ImmuneScore_estimate))]
meta[, stromal_score_z := as.numeric(scale(estimate$StromalScore_estimate))]
meta[, IFNG_expression_z := as.numeric(scale(log_expression_all["IFNG", sample_id]))]
meta[, IFNG_response_core_z := response_score$score[sample_id]]
stopifnot(!anyNA(meta[, .(IFNG_expression_z, IFNG_response_core_z)]))

for (key in names(upstream_score_objects)) {
  score_name <- paste0("module__", make.names(key))
  meta[[score_name]] <- upstream_score_objects[[key]]$score[meta$sample_id]
}

sample_scores <- meta[, .(
  sample_id, patient_name, degradation_group, Degradation_score,
  IFNG_log2_FPKM_plus1 = log_expression_all["IFNG", sample_id],
  IFNG_expression_z, IFNG_response_core_z,
  ESTIMATE_immune_score = estimate$ImmuneScore_estimate,
  ESTIMATE_stromal_score = estimate$StromalScore_estimate
)]
for (key in names(upstream_score_objects)) {
  score_name <- paste0("module__", make.names(key))
  sample_scores[[score_name]] <- meta[[score_name]]
}
fwrite(sample_scores, file.path(table_dir, "05_IFNG_axis_sample_level_scores.csv"))

fit_single_association <- function(outcome, predictor, model_type, module_label, outcome_label, genes_n) {
  base_terms <- c("predictor", "therapy_short", "response_group", "age_z", "gender")
  if (model_type == "Clinical + ESTIMATE immune/stromal adjusted") {
    base_terms <- c(base_terms, "immune_score_z", "stromal_score_z")
  }
  dat <- copy(meta)
  dat[, outcome := outcome]
  dat[, predictor := predictor]
  formula_text <- paste("outcome ~", paste(base_terms, collapse = " + "))
  fit <- lm(as.formula(formula_text), data = dat)
  if (qr(model.matrix(fit))$rank != ncol(model.matrix(fit))) stop("Rank-deficient association model")
  co <- summary(fit)$coefficients["predictor", ]
  critical <- qt(0.975, df = df.residual(fit))
  data.table(
    outcome = outcome_label, module_label = module_label, model = model_type,
    predictor_genes_n = genes_n, n = nobs(fit), standardized_beta = unname(co["Estimate"]),
    standard_error = unname(co["Std. Error"]),
    CI_low = unname(co["Estimate"] - critical * co["Std. Error"]),
    CI_high = unname(co["Estimate"] + critical * co["Std. Error"]),
    t = unname(co["t value"]), p_value = unname(co["Pr(>|t|)"]),
    adjusted_R_squared = summary(fit)$adj.r.squared
  )
}

association_results <- list()
for (i in seq_len(nrow(upstream_catalog))) {
  label <- upstream_catalog$module_label[i]
  for (outcome_key in c("IFNG", "Response")) {
    key <- paste(label, outcome_key, sep = "||")
    if (!key %in% names(upstream_score_objects)) next
    predictor <- upstream_score_objects[[key]]$score[meta$sample_id]
    outcome <- if (outcome_key == "IFNG") meta$IFNG_expression_z else meta$IFNG_response_core_z
    outcome_label <- if (outcome_key == "IFNG") "IFNG expression" else "Cross-collection IFN-gamma response core"
    for (model_type in c("Clinical-covariate adjusted", "Clinical + ESTIMATE immune/stromal adjusted")) {
      association_results[[length(association_results) + 1L]] <- fit_single_association(
        outcome, predictor, model_type, label, outcome_label,
        upstream_score_objects[[key]]$n
      )
    }
  }
}
association_results <- rbindlist(association_results)
association_results[, BH_q := p.adjust(p_value, method = "BH"), by = .(outcome, model)]
association_results[, direction := fifelse(standardized_beta > 0, "Positive", "Negative")]
setorder(association_results, outcome, model, BH_q, -standardized_beta)
fwrite(association_results, file.path(table_dir, "06_upstream_module_associations_with_IFNG_and_response.csv"))

# Relate IFNG outcomes to cell estimates while retaining the Degradation group in
# the model. Both axes derive from bulk RNA and therefore remain non-independent.
deconv_core <- fread(deconv_core_path, check.names = FALSE)
cell_features <- c("CD8 T cells", "NK cells", "Cytotoxic lymphocytes", "Pan T cells")
cell_long <- deconv_core[analysis_feature %in% cell_features]
stopifnot(all(cell_long$sample_id %in% meta$sample_id))

fit_cell_association <- function(x, outcome_name, outcome_values) {
  dat <- merge(
    meta[, .(sample_id, therapy_short, response_group, age_z, gender, degradation_group)],
    x[, .(sample_id, predictor = analysis_value)], by = "sample_id", all = FALSE
  )
  dat[, outcome := outcome_values[match(sample_id, meta$sample_id)]]
  fit <- lm(
    outcome ~ predictor + therapy_short + response_group + age_z + gender + degradation_group,
    data = dat
  )
  co <- summary(fit)$coefficients["predictor", ]
  critical <- qt(0.975, df = df.residual(fit))
  data.table(
    outcome = outcome_name, method = unique(x$method),
    cell_estimate = unique(x$analysis_feature), measure_type = unique(x$measure_type),
    n = nobs(fit), standardized_beta = unname(co["Estimate"]),
    standard_error = unname(co["Std. Error"]),
    CI_low = unname(co["Estimate"] - critical * co["Std. Error"]),
    CI_high = unname(co["Estimate"] + critical * co["Std. Error"]),
    t = unname(co["t value"]), p_value = unname(co["Pr(>|t|)"])
  )
}

cell_associations <- rbindlist(lapply(
  split(cell_long, interaction(cell_long$method, cell_long$analysis_feature, drop = TRUE)),
  function(x) rbindlist(list(
    fit_cell_association(x, "IFNG expression", meta$IFNG_expression_z),
    fit_cell_association(x, "Cross-collection IFN-gamma response core", meta$IFNG_response_core_z)
  ))
))
cell_associations[, BH_q := p.adjust(p_value, method = "BH"), by = outcome]
cell_associations[, estimate_label := paste(method, cell_estimate, sep = " — ")]
setorder(cell_associations, outcome, BH_q, -standardized_beta)
fwrite(cell_associations, file.path(table_dir, "07_cell_estimate_associations_with_IFNG_and_response.csv"))

# Orthogonal transcription-factor-target enrichment. MSigDB C3 legacy motif sets
# are promoter motif collections, not direct TF activity measurements.
tft_pathways <- lapply(split(tft_membership$gene_symbol, tft_membership$gs_name), unique)
run_tft_fgsea <- function(stats, model_label) {
  set.seed(42)
  out <- as.data.table(fgseaMultilevel(
    pathways = tft_pathways, stats = stats, minSize = 15L, maxSize = 500L,
    eps = 0, scoreType = "std", nproc = 1,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  ))
  out[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
  out[, `:=`(
    model = model_label,
    direction = fifelse(NES > 0, "Degradation-High", "Degradation-Low")
  )]
  out[, abs_NES_for_sort := abs(NES)]
  setorder(out, padj, -abs_NES_for_sort, pathway)
  out[, abs_NES_for_sort := NULL]
  out[]
}
tft_primary <- run_tft_fgsea(rank_primary, "Clinical-covariate adjusted")
tft_tme <- run_tft_fgsea(rank_tme, "Clinical + ESTIMATE immune/stromal adjusted")
tft_all <- rbindlist(list(tft_primary, tft_tme), use.names = TRUE)
fwrite(tft_all, file.path(table_dir, "08_C3_TFT_legacy_all_results.csv"))
focal_tft_patterns <- c("^IRF1_", "^STAT1_", "^STAT4_", "^NFAT", "NFAT_Q")
focal_tft <- tft_all[Reduce(`|`, lapply(focal_tft_patterns, grepl, x = pathway))]
focal_tft[, TF_family := fifelse(
  grepl("^IRF1", pathway), "IRF1",
  fifelse(grepl("^STAT1", pathway), "STAT1",
          fifelse(grepl("^STAT4", pathway), "STAT4", "NFAT"))
)]
setorder(focal_tft, TF_family, model, padj)
fwrite(focal_tft, file.path(table_dir, "09_IFNG_focal_TF_motif_enrichment.csv"))

# Figures ------------------------------------------------------------------
pathway_plot <- copy(pathway_evidence)
pathway_plot[, model_short := fifelse(
  model == "Clinical-covariate adjusted", "Primary clinical model", "Plus immune/stromal scores"
)]
pathway_plot[, model_short := factor(
  model_short, levels = c("Primary clinical model", "Plus immune/stromal scores")
)]
pathway_plot[, display_label := factor(label, levels = rev(unique(pathway_catalog$label)))]
pathway_plot[, significance := fifelse(
  is.na(padj), "Not formally tested",
  fifelse(padj < 0.05, "BH q<0.05", "BH q>=0.05")
)]
pathway_plot[, point_size := fifelse(is.na(padj), 1.8, pmin(-log10(pmax(padj, 1e-300)), 50))]
pathway_plot[, x_plot := fifelse(is.na(NES), 0, NES)]

p_pathway <- ggplot(pathway_plot, aes(x = x_plot, y = display_label, color = NES, size = point_size, shape = significance)) +
  geom_vline(xintercept = 0, color = "grey65", linewidth = 0.45) +
  geom_point(alpha = 0.92, na.rm = TRUE) +
  facet_wrap(~ model_short, ncol = 2) +
  scale_color_gradient2(low = "#1B9E8F", mid = "grey75", high = "#C43C3C", midpoint = 0, na.value = "grey45", name = "NES") +
  scale_size_continuous(range = c(2.4, 7.2), name = expression(-log[10](BH~q))) +
  scale_shape_manual(values = c("BH q<0.05" = 16, "BH q>=0.05" = 1, "Not formally tested" = 4)) +
  labs(
    title = "IFN-gamma production context and response-axis pathway evidence",
    subtitle = "Positive NES = Degradation-High; x at zero = mapped set below the original GSEA size threshold",
    x = "Normalized enrichment score (High vs Low)", y = NULL, shape = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"),
    legend.position = "bottom", panel.grid = element_blank()
  )
ggsave(file.path(figure_dir, "Fig_IFNG_Axis_Pathway_Evidence.png"), p_pathway,
       width = 13.2, height = 9.2, dpi = 600, bg = "white")
ggsave(file.path(figure_dir, "Fig_IFNG_Axis_Pathway_Evidence.pdf"), p_pathway,
       width = 13.2, height = 9.2, device = cairo_pdf)

focal_plot <- focal_genes[retained_by_expression_filter == TRUE]
focal_plot <- rbindlist(list(
  focal_plot[, .(gene_symbol, role, model = "Primary clinical model", effect = primary_log2FC, BH_q = primary_BH_q)],
  focal_plot[, .(gene_symbol, role, model = "Plus immune/stromal scores", effect = TME_adjusted_log2FC, BH_q = TME_adjusted_BH_q)]
))
focal_plot[, model := factor(model, levels = c("Primary clinical model", "Plus immune/stromal scores"))]
focal_plot[, gene_symbol := factor(gene_symbol, levels = rev(unique(focal_gene_map$gene_symbol[focal_gene_map$gene_symbol %in% focal_plot$gene_symbol])))]
focal_plot[, significant := BH_q < 0.05]
p_gene <- ggplot(focal_plot, aes(x = effect, y = gene_symbol, color = effect, shape = significant)) +
  geom_vline(xintercept = 0, color = "grey65", linewidth = 0.45) +
  geom_point(size = 3.2, alpha = 0.92) +
  facet_wrap(~ model, ncol = 2) +
  scale_color_gradient2(low = "#1B9E8F", mid = "grey75", high = "#C43C3C", midpoint = 0, name = expression(log[2]*FC)) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "BH q<0.05", `FALSE` = "BH q>=0.05"), name = NULL) +
  labs(
    title = "Focal IFN-gamma-axis genes",
    subtitle = "Bulk expression effects; immune/stromal adjustment is a sensitivity model, not purified-cell evidence",
    x = expression(log[2]*" fold change (High minus Low)"), y = NULL
  ) +
  theme_classic(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"), legend.position = "bottom")
ggsave(file.path(figure_dir, "Fig_IFNG_Axis_Focal_Genes.png"), p_gene,
       width = 10.8, height = 12.8, dpi = 600, bg = "white")
ggsave(file.path(figure_dir, "Fig_IFNG_Axis_Focal_Genes.pdf"), p_gene,
       width = 10.8, height = 12.8, device = cairo_pdf)

assoc_plot <- association_results[
  model == "Clinical + ESTIMATE immune/stromal adjusted" &
    outcome == "Cross-collection IFN-gamma response core"
]
assoc_plot[, module_label := factor(module_label, levels = rev(module_label[order(standardized_beta)]))]
p_assoc <- ggplot(assoc_plot, aes(x = standardized_beta, y = module_label, color = BH_q < 0.05)) +
  geom_vline(xintercept = 0, color = "grey65", linewidth = 0.45) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0.14, linewidth = 0.7) +
  geom_point(size = 3.6) +
  scale_color_manual(values = c(`TRUE` = "#C43C3C", `FALSE` = "grey55"),
                     labels = c(`TRUE` = "BH q<0.05", `FALSE` = "BH q>=0.05"), name = NULL) +
  labs(
    title = "Non-overlapping upstream modules associated with the IFN-gamma response core",
    subtitle = "Clinical + ESTIMATE immune/stromal adjusted; response-set genes removed from each predictor",
    x = "Standardized association coefficient", y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
ggsave(file.path(figure_dir, "Fig_IFNG_Upstream_to_Response_Associations.png"), p_assoc,
       width = 10.8, height = 7.4, dpi = 600, bg = "white")
ggsave(file.path(figure_dir, "Fig_IFNG_Upstream_to_Response_Associations.pdf"), p_assoc,
       width = 10.8, height = 7.4, device = cairo_pdf)

cell_plot <- copy(cell_associations)
cell_plot[, estimate_label := factor(estimate_label, levels = rev(unique(estimate_label[order(cell_estimate, method)])))]
p_cell <- ggplot(cell_plot, aes(x = standardized_beta, y = estimate_label, color = BH_q < 0.05)) +
  geom_vline(xintercept = 0, color = "grey65", linewidth = 0.45) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0.14, linewidth = 0.65) +
  geom_point(size = 3.2) +
  facet_wrap(~ outcome, ncol = 2) +
  scale_color_manual(values = c(`TRUE` = "#C43C3C", `FALSE` = "grey55"),
                     labels = c(`TRUE` = "BH q<0.05", `FALSE` = "BH q>=0.05"), name = NULL) +
  labs(
    title = "Cell-estimate associations with IFNG expression and response",
    subtitle = "Adjusted for clinical covariates and Degradation group; both axes derive from bulk RNA",
    x = "Standardized association coefficient", y = NULL
  ) +
  theme_classic(base_size = 10.5) +
  theme(plot.title = element_text(face = "bold"), strip.text = element_text(face = "bold"), legend.position = "bottom")
ggsave(file.path(figure_dir, "Fig_IFNG_Cell_Estimate_Associations.png"), p_cell,
       width = 13.5, height = 9.2, dpi = 600, bg = "white")
ggsave(file.path(figure_dir, "Fig_IFNG_Cell_Estimate_Associations.pdf"), p_cell,
       width = 13.5, height = 9.2, device = cairo_pdf)

# Candidate tiering is evidence-based and deliberately avoids a numerical causal score.
primary_lookup <- pathway_evidence[model == "Clinical-covariate adjusted"]
tme_lookup <- pathway_evidence[model == "Clinical + ESTIMATE immune/stromal adjusted"]
candidate_summary <- data.table(
  candidate = c(
    "T/NK activation-cytotoxicity -> type II IFN production",
    "IL-18 co-stimulation",
    "IL-12 signaling",
    "IFNGR-JAK1/2-STAT1-IRF1 response",
    "IL6-JAK-STAT3 / TNFA-NFKB inflammatory context"
  ),
  role = c("Candidate upstream context", "Candidate upstream co-stimulus", "Candidate upstream co-stimulus", "Downstream response route", "Context/comparator"),
  evidence_tier = c("High-priority", "Plausible but underpowered", "Incomplete", "High-priority", "Co-enriched; non-specific"),
  key_limit = c(
    "Bulk signatures and cell estimates cannot prove the IFNG-producing cell or direction of causality",
    "IL-18 gene sets mapped to fewer than 15 retained genes and were not formally tested by the predeclared GSEA",
    "Reactome IL-12 family/signaling did not pass BH q<0.05; IL12A/B were not retained by the expression filter",
    "Strong downstream response does not identify which upstream stimulus produced IFNG",
    "These pathways can accompany inflammation without being the direct IFNG trigger"
  )
)
fwrite(candidate_summary, file.path(table_dir, "10_IFNG_candidate_evidence_tiers.csv"))

method_contract <- data.table(
  field = c(
    "cohort", "contrast", "direction", "pathway_sources", "response_score",
    "upstream_score", "association_models", "cell_association_model", "TF_analysis",
    "multiple_testing", "causal_boundary"
  ),
  value = c(
    "TIGER melanoma PRE-only n=73 unique patients",
    "Sarcosine-degradation transcriptional score High (n=36) versus Low (n=37)",
    "Positive log2FC/NES/beta denotes Degradation-High or a positive association; it is not measured sarcosine concentration",
    paste0("MSigDB ", msigdb_version, " Hallmark, Reactome and GO:BP; validated fgsea results from analyses 104/106"),
    "Mean gene-wise z score for genes present in at least 2 of Hallmark IFNG response, Reactome IFNG signaling, and GO response to type II IFN; IFNG excluded",
    "Mean gene-wise z score; IFNG removed for IFNG-expression models; all IFNG-response-union genes removed for response models to reduce circularity",
    "Standardized OLS associations with therapy, response, age and sex; sensitivity additionally includes standardized ESTIMATE immune/stromal scores",
    "Standardized OLS with therapy, response, age, sex and Degradation group",
    "MSigDB C3:TFT legacy promoter-motif GSEA; descriptive TF-target enrichment, not direct TF activity",
    "BH within each outcome/model family or TF model",
    "Observational bulk-RNA results prioritize mechanisms but do not establish sarcosine flux, cellular origin, temporal mediation or causality"
  )
)
fwrite(method_contract, file.path(table_dir, "00_method_contract.csv"))

input_manifest <- data.table(
  input = basename(required_inputs), path = required_inputs,
  bytes = file.info(required_inputs)$size,
  md5 = unname(tools::md5sum(required_inputs))
)
fwrite(input_manifest, file.path(table_dir, "00_input_manifest.csv"))

package_versions <- data.table(
  package = c("R", required_packages),
  version = c(as.character(getRversion()), vapply(required_packages, function(x) as.character(packageVersion(x)), character(1)))
)
fwrite(package_versions, file.path(table_dir, "00_package_versions.csv"))

validation <- data.table(
  check = c(
    "PRE n=73", "High/Low counts 36/37", "IFNG retained", "response core >=10 genes",
    "response score finite", "upstream scores finite", "pathway directions valid",
    "association p/q valid", "cell association p/q valid", "TFT q valid",
    "IL18 GSEA ineligible documented", "all figures created"
  ),
  pass = c(
    nrow(meta) == 73L,
    sum(meta$degradation_group == "High") == 36L && sum(meta$degradation_group == "Low") == 37L,
    "IFNG" %in% primary_genes$gene_symbol,
    response_score$n >= 10L,
    all(is.finite(meta$IFNG_response_core_z)),
    all(vapply(upstream_score_objects, function(x) all(is.finite(x$score)), logical(1))),
    all(is.na(pathway_evidence$NES) | is.finite(pathway_evidence$NES)),
    all(association_results$p_value >= 0 & association_results$p_value <= 1) &&
      all(association_results$BH_q >= 0 & association_results$BH_q <= 1),
    all(cell_associations$p_value >= 0 & cell_associations$p_value <= 1) &&
      all(cell_associations$BH_q >= 0 & cell_associations$BH_q <= 1),
    all(is.na(tft_all$padj) | tft_all$padj >= 0 & tft_all$padj <= 1),
    pathway_evidence[pathway == "REACTOME_INTERLEUKIN_18_SIGNALING" & model == "Clinical-covariate adjusted", formally_tested] == FALSE &&
      pathway_evidence[pathway == "REACTOME_INTERLEUKIN_18_SIGNALING" & model == "Clinical-covariate adjusted", mapped_gene_count] < 15L,
    all(file.exists(file.path(figure_dir, c(
      "Fig_IFNG_Axis_Pathway_Evidence.png", "Fig_IFNG_Axis_Focal_Genes.png",
      "Fig_IFNG_Upstream_to_Response_Associations.png", "Fig_IFNG_Cell_Estimate_Associations.png"
    ))))
  )
)
fwrite(validation, file.path(table_dir, "11_validation_checks.csv"))
if (!all(validation$pass)) stop("Validation failure; inspect 11_validation_checks.csv")

primary_ifng_pathways <- pathway_evidence[
  model == "Clinical-covariate adjusted" &
    pathway %in% c(
      "GOBP_TYPE_II_INTERFERON_PRODUCTION", "GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
      "GOBP_NATURAL_KILLER_CELL_ACTIVATION", "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY",
      "REACTOME_INTERFERON_GAMMA_SIGNALING", "HALLMARK_INTERFERON_GAMMA_RESPONSE"
    ),
  .(label, NES, BH_q = padj)
]
top_upstream_response <- association_results[
  outcome == "Cross-collection IFN-gamma response core" &
    model == "Clinical + ESTIMATE immune/stromal adjusted"
][order(BH_q, -standardized_beta)][1:min(.N, 8L)]
top_cell <- cell_associations[order(BH_q, -standardized_beta)][1:min(.N, 12L)]
focal_gene_report <- focal_genes[gene_symbol %in% c(
  "IFNG", "IL18", "IL18R1", "IL18RAP", "STAT4", "TBX21", "STAT1", "IRF1",
  "CXCL9", "CXCL10", "CXCL11", "NLRC5", "CIITA"
), .(gene_symbol, primary_log2FC, primary_BH_q, TME_adjusted_log2FC, TME_adjusted_BH_q)]

report <- c(
  "# TIGER PRE73 focused IFN-gamma-axis analysis",
  "",
  "## Scope and interpretation boundary",
  "- PRE-only n=73; Degradation-High n=36 and Low n=37.",
  "- Positive effects denote Degradation-High. This is a transcriptional degradation-state comparison, not measured sarcosine concentration.",
  "- The analysis distinguishes candidate IFNG-production context from the downstream IFNGR response.",
  "- All results are observational bulk RNA. They prioritize experiments but do not establish cellular origin, mediation or causality.",
  "",
  "## Key pathway evidence",
  paste(capture.output(print(primary_ifng_pathways)), collapse = "\n"),
  "",
  "## Focal gene evidence",
  paste(capture.output(print(focal_gene_report)), collapse = "\n"),
  "",
  "## Upstream-module associations with the non-overlapping IFN-gamma response core",
  paste(capture.output(print(top_upstream_response[, .(
    module_label, predictor_genes_n, standardized_beta, CI_low, CI_high, p_value, BH_q
  )])), collapse = "\n"),
  "",
  "## Cell-estimate associations",
  paste(capture.output(print(top_cell[, .(
    outcome, method, cell_estimate, standardized_beta, CI_low, CI_high, p_value, BH_q
  )])), collapse = "\n"),
  "",
  "## TF-motif note",
  "C3:TFT legacy results are promoter-motif target enrichment and must not be read as direct phospho-protein or TF-activity measurements.",
  paste(capture.output(print(focal_tft[, .(TF_family, pathway, model, NES, padj)])), collapse = "\n"),
  "",
  "## Experimental resolution",
  "The highest-priority test is a two-arm validation: quantify the source arm (CD8/NK activation, IFNG, STAT4/TBX21) and the response arm (pSTAT1-Y701, IRF1, CXCL9/10/11, HLA-I/II). Perturb IFNG/IFNGR or JAK1/2 and repeat in tumour-only and immune co-culture settings. IL-18 is a plausible co-stimulus requiring targeted testing; the current IL-18 gene sets were below the predeclared GSEA size threshold. IL-12 is not established because Reactome IL-12 signaling did not pass BH correction."
)
writeLines(report, file.path(output_root, "FINAL_IFNG_AXIS_REPORT.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))

cat("Focused IFN-gamma-axis analysis complete\n")
cat("Output:", output_root, "\n")
print(primary_ifng_pathways)
print(top_upstream_response[, .(module_label, standardized_beta, BH_q)])
print(focal_tft[, .(TF_family, pathway, model, NES, padj)])
