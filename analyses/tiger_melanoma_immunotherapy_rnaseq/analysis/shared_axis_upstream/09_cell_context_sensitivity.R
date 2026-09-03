#!/usr/bin/env Rscript

# TIGER melanoma PRE73: cell-lineage attribution and tumor-microenvironment
# context analysis for the sarcosine-degradation / CD28 / NF-kappaB / IFNG axis.
#
# Integrity boundaries
# --------------------
# - The samples are bulk PRE-treatment tumor RNA-seq. Deconvolution estimates
#   support lineage-context attribution but cannot prove the cell of origin.
# - The exposure is a SARDH/PIPOX transcriptional degradation score, not a
#   measured sarcosine concentration or metabolic flux.
# - Candidate TME programs are association contexts. Cross-sectional models do
#   not establish that a program induces CD28, NF-kappaB, or IFNG signaling.
# - Candidate-context gene sets are made disjoint from the downstream CD28,
#   NFKB1/RELA-regulon, IFNG-response, and IFNG genes before scoring.
# - Reader-facing figures use biological labels; no opaque "+ESTIMATE" panels.

options(stringsAsFactors = FALSE, width = 180)
set.seed(260902)

analysis_root <- normalizePath(getwd())
if (basename(analysis_root) != "Melanoma-PRJEB23709") {
  stop("Run from the Melanoma-PRJEB23709 analysis root: ", analysis_root)
}

gsea_search_root <- file.path(dirname(dirname(dirname(analysis_root))), "2024_Drug_Res_Updates_NSCLC")
gsea_candidates <- list.dirs(gsea_search_root, recursive = TRUE, full.names = TRUE)
gsea_candidates <- gsea_candidates[
  basename(gsea_candidates) == "R_libs" &
    grepl("hallmark_GSEA_Q4_vs_Q1_26.08.26", gsea_candidates, fixed = TRUE)
]
if (length(gsea_candidates) != 1L) stop("Expected exactly one bundled GSEA R library; found ", length(gsea_candidates))
default_gsea_lib <- gsea_candidates[[1]]
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

project_parent <- dirname(analysis_root)
find_unique_file <- function(root, filename) {
  hits <- list.files(root, pattern = paste0("^", filename, "$"), recursive = TRUE, full.names = TRUE)
  if (length(hits) != 1L) stop("Expected exactly one ", filename, "; found ", length(hits))
  normalizePath(hits[[1]])
}
expression_path <- find_unique_file(project_parent, "TIGER_PRJEB23709_expression_FPKM_gene_by_sample.csv")
metadata_path <- find_unique_file(project_parent, "TIGER_PRE73_Fig4bc_analysis_data.csv")
route_root <- file.path(analysis_root, "results", "TIGER_PRE73_CD28_to_IFNG_Route_Narrowing_26.09.02")
deconv_root <- file.path(analysis_root, "results", "TIGER_PRE73_Degradation_HighLow_Deconvolution_7methods_26.09.02")
tf_root <- file.path(analysis_root, "results", "TIGER_PRE73_CD28_Upstream_TF_26.09.02")

route_score_path <- file.path(route_root, "tables", "06_sample_level_route_scores.csv")
score_manifest_path <- file.path(route_root, "tables", "03_score_gene_manifest.csv")
source_membership_path <- file.path(route_root, "tables", "02_MSigDB_Reactome_membership.csv")
deconv_path <- file.path(deconv_root, "tables", "05_core_celltype_values_long.csv")
regulon_path <- file.path(tf_root, "tables", "01_DoRothEA_AB_signed_regulon_all.csv")
required_inputs <- c(expression_path, metadata_path, route_score_path, score_manifest_path,
                     source_membership_path, deconv_path, regulon_path)
for (path in required_inputs) if (!file.exists(path)) stop("Missing input: ", path)

output_root <- file.path(analysis_root, "results", "TIGER_PRE73_CD28_IFNG_CellOrigin_TME_Context_26.09.02")
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

extract_lm_term <- function(fit, term, outcome, predictor, model, family = NA_character_) {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) stop("Term absent from model: ", term)
  estimate <- unname(sm[term, "Estimate"])
  se <- unname(sm[term, "Std. Error"])
  df <- df.residual(fit)
  crit <- qt(0.975, df)
  data.table(
    family = family, outcome = outcome, predictor = predictor, model = model,
    standardized_beta = estimate, standard_error = se,
    CI_low = estimate - crit * se, CI_high = estimate + crit * se,
    t_value = unname(sm[term, "t value"]), residual_df = df,
    p_value = unname(sm[term, "Pr(>|t|)"]),
    partial_R2 = unname(sm[term, "t value"])^2 / (unname(sm[term, "t value"])^2 + df)
  )
}

fit_group_matrix <- function(score_mat, data, extra_terms = character(), model_label) {
  rhs <- c("therapy_short", "response_group", "age_z", "gender", extra_terms, "degradation_group")
  design <- model.matrix(as.formula(paste("~", paste(rhs, collapse = " + "))), data = data)
  if (qr(design)$rank != ncol(design)) stop("Rank-deficient group model: ", model_label)
  fit <- eBayes(lmFit(score_mat[, data$sample_id, drop = FALSE], design), robust = TRUE)
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

vif_one <- function(term, others, data) {
  fit <- lm(as.formula(paste(term, "~", paste(others, collapse = " + "))), data = data)
  1 / (1 - summary(fit)$r.squared)
}

# -------------------------------------------------------------------------
# Inputs and PRE73 expression
# -------------------------------------------------------------------------
expr_dt <- fread(expression_path, check.names = FALSE)
meta <- fread(metadata_path, check.names = FALSE)
route_scores <- fread(route_score_path, check.names = FALSE)
score_manifest <- fread(score_manifest_path, check.names = FALSE)
source_membership <- fread(source_membership_path, check.names = FALSE)
deconv <- fread(deconv_path, check.names = FALSE)
regulon <- fread(regulon_path, check.names = FALSE)

stopifnot(
  nrow(meta) == 73L, uniqueN(meta$sample_id) == 73L,
  nrow(route_scores) == 73L, uniqueN(route_scores$sample_id) == 73L,
  uniqueN(deconv$sample_id) == 73L,
  all(c("sample_id", "therapy_short", "response_group", "age", "gender", "Degradation_score") %in% names(meta)),
  all(c("method", "analysis_feature", "sample_id", "analysis_value") %in% names(deconv))
)

sample_ids <- meta$sample_id
gene_symbols <- trimws(expr_dt[[1]])
if (anyNA(gene_symbols) || any(!nzchar(gene_symbols)) || anyDuplicated(gene_symbols)) stop("Invalid gene-symbol column")
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

analysis_data <- merge(meta, route_scores[, !c("degradation_group", "Degradation_score", "degradation_score_z")], by = "sample_id", all.x = TRUE, sort = FALSE)
analysis_data <- analysis_data[match(sample_ids, sample_id)]
if (anyNA(analysis_data)) stop("Analysis data contain missing values")

# -------------------------------------------------------------------------
# Cross-method deconvolution consensus scores
# -------------------------------------------------------------------------
lineage_map <- data.table(
  analysis_feature = c(
    "B cells", "CD4 T cells", "CD8 T cells", "Dendritic cells",
    "Macrophages", "Macrophages M1", "Macrophages M2", "Monocytes",
    "NK cells", "Neutrophils", "Tregs", "Pan T cells",
    "Cytotoxic lymphocytes", "CAFs", "Fibroblasts"
  ),
  lineage = c(
    "B cells", "CD4 T cells", "CD8 T cells", "Dendritic cells",
    "Macrophages", "Macrophages", "Macrophages", "Monocytes",
    "NK cells", "Neutrophils", "Tregs", "Pan T cells",
    "Cytotoxic lymphocytes", "Stromal cells", "Stromal cells"
  )
)
deconv_lineage <- merge(deconv, lineage_map, by = "analysis_feature", all = FALSE)
deconv_method_lineage <- deconv_lineage[, .(method_lineage_z = mean(analysis_value)), by = .(method, lineage, sample_id)]
deconv_method_lineage[, method_lineage_z := zscore(method_lineage_z), by = .(method, lineage)]
lineage_methods <- deconv_method_lineage[, .(
  contributing_methods_n = uniqueN(method),
  contributing_methods = paste(sort(unique(method)), collapse = ";")
), by = lineage]
lineage_consensus <- deconv_method_lineage[, .(lineage_consensus_raw = mean(method_lineage_z)), by = .(lineage, sample_id)]
lineage_consensus[, lineage_score_z := zscore(lineage_consensus_raw), by = lineage]
stopifnot(all(lineage_consensus[, uniqueN(sample_id), by = lineage]$V1 == 73L))

lineage_wide <- dcast(lineage_consensus, sample_id ~ lineage, value.var = "lineage_score_z")
safe_lineage_col <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
setnames(lineage_wide, setdiff(names(lineage_wide), "sample_id"), safe_lineage_col(setdiff(names(lineage_wide), "sample_id")))
analysis_data <- merge(analysis_data, lineage_wide, by = "sample_id", all.x = TRUE, sort = FALSE)
analysis_data <- analysis_data[match(sample_ids, sample_id)]

analysis_data[, T_cell_context_z := zscore(rowMeans(cbind(CD4_T_cells, CD8_T_cells, Pan_T_cells)))]
analysis_data[, APC_myeloid_context_z := zscore(rowMeans(cbind(B_cells, Dendritic_cells, Macrophages, Monocytes)))]
lineage_composite_manifest <- data.table(
  composite = c("T-cell context", "APC/myeloid-cell context"),
  component_consensus_scores = c("CD4 T cells;CD8 T cells;Pan T cells", "B cells;Dendritic cells;Macrophages;Monocytes"),
  interpretation = c(
    "Mean of standardized cross-method lineage estimates; context/sensitivity variable, not a measured fraction",
    "Mean of standardized cross-method lineage estimates; context/sensitivity variable, not a measured fraction"
  )
)

# -------------------------------------------------------------------------
# Cell-lineage associations with route nodes
# -------------------------------------------------------------------------
node_catalog <- data.table(
  node_order = 1:6,
  node = c("SPI1/PU.1 activity", "CD80/CD86 expression", "CD28 pathway expression",
           "NFKB1/RELA activity", "IFNG expression", "IFNG-response core"),
  column = c("SPI1_lro", "APC_ligand_z", "SPI1_chain_CD28_z", "branch_nfkb", "IFNG_expression_z", "IFNG_response_core_z")
)
lineage_display <- c("Dendritic cells", "Monocytes", "Macrophages", "B cells", "Pan T cells",
                     "CD8 T cells", "CD4 T cells", "NK cells", "Cytotoxic lymphocytes", "Stromal cells")
lineage_catalog <- merge(
  data.table(lineage = lineage_display, lineage_order = seq_along(lineage_display)),
  lineage_methods, by = "lineage", all.x = TRUE, sort = FALSE
)
lineage_catalog[, column := safe_lineage_col(lineage)]
if (anyNA(lineage_catalog)) stop("Missing lineage consensus metadata")

clinical_covars <- c("therapy_short", "response_group", "age_z", "gender")
cell_node_associations <- rbindlist(lapply(seq_len(nrow(node_catalog)), function(i) {
  rbindlist(lapply(seq_len(nrow(lineage_catalog)), function(j) {
    outcome_col <- node_catalog$column[i]
    lineage_col <- lineage_catalog$column[j]
    rhs <- c("degradation_score_z", lineage_col, clinical_covars)
    fit <- lm(as.formula(paste(outcome_col, "~", paste(rhs, collapse = " + "))), data = analysis_data)
    extract_lm_term(fit, lineage_col, node_catalog$node[i], lineage_catalog$lineage[j],
                    "Clinical covariates + degradation score", "Lineage-node association")
  }))
}))
cell_node_associations[, BH_q_within_node := p.adjust(p_value, method = "BH"), by = outcome]
cell_node_associations[, BH_q_global := p.adjust(p_value, method = "BH")]
cell_node_associations <- merge(cell_node_associations, lineage_catalog[, .(predictor = lineage, contributing_methods_n, contributing_methods)], by = "predictor", all.x = TRUE)

# Method-level audit: consensus findings should not be carried by one algorithm.
method_node_associations <- rbindlist(lapply(seq_len(nrow(node_catalog)), function(i) {
  rbindlist(lapply(seq_len(nrow(lineage_catalog)), function(j) {
    pieces <- deconv_method_lineage[lineage == lineage_catalog$lineage[j]]
    rbindlist(lapply(unique(pieces$method), function(this_method) {
      values <- pieces[method == this_method, .(sample_id, method_lineage_z)]
      d <- merge(analysis_data, values, by = "sample_id", all.x = TRUE, sort = FALSE)
      d <- d[match(sample_ids, sample_id)]
      fit <- lm(as.formula(paste(node_catalog$column[i], "~ degradation_score_z + method_lineage_z +",
                                 paste(clinical_covars, collapse = " + "))), data = d)
      out <- extract_lm_term(fit, "method_lineage_z", node_catalog$node[i], lineage_catalog$lineage[j],
                             this_method, "Method-level lineage-node association")
      out[, method := this_method]
      out
    }))
  }))
}))
method_node_associations[, BH_q_within_node := p.adjust(p_value, method = "BH"), by = outcome]
method_node_associations[, BH_q_global := p.adjust(p_value, method = "BH")]
method_node_summary <- method_node_associations[, .(
  methods_n = .N,
  positive_methods_n = sum(standardized_beta > 0),
  positive_direction_fraction = mean(standardized_beta > 0),
  nominal_positive_methods_n = sum(standardized_beta > 0 & p_value < 0.05),
  BH_positive_methods_n = sum(standardized_beta > 0 & BH_q_within_node < 0.05),
  median_standardized_beta = median(standardized_beta),
  minimum_BH_q_within_node = min(BH_q_within_node)
), by = .(outcome, lineage = predictor)]

node_score_mat <- do.call(rbind, lapply(seq_len(nrow(node_catalog)), function(i) analysis_data[[node_catalog$column[i]]]))
rownames(node_score_mat) <- node_catalog$node; colnames(node_score_mat) <- sample_ids
node_group_effects <- rbindlist(list(
  fit_group_matrix(node_score_mat, analysis_data, model_label = "Clinical covariates"),
  fit_group_matrix(node_score_mat, analysis_data, extra_terms = "T_cell_context_z", model_label = "After T-cell context"),
  fit_group_matrix(node_score_mat, analysis_data, extra_terms = "APC_myeloid_context_z", model_label = "After APC/myeloid context"),
  fit_group_matrix(node_score_mat, analysis_data, extra_terms = c("T_cell_context_z", "APC_myeloid_context_z"), model_label = "After both lineage contexts")
))
baseline_effect <- node_group_effects[model == "Clinical covariates", .(feature, baseline_effect = standardized_effect)]
node_group_effects <- merge(node_group_effects, baseline_effect, by = "feature", all.x = TRUE)
node_group_effects[, attenuation_fraction := fifelse(abs(baseline_effect) > 1e-8,
                                                     1 - standardized_effect / baseline_effect, NA_real_)]

# -------------------------------------------------------------------------
# Candidate TME programs, with downstream-route genes removed before scoring
# -------------------------------------------------------------------------
hallmark <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
reactome <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))
gobp <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"))
stopifnot(uniqueN(hallmark$db_version) == 1L,
          identical(unique(hallmark$db_version), unique(reactome$db_version)),
          identical(unique(hallmark$db_version), unique(gobp$db_version)))
msigdb_version <- unique(hallmark$db_version)

tme_catalog <- data.table(
  program_order = 1:17,
  category = c(rep("Immune-supportive context", 7), rep("Inflammatory context", 4), rep("Tumor/stromal context", 6)),
  program = c(
    "Antigen cross-presentation", "MHC-II antigen presentation", "Myeloid dendritic-cell activation",
    "IL-12 production", "T-cell receptor signaling",
    "NK-cell activation", "Leukocyte cytotoxicity",
    "Inflammatory response", "TNF-alpha/NF-kappaB", "IL-6/JAK/STAT3", "Complement",
    "TGF-beta signaling", "Hypoxia", "Glycolysis", "Oxidative phosphorylation",
    "Epithelial-mesenchymal transition", "Angiogenesis"
  ),
  collection = c(
    "Reactome", "Reactome", "GO:BP", "GO:BP", "GO:BP", "GO:BP", "GO:BP",
    rep("Hallmark", 10)
  ),
  source_id = c(
    "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
    "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
    "GOBP_MYELOID_DENDRITIC_CELL_ACTIVATION",
    "GOBP_INTERLEUKIN_12_PRODUCTION",
    "GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY",
    "GOBP_NATURAL_KILLER_CELL_ACTIVATION",
    "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY",
    "HALLMARK_INFLAMMATORY_RESPONSE",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_COMPLEMENT",
    "HALLMARK_TGF_BETA_SIGNALING",
    "HALLMARK_HYPOXIA",
    "HALLMARK_GLYCOLYSIS",
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_ANGIOGENESIS"
  )
)

collection_tables <- list(Hallmark = hallmark, Reactome = reactome, `GO:BP` = gobp)
all_available <- unique(c(hallmark$gs_name, reactome$gs_name, gobp$gs_name))
if (!all(tme_catalog$source_id %in% all_available)) {
  stop("Missing candidate TME gene sets: ", paste(setdiff(tme_catalog$source_id, all_available), collapse = ", "))
}

get_manifest_genes <- function(score_name) {
  x <- score_manifest[score == score_name, scored_genes]
  if (length(x) != 1L) stop("Score manifest row not found: ", score_name)
  strsplit(x, ";", fixed = TRUE)[[1]]
}
cd28_genes <- unique(source_membership[source_key == "CD28" & expressed_in_PRE73 == TRUE, gene_symbol])
nfkb_regulon_genes <- unique(regulon[TF %in% c("NFKB1", "RELA") & expressed == TRUE & confidence %in% c("A", "B"), target])
ifng_core_genes <- get_manifest_genes("IFNG response core")
downstream_union <- unique(c(cd28_genes, nfkb_regulon_genes, ifng_core_genes, "IFNG"))

tme_scores <- list()
tme_gene_manifest <- rbindlist(lapply(seq_len(nrow(tme_catalog)), function(i) {
  row <- tme_catalog[i]
  tbl <- collection_tables[[row$collection]]
  original <- unique(tbl[gs_name == row$source_id, gene_symbol])
  expressed <- intersect(original, rownames(log_expression))
  scored <- setdiff(expressed, downstream_union)
  obj <- score_gene_set(scored, row$program, log_expression, minimum_genes = 5L)
  tme_scores[[row$program]] <<- obj$score
  data.table(
    program = row$program, category = row$category, collection = row$collection,
    source_id = row$source_id, original_genes_n = length(original),
    expressed_genes_n = length(expressed), removed_downstream_overlap_n = length(intersect(expressed, downstream_union)),
    scored_genes_n = length(obj$genes), scored_genes = paste(obj$genes, collapse = ";"),
    MSigDB_version = msigdb_version
  )
}))

tme_score_mat <- do.call(rbind, tme_scores)
rownames(tme_score_mat) <- names(tme_scores); colnames(tme_score_mat) <- sample_ids
for (nm in names(tme_scores)) analysis_data[, (paste0("tme_", match(nm, names(tme_scores)))) := tme_scores[[nm]]]
tme_group_effects <- fit_group_matrix(tme_score_mat, analysis_data, model_label = "Clinical covariates")
setnames(tme_group_effects, "feature", "program")
tme_group_effects[, edge := "Degradation High - Low"]

tme_route_associations <- rbindlist(lapply(seq_len(nrow(tme_catalog)), function(i) {
  program_col <- paste0("tme_", i)
  program_label <- tme_catalog$program[i]
  fit_cd28 <- lm(as.formula(paste("cd28_route_1 ~ degradation_score_z +", program_col, "+", paste(clinical_covars, collapse = " + "))), data = analysis_data)
  fit_nfkb <- lm(as.formula(paste("branch_nfkb ~ degradation_score_z + cd28_route_1 +", program_col, "+", paste(clinical_covars, collapse = " + "))), data = analysis_data)
  fit_ifng <- lm(as.formula(paste("IFNG_expression_z ~ degradation_score_z + cd28_route_1 + branch_nfkb +", program_col, "+", paste(clinical_covars, collapse = " + "))), data = analysis_data)
  fit_resp <- lm(as.formula(paste("IFNG_response_core_z ~ degradation_score_z + cd28_route_1 + branch_nfkb +", program_col, "+", paste(clinical_covars, collapse = " + "))), data = analysis_data)
  primary <- rbindlist(list(
    extract_lm_term(fit_cd28, program_col, "CD28 pathway expression", program_label, "Clinical covariates + degradation score", "TME-route association"),
    extract_lm_term(fit_nfkb, program_col, "NFKB1/RELA activity", program_label, "Clinical covariates + degradation score + CD28", "TME-route association"),
    extract_lm_term(fit_ifng, program_col, "IFNG expression", program_label, "Clinical covariates + degradation score + CD28 + NF-kappaB", "TME-route association"),
    extract_lm_term(fit_resp, program_col, "IFNG-response core", program_label, "Clinical covariates + degradation score + CD28 + NF-kappaB", "TME-route association")
  ))
  fit_cd28_t <- lm(as.formula(paste("cd28_route_1 ~ degradation_score_z +", program_col, "+ T_cell_context_z +", paste(clinical_covars, collapse = " + "))), data = analysis_data)
  fit_nfkb_t <- lm(as.formula(paste("branch_nfkb ~ degradation_score_z + cd28_route_1 +", program_col, "+ T_cell_context_z +", paste(clinical_covars, collapse = " + "))), data = analysis_data)
  fit_ifng_t <- lm(as.formula(paste("IFNG_expression_z ~ degradation_score_z + cd28_route_1 + branch_nfkb +", program_col, "+ T_cell_context_z +", paste(clinical_covars, collapse = " + "))), data = analysis_data)
  fit_resp_t <- lm(as.formula(paste("IFNG_response_core_z ~ degradation_score_z + cd28_route_1 + branch_nfkb +", program_col, "+ T_cell_context_z +", paste(clinical_covars, collapse = " + "))), data = analysis_data)
  sensitivity <- rbindlist(list(
    extract_lm_term(fit_cd28_t, program_col, "CD28 pathway expression", program_label, "After T-cell context", "TME-route association"),
    extract_lm_term(fit_nfkb_t, program_col, "NFKB1/RELA activity", program_label, "After T-cell context", "TME-route association"),
    extract_lm_term(fit_ifng_t, program_col, "IFNG expression", program_label, "After T-cell context", "TME-route association"),
    extract_lm_term(fit_resp_t, program_col, "IFNG-response core", program_label, "After T-cell context", "TME-route association")
  ))
  rbind(primary, sensitivity)
}))
tme_route_associations[, BH_q := p.adjust(p_value, method = "BH"), by = .(model, outcome)]

tme_group_long <- tme_group_effects[, .(
  program, outcome = "Degradation High - Low", predictor = "Degradation group",
  model, standardized_beta = standardized_effect, standard_error, CI_low, CI_high,
  p_value, BH_q
)]
tme_primary_long <- tme_route_associations[model != "After T-cell context", .(
  program = predictor, outcome, predictor = "TME program", model,
  standardized_beta, standard_error, CI_low, CI_high, p_value, BH_q
)]
tme_context_evidence <- rbind(tme_group_long, tme_primary_long, fill = TRUE)
tme_context_evidence <- merge(tme_context_evidence, tme_catalog[, .(program, category, program_order)], by = "program", all.x = TRUE)

# Strict descriptive prioritization: a candidate must be higher in
# Degradation-High and positively track CD28, NFKB1/RELA, and the IFNG-response
# core. A second count records which links remain after T-cell-context
# sensitivity. This is ranking, not a new causal test.
get_tme_link <- function(program_name, outcome_name, model_filter) {
  x <- tme_route_associations[predictor == program_name & outcome == outcome_name & model == model_filter]
  if (nrow(x) != 1L) stop("Unexpected TME link multiplicity")
  x
}
tme_candidate_ranking <- rbindlist(lapply(tme_catalog$program, function(program_name) {
  g <- tme_group_effects[program == program_name]
  p_cd28 <- get_tme_link(program_name, "CD28 pathway expression", "Clinical covariates + degradation score")
  p_nfkb <- get_tme_link(program_name, "NFKB1/RELA activity", "Clinical covariates + degradation score + CD28")
  p_resp <- get_tme_link(program_name, "IFNG-response core", "Clinical covariates + degradation score + CD28 + NF-kappaB")
  t_cd28 <- get_tme_link(program_name, "CD28 pathway expression", "After T-cell context")
  t_nfkb <- get_tme_link(program_name, "NFKB1/RELA activity", "After T-cell context")
  t_resp <- get_tme_link(program_name, "IFNG-response core", "After T-cell context")
  primary_pass <- c(g$standardized_effect > 0 & g$BH_q < 0.05,
                    p_cd28$standardized_beta > 0 & p_cd28$BH_q < 0.05,
                    p_nfkb$standardized_beta > 0 & p_nfkb$BH_q < 0.05,
                    p_resp$standardized_beta > 0 & p_resp$BH_q < 0.05)
  tcell_pass <- c(t_cd28$standardized_beta > 0 & t_cd28$BH_q < 0.05,
                  t_nfkb$standardized_beta > 0 & t_nfkb$BH_q < 0.05,
                  t_resp$standardized_beta > 0 & t_resp$BH_q < 0.05)
  data.table(
    program = program_name,
    High_vs_Low_effect = g$standardized_effect, High_vs_Low_BH_q = g$BH_q,
    CD28_beta = p_cd28$standardized_beta, CD28_BH_q = p_cd28$BH_q,
    NFKB_beta = p_nfkb$standardized_beta, NFKB_BH_q = p_nfkb$BH_q,
    IFNG_response_beta = p_resp$standardized_beta, IFNG_response_BH_q = p_resp$BH_q,
    primary_positive_evidence_count = sum(primary_pass),
    Tcell_context_positive_link_count = sum(tcell_pass),
    strict_4_of_4_primary = all(primary_pass),
    strict_3_of_3_after_Tcell_context = all(tcell_pass)
  )
}))
tme_candidate_ranking <- merge(tme_candidate_ranking, tme_catalog[, .(program, category, program_order)], by = "program", all.x = TRUE)
tme_candidate_ranking <- tme_candidate_ranking[order(-strict_4_of_4_primary, -strict_3_of_3_after_Tcell_context,
                                                     -primary_positive_evidence_count, High_vs_Low_BH_q, program_order)]

vif_audit <- rbindlist(list(
  data.table(
    model = "Both lineage contexts", term = c("degradation_score_z", "T_cell_context_z", "APC_myeloid_context_z"),
    VIF = c(
      vif_one("degradation_score_z", c("T_cell_context_z", "APC_myeloid_context_z", clinical_covars), analysis_data),
      vif_one("T_cell_context_z", c("degradation_score_z", "APC_myeloid_context_z", clinical_covars), analysis_data),
      vif_one("APC_myeloid_context_z", c("degradation_score_z", "T_cell_context_z", clinical_covars), analysis_data)
    )
  ),
  rbindlist(lapply(seq_len(nrow(tme_catalog)), function(i) {
    term <- paste0("tme_", i)
    data.table(
      model = "TME -> IFNG response after T-cell context",
      term = tme_catalog$program[i],
      VIF = vif_one(term, c("degradation_score_z", "cd28_route_1", "branch_nfkb", "T_cell_context_z", clinical_covars), analysis_data)
    )
  }))
))
vif_audit[, high_collinearity_flag := VIF >= 5]

# -------------------------------------------------------------------------
# Tables, validation, and report
# -------------------------------------------------------------------------
input_manifest <- data.table(
  path = required_inputs,
  role = c("FPKM expression", "PRE73 metadata", "route sample scores", "route score genes",
           "Reactome source membership", "seven-method deconvolution estimates", "signed DoRothEA A/B regulon"),
  sha256 = vapply(required_inputs, sha256_file, character(1))
)
method_contract <- data.table(
  item = c(
    "Cohort", "Exposure", "Cell-attribution input", "Lineage consensus",
    "Candidate TME programs", "Overlap control", "Primary covariates",
    "Multiple testing", "Causal boundary"
  ),
  specification = c(
    "73 PRE-treatment TIGER melanoma bulk-RNA tumors",
    "Existing SARDH/PIPOX degradation transcriptional score; median Low n=37 / High n=36",
    "CIBERSORT, EPIC, ESTIMATE, MCPcounter, TIMER, quanTIseq and xCell output; ESTIMATE/xCell generic immune-score rows excluded from lineage attribution",
    "Within lineage, average duplicate subtypes per method, z-score per method, average across methods, z-score again",
    "17 predeclared Hallmark/Reactome/GO:BP programs",
    "Remove expressed CD28-set, NFKB1/RELA DoRothEA A/B targets, IFNG-response-core genes and IFNG from every TME program before scoring",
    "Therapy + clinical response + age + sex; degradation score retained in association models",
    "BH within each outcome/model family; global BH additionally reported for lineage-node screen",
    "Bulk/transcript-level cross-sectional associations; no cell-of-origin or induction/causation claim"
  )
)
package_versions <- data.table(
  package = c("R", required_packages),
  version = c(as.character(getRversion()), vapply(required_packages, function(p) as.character(packageVersion(p)), character(1)))
)

validation <- data.table(
  check = c(
    "PRE sample count is 73", "Degradation groups are 37 Low / 36 High",
    "Every lineage has 73 samples", "At least six independent methods contribute to CD8 consensus",
    "No generic ESTIMATE score is used in lineage consensus", "All six route nodes are complete",
    "All candidate TME programs retain at least five genes", "TME scored genes exclude downstream union",
    "Lineage-node screen has expected 60 rows", "Method-level lineage audit is complete",
    "TME candidate ranking has expected 17 rows", "TME primary evidence has expected 85 rows",
    "TME sensitivity screen has expected 68 rows", "All reported coefficients are finite",
    "Most parameterized VIF is below 10", "Reader-facing labels contain no +ESTIMATE",
    "MSigDB version is singular"
  ),
  passed = c(
    nrow(analysis_data) == 73L,
    sum(analysis_data$degradation_group == "Low") == 37L && sum(analysis_data$degradation_group == "High") == 36L,
    all(lineage_consensus[, uniqueN(sample_id), by = lineage]$V1 == 73L),
    lineage_methods[lineage == "CD8 T cells", contributing_methods_n] >= 6L,
    !any(deconv_lineage$method == "ESTIMATE"),
    all(complete.cases(analysis_data[, ..node_catalog$column])),
    all(tme_gene_manifest$scored_genes_n >= 5L),
    all(vapply(strsplit(tme_gene_manifest$scored_genes, ";", fixed = TRUE), function(x) length(intersect(x, downstream_union)) == 0L, logical(1))),
    nrow(cell_node_associations) == 60L,
    nrow(method_node_associations) == sum(lineage_catalog$contributing_methods_n) * nrow(node_catalog),
    nrow(tme_candidate_ranking) == 17L,
    nrow(tme_context_evidence) == 85L,
    nrow(tme_route_associations[model == "After T-cell context"]) == 68L,
    all(is.finite(c(cell_node_associations$standardized_beta, node_group_effects$standardized_effect,
                    tme_route_associations$standardized_beta, tme_group_effects$standardized_effect))),
    max(vif_audit$VIF) < 10,
    !any(grepl("\\+ESTIMATE", c(node_catalog$node, lineage_catalog$lineage, tme_catalog$program))),
    length(msigdb_version) == 1L
  )
)
if (!all(validation$passed)) stop("Validation failed: ", paste(validation[passed == FALSE, check], collapse = "; "),
                                  "; max VIF = ", sprintf("%.3f", max(vif_audit$VIF)))

fwrite(input_manifest, file.path(table_dir, "00_input_manifest.csv"))
fwrite(method_contract, file.path(table_dir, "00_method_contract.csv"))
fwrite(package_versions, file.path(table_dir, "00_package_versions.csv"))
fwrite(lineage_methods, file.path(table_dir, "01_lineage_consensus_method_manifest.csv"))
fwrite(lineage_composite_manifest, file.path(table_dir, "02_lineage_composite_manifest.csv"))
fwrite(deconv_method_lineage, file.path(table_dir, "03_method_level_lineage_values.csv"))
fwrite(lineage_consensus, file.path(table_dir, "04_cross_method_lineage_consensus_scores.csv"))
fwrite(cell_node_associations, file.path(table_dir, "05_lineage_associations_with_route_nodes.csv"))
fwrite(method_node_associations, file.path(table_dir, "05b_method_level_lineage_node_associations.csv"))
fwrite(method_node_summary, file.path(table_dir, "05c_method_agreement_summary.csv"))
fwrite(node_group_effects, file.path(table_dir, "06_route_node_High_vs_Low_lineage_sensitivity.csv"))
fwrite(tme_catalog, file.path(table_dir, "07_candidate_TME_program_catalog.csv"))
fwrite(tme_gene_manifest, file.path(table_dir, "08_candidate_TME_program_gene_manifest_disjoint.csv"))
fwrite(tme_group_effects, file.path(table_dir, "09_candidate_TME_program_High_vs_Low.csv"))
fwrite(tme_route_associations, file.path(table_dir, "10_candidate_TME_program_route_associations.csv"))
fwrite(tme_context_evidence, file.path(table_dir, "11_candidate_TME_context_evidence_primary.csv"))
fwrite(tme_candidate_ranking, file.path(table_dir, "11b_candidate_TME_context_strict_ranking.csv"))
fwrite(vif_audit, file.path(table_dir, "12_VIF_audit.csv"))
fwrite(validation, file.path(table_dir, "13_validation_checks.csv"))
fwrite(analysis_data[, c("sample_id", "degradation_group", "Degradation_score", node_catalog$column,
                         "T_cell_context_z", "APC_myeloid_context_z", paste0("tme_", seq_len(nrow(tme_catalog)))), with = FALSE],
       file.path(table_dir, "14_sample_level_scores.csv"))

# -------------------------------------------------------------------------
# Reader-facing figures
# -------------------------------------------------------------------------
red <- "#C83B3D"; teal <- "#1D9E92"
theme_reader <- theme_classic(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(color = "#4D4D4D", size = 12),
        axis.title = element_text(face = "bold"),
        strip.background = element_rect(fill = "white", color = "#333333"),
        strip.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"))

cell_plot_dt <- copy(cell_node_associations)
cell_plot_dt <- merge(cell_plot_dt, node_catalog[, .(outcome = node, node_order)], by = "outcome", all.x = TRUE)
cell_plot_dt <- merge(cell_plot_dt, lineage_catalog[, .(predictor = lineage, lineage_order)], by = "predictor", all.x = TRUE)
cell_plot_dt[, lineage_label := paste0(predictor, " (", contributing_methods_n, " method", ifelse(contributing_methods_n == 1, "", "s"), ")")]
cell_plot_dt[, node_label := factor(outcome, levels = node_catalog$node)]
lineage_label_order <- cell_plot_dt[order(lineage_order), unique(lineage_label)]
cell_plot_dt[, lineage_label := factor(lineage_label, levels = rev(lineage_label_order))]

p_cell <- ggplot(cell_plot_dt, aes(node_label, lineage_label, fill = standardized_beta)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_point(data = cell_plot_dt[BH_q_global < 0.05], shape = 21, size = 2.6, fill = "black", color = "black") +
  scale_fill_gradient2(low = teal, mid = "white", high = red, midpoint = 0,
                       name = "Adjusted association\ncoefficient") +
  labs(
    title = "Which inferred immune-cell contexts track the CD28–IFNG axis?",
    subtitle = "Cross-method deconvolution consensus; models adjust for degradation score and clinical covariates. Black dot: global BH q < 0.05.",
    x = NULL, y = NULL,
    caption = "Associations in bulk RNA support lineage context but do not prove the cellular source of a transcript or pathway."
  ) + theme_reader +
  theme(axis.text.x = element_text(angle = 25, hjust = 1), panel.grid = element_blank(),
        plot.caption = element_text(hjust = 0, color = "#555555"))
ggsave(file.path(figure_dir, "Fig_Cell_Lineage_Associations_with_CD28_IFNG_Nodes.png"), p_cell, width = 15, height = 9, dpi = 300)
ggsave(file.path(figure_dir, "Fig_Cell_Lineage_Associations_with_CD28_IFNG_Nodes.pdf"), p_cell, width = 15, height = 9)

model_levels <- c("Clinical covariates", "After T-cell context", "After APC/myeloid context", "After both lineage contexts")
model_labels <- c(
  "High vs Low\n(clinical covariates)", "After T-cell\nabundance context",
  "After APC/myeloid\nabundance context", "After both\nlineage contexts"
)
node_plot_dt <- copy(node_group_effects)
node_plot_dt[, feature := factor(feature, levels = rev(node_catalog$node))]
node_plot_dt[, model := factor(model, levels = model_levels, labels = model_labels)]
p_node <- ggplot(node_plot_dt, aes(model, feature, fill = standardized_effect)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.2f%s", standardized_effect, ifelse(BH_q < 0.05, "*", ""))), size = 4) +
  scale_fill_gradient2(low = teal, mid = "white", high = red, midpoint = 0,
                       name = "Degradation High-minus-Low\neffect (SD)") +
  labs(
    title = "How much of the Degradation-High signal follows inferred lineage abundance?",
    subtitle = "Asterisk: BH q < 0.05 within each model. Context adjustment is a sensitivity analysis, not proof of confounding or mediation.",
    x = NULL, y = NULL
  ) + theme_reader +
  theme(axis.text.x = element_text(size = 11), panel.grid = element_blank())
ggsave(file.path(figure_dir, "Fig_Degradation_HighLow_Node_Effects_After_Lineage_Context.png"), p_node, width = 13, height = 7.8, dpi = 300)
ggsave(file.path(figure_dir, "Fig_Degradation_HighLow_Node_Effects_After_Lineage_Context.pdf"), p_node, width = 13, height = 7.8)

tme_plot_dt <- copy(tme_context_evidence)
outcome_levels <- c("Degradation High - Low", "CD28 pathway expression", "NFKB1/RELA activity", "IFNG expression", "IFNG-response core")
outcome_labels <- c("Higher in\nDegradation-High", "Tracks CD28\npathway", "Tracks NFKB1/RELA\nactivity", "Tracks IFNG\nexpression", "Tracks IFNG-response\ncore")
tme_plot_dt[, outcome := factor(outcome, levels = outcome_levels, labels = outcome_labels)]
tme_plot_dt[, program := factor(program, levels = rev(tme_catalog$program))]
tme_plot_dt[, neglog_q := pmin(-log10(pmax(BH_q, .Machine$double.xmin)), 20)]
tme_plot_dt[, significant := BH_q < 0.05]
p_tme <- ggplot(tme_plot_dt, aes(outcome, program)) +
  geom_point(aes(fill = standardized_beta, size = neglog_q, alpha = significant), shape = 21, color = "#333333", stroke = 0.6) +
  scale_fill_gradient2(low = teal, mid = "white", high = red, midpoint = 0,
                       name = "Adjusted\ncoefficient") +
  scale_size_continuous(name = expression(-log[10](BH~q)), range = c(1.8, 7)) +
  scale_alpha_manual(values = c(`FALSE` = 0.35, `TRUE` = 1), guide = "none") +
  labs(
    title = "Tumor-microenvironment programs associated with the CD28–NF-kappaB–IFNG axis",
    subtitle = "Every program score excludes downstream-route genes. Filled/opaque points pass BH q < 0.05 within the displayed outcome family.",
    x = NULL, y = NULL,
    caption = "These are candidate associated contexts in bulk PRE-treatment tumors; they do not establish an inducing mechanism."
  ) + theme_reader +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), panel.grid = element_blank(),
        plot.caption = element_text(hjust = 0, color = "#555555"))
ggsave(file.path(figure_dir, "Fig_TME_Context_Candidates_for_CD28_IFNG_Axis.png"), p_tme, width = 15, height = 11, dpi = 300)
ggsave(file.path(figure_dir, "Fig_TME_Context_Candidates_for_CD28_IFNG_Axis.pdf"), p_tme, width = 15, height = 11)

sig_cell <- cell_node_associations[BH_q_global < 0.05][order(outcome, -abs(standardized_beta))]
sig_tme_group <- tme_group_effects[BH_q < 0.05][order(BH_q)]
sig_tme_route <- tme_route_associations[model != "After T-cell context" & BH_q < 0.05][order(outcome, BH_q)]
sig_tme_route_t <- tme_route_associations[model == "After T-cell context" & BH_q < 0.05][order(outcome, BH_q)]
strict_tme_candidates <- tme_candidate_ranking[strict_4_of_4_primary == TRUE][order(High_vs_Low_BH_q)]
key_method_agreement <- method_node_summary[
  (outcome == "SPI1/PU.1 activity" & lineage %in% c("Dendritic cells", "Macrophages")) |
    (outcome == "CD80/CD86 expression" & lineage %in% c("Dendritic cells", "Macrophages")) |
    (outcome == "CD28 pathway expression" & lineage %in% c("Pan T cells", "CD8 T cells", "CD4 T cells")) |
    (outcome %in% c("IFNG expression", "IFNG-response core") & lineage %in% c("CD8 T cells", "NK cells", "Cytotoxic lymphocytes"))
][order(outcome, lineage)]

fmt_lines <- function(dt, formatter, empty_text) {
  if (!nrow(dt)) return(paste0("- ", empty_text))
  paste(vapply(seq_len(nrow(dt)), function(i) paste0("- ", formatter(dt[i])), character(1)), collapse = "\n")
}
report_lines <- c(
  "# TIGER PRE73 CD28–IFNG cell-origin and TME-context analysis", "",
  "## Integrity boundary", "",
  "This is a bulk PRE-treatment tumor RNA-seq analysis. The results can identify immune-lineage and tumor-microenvironment contexts that covary with the transcript programs, but cannot assign a transcript to a specific cell or prove that a TME program induces the route. Degradation-High, not Degradation-Low, is the direction of the previously observed CD28/NF-kappaB/IFNG-associated signal.", "",
  "## Cross-method lineage associations", "",
  fmt_lines(sig_cell,
            function(x) sprintf("%s with %s: beta = %.3f, global BH q = %.3g (%d contributing method%s).",
                                x$predictor, x$outcome, x$standardized_beta, x$BH_q_global,
                                x$contributing_methods_n, ifelse(x$contributing_methods_n == 1, "", "s")),
            "No lineage-node pair passed the global BH q < 0.05 screen."), "",
  "Method-level direction audit for key lineage hypotheses:", "",
  fmt_lines(key_method_agreement,
            function(x) sprintf("%s with %s: %d/%d methods positive; median beta = %.3f; %d methods pass node-wise BH q < 0.05.",
                                x$lineage, x$outcome, x$positive_methods_n, x$methods_n,
                                x$median_standardized_beta, x$BH_positive_methods_n),
            "No key lineage hypotheses were available for the method-level audit."), "",
  "## Degradation-High effects after lineage-context sensitivity", "",
  paste0("Exact effects and q values are in `tables/06_route_node_High_vs_Low_lineage_sensitivity.csv`. Adjustment for a lineage context cannot distinguish confounding from mediation. Maximum VIF across audited models was ", sprintf("%.2f", max(vif_audit$VIF)), "; VIF values above 5 are flagged as potentially unstable rather than treated as confirmatory."), "",
  "## Candidate TME programs", "",
  "Programs significantly different between Degradation-High and Low:", "",
  fmt_lines(sig_tme_group,
            function(x) sprintf("%s: High-minus-Low effect = %.3f, BH q = %.3g.", x$program, x$standardized_effect, x$BH_q),
            "No predeclared program passed BH q < 0.05."), "",
  "Strict context candidates (positive High-vs-Low effect plus positive links to CD28, NFKB1/RELA, and the IFNG-response core; all four BH q < 0.05):", "",
  fmt_lines(strict_tme_candidates,
            function(x) sprintf("%s: primary evidence %d/4; after T-cell-context sensitivity %d/3 positive links remain.",
                                x$program, x$primary_positive_evidence_count, x$Tcell_context_positive_link_count),
            "No candidate met the strict 4-of-4 descriptive rule."), "",
  "Gene-disjoint programs significantly associated with route nodes in the primary model:", "",
  fmt_lines(sig_tme_route,
            function(x) sprintf("%s -> %s: beta = %.3f, BH q = %.3g.", x$predictor, x$outcome, x$standardized_beta, x$BH_q),
            "No program-to-route association passed BH q < 0.05."), "",
  "Associations remaining significant after the named T-cell-context sensitivity:", "",
  fmt_lines(sig_tme_route_t,
            function(x) sprintf("%s -> %s: beta = %.3f, BH q = %.3g.", x$predictor, x$outcome, x$standardized_beta, x$BH_q),
            "No association remained significant after T-cell-context sensitivity."), "",
  "## Interpretation", "",
  "- CD28 transcript/pathway and IFNG signals should be described as immune-context signals in bulk tumor, not as tumor-cell-autonomous activity unless orthogonal cell-resolved data are added.",
  "- SPI1/PU.1 and CD80/CD86 plausibly mark an APC/myeloid side of the association, whereas CD28 and much of IFNG production plausibly mark T/NK-cell biology; the exact cell assignment remains unproven here.",
  "- Candidate TME programs are prioritized by association only. A mechanistic statement requires spatial, single-cell, sorted-cell, perturbation, or protein/phospho-level validation.", "",
  "## Reproducibility", "",
  paste0("- MSigDB release: ", msigdb_version, "."),
  "- All candidate TME scores exclude the union of expressed CD28-set genes, NFKB1/RELA signed DoRothEA A/B targets, the IFNG-response core, and IFNG itself.",
  paste0("- Validation: ", sum(validation$passed), "/", nrow(validation), " checks passed."),
  "- Input SHA-256 hashes are recorded in `tables/00_input_manifest.csv`."
)
writeLines(report_lines, file.path(output_root, "REPORT.md"))

figure_files <- list.files(figure_dir, full.names = TRUE)
figure_hashes <- data.table(file = basename(figure_files), sha256 = vapply(figure_files, sha256_file, character(1)))
fwrite(figure_hashes, file.path(table_dir, "15_figure_hashes.csv"))

run_summary <- c(
  paste0("Completed: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Analysis: TIGER PRE73 CD28-IFNG cell-origin and TME context",
  paste0("Samples: ", nrow(analysis_data), " (Low ", sum(analysis_data$degradation_group == "Low"), ", High ", sum(analysis_data$degradation_group == "High"), ")"),
  paste0("Lineages screened: ", nrow(lineage_catalog)),
  paste0("Candidate TME programs: ", nrow(tme_catalog)),
  paste0("Significant lineage-node pairs (global BH q < 0.05): ", nrow(sig_cell)),
  paste0("Significant TME High-vs-Low programs: ", nrow(sig_tme_group)),
  paste0("Significant primary TME-route associations: ", nrow(sig_tme_route)),
  paste0("Strict 4-of-4 TME context candidates: ", nrow(strict_tme_candidates)),
  paste0("Maximum audited VIF: ", sprintf("%.3f", max(vif_audit$VIF))),
  paste0("Validation passed: ", sum(validation$passed), "/", nrow(validation))
)
writeLines(run_summary, file.path(log_dir, "run_summary.txt"))
cat(paste(run_summary, collapse = "\n"), "\n")
