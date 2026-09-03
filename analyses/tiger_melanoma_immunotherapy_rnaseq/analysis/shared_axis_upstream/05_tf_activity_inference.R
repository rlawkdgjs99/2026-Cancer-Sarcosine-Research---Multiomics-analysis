#!/usr/bin/env Rscript

# Candidate transcriptional regulators of the CD28 costimulation program in
# TIGER melanoma PRE-only bulk RNA-seq.
#
# Scientific scope
# ----------------
# 1. DoRothEA A/B signed TF-target interactions are used to infer sample-level
#    TF regulon activities by the univariate linear-model (ULM) statistic used
#    by decoupleR (correlation converted to a t statistic).
# 2. Candidate regulators are prioritized only when their curated targets
#    overlap the *full* Reactome CD28 pathway memberships. This structural
#    overlap is independent of the observed High-vs-Low leading edge.
# 3. High vs Low is the only biological grouping displayed in figures.
#    Clinical covariates are handled in the model and described in captions;
#    ESTIMATE immune/stromal adjustment is retained only as a sensitivity table.
#
# This observational bulk-RNA analysis identifies TF-regulated expression
# programs associated with the CD28 pathway. It does not establish a causal
# TF -> CD28 -> IFNG sequence or a T-cell-intrinsic mechanism.

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

prior_root <- file.path(
  analysis_root, "results", "TIGER_PRE73_Degradation_HighLow_GSEA_GO_BP_26.09.02"
)
integrated_root <- file.path(
  analysis_root, "results", "TIGER_PRE73_ThreeCollection_Immune_TME_Integrated_26.09.02"
)
deconv_root <- file.path(
  analysis_root, "results", "TIGER_PRE73_Degradation_HighLow_Deconvolution_7methods_26.09.02"
)
ifng_root <- file.path(
  analysis_root, "results", "TIGER_PRE73_IFNG_Axis_Focused_26.09.02"
)
output_root <- file.path(
  analysis_root, "results", "TIGER_PRE73_CD28_Upstream_TF_26.09.02"
)
table_dir <- file.path(output_root, "tables")
figure_dir <- file.path(output_root, "figures")
log_dir <- file.path(output_root, "logs")
resource_dir <- file.path(output_root, "resources")
for (d in c(table_dir, figure_dir, log_dir, resource_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

used_data_root <- normalizePath(file.path(analysis_root, "..", "사용데이터_모음"))
expression_path <- file.path(
  used_data_root, "Source_Input", "TIGER_PRJEB23709_expression_FPKM_gene_by_sample.csv"
)
metadata_path <- file.path(
  used_data_root, "Analysis_Ready", "TIGER_PRE73_Fig4bc_analysis_data.csv"
)
primary_gene_path <- file.path(prior_root, "tables", "02_limma_adjusted_High_vs_Low_all_genes.csv")
primary_reactome_path <- file.path(integrated_root, "tables", "01_fgsea_Reactome_adjusted_High_vs_Low.csv")
tme_gene_path <- file.path(integrated_root, "tables", "08_limma_TME_adjusted_High_vs_Low_all_genes.csv")
estimate_path <- file.path(deconv_root, "raw_deconvolution", "deconvolution_ESTIMATE.csv")
c3_tft_path <- file.path(ifng_root, "tables", "08_C3_TFT_legacy_all_results.csv")
regulon_path <- file.path(resource_dir, "OmniPath_DoRothEA_ABC_human_2026-09-02.tsv")

required_inputs <- c(
  expression_path, metadata_path, primary_gene_path, primary_reactome_path,
  tme_gene_path, estimate_path, c3_tft_path, regulon_path
)
for (path in required_inputs) if (!file.exists(path)) stop("Missing input: ", path)

sha256_file <- function(path) {
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  strsplit(out, "[[:space:]]+")[[1]][1]
}

expr_dt <- fread(expression_path, check.names = FALSE)
meta <- fread(metadata_path, check.names = FALSE)
primary_genes <- fread(primary_gene_path, check.names = FALSE)
tme_genes <- fread(tme_gene_path, check.names = FALSE)
primary_reactome <- fread(primary_reactome_path, check.names = FALSE)
estimate <- fread(estimate_path, check.names = FALSE)
c3_tft <- fread(c3_tft_path, check.names = FALSE)
dorothea_raw <- fread(regulon_path, sep = "\t", quote = "", check.names = FALSE)

stopifnot(
  nrow(meta) == 73L,
  uniqueN(meta$sample_id) == 73L,
  all(meta$timepoint == "PRE"),
  all(c(
    "sample_id", "therapy_short", "response_group", "age", "gender",
    "Degradation_score"
  ) %in% names(meta)),
  all(c("gene_symbol", "moderated_t") %in% names(primary_genes)),
  all(c("gene_symbol", "moderated_t") %in% names(tme_genes)),
  all(c("pathway", "NES", "padj", "leadingEdge") %in% names(primary_reactome)),
  all(c("ID", "ImmuneScore_estimate", "StromalScore_estimate") %in% names(estimate)),
  all(c(
    "source_genesymbol", "target_genesymbol", "is_stimulation",
    "is_inhibition", "dorothea_level"
  ) %in% names(dorothea_raw))
)

gene_symbols <- trimws(expr_dt[[1]])
stopifnot(!anyNA(gene_symbols), all(nzchar(gene_symbols)), !anyDuplicated(gene_symbols))
stopifnot(all(meta$sample_id %in% names(expr_dt)))
sample_ids <- meta$sample_id
fpkm <- as.matrix(expr_dt[, ..sample_ids])
storage.mode(fpkm) <- "double"
rownames(fpkm) <- gene_symbols
colnames(fpkm) <- sample_ids
stopifnot(ncol(fpkm) == 73L, !anyNA(fpkm), all(is.finite(fpkm)), all(fpkm >= 0))

# Match the validated preceding analysis filter.
minimum_samples <- ceiling(0.10 * ncol(fpkm))
keep_expression <- rowSums(fpkm >= 1) >= minimum_samples
log_expression <- log2(fpkm[keep_expression, , drop = FALSE] + 1)
keep_variance <- apply(log_expression, 1, var) > 0
log_expression <- log_expression[keep_variance, , drop = FALSE]
stopifnot(setequal(rownames(log_expression), primary_genes$gene_symbol))

degradation_median <- median(meta$Degradation_score)
meta[, degradation_group := factor(
  ifelse(Degradation_score > degradation_median, "High", "Low"),
  levels = c("Low", "High")
)]
stopifnot(sum(meta$degradation_group == "Low") == 37L)
stopifnot(sum(meta$degradation_group == "High") == 36L)

estimate <- estimate[match(meta$sample_id, ID)]
stopifnot(identical(estimate$ID, meta$sample_id))
meta[, therapy_short := relevel(factor(therapy_short), ref = "antiPD1")]
meta[, response_group := relevel(factor(response_group), ref = "NR")]
meta[, gender := relevel(factor(gender), ref = "Male")]
meta[, age_z := as.numeric(scale(age))]
meta[, immune_score_z := as.numeric(scale(estimate$ImmuneScore_estimate))]
meta[, stromal_score_z := as.numeric(scale(estimate$StromalScore_estimate))]

# -------------------------------------------------------------------------
# High-confidence signed DoRothEA regulon
# -------------------------------------------------------------------------
as_flag <- function(x) tolower(as.character(x)) == "true"
dorothea_raw[, `:=`(
  stimulation_flag = as_flag(is_stimulation),
  inhibition_flag = as_flag(is_inhibition),
  confidence_rank = fifelse(grepl("(^|;)A($|;)", dorothea_level), 1L,
                            fifelse(grepl("(^|;)B($|;)", dorothea_level), 2L, NA_integer_))
)]

# Unknown or contradictory directions are not suitable for signed activity.
dorothea_signed <- dorothea_raw[
  !is.na(confidence_rank) & stimulation_flag != inhibition_flag,
  .(
    TF = trimws(source_genesymbol), target = trimws(target_genesymbol),
    mor = fifelse(stimulation_flag, 1, -1),
    confidence_rank, dorothea_level, sources, references
  )
]
dorothea_signed <- dorothea_signed[nzchar(TF) & nzchar(target)]

# Resolve duplicate TF-target rows using the best confidence. Conflicting signs
# at the same best confidence are excluded rather than guessed.
dorothea_signed[, best_rank := min(confidence_rank), by = .(TF, target)]
dorothea_signed <- dorothea_signed[confidence_rank == best_rank]
regulon <- dorothea_signed[, .(
  mor_values_n = uniqueN(mor),
  mor = if (uniqueN(mor) == 1L) unique(mor) else NA_real_,
  confidence = if (min(confidence_rank) == 1L) "A" else "B",
  sources = paste(sort(unique(unlist(strsplit(sources, ";", fixed = TRUE)))), collapse = ";"),
  references = paste(sort(unique(unlist(strsplit(references, ";", fixed = TRUE)))), collapse = ";")
), by = .(TF, target)]
regulon <- regulon[!is.na(mor)]
regulon <- unique(regulon, by = c("TF", "target"))
regulon[, expressed := target %in% rownames(log_expression)]

min_targets <- 10L
tf_coverage <- regulon[expressed == TRUE, .(
  expressed_targets_n = uniqueN(target),
  activating_targets_n = uniqueN(target[mor > 0]),
  repressing_targets_n = uniqueN(target[mor < 0]),
  confidence_A_targets_n = uniqueN(target[confidence == "A"]),
  confidence_B_targets_n = uniqueN(target[confidence == "B"])
), by = TF]
eligible_tfs <- tf_coverage[expressed_targets_n >= min_targets, TF]
regulon_tested <- regulon[expressed == TRUE & TF %in% eligible_tfs]
stopifnot(!anyDuplicated(regulon_tested[, .(TF, target)]), length(eligible_tfs) > 50L)

fwrite(regulon, file.path(table_dir, "01_DoRothEA_AB_signed_regulon_all.csv"))
fwrite(tf_coverage[order(-expressed_targets_n)], file.path(table_dir, "02_TF_regulon_coverage.csv"))

# -------------------------------------------------------------------------
# Reactome CD28 pathway memberships and regulator-target overlap
# -------------------------------------------------------------------------
cd28_catalog <- data.table(
  pathway = c(
    "REACTOME_CO_STIMULATION_BY_CD28",
    "REACTOME_CD28_DEPENDENT_PI3K_AKT_SIGNALING"
  ),
  pathway_label = c("CD28 co-stimulation", "CD28-dependent PI3K/AKT signaling"),
  pathway_order = 1:2
)
reactome_membership <- as.data.table(msigdbr(
  species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"
))
stopifnot(uniqueN(reactome_membership$db_version) == 1L)
msigdb_version <- unique(reactome_membership$db_version)
cd28_membership <- unique(reactome_membership[
  gs_name %in% cd28_catalog$pathway,
  .(
    pathway = gs_name, gene_symbol,
    Reactome_ID = gs_exact_source, pathway_url = gs_url,
    description = gs_description, MSigDB_version = db_version
  )
])
cd28_membership <- merge(cd28_membership, cd28_catalog, by = "pathway", all.x = TRUE)
cd28_union <- unique(cd28_membership$gene_symbol)
cd28_union_expressed <- intersect(cd28_union, rownames(log_expression))
stopifnot(length(cd28_union_expressed) >= 20L)

primary_cd28 <- merge(
  cd28_catalog,
  primary_reactome[, .(pathway, NES, pval, BH_q = padj, size, leadingEdge)],
  by = "pathway", all.x = TRUE
)
stopifnot(!anyNA(primary_cd28$NES))

membership_manifest <- copy(cd28_membership)
membership_manifest[, expressed_in_PRE73 := gene_symbol %in% rownames(log_expression)]
fwrite(membership_manifest, file.path(table_dir, "03_CD28_pathway_gene_manifest.csv"))
fwrite(primary_cd28, file.path(table_dir, "04_CD28_pathway_primary_GSEA_summary.csv"))

ora_universe <- rownames(log_expression)
ora_one <- function(tf, pathway_name, pathway_genes) {
  tf_targets <- unique(regulon_tested[TF == tf, target])
  pathway_genes <- intersect(unique(pathway_genes), ora_universe)
  overlap <- intersect(tf_targets, pathway_genes)
  a <- length(overlap)
  b <- length(setdiff(tf_targets, pathway_genes))
  c <- length(setdiff(pathway_genes, tf_targets))
  d <- length(ora_universe) - a - b - c
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2L), alternative = "greater")
  data.table(
    TF = tf, pathway = pathway_name,
    TF_targets_n = length(tf_targets), pathway_genes_n = length(pathway_genes),
    overlap_n = a, overlap_genes = paste(sort(overlap), collapse = ";"),
    odds_ratio = unname(ft$estimate), p_value = ft$p.value
  )
}

ora_results <- rbindlist(lapply(eligible_tfs, function(tf) {
  rbindlist(list(
    ora_one(tf, cd28_catalog$pathway[1], cd28_membership[pathway == cd28_catalog$pathway[1], gene_symbol]),
    ora_one(tf, cd28_catalog$pathway[2], cd28_membership[pathway == cd28_catalog$pathway[2], gene_symbol]),
    ora_one(tf, "CD28_CORE_UNION", cd28_union)
  ))
}))
ora_results[, BH_q := p.adjust(p_value, method = "BH"), by = pathway]
ora_results <- merge(
  ora_results,
  rbind(
    cd28_catalog[, .(pathway, pathway_label)],
    data.table(pathway = "CD28_CORE_UNION", pathway_label = "CD28 core union")
  ),
  by = "pathway", all.x = TRUE
)
setorder(ora_results, pathway, BH_q, -overlap_n, -odds_ratio)
fwrite(ora_results, file.path(table_dir, "05_TF_target_overlap_with_CD28_pathways.csv"))

# TFs directly regulating the CD28 gene in the high-confidence signed network.
direct_cd28 <- regulon[target == "CD28" & confidence %in% c("A", "B")]
direct_cd28 <- merge(direct_cd28, tf_coverage, by = "TF", all.x = TRUE)
fwrite(direct_cd28, file.path(table_dir, "06_direct_CD28_gene_regulators.csv"))

# -------------------------------------------------------------------------
# ULM activity inference, reproducing the documented decoupleR statistic.
# -------------------------------------------------------------------------
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
  i <- match(net$target, shared_targets)
  j <- match(net$TF, tfs)
  mor_mat[cbind(i, j)] <- net$mor
  y <- mat[shared_targets, , drop = FALSE]
  if (center_rows) y <- sweep(y, 1L, rowMeans(y), FUN = "-")
  r <- cor(mor_mat, y)
  r[!is.finite(r)] <- NA_real_
  r <- pmax(pmin(r, 1 - 1e-12), -1 + 1e-12)
  df <- nrow(mor_mat) - 2L
  score <- r * sqrt(df / ((1 - r + 1e-20) * (1 + r + 1e-20)))
  attr(score, "n_shared_targets") <- nrow(mor_mat)
  attr(score, "df") <- df
  score
}

activity_tf_by_sample <- ulm_from_matrix(
  log_expression, regulon_tested, minsize = min_targets, center_rows = TRUE
)
stopifnot(
  identical(colnames(activity_tf_by_sample), meta$sample_id),
  !anyNA(activity_tf_by_sample),
  all(is.finite(activity_tf_by_sample))
)

fit_activity_model <- function(activity_tf_by_sample, meta, add_tme = FALSE) {
  design_formula <- if (add_tme) {
    ~ therapy_short + response_group + age_z + gender + immune_score_z + stromal_score_z + degradation_group
  } else {
    ~ therapy_short + response_group + age_z + gender + degradation_group
  }
  design <- model.matrix(design_formula, data = meta)
  coef_name <- "degradation_groupHigh"
  if (!coef_name %in% colnames(design)) stop("Group coefficient absent from model")
  y <- activity_tf_by_sample
  fit <- eBayes(lmFit(y, design), trend = FALSE, robust = TRUE)
  tt <- as.data.table(topTable(
    fit, coef = coef_name, number = Inf, sort.by = "none", adjust.method = "BH"
  ), keep.rownames = "TF")
  outcome_sd <- apply(y, 1L, sd)
  se <- fit$stdev.unscaled[, coef_name] * sqrt(fit$s2.post)
  out <- tt[, .(
    TF, activity_difference = logFC, average_activity = AveExpr,
    moderated_t = t, p_value = P.Value, BH_q = adj.P.Val
  )]
  out[, standard_error := se[TF]]
  out[, `:=`(
    CI_low = activity_difference - qt(0.975, df = fit$df.total[TF]) * standard_error,
    CI_high = activity_difference + qt(0.975, df = fit$df.total[TF]) * standard_error
  )]
  out[, `:=`(
    standardized_effect = activity_difference / outcome_sd[TF],
    standardized_CI_low = CI_low / outcome_sd[TF],
    standardized_CI_high = CI_high / outcome_sd[TF],
    direction = fifelse(activity_difference > 0, "Higher in Degradation-High", "Higher in Degradation-Low"),
    model = if (add_tme) "Clinical + ESTIMATE sensitivity" else "Clinical-covariate adjusted"
  )]
  list(results = out, design = design, fit = fit)
}

primary_fit <- fit_activity_model(activity_tf_by_sample, meta, add_tme = FALSE)
tme_fit <- fit_activity_model(activity_tf_by_sample, meta, add_tme = TRUE)
primary_activity <- merge(primary_fit$results, tf_coverage, by = "TF", all.x = TRUE)
tme_activity <- merge(tme_fit$results, tf_coverage, by = "TF", all.x = TRUE)
primary_activity[, absolute_standardized_effect := abs(standardized_effect)]
tme_activity[, absolute_standardized_effect := abs(standardized_effect)]
setorder(primary_activity, BH_q, -absolute_standardized_effect)
setorder(tme_activity, BH_q, -absolute_standardized_effect)
fwrite(primary_activity, file.path(table_dir, "07_TF_ULM_activity_High_vs_Low.csv"))
fwrite(tme_activity, file.path(table_dir, "08_TF_ULM_activity_ESTIMATE_sensitivity.csv"))

activity_long <- as.data.table(as.table(activity_tf_by_sample))
setnames(activity_long, c("TF", "sample_id", "ULM_activity"))
activity_long <- merge(
  activity_long,
  meta[, .(sample_id, degradation_group, Degradation_score)],
  by = "sample_id", all.x = TRUE
)
stopifnot(nrow(activity_long) == nrow(activity_tf_by_sample) * ncol(activity_tf_by_sample))
fwrite(activity_long, file.path(table_dir, "09_sample_level_TF_ULM_activities_long.csv"))

# Contrast-level ULM is an independent implementation check using the prior
# covariate-adjusted moderated-t gene ranking rather than sample activities.
ulm_from_vector <- function(values, net, minsize = 10L) {
  stopifnot(!is.null(names(values)), !anyDuplicated(names(values)))
  net <- net[target %in% names(values)]
  keep_tf <- net[, uniqueN(target), by = TF][V1 >= minsize, TF]
  net <- net[TF %in% keep_tf]
  shared_targets <- sort(unique(net$target))
  tfs <- sort(unique(net$TF))
  mor_mat <- matrix(0, nrow = length(shared_targets), ncol = length(tfs),
                    dimnames = list(shared_targets, tfs))
  mor_mat[cbind(match(net$target, shared_targets), match(net$TF, tfs))] <- net$mor
  y <- matrix(values[shared_targets], ncol = 1L, dimnames = list(shared_targets, "High_vs_Low"))
  r <- cor(mor_mat, y)
  df <- nrow(mor_mat) - 2L
  score <- r * sqrt(df / ((1 - r + 1e-20) * (1 + r + 1e-20)))
  data.table(TF = rownames(score), contrast_ULM_score = score[, 1])
}

primary_t <- primary_genes$moderated_t
names(primary_t) <- primary_genes$gene_symbol
contrast_ulm <- ulm_from_vector(primary_t, regulon_tested, minsize = min_targets)
contrast_check <- merge(
  primary_activity[, .(TF, sample_model_moderated_t = moderated_t, sample_model_BH_q = BH_q)],
  contrast_ulm, by = "TF"
)
fwrite(contrast_check, file.path(table_dir, "10_contrast_vs_sample_ULM_validation.csv"))

# -------------------------------------------------------------------------
# Candidate prioritization and reader-facing figures
# -------------------------------------------------------------------------
union_ora <- ora_results[pathway == "CD28_CORE_UNION"]
candidate_evidence <- merge(primary_activity, union_ora, by = "TF", all.x = TRUE)
candidate_evidence[, direct_CD28_regulator := TF %in% direct_cd28$TF]
candidate_evidence[, structural_link := overlap_n >= 2L & p_value.y < 0.05]
candidate_evidence[, activity_FDR_significant := BH_q.x < 0.05]
candidate_evidence[, candidate_class := fifelse(
  direct_CD28_regulator & activity_FDR_significant,
  "Direct CD28 regulator + differential TF activity",
  fifelse(
    structural_link & activity_FDR_significant,
    "CD28-pathway target overlap + differential TF activity",
    fifelse(direct_CD28_regulator, "Direct CD28 regulator", "Other tested TF")
  )
)]

# Rename colliding columns created by the merge.
setnames(candidate_evidence,
         old = c("p_value.x", "BH_q.x", "p_value.y", "BH_q.y"),
         new = c("activity_p", "activity_BH_q", "overlap_p", "overlap_BH_q"))

# Activity-significant TFs with a direct/structural CD28 link are the main
# candidates. Direct CD28 regulators are included even if their activity is not
# significant, so the figure cannot silently omit a biologically direct edge.
main_candidates <- candidate_evidence[
  (activity_BH_q < 0.05 & overlap_n >= 2L & overlap_p < 0.05) |
    direct_CD28_regulator == TRUE
]
main_candidates[, priority_group := fifelse(
  direct_CD28_regulator, 1L,
  fifelse(overlap_BH_q < 0.10, 2L, 3L)
)]
setorder(main_candidates, priority_group, activity_BH_q, overlap_p, -overlap_n)

# If fewer than six candidates pass the predeclared intersection, supplement
# with the most significant differentially active TFs having >=2 pathway
# targets. These rows are explicitly labeled as exploratory in the table.
if (nrow(main_candidates) < 6L) {
  fillers <- candidate_evidence[
    activity_BH_q < 0.05 & overlap_n >= 2L & !TF %in% main_candidates$TF
  ]
  setorder(fillers, activity_BH_q, overlap_p, -overlap_n)
  fillers[, candidate_class := "Exploratory: differential activity + >=2 CD28-pathway targets"]
  main_candidates <- rbindlist(list(main_candidates, head(fillers, 6L - nrow(main_candidates))), fill = TRUE)
}
main_candidates <- unique(main_candidates, by = "TF")
if (nrow(main_candidates) > 8L) main_candidates <- main_candidates[1:8]
if (!nrow(main_candidates)) stop("No TF candidates satisfied the transparent selection rule")

candidate_tfs <- main_candidates$TF
candidate_evidence[, selected_for_figure := TF %in% candidate_tfs]
setorder(candidate_evidence, -selected_for_figure, activity_BH_q, overlap_p)
fwrite(candidate_evidence, file.path(table_dir, "11_CD28_candidate_TF_evidence.csv"))
fwrite(main_candidates, file.path(table_dir, "12_CD28_candidate_TFs_selected_for_figure.csv"))

candidate_sample <- activity_long[TF %in% candidate_tfs]
fwrite(candidate_sample, file.path(table_dir, "13_candidate_TF_sample_activities.csv"))

# Cross-check candidate names against existing MSigDB C3:TFT legacy terms.
# C3 legacy motifs often use family/alias labels rather than HGNC symbols:
# NFKB1/RELA -> NFKB and SPI1 -> PU1. These are explicitly recorded as
# family/alias cross-checks, not treated as one-to-one TF evidence.
c3_alias <- data.table(
  TF = c("EGR1", "SP1", "SPI1", "NFKB1", "RELA"),
  c3_pattern = c("EGR1", "SP1", "PU1", "NFKB", "NFKB"),
  mapping_type = c("exact symbol", "exact symbol", "HGNC alias", "TF-family motif", "TF-family motif")
)
c3_alias <- c3_alias[TF %in% candidate_tfs]
c3_crosscheck <- rbindlist(lapply(seq_len(nrow(c3_alias)), function(i) {
  hit <- copy(c3_tft[grepl(c3_alias$c3_pattern[i], pathway, ignore.case = TRUE)])
  if (nrow(hit)) hit[, `:=`(
    candidate_TF = c3_alias$TF[i],
    motif_query = c3_alias$c3_pattern[i],
    mapping_type = c3_alias$mapping_type[i]
  )]
  hit
}), fill = TRUE)
if (nrow(c3_crosscheck)) {
  setcolorder(c3_crosscheck, c("candidate_TF", "motif_query", "mapping_type", setdiff(names(c3_crosscheck), c("candidate_TF", "motif_query", "mapping_type"))))
}
fwrite(c3_crosscheck, file.path(table_dir, "14_candidate_TF_C3_TFT_crosscheck.csv"))

high_color <- "#C93C3C"
low_color <- "#1F9E93"
neutral_color <- "#8E8E8E"

plot_order <- main_candidates[order(standardized_effect), TF]
forest_dt <- copy(main_candidates)
forest_dt[, TF_plot := factor(TF, levels = plot_order)]
forest_dt[, significant := activity_BH_q < 0.05]
forest_dt[, q_label := paste0("q = ", formatC(activity_BH_q, format = "g", digits = 2))]

p_forest <- ggplot(forest_dt, aes(x = standardized_effect, y = TF_plot)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.5, color = "#666666") +
  geom_errorbar(
    aes(xmin = standardized_CI_low, xmax = standardized_CI_high),
    orientation = "y", width = 0, linewidth = 0.85, color = "#555555"
  ) +
  geom_point(aes(color = significant, size = overlap_n), alpha = 0.95) +
  scale_color_manual(values = c(`TRUE` = high_color, `FALSE` = neutral_color), guide = "none") +
  scale_size_continuous(name = "CD28-pathway\ntargets", range = c(3.2, 8)) +
  annotate("text", x = -Inf, y = Inf, label = "Degradation-Low higher",
           hjust = -0.02, vjust = 1.6, color = low_color, fontface = "bold", size = 4.2) +
  annotate("text", x = Inf, y = Inf, label = "Degradation-High higher",
           hjust = 1.02, vjust = 1.6, color = high_color, fontface = "bold", size = 4.2) +
  labs(
    title = "TF programs linked to CD28-pathway genes",
    subtitle = "Signed DoRothEA A/B target evidence; dots show adjusted High-minus-Low inferred TF-activity effects",
    x = "Standardized TF-activity difference (High - Low)", y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 11, color = "#4D4D4D"),
    axis.text.y = element_text(face = "bold"),
    legend.position = "right"
  )

box_dt <- candidate_sample[TF %in% candidate_tfs]
box_dt[, TF := factor(TF, levels = candidate_tfs)]
box_dt[, degradation_group := factor(degradation_group, levels = c("Low", "High"))]
q_annot <- main_candidates[, .(
  TF, label = paste0("q = ", formatC(activity_BH_q, format = "g", digits = 2))
)]
box_ranges <- box_dt[, .(y = max(ULM_activity) + 0.10 * diff(range(ULM_activity))), by = TF]
q_annot <- merge(q_annot, box_ranges, by = "TF", all.x = TRUE)
q_annot[, TF := factor(TF, levels = candidate_tfs)]

p_box <- ggplot(box_dt, aes(x = degradation_group, y = ULM_activity, color = degradation_group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, linewidth = 0.65, color = "#222222", fill = "white") +
  geom_jitter(width = 0.12, height = 0, size = 1.25, alpha = 0.72) +
  geom_text(
    data = q_annot, aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE, size = 3.2, color = "#333333"
  ) +
  facet_wrap(~ TF, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c(Low = low_color, High = high_color), guide = "none") +
  labs(
    title = "Tumor-level inferred TF regulon activities",
    subtitle = "Low n = 37; High n = 36. Points represent PRE-treatment tumors.",
    x = "Sarcosine-degradation group", y = "ULM TF-activity score"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 10, color = "#4D4D4D"),
    strip.text = element_text(face = "bold", size = 11),
    axis.text.x = element_text(face = "bold")
  )

combined <- p_forest / p_box +
  plot_annotation(
    title = "Candidate transcription-factor programs linked to CD28-pathway expression",
    subtitle = paste0(
      "TIGER melanoma PRE-only; Sarcosine Degradation High vs Low. ",
      "Group effects adjusted for therapy, clinical response, age and sex."
    ),
    tag_levels = "a",
    theme = theme(
      plot.title = element_text(face = "bold", size = 22),
      plot.subtitle = element_text(size = 12, color = "#4D4D4D")
    )
  ) + plot_layout(heights = c(0.9, 1.45))

ggsave(
  file.path(figure_dir, "Fig_CD28_Candidate_TF_Activity_High_vs_Low.png"),
  combined, width = 14.5, height = 11.5, dpi = 400, bg = "white"
)
ggsave(
  file.path(figure_dir, "Fig_CD28_Candidate_TF_Activity_High_vs_Low.pdf"),
  combined, width = 14.5, height = 11.5, device = cairo_pdf, bg = "white"
)

# Sample activity heatmap with only High/Low group annotation.
heat_mat <- activity_tf_by_sample[candidate_tfs, , drop = FALSE]
heat_mat <- t(scale(t(heat_mat)))
heat_mat[heat_mat > 2.5] <- 2.5
heat_mat[heat_mat < -2.5] <- -2.5
sample_order <- meta[order(degradation_group, Degradation_score), sample_id]
heat_mat <- heat_mat[, sample_order, drop = FALSE]
annotation_col <- data.frame(
  `Sarcosine degradation` = meta[match(sample_order, sample_id), as.character(degradation_group)],
  row.names = sample_order,
  check.names = FALSE
)
annotation_colors <- list(`Sarcosine degradation` = c(Low = low_color, High = high_color))
heat_png <- file.path(figure_dir, "Fig_CD28_Candidate_TF_Activity_Heatmap_High_vs_Low.png")
png(heat_png, width = 3600, height = 1900, res = 320, bg = "white")
pheatmap(
  heat_mat,
  cluster_rows = TRUE, cluster_cols = FALSE,
  show_colnames = FALSE, show_rownames = TRUE,
  annotation_col = annotation_col, annotation_colors = annotation_colors,
  annotation_legend = TRUE, border_color = NA,
  gaps_col = sum(meta$degradation_group == "Low"),
  color = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(101),
  breaks = seq(-2.5, 2.5, length.out = 102),
  fontsize = 13, fontsize_row = 13,
  main = "Candidate CD28-regulatory TF activities across PRE-treatment tumors"
)
dev.off()

# TF-to-pathway target map: structural resource evidence, not expression effect.
target_edges <- regulon_tested[
  TF %in% candidate_tfs & target %in% cd28_union_expressed,
  .(TF, target, mor, confidence)
]
fwrite(target_edges, file.path(table_dir, "15_candidate_TF_to_CD28_pathway_edges.csv"))
if (nrow(target_edges)) {
  target_order <- target_edges[, .N, by = target][order(-N, target), target]
  edge_plot <- ggplot(target_edges, aes(x = factor(target, levels = target_order), y = factor(TF, levels = rev(candidate_tfs)))) +
    geom_point(aes(fill = factor(mor), size = confidence), shape = 21, color = "#222222", stroke = 0.35) +
    scale_fill_manual(
      values = c(`-1` = "#2166AC", `1` = "#B2182B"),
      labels = c(`-1` = "Represses target", `1` = "Activates target"),
      name = "Regulatory sign"
    ) +
    scale_size_manual(values = c(A = 5.5, B = 4), name = "DoRothEA confidence") +
    labs(
      title = "Curated TF-to-CD28-pathway regulatory links",
      subtitle = "Full Reactome CD28 core membership; absence of a dot means no signed A/B DoRothEA edge",
      x = "CD28-pathway gene", y = "Candidate TF"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", size = 18),
      axis.text.x = element_text(angle = 50, hjust = 1, vjust = 1, size = 9),
      axis.text.y = element_text(face = "bold")
    )
  ggsave(
    file.path(figure_dir, "Fig_CD28_Candidate_TF_to_Pathway_Targets.png"),
    edge_plot, width = max(10, 0.42 * length(target_order) + 5), height = 6.8,
    dpi = 400, limitsize = FALSE, bg = "white"
  )
}

# -------------------------------------------------------------------------
# Reproducibility and validation records
# -------------------------------------------------------------------------
input_manifest <- data.table(
  role = c(
    "FPKM expression", "PRE73 metadata", "primary gene model",
    "primary Reactome GSEA", "TME gene sensitivity", "ESTIMATE scores",
    "existing C3:TFT cross-check", "DoRothEA TF-target resource"
  ),
  path = normalizePath(required_inputs),
  sha256 = vapply(required_inputs, sha256_file, character(1)),
  size_bytes = file.info(required_inputs)$size
)
fwrite(input_manifest, file.path(table_dir, "00_input_manifest.csv"))

method_contract <- data.table(
  field = c(
    "analysis_unit", "population", "group_definition", "primary_estimand",
    "primary_covariates", "sensitivity_covariates", "expression_space",
    "expression_filter", "regulon_resource", "regulon_confidence",
    "activity_method", "activity_input", "minimum_targets", "multiple_testing",
    "structural_overlap", "figure_grouping", "interpretation_limit"
  ),
  value = c(
    "one PRE-treatment tumor per patient",
    "TIGER melanoma PRE-only n=73",
    "median split of z-scored SARDH/PIPOX degradation module: High n=36, Low n=37",
    "difference in signed TF regulon activity: High minus Low",
    "therapy + clinical response + age + sex",
    "primary covariates + ESTIMATE immune and stromal scores",
    "log2(FPKM+1), gene-centered for sample-level ULM",
    "FPKM >=1 in at least ceiling(10% of 73) samples; non-zero variance",
    paste0("OmniPath DoRothEA human query downloaded 2026-09-02; MSigDB ", msigdb_version, " Reactome"),
    "signed A/B only; unknown/contradictory directions excluded",
    "ULM t-statistic matching the documented decoupleR correlation-to-t calculation",
    "sample-level centered expression; contrast-level moderated-t used only for validation",
    as.character(min_targets),
    "Benjamini-Hochberg across all tested TFs; separately within each structural-overlap pathway",
    "one-sided Fisher exact test against all expression-filtered genes; full pathway membership",
    "Sarcosine Degradation Low versus High only; statistical model names omitted from reader-facing figures",
    "bulk-RNA association; not causal or T-cell intrinsic"
  )
)
fwrite(method_contract, file.path(table_dir, "00_method_contract.csv"))

pkg_versions <- data.table(
  package = required_packages,
  version = vapply(required_packages, function(x) as.character(packageVersion(x)), character(1))
)
pkg_versions <- rbind(
  data.table(package = "R", version = R.version.string),
  pkg_versions
)
fwrite(pkg_versions, file.path(table_dir, "00_package_versions.csv"))

validation <- data.table(
  check = c(
    "PRE sample count", "unique sample IDs", "group sizes 37/36",
    "expression genes match prior primary table", "signed regulon unique TF-target edges",
    "tested TFs have >=10 targets", "ULM values finite",
    "sample-model versus contrast-ULM direction agreement",
    "direct CD28 regulators retained in candidate audit",
    "candidate figure contains no model-group labels",
    "primary design full rank", "TME sensitivity design full rank"
  ),
  value = c(
    nrow(meta), uniqueN(meta$sample_id),
    paste(sum(meta$degradation_group == "Low"), sum(meta$degradation_group == "High"), sep = "/"),
    setequal(rownames(log_expression), primary_genes$gene_symbol),
    !anyDuplicated(regulon_tested[, .(TF, target)]),
    min(tf_coverage[TF %in% eligible_tfs, expressed_targets_n]),
    all(is.finite(activity_tf_by_sample)),
    mean(sign(contrast_check$sample_model_moderated_t) == sign(contrast_check$contrast_ULM_score)),
    all(direct_cd28$TF %in% candidate_evidence[direct_CD28_regulator == TRUE, TF]),
    !any(grepl("Primary clinical model|ESTIMATE immune/stromal", c(
      p_forest$labels$title, p_forest$labels$subtitle,
      p_box$labels$title, p_box$labels$subtitle
    ), ignore.case = TRUE)),
    qr(primary_fit$design)$rank == ncol(primary_fit$design),
    qr(tme_fit$design)$rank == ncol(tme_fit$design)
  ),
  expected = c(
    "73", "73", "37/36", "TRUE", "TRUE", ">=10", "TRUE", ">=0.90",
    "TRUE", "TRUE", "TRUE", "TRUE"
  )
)
validation[, pass := c(
  value[1] == 73, value[2] == 73, value[3] == "37/36",
  value[4] == TRUE, value[5] == TRUE, as.numeric(value[6]) >= 10,
  value[7] == TRUE, as.numeric(value[8]) >= 0.90,
  value[9] == TRUE, value[10] == TRUE, value[11] == TRUE, value[12] == TRUE
)]
fwrite(validation, file.path(table_dir, "16_validation_checks.csv"))
if (!all(validation$pass)) {
  print(validation)
  stop("One or more validation checks failed")
}

summary_lines <- c(
  "TIGER PRE73 CD28 upstream-TF analysis completed.",
  paste0("Expression-filtered genes: ", nrow(log_expression)),
  paste0("Signed DoRothEA A/B edges tested: ", nrow(regulon_tested)),
  paste0("TFs tested (>=", min_targets, " expressed targets): ", length(eligible_tfs)),
  paste0("Figure candidates: ", paste(candidate_tfs, collapse = ", ")),
  paste0("Activity-FDR-significant tested TFs: ", sum(primary_activity$BH_q < 0.05)),
  paste0("Contrast/sample ULM direction agreement: ",
         sprintf("%.3f", mean(sign(contrast_check$sample_model_moderated_t) == sign(contrast_check$contrast_ULM_score)))),
  "Reader-facing figures show only Sarcosine Degradation High and Low groups.",
  "ESTIMATE adjustment is preserved only in table 08 as a sensitivity analysis.",
  "Interpretation: candidate TF-regulated programs, not causal upstream regulators."
)
writeLines(summary_lines, file.path(log_dir, "analysis_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
