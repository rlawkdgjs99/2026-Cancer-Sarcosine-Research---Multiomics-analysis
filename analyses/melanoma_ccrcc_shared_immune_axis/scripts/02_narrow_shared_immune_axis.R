#!/usr/bin/env Rscript

# Cross-cohort exploratory narrowing of immune mechanisms shared by:
#   1) TIGER melanoma sarcosine-degradation High tumors, and
#   2) matched ccRCC measured-tumor-Sarcosine Low tumors.
#
# The analysis starts from a broad, biologically motivated candidate catalog,
# applies exact-pathway replication gates, then compares gene-disjoint sample-
# level route scores. It ranks observational routes; it does not establish
# temporal causality, same-cell signaling, metabolite flux, or protein activity.

options(stringsAsFactors = FALSE, width = 200)
set.seed(260902)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
module_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
input_dir <- Sys.getenv("SARCO_SHARED_AXIS_INPUT_DIR", unset = file.path(module_dir, "inputs"))

gsea_lib <- Sys.getenv("SARCO_GSEA_R_LIB", unset = "")
if (nzchar(gsea_lib)) {
  if (!dir.exists(gsea_lib)) stop("SARCO_GSEA_R_LIB does not exist: ", gsea_lib)
  .libPaths(unique(c(normalizePath(gsea_lib), .libPaths())))
}

required_packages <- c("data.table", "msigdbr", "ggplot2", "patchwork", "pheatmap")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
  library(pheatmap)
})

find_unique_dir <- function(root, basename_target) {
  hits <- list.dirs(root, recursive = TRUE, full.names = TRUE)
  hits <- hits[basename(hits) == basename_target]
  if (length(hits) != 1L) stop("Expected one directory named ", basename_target, "; found ", length(hits))
  normalizePath(hits[[1]])
}

find_unique_file <- function(root, filename) {
  hits <- list.files(root, pattern = paste0("^", filename, "$"), recursive = TRUE, full.names = TRUE)
  if (length(hits) != 1L) stop("Expected one file named ", filename, "; found ", length(hits))
  normalizePath(hits[[1]])
}

sha256_file <- function(path) {
  out <- system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE)
  strsplit(out, "[[:space:]]+")[[1]][1]
}

zscore <- function(x) {
  ans <- as.numeric(scale(as.numeric(x)))
  if (anyNA(ans) || any(!is.finite(ans))) stop("Non-finite z score")
  ans
}

score_gene_set <- function(genes, label, expression_matrix, minimum_genes = 10L) {
  mapped <- sort(intersect(unique(genes), rownames(expression_matrix)))
  mapped <- mapped[apply(expression_matrix[mapped, , drop = FALSE], 1L, var) > 0]
  if (length(mapped) < minimum_genes) {
    stop(label, " has only ", length(mapped), " usable genes; minimum is ", minimum_genes)
  }
  gene_z <- t(scale(t(expression_matrix[mapped, , drop = FALSE])))
  if (anyNA(gene_z) || any(!is.finite(gene_z))) stop("Non-finite row z score in ", label)
  score <- zscore(colMeans(gene_z))
  names(score) <- colnames(expression_matrix)
  list(score = score, genes = mapped)
}

ulm_from_matrix <- function(mat, net, minsize = 10L) {
  stopifnot(is.matrix(mat), all(c("TF", "target", "mor") %in% names(net)))
  net <- unique(net[, .(TF, target, mor)])
  shared_targets <- sort(intersect(rownames(mat), unique(net$target)))
  net <- net[target %in% shared_targets]
  keep_tf <- net[, uniqueN(target), by = TF][V1 >= minsize, TF]
  net <- net[TF %in% keep_tf]
  shared_targets <- sort(unique(net$target))
  tfs <- sort(unique(net$TF))
  mor_mat <- matrix(0, nrow = length(shared_targets), ncol = length(tfs),
                    dimnames = list(shared_targets, tfs))
  mor_mat[cbind(match(net$target, shared_targets), match(net$TF, tfs))] <- net$mor
  y <- sweep(mat[shared_targets, , drop = FALSE], 1L,
             rowMeans(mat[shared_targets, , drop = FALSE]), FUN = "-")
  r <- cor(mor_mat, y)
  r <- pmax(pmin(r, 1 - 1e-12), -1 + 1e-12)
  df <- nrow(mor_mat) - 2L
  score <- r * sqrt(df / ((1 - r + 1e-20) * (1 + r + 1e-20)))
  if (anyNA(score) || any(!is.finite(score))) stop("Non-finite ULM score")
  score
}

extract_term <- function(fit, term, cohort, model, outcome, predictor, route = NA_character_, edge = NA_character_) {
  sm <- summary(fit)$coefficients
  if (!term %in% rownames(sm)) stop("Term absent: ", term, " in ", cohort, " / ", model, " / ", outcome)
  estimate <- unname(sm[term, "Estimate"])
  se <- unname(sm[term, "Std. Error"])
  df <- df.residual(fit)
  crit <- qt(0.975, df)
  data.table(
    cohort = cohort, model = model, route = route, edge = edge,
    outcome = outcome, predictor = predictor,
    standardized_beta = estimate, standard_error = se,
    CI_low = estimate - crit * se, CI_high = estimate + crit * se,
    t_value = unname(sm[term, "t value"]),
    p_value = unname(sm[term, "Pr(>|t|)"]), residual_df = df,
    partial_R2 = unname(sm[term, "t value"])^2 /
      (unname(sm[term, "t value"])^2 + df)
  )
}

fit_term <- function(data, outcome, predictor, exposure, covars, extras, cohort, model,
                     route = NA_character_, edge = NA_character_) {
  data <- droplevels(copy(data))
  rhs <- unique(c(exposure, predictor, covars, extras))
  form <- as.formula(paste(outcome, "~", paste(rhs, collapse = " + ")))
  x <- model.matrix(form, data = data)
  if (qr(x)$rank != ncol(x)) stop("Rank-deficient model: ", cohort, " / ", model, " / ", outcome)
  fit <- lm(form, data = data)
  extract_term(fit, predictor, cohort, model, outcome, predictor, route, edge)
}

serial_design <- function(data, nodes, exposure, covars, extras) {
  data <- droplevels(copy(data))
  xs <- vector("list", length(nodes)); ys <- vector("list", length(nodes)); terms <- character(length(nodes))
  for (i in seq_along(nodes)) {
    if (i == 1L) {
      rhs <- unique(c(exposure, covars, extras)); term <- exposure
    } else {
      rhs <- unique(c(exposure, nodes[seq_len(i - 1L)], covars, extras)); term <- nodes[i - 1L]
    }
    form <- as.formula(paste("~", paste(rhs, collapse = " + ")))
    x <- model.matrix(form, data = data)
    if (qr(x)$rank != ncol(x)) stop("Rank-deficient serial design at ", nodes[i])
    xs[[i]] <- x; ys[[i]] <- data[[nodes[i]]]; terms[i] <- term
  }
  list(x = xs, y = ys, term = terms, nodes = nodes)
}

serial_edges <- function(data, nodes, exposure, covars, extras, cohort, model, route, direction) {
  obj <- serial_design(data, nodes, exposure, covars, extras)
  rbindlist(lapply(seq_along(nodes), function(i) {
    fit <- lm(as.formula(paste(nodes[i], "~", paste(unique(c(
      exposure,
      if (i > 1L) nodes[seq_len(i - 1L)] else character(),
      covars, extras
    )), collapse = " + "))), data = data)
    edge_label <- if (i == 1L) paste("Target state ->", nodes[i]) else paste(nodes[i - 1L], "->", nodes[i])
    extract_term(fit, obj$term[i], cohort, model, nodes[i], obj$term[i], route,
                 paste(direction, edge_label))[, `:=`(direction = direction, edge_order = i)]
  }))
}

serial_bootstrap <- function(data, nodes, exposure, covars, extras, cohort, model, route,
                             direction, B = 5000L, seed = 260902L) {
  obj <- serial_design(data, nodes, exposure, covars, extras)
  indices <- mapply(function(x, term) match(term, colnames(x)), obj$x, obj$term)
  observed_edges <- mapply(function(x, y, idx) coef(lm.fit(x, y))[idx], obj$x, obj$y, indices)
  values <- rep(NA_real_, B)
  set.seed(seed)
  for (b in seq_len(B)) {
    ids <- sample.int(nrow(data), nrow(data), replace = TRUE)
    cf <- mapply(function(x, y, idx) {
      tryCatch(coef(lm.fit(x[ids, , drop = FALSE], y[ids]))[idx], error = function(e) NA_real_)
    }, obj$x, obj$y, indices)
    values[b] <- prod(cf)
  }
  values <- values[is.finite(values)]
  if (length(values) < 0.95 * B) stop("Too many invalid bootstrap replicates: ", cohort, " / ", route)
  ci <- quantile(values, c(0.025, 0.975), names = FALSE)
  data.table(
    cohort = cohort, model = model, route = route, direction = direction,
    nodes = paste(nodes, collapse = " -> "),
    serial_product = prod(observed_edges), bootstrap_CI_low = ci[1], bootstrap_CI_high = ci[2],
    positive_fraction = mean(values > 0), negative_fraction = mean(values < 0),
    CI_excludes_zero = ci[1] > 0 | ci[2] < 0,
    bootstrap_replicates_requested = B, bootstrap_replicates_valid = length(values)
  )
}

# -------------------------------------------------------------------------
# Input contracts
# -------------------------------------------------------------------------
input_path <- function(env_name, default_name) {
  value <- Sys.getenv(env_name, unset = file.path(input_dir, default_name))
  if (!file.exists(value)) stop("Missing input ", env_name, ": ", value)
  normalizePath(value, mustWork = TRUE)
}

cc_expression_path <- input_path("SARCO_CCRCC_EXPRESSION", "ccrcc_expression_TPM.csv")
cc_metadata_path <- input_path("SARCO_CCRCC_METADATA", "ccrcc_sarcosine_group_map.csv")
cc_deconv_path <- input_path("SARCO_CCRCC_DECONVOLUTION", "ccrcc_deconvolution.csv")
cc_regulon_path <- input_path("SARCO_CCRCC_REGULON", "ccrcc_signed_regulon.csv")
common_pathway_path <- input_path("SARCO_COMMON_PATHWAY_JOIN", "common_exact_pathway_join.csv")
tiger_expression_path <- input_path("SARCO_TIGER_EXPRESSION", "tiger_expression_FPKM.csv")
tiger_metadata_path <- input_path("SARCO_TIGER_METADATA", "tiger_PRE73_metadata.csv")
tiger_regulon_path <- input_path("SARCO_TIGER_REGULON", "tiger_signed_regulon.csv")
tiger_context_path <- input_path("SARCO_TIGER_CELL_CONTEXT", "tiger_cell_context.csv")

input_paths <- c(
  cc_expression_path, cc_metadata_path, cc_deconv_path, cc_regulon_path,
  common_pathway_path, tiger_expression_path, tiger_metadata_path,
  tiger_regulon_path, tiger_context_path
)
for (path in input_paths) if (!file.exists(path)) stop("Missing input: ", path)

output_root <- file.path(module_dir, "results", "mechanism_narrowing")
table_dir <- output_root
figure_dir <- file.path(output_root, "figures")
reader_figure_dir <- file.path(output_root, "figures_reader")
log_dir <- file.path(output_root, "logs")
for (d in c(table_dir, figure_dir, reader_figure_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Red/teal is reserved for target/reference states within a cohort. Cross-cohort
# comparisons use purple/blue so that cohort identity cannot be confused with
# biological direction.
cohort_palette <- c(
  "Melanoma: Degradation-High" = "#7B2CBF",
  "ccRCC: Sarcosine-Low" = "#1976D2"
)

# -------------------------------------------------------------------------
# Expression and metadata
# -------------------------------------------------------------------------
tiger_meta <- fread(tiger_metadata_path, check.names = FALSE)
stopifnot(nrow(tiger_meta) == 73L, uniqueN(tiger_meta$sample_id) == 73L)
tiger_ids <- tiger_meta$sample_id
tiger_expr_dt <- fread(tiger_expression_path, check.names = FALSE)
tiger_genes <- trimws(tiger_expr_dt[[1]])
if (anyNA(tiger_genes) || any(!nzchar(tiger_genes)) || anyDuplicated(tiger_genes)) stop("Invalid TIGER genes")
if (!all(tiger_ids %in% names(tiger_expr_dt))) stop("TIGER expression lacks PRE73 samples")
tiger_fpkm <- as.matrix(tiger_expr_dt[, ..tiger_ids]); storage.mode(tiger_fpkm) <- "double"
rownames(tiger_fpkm) <- tiger_genes; colnames(tiger_fpkm) <- tiger_ids
if (anyNA(tiger_fpkm) || any(!is.finite(tiger_fpkm)) || any(tiger_fpkm < 0)) stop("Invalid TIGER FPKM")
tiger_log_all <- log2(tiger_fpkm + 1)
tiger_keep <- rowSums(tiger_fpkm >= 1) >= ceiling(0.10 * ncol(tiger_fpkm)) &
  apply(tiger_log_all, 1L, var) > 0
tiger_log <- tiger_log_all[tiger_keep, , drop = FALSE]

tiger_meta[, `:=`(
  target_group = factor(ifelse(Degradation_score > median(Degradation_score), "Target", "Reference"),
                        levels = c("Reference", "Target")),
  target_group_z = zscore(as.integer(Degradation_score > median(Degradation_score))),
  target_continuous_z = zscore(Degradation_score),
  therapy_short = relevel(factor(therapy_short), ref = "antiPD1"),
  response_group = relevel(factor(response_group), ref = "NR"),
  gender = relevel(factor(gender), ref = "Male"),
  age_z = zscore(age)
)]
stopifnot(sum(tiger_meta$target_group == "Target") == 36L, sum(tiger_meta$target_group == "Reference") == 37L)
tiger_context <- fread(tiger_context_path, check.names = FALSE)[,
  .(sample_id, T_cell_context_z, APC_myeloid_context_z)]
tiger_meta <- merge(tiger_meta, tiger_context, by = "sample_id", all.x = TRUE, sort = FALSE)
tiger_meta <- tiger_meta[match(tiger_ids, sample_id)]
if (anyNA(tiger_meta[, .(T_cell_context_z, APC_myeloid_context_z)])) stop("Missing TIGER context")

cc_meta <- fread(cc_metadata_path, check.names = FALSE)
stopifnot(nrow(cc_meta) == 100L, uniqueN(cc_meta$sample_id) == 100L)
cc_ids <- cc_meta$sample_id
cc_expr_dt <- fread(cc_expression_path, check.names = FALSE)
cc_genes <- trimws(cc_expr_dt[[1]])
if (anyNA(cc_genes) || any(!nzchar(cc_genes)) || anyDuplicated(cc_genes)) stop("Invalid ccRCC genes")
if (!all(cc_ids %in% names(cc_expr_dt))) stop("ccRCC expression lacks matched samples")
cc_log_all <- as.matrix(cc_expr_dt[, ..cc_ids]); storage.mode(cc_log_all) <- "double"
rownames(cc_log_all) <- cc_genes; colnames(cc_log_all) <- cc_ids
if (anyNA(cc_log_all) || any(!is.finite(cc_log_all)) || any(cc_log_all < 0)) stop("Invalid ccRCC expression")
cc_keep <- rowSums(cc_log_all >= 1) >= ceiling(0.10 * ncol(cc_log_all)) & apply(cc_log_all, 1L, var) > 0
cc_log <- cc_log_all[cc_keep, , drop = FALSE]
if (nrow(cc_log) != 15119L) stop("Unexpected ccRCC expression universe: ", nrow(cc_log))

cc_meta[, `:=`(
  target_group = factor(ifelse(sarcosine_group == "Low", "Target", "Reference"),
                        levels = c("Reference", "Target")),
  target_group_z = zscore(as.integer(sarcosine_group == "Low")),
  target_continuous_z = zscore(-log2_sarcosine_intensity),
  batch = factor(batch),
  sex = relevel(factor(sex), ref = "female"),
  age = factor(age, levels = c("40-60", "<40", ">60")),
  grade = factor(grade, levels = c("1", "2", "3", "4"))
)]
stopifnot(sum(cc_meta$target_group == "Target") == 50L, sum(cc_meta$target_group == "Reference") == 50L)

cc_deconv <- fread(cc_deconv_path, check.names = FALSE)
t_features <- cc_deconv[
  (cell_type %in% c("CD4 T cells", "CD8 T cells") | (method == "MCPcounter" & cell_type == "Pan T cells")) &
    method != "ESTIMATE"
]
apc_features <- cc_deconv[
  cell_type %in% c("Dendritic cells", "Monocytes", "Macrophages", "Macrophages M1", "Macrophages M2") &
    method != "ESTIMATE"
]
if (uniqueN(t_features$feature_id) < 8L || uniqueN(apc_features$feature_id) < 8L) stop("Insufficient ccRCC named cell context")
cc_t <- t_features[, .(T_cell_context = mean(analysis_value)), by = sample_id]
cc_t[, T_cell_context_z := zscore(T_cell_context)]
cc_t[, T_cell_context := NULL]
cc_apc <- apc_features[, .(APC_myeloid_context = mean(analysis_value)), by = sample_id]
cc_apc[, APC_myeloid_context_z := zscore(APC_myeloid_context)]
cc_apc[, APC_myeloid_context := NULL]
cc_context <- merge(cc_t, cc_apc, by = "sample_id", all = TRUE)
cc_meta <- merge(cc_meta, cc_context, by = "sample_id", all.x = TRUE, sort = FALSE)
cc_meta <- cc_meta[match(cc_ids, sample_id)]
if (anyNA(cc_meta[, .(T_cell_context_z, APC_myeloid_context_z)])) stop("Missing ccRCC context")

# -------------------------------------------------------------------------
# Broad candidate screen and common gene-disjoint modules
# -------------------------------------------------------------------------
candidate_catalog <- data.table(
  candidate_order = 1:11,
  candidate = c(
    "APC cross-presentation", "TCR signaling", "CD28-family regulation (broad)",
    "CD28 co-stimulation (specific)", "Canonical NF-kappaB", "Noncanonical NF-kappaB",
    "CD28-PI3K-AKT", "MAPK/AP-1", "Calcineurin-NFAT", "IL-12/STAT4", "IFNG response"
  ),
  pathway = c(
    "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION", "REACTOME_TCR_SIGNALING",
    "REACTOME_REGULATION_OF_T_CELL_ACTIVATION_BY_CD28_FAMILY", "REACTOME_CO_STIMULATION_BY_CD28",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY",
    "REACTOME_CD28_DEPENDENT_PI3K_AKT_SIGNALING", "REACTOME_MAPK_FAMILY_SIGNALING_CASCADES",
    "GOBP_CALCINEURIN_MEDIATED_SIGNALING", "REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE"
  ),
  candidate_class = c(
    "Upstream immune context", "T-cell receptor", "Co-stimulatory/checkpoint mixture",
    "Specific co-stimulation", "Downstream branch", "Downstream branch",
    "Downstream branch", "Downstream branch", "Downstream branch", "Parallel cytokine branch", "Outcome"
  )
)

common_pathways <- fread(common_pathway_path, check.names = FALSE)
screen <- merge(candidate_catalog, common_pathways, by = "pathway", all.x = TRUE, sort = FALSE)
if (nrow(screen) != nrow(candidate_catalog) || anyNA(screen$tiger_q) || anyNA(screen$ccrcc_q)) {
  stop("Candidate pathway screen did not map completely")
}
screen <- screen[order(candidate_order)]
screen[, exact_replication_pass := target_direction_aligned & tiger_q < 0.05 & ccrcc_q < 0.05]
screen[, screen_decision := fifelse(
  exact_replication_pass & candidate == "CD28-family regulation (broad)",
  "Retain only as a mixed regulatory sensitivity node",
  fifelse(exact_replication_pass, "Advance", "Do not advance as a shared core node")
)]
if (screen[candidate == "CD28 co-stimulation (specific)", exact_replication_pass]) {
  stop("CD28-specific gate unexpectedly passed; revise prespecified route logic")
}

msig <- as.data.table(msigdbr(species = "Homo sapiens", db_species = "HS"))
msigdb_version <- unique(msig$db_version)
if (length(msigdb_version) != 1L) stop("Non-unique MSigDB version")
get_genes <- function(set_name) unique(msig[gs_name == set_name, gene_symbol])

response_ids <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "GOBP_RESPONSE_TO_TYPE_II_INTERFERON"
)
if (!all(response_ids %in% msig$gs_name)) stop("Missing IFNG response source set")
ifng_core <- msig[gs_name %in% response_ids, .N, by = gene_symbol][N >= 2L, gene_symbol]
ifng_core <- setdiff(ifng_core, "IFNG")

predictor_set_ids <- c(
  APC = "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
  TCR = "REACTOME_TCR_SIGNALING",
  CD28_family = "REACTOME_REGULATION_OF_T_CELL_ACTIVATION_BY_CD28_FAMILY",
  Canonical_NFkB = "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  Noncanonical_NFkB = "REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY"
)
raw_predictors <- lapply(predictor_set_ids, get_genes)
without_outcome <- lapply(raw_predictors, setdiff, y = ifng_core)
membership_counts <- table(unlist(without_outcome, use.names = FALSE))
shared_predictor_genes <- names(membership_counts[membership_counts > 1L])
strict_predictors <- lapply(without_outcome, setdiff, y = shared_predictor_genes)
if (any(lengths(strict_predictors) < 20L)) stop("A strict predictor module has fewer than 20 genes")

tiger_score_objects <- lapply(names(strict_predictors), function(nm) {
  score_gene_set(strict_predictors[[nm]], paste("TIGER", nm), tiger_log, 15L)
})
names(tiger_score_objects) <- names(strict_predictors)
cc_score_objects <- lapply(names(strict_predictors), function(nm) {
  score_gene_set(strict_predictors[[nm]], paste("ccRCC", nm), cc_log, 15L)
})
names(cc_score_objects) <- names(strict_predictors)
tiger_ifng <- score_gene_set(ifng_core, "TIGER strict IFNG response", tiger_log, 30L)
cc_ifng <- score_gene_set(ifng_core, "ccRCC strict IFNG response", cc_log, 30L)

score_gene_lists <- c(strict_predictors, list(IFNG_response = ifng_core))
overlap_audit <- rbindlist(lapply(names(score_gene_lists), function(a) {
  rbindlist(lapply(names(score_gene_lists), function(b) {
    data.table(score_1 = a, score_2 = b,
               overlap_n = length(intersect(score_gene_lists[[a]], score_gene_lists[[b]])))
  }))
}))
if (overlap_audit[score_1 != score_2, max(overlap_n)] != 0L) stop("Strict score overlap remains")

score_manifest <- rbindlist(lapply(names(score_gene_lists), function(nm) {
  genes <- score_gene_lists[[nm]]
  data.table(
    score = nm,
    source_pathway = if (nm == "IFNG_response") paste(response_ids, collapse = ";") else predictor_set_ids[[nm]],
    original_genes_n = if (nm == "IFNG_response") length(ifng_core) else length(raw_predictors[[nm]]),
    removed_IFNG_overlap_n = if (nm == "IFNG_response") 0L else length(intersect(raw_predictors[[nm]], ifng_core)),
    removed_predictor_overlap_n = if (nm == "IFNG_response") 0L else length(intersect(without_outcome[[nm]], shared_predictor_genes)),
    strict_gene_universe_n = length(genes),
    TIGER_mapped_n = if (nm == "IFNG_response") length(tiger_ifng$genes) else length(tiger_score_objects[[nm]]$genes),
    ccRCC_mapped_n = if (nm == "IFNG_response") length(cc_ifng$genes) else length(cc_score_objects[[nm]]$genes),
    strict_genes = paste(sort(genes), collapse = ";"), MSigDB_version = msigdb_version
  )
}))

add_scores <- function(meta, score_objects, ifng_object) {
  for (nm in names(score_objects)) meta[, (nm) := score_objects[[nm]]$score[sample_id]]
  meta[, IFNG_response := ifng_object$score[sample_id]]
  if (anyNA(meta[, c(names(score_objects), "IFNG_response"), with = FALSE])) stop("Missing score")
  meta
}
tiger_data <- add_scores(tiger_meta, tiger_score_objects, tiger_ifng)
cc_data <- add_scores(cc_meta, cc_score_objects, cc_ifng)

# -------------------------------------------------------------------------
# Cross-cohort-consensus TF activities, independent of route-score genes
# -------------------------------------------------------------------------
tiger_reg <- fread(tiger_regulon_path, check.names = FALSE)
cc_reg <- fread(cc_regulon_path, check.names = FALSE)
focal_tfs_attempted <- c("SPI1", "NFKB1", "RELA", "NFKB2", "RELB")
required_inferable_tfs <- c("SPI1", "NFKB1", "RELA", "NFKB2")
tiger_reg <- tiger_reg[TF %in% focal_tfs_attempted & confidence %in% c("A", "B") & expressed == TRUE,
                       .(TF, target, tiger_mor = mor)]
cc_reg <- cc_reg[TF %in% focal_tfs_attempted & confidence %in% c("A", "B") & expressed == TRUE,
                 .(TF, target, cc_mor = mor)]
consensus_reg <- merge(tiger_reg, cc_reg, by = c("TF", "target"), allow.cartesian = FALSE)
consensus_reg <- consensus_reg[sign(tiger_mor) == sign(cc_mor)]
excluded_tf_targets <- unique(c(unlist(score_gene_lists, use.names = FALSE), "IFNG"))
consensus_reg <- consensus_reg[!target %in% excluded_tf_targets,
                               .(TF, target, mor = sign(tiger_mor))]
tf_coverage_observed <- consensus_reg[, .(retained_targets_n = uniqueN(target),
                                          targets = paste(sort(unique(target)), collapse = ";")), by = TF]
tf_coverage <- merge(data.table(TF = focal_tfs_attempted), tf_coverage_observed, by = "TF", all.x = TRUE)
tf_coverage[is.na(retained_targets_n), `:=`(retained_targets_n = 0L, targets = "")]
tf_coverage[, estimable_at_minsize_10 := retained_targets_n >= 10L]
if (!all(required_inferable_tfs %in% tf_coverage[estimable_at_minsize_10 == TRUE, TF])) {
  stop("Insufficient consensus TF coverage for a required TF")
}
inferable_tfs <- tf_coverage[estimable_at_minsize_10 == TRUE, TF]
consensus_reg_inferable <- consensus_reg[TF %in% inferable_tfs]
tiger_tf <- ulm_from_matrix(tiger_log, consensus_reg_inferable, 10L)[inferable_tfs, tiger_ids, drop = FALSE]
cc_tf <- ulm_from_matrix(cc_log, consensus_reg_inferable, 10L)[inferable_tfs, cc_ids, drop = FALSE]
tiger_tf <- t(apply(tiger_tf, 1L, zscore)); cc_tf <- t(apply(cc_tf, 1L, zscore))
rownames(tiger_tf) <- inferable_tfs; colnames(tiger_tf) <- tiger_ids
rownames(cc_tf) <- inferable_tfs; colnames(cc_tf) <- cc_ids
tiger_data[, `:=`(
  SPI1_activity = as.numeric(tiger_tf["SPI1", sample_id]),
  NFKB1_activity = as.numeric(tiger_tf["NFKB1", sample_id]),
  RELA_activity = as.numeric(tiger_tf["RELA", sample_id]),
  NFKB2_activity = as.numeric(tiger_tf["NFKB2", sample_id]),
  NFKB2_expression = zscore(tiger_log["NFKB2", sample_id]),
  RELB_expression = zscore(tiger_log["RELB", sample_id])
)]
cc_data[, `:=`(
  SPI1_activity = as.numeric(cc_tf["SPI1", sample_id]),
  NFKB1_activity = as.numeric(cc_tf["NFKB1", sample_id]),
  RELA_activity = as.numeric(cc_tf["RELA", sample_id]),
  NFKB2_activity = as.numeric(cc_tf["NFKB2", sample_id]),
  NFKB2_expression = zscore(cc_log["NFKB2", sample_id]),
  RELB_expression = zscore(cc_log["RELB", sample_id])
)]
tiger_data[, NFKB1_RELA_activity := zscore((NFKB1_activity + RELA_activity) / 2)]
cc_data[, NFKB1_RELA_activity := zscore((NFKB1_activity + RELA_activity) / 2)]

# -------------------------------------------------------------------------
# Cohort contracts, exposure effects, route edges, and competing branches
# -------------------------------------------------------------------------
cohorts <- list(
  "Melanoma: Degradation-High" = list(
    data = tiger_data,
    covars = c("therapy_short", "response_group", "age_z", "gender"),
    group_exposure = "target_group_z", continuous_exposure = "target_continuous_z"
  ),
  "ccRCC: Sarcosine-Low" = list(
    data = cc_data,
    covars = c("batch", "sex", "age", "grade"),
    group_exposure = "target_group_z", continuous_exposure = "target_continuous_z"
  )
)

models <- data.table(
  model_order = 1:3,
  model = c("Grouped primary", "Continuous sensitivity", "Grouped + named T-cell/APC context"),
  exposure_type = c("group", "continuous", "group"),
  extras = c("", "", "T_cell_context_z;APC_myeloid_context_z"),
  reader_facing = c(TRUE, FALSE, FALSE)
)

node_columns <- c("APC", "TCR", "CD28_family", "Canonical_NFkB", "Noncanonical_NFkB", "IFNG_response")
node_effects <- rbindlist(lapply(names(cohorts), function(cohort_name) {
  spec <- cohorts[[cohort_name]]
  rbindlist(lapply(seq_len(nrow(models)), function(i) {
    exposure <- if (models$exposure_type[i] == "group") spec$group_exposure else spec$continuous_exposure
    extras <- if (nzchar(models$extras[i])) strsplit(models$extras[i], ";", fixed = TRUE)[[1]] else character()
    out <- rbindlist(lapply(node_columns, function(node) {
      fit_term(spec$data, node, exposure, exposure, spec$covars, extras,
               cohort_name, models$model[i])
    }))
    out[, BH_q := p.adjust(p_value, method = "BH")]
    out
  }))
}))

route_catalog <- data.table(
  route_order = 1:4,
  route = c(
    "APC -> TCR -> canonical NF-kappaB -> IFNG",
    "APC -> TCR -> noncanonical NF-kappaB -> IFNG",
    "APC -> broad CD28-family regulation -> canonical NF-kappaB -> IFNG",
    "APC -> broad CD28-family regulation -> noncanonical NF-kappaB -> IFNG"
  ),
  middle = c("TCR", "TCR", "CD28_family", "CD28_family"),
  branch = c("Canonical_NFkB", "Noncanonical_NFkB", "Canonical_NFkB", "Noncanonical_NFkB"),
  role = c("Shared TCR route", "Shared TCR route", "Mixed CD28-family sensitivity", "Mixed CD28-family sensitivity")
)

route_edges <- list(); route_boot <- list(); counter <- 0L
for (cohort_name in names(cohorts)) {
  spec <- cohorts[[cohort_name]]
  for (m in seq_len(nrow(models))) {
    exposure <- if (models$exposure_type[m] == "group") spec$group_exposure else spec$continuous_exposure
    extras <- if (nzchar(models$extras[m])) strsplit(models$extras[m], ";", fixed = TRUE)[[1]] else character()
    for (r in seq_len(nrow(route_catalog))) {
      forward_nodes <- c("APC", route_catalog$middle[r], route_catalog$branch[r], "IFNG_response")
      reverse_nodes <- rev(forward_nodes)
      counter <- counter + 1L
      route_edges[[counter]] <- rbindlist(list(
        serial_edges(spec$data, forward_nodes, exposure, spec$covars, extras,
                     cohort_name, models$model[m], route_catalog$route[r], "Forward"),
        serial_edges(spec$data, reverse_nodes, exposure, spec$covars, extras,
                     cohort_name, models$model[m], route_catalog$route[r], "Reverse")
      ))
      route_boot[[counter]] <- rbindlist(list(
        serial_bootstrap(spec$data, forward_nodes, exposure, spec$covars, extras,
                         cohort_name, models$model[m], route_catalog$route[r], "Forward",
                         B = 5000L, seed = 260902L + 1000L * counter + 1L),
        serial_bootstrap(spec$data, reverse_nodes, exposure, spec$covars, extras,
                         cohort_name, models$model[m], route_catalog$route[r], "Reverse",
                         B = 5000L, seed = 260902L + 1000L * counter + 2L)
      ))
    }
  }
}
route_edges <- rbindlist(route_edges)
route_edges[, BH_q := p.adjust(p_value, method = "BH"), by = .(cohort, model, direction, edge_order)]
route_bootstrap <- rbindlist(route_boot)

joint_results <- rbindlist(lapply(names(cohorts), function(cohort_name) {
  spec <- cohorts[[cohort_name]]
  rbindlist(lapply(seq_len(nrow(models)), function(i) {
    exposure <- if (models$exposure_type[i] == "group") spec$group_exposure else spec$continuous_exposure
    extras <- if (nzchar(models$extras[i])) strsplit(models$extras[i], ";", fixed = TRUE)[[1]] else character()
    predictors <- c("TCR", "CD28_family", "Canonical_NFkB", "Noncanonical_NFkB")
    rhs <- unique(c(exposure, "APC", predictors, spec$covars, extras))
    form <- as.formula(paste("IFNG_response ~", paste(rhs, collapse = " + ")))
    x <- model.matrix(form, data = spec$data)
    if (qr(x)$rank != ncol(x)) stop("Rank-deficient joint model: ", cohort_name, " / ", models$model[i])
    fit <- lm(form, data = spec$data)
    out <- rbindlist(lapply(predictors, function(term) {
      extract_term(fit, term, cohort_name, models$model[i], "IFNG_response", term,
                   route = "Joint competition", edge = paste(term, "-> IFNG | all branches"))
    }))
    out[, BH_q := p.adjust(p_value, method = "BH")]
    out
  }))
}))

tf_effects <- rbindlist(lapply(names(cohorts), function(cohort_name) {
  spec <- cohorts[[cohort_name]]
  rbindlist(lapply(seq_len(nrow(models)), function(i) {
    exposure <- if (models$exposure_type[i] == "group") spec$group_exposure else spec$continuous_exposure
    extras <- if (nzchar(models$extras[i])) strsplit(models$extras[i], ";", fixed = TRUE)[[1]] else character()
    out <- rbindlist(lapply(c("SPI1_activity", "NFKB1_activity", "RELA_activity", "NFKB1_RELA_activity", "NFKB2_activity"), function(node) {
      fit_term(spec$data, node, exposure, exposure, spec$covars, extras,
               cohort_name, models$model[i])
    }))
    out[, BH_q := p.adjust(p_value, method = "BH")]
    out
  }))
}))

tf_transcript_effects <- rbindlist(lapply(names(cohorts), function(cohort_name) {
  spec <- cohorts[[cohort_name]]
  rbindlist(lapply(seq_len(nrow(models)), function(i) {
    exposure <- if (models$exposure_type[i] == "group") spec$group_exposure else spec$continuous_exposure
    extras <- if (nzchar(models$extras[i])) strsplit(models$extras[i], ";", fixed = TRUE)[[1]] else character()
    out <- rbindlist(lapply(c("NFKB2_expression", "RELB_expression"), function(node) {
      fit_term(spec$data, node, exposure, exposure, spec$covars, extras,
               cohort_name, models$model[i])
    }))
    out[, BH_q := p.adjust(p_value, method = "BH")]
    out
  }))
}))

tf_links <- rbindlist(lapply(names(cohorts), function(cohort_name) {
  spec <- cohorts[[cohort_name]]
  exposure <- spec$group_exposure
  out <- rbindlist(list(
    fit_term(spec$data, "APC", "SPI1_activity", exposure, spec$covars, character(),
             cohort_name, "Grouped primary", edge = "SPI1 activity -> APC cross-presentation"),
    fit_term(spec$data, "Canonical_NFkB", "NFKB1_RELA_activity", exposure, spec$covars, character(),
             cohort_name, "Grouped primary", edge = "NFKB1/RELA activity -> canonical NF-kappaB program"),
    fit_term(spec$data, "Noncanonical_NFkB", "NFKB2_activity", exposure, spec$covars, character(),
             cohort_name, "Grouped primary", edge = "NFKB2 activity -> noncanonical NF-kappaB program")
  ))
  out[, BH_q := p.adjust(p_value, method = "BH")]
  out
}))
tf_support <- rbindlist(list(
  copy(tf_effects)[, analysis := "Exposure effect"],
  copy(tf_transcript_effects)[, analysis := "Focal transcript exposure effect"],
  copy(tf_links)[, analysis := "Program link"]
), fill = TRUE)

# -------------------------------------------------------------------------
# Transparent evidence score and final leading candidate
# -------------------------------------------------------------------------
get_screen <- function(label) screen[candidate == label]
cohort_names <- names(cohorts)
primary_boot <- route_bootstrap[model == "Grouped primary"]
continuous_boot <- route_bootstrap[model == "Continuous sensitivity"]
primary_edges <- route_edges[model == "Grouped primary" & direction == "Forward"]
primary_joint <- joint_results[model == "Grouped primary"]

fixed_effect_meta <- function(x, branch) {
  w <- 1 / x$standard_error^2
  estimate <- sum(w * x$standardized_beta) / sum(w)
  se <- sqrt(1 / sum(w))
  z <- estimate / se
  Q <- sum(w * (x$standardized_beta - estimate)^2)
  Q_df <- nrow(x) - 1L
  data.table(
    branch = branch, cohorts_n = nrow(x), fixed_effect_beta = estimate,
    fixed_effect_SE = se, fixed_effect_CI_low = estimate - qnorm(0.975) * se,
    fixed_effect_CI_high = estimate + qnorm(0.975) * se,
    fixed_effect_z = z, fixed_effect_p = 2 * pnorm(abs(z), lower.tail = FALSE),
    Cochran_Q = Q, Q_df = Q_df,
    heterogeneity_p = pchisq(Q, df = Q_df, lower.tail = FALSE),
    I2_percent = if (Q > 0) max(0, (Q - Q_df) / Q) * 100 else 0,
    direction_consistent_positive = all(x$standardized_beta > 0),
    cohort_betas = paste(paste(x$cohort, formatC(x$standardized_beta, digits = 5, format = "fg"), sep = ":"), collapse = ";")
  )
}
branch_meta <- rbindlist(lapply(c("Canonical_NFkB", "Noncanonical_NFkB"), function(branch) {
  fixed_effect_meta(primary_joint[predictor == branch], branch)
}))

route_evidence <- rbindlist(lapply(seq_len(nrow(route_catalog)), function(i) {
  route_name <- route_catalog$route[i]
  middle_label <- if (route_catalog$middle[i] == "TCR") "TCR signaling" else "CD28-family regulation (broad)"
  branch_label <- if (route_catalog$branch[i] == "Canonical_NFkB") "Canonical NF-kappaB" else "Noncanonical NF-kappaB"
  anchors <- c("APC cross-presentation", middle_label, branch_label, "IFNG response")
  anchor_pass <- all(screen[candidate %in% anchors, exact_replication_pass])
  specific_middle_pass <- if (route_catalog$middle[i] == "TCR") {
    get_screen("TCR signaling")$exact_replication_pass
  } else {
    get_screen("CD28 co-stimulation (specific)")$exact_replication_pass
  }
  edge_ok <- sapply(cohort_names, function(cn) {
    z <- primary_edges[cohort == cn & route == route_name]
    nrow(z) == 4L && all(z$standardized_beta > 0 & z$p_value < 0.05)
  })
  forward_ok <- sapply(cohort_names, function(cn) {
    z <- primary_boot[cohort == cn & route == route_name & direction == "Forward"]
    nrow(z) == 1L && z$bootstrap_CI_low > 0
  })
  reverse_not_supported <- sapply(cohort_names, function(cn) {
    z <- primary_boot[cohort == cn & route == route_name & direction == "Reverse"]
    nrow(z) == 1L && !z$CI_excludes_zero
  })
  joint_ok <- sapply(cohort_names, function(cn) {
    z <- primary_joint[cohort == cn & predictor == route_catalog$branch[i]]
    nrow(z) == 1L && z$standardized_beta > 0 & z$BH_q < 0.05
  })
  continuous_positive <- sapply(cohort_names, function(cn) {
    z <- continuous_boot[cohort == cn & route == route_name & direction == "Forward"]
    nrow(z) == 1L && z$serial_product > 0
  })
  shared_le <- get_screen(branch_label)$shared_leading_edge_n
  criteria <- c(
    anchor_pass,
    specific_middle_pass,
    edge_ok[[1]], edge_ok[[2]],
    forward_ok[[1]], forward_ok[[2]],
    reverse_not_supported[[1]], reverse_not_supported[[2]],
    joint_ok[[1]], joint_ok[[2]],
    continuous_positive[[1]], continuous_positive[[2]],
    shared_le >= 10L
  )
  data.table(
    route_order = route_catalog$route_order[i], route = route_name, role = route_catalog$role[i],
    pathway_anchors_replicated = anchor_pass,
    receptor_specificity_pass = specific_middle_pass,
    TIGER_all_forward_edges_positive_p_lt_0_05 = edge_ok[[1]],
    ccRCC_all_forward_edges_positive_p_lt_0_05 = edge_ok[[2]],
    TIGER_forward_serial_CI_positive = forward_ok[[1]],
    ccRCC_forward_serial_CI_positive = forward_ok[[2]],
    TIGER_reverse_serial_not_supported = reverse_not_supported[[1]],
    ccRCC_reverse_serial_not_supported = reverse_not_supported[[2]],
    TIGER_branch_joint_positive_p_lt_0_05 = joint_ok[[1]],
    ccRCC_branch_joint_positive_p_lt_0_05 = joint_ok[[2]],
    TIGER_continuous_serial_positive = continuous_positive[[1]],
    ccRCC_continuous_serial_positive = continuous_positive[[2]],
    branch_shared_leading_edge_n = shared_le,
    branch_shared_leading_edge_ge_10 = shared_le >= 10L,
    branch_direction_consistent_positive = branch_meta[branch == route_catalog$branch[i], direction_consistent_positive],
    branch_fixed_effect_beta = branch_meta[branch == route_catalog$branch[i], fixed_effect_beta],
    branch_fixed_effect_p = branch_meta[branch == route_catalog$branch[i], fixed_effect_p],
    branch_I2_percent = branch_meta[branch == route_catalog$branch[i], I2_percent],
    criteria_passed = sum(criteria), criteria_total = length(criteria)
  )
}))

setorder(route_evidence, -criteria_passed, route_order)
eligible_branches <- branch_meta[direction_consistent_positive == TRUE & fixed_effect_p < 0.05]
if (nrow(eligible_branches) == 0L) stop("No direction-consistent common NF-kappaB branch")
leading_branch <- eligible_branches[order(fixed_effect_p), branch][1]
leading_route <- route_catalog[middle == "TCR" & branch == leading_branch, route]
if (length(leading_route) != 1L) stop("Could not select one TCR-based leading route")
leading_catalog <- route_catalog[route == leading_route]
complete_serial_replicated <- all(
  primary_boot[route == leading_route & direction == "Forward", bootstrap_CI_low] > 0
)
final_decision <- data.table(
  rank = frank(-as.integer(route_evidence$route == leading_route), ties.method = "first"),
  route = route_evidence$route,
  criteria_passed = route_evidence$criteria_passed,
  criteria_total = route_evidence$criteria_total,
  designation = fifelse(
    route_evidence$route == leading_route,
    "Direction-consistent branch association inside the shared co-enrichment program; complete serial route not replicated",
    fifelse(grepl("broad CD28", route_evidence$route),
            "Excluded as a specific CD28 route; broad mixed-program sensitivity only",
            "Lower-priority comparator")
  ),
  complete_serial_replicated_in_both_cohorts = route_evidence$route == leading_route & complete_serial_replicated,
  claim_boundary = "Exploratory cross-sectional bulk-RNA route; not a proven temporal or same-cell mechanism"
)
setorder(final_decision, rank)

# -------------------------------------------------------------------------
# Reader-facing figures
# -------------------------------------------------------------------------
plot_screen <- rbindlist(list(
  screen[, .(candidate_order, candidate, cohort = "Melanoma\nDegradation-High",
             NES = tiger_harmonized_NES, BH_q = tiger_q)],
  screen[, .(candidate_order, candidate, cohort = "ccRCC\nSarcosine-Low",
             NES = ccrcc_harmonized_NES, BH_q = ccrcc_q)]
))
plot_screen[, candidate := factor(candidate, levels = rev(candidate_catalog$candidate))]
plot_screen[, significant := BH_q < 0.05]
plot_screen[, q_plot := pmin(-log10(pmax(BH_q, 1e-50)), 20)]

p_a <- ggplot(plot_screen, aes(cohort, candidate)) +
  geom_tile(aes(fill = NES), color = "white", linewidth = 0.8) +
  geom_point(aes(size = q_plot, shape = significant), fill = "black", color = "black") +
  scale_shape_manual(values = c(`TRUE` = 21, `FALSE` = 1),
                     labels = c(`TRUE` = "BH q < 0.05", `FALSE` = "BH q >= 0.05")) +
  scale_size_continuous(name = expression(-log[10]("BH q")), range = c(1.5, 6), breaks = c(2, 5, 10, 20)) +
  scale_fill_gradient2(low = "#2A9D8F", mid = "white", high = "#C73E3A", midpoint = 0,
                       name = "Harmonized\nNES") +
  labs(title = "a  Broad candidate screen", subtitle = "Exact MSigDB pathways; red denotes the target state",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(size = 9), plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(color = "grey35"), legend.position = "right")

joint_branch_plot <- joint_results[
  model == "Grouped primary" & predictor %in% c("Canonical_NFkB", "Noncanonical_NFkB")
]
branch_labels <- c(Canonical_NFkB = "Canonical NF-kappaB", Noncanonical_NFkB = "Noncanonical NF-kappaB")
joint_branch_plot[, branch_label := factor(branch_labels[predictor], levels = rev(unname(branch_labels)))]
p_b <- ggplot(joint_branch_plot, aes(standardized_beta, branch_label, color = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0, linewidth = 0.8,
                 position = position_dodge(width = 0.42)) +
  geom_point(aes(shape = BH_q < 0.05), size = 3.1, position = position_dodge(width = 0.42)) +
  scale_color_manual(values = cohort_palette) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(title = "b  Joint branch-to-IFNG competition",
       subtitle = "Both branches entered together; filled points pass within-panel BH q < 0.05",
       x = "Adjusted association with IFNG response", y = NULL, color = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(color = "grey35"), legend.position = "bottom")

primary_serial_plot <- primary_boot[direction == "Forward"]
primary_serial_plot[, short_route := fifelse(
  grepl("TCR.*noncanonical", route), "TCR -> noncanonical NF-kappaB",
  fifelse(grepl("TCR.*canonical", route), "TCR -> canonical NF-kappaB",
          fifelse(grepl("CD28.*noncanonical", route), "Broad CD28-family -> noncanonical NF-kappaB",
                  "Broad CD28-family -> canonical NF-kappaB"))
)]
primary_serial_plot[, short_route := factor(short_route, levels = rev(c(
  "TCR -> canonical NF-kappaB", "TCR -> noncanonical NF-kappaB",
  "Broad CD28-family -> canonical NF-kappaB", "Broad CD28-family -> noncanonical NF-kappaB"
)))]
p_c <- ggplot(primary_serial_plot, aes(serial_product, short_route, color = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbar(aes(xmin = bootstrap_CI_low, xmax = bootstrap_CI_high), orientation = "y", width = 0,
                linewidth = 0.9, position = position_dodge(width = 0.45)) +
  geom_point(aes(shape = CI_excludes_zero), size = 3.2, position = position_dodge(width = 0.45)) +
  scale_color_manual(values = cohort_palette) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(title = "c  Forward serial association", subtitle = "Exposure -> APC -> receptor-family node -> branch -> IFNG; 5,000-bootstrap 95% CI",
       x = "Product of standardized path coefficients", y = NULL, color = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(color = "grey35"), legend.position = "bottom")

evidence_display <- copy(route_evidence)
evidence_display[, `:=`(
  both_reverse_orders_not_supported = TIGER_reverse_serial_not_supported & ccRCC_reverse_serial_not_supported,
  both_branch_joint_BH_positive = TIGER_branch_joint_positive_p_lt_0_05 & ccRCC_branch_joint_positive_p_lt_0_05,
  both_continuous_serial_positive = TIGER_continuous_serial_positive & ccRCC_continuous_serial_positive
)]
evidence_long <- melt(
  evidence_display,
  id.vars = c("route_order", "route", "role", "criteria_passed", "criteria_total", "branch_shared_leading_edge_n"),
  measure.vars = c(
    "pathway_anchors_replicated", "receptor_specificity_pass",
    "TIGER_forward_serial_CI_positive", "ccRCC_forward_serial_CI_positive",
    "branch_direction_consistent_positive", "both_branch_joint_BH_positive",
    "both_continuous_serial_positive", "both_reverse_orders_not_supported",
    "branch_shared_leading_edge_ge_10"
  ),
  variable.name = "criterion", value.name = "pass"
)
criterion_labels <- c(
  pathway_anchors_replicated = "Exact pathways replicated",
  receptor_specificity_pass = "Receptor-specific gate",
  TIGER_forward_serial_CI_positive = "Melanoma: forward serial CI > 0",
  ccRCC_forward_serial_CI_positive = "ccRCC: forward serial CI > 0",
  branch_direction_consistent_positive = "Branch direction agrees",
  both_branch_joint_BH_positive = "Branch survives both cohorts",
  both_continuous_serial_positive = "Continuous direction agrees",
  both_reverse_orders_not_supported = "Reverse unsupported in both",
  branch_shared_leading_edge_ge_10 = "Shared branch leading edge >= 10"
)
route_short_labels <- c(
  "APC -> TCR -> canonical NF-kappaB -> IFNG" = "TCR / canonical NF-kappaB",
  "APC -> TCR -> noncanonical NF-kappaB -> IFNG" = "TCR / noncanonical NF-kappaB",
  "APC -> broad CD28-family regulation -> canonical NF-kappaB -> IFNG" = "Broad CD28-family / canonical NF-kappaB",
  "APC -> broad CD28-family regulation -> noncanonical NF-kappaB -> IFNG" = "Broad CD28-family / noncanonical NF-kappaB"
)
evidence_long[, criterion_label := factor(criterion_labels[as.character(criterion)], levels = rev(criterion_labels))]
evidence_long[, route_label := factor(route_short_labels[route], levels = unname(route_short_labels))]
p_d <- ggplot(evidence_long, aes(criterion_label, route_label, fill = pass)) +
  geom_tile(color = "white", linewidth = 0.7) +
  scale_fill_manual(values = c(`TRUE` = "#C73E3A", `FALSE` = "#E6E6E6"), guide = "none") +
  coord_cartesian(clip = "off") +
  labs(title = "d  Hierarchical evidence gates", subtitle = "Red = criterion passed; no complete serial route passed both cohorts",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(face = "bold"), plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(color = "grey35"))

combined <- (p_a | p_b) / (p_c | p_d) +
  plot_annotation(
    title = "Cross-cohort pathway screen and mechanism-replication diagnostics",
    subtitle = "Target states: melanoma sarcosine-degradation High and ccRCC measured Sarcosine Low",
    caption = "Exploratory bulk-RNA associations. Pathway enrichment and serial models do not establish temporal causality or same-cell signaling.",
    theme = theme(plot.title = element_text(face = "bold", size = 21),
                  plot.subtitle = element_text(size = 13, color = "grey25"),
                  plot.caption = element_text(size = 9, color = "grey40"))
  )
ggsave(file.path(figure_dir, "Fig_Common_Immune_Mechanism_Broad_to_Final_Narrowing.png"), combined,
       width = 18, height = 13.5, dpi = 320, bg = "white")
ggsave(file.path(figure_dir, "Fig_Common_Immune_Mechanism_Broad_to_Final_Narrowing.pdf"), combined,
       width = 18, height = 13.5, bg = "white")

# Dedicated reader-facing evidence figure explaining why the selected downstream
# axis is called shared across the two cohorts. This separates three questions:
# exact-set enrichment, gene-level leading-edge overlap, and branch competition.
shared_axis_candidates <- c(
  "APC cross-presentation", "TCR signaling",
  "Noncanonical NF-kappaB", "IFNG response"
)
shared_axis_labels <- c(
  "APC cross-presentation" = "APC cross-presentation",
  "TCR signaling" = "TCR signaling",
  "Noncanonical NF-kappaB" = "TNFR2-related\nnoncanonical NF-kappaB",
  "IFNG response" = "IFNG response"
)
shared_axis_screen <- screen[candidate %in% shared_axis_candidates]
if (nrow(shared_axis_screen) != length(shared_axis_candidates) ||
    !all(shared_axis_screen$exact_replication_pass)) {
  stop("The four selected shared-axis sets did not all pass the exact two-cohort gate")
}
format_reader_q <- function(x) {
  if (x < 0.001) formatC(x, format = "e", digits = 1) else formatC(x, format = "f", digits = 3)
}
shared_axis_long <- rbindlist(list(
  shared_axis_screen[, .(
    candidate, cohort = "Melanoma\nDegradation-High",
    NES = tiger_harmonized_NES, BH_q = tiger_q
  )],
  shared_axis_screen[, .(
    candidate, cohort = "ccRCC\nSarcosine-Low",
    NES = ccrcc_harmonized_NES, BH_q = ccrcc_q
  )]
))
shared_axis_long[, axis_label := shared_axis_labels[candidate]]
shared_axis_long[, axis_label := factor(axis_label, levels = rev(unname(shared_axis_labels[shared_axis_candidates])))]
shared_axis_long[, cohort := factor(cohort, levels = c("Melanoma\nDegradation-High", "ccRCC\nSarcosine-Low"))]
shared_axis_long[, cell_label := paste0(
  "NES ", sprintf("%.2f", NES), "\nq ", vapply(BH_q, format_reader_q, character(1))
)]

p_shared_a <- ggplot(shared_axis_long, aes(cohort, axis_label)) +
  geom_tile(aes(fill = NES), color = "black", linewidth = 0.9) +
  geom_text(aes(label = cell_label), size = 3.8, lineheight = 0.95) +
  scale_fill_gradient2(low = "#2A9D8F", mid = "white", high = "#C73E3A", midpoint = 0,
                       limits = c(-3.3, 3.3), name = "Target-direction\nNES") +
  labs(
    title = "a  Independent exact-set replication",
    subtitle = "Every displayed MSigDB set has positive NES and BH q < 0.05 in each cohort",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(), axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 10), plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey35"), legend.position = "right"
  )

shared_axis_overlap <- copy(shared_axis_screen)
shared_axis_overlap[, axis_label := shared_axis_labels[candidate]]
shared_axis_overlap[, axis_label := factor(axis_label, levels = rev(unname(shared_axis_labels[shared_axis_candidates])))]
shared_axis_overlap[, overlap_label := fifelse(
  candidate == "Noncanonical NF-kappaB",
  paste0("n=", shared_leading_edge_n, ": NFKB2, RELB, TRAF3"),
  paste0("n=", shared_leading_edge_n)
)]
p_shared_b <- ggplot(shared_axis_overlap, aes(leading_edge_jaccard, axis_label)) +
  geom_segment(aes(x = 0, xend = leading_edge_jaccard, yend = axis_label),
               linewidth = 0.8, color = "grey70") +
  geom_point(aes(size = shared_leading_edge_n), color = "#6A51A3", alpha = 0.9) +
  geom_text(aes(x = pmin(leading_edge_jaccard + 0.018, 0.36), label = overlap_label),
            hjust = 0, size = 3.6, color = "grey20") +
  scale_x_continuous(
    labels = function(x) paste0(round(100 * x), "%"),
    limits = c(0, 0.43), breaks = seq(0, 0.4, 0.1), expand = expansion(mult = c(0, 0.02))
  ) +
  scale_size_continuous(name = "Shared leading-\nedge genes", range = c(3, 9)) +
  labs(
    title = "b  Cohort leading-edge overlap",
    subtitle = "Intersection of the two cohort-specific leading edges for the same exact set",
    x = "Leading-edge Jaccard overlap", y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(), axis.text.y = element_text(size = 9),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey35"), legend.position = "right"
  )

shared_branch_forest <- copy(joint_branch_plot)
shared_branch_meta <- branch_meta[branch %in% c("Canonical_NFkB", "Noncanonical_NFkB")]
shared_branch_meta[, branch_label := factor(branch_labels[branch], levels = rev(unname(branch_labels)))]
shared_branch_meta[, meta_label := paste0(
  "2-cohort summary: beta ", sprintf("%.3f", fixed_effect_beta),
  "; p ", vapply(fixed_effect_p, format_reader_q, character(1)),
  "; I2 ", sprintf("%.1f%%", I2_percent)
)]
p_shared_c <- ggplot(shared_branch_forest, aes(standardized_beta, branch_label, color = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0,
                linewidth = 0.9, position = position_dodge(width = 0.42)) +
  geom_point(aes(shape = BH_q < 0.05), size = 3.3, position = position_dodge(width = 0.42)) +
  geom_text(
    data = shared_branch_meta,
    aes(x = 0.60, y = branch_label, label = meta_label),
    inherit.aes = FALSE, hjust = 0, size = 3.5, color = "grey20"
  ) +
  scale_color_manual(values = cohort_palette) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "BH q < 0.05", `FALSE` = "BH q >= 0.05")) +
  coord_cartesian(xlim = c(-0.38, 1.02), clip = "off") +
  labs(
    title = "c  Why the noncanonical branch was retained",
    subtitle = "Both branches entered together: noncanonical is positive in both cohorts; canonical reverses direction",
    x = "Adjusted branch association with IFNG response", y = NULL, color = NULL, shape = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey35"), legend.position = "bottom",
    plot.margin = margin(5.5, 55, 5.5, 5.5)
  )

shared_evidence_figure <- (p_shared_a | p_shared_b) / p_shared_c +
  plot_layout(heights = c(1.15, 0.85)) +
  plot_annotation(
    title = "Cross-cohort evidence for the shared co-enriched immune program",
    subtitle = "Target states: melanoma Degradation-High and ccRCC measured Sarcosine-Low; cohorts were tested independently",
    caption = paste0(
      "Conclusion: APC cross-presentation, TCR signaling, the TNFR2-related noncanonical NF-kappaB expression program, and IFNG response are shared pathway-level signals.\n",
      "Limits: the noncanonical leading-edge intersection is small (3 genes), common NFKB2/RELB TF activation was not demonstrated, and the complete exposure-linked serial route did not replicate in ccRCC."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 21),
      plot.subtitle = element_text(size = 12.5, color = "grey25"),
      plot.caption = element_text(size = 9, color = "grey35", hjust = 0)
    )
  )
ggsave(file.path(figure_dir, "Fig_Cross_Cohort_Evidence_for_Selected_Immune_Axis.png"), shared_evidence_figure,
       width = 18, height = 12.5, dpi = 320, bg = "white")
ggsave(file.path(figure_dir, "Fig_Cross_Cohort_Evidence_for_Selected_Immune_Axis.pdf"), shared_evidence_figure,
       width = 18, height = 12.5, bg = "white")

# Reader-facing final schematic: deliberately concise.
schematic_nodes <- data.table(
  x = 1:4, y = 1,
  label = c("APC cross-\npresentation", "TCR signaling",
            if (leading_catalog$branch == "Canonical_NFkB") "Canonical NF-kappaB\n(NFKB1 / RELA)" else "TNFR2-related\nnoncanonical NF-kappaB\nexpression program",
            "IFNG-response\nprogram")
)
schematic_plus <- data.table(x = c(1.5, 2.5, 3.5), y = 1, label = "+")
p_schematic <- ggplot() +
  geom_text(data = schematic_plus, aes(x, y, label = label), color = "#6A51A3",
            size = 9, fontface = "bold") +
  geom_label(data = schematic_nodes, aes(x, y, label = label), size = 5.2,
             linewidth = 0.8, label.padding = unit(0.28, "lines"), fontface = "bold") +
  annotate("text", x = 2.5, y = 0.50,
           label = "Co-enriched gene-set program only: no temporal order and no shared exposure-linked serial mechanism inferred",
           size = 4.3, color = "grey30") +
  coord_cartesian(xlim = c(0.55, 4.45), ylim = c(0.30, 1.45), clip = "off") +
  labs(title = "Shared co-enriched immune program across both cohorts",
       subtitle = "Plus signs denote joint membership in the replicated GSEA signal; they are not signaling arrows") +
  theme_void(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 21),
        plot.subtitle = element_text(size = 12.5, color = "grey35"),
        plot.margin = margin(25, 35, 25, 35))
ggsave(file.path(figure_dir, "Fig_Leading_Common_Immune_Candidate_Route.png"), p_schematic,
       width = 17, height = 5.8, dpi = 320, bg = "white")
ggsave(file.path(figure_dir, "Fig_Leading_Common_Immune_Candidate_Route.pdf"), p_schematic,
       width = 17, height = 5.8, bg = "white")

# Targeted TF follow-up for the selected noncanonical branch.
tf_reader <- rbindlist(list(
  tf_effects[model == "Grouped primary" & outcome == "NFKB2_activity",
             .(cohort, evidence = "NFKB2 inferred TF activity", standardized_beta, CI_low, CI_high, BH_q)],
  tf_transcript_effects[model == "Grouped primary",
                        .(cohort,
                          evidence = fifelse(outcome == "NFKB2_expression", "NFKB2 transcript", "RELB transcript"),
                          standardized_beta, CI_low, CI_high, BH_q)]
), use.names = TRUE)
tf_reader[, evidence := factor(evidence, levels = rev(c(
  "NFKB2 inferred TF activity", "NFKB2 transcript", "RELB transcript"
)))]
p_tf_a <- ggplot(tf_reader, aes(standardized_beta, evidence, color = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0,
                linewidth = 0.8, position = position_dodge(width = 0.42)) +
  geom_point(aes(shape = BH_q < 0.05), size = 3.1, position = position_dodge(width = 0.42)) +
  scale_color_manual(values = cohort_palette) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "BH q < 0.05", `FALSE` = "BH q >= 0.05")) +
  labs(title = "a  Target-state TF/transcript contrasts",
       subtitle = "Adjusted target-minus-reference effects",
       x = "Standardized effect", y = NULL, color = NULL, shape = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(color = "grey35"), legend.position = "bottom")

nfkb2_link_plot <- tf_links[edge == "NFKB2 activity -> noncanonical NF-kappaB program"]
p_tf_b <- ggplot(nfkb2_link_plot, aes(standardized_beta, cohort, color = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0, linewidth = 0.8) +
  geom_point(aes(shape = BH_q < 0.05), size = 3.2) +
  scale_color_manual(values = cohort_palette) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), guide = "none") +
  labs(title = "b  NFKB2 activity-to-branch association",
       subtitle = "Exposure-adjusted link to the gene-disjoint noncanonical program",
       x = "Standardized association", y = NULL, color = NULL) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 15),
        plot.subtitle = element_text(color = "grey35"), legend.position = "none")

tf_combined <- (p_tf_a | p_tf_b) +
  plot_annotation(
    title = "Targeted TF check: common NFKB2/RELB activation is not supported",
    subtitle = "Strict cross-cohort-consensus signed DoRothEA A/B regulons exclude every route-score and IFNG-response gene",
    caption = "NFKB2 is estimable from 10 retained targets. RELB activity is not estimable at the prespecified minimum: 5 retained targets (<10). Transcript expression is not TF activity.",
    theme = theme(plot.title = element_text(face = "bold", size = 20),
                  plot.subtitle = element_text(size = 12.5, color = "grey30"),
                  plot.caption = element_text(size = 9, color = "grey40"))
  )
ggsave(file.path(figure_dir, "Fig_Common_Noncanonical_NFkB_TF_Check.png"), tf_combined,
       width = 16, height = 7.5, dpi = 320, bg = "white")
ggsave(file.path(figure_dir, "Fig_Common_Noncanonical_NFkB_TF_Check.pdf"), tf_combined,
       width = 16, height = 7.5, bg = "white")

# Direct visualization of group separation for the gene-disjoint sample-level
# route scores. This is deliberately shown separately from the exact-set GSEA:
# the two analyses answer different questions and should not be conflated.
route_score_values <- rbindlist(list(
  tiger_data[, .(
    cohort = "Melanoma", sample_id,
    group = fifelse(target_group == "Target", "Degradation-High", "Degradation-Low"),
    APC, TCR, Noncanonical_NFkB, IFNG_response
  )],
  cc_data[, .(
    cohort = "ccRCC", sample_id,
    group = fifelse(target_group == "Target", "Sarcosine-Low", "Sarcosine-High"),
    APC, TCR, Noncanonical_NFkB, IFNG_response
  )]
), use.names = TRUE)
route_score_values <- melt(
  route_score_values,
  id.vars = c("cohort", "sample_id", "group"),
  variable.name = "node", value.name = "score"
)
route_score_labels <- c(
  APC = "APC cross-presentation", TCR = "TCR signaling",
  Noncanonical_NFkB = "TNFR2-related\nnoncanonical NF-kappaB",
  IFNG_response = "IFNG response"
)
route_score_values[, node_label := factor(route_score_labels[as.character(node)],
                                           levels = unname(route_score_labels))]

route_group_effects <- node_effects[
  model == "Grouped primary" & outcome %in% names(route_score_labels)
]
route_group_effects[, node_label := factor(route_score_labels[outcome], levels = rev(unname(route_score_labels)))]
route_group_effects[, q_label := paste0("q ", vapply(BH_q, format_reader_q, character(1)))]

score_limits <- range(route_score_values$score, finite = TRUE)
score_padding <- 0.10 * diff(score_limits)
score_limits <- score_limits + c(-score_padding, score_padding)

make_group_distribution <- function(cohort_name, low_reference, high_target, panel_title) {
  plot_data <- route_score_values[cohort == cohort_name]
  plot_data[, group := factor(as.character(group), levels = c(low_reference, high_target))]
  group_palette <- setNames(c("#2A9D8F", "#C73E3A"), c(low_reference, high_target))
  ggplot(plot_data, aes(group, score, color = group)) +
    geom_boxplot(width = 0.48, outlier.shape = NA, linewidth = 0.7) +
    geom_jitter(width = 0.13, height = 0, alpha = 0.55, size = 1.25) +
    facet_wrap(~node_label, nrow = 1) +
    scale_color_manual(values = group_palette, guide = "none") +
    coord_cartesian(ylim = score_limits) +
    labs(title = panel_title, x = NULL, y = "Standardized gene-disjoint score") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.minor = element_blank(), strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 20, hjust = 1),
      plot.title = element_text(face = "bold", size = 14)
    )
}

p_group_melanoma <- make_group_distribution(
  "Melanoma", "Degradation-Low", "Degradation-High",
  "a  Melanoma: 4/4 adjusted node effects are significant (all q < 0.002)"
)
p_group_ccrcc <- make_group_distribution(
  "ccRCC", "Sarcosine-High", "Sarcosine-Low",
  "b  ccRCC: 0/4 adjusted node effects are significant (all q = 0.947)"
)

p_group_effect <- ggplot(route_group_effects, aes(standardized_beta, node_label, color = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0,
                linewidth = 0.85, position = position_dodge(width = 0.42)) +
  geom_point(aes(shape = BH_q < 0.05), size = 3.1, position = position_dodge(width = 0.42)) +
  geom_text(aes(x = CI_high + 0.035, label = q_label), hjust = 0, size = 3.2,
            position = position_dodge(width = 0.42), show.legend = FALSE) +
  scale_color_manual(values = cohort_palette) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1),
                     labels = c(`TRUE` = "BH q < 0.05", `FALSE` = "BH q >= 0.05")) +
  coord_cartesian(xlim = c(-0.38, 1.02), clip = "off") +
  labs(
    title = "c  Adjusted target-state effects",
    subtitle = "Positive values denote Degradation-High or Sarcosine-Low",
    x = "Standardized target-minus-reference effect", y = NULL,
    color = NULL, shape = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "grey35"), legend.position = "bottom",
    plot.margin = margin(5.5, 45, 5.5, 5.5)
  )

route_group_contrast_figure <- ((p_group_melanoma / p_group_ccrcc) | p_group_effect) +
  plot_layout(widths = c(2.15, 1.0)) +
  plot_annotation(
    title = "Sample-level route-score contrasts support melanoma but not ccRCC",
    subtitle = "Within boxplots: red = target state, teal = reference state; forest colors identify cohorts (purple = melanoma, blue = ccRCC)",
    caption = "Interpretation: the shared exact-set enrichment signal does not reproduce as a group-separating sample-level route in ccRCC. Heatmaps remain descriptive views of within-group covariance.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 20),
      plot.subtitle = element_text(size = 12, color = "grey30"),
      plot.caption = element_text(size = 9, color = "grey35", hjust = 0)
    )
  )
ggsave(file.path(figure_dir, "Fig_Route_Score_Group_Contrasts.png"), route_group_contrast_figure,
       width = 18, height = 10.5, dpi = 320, bg = "white")
ggsave(file.path(figure_dir, "Fig_Route_Score_Group_Contrasts.pdf"), route_group_contrast_figure,
       width = 18, height = 10.5, bg = "white")

# -------------------------------------------------------------------------
# Authoritative two-figure reader set
# -------------------------------------------------------------------------
# These two figures are the only standalone reader-facing summary. They make
# the evidence hierarchy explicit so that exact-set GSEA replication cannot be
# mistaken for replication of a sample-level serial mechanism.

selected_serial <- primary_boot[route == leading_route & direction == "Forward"]
selected_noncanonical <- joint_branch_plot[predictor == "Noncanonical_NFkB"]
mel_noncanonical <- selected_noncanonical[cohort == "Melanoma: Degradation-High"]
cc_noncanonical <- selected_noncanonical[cohort == "ccRCC: Sarcosine-Low"]

evidence_ladder <- data.table(
  evidence_order = 1:5,
  evidence_level = c(
    "Exact-set ranked-list GSEA",
    "Gene-disjoint node group effects",
    "Noncanonical branch-IFNG association",
    "Full exposure-linked forward serial model",
    "NFKB2/RELB TF activation"
  ),
  melanoma_text = c(
    "4/4 positive\nall BH q < 0.05",
    "4/4 positive\nall BH q < 0.002",
    paste0("beta ", sprintf("%.3f", mel_noncanonical$standardized_beta),
           "\nq ", format_reader_q(mel_noncanonical$BH_q)),
    "95% CI excludes 0",
    "Not concordant"
  ),
  ccrcc_text = c(
    "4/4 positive\nall BH q < 0.05",
    "0/4 significant\nall q = 0.947",
    paste0("beta ", sprintf("%.3f", cc_noncanonical$standardized_beta),
           "\nq ", format_reader_q(cc_noncanonical$BH_q)),
    "95% CI crosses 0",
    "Not significant"
  ),
  conclusion_text = c(
    "SHARED\nCO-ENRICHMENT",
    "NOT REPLICATED",
    "Direction-consistent\nexploratory association",
    "NOT REPLICATED",
    "COMMON TF ACTIVATION\nNOT SUPPORTED"
  ),
  melanoma_status = c("shared", "shared", "partial", "partial", "not_supported"),
  ccrcc_status = c("shared", "not_supported", "shared", "not_supported", "not_supported"),
  conclusion_status = c("shared", "not_supported", "partial", "not_supported", "not_supported")
)

evidence_ladder_long <- rbindlist(list(
  evidence_ladder[, .(
    evidence_order, evidence_level,
    column = "Melanoma\nDegradation-High", label = melanoma_text, status = melanoma_status
  )],
  evidence_ladder[, .(
    evidence_order, evidence_level,
    column = "ccRCC\nSarcosine-Low", label = ccrcc_text, status = ccrcc_status
  )],
  evidence_ladder[, .(
    evidence_order, evidence_level,
    column = "Cross-cohort\ninterpretation", label = conclusion_text, status = conclusion_status
  )]
))
evidence_ladder_long[, evidence_level := factor(
  evidence_level, levels = rev(evidence_ladder$evidence_level)
)]
evidence_ladder_long[, column := factor(
  column, levels = c("Melanoma\nDegradation-High", "ccRCC\nSarcosine-Low", "Cross-cohort\ninterpretation")
)]

p_reader_ladder <- ggplot(evidence_ladder_long, aes(column, evidence_level)) +
  geom_tile(aes(fill = status), color = "white", linewidth = 1.0) +
  geom_text(aes(label = label), size = 3.5, lineheight = 0.95, color = "grey10") +
  scale_fill_manual(
    values = c(shared = "#74C69D", partial = "#F2CC8F", not_supported = "#D9D9D9"),
    breaks = c("shared", "partial", "not_supported"),
    labels = c("Replicated/shared", "Partial/exploratory", "Not replicated/supported"),
    name = "Evidence status"
  ) +
  labs(
    title = "c  Evidence ladder resolves the apparent contradiction",
    subtitle = "A result can replicate at the gene-set level without replicating as a sample-level mechanism",
    x = NULL, y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(), axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 10), legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey35")
  )

p_reader_program <- ggplot() +
  annotate("rect", xmin = 0.10, xmax = 1.30, ymin = 0.62, ymax = 1.22,
           fill = "#F1E8FA", color = "#51247A", linewidth = 0.9) +
  annotate("rect", xmin = 1.52, xmax = 4.48, ymin = 0.52, ymax = 1.32,
           fill = "#E9F5EE", color = "#1B4332", linewidth = 1.1) +
  annotate("rect", xmin = 4.70, xmax = 5.90, ymin = 0.62, ymax = 1.22,
           fill = "#E7F1FB", color = "#125A93", linewidth = 0.9) +
  geom_segment(aes(x = 1.30, y = 0.92, xend = 1.52, yend = 0.92),
               linewidth = 1.1, color = "#7B2CBF") +
  geom_segment(aes(x = 4.48, y = 0.92, xend = 4.70, yend = 0.92),
               linewidth = 1.1, color = "#1976D2") +
  annotate("text", x = 0.70, y = 0.92,
           label = "Melanoma\nDegradation-High\n4/4 exact sets",
           size = 3.25, fontface = "bold", color = "#51247A", lineheight = 1.05) +
  annotate("text", x = 5.30, y = 0.92,
           label = "ccRCC\nSarcosine-Low\n4/4 exact sets",
           size = 3.25, fontface = "bold", color = "#125A93", lineheight = 1.05) +
  annotate("text", x = 3.00, y = 1.08,
           label = "SHARED INTERSECTION: 4 exact gene sets",
           size = 4.05, fontface = "bold", color = "#1B4332") +
  annotate("text", x = 3.00, y = 0.77,
           label = paste0(
             "APC cross-presentation  +  TCR signaling\n",
             "TNFR2-related noncanonical NF-kappaB  +  IFNG response"
           ),
           size = 3.35, fontface = "bold", color = "#1B4332", lineheight = 1.05) +
  annotate(
    "text", x = 3.0, y = 0.25,
    label = "The intersection is replicated enrichment - not a proven temporal or same-cell signaling chain",
    size = 3.4, color = "grey35"
  ) +
  coord_cartesian(xlim = c(0, 6), ylim = c(0.05, 1.45), clip = "off") +
  labs(
    title = "c  The direct two-cohort intersection",
    subtitle = "The same four exact MSigDB sets pass the target-direction gate independently in both cohorts"
  ) +
  theme_void(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(color = "grey35"),
    plot.margin = margin(10, 18, 10, 18)
  )

# Clean main figure: the selected common candidate axis plus only the evidence
# that is shared between cohorts. Broad screens and failed diagnostic gates are
# deliberately excluded from this figure.
route_candidate_order <- c(
  "APC cross-presentation", "TCR signaling",
  "Noncanonical NF-kappaB", "IFNG response"
)
reader_route <- shared_axis_screen[match(route_candidate_order, candidate)]
stopifnot(identical(reader_route$candidate, route_candidate_order))
reader_route[, `:=`(
  x = c(1.0, 2.7, 4.4, 6.1),
  node_label = c(
    "APC cross-\npresentation", "TCR signaling",
    "TNFR2-related\nnoncanonical NF-kappaB", "IFNG response"
  ),
  melanoma_label = paste0(
    "Melanoma: NES ", sprintf("%.2f", tiger_harmonized_NES),
    " | q ", vapply(tiger_q, format_reader_q, character(1))
  ),
  ccrcc_label = paste0(
    "ccRCC: NES ", sprintf("%.2f", ccrcc_harmonized_NES),
    " | q ", vapply(ccrcc_q, format_reader_q, character(1))
  ),
  overlap_label = paste0(
    "Shared leading edge: n=", shared_leading_edge_n,
    " | Jaccard ", sprintf("%.1f%%", 100 * leading_edge_jaccard)
  )
)]

reader_gene_examples <- list(
  "APC cross-presentation" = c("B2M", "HLA-A", "HLA-B", "HLA-C", "PSMB8", "PSMB9", "TAP1", "TAPBP"),
  "TCR signaling" = c("CARD11", "ITK", "LCP2", "NFKB1", "NFKBIA", "PTPN22", "HLA-DRA", "HLA-DRB1"),
  "Noncanonical NF-kappaB" = c("NFKB2", "RELB", "TRAF3"),
  "IFNG response" = c("CXCL9", "CXCL10", "CXCL11", "IRF1", "IRF8", "STAT1", "NFKB1", "CD86", "TAP1", "TAPBP")
)
reader_route[, gene_examples := vapply(candidate, function(candidate_name) {
  shared_genes <- strsplit(shared_leading_edge_genes[candidate == candidate_name], ";", fixed = TRUE)[[1]]
  examples <- reader_gene_examples[[candidate_name]]
  if (!all(examples %in% shared_genes)) stop("Reader gene example is not in the exact shared leading edge: ", candidate_name)
  prefix <- if (candidate_name == "Noncanonical NF-kappaB") "All shared leading-edge genes: " else "Shared leading-edge examples: "
  paste(strwrap(paste0(prefix, paste(examples, collapse = ", ")), width = 35), collapse = "\n")
}, character(1))]

reader_route[, `:=`(
  node_xmin = x - c(0.57, 0.55, 0.76, 0.55),
  node_xmax = x + c(0.57, 0.55, 0.76, 0.55),
  card_xmin = x - 0.78,
  card_xmax = x + 0.78
)]
reader_edges <- data.table(
  x = reader_route$node_xmax[-nrow(reader_route)] + 0.05,
  xend = reader_route$node_xmin[-1] - 0.05,
  y = 3.43, yend = 3.43
)

p_shared_route_main <- ggplot() +
  annotate("label", x = 3.55, y = 4.28,
           label = "REPLICATED IN BOTH COHORTS: 4/4 EXACT PATHWAY SETS",
           size = 4.7, fontface = "bold", fill = "#CDEED8", color = "#1B4332",
           label.padding = unit(0.28, "lines"), linewidth = 0.9) +
  geom_segment(
    data = reader_edges, aes(x, y, xend = xend, yend = yend),
    color = "#C73E3A", linewidth = 1.1, linetype = 2,
    arrow = arrow(length = unit(0.15, "inches"), type = "closed")
  ) +
  geom_rect(
    data = reader_route,
    aes(xmin = node_xmin, xmax = node_xmax, ymin = 3.08, ymax = 3.78),
    fill = "#E9F5EE", color = "#1B4332", linewidth = 0.9
  ) +
  geom_text(data = reader_route, aes(x, y = 3.43, label = node_label),
            size = 4.2, fontface = "bold", color = "#173F35", lineheight = 1.02) +
  geom_rect(
    data = reader_route,
    aes(xmin = card_xmin, xmax = card_xmax, ymin = 0.48, ymax = 2.55),
    fill = "white", color = "#74A892", linewidth = 0.8
  ) +
  geom_text(data = reader_route, aes(x, y = 2.27),
            label = "BOTH COHORTS: positive NES and BH q < 0.05",
            size = 3.15, fontface = "bold", color = "#1B4332") +
  geom_text(data = reader_route, aes(x, y = 1.89, label = melanoma_label),
            size = 3.05, color = "#7B2CBF", fontface = "bold") +
  geom_text(data = reader_route, aes(x, y = 1.55, label = ccrcc_label),
            size = 3.05, color = "#1976D2", fontface = "bold") +
  geom_text(data = reader_route, aes(x, y = 1.16, label = overlap_label),
            size = 3.05, color = "grey15", fontface = "bold") +
  geom_text(data = reader_route, aes(x, y = 0.76, label = gene_examples),
            size = 2.55, color = "grey25", lineheight = 1.0) +
  annotate("text", x = 3.55, y = 2.82,
           label = "Dashed arrows show the selected candidate ordering; shared enrichment and leading-edge overlap support the axis, not temporal causality.",
           size = 3.35, color = "grey35") +
  coord_cartesian(xlim = c(0.10, 7.00), ylim = c(0.20, 4.62), clip = "off") +
  labs(
    title = "Selected immune axis shared by melanoma and ccRCC",
    subtitle = "Target states: melanoma Sarcosine-degradation High and ccRCC measured Sarcosine-Low",
    caption = "Each node is the same exact MSigDB pathway in both cohorts. Gene examples are drawn only from the exact cross-cohort leading-edge intersection."
  ) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 23),
    plot.subtitle = element_text(size = 13, color = "grey25"),
    plot.caption = element_text(size = 9, color = "grey40", hjust = 0),
    plot.margin = margin(25, 35, 25, 35)
  )

reader_shared_figure <- p_shared_route_main

ggsave(file.path(reader_figure_dir, "Fig_01_What_Is_Shared_Across_Cohorts.png"), reader_shared_figure,
       width = 18, height = 8.8, dpi = 320, bg = "white")
ggsave(file.path(reader_figure_dir, "Fig_01_What_Is_Shared_Across_Cohorts.pdf"), reader_shared_figure,
       width = 18, height = 8.8, bg = "white")

# Replace the previously broad/diagnostic overview at the exact legacy path
# referenced by the reader. The filename is retained for compatibility, but
# its content now contains only the cross-cohort overlap result.
ggsave(file.path(figure_dir, "Fig_Common_Immune_Mechanism_Broad_to_Final_Narrowing.png"), reader_shared_figure,
       width = 18, height = 8.8, dpi = 320, bg = "white")
ggsave(file.path(figure_dir, "Fig_Common_Immune_Mechanism_Broad_to_Final_Narrowing.pdf"), reader_shared_figure,
       width = 18, height = 8.8, bg = "white")
ggsave(file.path(figure_dir, "Fig_Leading_Common_Immune_Candidate_Route.png"), reader_shared_figure,
       width = 18, height = 8.8, dpi = 320, bg = "white")
ggsave(file.path(figure_dir, "Fig_Leading_Common_Immune_Candidate_Route.pdf"), reader_shared_figure,
       width = 18, height = 8.8, bg = "white")

shared_overlap_support_figure <- (p_shared_a | p_shared_b) +
  plot_annotation(
    title = "Cross-cohort overlap supporting the selected immune axis",
    subtitle = "The same four exact pathways are enriched in both target states and share cohort-specific leading-edge genes",
    caption = "These are independent cohort tests of the same MSigDB sets; shared leading-edge genes are their exact intersections.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 21),
      plot.subtitle = element_text(size = 12.5, color = "grey25"),
      plot.caption = element_text(size = 9, color = "grey40", hjust = 0)
    )
  )
ggsave(file.path(figure_dir, "Fig_Cross_Cohort_Evidence_for_Selected_Immune_Axis.png"), shared_overlap_support_figure,
       width = 18, height = 7.0, dpi = 320, bg = "white")
ggsave(file.path(figure_dir, "Fig_Cross_Cohort_Evidence_for_Selected_Immune_Axis.pdf"), shared_overlap_support_figure,
       width = 18, height = 7.0, bg = "white")

reader_group_contrast_figure <- ((p_group_melanoma / p_group_ccrcc) | p_group_effect) +
  plot_layout(widths = c(2.15, 1.0)) +
  plot_annotation(
    title = "Boundary of the shared four-program intersection",
    subtitle = paste0(
      "The four exact gene sets remain the cross-cohort shared result. Within boxplots: red = target state and teal = reference state.  ",
      "Forest: purple = melanoma and blue = ccRCC."
    ),
    caption = paste0(
      "Melanoma shows separation for all four gene-disjoint scores; ccRCC shows none (all BH q = 0.947). ",
      "These results do not negate exact-set GSEA replication, but they preclude a shared serial-mechanism claim."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 20),
      plot.subtitle = element_text(size = 12, color = "grey30"),
      plot.caption = element_text(size = 9, color = "grey35", hjust = 0)
    )
  )

ggsave(file.path(reader_figure_dir, "Fig_02_Cohort_Specific_Group_Contrasts.png"), reader_group_contrast_figure,
       width = 18, height = 10.5, dpi = 320, bg = "white")
ggsave(file.path(reader_figure_dir, "Fig_02_Cohort_Specific_Group_Contrasts.pdf"), reader_group_contrast_figure,
       width = 18, height = 10.5, bg = "white")

fwrite(evidence_ladder, file.path(table_dir, "16_reader_evidence_ladder.csv"))
fwrite(data.table(
  order = 1:2,
  figure = c("Fig_01_What_Is_Shared_Across_Cohorts", "Fig_02_Cohort_Specific_Group_Contrasts"),
  purpose = c(
    "Presents the selected four-node common candidate axis with paired cohort NES/q values and exact leading-edge intersections",
    "Shows why the gene-disjoint sample-level route supports melanoma but not ccRCC"
  ),
  standalone_claim = c(
    "Shared co-enriched immune program; no shared serial mechanism established",
    "Sample-level group separation replicates in melanoma only"
  )
), file.path(reader_figure_dir, "FIGURE_MANIFEST.csv"))

writeLines(c(
  "# Authoritative reader-facing cross-cohort figure set",
  "",
  "Use these figures in order:",
  "",
  "1. `Fig_01_What_Is_Shared_Across_Cohorts`: the four exact MSigDB gene sets are enriched in the target direction in both cohorts. This is the shared result.",
  "2. `Fig_02_Cohort_Specific_Group_Contrasts`: the gene-disjoint sample-level scores separate melanoma groups but not ccRCC groups. Therefore a shared exposure-linked serial mechanism is not claimed.",
  "",
  "Terminology:",
  "",
  "- Supported: **shared co-enriched immune program**.",
  "- Not supported: **shared causal/serial pathway**, **same-cell signaling chain**, or **common NFKB2/RELB TF activation**.",
  "",
  "Color contract:",
  "",
  "- Within-cohort group plots: red = target state; teal = reference state.",
  "- Cross-cohort comparison plots: purple = melanoma; blue = ccRCC.",
  "",
  "The files in `../figures` are extended diagnostics/provenance and should not be presented as a standalone figure sequence."
), file.path(reader_figure_dir, "README.md"), useBytes = TRUE)

# Descriptive heatmap of the final-route nodes, separately ordered by group.
heatmap_one <- function(data, cohort_label, target_label, reference_label, path) {
  branch <- leading_catalog$branch
  mat <- t(as.matrix(data[, .(APC, TCR, branch_value = get(branch), IFNG_response)]))
  rownames(mat) <- c("APC cross-presentation", "TCR signaling",
                    if (branch == "Canonical_NFkB") "Canonical NF-kappaB" else "TNFR2-related noncanonical NF-kappaB",
                    "IFNG response")
  colnames(mat) <- data$sample_id
  ord <- order(data$target_group, data$IFNG_response)
  mat <- mat[, ord, drop = FALSE]
  ann <- data.frame("Exposure group" = ifelse(data$target_group[ord] == "Target", target_label, reference_label),
                    row.names = data$sample_id[ord], check.names = FALSE)
  group_colors <- setNames(c("#C73E3A", "#2A9D8F"), c(target_label, reference_label))
  pheatmap(
    mat, cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE,
    color = colorRampPalette(c("#2A9D8F", "white", "#C73E3A"))(101),
    breaks = seq(-2.5, 2.5, length.out = 102), annotation_col = ann,
    annotation_colors = setNames(list(group_colors), "Exposure group"),
    main = cohort_label, fontsize_row = 11, border_color = NA,
    filename = path, width = 12, height = 4.8
  )
}
heatmap_one(tiger_data, "Melanoma: descriptive route-score heatmap\nOrdered by IFNG within group; group effects are shown separately",
            "Degradation-High", "Degradation-Low",
            file.path(figure_dir, "Fig_Leading_Common_Route_Heatmap_Melanoma.png"))
heatmap_one(cc_data, "ccRCC: descriptive route-score heatmap\nOrdered by IFNG within group; group effects are shown separately",
            "Sarcosine-Low", "Sarcosine-High",
            file.path(figure_dir, "Fig_Leading_Common_Route_Heatmap_ccRCC.png"))

# -------------------------------------------------------------------------
# Outputs and report
# -------------------------------------------------------------------------
sample_scores <- rbindlist(list(
  tiger_data[, .(cohort = "Melanoma: Degradation-High", sample_id, target_group,
                 target_continuous_z, APC, TCR, CD28_family, Canonical_NFkB,
                 Noncanonical_NFkB, IFNG_response, SPI1_activity, NFKB1_activity,
                 RELA_activity, NFKB1_RELA_activity, NFKB2_activity, NFKB2_expression,
                 RELB_expression, T_cell_context_z, APC_myeloid_context_z)],
  cc_data[, .(cohort = "ccRCC: Sarcosine-Low", sample_id, target_group,
              target_continuous_z, APC, TCR, CD28_family, Canonical_NFkB,
              Noncanonical_NFkB, IFNG_response, SPI1_activity, NFKB1_activity,
              RELA_activity, NFKB1_RELA_activity, NFKB2_activity, NFKB2_expression,
              RELB_expression, T_cell_context_z, APC_myeloid_context_z)]
), use.names = TRUE)

correlation_variables <- c("APC", "TCR", "CD28_family", "Canonical_NFkB", "Noncanonical_NFkB", "IFNG_response")
correlation_audit <- rbindlist(lapply(unique(sample_scores$cohort), function(cohort_name) {
  mat <- cor(as.matrix(sample_scores[cohort == cohort_name, ..correlation_variables]))
  out <- as.data.table(as.table(mat))
  setnames(out, c("score_1", "score_2", "Pearson_r"))
  out[, cohort := cohort_name]
  out
}))
vif_audit <- rbindlist(lapply(unique(sample_scores$cohort), function(cohort_name) {
  d <- sample_scores[cohort == cohort_name]
  predictors <- setdiff(correlation_variables, "IFNG_response")
  rbindlist(lapply(predictors, function(term) {
    others <- setdiff(predictors, term)
    fit <- lm(as.formula(paste(term, "~", paste(others, collapse = " + "))), data = d)
    data.table(cohort = cohort_name, predictor = term, VIF = 1 / (1 - summary(fit)$r.squared))
  }))
}))

input_manifest <- data.table(
  input_role = c(
    "ccRCC expression", "ccRCC metadata", "ccRCC named-cell deconvolution", "ccRCC signed TF regulon",
    "Prior exact common-pathway table", "TIGER expression", "TIGER PRE73 metadata",
    "TIGER signed TF regulon", "TIGER named-cell context"
  ),
  path = input_paths,
  SHA256 = vapply(input_paths, sha256_file, character(1))
)

fwrite(input_manifest, file.path(table_dir, "01_input_manifest_sha256.csv"))
fwrite(screen, file.path(table_dir, "02_broad_candidate_pathway_screen.csv"))
fwrite(score_manifest, file.path(table_dir, "03_gene_disjoint_score_manifest.csv"))
fwrite(overlap_audit, file.path(table_dir, "04_score_overlap_audit.csv"))
fwrite(sample_scores, file.path(table_dir, "05_sample_level_scores.csv"))
fwrite(node_effects, file.path(table_dir, "06_node_exposure_effects.csv"))
fwrite(route_edges, file.path(table_dir, "07_route_edge_models.csv"))
fwrite(route_bootstrap, file.path(table_dir, "08_route_serial_bootstrap.csv"))
fwrite(joint_results, file.path(table_dir, "09_joint_branch_competition.csv"))
fwrite(branch_meta, file.path(table_dir, "09b_two_cohort_branch_fixed_effect_meta.csv"))
fwrite(tf_coverage, file.path(table_dir, "10_consensus_TF_regulon_coverage.csv"))
fwrite(tf_support, file.path(table_dir, "11_TF_support_models.csv"))
fwrite(route_evidence, file.path(table_dir, "12_route_evidence_summary.csv"))
fwrite(final_decision, file.path(table_dir, "13_final_route_decision.csv"))
fwrite(correlation_audit, file.path(table_dir, "14_score_correlation_audit.csv"))
fwrite(vif_audit, file.path(table_dir, "15_score_VIF_audit.csv"))

lead_primary <- primary_boot[route == leading_route & direction == "Forward"]
lead_reverse <- primary_boot[route == leading_route & direction == "Reverse"]
lead_joint <- primary_joint[predictor == leading_catalog$branch]
lead_meta <- branch_meta[branch == leading_catalog$branch]
canonical_meta <- branch_meta[branch == "Canonical_NFkB"]
canonical_screen <- get_screen("Canonical NF-kappaB")
noncanonical_screen <- get_screen("Noncanonical NF-kappaB")
cd28_specific_screen <- get_screen("CD28 co-stimulation (specific)")
il12_screen <- get_screen("IL-12/STAT4")
nfkb2_coverage <- tf_coverage[TF == "NFKB2"]
relb_coverage <- tf_coverage[TF == "RELB"]
nfkb2_group <- tf_effects[model == "Grouped primary" & outcome == "NFKB2_activity"]
nfkb2_link <- tf_links[edge == "NFKB2 activity -> noncanonical NF-kappaB program"]
nfkb2_transcript_group <- tf_transcript_effects[model == "Grouped primary" & outcome == "NFKB2_expression"]
relb_transcript_group <- tf_transcript_effects[model == "Grouped primary" & outcome == "RELB_expression"]

fmt <- function(x, digits = 3) formatC(x, digits = digits, format = "fg")
report_lines <- c(
  "# Cross-cohort common immune-program final narrowing",
  "",
  "## Bottom line",
  "",
  "The reproducible cross-cohort result is a **shared four-program enrichment intersection**: APC cross-presentation, TCR signaling, the TNFR2-related noncanonical NF-kappaB expression program, and IFNG response are each enriched in the target direction in both cohorts.",
  paste0("Within that shared program, **", gsub(" -> ", " → ", leading_route), "** is a direction-consistent exploratory branch association, not a replicated serial mechanism."),
  "No complete exposure-linked serial route was reproduced in both cohorts: the first target-state-to-APC link was absent in ccRCC, and common NFKB2/RELB TF activation was not demonstrated.",
  "",
  "## Target-state definitions",
  "",
  "- TIGER melanoma: sarcosine-degradation High versus Low among 73 PRE-treatment tumors.",
  "- ccRCC: measured tumor Sarcosine-Low versus High by the prespecified median split among 100 matched tumors.",
  "- Positive harmonized NES and positive model effects denote those target states.",
  "",
  "## Broad screen",
  "",
  paste0("- Exact common-pathway gates were evaluated for ", nrow(screen), " prespecified candidates."),
  paste0("- ", sum(screen$exact_replication_pass), " candidates were direction-aligned and BH q<0.05 in both cohorts."),
  paste0("- The displayed shared downstream axis passed the exact-set gate at all four nodes: APC cross-presentation (shared leading edge n=", get_screen("APC cross-presentation")$shared_leading_edge_n,
         "), TCR signaling (n=", get_screen("TCR signaling")$shared_leading_edge_n,
         "), noncanonical NF-kappaB (n=", noncanonical_screen$shared_leading_edge_n,
         "; NFKB2, RELB, TRAF3), and IFNG response (n=", get_screen("IFNG response")$shared_leading_edge_n, ")."),
  paste0("- Canonical NF-kappaB: TIGER NES ", fmt(canonical_screen$tiger_harmonized_NES), ", q ", fmt(canonical_screen$tiger_q),
         "; ccRCC NES ", fmt(canonical_screen$ccrcc_harmonized_NES), ", q ", fmt(canonical_screen$ccrcc_q),
         "; shared leading-edge genes ", canonical_screen$shared_leading_edge_n, "."),
  paste0("- Noncanonical NF-kappaB: TIGER NES ", fmt(noncanonical_screen$tiger_harmonized_NES), ", q ", fmt(noncanonical_screen$tiger_q),
         "; ccRCC NES ", fmt(noncanonical_screen$ccrcc_harmonized_NES), ", q ", fmt(noncanonical_screen$ccrcc_q),
         "; shared leading-edge genes ", noncanonical_screen$shared_leading_edge_n, "."),
  paste0("- Specific CD28 co-stimulation failed the two-cohort gate (TIGER q ", fmt(cd28_specific_screen$tiger_q),
         "; ccRCC q ", fmt(cd28_specific_screen$ccrcc_q), "). The broad CD28-family set was retained only as a mixed regulatory sensitivity node."),
  paste0("- IL-12/STAT4 failed the two-cohort gate (TIGER q ", fmt(il12_screen$tiger_q),
         "; ccRCC q ", fmt(il12_screen$ccrcc_q), ")."),
  "",
  "## Gene-disjoint route design",
  "",
  "APC cross-presentation, TCR, broad CD28-family regulation, canonical NF-kappaB, noncanonical NF-kappaB, and the IFNG-response outcome were scored from MSigDB 2026.1.Hs. IFNG-response genes were removed from every predictor, and genes shared by two or more predictor modules were removed from every predictor. All retained predictor/outcome scores therefore have zero pairwise gene overlap.",
  "",
  "## Direct sample-level group separation",
  "",
  paste0("- Melanoma gene-disjoint route scores all showed positive adjusted Degradation-High effects: APC beta ",
         fmt(route_group_effects[cohort == "Melanoma: Degradation-High" & outcome == "APC", standardized_beta]),
         ", TCR beta ", fmt(route_group_effects[cohort == "Melanoma: Degradation-High" & outcome == "TCR", standardized_beta]),
         ", noncanonical NF-kappaB beta ", fmt(route_group_effects[cohort == "Melanoma: Degradation-High" & outcome == "Noncanonical_NFkB", standardized_beta]),
         ", and IFNG-response beta ", fmt(route_group_effects[cohort == "Melanoma: Degradation-High" & outcome == "IFNG_response", standardized_beta]),
         "; all BH q<0.002."),
  paste0("- ccRCC gene-disjoint route scores did not separate Sarcosine-Low from High: APC beta ",
         fmt(route_group_effects[cohort == "ccRCC: Sarcosine-Low" & outcome == "APC", standardized_beta]),
         ", TCR beta ", fmt(route_group_effects[cohort == "ccRCC: Sarcosine-Low" & outcome == "TCR", standardized_beta]),
         ", noncanonical NF-kappaB beta ", fmt(route_group_effects[cohort == "ccRCC: Sarcosine-Low" & outcome == "Noncanonical_NFkB", standardized_beta]),
         ", and IFNG-response beta ", fmt(route_group_effects[cohort == "ccRCC: Sarcosine-Low" & outcome == "IFNG_response", standardized_beta]),
         "; all BH q=0.947."),
  "This explains the weak visual group separation in the ccRCC heatmap. The cross-cohort overlap is restricted to full exact-set ranked-list GSEA and downstream branch-level covariance; it is not a replicated ccRCC sample-level route contrast.",
  "",
  "## Primary grouped serial products for the leading route",
  "",
  paste0("- Melanoma: ", fmt(lead_primary[cohort == "Melanoma: Degradation-High", serial_product]), " [",
         fmt(lead_primary[cohort == "Melanoma: Degradation-High", bootstrap_CI_low]), ", ",
         fmt(lead_primary[cohort == "Melanoma: Degradation-High", bootstrap_CI_high]), "]."),
  paste0("- ccRCC: ", fmt(lead_primary[cohort == "ccRCC: Sarcosine-Low", serial_product]), " [",
         fmt(lead_primary[cohort == "ccRCC: Sarcosine-Low", bootstrap_CI_low]), ", ",
         fmt(lead_primary[cohort == "ccRCC: Sarcosine-Low", bootstrap_CI_high]), "]."),
  paste0("- Complete forward serial replication in both cohorts: **", ifelse(complete_serial_replicated, "PASS", "FAIL"), "**."),
  "- Reverse-order products are retained in the source table and are treated as an ambiguity diagnostic, not as proof of reverse biology.",
  "",
  "## Branch competition",
  "",
  paste0("- Melanoma leading-branch beta: ", fmt(lead_joint[cohort == "Melanoma: Degradation-High", standardized_beta]),
         ", BH q ", fmt(lead_joint[cohort == "Melanoma: Degradation-High", BH_q]), "."),
  paste0("- ccRCC leading-branch beta: ", fmt(lead_joint[cohort == "ccRCC: Sarcosine-Low", standardized_beta]),
         ", BH q ", fmt(lead_joint[cohort == "ccRCC: Sarcosine-Low", BH_q]), "."),
  paste0("- Two-cohort fixed-effect branch estimate: beta ", fmt(lead_meta$fixed_effect_beta), " [",
         fmt(lead_meta$fixed_effect_CI_low), ", ", fmt(lead_meta$fixed_effect_CI_high),
         "], p ", fmt(lead_meta$fixed_effect_p), ", I2 ", fmt(lead_meta$I2_percent), "%."),
  paste0("- Canonical NF-kappaB was direction-discordant after joint competition; its fixed-effect beta was ",
         fmt(canonical_meta$fixed_effect_beta), ", p ", fmt(canonical_meta$fixed_effect_p),
         ", I2 ", fmt(canonical_meta$I2_percent), "%. This is why the strong canonical enrichment was not promoted to the final common branch."),
  "",
  "## TF evidence",
  "",
  paste0("SPI1, NFKB1, RELA, and NFKB2 activities were inferred from cross-cohort-consensus signed DoRothEA A/B target edges after excluding every route-score and IFNG-response gene. NFKB2 retained ",
         nfkb2_coverage$retained_targets_n, " targets and met the prespecified minimum of 10; RELB retained only ",
         relb_coverage$retained_targets_n, " targets and its activity was therefore not estimated."),
  paste0("- NFKB2 inferred-activity target-state effect: melanoma beta ",
         fmt(nfkb2_group[cohort == "Melanoma: Degradation-High", standardized_beta]), ", q ",
         fmt(nfkb2_group[cohort == "Melanoma: Degradation-High", BH_q]), "; ccRCC beta ",
         fmt(nfkb2_group[cohort == "ccRCC: Sarcosine-Low", standardized_beta]), ", q ",
         fmt(nfkb2_group[cohort == "ccRCC: Sarcosine-Low", BH_q]), "."),
  paste0("- NFKB2 activity-to-noncanonical-program association: melanoma beta ",
         fmt(nfkb2_link[cohort == "Melanoma: Degradation-High", standardized_beta]), ", q ",
         fmt(nfkb2_link[cohort == "Melanoma: Degradation-High", BH_q]), "; ccRCC beta ",
         fmt(nfkb2_link[cohort == "ccRCC: Sarcosine-Low", standardized_beta]), ", q ",
         fmt(nfkb2_link[cohort == "ccRCC: Sarcosine-Low", BH_q]), "."),
  paste0("- NFKB2 transcript target-state effect: melanoma beta ",
         fmt(nfkb2_transcript_group[cohort == "Melanoma: Degradation-High", standardized_beta]), ", q ",
         fmt(nfkb2_transcript_group[cohort == "Melanoma: Degradation-High", BH_q]), "; ccRCC beta ",
         fmt(nfkb2_transcript_group[cohort == "ccRCC: Sarcosine-Low", standardized_beta]), ", q ",
         fmt(nfkb2_transcript_group[cohort == "ccRCC: Sarcosine-Low", BH_q]), "."),
  paste0("- RELB transcript target-state effect: melanoma beta ",
         fmt(relb_transcript_group[cohort == "Melanoma: Degradation-High", standardized_beta]), ", q ",
         fmt(relb_transcript_group[cohort == "Melanoma: Degradation-High", BH_q]), "; ccRCC beta ",
         fmt(relb_transcript_group[cohort == "ccRCC: Sarcosine-Low", standardized_beta]), ", q ",
         fmt(relb_transcript_group[cohort == "ccRCC: Sarcosine-Low", BH_q]), "."),
  "SPI1 remains an APC-context check, while NFKB1/RELA remain the prespecified canonical-branch check. TF activities and focal transcripts are corroborative diagnostics and are not inserted as causal nodes.",
  "",
  "## Interpretation boundary",
  "",
  "The two exposures are not identical measurements: TIGER uses a SARDH/PIPOX transcriptional degradation score, whereas ccRCC uses measured tumor Sarcosine. Bulk RNA cannot identify ligand secretion, receptor engagement, phosphorylation, cell of origin, or temporal ordering. Named T-cell/APC context-adjusted models are sensitivity analyses because adjustment can remove either confounding or true immune-composition mediation.",
  "",
  "## Reader-facing outputs",
  "",
  "Use only the following two figures as the connected standalone reader sequence:",
  "",
  "1. `figures_final_reader/Fig_01_What_Is_Shared_Across_Cohorts.png`: directly presents the selected common candidate axis with paired cohort NES/q values and shared leading-edge evidence at every node.",
  "2. `figures_final_reader/Fig_02_Cohort_Specific_Group_Contrasts.png`: shows why the shared enrichment cannot be promoted to a replicated sample-level serial mechanism.",
  "",
  "The legacy `figures/` directory contains extended diagnostics and provenance. It is not a standalone narrative sequence. Full clinical, continuous, reverse-order, and named-cell-context results remain in the tables."
)
writeLines(report_lines, file.path(output_root, "REPORT.md"), useBytes = TRUE)

capture.output(sessionInfo(), file = file.path(output_root, "sessionInfo.txt"))
writeLines(c(
  paste0("Script SHA256: ", sha256_file(script_path)),
  paste0("MSigDB version: ", msigdb_version),
  paste0("Bootstrap replicates per route/direction/model/cohort: 5000"),
  paste0("Leading route: ", leading_route)
), file.path(log_dir, "run_summary.txt"))

cat("Completed cross-cohort mechanism narrowing.\n")
cat("Leading route:", leading_route, "\n")
cat("Output:", output_root, "\n")
