#!/usr/bin/env Rscript

## TCGA NSCLC validation of the prespecified shared immune axis.
## NSCLC = TCGA LUAD + LUSC primary tumors.
## Exposure definition exactly follows the melanoma analysis:
## degradation score = mean of within-cohort SARDH and PIPOX z scores;
## High if score > within-cohort median, otherwise Low.
##
## Primary inference:
##   limma group contrast (High - Low), adjusted for histology, age and sex;
##   exact-set preranked GSEA for four pathways frozen before this analysis.
## Supporting inference:
##   continuous degradation-score GSEA;
##   gene-disjoint sample-level node scores and adjacent-node associations.

options(stringsAsFactors = FALSE, warn = 1)
set.seed(20260903)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

local_lib <- Sys.getenv("SARCO_GSEA_R_LIB", unset = "")
if (nzchar(local_lib)) {
  if (!dir.exists(local_lib)) stop("SARCO_GSEA_R_LIB does not exist: ", local_lib)
  .libPaths(c(normalizePath(local_lib), .libPaths()))
}

suppressPackageStartupMessages({
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
})

dir.create(file.path(project_dir, "results", "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_dir, "results", "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(project_dir, "results", "logs"), recursive = TRUE, showWarnings = FALSE)
table_dir <- file.path(project_dir, "results", "tables")
figure_dir <- file.path(project_dir, "results", "figures")
log_dir <- file.path(project_dir, "results", "logs")

input_raw <- file.path(project_dir, "data", "01_xena_raw.rds")
input_scores <- file.path(project_dir, "data", "02_scores.rds")
input_strict_manifest <- file.path(project_dir, "config", "gene_disjoint_score_manifest.csv")
stopifnot(file.exists(input_raw), file.exists(input_scores), file.exists(input_strict_manifest))

raw <- readRDS(input_raw)
precomputed_scores <- readRDS(input_scores)
strict_manifest <- read.csv(input_strict_manifest, check.names = FALSE)

pathway_meta <- data.frame(
  node_order = 1:4,
  node = c(
    "APC cross-presentation",
    "TCR signaling",
    "TNFR2-related\nnoncanonical NF-kappaB",
    "IFNG response"
  ),
  short_node = c("APC", "TCR", "Noncanonical NF-kappaB", "IFNG"),
  pathway = c(
    "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION",
    "REACTOME_TCR_SIGNALING",
    "REACTOME_TNFR2_NON_CANONICAL_NF_KB_PATHWAY",
    "HALLMARK_INTERFERON_GAMMA_RESPONSE"
  ),
  collection = c("Reactome", "Reactome", "Reactome", "Hallmark"),
  stringsAsFactors = FALSE
)

cat("Loading local MSigDB 2026.1.Hs gene sets...\n")
msig_h <- suppressMessages(msigdbr(
  species = "Homo sapiens", db_species = "HS", collection = "H"
))
msig_r <- suppressMessages(msigdbr(
  species = "Homo sapiens", db_species = "HS",
  collection = "C2", subcollection = "CP:REACTOME"
))
msig_version <- unique(c(msig_h$db_version, msig_r$db_version))
if (length(msig_version) != 1L || msig_version != "2026.1.Hs") {
  stop("Expected one local MSigDB version (2026.1.Hs); found: ",
       paste(msig_version, collapse = ", "))
}

split_unique <- function(gene, set) {
  out <- split(gene, set)
  lapply(out, unique)
}
hallmark_sets <- split_unique(msig_h$gene_symbol, msig_h$gs_name)
reactome_sets <- split_unique(msig_r$gene_symbol, msig_r$gs_name)
exact_sets <- lapply(pathway_meta$pathway, function(pw) {
  if (pw %in% names(hallmark_sets)) hallmark_sets[[pw]] else reactome_sets[[pw]]
})
names(exact_sets) <- pathway_meta$pathway
if (any(lengths(exact_sets) == 0L)) stop("At least one prespecified pathway was not found.")
exact_score_sets <- exact_sets
names(exact_score_sets) <- c("APC", "TCR", "Noncanonical_NFkB", "IFNG")

strict_wanted <- c(
  "APC" = "APC",
  "TCR" = "TCR",
  "Noncanonical_NFkB" = "Noncanonical_NFkB",
  "IFNG" = "IFNG_response"
)
strict_display <- c(
  APC = "APC cross-presentation",
  TCR = "TCR signaling",
  Noncanonical_NFkB = "TNFR2-related\nnoncanonical NF-kappaB",
  IFNG = "IFNG response"
)
strict_sets <- lapply(strict_wanted, function(key) {
  row <- strict_manifest[strict_manifest$score == key, , drop = FALSE]
  if (nrow(row) != 1L) stop("Strict score definition not unique for: ", key)
  unique(strsplit(row$strict_genes, ";", fixed = TRUE)[[1]])
})

fmt_q <- function(x) {
  ifelse(is.na(x), "NA",
         ifelse(x < 1e-3, formatC(x, format = "e", digits = 1),
                formatC(x, format = "f", digits = 3)))
}

as_numeric_clean <- function(x) suppressWarnings(as.numeric(as.character(x)))

prepare_cohort <- function(cohort) {
  expr_df <- raw[[cohort]]$expr
  genes <- as.character(expr_df[[1]])
  expr <- as.matrix(expr_df[, -1, drop = FALSE])
  storage.mode(expr) <- "double"
  rownames(expr) <- genes
  expr <- expr[!duplicated(rownames(expr)), , drop = FALSE]

  score_df <- precomputed_scores[[cohort]]
  sample_idx <- match(score_df$sample, colnames(expr))
  if (anyNA(sample_idx)) stop(cohort, ": score samples missing from expression.")
  expr <- expr[, sample_idx, drop = FALSE]
  if (!identical(colnames(expr), score_df$sample)) stop(cohort, ": sample-order mismatch.")

  clin <- raw[[cohort]]$clin
  clin_idx <- match(score_df$sample, clin$sampleID)
  if (anyNA(clin_idx)) stop(cohort, ": samples missing from clinical metadata.")
  clin <- clin[clin_idx, , drop = FALSE]

  degradation_recomputed <- rowMeans(cbind(
    as.numeric(scale(score_df$SARDH)),
    as.numeric(scale(score_df$PIPOX))
  ))
  max_delta <- max(abs(degradation_recomputed - score_df$degradation_z), na.rm = TRUE)
  if (max_delta > 1e-10) stop(cohort, ": degradation-score reconstruction mismatch.")

  cutoff <- median(score_df$degradation_z, na.rm = TRUE)
  group <- ifelse(score_df$degradation_z > cutoff, "High", "Low")
  group <- factor(group, levels = c("Low", "High"))

  disease <- as.character(clin$`_primary_disease`)
  subtype <- ifelse(grepl("adenocarcinoma", disease, ignore.case = TRUE), "LUAD",
                    ifelse(grepl("squamous", disease, ignore.case = TRUE), "LUSC", NA_character_))
  age <- as_numeric_clean(clin$age_at_initial_pathologic_diagnosis)
  sex <- toupper(trimws(as.character(clin$gender)))
  sex[!sex %in% c("FEMALE", "MALE")] <- NA_character_

  meta <- data.frame(
    sample = score_df$sample,
    cohort = cohort,
    subtype = factor(subtype),
    age = age,
    age_z = as.numeric(scale(age)),
    sex = factor(sex),
    SARDH = score_df$SARDH,
    PIPOX = score_df$PIPOX,
    degradation_z = score_df$degradation_z,
    cutoff = cutoff,
    degradation_group = group,
    stringsAsFactors = FALSE
  )
  complete <- complete.cases(meta[, c("subtype", "age_z", "sex", "degradation_group")])
  if (sum(complete) < 0.9 * nrow(meta)) stop(cohort, ": >10% excluded by primary covariates.")
  if (any(table(meta$degradation_group[complete]) < 20L)) stop(cohort, ": small primary group.")

  list(expr = expr, meta = meta, complete = complete, max_delta = max_delta)
}

if (!identical(names(raw), "NSCLC") || !identical(names(precomputed_scores), "NSCLC")) {
  stop("The public validation input must contain the NSCLC cohort only.")
}
prepared <- list(NSCLC = prepare_cohort("NSCLC"))

cohort_summary <- do.call(rbind, lapply(names(prepared), function(cohort) {
  z <- prepared[[cohort]]
  all_tab <- table(z$meta$degradation_group)
  model_tab <- table(z$meta$degradation_group[z$complete])
  data.frame(
    cohort = cohort,
    n_primary = nrow(z$meta),
    median_cutoff = unique(z$meta$cutoff),
    low_n = unname(all_tab["Low"]),
    high_n = unname(all_tab["High"]),
    model_n = sum(z$complete),
    model_low_n = unname(model_tab["Low"]),
    model_high_n = unname(model_tab["High"]),
    subtypes = paste(levels(z$meta$subtype), collapse = ";"),
    score_reconstruction_max_abs_delta = z$max_delta,
    stringsAsFactors = FALSE
  )
}))
write.csv(cohort_summary, file.path(table_dir, "01_cohort_and_group_definition.csv"), row.names = FALSE)

fit_limma_ranks <- function(z, predictor = c("group", "continuous"), subset = NULL,
                            adjust_subtype = TRUE) {
  predictor <- match.arg(predictor)
  if (is.null(subset)) subset <- rep(TRUE, nrow(z$meta))
  keep <- z$complete & subset
  md <- droplevels(z$meta[keep, , drop = FALSE])
  if (nlevels(md$degradation_group) != 2L) stop("Both groups required.")
  if (adjust_subtype && nlevels(md$subtype) > 1L) {
    design <- if (predictor == "group") {
      model.matrix(~ degradation_group + subtype + age_z + sex, data = md)
    } else {
      model.matrix(~ degradation_z + subtype + age_z + sex, data = md)
    }
  } else {
    design <- if (predictor == "group") {
      model.matrix(~ degradation_group + age_z + sex, data = md)
    } else {
      model.matrix(~ degradation_z + age_z + sex, data = md)
    }
  }
  coef_name <- if (predictor == "group") "degradation_groupHigh" else "degradation_z"
  if (!coef_name %in% colnames(design)) stop("Predictor coefficient absent: ", coef_name)
  fit <- eBayes(lmFit(z$expr[, keep, drop = FALSE], design))
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  ranks <- fit$t[, coef_name]
  names(ranks) <- rownames(fit$t)
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  list(ranks = ranks, table = tt, n = nrow(md), meta = md)
}

run_full_gsea <- function(ranks) {
  serial_param <- BiocParallel::SerialParam()
  gh <- as.data.frame(fgseaMultilevel(
    pathways = hallmark_sets, stats = ranks,
    minSize = 10, maxSize = 500, eps = 0, nPermSimple = 10000,
    BPPARAM = serial_param
  ))
  gr <- as.data.frame(fgseaMultilevel(
    pathways = reactome_sets, stats = ranks,
    minSize = 10, maxSize = 500, eps = 0, nPermSimple = 10000,
    BPPARAM = serial_param
  ))
  gh$collection <- "Hallmark"
  gr$collection <- "Reactome"
  rbind(gh, gr)
}

run_targeted_gsea <- function(ranks) {
  serial_param <- BiocParallel::SerialParam()
  as.data.frame(fgseaMultilevel(
    pathways = exact_sets, stats = ranks,
    minSize = 10, maxSize = 500, eps = 0, nPermSimple = 10000,
    BPPARAM = serial_param
  ))
}

main_gsea <- list()
continuous_gsea <- list()
full_gsea_outputs <- list()
limma_outputs <- list()

for (cohort in names(prepared)) {
  cat("Fitting ", cohort, " High-versus-Low gene-level model...\n", sep = "")
  fg <- fit_limma_ranks(prepared[[cohort]], "group")
  limma_outputs[[paste0(cohort, "_group")]] <- fg
  full <- run_full_gsea(fg$ranks)
  full_gsea_outputs[[cohort]] <- full
  sel <- full[match(pathway_meta$pathway, full$pathway), , drop = FALSE]
  if (anyNA(sel$pathway)) stop(cohort, ": exact pathway missing from full GSEA.")
  sel$q_prespecified_4 <- p.adjust(sel$pval, method = "BH")
  sel$cohort <- cohort
  sel$model_n <- fg$n
  main_gsea[[cohort]] <- sel

  cat("Fitting ", cohort, " continuous degradation-score model...\n", sep = "")
  fc <- fit_limma_ranks(prepared[[cohort]], "continuous")
  limma_outputs[[paste0(cohort, "_continuous")]] <- fc
  cg <- run_targeted_gsea(fc$ranks)
  cg$q_prespecified_4 <- p.adjust(cg$pval, method = "BH")
  cg$cohort <- cohort
  cg$model_n <- fc$n
  continuous_gsea[[cohort]] <- cg

  full_write <- full
  full_write$leadingEdge <- vapply(full_write$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  write.csv(full_write,
            gzfile(file.path(table_dir, paste0("02_full_collection_GSEA_", cohort, ".csv.gz"))),
            row.names = FALSE)

  gene_group <- fg$table[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
  write.csv(gene_group,
            gzfile(file.path(table_dir, paste0("03_gene_level_limma_High_vs_Low_", cohort, ".csv.gz"))),
            row.names = FALSE)
}

main_gsea_df <- do.call(rbind, main_gsea)
main_gsea_df <- merge(pathway_meta, main_gsea_df, by = "pathway", all.x = TRUE, sort = FALSE)
main_gsea_df <- main_gsea_df[order(main_gsea_df$node_order, main_gsea_df$cohort), ]
main_gsea_df$direction <- ifelse(main_gsea_df$NES > 0, "High", "Low")
main_gsea_df$replication_pass <- main_gsea_df$NES > 0 & main_gsea_df$padj < 0.05
main_gsea_df$leading_edge_genes <- vapply(main_gsea_df$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
main_gsea_df$leadingEdge <- NULL
names(main_gsea_df)[names(main_gsea_df) == "padj"] <- "q_within_collection"
write.csv(main_gsea_df, file.path(table_dir, "04_prespecified_exact_pathway_GSEA_High_vs_Low.csv"), row.names = FALSE)

continuous_gsea_df <- do.call(rbind, continuous_gsea)
continuous_gsea_df <- merge(pathway_meta, continuous_gsea_df, by = "pathway", all.x = TRUE, sort = FALSE)
continuous_gsea_df <- continuous_gsea_df[order(continuous_gsea_df$node_order, continuous_gsea_df$cohort), ]
continuous_gsea_df$direction <- ifelse(continuous_gsea_df$NES > 0, "Higher degradation", "Lower degradation")
continuous_gsea_df$replication_pass <- continuous_gsea_df$NES > 0 & continuous_gsea_df$q_prespecified_4 < 0.05
continuous_gsea_df$leading_edge_genes <- vapply(continuous_gsea_df$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
continuous_gsea_df$leadingEdge <- NULL
names(continuous_gsea_df)[names(continuous_gsea_df) == "padj"] <- "q_fgsea_targeted"
write.csv(continuous_gsea_df, file.path(table_dir, "05_prespecified_exact_pathway_GSEA_continuous.csv"), row.names = FALSE)

## Strict, mutually gene-disjoint sample-level node scores.
score_long <- list()
score_effects <- list()
edge_effects <- list()
score_mapping <- list()
exact_score_long <- list()
exact_score_effects <- list()
exact_score_mapping <- list()

make_strict_scores <- function(z, cohort) {
  out <- data.frame(sample = z$meta$sample, stringsAsFactors = FALSE)
  for (nm in names(strict_sets)) {
    mapped <- intersect(strict_sets[[nm]], rownames(z$expr))
    sds <- apply(z$expr[mapped, , drop = FALSE], 1, sd, na.rm = TRUE)
    mapped <- mapped[is.finite(sds) & sds > 0]
    if (length(mapped) < 10L) stop(cohort, ": too few genes for strict score ", nm)
    zz <- t(scale(t(z$expr[mapped, , drop = FALSE])))
    raw_score <- colMeans(zz, na.rm = TRUE)
    out[[nm]] <- as.numeric(scale(raw_score))
    score_mapping[[length(score_mapping) + 1L]] <<- data.frame(
      cohort = cohort, score = nm,
      prespecified_n = length(strict_sets[[nm]]), mapped_nonconstant_n = length(mapped),
      mapped_genes = paste(mapped, collapse = ";"), stringsAsFactors = FALSE
    )
  }
  out
}

make_exact_scores <- function(z, cohort) {
  out <- data.frame(sample = z$meta$sample, stringsAsFactors = FALSE)
  for (nm in names(exact_score_sets)) {
    mapped <- intersect(exact_score_sets[[nm]], rownames(z$expr))
    sds <- apply(z$expr[mapped, , drop = FALSE], 1, sd, na.rm = TRUE)
    mapped <- mapped[is.finite(sds) & sds > 0]
    if (length(mapped) < 10L) stop(cohort, ": too few genes for exact score ", nm)
    zz <- t(scale(t(z$expr[mapped, , drop = FALSE])))
    raw_score <- colMeans(zz, na.rm = TRUE)
    out[[nm]] <- as.numeric(scale(raw_score))
    exact_score_mapping[[length(exact_score_mapping) + 1L]] <<- data.frame(
      cohort = cohort, score = nm,
      prespecified_n = length(exact_score_sets[[nm]]), mapped_nonconstant_n = length(mapped),
      mapped_genes = paste(mapped, collapse = ";"), stringsAsFactors = FALSE
    )
  }
  out
}

extract_lm <- function(fit, term) {
  co <- summary(fit)$coefficients
  if (!term %in% rownames(co)) stop("Term missing from lm: ", term)
  ci <- confint(fit, term, level = 0.95)
  c(beta = unname(co[term, "Estimate"]), se = unname(co[term, "Std. Error"]),
    p = unname(co[term, "Pr(>|t|)"]), lo = unname(ci[1]), hi = unname(ci[2]))
}

for (cohort in names(prepared)) {
  z <- prepared[[cohort]]
  ss <- make_strict_scores(z, cohort)
  es <- make_exact_scores(z, cohort)
  md <- cbind(z$meta, ss[match(z$meta$sample, ss$sample), setdiff(names(ss), "sample"), drop = FALSE])
  md_exact <- cbind(z$meta, es[match(z$meta$sample, es$sample), setdiff(names(es), "sample"), drop = FALSE])
  md$model_complete <- z$complete
  md_exact$model_complete <- z$complete
  long <- reshape(
    md[, c("sample", "cohort", "subtype", "degradation_group", names(strict_sets))],
    varying = names(strict_sets), v.names = "score_value", timevar = "node",
    times = names(strict_sets), direction = "long"
  )
  rownames(long) <- NULL
  score_long[[cohort]] <- long

  exact_long <- reshape(
    md_exact[, c("sample", "cohort", "subtype", "degradation_group", names(exact_score_sets))],
    varying = names(exact_score_sets), v.names = "score_value", timevar = "node",
    times = names(exact_score_sets), direction = "long"
  )
  rownames(exact_long) <- NULL
  exact_score_long[[cohort]] <- exact_long

  eff <- list()
  for (node in names(strict_sets)) {
    dat <- md[z$complete, , drop = FALSE]
    f_group <- lm(reformulate(c("degradation_group", "subtype", "age_z", "sex"), response = node), data = dat)
    xg <- extract_lm(f_group, "degradation_groupHigh")
    f_cont <- lm(reformulate(c("degradation_z", "subtype", "age_z", "sex"), response = node), data = dat)
    xc <- extract_lm(f_cont, "degradation_z")
    eff[[node]] <- rbind(
      data.frame(cohort = cohort, node = node, model = "High vs Low", t(xg), n = nobs(f_group)),
      data.frame(cohort = cohort, node = node, model = "Continuous", t(xc), n = nobs(f_cont))
    )
  }
  eff <- do.call(rbind, eff)
  eff$q <- ave(eff$p, eff$model, FUN = function(x) p.adjust(x, method = "BH"))
  score_effects[[cohort]] <- eff

  exact_eff <- list()
  for (node in names(exact_score_sets)) {
    dat <- md_exact[z$complete, , drop = FALSE]
    f_group <- lm(reformulate(c("degradation_group", "subtype", "age_z", "sex"), response = node), data = dat)
    xg <- extract_lm(f_group, "degradation_groupHigh")
    f_cont <- lm(reformulate(c("degradation_z", "subtype", "age_z", "sex"), response = node), data = dat)
    xc <- extract_lm(f_cont, "degradation_z")
    exact_eff[[node]] <- rbind(
      data.frame(cohort = cohort, node = node, model = "High vs Low", t(xg), n = nobs(f_group)),
      data.frame(cohort = cohort, node = node, model = "Continuous", t(xc), n = nobs(f_cont))
    )
  }
  exact_eff <- do.call(rbind, exact_eff)
  exact_eff$q <- ave(exact_eff$p, exact_eff$model, FUN = function(x) p.adjust(x, method = "BH"))
  exact_score_effects[[cohort]] <- exact_eff

  edge_defs <- list(
    "APC -> TCR" = c("APC", "TCR"),
    "TCR -> Noncanonical NF-kappaB" = c("TCR", "Noncanonical_NFkB"),
    "Noncanonical NF-kappaB -> IFNG" = c("Noncanonical_NFkB", "IFNG")
  )
  ee <- lapply(names(edge_defs), function(edge) {
    upstream <- edge_defs[[edge]][1]
    downstream <- edge_defs[[edge]][2]
    dat <- md[z$complete, , drop = FALSE]
    fit <- lm(reformulate(c(upstream, "degradation_group", "subtype", "age_z", "sex"),
                          response = downstream), data = dat)
    x <- extract_lm(fit, upstream)
    data.frame(cohort = cohort, edge = edge, upstream = upstream, downstream = downstream,
               t(x), n = nobs(fit), stringsAsFactors = FALSE)
  })
  ee <- do.call(rbind, ee)
  ee$q <- p.adjust(ee$p, method = "BH")
  edge_effects[[cohort]] <- ee
}

score_long_df <- do.call(rbind, score_long)
score_effects_df <- do.call(rbind, score_effects)
exact_score_long_df <- do.call(rbind, exact_score_long)
exact_score_effects_df <- do.call(rbind, exact_score_effects)
edge_effects_df <- do.call(rbind, edge_effects)
score_mapping_df <- do.call(rbind, score_mapping)
exact_score_mapping_df <- do.call(rbind, exact_score_mapping)
write.csv(score_long_df, file.path(table_dir, "06_gene_disjoint_sample_scores.csv"), row.names = FALSE)
write.csv(score_effects_df, file.path(table_dir, "07_adjusted_gene_disjoint_score_effects.csv"), row.names = FALSE)
write.csv(edge_effects_df, file.path(table_dir, "08_adjusted_adjacent_node_associations.csv"), row.names = FALSE)
write.csv(score_mapping_df, file.path(table_dir, "09_strict_score_mapping.csv"), row.names = FALSE)
write.csv(exact_score_long_df, file.path(table_dir, "09b_exact_pathway_sample_scores.csv"), row.names = FALSE)
write.csv(exact_score_effects_df, file.path(table_dir, "09c_adjusted_exact_pathway_score_effects.csv"), row.names = FALSE)
write.csv(exact_score_mapping_df, file.path(table_dir, "09d_exact_pathway_score_mapping.csv"), row.names = FALSE)

summary_lines <- c(
  "TCGA NSCLC sarcosine-degradation immune-axis validation",
  paste0("Analysis date: ", Sys.Date()),
  "Definition: degradation_z = mean(z(SARDH), z(PIPOX)); High if above the NSCLC median.",
  paste0("MSigDB version: ", msig_version),
  paste0("Primary tumors: ", cohort_summary$n_primary,
         "; complete adjusted model: ", cohort_summary$model_n),
  paste0("Exact GSEA pathways positive and collection-wide BH q<0.05: ",
         sum(main_gsea_df$replication_pass), "/4"),
  paste0("Exact pathway-score High-vs-Low effects positive and BH q<0.05: ",
         sum(exact_score_effects_df$model == "High vs Low" &
               exact_score_effects_df$beta > 0 & exact_score_effects_df$q < 0.05), "/4"),
  paste0("Gene-disjoint adjacent-node associations positive and BH q<0.05: ",
         sum(edge_effects_df$beta > 0 & edge_effects_df$q < 0.05), "/3"),
  "Interpretation: coordinated bulk-tumor program; temporal causality and same-cell signaling are not established."
)
writeLines(summary_lines, file.path(log_dir, "run_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
