#!/usr/bin/env Rscript

## Independent checks for the NSCLC-only public validation output.

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
project_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
table_dir <- file.path(project_dir, "results", "tables")
log_dir <- file.path(project_dir, "results", "logs")

raw <- readRDS(file.path(project_dir, "data", "01_xena_raw.rds"))
scores <- readRDS(file.path(project_dir, "data", "02_scores.rds"))
groups <- read.csv(file.path(table_dir, "01_cohort_and_group_definition.csv"))
gsea <- read.csv(file.path(table_dir, "04_prespecified_exact_pathway_GSEA_High_vs_Low.csv"), check.names = FALSE)
exact_long <- read.csv(file.path(table_dir, "09b_exact_pathway_sample_scores.csv"), check.names = FALSE)
exact_eff <- read.csv(file.path(table_dir, "09c_adjusted_exact_pathway_score_effects.csv"), check.names = FALSE)
strict_map <- read.csv(file.path(table_dir, "09_strict_score_mapping.csv"), check.names = FALSE)
edges <- read.csv(file.path(table_dir, "08_adjusted_adjacent_node_associations.csv"), check.names = FALSE)

checks <- character()
check <- function(condition, label) {
  if (!isTRUE(condition)) stop("QC FAIL: ", label, call. = FALSE)
  checks <<- c(checks, paste0("PASS | ", label))
}
near <- function(x, y, tol = 1e-8) {
  isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol))
}

check(identical(names(raw), "NSCLC") && identical(names(scores), "NSCLC"),
      "public input contains NSCLC only")
s <- scores$NSCLC
check(nrow(s) == 1017L, "expected 1,017 primary tumors before covariate completeness")
check(!anyDuplicated(s$sample), "unique sample barcodes")
check(!anyDuplicated(substr(s$sample, 1, 12)), "one primary tumor per patient")
check(all(substr(s$sample, 14, 15) == "01"), "TCGA primary-tumor sample type")

reconstructed <- rowMeans(cbind(as.numeric(scale(s$SARDH)), as.numeric(scale(s$PIPOX))))
check(max(abs(reconstructed - s$degradation_z)) < 1e-10,
      "degradation score reconstructs from SARDH and PIPOX")
check(nrow(groups) == 1L && groups$cohort == "NSCLC" && groups$model_n == 988L,
      "NSCLC-only cohort summary and expected adjusted-model size")
check(nrow(gsea) == 4L && all(gsea$cohort == "NSCLC"), "four prespecified exact pathways")
check(all(gsea$NES > 0 & gsea$q_within_collection < 0.05),
      "four exact pathways favor Degradation-High at collection-wide BH q<0.05")
check(near(gsea$q_prespecified_4, p.adjust(gsea$pval, method = "BH")),
      "prespecified four-set BH values reproduce")

check(nrow(exact_eff) == 8L && all(exact_eff$cohort == "NSCLC"),
      "four exact scores under group and continuous models")
check(all(exact_eff$beta > 0 & exact_eff$q < 0.05),
      "all exact-score effects are positive and BH-significant")
check(nrow(edges) == 3L && all(edges$cohort == "NSCLC"),
      "three adjacent-node associations")
check(all(edges$beta > 0 & edges$q < 0.05),
      "all adjacent-node associations are positive and BH-significant")

gene_lists <- lapply(strict_map$mapped_genes, function(x) strsplit(x, ";", fixed = TRUE)[[1]])
check(all(strict_map$cohort == "NSCLC"), "gene-disjoint manifest contains NSCLC only")
check(all(vapply(gene_lists, function(x) !any(c("SARDH", "PIPOX") %in% x), logical(1))),
      "immune-node scores exclude the exposure genes")
for (pair in combn(seq_along(gene_lists), 2, simplify = FALSE)) {
  check(length(intersect(gene_lists[[pair[1]]], gene_lists[[pair[2]]])) == 0L,
        paste("gene-disjoint nodes", strict_map$score[pair[1]], "and", strict_map$score[pair[2]]))
}

cl <- raw$NSCLC$clin
ci <- match(s$sample, cl$sampleID)
check(!anyNA(ci), "clinical rows align to score samples")
cl <- cl[ci, , drop = FALSE]
disease <- as.character(cl$`_primary_disease`)
subtype <- factor(ifelse(grepl("adenocarcinoma", disease, ignore.case = TRUE), "LUAD", "LUSC"))
age <- suppressWarnings(as.numeric(as.character(cl$age_at_initial_pathologic_diagnosis)))
sex <- toupper(trimws(as.character(cl$gender)))
sex[!sex %in% c("FEMALE", "MALE")] <- NA_character_
base <- data.frame(
  sample = s$sample,
  degradation_group = factor(ifelse(s$degradation_z > median(s$degradation_z), "High", "Low"),
                             levels = c("Low", "High")),
  subtype = subtype,
  age_z = as.numeric(scale(age)),
  sex = factor(sex)
)
for (node in unique(exact_long$node)) {
  nd <- exact_long[exact_long$node == node, c("sample", "score_value")]
  dat <- merge(base, nd, by = "sample", sort = FALSE)
  fit <- lm(score_value ~ degradation_group + subtype + age_z + sex, data = dat)
  beta <- coef(fit)["degradation_groupHigh"]
  recorded <- exact_eff[exact_eff$node == node & exact_eff$model == "High vs Low", ]
  check(nrow(recorded) == 1L && near(beta, recorded$beta),
        paste(node, "adjusted score effect refits exactly"))
}

dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(c(
  "TCGA NSCLC immune-axis verification: PASS",
  paste0("Checks passed: ", length(checks)),
  checks
), file.path(log_dir, "verification_PASS.txt"))
cat(paste(c("VERIFICATION PASS", checks), collapse = "\n"), "\n")
