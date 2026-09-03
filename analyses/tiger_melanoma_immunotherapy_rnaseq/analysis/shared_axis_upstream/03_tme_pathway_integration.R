#!/usr/bin/env Rscript

# TIGER PRJEB23709 PRE-only integrated pathway/TME analysis.
# Primary pathway comparison reuses the validated adjusted limma rank from analysis 104.
# Positive NES denotes enrichment in sarcosine-degradation Degradation-High.

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

required_packages <- c(
  "data.table", "limma", "fgsea", "msigdbr", "ggplot2", "BiocParallel", "statmod"
)
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
deconv_root <- file.path(
  analysis_root, "results", "TIGER_PRE73_Degradation_HighLow_Deconvolution_7methods_26.09.02"
)
output_root <- file.path(
  analysis_root, "results", "TIGER_PRE73_ThreeCollection_Immune_TME_Integrated_26.09.02"
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
primary_gene_table_path <- file.path(prior_root, "tables", "02_limma_adjusted_High_vs_Low_all_genes.csv")
continuous_gene_table_path <- file.path(prior_root, "tables", "04_limma_adjusted_continuous_score_all_genes.csv")
hallmark_table_path <- file.path(prior_root, "tables", "05_fgsea_Hallmark_adjusted_High_vs_Low.csv")
gobp_table_path <- file.path(prior_root, "tables", "06_fgsea_GO_BP_adjusted_High_vs_Low.csv")
estimate_path <- file.path(deconv_root, "raw_deconvolution", "deconvolution_ESTIMATE.csv")

required_inputs <- c(
  expression_path, metadata_path, primary_gene_table_path, continuous_gene_table_path,
  hallmark_table_path, gobp_table_path, estimate_path
)
for (path in required_inputs) if (!file.exists(path)) stop("Missing input: ", path)

immune_pattern <- paste(c(
  "IMMUN", "(^|_)T_CELL($|_)", "(^|_)B_CELL($|_)", "LEUKOCYTE", "LYMPHOCYTE",
  "INTERFERON", "CYTOKINE(S)?($|_|\\b)", "ANTIGEN", "NATURAL_KILLER", "NK_CELL",
  "MACROPHAGE", "MONOCYTE", "NEUTROPHIL", "DENDRITIC_CELL", "INFLAM",
  "(^|_)COMPLEMENT($|_)", "CHEMOKINE", "TOLL_LIKE", "TNF", "INTERLEUKIN",
  "MHC", "CYTOTOX", "PHAGOCYT", "GRANULOCYTE", "MYELOID", "ALLOGRAFT",
  "(^|_)IL[0-9]+", "JAK_STAT", "NF.?KB", "(^|_)TCR($|_)", "(^|_)BCR($|_)"
), collapse = "|")

add_immune_flag <- function(x) {
  x[, immune_related := grepl(immune_pattern, pathway, ignore.case = TRUE)]
  x[]
}

rank_from_table <- function(x) {
  stopifnot(all(c("gene_symbol", "moderated_t") %in% names(x)))
  stats <- x$moderated_t
  names(stats) <- x$gene_symbol
  stats <- sort(stats, decreasing = TRUE)
  stopifnot(!anyNA(stats), all(is.finite(stats)), !anyDuplicated(names(stats)))
  stats
}

primary_genes <- fread(primary_gene_table_path, check.names = FALSE)
continuous_genes <- fread(continuous_gene_table_path, check.names = FALSE)
rank_primary <- rank_from_table(primary_genes)
rank_continuous <- rank_from_table(continuous_genes)
hallmark_primary <- add_immune_flag(fread(hallmark_table_path, check.names = FALSE))
gobp_primary <- add_immune_flag(fread(gobp_table_path, check.names = FALSE))

reactome_membership <- as.data.table(msigdbr(
  species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"
))
hallmark_membership <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
gobp_membership <- as.data.table(msigdbr(
  species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"
))
stopifnot(
  uniqueN(reactome_membership$db_version) == 1L,
  identical(unique(reactome_membership$db_version), unique(hallmark_membership$db_version)),
  identical(unique(reactome_membership$db_version), unique(gobp_membership$db_version))
)
msigdb_version <- unique(reactome_membership$db_version)

reactome_pathways <- lapply(split(reactome_membership$gene_symbol, reactome_membership$gs_name), unique)
hallmark_pathways <- lapply(split(hallmark_membership$gene_symbol, hallmark_membership$gs_name), unique)
gobp_pathways <- lapply(split(gobp_membership$gene_symbol, gobp_membership$gs_name), unique)

metadata_all <- rbindlist(list(
  unique(hallmark_membership[, .(
    collection = "Hallmark", pathway = gs_name, description = gs_description,
    database_id = gs_exact_source, url = gs_url
  )]),
  unique(reactome_membership[, .(
    collection = "Reactome", pathway = gs_name, description = gs_description,
    database_id = gs_exact_source, url = gs_url
  )]),
  unique(gobp_membership[, .(
    collection = "GO:BP", pathway = gs_name, description = gs_description,
    database_id = gs_exact_source, url = gs_url
  )])
), use.names = TRUE)
metadata_all[, MSigDB_version := msigdb_version]

run_fgsea <- function(pathways, stats, collection_label, min_size, max_size, model_label) {
  set.seed(42)
  result <- as.data.table(fgseaMultilevel(
    pathways = pathways, stats = stats, minSize = min_size, maxSize = max_size,
    eps = 0, scoreType = "std", nproc = 1,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  ))
  result[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
  result[, `:=`(
    collection = collection_label,
    model = model_label,
    direction = fifelse(NES > 0, "Higher in Degradation-High", "Higher in Degradation-Low")
  )]
  result <- merge(result, metadata_all[collection == collection_label],
                  by = c("collection", "pathway"), all.x = TRUE)
  result <- add_immune_flag(result)
  result[, abs_NES_for_sort := abs(NES)]
  setorder(result, padj, -abs_NES_for_sort, pathway)
  result[, abs_NES_for_sort := NULL]
  result[]
}

run_camera <- function(pathways, stats, collection_label, min_size, max_size, model_label) {
  index <- ids2indices(pathways, names(stats), remove.empty = TRUE)
  index <- index[lengths(index) >= min_size & lengths(index) <= max_size]
  result <- as.data.table(cameraPR(
    statistic = stats, index = index, use.ranks = FALSE,
    inter.gene.cor = 0.01, sort = FALSE, directional = TRUE
  ), keep.rownames = "pathway")
  result[, `:=`(collection = collection_label, model = model_label)]
  result <- merge(result, metadata_all[collection == collection_label],
                  by = c("collection", "pathway"), all.x = TRUE)
  result <- add_immune_flag(result)
  setorder(result, FDR, PValue, pathway)
  result[]
}

reactome_primary <- run_fgsea(
  reactome_pathways, rank_primary, "Reactome", 15L, 500L,
  "Adjusted Degradation High vs Low"
)
reactome_no_score <- run_fgsea(
  reactome_pathways, rank_primary[!names(rank_primary) %in% c("SARDH", "PIPOX")],
  "Reactome", 15L, 500L, "Adjusted High vs Low; SARDH/PIPOX excluded"
)
reactome_continuous_camera <- run_camera(
  reactome_pathways, rank_continuous, "Reactome", 15L, 500L,
  "Adjusted continuous Degradation score"
)
reactome_primary_camera <- run_camera(
  reactome_pathways, rank_primary, "Reactome", 15L, 500L,
  "Adjusted Degradation High vs Low"
)

fwrite(reactome_primary, file.path(table_dir, "01_fgsea_Reactome_adjusted_High_vs_Low.csv"))
fwrite(reactome_primary[immune_related == TRUE], file.path(table_dir, "02_fgsea_Reactome_immune_subset.csv"))
fwrite(reactome_no_score, file.path(table_dir, "03_fgsea_Reactome_score_genes_excluded.csv"))
fwrite(reactome_primary_camera, file.path(table_dir, "04_cameraPR_Reactome_adjusted.csv"))
fwrite(reactome_continuous_camera, file.path(table_dir, "05_cameraPR_Reactome_continuous.csv"))

# Cross-collection correspondence is descriptive. It does not treat ontologies as independent studies.
parse_leading_edge <- function(x) {
  if (is.na(x) || !nzchar(x)) character() else unique(strsplit(x, ";", fixed = TRUE)[[1]])
}

pair_overlap <- function(source_genes, target_row) {
  target_genes <- parse_leading_edge(target_row$leadingEdge)
  shared <- intersect(source_genes, target_genes)
  union_genes <- union(source_genes, target_genes)
  data.table(
    target_pathway = target_row$pathway,
    target_NES = target_row$NES,
    target_padj = target_row$padj,
    target_size = length(target_genes),
    shared_genes_n = length(shared),
    overlap_coefficient = if (min(length(source_genes), length(target_genes)) > 0) {
      length(shared) / min(length(source_genes), length(target_genes))
    } else 0,
    jaccard = if (length(union_genes) > 0) length(shared) / length(union_genes) else 0,
    shared_genes = paste(sort(shared), collapse = ";")
  )
}

positive_hallmark <- hallmark_primary[
  immune_related == TRUE & NES > 0 & !is.na(padj) & padj < 0.05
]
positive_reactome <- reactome_primary[
  immune_related == TRUE & NES > 0 & !is.na(padj) & padj < 0.05
]
positive_gobp <- gobp_primary[
  immune_related == TRUE & NES > 0 & !is.na(padj) & padj < 0.05
]

best_match <- function(source_row, target_dt, target_collection) {
  source_genes <- parse_leading_edge(source_row$leadingEdge)
  candidates <- rbindlist(lapply(seq_len(nrow(target_dt)), function(i) {
    pair_overlap(source_genes, target_dt[i])
  }))
  candidates[, target_collection := target_collection]
  setorder(candidates, -shared_genes_n, -overlap_coefficient, target_padj, -target_NES, target_pathway)
  candidates[1]
}

cross_collection <- rbindlist(lapply(seq_len(nrow(positive_hallmark)), function(i) {
  h <- positive_hallmark[i]
  h_genes <- parse_leading_edge(h$leadingEdge)
  r <- best_match(h, positive_reactome, "Reactome")
  g <- best_match(h, positive_gobp, "GO:BP")
  r_genes <- parse_leading_edge(positive_reactome[pathway == r$target_pathway]$leadingEdge[[1]])
  g_genes <- parse_leading_edge(positive_gobp[pathway == g$target_pathway]$leadingEdge[[1]])
  triple <- Reduce(intersect, list(h_genes, r_genes, g_genes))
  rg_shared <- intersect(r_genes, g_genes)
  data.table(
    program = sub("^HALLMARK_", "", h$pathway),
    Hallmark_pathway = h$pathway, Hallmark_NES = h$NES, Hallmark_q = h$padj,
    Reactome_pathway = r$target_pathway, Reactome_NES = r$target_NES,
    Reactome_q = r$target_padj, Hallmark_Reactome_shared_n = r$shared_genes_n,
    Hallmark_Reactome_overlap = r$overlap_coefficient,
    GO_BP_pathway = g$target_pathway, GO_BP_NES = g$target_NES,
    GO_BP_q = g$target_padj, Hallmark_GO_BP_shared_n = g$shared_genes_n,
    Hallmark_GO_BP_overlap = g$overlap_coefficient,
    Reactome_GO_BP_shared_n = length(rg_shared),
    triple_shared_genes_n = length(triple),
    triple_shared_genes = paste(sort(triple), collapse = ";")
  )
}), fill = TRUE)

cross_collection[, robust_correspondence :=
  Hallmark_Reactome_shared_n >= 5L & Hallmark_GO_BP_shared_n >= 5L &
  Hallmark_Reactome_overlap >= 0.20 & Hallmark_GO_BP_overlap >= 0.20 &
  Reactome_GO_BP_shared_n >= 3L & triple_shared_genes_n >= 2L
]
setorder(cross_collection, -robust_correspondence, -triple_shared_genes_n,
         -Hallmark_Reactome_shared_n, -Hallmark_GO_BP_shared_n)
fwrite(cross_collection, file.path(table_dir, "06_cross_collection_immune_program_correspondence.csv"))

plot_programs <- cross_collection[robust_correspondence == TRUE]
if (nrow(plot_programs) > 8L) plot_programs <- plot_programs[1:8]
if (!nrow(plot_programs)) plot_programs <- cross_collection[1:min(.N, 8L)]
program_long <- rbindlist(list(
  plot_programs[, .(program, collection = "Hallmark", NES = Hallmark_NES, q = Hallmark_q)],
  plot_programs[, .(program, collection = "Reactome", NES = Reactome_NES, q = Reactome_q)],
  plot_programs[, .(program, collection = "GO:BP", NES = GO_BP_NES, q = GO_BP_q)]
))
program_long[, program_label := gsub("_", " ", program)]
program_order <- rev(gsub("_", " ", plot_programs$program))
program_long[, program_label := factor(program_label, levels = program_order)]
program_long[, collection := factor(collection, levels = c("Hallmark", "Reactome", "GO:BP"))]
program_long[, star := fifelse(q < 0.001, "***", fifelse(q < 0.01, "**", fifelse(q < 0.05, "*", "")))]
fwrite(program_long, file.path(table_dir, "07_cross_collection_program_NES_long.csv"))

p_program <- ggplot(program_long, aes(x = collection, y = program_label, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = star), fontface = "bold", size = 4.2) +
  scale_fill_gradient2(low = "#1B9E8F", mid = "white", high = "#C43C3C", midpoint = 0,
                       name = "NES") +
  labs(
    title = "Immune programs concordantly enriched across three gene-set collections",
    subtitle = "Positive NES = Degradation-High; terms matched by leading-edge overlap; * BH q<0.05",
    x = NULL, y = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(face = "bold"), panel.grid = element_blank(), legend.position = "right")
ggsave(file.path(figure_dir, "Fig_CrossCollection_ImmuneProgram_NES_heatmap.png"),
       p_program, width = 9.4, height = max(5.5, 2.8 + 0.48 * nrow(plot_programs)), dpi = 600, bg = "white")
ggsave(file.path(figure_dir, "Fig_CrossCollection_ImmuneProgram_NES_heatmap.pdf"),
       p_program, width = 9.4, height = max(5.5, 2.8 + 0.48 * nrow(plot_programs)), device = cairo_pdf)

# TME-adjusted bulk-tumour sensitivity model. This is not a purified cancer-cell analysis.
expr_dt <- fread(expression_path, check.names = FALSE)
meta <- fread(metadata_path, check.names = FALSE)
estimate <- fread(estimate_path, check.names = FALSE)
stopifnot(nrow(meta) == 73L, uniqueN(meta$sample_id) == 73L, all(meta$timepoint == "PRE"))
stopifnot(setequal(estimate$ID, meta$sample_id))
estimate <- estimate[match(meta$sample_id, ID)]
stopifnot(identical(estimate$ID, meta$sample_id))

sample_ids <- meta$sample_id
gene_symbols <- trimws(expr_dt[[1]])
fpkm <- as.matrix(expr_dt[, ..sample_ids])
storage.mode(fpkm) <- "double"
rownames(fpkm) <- gene_symbols
colnames(fpkm) <- sample_ids
minimum_samples <- ceiling(0.10 * ncol(fpkm))
keep <- rowSums(fpkm >= 1) >= minimum_samples
log_expression <- log2(fpkm[keep, , drop = FALSE] + 1)
log_expression <- log_expression[apply(log_expression, 1, var) > 0, , drop = FALSE]
stopifnot(identical(rownames(log_expression), primary_genes$gene_symbol))

degradation_median <- median(meta$Degradation_score)
meta[, degradation_group := factor(
  ifelse(Degradation_score > degradation_median, "High", "Low"), levels = c("Low", "High")
)]
meta[, therapy_short := relevel(factor(therapy_short), ref = "antiPD1")]
meta[, response_group := relevel(factor(response_group), ref = "NR")]
meta[, gender := relevel(factor(gender), ref = "Male")]
meta[, age_z := as.numeric(scale(age))]
meta[, immune_score_z := as.numeric(scale(estimate$ImmuneScore_estimate))]
meta[, stromal_score_z := as.numeric(scale(estimate$StromalScore_estimate))]

design_tme_adjusted <- model.matrix(
  ~ therapy_short + response_group + age_z + gender + immune_score_z + stromal_score_z + degradation_group,
  data = meta
)
stopifnot(qr(design_tme_adjusted)$rank == ncol(design_tme_adjusted))
tme_coef <- match("degradation_groupHigh", colnames(design_tme_adjusted))
fit_tme <- eBayes(lmFit(log_expression, design_tme_adjusted), trend = TRUE, robust = TRUE)
tme_gene_table <- as.data.table(topTable(fit_tme, coef = tme_coef, number = Inf, sort.by = "none"), keep.rownames = "gene_symbol")
setnames(tme_gene_table, c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"),
         c("log2FC", "average_log2_FPKM_plus1", "moderated_t", "p_value", "BH_q", "B_statistic"))
tme_gene_table[, estimand := "Degradation High minus Low; adjusted for clinical covariates plus ESTIMATE immune/stromal scores"]
rank_tme <- rank_from_table(tme_gene_table)
fwrite(tme_gene_table, file.path(table_dir, "08_limma_TME_adjusted_High_vs_Low_all_genes.csv"))

hallmark_tme <- run_fgsea(hallmark_pathways, rank_tme, "Hallmark", 10L, 500L,
                          "Clinical plus ESTIMATE immune/stromal adjusted")
reactome_tme <- run_fgsea(reactome_pathways, rank_tme, "Reactome", 15L, 500L,
                          "Clinical plus ESTIMATE immune/stromal adjusted")
gobp_tme <- run_fgsea(gobp_pathways, rank_tme, "GO:BP", 15L, 500L,
                      "Clinical plus ESTIMATE immune/stromal adjusted")
fwrite(hallmark_tme, file.path(table_dir, "09_fgsea_Hallmark_TME_adjusted.csv"))
fwrite(reactome_tme, file.path(table_dir, "10_fgsea_Reactome_TME_adjusted.csv"))
fwrite(gobp_tme, file.path(table_dir, "11_fgsea_GO_BP_TME_adjusted.csv"))

select_tme_terms <- function(x, n_each = 5L) {
  sig <- x[!is.na(padj) & padj < 0.05]
  pos <- sig[NES > 0][order(padj, -NES)][1:min(.N, n_each)]
  neg <- sig[NES < 0][order(padj, NES)][1:min(.N, n_each)]
  rbind(pos, neg)
}
tme_plot <- rbindlist(list(
  select_tme_terms(hallmark_tme)[, .(collection, pathway, NES, padj, immune_related)],
  select_tme_terms(reactome_tme)[, .(collection, pathway, NES, padj, immune_related)],
  select_tme_terms(gobp_tme)[, .(collection, pathway, NES, padj, immune_related)]
), use.names = TRUE)
tme_plot[, pathway_label := sub("^(HALLMARK_|REACTOME_|GOBP_)", "", pathway)]
tme_plot[, pathway_label := gsub("_", " ", pathway_label)]
tme_plot[, pathway_label := vapply(pathway_label, function(z) paste(strwrap(z, width = 42), collapse = "\n"), character(1))]
tme_plot[, display_label := paste(pathway_label, collection, sep = " — ")]
tme_plot <- tme_plot[order(NES)]
tme_plot[, display_label := factor(display_label, levels = display_label)]
fwrite(tme_plot, file.path(table_dir, "12_top_TME_adjusted_pathways_for_plot.csv"))

p_tme_pathway <- ggplot(tme_plot, aes(x = NES, y = display_label, color = NES > 0, size = -log10(pmax(padj, 1e-300)))) +
  geom_vline(xintercept = 0, color = "grey60", linewidth = 0.45) +
  geom_point(alpha = 0.9) +
  scale_color_manual(values = c(`TRUE` = "#C43C3C", `FALSE` = "#1B9E8F"),
                     labels = c(`TRUE` = "Degradation-High", `FALSE` = "Degradation-Low"), name = NULL) +
  scale_size_continuous(range = c(2.8, 6.5), name = expression(-log[10](BH~q))) +
  labs(
    title = "Pathway associations after adjusting for broad immune and stromal scores",
    subtitle = "Bulk-tumour sensitivity analysis; persistence does not establish cancer-cell specificity",
    x = "Normalized enrichment score (High vs Low)", y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
ggsave(file.path(figure_dir, "Fig_TME_Adjusted_ThreeCollection_Pathways.png"),
       p_tme_pathway, width = 11.5, height = max(8, 3.5 + 0.31 * nrow(tme_plot)), dpi = 600, bg = "white")
ggsave(file.path(figure_dir, "Fig_TME_Adjusted_ThreeCollection_Pathways.pdf"),
       p_tme_pathway, width = 11.5, height = max(8, 3.5 + 0.31 * nrow(tme_plot)), device = cairo_pdf)

count_summary <- data.table(
  analysis = c(
    "Hallmark primary", "Reactome primary", "GO:BP primary",
    "Hallmark TME-adjusted", "Reactome TME-adjusted", "GO:BP TME-adjusted"
  ),
  tested = c(nrow(hallmark_primary), nrow(reactome_primary), nrow(gobp_primary),
             nrow(hallmark_tme), nrow(reactome_tme), nrow(gobp_tme)),
  BH_q_lt_0_05 = c(
    sum(hallmark_primary$padj < 0.05, na.rm = TRUE),
    sum(reactome_primary$padj < 0.05, na.rm = TRUE),
    sum(gobp_primary$padj < 0.05, na.rm = TRUE),
    sum(hallmark_tme$padj < 0.05, na.rm = TRUE),
    sum(reactome_tme$padj < 0.05, na.rm = TRUE),
    sum(gobp_tme$padj < 0.05, na.rm = TRUE)
  ),
  immune_positive_BH_q_lt_0_05 = c(
    sum(hallmark_primary$padj < 0.05 & hallmark_primary$NES > 0 & hallmark_primary$immune_related, na.rm = TRUE),
    sum(reactome_primary$padj < 0.05 & reactome_primary$NES > 0 & reactome_primary$immune_related, na.rm = TRUE),
    sum(gobp_primary$padj < 0.05 & gobp_primary$NES > 0 & gobp_primary$immune_related, na.rm = TRUE),
    sum(hallmark_tme$padj < 0.05 & hallmark_tme$NES > 0 & hallmark_tme$immune_related, na.rm = TRUE),
    sum(reactome_tme$padj < 0.05 & reactome_tme$NES > 0 & reactome_tme$immune_related, na.rm = TRUE),
    sum(gobp_tme$padj < 0.05 & gobp_tme$NES > 0 & gobp_tme$immune_related, na.rm = TRUE)
  )
)
fwrite(count_summary, file.path(table_dir, "13_enrichment_count_summary.csv"))

model_diagnostics <- data.table(
  metric = c(
    "PRE samples", "Degradation High", "Degradation Low",
    "ESTIMATE immune-stromal correlation", "TME-adjusted design rank",
    "TME-adjusted design columns", "TME-adjusted design kappa"
  ),
  value = c(
    nrow(meta), sum(meta$degradation_group == "High"), sum(meta$degradation_group == "Low"),
    cor(meta$immune_score_z, meta$stromal_score_z), qr(design_tme_adjusted)$rank,
    ncol(design_tme_adjusted), kappa(design_tme_adjusted)
  )
)
fwrite(model_diagnostics, file.path(table_dir, "00_model_diagnostics.csv"))

method_contract <- data.table(
  field = c(
    "primary_rank", "Reactome_database", "primary_enrichment", "multiple_testing",
    "cross_collection_mapping", "robust_correspondence_rule", "TME_sensitivity_model",
    "direction", "interpretive_boundary"
  ),
  value = c(
    "Validated analysis-104 clinical-covariate-adjusted limma moderated-t rank",
    paste0("MSigDB ", msigdb_version, " C2:CP:REACTOME via msigdbr ", packageVersion("msigdbr")),
    "fgseaMultilevel preranked GSEA; Reactome cameraPR and score-gene-removal sensitivities",
    "BH FDR within each collection/model",
    "Best positive significant immune Reactome and GO:BP term for each positive significant immune Hallmark, ranked by leading-edge overlap",
    "Each Hallmark-target overlap >=0.20 and >=5 genes; Reactome-GO shared >=3; triple shared >=2",
    "Clinical covariates plus standardized ESTIMATE immune and stromal scores; sensitivity only",
    "Positive NES = higher in Degradation-High",
    "Bulk tumour observational associations; neither deconvolution nor TME adjustment proves cell counts, sarcosine flux, causality, or cancer-cell-specific activity"
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
    "Primary rank finite", "Reactome sets tested", "Reactome q in range",
    "Reactome score genes removed", "cross-collection rows equal positive immune Hallmarks",
    "TME design full rank", "TME rank finite", "TME q in range", "figures created",
    "same retained genes as analysis 104"
  ),
  pass = c(
    all(is.finite(rank_primary)), nrow(reactome_primary) > 1000L,
    all(is.na(reactome_primary$padj) | reactome_primary$padj >= 0 & reactome_primary$padj <= 1),
    !any(c("SARDH", "PIPOX") %in% names(rank_primary[!names(rank_primary) %in% c("SARDH", "PIPOX")])),
    nrow(cross_collection) == nrow(positive_hallmark),
    qr(design_tme_adjusted)$rank == ncol(design_tme_adjusted), all(is.finite(rank_tme)),
    all(is.na(hallmark_tme$padj) | hallmark_tme$padj >= 0 & hallmark_tme$padj <= 1) &&
      all(is.na(reactome_tme$padj) | reactome_tme$padj >= 0 & reactome_tme$padj <= 1) &&
      all(is.na(gobp_tme$padj) | gobp_tme$padj >= 0 & gobp_tme$padj <= 1),
    all(file.exists(file.path(figure_dir, c(
      "Fig_CrossCollection_ImmuneProgram_NES_heatmap.png",
      "Fig_TME_Adjusted_ThreeCollection_Pathways.png"
    )))),
    identical(rownames(log_expression), primary_genes$gene_symbol)
  )
)
fwrite(validation, file.path(table_dir, "14_validation_checks.csv"))
if (!all(validation$pass)) stop("Validation failure; inspect 14_validation_checks.csv")

report <- c(
  "# TIGER PRE73 three-collection pathway and TME-integrated analysis",
  "",
  "## Scope",
  "- PRE-only n=73; Degradation-High n=36 and Low n=37.",
  "- Primary pathway rank is the validated clinical-covariate-adjusted limma rank from analysis 104.",
  paste0("- Gene sets: MSigDB ", msigdb_version, " Hallmark, Reactome and GO:BP; fgseaMultilevel; BH within each collection."),
  "- Positive NES denotes Degradation-High.",
  "",
  "## Result counts",
  paste(capture.output(print(count_summary)), collapse = "\n"),
  "",
  "## Cross-collection immune program correspondence",
  paste(capture.output(print(cross_collection[, .(
    program, Hallmark_NES, Reactome_pathway, Reactome_NES, GO_BP_pathway, GO_BP_NES,
    triple_shared_genes_n, robust_correspondence
  )])), collapse = "\n"),
  "",
  "## Interpretation boundary",
  "Cross-collection correspondence is a descriptive leading-edge overlap analysis; ontologies are not independent studies. TME-adjusted bulk-RNA associations are tumour-intrinsic candidates only, not purified cancer-cell evidence. Deconvolution and ESTIMATE scores can reflect expression state, tumour purity and reference mismatch as well as cell abundance."
)
writeLines(report, file.path(output_root, "FINAL_INTEGRATED_REPORT.md"), useBytes = TRUE)
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))

cat("Integrated three-collection/TME analysis complete\n")
cat("Output:", output_root, "\n")
print(count_summary)
print(cross_collection[, .(program, triple_shared_genes_n, robust_correspondence)])
