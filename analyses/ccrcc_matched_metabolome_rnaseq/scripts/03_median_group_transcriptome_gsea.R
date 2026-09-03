#!/usr/bin/env Rscript

# Matched TJ-RCC tumour tissue Sarcosine High/Low transcriptome and GSEA analysis.
# Biological comparison shown to readers: normalized GC-MS Sarcosine High versus Low.
# Group definition: tumour median 134.4850986; High n=50, Low n=50.
# Positive effects/NES mean higher expression or enrichment in Sarcosine-High.

options(stringsAsFactors = FALSE, width = 180)
set.seed(260902)

analysis_root <- normalizePath(getwd())
if (!all(file.exists(file.path(analysis_root, c("RNAseq_Data", "Metabolomics_Data", "analysis"))))) {
  stop("Run from the ccRCC matched multi-omics analysis root: ", analysis_root)
}

project_root <- normalizePath(file.path(analysis_root, ".."))
gsea_lib <- file.path(
  project_root,
  "2024_Drug_Res_Updates_NSCLC", "RNA-seq공공데이터_GSE207422",
  "analysis_sarcosine_FINAL_26.08.25", "03_lee_fig3_style_FINAL",
  "hallmark_GSEA_Q4_vs_Q1_26.08.26", "R_libs"
)
if (!dir.exists(gsea_lib)) stop("Missing frozen local GSEA R library: ", gsea_lib)
.libPaths(unique(c(normalizePath(gsea_lib), .libPaths())))

required_packages <- c(
  "data.table", "limma", "fgsea", "msigdbr", "ggplot2", "pheatmap",
  "BiocParallel", "grid"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(pheatmap)
})

expression_path <- file.path(analysis_root, "RNAseq_Data", "bulkRNA_matrix_TPM.csv")
group_path <- file.path(
  analysis_root, "results", "00_input_audit_26.09.02", "tumor_sarcosine_group_map.csv"
)
for (path in c(expression_path, group_path)) if (!file.exists(path)) stop("Missing input: ", path)

output_root <- file.path(
  analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_HighLow_Transcriptome_GSEA_26.09.02"
)
table_dir <- file.path(output_root, "tables")
figure_dir <- file.path(output_root, "figures")
log_dir <- file.path(output_root, "logs")
for (d in c(table_dir, figure_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

expr_dt <- fread(expression_path, check.names = FALSE)
meta <- fread(group_path, check.names = FALSE)
stopifnot(nrow(meta) == 100L, uniqueN(meta$sample_id) == 100L)
stopifnot(all(c(
  "sample_id", "sarcosine_normalized_intensity", "log2_sarcosine_intensity",
  "tumor_median_cutpoint", "sarcosine_group", "batch", "sex", "age", "stage", "grade"
) %in% names(meta)))
stopifnot(!anyNA(meta[, .(
  sample_id, sarcosine_normalized_intensity, log2_sarcosine_intensity,
  sarcosine_group, batch, sex, age, stage, grade
)]))

cutpoint <- unique(meta$tumor_median_cutpoint)
stopifnot(length(cutpoint) == 1L, abs(cutpoint - 134.4850986) < 1e-8)
stopifnot(all(meta$sarcosine_group == ifelse(meta$sarcosine_normalized_intensity > cutpoint, "High", "Low")))
meta[, sarcosine_group := factor(sarcosine_group, levels = c("Low", "High"))]
stopifnot(sum(meta$sarcosine_group == "Low") == 50L, sum(meta$sarcosine_group == "High") == 50L)
stopifnot(sum(meta$sarcosine_normalized_intensity == cutpoint) == 0L)

gene_symbols <- trimws(expr_dt[[1]])
stopifnot(!anyNA(gene_symbols), all(nzchar(gene_symbols)), !anyDuplicated(gene_symbols))
sample_columns <- names(expr_dt)[-1]
stopifnot(all(meta$sample_id %in% sample_columns))

sample_ids <- meta$sample_id
log_expression_all <- as.matrix(expr_dt[, ..sample_ids])
storage.mode(log_expression_all) <- "double"
rownames(log_expression_all) <- gene_symbols
colnames(log_expression_all) <- sample_ids
stopifnot(identical(colnames(log_expression_all), meta$sample_id))
stopifnot(!anyNA(log_expression_all), all(is.finite(log_expression_all)), all(log_expression_all >= 0))

# The official file is named TPM and contains non-negative, already log-scale-like values.
# A value >=1 corresponds to an inferred linear value >=1 under the observed log2(x+1) pattern.
minimum_samples <- ceiling(0.10 * ncol(log_expression_all))
keep_expression <- rowSums(log_expression_all >= 1) >= minimum_samples
log_expression <- log_expression_all[keep_expression, , drop = FALSE]
keep_variance <- apply(log_expression, 1L, var) > 0
log_expression <- log_expression[keep_variance, , drop = FALSE]
stopifnot(nrow(log_expression) > 5000L)

meta[, batch := factor(batch)]
meta[, sex := relevel(factor(sex), ref = "female")]
meta[, age := factor(age, levels = c("40-60", "<40", ">60"))]
meta[, grade := factor(grade, levels = c("1", "2", "3", "4"))]
meta[, stage := factor(stage, levels = c("stage I", "stage II", "stage III", "stage IV"))]
meta[, sarcosine_z := as.numeric(scale(log2_sarcosine_intensity))]

design_primary <- model.matrix(~ batch + sex + age + grade + sarcosine_group, data = meta)
design_clinical_no_batch <- model.matrix(~ sex + age + grade + sarcosine_group, data = meta)
design_unadjusted <- model.matrix(~ sarcosine_group, data = meta)
design_continuous <- model.matrix(~ batch + sex + age + grade + sarcosine_z, data = meta)
design_stage <- model.matrix(~ batch + sex + age + stage + sarcosine_group, data = meta)

designs <- list(
  primary_adjusted_High_vs_Low = design_primary,
  clinical_no_batch_High_vs_Low = design_clinical_no_batch,
  unadjusted_High_vs_Low = design_unadjusted,
  primary_continuous_log2_Sarcosine = design_continuous,
  stage_instead_of_grade_High_vs_Low = design_stage
)
for (nm in names(designs)) {
  if (qr(designs[[nm]])$rank != ncol(designs[[nm]])) stop("Rank-deficient design: ", nm)
}

coefficients <- c(
  primary_adjusted_High_vs_Low = "sarcosine_groupHigh",
  clinical_no_batch_High_vs_Low = "sarcosine_groupHigh",
  unadjusted_High_vs_Low = "sarcosine_groupHigh",
  primary_continuous_log2_Sarcosine = "sarcosine_z",
  stage_instead_of_grade_High_vs_Low = "sarcosine_groupHigh"
)

fit_limma <- function(expression_matrix, design, coefficient, model_label) {
  stopifnot(coefficient %in% colnames(design))
  fit <- eBayes(lmFit(expression_matrix, design), trend = TRUE, robust = TRUE)
  tt <- as.data.table(topTable(
    fit, coef = coefficient, number = Inf, sort.by = "none", adjust.method = "BH"
  ), keep.rownames = "gene_symbol")
  se <- fit$stdev.unscaled[, coefficient] * sqrt(fit$s2.post)
  crit <- qt(0.975, df = fit$df.total)
  if (length(crit) == 1L) crit <- rep(crit, nrow(tt))
  tt[, standard_error := se[gene_symbol]]
  tt[, `:=`(
    CI_low = logFC - crit * standard_error,
    CI_high = logFC + crit * standard_error,
    model = model_label
  )]
  setnames(
    tt,
    c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"),
    c("log2FC_or_SD_effect", "average_log_expression", "moderated_t", "p_value", "BH_q", "B_statistic")
  )
  stopifnot(nrow(tt) == nrow(expression_matrix), !anyNA(tt$moderated_t), all(is.finite(tt$moderated_t)))
  list(fit = fit, table = tt)
}

fits <- lapply(names(designs), function(nm) {
  fit_limma(log_expression, designs[[nm]], coefficients[[nm]], nm)
})
names(fits) <- names(designs)

for (nm in names(fits)) {
  fwrite(fits[[nm]]$table, file.path(table_dir, paste0("limma_", nm, "_all_genes.csv")))
}

# Exact unadjusted-coefficient reproduction check.
observed_difference <- rowMeans(log_expression[, meta$sarcosine_group == "High", drop = FALSE]) -
  rowMeans(log_expression[, meta$sarcosine_group == "Low", drop = FALSE])
unadjusted_lookup <- setNames(
  fits$unadjusted_High_vs_Low$table$log2FC_or_SD_effect,
  fits$unadjusted_High_vs_Low$table$gene_symbol
)
stopifnot(max(abs(observed_difference - unadjusted_lookup[names(observed_difference)])) < 1e-10)

model_summary <- rbindlist(lapply(names(fits), function(nm) {
  x <- fits[[nm]]$table
  data.table(
    model = nm,
    tested_genes = nrow(x),
    positive_BH_q_lt_0_05 = sum(x$BH_q < 0.05 & x$log2FC_or_SD_effect > 0),
    negative_BH_q_lt_0_05 = sum(x$BH_q < 0.05 & x$log2FC_or_SD_effect < 0),
    positive_BH_q_lt_0_05_abs_effect_ge_0_5 = sum(x$BH_q < 0.05 & x$log2FC_or_SD_effect >= 0.5),
    negative_BH_q_lt_0_05_abs_effect_ge_0_5 = sum(x$BH_q < 0.05 & x$log2FC_or_SD_effect <= -0.5)
  )
}))
fwrite(model_summary, file.path(table_dir, "limma_model_summary.csv"))

rank_from_table <- function(x) {
  stats <- x$moderated_t
  names(stats) <- x$gene_symbol
  stats <- sort(stats, decreasing = TRUE)
  stopifnot(!anyNA(stats), all(is.finite(stats)), !anyDuplicated(names(stats)))
  stats
}

rank_primary <- rank_from_table(fits$primary_adjusted_High_vs_Low$table)
rank_continuous <- rank_from_table(fits$primary_continuous_log2_Sarcosine$table)

hallmark_membership <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
reactome_membership <- as.data.table(msigdbr(
  species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"
))
gobp_membership <- as.data.table(msigdbr(
  species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"
))
stopifnot(uniqueN(c(
  hallmark_membership$db_version,
  reactome_membership$db_version,
  gobp_membership$db_version
)) == 1L)
msigdb_version <- unique(hallmark_membership$db_version)

make_pathways <- function(x) lapply(split(x$gene_symbol, x$gs_name), unique)
pathway_lists <- list(
  Hallmark = make_pathways(hallmark_membership),
  Reactome = make_pathways(reactome_membership),
  `GO:BP` = make_pathways(gobp_membership)
)
term_metadata <- rbindlist(list(
  unique(hallmark_membership[, .(
    collection = "Hallmark", pathway = gs_name, description = gs_description,
    exact_source = gs_exact_source, url = gs_url, MSigDB_version = db_version
  )]),
  unique(reactome_membership[, .(
    collection = "Reactome", pathway = gs_name, description = gs_description,
    exact_source = gs_exact_source, url = gs_url, MSigDB_version = db_version
  )]),
  unique(gobp_membership[, .(
    collection = "GO:BP", pathway = gs_name, description = gs_description,
    exact_source = gs_exact_source, url = gs_url, MSigDB_version = db_version
  )])
), use.names = TRUE)

run_fgsea <- function(pathways, stats, collection_label, model_label) {
  min_size <- if (collection_label == "Hallmark") 10L else 15L
  set.seed(260902)
  ans <- as.data.table(fgseaMultilevel(
    pathways = pathways,
    stats = stats,
    minSize = min_size,
    maxSize = 500L,
    eps = 0,
    scoreType = "std",
    nPermSimple = 100000,
    nproc = 1,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  ))
  ans[, `:=`(
    collection = collection_label,
    model = model_label,
    direction = fifelse(NES > 0, "Sarcosine-High", "Sarcosine-Low")
  )]
  merge(ans, term_metadata[collection == collection_label], by = c("pathway", "collection"), all.x = TRUE)
}

gsea_primary_raw <- rbindlist(lapply(names(pathway_lists), function(collection_label) {
  run_fgsea(pathway_lists[[collection_label]], rank_primary, collection_label, "Adjusted Sarcosine High vs Low")
}), use.names = TRUE)
gsea_continuous_raw <- rbindlist(lapply(names(pathway_lists), function(collection_label) {
  run_fgsea(pathway_lists[[collection_label]], rank_continuous, collection_label, "Per 1 SD higher log2 Sarcosine")
}), use.names = TRUE)

write_gsea <- function(x, filename) {
  out <- copy(x)
  out[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
  out[, abs_NES_sort := abs(NES)]
  setorder(out, padj, -abs_NES_sort, pathway)
  out[, abs_NES_sort := NULL]
  fwrite(out, file.path(table_dir, filename))
}
write_gsea(gsea_primary_raw, "GSEA_Hallmark_Reactome_GOBP_adjusted_High_vs_Low.csv")
write_gsea(gsea_continuous_raw, "GSEA_Hallmark_Reactome_GOBP_adjusted_continuous_log2_Sarcosine.csv")

gsea_summary <- rbindlist(lapply(list(
  adjusted_High_vs_Low = gsea_primary_raw,
  continuous_log2_Sarcosine = gsea_continuous_raw
), function(x) {
  x[, .(
    tested_pathways = .N,
    positive_BH_q_lt_0_05 = sum(padj < 0.05 & NES > 0, na.rm = TRUE),
    negative_BH_q_lt_0_05 = sum(padj < 0.05 & NES < 0, na.rm = TRUE),
    NA_pathways = sum(is.na(padj) | is.na(NES))
  ), by = collection]
}), idcol = "model")
fwrite(gsea_summary, file.path(table_dir, "GSEA_collection_summary.csv"))

focus_pathways <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "REACTOME_CO_STIMULATION_BY_CD28",
  "REACTOME_CD28_DEPENDENT_PI3K_AKT_SIGNALING",
  "REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING",
  "GOBP_INTERFERON_GAMMA_PRODUCTION",
  "GOBP_POSITIVE_REGULATION_OF_INTERFERON_GAMMA_PRODUCTION",
  "GOBP_RESPONSE_TO_INTERFERON_GAMMA",
  "GOBP_CELLULAR_RESPONSE_TO_INTERFERON_GAMMA",
  "GOBP_ANTIGEN_RECEPTOR_MEDIATED_SIGNALING_PATHWAY",
  "GOBP_T_CELL_ACTIVATION",
  "GOBP_T_CELL_PROLIFERATION",
  "GOBP_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_T_CELL_MEDIATED_CYTOTOXICITY",
  "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY",
  "GOBP_NATURAL_KILLER_CELL_ACTIVATION",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I"
)
focus_gsea <- gsea_primary_raw[pathway %in% focus_pathways]
focus_out <- copy(focus_gsea)
focus_out[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
focus_out[, pathway_order := match(pathway, focus_pathways)]
setorder(focus_out, pathway_order)
focus_out[, pathway_order := NULL]
fwrite(focus_out, file.path(table_dir, "GSEA_IFNG_CD28_IL12_immune_focus.csv"))

# Complete source membership for focused sets.
focus_membership <- rbindlist(lapply(names(pathway_lists), function(collection_label) {
  x <- if (collection_label == "Hallmark") hallmark_membership else if (collection_label == "Reactome") reactome_membership else gobp_membership
  unique(x[gs_name %in% focus_pathways, .(
    collection = collection_label, pathway = gs_name, gene_symbol,
    MSigDB_version = db_version, exact_source = gs_exact_source, url = gs_url
  )])
}), use.names = TRUE)
fwrite(focus_membership, file.path(table_dir, "GSEA_IFNG_CD28_IL12_focus_pathway_membership.csv"))

# -------------------------------------------------------------------------
# Figures
# -------------------------------------------------------------------------
COL_LOW <- "#1B9E8F"
COL_HIGH <- "#C43C3C"
COL_NS <- "#8F8F8F"

theme_publication <- theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 17),
    plot.subtitle = element_text(size = 11, colour = "#444444"),
    axis.title = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 12),
    legend.position = "bottom"
  )

primary_tt <- copy(fits$primary_adjusted_High_vs_Low$table)
primary_tt[, significance := fifelse(
  BH_q < 0.05 & log2FC_or_SD_effect > 0, "Sarcosine-High",
  fifelse(BH_q < 0.05 & log2FC_or_SD_effect < 0, "Sarcosine-Low", "BH q>=0.05")
)]
primary_tt[, minus_log10_q := -log10(pmax(BH_q, .Machine$double.xmin))]
label_genes <- unique(c(
  primary_tt[order(BH_q)][1:min(12L, .N), gene_symbol],
  c("SARDH", "PIPOX", "DMGDH", "GNMT", "CD28", "IFNG", "IL12RB1", "IL12RB2", "STAT4", "NFKB1", "RELA")
))
label_dt <- primary_tt[gene_symbol %in% label_genes]

p_volcano <- ggplot(primary_tt, aes(log2FC_or_SD_effect, minus_log10_q, colour = significance)) +
  geom_point(alpha = 0.55, size = 1.2) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#777777") +
  scale_colour_manual(values = c(
    "Sarcosine-Low" = COL_LOW, "BH q>=0.05" = "#C8C8C8", "Sarcosine-High" = COL_HIGH
  )) +
  geom_text(
    data = label_dt, aes(label = gene_symbol), check_overlap = TRUE,
    colour = "black", size = 3.2, vjust = -0.5
  ) +
  labs(
    title = "Tumour transcriptome associated with measured tissue Sarcosine",
    subtitle = "Sarcosine-High versus -Low; adjusted for RNA batch, sex, age and grade",
    x = "Adjusted High - Low expression difference", y = expression(-log[10]("BH q")), colour = NULL
  ) + theme_publication
ggsave(file.path(figure_dir, "Fig_Transcriptome_Volcano_Sarcosine_High_vs_Low.png"), p_volcano, width = 10, height = 7.5, dpi = 400)

sig_gsea <- gsea_primary_raw[padj < 0.05]
overview_gsea <- sig_gsea[order(padj, -abs(NES)), head(.SD, 6L), by = .(collection, direction)]
if (!nrow(overview_gsea)) overview_gsea <- gsea_primary_raw[order(padj, -abs(NES)), head(.SD, 10L), by = collection]
overview_gsea[, pathway_label := gsub("^(HALLMARK_|REACTOME_|GOBP_)", "", pathway)]
overview_gsea[, pathway_label := gsub("_", " ", pathway_label)]
overview_gsea[, minus_log10_q := pmin(-log10(pmax(padj, .Machine$double.xmin)), 50)]
overview_export <- copy(overview_gsea)
overview_export[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
fwrite(overview_export, file.path(table_dir, "Figure_Source_GSEA_Overview.csv"))

p_gsea <- ggplot(overview_gsea, aes(NES, reorder(pathway_label, NES), colour = direction, size = minus_log10_q)) +
  geom_vline(xintercept = 0, colour = "#888888") +
  geom_point(alpha = 0.9) +
  facet_wrap(~ collection, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = c("Sarcosine-Low" = COL_LOW, "Sarcosine-High" = COL_HIGH)) +
  scale_size_continuous(range = c(2.5, 8)) +
  labs(
    title = "Hallmark, Reactome and GO:BP programs by tissue Sarcosine state",
    subtitle = "Positive NES = Sarcosine-High; displayed terms have BH q<0.05 when available",
    x = "Normalized enrichment score (High vs Low)", y = NULL,
    colour = NULL, size = expression(-log[10]("BH q"))
  ) + theme_publication + theme(strip.background = element_blank())
ggsave(file.path(figure_dir, "Fig_ThreeCollection_GSEA_Sarcosine_High_vs_Low.png"), p_gsea, width = 12, height = 16, dpi = 400)

focus_plot <- copy(focus_gsea)
focus_plot[, pathway_label := gsub("^(HALLMARK_|REACTOME_|GOBP_)", "", pathway)]
focus_plot[, pathway_label := gsub("_", " ", pathway_label)]
focus_plot[, significance := fifelse(padj < 0.05, "BH q<0.05", "BH q>=0.05")]
focus_plot[, minus_log10_q := pmin(-log10(pmax(padj, .Machine$double.xmin)), 50)]
focus_plot[, direction := fifelse(NES > 0, "Sarcosine-High", "Sarcosine-Low")]

p_focus <- ggplot(focus_plot, aes(NES, factor(pathway_label, levels = rev(pathway_label)))) +
  geom_vline(xintercept = 0, colour = "#888888") +
  geom_point(aes(colour = direction, size = minus_log10_q, shape = significance), alpha = 0.9) +
  scale_colour_manual(values = c("Sarcosine-Low" = COL_LOW, "Sarcosine-High" = COL_HIGH)) +
  scale_shape_manual(values = c("BH q<0.05" = 16, "BH q>=0.05" = 1)) +
  scale_size_continuous(range = c(2.5, 8)) +
  labs(
    title = "IFN-gamma, CD28, IL-12 and cytotoxic immune programs",
    subtitle = "Tumour Sarcosine-High versus -Low; covariate-adjusted ranked transcriptome",
    x = "Normalized enrichment score (High vs Low)", y = NULL,
    colour = NULL, shape = NULL, size = expression(-log[10]("BH q"))
  ) + theme_publication
ggsave(file.path(figure_dir, "Fig_IFNG_CD28_IL12_Focused_GSEA.png"), p_focus, width = 13, height = 10.5, dpi = 400)

# Focused leading-edge sample heatmap. Up to 15 strongest primary leading-edge genes per pathway.
heatmap_targets <- intersect(c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "REACTOME_CO_STIMULATION_BY_CD28",
  "REACTOME_CD28_DEPENDENT_PI3K_AKT_SIGNALING",
  "REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING"
), gsea_primary_raw$pathway)
heatmap_members <- rbindlist(lapply(heatmap_targets, function(path) {
  genes <- gsea_primary_raw[pathway == path, leadingEdge][[1]]
  gene_effects <- primary_tt[gene_symbol %in% genes]
  gene_effects[, abs_t_sort := abs(moderated_t)]
  setorder(gene_effects, -abs_t_sort)
  gene_effects[, abs_t_sort := NULL]
  gene_effects[1:min(15L, .N), .(pathway = path, gene_symbol, moderated_t, BH_q)]
}), use.names = TRUE)
heatmap_genes <- unique(heatmap_members$gene_symbol)
stopifnot(length(heatmap_genes) >= 10L)
heatmap_expression <- log_expression[heatmap_genes, , drop = FALSE]
heatmap_z <- t(scale(t(heatmap_expression)))
heatmap_z[!is.finite(heatmap_z)] <- 0
sample_order <- meta[order(sarcosine_group, log2_sarcosine_intensity), sample_id]
heatmap_z <- heatmap_z[, sample_order, drop = FALSE]
annotation_col <- as.data.frame(meta[match(sample_order, sample_id), .(
  Sarcosine = sarcosine_group,
  `log2 Sarcosine signal` = log2_sarcosine_intensity
)])
rownames(annotation_col) <- sample_order
heatmap_membership_label <- heatmap_members[, .(
  pathway_membership = paste(gsub("^(HALLMARK_|REACTOME_)", "", pathway), collapse = ";")
), by = gene_symbol]
annotation_row <- data.frame(
  Focus = heatmap_membership_label$pathway_membership[match(rownames(heatmap_z), heatmap_membership_label$gene_symbol)],
  row.names = rownames(heatmap_z), check.names = FALSE
)
fwrite(
  as.data.table(heatmap_z, keep.rownames = "gene_symbol"),
  file.path(table_dir, "Figure_Source_IFNG_CD28_IL12_leading_edge_heatmap_zscores.csv")
)
fwrite(heatmap_members, file.path(table_dir, "Figure_Source_IFNG_CD28_IL12_heatmap_gene_selection.csv"))

png(
  file.path(figure_dir, "Fig_IFNG_CD28_IL12_LeadingEdge_Sample_Heatmap.png"),
  width = 4800, height = max(3400, 115 * nrow(heatmap_z)), res = 350
)
pheatmap(
  heatmap_z,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  show_colnames = FALSE,
  border_color = NA,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks = seq(-2.5, 2.5, length.out = 102),
  annotation_col = annotation_col,
  annotation_colors = list(Sarcosine = c(Low = COL_LOW, High = COL_HIGH)),
  gaps_col = 50,
  fontsize_row = 9,
  main = "Core genes of IFN-gamma/CD28/IL-12 programs across measured Sarcosine groups"
)
dev.off()

target_genes <- c(
  "SARDH", "PIPOX", "DMGDH", "GNMT", "GLDC", "GCSH", "AMT",
  "CD28", "CD80", "CD86", "ICOS", "LCK", "VAV1", "PIK3CG",
  "IFNG", "IFNGR1", "IFNGR2", "JAK1", "JAK2", "STAT1", "IRF1",
  "IL12A", "IL12B", "IL12RB1", "IL12RB2", "STAT4", "NFKB1", "RELA",
  "CD3D", "CD8A", "GZMB", "PRF1"
)
target_effects <- primary_tt[gene_symbol %in% target_genes]
target_effects[, family := fifelse(
  gene_symbol %in% c("SARDH", "PIPOX", "DMGDH", "GNMT", "GLDC", "GCSH", "AMT"), "Sarcosine metabolism",
  fifelse(gene_symbol %in% c("CD28", "CD80", "CD86", "ICOS", "LCK", "VAV1", "PIK3CG", "CD3D", "CD8A"), "CD28/T-cell",
  fifelse(gene_symbol %in% c("IL12A", "IL12B", "IL12RB1", "IL12RB2", "STAT4"), "IL-12 family", "IFN-gamma/cytotoxicity"))
)]
target_effects[, significance := fifelse(BH_q < 0.05, "BH q<0.05", "BH q>=0.05")]
setorder(target_effects, family, log2FC_or_SD_effect)
fwrite(target_effects, file.path(table_dir, "Focused_gene_effects_Sarcosine_IFNG_CD28_IL12.csv"))

p_genes <- ggplot(target_effects, aes(log2FC_or_SD_effect, reorder(gene_symbol, log2FC_or_SD_effect))) +
  geom_vline(xintercept = 0, colour = "#888888") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high, colour = significance), width = 0, orientation = "y") +
  geom_point(aes(colour = significance), size = 2.8) +
  facet_wrap(~ family, scales = "free_y", ncol = 2) +
  scale_colour_manual(values = c("BH q<0.05" = COL_HIGH, "BH q>=0.05" = COL_NS)) +
  labs(
    title = "Focused gene-level associations with measured tissue Sarcosine",
    subtitle = "Positive effect = higher expression in Sarcosine-High",
    x = "Adjusted High - Low expression difference (95% CI)", y = NULL, colour = NULL
  ) + theme_publication + theme(strip.background = element_blank())
ggsave(file.path(figure_dir, "Fig_Focused_Gene_Effects_Sarcosine_IFNG_CD28_IL12.png"), p_genes, width = 13, height = 11, dpi = 400)

# Analysis manifest and executable validation.
input_manifest <- data.table(
  input = c(expression_path, group_path),
  md5 = unname(tools::md5sum(c(expression_path, group_path))),
  size_bytes = file.info(c(expression_path, group_path))$size
)
fwrite(input_manifest, file.path(log_dir, "input_manifest.csv"))

validation <- data.table(
  check = c(
    "100 unique tumour samples", "fixed median cutpoint", "High n=50", "Low n=50",
    "no value equal to median", "RNA IDs aligned by name", "no missing expression",
    "more than 5000 genes after filter", "all five designs full rank",
    "unadjusted coefficient reproduces mean difference", "one MSigDB version",
    "three collections tested", "focused pathways found", "heatmap has >=10 genes",
    "all required figures exist"
  ),
  passed = c(
    nrow(meta) == 100L && uniqueN(meta$sample_id) == 100L,
    abs(cutpoint - 134.4850986) < 1e-8,
    sum(meta$sarcosine_group == "High") == 50L,
    sum(meta$sarcosine_group == "Low") == 50L,
    sum(meta$sarcosine_normalized_intensity == cutpoint) == 0L,
    identical(colnames(log_expression_all), meta$sample_id),
    !anyNA(log_expression),
    nrow(log_expression) > 5000L,
    all(vapply(designs, function(x) qr(x)$rank == ncol(x), logical(1))),
    max(abs(observed_difference - unadjusted_lookup[names(observed_difference)])) < 1e-10,
    length(msigdb_version) == 1L,
    setequal(unique(gsea_primary_raw$collection), c("Hallmark", "Reactome", "GO:BP")),
    nrow(focus_gsea) >= 10L,
    length(heatmap_genes) >= 10L,
    all(file.exists(file.path(figure_dir, c(
      "Fig_Transcriptome_Volcano_Sarcosine_High_vs_Low.png",
      "Fig_ThreeCollection_GSEA_Sarcosine_High_vs_Low.png",
      "Fig_IFNG_CD28_IL12_Focused_GSEA.png",
      "Fig_IFNG_CD28_IL12_LeadingEdge_Sample_Heatmap.png",
      "Fig_Focused_Gene_Effects_Sarcosine_IFNG_CD28_IL12.png"
    ))))
  )
)
fwrite(validation, file.path(log_dir, "validation_checks.csv"))
if (!all(validation$passed)) stop("One or more validation checks failed")

capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))
cat("Analysis complete: ", output_root, "\n", sep = "")
cat("Expression-filtered genes: ", nrow(log_expression), "\n", sep = "")
print(model_summary)
print(gsea_summary)
print(focus_out[, .(collection, pathway, NES, pval, padj, size)])
