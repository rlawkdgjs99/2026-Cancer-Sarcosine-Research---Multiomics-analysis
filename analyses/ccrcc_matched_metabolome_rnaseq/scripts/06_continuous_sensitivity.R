#!/usr/bin/env Rscript

# Influence/robustness audit for the continuous tissue-Sarcosine sensitivity.
# The reader-facing primary comparison remains the frozen median split (50 Low/50 High).

options(stringsAsFactors = FALSE, width = 180)
set.seed(260902)

analysis_root <- normalizePath(getwd())
project_root <- normalizePath(file.path(analysis_root, ".."))
gsea_lib <- file.path(
  project_root, "2024_Drug_Res_Updates_NSCLC", "RNA-seq공공데이터_GSE207422",
  "analysis_sarcosine_FINAL_26.08.25", "03_lee_fig3_style_FINAL",
  "hallmark_GSEA_Q4_vs_Q1_26.08.26", "R_libs"
)
.libPaths(unique(c(normalizePath(gsea_lib), .libPaths())))

required <- c("data.table", "limma", "fgsea", "msigdbr", "ggplot2", "BiocParallel")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table); library(limma); library(fgsea); library(msigdbr); library(ggplot2)
})

expression_path <- file.path(analysis_root, "RNAseq_Data", "bulkRNA_matrix_TPM.csv")
group_path <- file.path(analysis_root, "results", "00_input_audit_26.09.02", "tumor_sarcosine_group_map.csv")
out_root <- file.path(analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_Continuous_Influence_Sensitivity_26.09.02")
table_dir <- file.path(out_root, "tables"); figure_dir <- file.path(out_root, "figures"); log_dir <- file.path(out_root, "logs")
for (d in c(table_dir, figure_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

expr_dt <- fread(expression_path, check.names = FALSE)
meta <- fread(group_path, check.names = FALSE)
gene <- trimws(expr_dt[[1]])
sample_ids <- meta$sample_id
expr <- as.matrix(expr_dt[, ..sample_ids]); storage.mode(expr) <- "double"
rownames(expr) <- gene; colnames(expr) <- meta$sample_id
keep <- rowSums(expr >= 1) >= 10 & apply(expr, 1L, var) > 0
expr <- expr[keep, , drop = FALSE]
stopifnot(nrow(expr) == 15119L, identical(colnames(expr), meta$sample_id))

prepare_meta <- function(x) {
  x[, batch := factor(batch)]
  x[, sex := relevel(factor(sex), ref = "female")]
  x[, age := factor(age, levels = c("40-60", "<40", ">60"))]
  x[, grade := factor(grade, levels = c("1", "2", "3", "4"))]
  x
}
meta <- prepare_meta(meta)

q1 <- quantile(meta$sarcosine_normalized_intensity, 0.25)
q3 <- quantile(meta$sarcosine_normalized_intensity, 0.75)
tukey_low <- q1 - 1.5 * IQR(meta$sarcosine_normalized_intensity)
outlier_ids <- meta[sarcosine_normalized_intensity < tukey_low, sample_id]
stopifnot(identical(sort(outlier_ids), sort(c("H46_T", "N39_T", "R25_T", "R82_T", "Z16_T"))))

models <- list()
meta_all <- copy(meta)
meta_all[, exposure := as.numeric(scale(log2_sarcosine_intensity))]
models[["All tumours"]] <- list(meta = meta_all, expr = expr)

meta_no <- copy(meta[!sample_id %in% outlier_ids])
meta_no[, exposure := as.numeric(scale(log2_sarcosine_intensity))]
models[["Exclude 5 Tukey-low outliers"]] <- list(meta = meta_no, expr = expr[, meta_no$sample_id, drop = FALSE])

meta_win <- copy(meta)
win_lo <- quantile(meta_win$log2_sarcosine_intensity, 0.05)
win_hi <- quantile(meta_win$log2_sarcosine_intensity, 0.95)
meta_win[, log2_winsorized := pmin(pmax(log2_sarcosine_intensity, win_lo), win_hi)]
meta_win[, exposure := as.numeric(scale(log2_winsorized))]
models[["Winsorize log2 Sarcosine at 5/95%"]] <- list(meta = meta_win, expr = expr)

fit_one <- function(item, label) {
  m <- item$meta
  d <- model.matrix(~ batch + sex + age + grade + exposure, data = m)
  if (qr(d)$rank != ncol(d)) stop("Rank-deficient design: ", label)
  fit <- eBayes(lmFit(item$expr, d), trend = TRUE, robust = TRUE)
  tt <- as.data.table(topTable(fit, coef = "exposure", number = Inf, sort.by = "none", adjust.method = "BH"), keep.rownames = "gene_symbol")
  setnames(tt, c("logFC", "t", "P.Value", "adj.P.Val"), c("effect_per_SD", "moderated_t", "p_value", "BH_q"))
  stopifnot(nrow(tt) == 15119L, all(is.finite(tt$moderated_t)))
  tt[, model := label]
  list(table = tt, rank = setNames(tt$moderated_t, tt$gene_symbol), design = d)
}
fits <- Map(fit_one, models, names(models)); names(fits) <- names(models)

gene_summary <- rbindlist(lapply(names(fits), function(nm) {
  x <- fits[[nm]]$table
  data.table(
    model = nm, n = nrow(models[[nm]]$meta),
    positive_BH_q_lt_0_05 = sum(x$BH_q < 0.05 & x$effect_per_SD > 0),
    negative_BH_q_lt_0_05 = sum(x$BH_q < 0.05 & x$effect_per_SD < 0)
  )
}))
fwrite(gene_summary, file.path(table_dir, "continuous_gene_level_sensitivity_summary.csv"))

msig <- as.data.table(msigdbr(species = "Homo sapiens", db_species = "HS"))
focus <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "REACTOME_INTERFERON_GAMMA_SIGNALING", "REACTOME_CO_STIMULATION_BY_CD28",
  "REACTOME_CD28_DEPENDENT_PI3K_AKT_SIGNALING", "REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING",
  "GOBP_T_CELL_PROLIFERATION", "GOBP_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_T_CELL_MEDIATED_CYTOTOXICITY", "GOBP_LEUKOCYTE_MEDIATED_CYTOTOXICITY"
)
msig_focus <- unique(msig[gs_name %in% focus, .(gs_name, gene_symbol)])
pathways <- split(msig_focus$gene_symbol, msig_focus$gs_name)
stopifnot(setequal(names(pathways), focus))

run_focus <- function(stats, label) {
  ans <- as.data.table(fgseaMultilevel(
    pathways = pathways, stats = stats, minSize = 10L, maxSize = 500L,
    eps = 0, nPermSimple = 100000L, nproc = 1,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  ))
  ans[, `:=`(model = label, direction = fifelse(NES > 0, "Higher Sarcosine", "Lower Sarcosine"))]
  ans
}
gsea <- rbindlist(Map(function(x, nm) run_focus(x$rank, nm), fits, names(fits)), use.names = TRUE)
gsea_out <- copy(gsea); gsea_out[, leadingEdge := vapply(leadingEdge, paste, character(1), collapse = ";")]
setorder(gsea_out, model, padj, pathway)
fwrite(gsea_out, file.path(table_dir, "focused_pathway_continuous_influence_sensitivity.csv"))
fwrite(
  gsea_out[is.na(NES) | is.na(padj), .(model, pathway, reason = "fgsea unbalanced gene-level statistic; NES/p/q not estimable")],
  file.path(table_dir, "focused_pathway_unbalanced_NA_audit.csv")
)

manifest <- data.table(
  analysis = names(models),
  n = vapply(models, function(x) nrow(x$meta), integer(1)),
  exposure = c("z(log2 Sarcosine)", "z(log2 Sarcosine), Tukey-low outliers removed", "z(winsorized log2 Sarcosine)"),
  removed_or_capped = c("none", paste(outlier_ids, collapse = ";"), sprintf("5th=%.4f;95th=%.4f", win_lo, win_hi))
)
fwrite(manifest, file.path(table_dir, "continuous_sensitivity_model_manifest.csv"))

plot_dt <- copy(gsea)
plot_dt[, label := gsub("^(HALLMARK_|REACTOME_|GOBP_)", "", pathway)]
plot_dt[, label := gsub("_", " ", label)]
label_levels <- plot_dt[, .(mean_NES = mean(NES, na.rm = TRUE)), by = label][order(mean_NES), label]
plot_dt[, label := factor(label, levels = label_levels)]
plot_dt[, plot_NES := fifelse(is.finite(NES), NES, 0)]
plot_dt[, significance := fifelse(!is.finite(padj), "Not estimable", fifelse(padj < 0.05, "BH q<0.05", "BH q>=0.05"))]
plot_dt[, plot_direction := fifelse(!is.finite(NES), "Not estimable", direction)]
p <- ggplot(plot_dt, aes(plot_NES, label, colour = plot_direction, shape = significance)) +
  geom_vline(xintercept = 0, colour = "#777777") + geom_point(size = 3) +
  facet_wrap(~ model, ncol = 1) +
  scale_colour_manual(values = c("Higher Sarcosine" = "#C43C3C", "Lower Sarcosine" = "#1B9E8F", "Not estimable" = "#777777")) +
  scale_shape_manual(values = c("BH q<0.05" = 16, "BH q>=0.05" = 1, "Not estimable" = 4)) +
  labs(
    title = "Continuous tissue-Sarcosine pathway robustness",
    subtitle = "Positive NES = higher with Sarcosine; adjusted for batch, sex, age and grade; unestimable unbalanced sets omitted",
    x = "Normalized enrichment score per higher Sarcosine", y = NULL, colour = NULL, shape = NULL
  ) + theme_classic(base_size = 13) + theme(
    plot.title = element_text(face = "bold", size = 17), strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )
ggsave(file.path(figure_dir, "Fig_Continuous_Sarcosine_Influence_Sensitivity.png"), p, width = 13, height = 14, dpi = 400)

validation <- data.table(
  check = c("three prespecified models", "five Tukey-low outliers", "all designs full rank", "all focused pathways returned", "robust-treatment GSEA finite", "unbalanced NA pathways audited", "figure exists"),
  passed = c(
    length(models) == 3L, length(outlier_ids) == 5L,
    all(vapply(fits, function(x) qr(x$design)$rank == ncol(x$design), logical(1))),
    nrow(gsea) == 3L * length(focus),
    gsea[model != "All tumours", all(is.finite(NES)) && all(is.finite(padj))],
    gsea[is.na(NES) | is.na(padj), .N] == 4L && gsea[is.na(NES) | is.na(padj), all(model == "All tumours")],
    file.exists(file.path(figure_dir, "Fig_Continuous_Sarcosine_Influence_Sensitivity.png"))
  )
)
fwrite(validation, file.path(log_dir, "validation_checks.csv"))
if (!all(validation$passed)) stop("Continuous sensitivity validation failed")
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))
fwrite(data.table(input = c(expression_path, group_path), md5 = unname(tools::md5sum(c(expression_path, group_path)))), file.path(log_dir, "input_manifest.csv"))

cat("Continuous influence sensitivity complete: ", out_root, "\n", sep = "")
print(gene_summary)
print(gsea[, .(model, pathway, NES, padj, direction)])
