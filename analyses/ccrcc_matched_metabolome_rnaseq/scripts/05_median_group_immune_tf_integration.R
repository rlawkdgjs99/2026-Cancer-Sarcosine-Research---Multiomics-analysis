#!/usr/bin/env Rscript

# Integrated immune-program, abundance-adjustment and CD28-linked TF analysis.
# Biological grouping displayed: measured tissue Sarcosine High versus Low.
# Positive effects mean higher in Sarcosine-High.

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

required_packages <- c("data.table", "limma", "msigdbr", "ggplot2", "pheatmap")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required R packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(msigdbr)
  library(ggplot2)
  library(pheatmap)
})

expression_path <- file.path(analysis_root, "RNAseq_Data", "bulkRNA_matrix_TPM.csv")
group_path <- file.path(analysis_root, "results", "00_input_audit_26.09.02", "tumor_sarcosine_group_map.csv")
deconv_root <- file.path(analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_HighLow_Deconvolution_26.09.02")
deconv_values_path <- file.path(deconv_root, "tables", "core_celltype_values_long.csv")
regulon_path <- file.path(
  project_root, "공공_BulkRNAseq_ICI반응성관련_분석모음",
  "Sarcosine-TIGER상_ICI반응성_RNAseq분석", "Melanoma-PRJEB23709",
  "results", "TIGER_PRE73_CD28_Upstream_TF_26.09.02", "tables",
  "01_DoRothEA_AB_signed_regulon_all.csv"
)
for (path in c(expression_path, group_path, deconv_values_path, regulon_path)) {
  if (!file.exists(path)) stop("Missing input: ", path)
}

output_root <- file.path(
  analysis_root, "results", "TJ_RCC_Tumor_Sarcosine_HighLow_Immune_TF_Integrated_26.09.02"
)
table_dir <- file.path(output_root, "tables")
figure_dir <- file.path(output_root, "figures")
log_dir <- file.path(output_root, "logs")
for (d in c(table_dir, figure_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

expr_dt <- fread(expression_path, check.names = FALSE)
meta <- fread(group_path, check.names = FALSE)
deconv <- fread(deconv_values_path, check.names = FALSE)
regulon <- fread(regulon_path, check.names = FALSE)
stopifnot(nrow(meta) == 100L, uniqueN(meta$sample_id) == 100L)
stopifnot(all(meta$sarcosine_group == ifelse(meta$sarcosine_normalized_intensity > 134.4850986, "High", "Low")))
meta[, sarcosine_group := factor(sarcosine_group, levels = c("Low", "High"))]
stopifnot(sum(meta$sarcosine_group == "Low") == 50L, sum(meta$sarcosine_group == "High") == 50L)

gene_symbols <- trimws(expr_dt[[1]])
sample_ids <- meta$sample_id
stopifnot(!anyDuplicated(gene_symbols), all(sample_ids %in% names(expr_dt)))
log_tpm_all <- as.matrix(expr_dt[, ..sample_ids])
storage.mode(log_tpm_all) <- "double"
rownames(log_tpm_all) <- gene_symbols
colnames(log_tpm_all) <- sample_ids
keep_expression <- rowSums(log_tpm_all >= 1) >= ceiling(0.10 * ncol(log_tpm_all))
log_tpm <- log_tpm_all[keep_expression, , drop = FALSE]
log_tpm <- log_tpm[apply(log_tpm, 1L, var) > 0, , drop = FALSE]
stopifnot(nrow(log_tpm) == 15119L, !anyNA(log_tpm), all(is.finite(log_tpm)))

meta[, batch := factor(batch)]
meta[, sex := relevel(factor(sex), ref = "female")]
meta[, age := factor(age, levels = c("40-60", "<40", ">60"))]
meta[, grade := factor(grade, levels = c("1", "2", "3", "4"))]
meta[, sarcosine_z := as.numeric(scale(log2_sarcosine_intensity))]

# -------------------------------------------------------------------------
# Frozen immune/cancer-relevant program definitions from MSigDB 2026.1.Hs.
# -------------------------------------------------------------------------
program_catalog <- data.table(
  program_id = c(
    "HALLMARK_INTERFERON_GAMMA_RESPONSE",
    "REACTOME_INTERFERON_GAMMA_SIGNALING",
    "REACTOME_CO_STIMULATION_BY_CD28",
    "REACTOME_CD28_DEPENDENT_PI3K_AKT_SIGNALING",
    "REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING",
    "GOBP_T_CELL_PROLIFERATION",
    "GOBP_T_CELL_DIFFERENTIATION_INVOLVED_IN_IMMUNE_RESPONSE",
    "GOBP_T_CELL_MEDIATED_CYTOTOXICITY",
    "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_I",
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_E2F_TARGETS"
  ),
  program_label = c(
    "IFN-gamma response", "IFN-gamma signaling", "CD28 co-stimulation",
    "CD28-PI3K/AKT", "IL-12-family signaling", "T-cell proliferation",
    "Effector-associated T-cell differentiation", "T-cell cytotoxicity",
    "MHC-I antigen presentation", "TNFalpha-NF-kB", "G2M checkpoint", "E2F targets"
  ),
  family = c(
    rep("IFN-gamma axis", 2), rep("CD28/IL-12 upstream", 3),
    rep("T-cell state", 4), "Inflammatory signaling", rep("Cancer-cell proliferation", 2)
  ),
  collection = c(
    "Hallmark", "Reactome", "Reactome", "Reactome", "Reactome",
    "GO:BP", "GO:BP", "GO:BP", "GO:BP", "Hallmark", "Hallmark", "Hallmark"
  ),
  display_order = 1:12
)

hallmark <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
reactome <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))
gobp <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"))
stopifnot(uniqueN(c(hallmark$db_version, reactome$db_version, gobp$db_version)) == 1L)
membership <- rbindlist(list(
  hallmark[gs_name %in% program_catalog$program_id, .(program_id = gs_name, gene_symbol, MSigDB_version = db_version, url = gs_url)],
  reactome[gs_name %in% program_catalog$program_id, .(program_id = gs_name, gene_symbol, MSigDB_version = db_version, url = gs_url)],
  gobp[gs_name %in% program_catalog$program_id, .(program_id = gs_name, gene_symbol, MSigDB_version = db_version, url = gs_url)]
), use.names = TRUE)
membership <- unique(membership, by = c("program_id", "gene_symbol"))
stopifnot(setequal(unique(membership$program_id), program_catalog$program_id))
membership[, expressed := gene_symbol %in% rownames(log_tpm)]

score_one <- function(genes) {
  genes <- intersect(unique(genes), rownames(log_tpm))
  if (length(genes) < 10L) stop("Program has fewer than 10 expressed genes")
  z <- t(scale(t(log_tpm[genes, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z)))
}
program_mat <- do.call(rbind, lapply(program_catalog$program_id, function(pid) {
  score_one(membership[program_id == pid, gene_symbol])
}))
rownames(program_mat) <- program_catalog$program_id
colnames(program_mat) <- sample_ids
stopifnot(!anyNA(program_mat), all(is.finite(program_mat)))

program_manifest <- merge(
  program_catalog,
  membership[, .(total_members = .N, expressed_members = sum(expressed)), by = program_id],
  by = "program_id", all.x = TRUE
)
setorder(program_manifest, display_order)
fwrite(program_manifest, file.path(table_dir, "immune_program_manifest.csv"))
fwrite(membership, file.path(table_dir, "immune_program_membership.csv"))

program_long <- as.data.table(as.table(program_mat))
setnames(program_long, c("program_id", "sample_id", "program_score_z"))
program_long <- merge(program_long, program_catalog, by = "program_id", all.x = TRUE)
program_long <- merge(
  program_long,
  meta[, .(sample_id, sarcosine_group, log2_sarcosine_intensity)],
  by = "sample_id", all.x = TRUE
)
fwrite(program_long, file.path(table_dir, "immune_program_scores_per_tumour.csv"))

# -------------------------------------------------------------------------
# Abundance and TME adjusters from the seven-method deconvolution.
# -------------------------------------------------------------------------
stopifnot(all(c("sample_id", "method", "cell_type", "analysis_value") %in% names(deconv)))
cd8_long <- deconv[cell_type == "CD8 T cells" & method %in% c(
  "CIBERSORT", "EPIC", "MCPcounter", "quanTIseq", "TIMER", "xCell"
), .(sample_id, method, analysis_value)]
stopifnot(nrow(cd8_long) == 600L)
cd8_wide <- dcast(cd8_long, sample_id ~ method, value.var = "analysis_value")
cd8_wide <- cd8_wide[match(sample_ids, sample_id)]
cd8_wide[, CD8_consensus_z := as.numeric(scale(apply(as.matrix(.SD), 1L, median))), .SDcols = setdiff(names(cd8_wide), "sample_id")]
cd8_wide[, CD8_consensus_no_CIBERSORT_z := as.numeric(scale(apply(
  as.matrix(.SD), 1L, median
))), .SDcols = setdiff(names(cd8_wide), c("sample_id", "CIBERSORT", "CD8_consensus_z"))]

pan_t <- deconv[method == "MCPcounter" & cell_type == "Pan T cells", .(sample_id, pan_T_z = analysis_value)]
estimate <- deconv[method == "ESTIMATE" & cell_type %in% c("Immune score", "Stromal score"), .(sample_id, cell_type, analysis_value)]
estimate_wide <- dcast(estimate, sample_id ~ cell_type, value.var = "analysis_value")
setnames(estimate_wide, c("Immune score", "Stromal score"), c("immune_score_z", "stromal_score_z"))

model_meta <- Reduce(function(x, y) merge(x, y, by = "sample_id", all.x = TRUE, sort = FALSE), list(
  meta,
  cd8_wide[, .(sample_id, CD8_consensus_z, CD8_consensus_no_CIBERSORT_z)],
  pan_t,
  estimate_wide
))
model_meta <- model_meta[match(sample_ids, sample_id)]
stopifnot(identical(model_meta$sample_id, sample_ids), !anyNA(model_meta[, .(
  CD8_consensus_z, CD8_consensus_no_CIBERSORT_z, pan_T_z, immune_score_z, stromal_score_z
)]))

model_formulas <- list(
  `High vs Low` = ~ batch + sex + age + grade + sarcosine_group,
  `Adjusted for CD8 consensus` = ~ batch + sex + age + grade + CD8_consensus_z + sarcosine_group,
  `Adjusted for CD8 consensus excluding CIBERSORT` = ~ batch + sex + age + grade + CD8_consensus_no_CIBERSORT_z + sarcosine_group,
  `Adjusted for pan-T estimate` = ~ batch + sex + age + grade + pan_T_z + sarcosine_group,
  `Adjusted for immune/stromal scores` = ~ batch + sex + age + grade + immune_score_z + stromal_score_z + sarcosine_group
)

fit_program_model <- function(label, formula) {
  design <- model.matrix(formula, data = model_meta)
  stopifnot(qr(design)$rank == ncol(design), "sarcosine_groupHigh" %in% colnames(design))
  fit <- eBayes(lmFit(program_mat, design), robust = TRUE)
  tt <- as.data.table(topTable(
    fit, coef = "sarcosine_groupHigh", number = Inf, sort.by = "none", adjust.method = "BH"
  ), keep.rownames = "program_id")
  se <- fit$stdev.unscaled[, "sarcosine_groupHigh"] * sqrt(fit$s2.post)
  crit <- qt(0.975, df = fit$df.total)
  if (length(crit) == 1L) crit <- rep(crit, nrow(tt))
  tt[, standard_error := se[program_id]]
  tt[, `:=`(CI_low = logFC - crit * standard_error, CI_high = logFC + crit * standard_error)]
  setnames(tt, c("logFC", "t", "P.Value", "adj.P.Val"), c("effect_High_minus_Low", "moderated_t", "p_value", "BH_q"))
  tt[, adjustment := label]
  merge(tt, program_catalog, by = "program_id", all.x = TRUE)
}
program_effects <- rbindlist(lapply(names(model_formulas), function(label) {
  fit_program_model(label, model_formulas[[label]])
}), use.names = TRUE)
program_effects[, adjustment := factor(adjustment, levels = names(model_formulas))]
setorder(program_effects, display_order, adjustment)
fwrite(program_effects, file.path(table_dir, "immune_program_High_vs_Low_abundance_adjusted_effects.csv"))

# -------------------------------------------------------------------------
# CD28-pathway-linked signed DoRothEA A/B TF-program inference.
# -------------------------------------------------------------------------
stopifnot(all(c("TF", "target", "mor", "confidence") %in% names(regulon)))
regulon <- unique(regulon[confidence %in% c("A", "B") & is.finite(mor), .(TF, target, mor, confidence)], by = c("TF", "target"))
regulon[, expressed := target %in% rownames(log_tpm)]
tf_coverage <- regulon[expressed == TRUE, .(
  expressed_targets_n = uniqueN(target),
  activating_targets_n = uniqueN(target[mor > 0]),
  repressing_targets_n = uniqueN(target[mor < 0])
), by = TF]
eligible_tfs <- tf_coverage[expressed_targets_n >= 10L, TF]
regulon_tested <- regulon[expressed == TRUE & TF %in% eligible_tfs]
stopifnot(length(eligible_tfs) > 50L)

ulm_from_matrix <- function(mat, net, minsize = 10L) {
  shared_targets <- sort(intersect(rownames(mat), unique(net$target)))
  net <- net[target %in% shared_targets]
  keep_tf <- net[, uniqueN(target), by = TF][V1 >= minsize, TF]
  net <- net[TF %in% keep_tf]
  shared_targets <- sort(unique(net$target))
  tfs <- sort(unique(net$TF))
  mor_mat <- matrix(0, nrow = length(shared_targets), ncol = length(tfs), dimnames = list(shared_targets, tfs))
  mor_mat[cbind(match(net$target, shared_targets), match(net$TF, tfs))] <- net$mor
  y <- mat[shared_targets, , drop = FALSE]
  y <- sweep(y, 1L, rowMeans(y), FUN = "-")
  r <- cor(mor_mat, y)
  r[!is.finite(r)] <- NA_real_
  r <- pmax(pmin(r, 1 - 1e-12), -1 + 1e-12)
  df <- nrow(mor_mat) - 2L
  score <- r * sqrt(df / ((1 - r + 1e-20) * (1 + r + 1e-20)))
  score
}

tf_activity <- ulm_from_matrix(log_tpm, regulon_tested, minsize = 10L)
stopifnot(identical(colnames(tf_activity), sample_ids), !anyNA(tf_activity), all(is.finite(tf_activity)))
design_tf <- model.matrix(~ batch + sex + age + grade + sarcosine_group, data = meta)
fit_tf <- eBayes(lmFit(tf_activity, design_tf), robust = TRUE)
tf_tt <- as.data.table(topTable(
  fit_tf, coef = "sarcosine_groupHigh", number = Inf, sort.by = "none", adjust.method = "BH"
), keep.rownames = "TF")
tf_sd <- apply(tf_activity, 1L, sd)
tf_se <- fit_tf$stdev.unscaled[, "sarcosine_groupHigh"] * sqrt(fit_tf$s2.post)
tf_crit <- qt(0.975, df = fit_tf$df.total)
if (length(tf_crit) == 1L) tf_crit <- rep(tf_crit, nrow(tf_tt))
tf_tt[, standard_error := tf_se[TF]]
tf_tt[, `:=`(CI_low = logFC - tf_crit * standard_error, CI_high = logFC + tf_crit * standard_error)]
tf_tt[, `:=`(
  standardized_effect = logFC / tf_sd[TF],
  standardized_CI_low = CI_low / tf_sd[TF],
  standardized_CI_high = CI_high / tf_sd[TF]
)]
setnames(tf_tt, c("logFC", "t", "P.Value", "adj.P.Val"), c("activity_difference", "moderated_t", "p_value", "BH_q"))
tf_tt <- merge(tf_tt, tf_coverage, by = "TF", all.x = TRUE)

cd28_programs <- c("REACTOME_CO_STIMULATION_BY_CD28", "REACTOME_CD28_DEPENDENT_PI3K_AKT_SIGNALING")
cd28_genes <- unique(membership[program_id %in% cd28_programs, gene_symbol])
tf_overlap <- regulon_tested[target %in% cd28_genes, .(
  CD28_pathway_targets_n = uniqueN(target),
  CD28_pathway_targets = paste(sort(unique(target)), collapse = ";")
), by = TF]
tf_tt <- merge(tf_tt, tf_overlap, by = "TF", all.x = TRUE)
tf_tt[is.na(CD28_pathway_targets_n), `:=`(CD28_pathway_targets_n = 0L, CD28_pathway_targets = "")]
tf_tt[, direct_CD28_gene_regulator := TF %in% regulon[target == "CD28", TF]]
tf_tt[, abs_effect_sort := abs(standardized_effect)]
setorder(tf_tt, BH_q, -abs_effect_sort)
tf_tt[, abs_effect_sort := NULL]
fwrite(regulon_tested, file.path(table_dir, "DoRothEA_AB_signed_regulon_tested.csv"))
fwrite(tf_tt, file.path(table_dir, "DoRothEA_TF_activity_all_eligible_High_vs_Low.csv"))
fwrite(
  tf_tt[CD28_pathway_targets_n > 0 | direct_CD28_gene_regulator == TRUE],
  file.path(table_dir, "CD28_linked_TF_activity_High_vs_Low.csv")
)

tf_long <- as.data.table(as.table(tf_activity))
setnames(tf_long, c("TF", "sample_id", "ULM_activity"))
tf_long <- merge(tf_long, meta[, .(sample_id, sarcosine_group)], by = "sample_id", all.x = TRUE)
fwrite(tf_long, file.path(table_dir, "sample_level_TF_activities_long.csv"))

# -------------------------------------------------------------------------
# Reader-facing figures.
# -------------------------------------------------------------------------
COL_LOW <- "#1B9E8F"
COL_HIGH <- "#C43C3C"
COL_NS <- "#8F8F8F"
theme_publication <- theme_classic(base_size = 14) + theme(
  plot.title = element_text(face = "bold", size = 17),
  plot.subtitle = element_text(size = 11, colour = "#444444"),
  axis.title = element_text(face = "bold"),
  strip.text = element_text(face = "bold"),
  legend.position = "bottom"
)

main_programs <- program_catalog[display_order <= 10, program_id]
program_plot <- program_long[program_id %in% main_programs]
program_plot[, program_label := factor(program_label, levels = program_catalog[program_id %in% main_programs, program_label])]
p_scores <- ggplot(program_plot, aes(sarcosine_group, program_score_z, colour = sarcosine_group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "black", fill = "white") +
  geom_jitter(width = 0.12, alpha = 0.45, size = 1.0) +
  facet_wrap(~ program_label, scales = "free_y", ncol = 3) +
  scale_colour_manual(values = c(Low = COL_LOW, High = COL_HIGH)) +
  labs(
    title = "Immune programs across measured tissue Sarcosine groups",
    subtitle = "Every panel shows the same 50 Sarcosine-Low and 50 Sarcosine-High tumours",
    x = "Tissue Sarcosine group", y = "Pathway score (z)", colour = NULL
  ) + theme_publication
ggsave(file.path(figure_dir, "Fig_Immune_Program_Scores_Sarcosine_High_vs_Low.png"), p_scores, width = 15, height = 12, dpi = 400)

effect_plot <- program_effects[program_id %in% main_programs]
effect_plot[, significant := fifelse(BH_q < 0.05, "BH q<0.05", "BH q>=0.05")]
p_adjust <- ggplot(effect_plot, aes(effect_High_minus_Low, reorder(program_label, effect_High_minus_Low))) +
  geom_vline(xintercept = 0, colour = "#777777") +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high, colour = significant), width = 0, orientation = "y") +
  geom_point(aes(colour = significant), size = 2.5) +
  facet_wrap(~ adjustment, ncol = 2) +
  scale_colour_manual(values = c("BH q<0.05" = COL_HIGH, "BH q>=0.05" = COL_NS)) +
  labs(
    title = "Do immune-program differences persist after T-cell/TME adjustment?",
    subtitle = "Positive effect = higher in Sarcosine-High; adjustments are sensitivity analyses",
    x = "Adjusted High - Low difference (program-score SD; 95% CI)", y = NULL, colour = NULL
  ) + theme_publication
ggsave(file.path(figure_dir, "Fig_Immune_Programs_Tcell_TME_Adjustment.png"), p_adjust, width = 15, height = 14, dpi = 400)

candidate_tf <- tf_tt[CD28_pathway_targets_n > 0 | direct_CD28_gene_regulator == TRUE]
candidate_tf <- unique(rbind(
  candidate_tf[order(BH_q, -abs(standardized_effect))][1:min(20L, .N)],
  tf_tt[TF %in% c("NFKB1", "RELA", "SPI1", "EGR1", "SP1", "NFATC1", "JUN", "FOS")]
), by = "TF")
candidate_tf[, significant := fifelse(BH_q < 0.05, "BH q<0.05", "BH q>=0.05")]
setorder(candidate_tf, standardized_effect)
fwrite(candidate_tf, file.path(table_dir, "Figure_Source_CD28_candidate_TFs.csv"))
p_tf <- ggplot(candidate_tf, aes(standardized_effect, reorder(TF, standardized_effect))) +
  geom_vline(xintercept = 0, colour = "#777777") +
  geom_errorbar(aes(xmin = standardized_CI_low, xmax = standardized_CI_high, colour = significant), width = 0, orientation = "y") +
  geom_point(aes(colour = significant, size = pmax(CD28_pathway_targets_n, 1L))) +
  scale_colour_manual(values = c("BH q<0.05" = COL_HIGH, "BH q>=0.05" = COL_NS)) +
  scale_size_continuous(range = c(2.5, 7)) +
  labs(
    title = "Candidate TF programs linked to the CD28 pathway",
    subtitle = "Signed DoRothEA A/B regulon activity; positive effect = higher in Sarcosine-High",
    x = "Adjusted High - Low TF-activity difference (SD; 95% CI)", y = NULL,
    colour = NULL, size = "CD28-pathway\ntargets"
  ) + theme_publication
ggsave(file.path(figure_dir, "Fig_CD28_Linked_TF_Activity_Sarcosine_High_vs_Low.png"), p_tf, width = 10, height = 10, dpi = 400)

tf_heatmap_names <- candidate_tf[order(BH_q, -abs(standardized_effect))][1:min(15L, .N), TF]
tf_heat <- tf_activity[tf_heatmap_names, sample_ids, drop = FALSE]
tf_heat <- t(scale(t(tf_heat)))
tf_heat[!is.finite(tf_heat)] <- 0
sample_order <- meta[order(sarcosine_group, log2_sarcosine_intensity), sample_id]
tf_heat <- tf_heat[, sample_order, drop = FALSE]
ann_col <- data.frame(
  Sarcosine = meta[match(sample_order, sample_id), sarcosine_group],
  row.names = sample_order
)
png(file.path(figure_dir, "Fig_CD28_Linked_TF_Activity_Heatmap.png"), width = 4200, height = 2500, res = 350)
pheatmap(
  tf_heat, cluster_rows = TRUE, cluster_cols = FALSE, show_colnames = FALSE,
  border_color = NA, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(101),
  breaks = seq(-2.5, 2.5, length.out = 102), annotation_col = ann_col,
  annotation_colors = list(Sarcosine = c(Low = COL_LOW, High = COL_HIGH)),
  gaps_col = 50, main = "CD28-linked TF-program activity across measured Sarcosine groups"
)
dev.off()

input_manifest <- data.table(
  input = c(expression_path, group_path, deconv_values_path, regulon_path),
  md5 = unname(tools::md5sum(c(expression_path, group_path, deconv_values_path, regulon_path))),
  size_bytes = file.info(c(expression_path, group_path, deconv_values_path, regulon_path))$size
)
fwrite(input_manifest, file.path(log_dir, "input_manifest.csv"))

validation <- data.table(
  check = c(
    "100 unique tumours", "High n=50", "Low n=50", "15119 filtered genes",
    "all 12 programs found", "all program scores finite", "six CD8 estimates per sample",
    "all five program models full rank", ">50 eligible TFs", "all TF activities finite",
    "CD28-linked candidates retained", "all four figures exist"
  ),
  passed = c(
    nrow(meta) == 100L && uniqueN(meta$sample_id) == 100L,
    sum(meta$sarcosine_group == "High") == 50L,
    sum(meta$sarcosine_group == "Low") == 50L,
    nrow(log_tpm) == 15119L,
    setequal(unique(membership$program_id), program_catalog$program_id),
    !anyNA(program_mat) && all(is.finite(program_mat)),
    nrow(cd8_long) == 600L,
    all(vapply(model_formulas, function(f) {d <- model.matrix(f, data = model_meta); qr(d)$rank == ncol(d)}, logical(1))),
    length(eligible_tfs) > 50L,
    !anyNA(tf_activity) && all(is.finite(tf_activity)),
    nrow(candidate_tf) >= 10L,
    all(file.exists(file.path(figure_dir, c(
      "Fig_Immune_Program_Scores_Sarcosine_High_vs_Low.png",
      "Fig_Immune_Programs_Tcell_TME_Adjustment.png",
      "Fig_CD28_Linked_TF_Activity_Sarcosine_High_vs_Low.png",
      "Fig_CD28_Linked_TF_Activity_Heatmap.png"
    ))))
  )
)
fwrite(validation, file.path(log_dir, "validation_checks.csv"))
if (!all(validation$passed)) stop("One or more validation checks failed")
capture.output(sessionInfo(), file = file.path(log_dir, "sessionInfo.txt"))

cat("Integrated immune/TF analysis complete: ", output_root, "\n", sep = "")
print(program_effects[, .(program_label, adjustment, effect_High_minus_Low, CI_low, CI_high, p_value, BH_q)])
print(candidate_tf[, .(TF, standardized_effect, standardized_CI_low, standardized_CI_high, p_value, BH_q, CD28_pathway_targets_n, direct_CD28_gene_regulator)])
