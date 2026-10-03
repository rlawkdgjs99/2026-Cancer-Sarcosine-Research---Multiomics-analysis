#!/usr/bin/env Rscript
# =============================================================================
# Figure 1b companion — ROC + median-split stacked bar for faecal sarcosine
#   CRC vs healthy controls, MetaboLights MTBLS10232
#
# Modelled on Lee et al., Drug Resist Updat 77:101159 (2024), Fig 1B / 1C:
#   1B = ROC on the CONTINUOUS metabolite level (AUC)
#   1C = stacked bars of Low/High groups split at the metabolite MEDIAN,
#        tested by chi-square.
# Here the contrast is DISEASE STATUS (CRC vs healthy), not ICI response.
#
# Run with the working directory set to CRC_Metabolomics_MTBLS10232_260428/
#   Rscript R_scripts/101_ROC_stackedbar_26.07.27.R
# Existing analyses and figures are read-only; output goes to a NEW folder.
#
# DESIGN
#   Unit of observation : one stool sample = one subject; groups are INDEPENDENT
#   n                   : 246 CTRL vs 311 CRC (557 total)
#   Values              : log2 quantile-normalized NEG-mode LC-MS intensity,
#                         re-derived here from analysis/intensity_matrices.rds
#                         exactly as in 99_publication_style_plots_26.07.23.R.
#                         Already log-transformed -> no further transformation.
#   ROC                 : direction FIXED A PRIORI (CRC > CTRL). pROC's default
#                         direction="auto" forces AUC >= 0.5 and would bias it.
#   AUC CI              : DeLong — valid, because observations are independent.
#   Split               : POOLED median, per Lee 2024 ("the median value of each
#                         metabolite"). A within-group median would give 50:50
#                         by construction and be meaningless.
#   2x2 test            : Pearson chi-square (independent samples), with Fisher
#                         exact reported alongside.
#   Covariate           : the upstream limma model adjusted for age_cohort_bin,
#                         so an age-stratified AUC is reported as a sensitivity
#                         analysis. The headline AUC is crude, matching Lee 2024.
#
# CAVEATS carried from the upstream analysis (do not drop when reporting):
#   - Sarcosine here is a PUTATIVE annotation (MSI level >= 2), isobaric with
#     alanine/beta-alanine, no MS/MS deposited.
#   - Its age-adjusted log2FC (+0.43) does NOT meet this pipeline's own
#     discovery effect-size cutoff of |log2FC| > 0.5; it is highly significant
#     by p-value but modest in magnitude. The AUC reflects that.
#
# NOTE ON FILE I/O: this project lives under a Korean-named directory. macOS
# stores those names in NFD while R emits NFC, so readr/readxl cannot open
# files here even via relative paths. Base-R I/O and ggsave DO work, so this
# script deliberately uses read.delim/readRDS/write.csv/ggsave only.
# =============================================================================

suppressPackageStartupMessages({
  library(limma); library(pROC); library(ggplot2)
})
set.seed(42)

out_dir <- file.path("analysis", "main", "ROC_stackedbar_26.07.27")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

COL_LOW <- "#0072B2"; COL_HIGH <- "#E69F00"   # project Low/High standard
font_family <- "sans"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  arial <- systemfonts::match_fonts("Arial")
  if (nrow(arial) == 1L && nzchar(arial$path[1]) && file.exists(arial$path[1])) {
    font_family <- "Arial"
  }
}

# Saving PNGs from inside this Korean-named directory: the ragg device (which
# ggsave uses by default, and which 99_publication_style_plots_26.07.23.R
# requests explicitly) FAILS here and only emits a warning -- it returns
# silently with NO file written. grDevices::png works. save_png() therefore
# pins the device AND asserts the file really appeared, so a silent write
# failure becomes a hard error instead of a missing figure.
save_png <- function(path, plot, width, height, dpi = 300) {
  suppressWarnings(ggsave(path, plot, width = width, height = height,
                          dpi = dpi, bg = "white", device = grDevices::png))
  if (!file.exists(path) || file.size(path) < 10000)
    stop("PNG was not written (or is truncated): ", path)
  invisible(path)
}

# Publication theme — identical sizing to 99_publication_style_plots_26.07.23.R
# so these panels sit alongside the existing figures without a font mismatch.
theme_publication <- function() {
  theme_classic(base_family = font_family, base_size = 15) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5,
                                margin = margin(b = 8)),
      axis.title = element_text(size = 17, colour = "black"),
      axis.text  = element_text(size = 15, colour = "black"),
      axis.line  = element_line(linewidth = 0.65, colour = "black"),
      axis.ticks = element_line(linewidth = 0.55, colour = "black"),
      legend.position = "none",
      plot.margin = margin(12, 14, 10, 12)
    )
}
# The stacked bar needs its Low/High key, so it re-enables the legend at the
# same 15 pt as the axis text.
theme_publication_legend <- function() {
  theme_publication() +
    theme(legend.position = "right",
          legend.title = element_blank(),
          legend.text = element_text(size = 15, colour = "black"),
          legend.key.size = unit(0.9, "lines"))
}
ANNOT_SIZE <- 5.2   # ggplot size units (mm) ~= 15 pt, matching axis.text

log_lines <- character(0)
say <- function(...) { s <- sprintf(...); cat(s, "\n"); log_lines <<- c(log_lines, s) }

# -----------------------------------------------------------------------------
# 1. Re-derive the sarcosine vector from source (same steps as 99_ script)
# -----------------------------------------------------------------------------
dat      <- readRDS("analysis/intensity_matrices.rds")
metadata <- read.delim("metadata_merged_619samples.tsv",
                       check.names = FALSE, stringsAsFactors = FALSE)
diff_neg <- read.delim("analysis/main/diff_results_NEG.tsv",
                       stringsAsFactors = FALSE, check.names = FALSE)

stopifnot(is.matrix(dat$intensity_neg),
          identical(colnames(dat$intensity_neg), dat$sample_cols_neg),
          !anyDuplicated(metadata$MTBLS_Sample_Name))

sarc_idx <- which(dat$feature_meta_neg$metabolite_identification == "Sarcosine" &
                    dat$feature_meta_neg$database_identifier == "HMDB0000271")
if (length(sarc_idx) != 1L)
  stop("Expected exactly one NEG-mode Sarcosine/HMDB0000271 feature; found ", length(sarc_idx))

imp <- dat$intensity_neg
pos <- imp[is.finite(imp) & imp > 0]
if (!length(pos)) stop("NEG intensity matrix has no positive values.")
imp[imp == 0] <- min(pos) / 2
qn_neg <- normalizeBetweenArrays(log2(imp), method = "quantile")

meta_idx <- match(colnames(qn_neg), metadata$MTBLS_Sample_Name)
if (anyNA(meta_idx)) stop("At least one NEG-mode sample is absent from metadata.")
sc_meta <- metadata[meta_idx, , drop = FALSE]

# The metadata has no ready-made age column; `cohort_subset` encodes both the
# group and the age stratum as older_CRC / older_CTRL / younger_CRC /
# younger_CTRL (plus QC, which is dropped with the NA group). The age bin is
# derived from it and the resulting cross-tabulation is asserted below.
df <- data.frame(
  sample = colnames(qn_neg),
  value  = as.numeric(qn_neg[sarc_idx, ]),
  group  = ifelse(grepl("CRC",  sc_meta$cohort_subset), "CRC",
           ifelse(grepl("CTRL", sc_meta$cohort_subset), "CTRL", NA_character_)),
  age    = ifelse(grepl("^older",   sc_meta$cohort_subset), "older",
           ifelse(grepl("^younger", sc_meta$cohort_subset), "younger", NA_character_)),
  stringsAsFactors = FALSE)
df <- df[!is.na(df$group), , drop = FALSE]
df$group <- factor(df$group, levels = c("CTRL", "CRC"))

xt <- table(df$group, df$age)
stopifnot(nrow(df) == 557L,
          sum(df$group == "CTRL") == 246L, sum(df$group == "CRC") == 311L,
          !anyNA(df$value), !anyNA(df$age), !anyDuplicated(df$sample),
          xt["CRC","older"] == 173L, xt["CTRL","older"] == 122L,
          xt["CRC","younger"] == 138L, xt["CTRL","younger"] == 124L)
say("Re-derived from source: n = %d (CTRL %d / CRC %d)",
    nrow(df), sum(df$group == "CTRL"), sum(df$group == "CRC"))
say("Age strata (derived from cohort_subset): CRC %d older / %d younger; CTRL %d older / %d younger",
    xt["CRC","older"], xt["CRC","younger"], xt["CTRL","older"], xt["CTRL","younger"])

# -----------------------------------------------------------------------------
# 2. Integrity: recomputation must reproduce the stored upstream statistics
# -----------------------------------------------------------------------------
v_ctrl <- df$value[df$group == "CTRL"]; v_crc <- df$value[df$group == "CRC"]
wt <- wilcox.test(v_crc, v_ctrl)
up_w <- diff_neg[diff_neg$method == "wilcoxon" & diff_neg$metabolite == "Sarcosine", ]
up_l <- diff_neg[diff_neg$method == "limma"    & diff_neg$metabolite == "Sarcosine", ]
stopifnot(nrow(up_w) == 1L, nrow(up_l) == 1L)
stopifnot(isTRUE(all.equal(wt$p.value, up_w$pvalue, tolerance = 1e-6)),
          isTRUE(all.equal(mean(v_crc) - mean(v_ctrl), up_w$log2FC, tolerance = 1e-6)))
say("Integrity OK: recomputed Wilcoxon p = %.6e matches stored %.6e", wt$p.value, up_w$pvalue)
say("Integrity OK: mean difference %.6f matches stored wilcoxon log2FC %.6f",
    mean(v_crc) - mean(v_ctrl), up_w$log2FC)
say("Upstream age-adjusted limma log2FC = %.4f (padj %.3e) — reported in the manuscript as +0.43",
    up_l$log2FC, up_l$padj)
say("Direction check: mean CRC %.4f > mean CTRL %.4f", mean(v_crc), mean(v_ctrl))

# -----------------------------------------------------------------------------
# 3. ROC  (Lee 2024 Fig 1B analogue)
# -----------------------------------------------------------------------------
# pROC semantics, verified empirically:
#   levels=c(control,case) + direction="<"  ->  AUC = P(case > control)
# Flipping BOTH levels and direction cancels out; flip only one for the complement.
roc_obj <- roc(response = df$group, predictor = df$value,
               levels = c("CTRL", "CRC"), direction = "<", quiet = TRUE)
auc_pt  <- as.numeric(auc(roc_obj))
ci_d    <- ci.auc(roc_obj, method = "delong")

# Cross-check with a brute-force concordance count that uses no pROC convention.
auc_brute <- mean(outer(v_crc, v_ctrl, function(x, y) (x > y) + 0.5 * (x == y)))
stopifnot(isTRUE(all.equal(auc_pt, auc_brute, tolerance = 1e-10)))
say("ROC  AUC = %.4f (95%% CI DeLong %.4f-%.4f); brute-force cross-check %.4f OK",
    auc_pt, ci_d[1], ci_d[3], auc_brute)
say("     Wilcoxon rank-sum p = %.3e", wt$p.value)

say("Age-stratified AUC (sensitivity; upstream model adjusted for age_cohort_bin):")
for (a in sort(unique(df$age))) {
  s <- df[df$age == a, ]
  r <- roc(s$group, s$value, levels = c("CTRL", "CRC"), direction = "<", quiet = TRUE)
  say("   %-8s n = %3d CTRL / %3d CRC   AUC = %.4f (95%% CI %.4f-%.4f)",
      a, sum(s$group == "CTRL"), sum(s$group == "CRC"),
      as.numeric(auc(r)), ci.auc(r)[1], ci.auc(r)[3])
}

# -----------------------------------------------------------------------------
# 4. Median split + stacked bar  (Lee 2024 Fig 1C analogue)
# -----------------------------------------------------------------------------
med <- median(df$value)                                   # POOLED median
df$level <- factor(ifelse(df$value > med, "High", "Low"), levels = c("Low", "High"))
tab <- table(df$group, df$level)
chi <- chisq.test(tab, correct = FALSE)
fis <- fisher.test(tab)
stopifnot(min(chi$expected) > 5)                          # chi-square validity
say("Median split at pooled median %.4f (%d value(s) exactly on the split point)",
    med, sum(df$value == med))
say("  CTRL  Low %3d / High %3d   (%.1f%% / %.1f%%)",
    tab["CTRL","Low"], tab["CTRL","High"],
    100*tab["CTRL","Low"]/sum(tab["CTRL",]), 100*tab["CTRL","High"]/sum(tab["CTRL",]))
say("  CRC   Low %3d / High %3d   (%.1f%% / %.1f%%)",
    tab["CRC","Low"], tab["CRC","High"],
    100*tab["CRC","Low"]/sum(tab["CRC",]), 100*tab["CRC","High"]/sum(tab["CRC",]))
say("  Pearson chi-square X2 = %.2f, df = %d, p = %.3e (min expected %.1f)",
    chi$statistic, chi$parameter, chi$p.value, min(chi$expected))
say("  Fisher exact p = %.3e", fis$p.value)

# -----------------------------------------------------------------------------
# 5. Figures
# -----------------------------------------------------------------------------

p_roc <- ggplot(data.frame(x = 100 * (1 - roc_obj$specificities),
                           y = 100 * roc_obj$sensitivities), aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed", linewidth = .4) +
  geom_step(direction = "hv", linewidth = .7) +
  annotate("text", x = 98, y = 10, hjust = 1, size = ANNOT_SIZE,
           label = sprintf("AUC = %.4f\n95%% CI %.3f-%.3f\np = %.1e", auc_pt, ci_d[1], ci_d[3], wt$p.value)) +
  coord_equal() + scale_x_continuous(limits = c(0, 100)) + scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Fecal Sarcosine",
       x = "100% - Specificity, %", y = "Sensitivity, %") + theme_publication()

bar_df <- as.data.frame(tab); names(bar_df) <- c("grp", "level", "n")
bar_df$pct <- 100 * bar_df$n / ave(bar_df$n, bar_df$grp, FUN = sum)
bar_df$level <- factor(bar_df$level, levels = c("High", "Low"))
p_bar <- ggplot(bar_df, aes(grp, pct, fill = level)) +
  geom_col(width = .6) +
  scale_fill_manual(values = c(Low = COL_LOW, High = COL_HIGH)) +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  annotate("text", x = 1.5, y = 106, size = ANNOT_SIZE,
           label = sprintf("chi-square p = %.1e", chi$p.value)) +
  labs(title = "Fecal Sarcosine", x = NULL, y = "% within group") +
  theme_publication_legend()

save_png(file.path(out_dir, "Fig1b_ROC_CRC_vs_healthy.png"),      p_roc, 6.4, 6.6)
save_png(file.path(out_dir, "Fig1b_stackedbar_median_split.png"), p_bar, 5.8, 6.6)

# -----------------------------------------------------------------------------
# 6. Machine-readable results + panel source data + validation log
# -----------------------------------------------------------------------------
res <- data.frame(
  metric = c("n_CRC", "n_CTRL", "AUC", "AUC_95CI_low_DeLong", "AUC_95CI_high_DeLong",
             "wilcoxon_p", "upstream_wilcoxon_p", "upstream_limma_log2FC", "upstream_limma_padj",
             "pooled_median", "CTRL_low", "CTRL_high", "CRC_low", "CRC_high",
             "chisq_X2", "chisq_p", "fisher_p", "min_expected_count"),
  value = as.character(c(length(v_crc), length(v_ctrl), round(auc_pt, 4), round(ci_d[1], 4), round(ci_d[3], 4),
            signif(wt$p.value, 4), signif(up_w$pvalue, 4), round(up_l$log2FC, 4), signif(up_l$padj, 4),
            round(med, 4), tab["CTRL","Low"], tab["CTRL","High"], tab["CRC","Low"], tab["CRC","High"],
            round(chi$statistic, 3), signif(chi$p.value, 4), signif(fis$p.value, 4),
            round(min(chi$expected), 1))),
  stringsAsFactors = FALSE)
write.csv(res, file.path(out_dir, "Fig1b_ROC_stackedbar_stats.csv"), row.names = FALSE)
write.csv(df[, c("sample", "group", "age", "value", "level")],
          file.path(out_dir, "Fig1b_ROC_stackedbar_source_data.csv"), row.names = FALSE)

writeLines(c("VALIDATION LOG - Fig1b ROC + median-split stacked bar",
             format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), strrep("-", 70), log_lines,
             strrep("-", 70),
             "Assertions enforced: n = 557 (246 CTRL / 311 CRC); exactly one",
             "NEG-mode Sarcosine/HMDB0000271 feature; recomputed Wilcoxon p and",
             "mean difference equal the stored diff_results_NEG.tsv values;",
             "pROC AUC equals a brute-force concordance count; chi-square minimum",
             "expected cell count > 5. Upstream results were not overwritten."),
           file.path(out_dir, "validation_report.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

cat("\nWrote ->", out_dir, "\n")
print(res, row.names = FALSE)
