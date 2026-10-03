#!/usr/bin/env Rscript
# =============================================================================
# Figure 1c companion — ROC + median-split stacked bar for tissue sarcosine
#   ccRCC tumour vs matched adjacent-normal kidney, Hakimi et al. Cancer Cell 2016
#
# Modelled on Lee et al., Drug Resist Updat 77:101159 (2024), Fig 1B / 1C.
# BUT the design differs from Lee 2024 in one decisive way, and the statistics
# are changed accordingly:
#
#   Lee 2024        : responders vs non-responders = DIFFERENT patients (independent)
#   Hakimi / Fig 1c : tumour vs normal = THE SAME 138 patients (PAIRED)
#
#   -> the 276 observations are 138 independent units, not 276.
#   -> AUC point estimate is unbiased, but its CI must account for clustering:
#      DeLong assumes independence, so a PATIENT-LEVEL CLUSTER BOOTSTRAP is used.
#   -> the 2x2 table has the same patient in both rows, so the correct test is
#      McNEMAR, not chi-square. Chi-square is reported only to document that it
#      is anticonservative here.
#
# Run with the working directory set to RCC_Cancer_Cell_2016/
#   Rscript R_scripts/101_ROC_stackedbar_26.07.27.R
# Existing analyses and figures are read-only; output goes to a NEW folder.
#
# DESIGN
#   Unit of observation : patient (n = 138), each contributing tumour + normal
#   Values              : log2 of the authors' median-normalized Metabolon
#                         abundance (GC/MS). Already normalized upstream.
#   ROC direction       : FIXED A PRIORI (tumour > normal). pROC's default
#                         direction="auto" forces AUC >= 0.5 and would bias it.
#   Split               : POOLED median over all 276 values, per Lee 2024.
#   Sensitivity         : complete-case analysis. In the "Median Normalized"
#                         sheet, below-detection values are imputed with the
#                         per-metabolite minimum. Missingness is DIFFERENTIAL
#                         (4 tumour vs 16 normal), it sits at the very bottom of
#                         the distribution, and therefore acts directly on a
#                         median split and on the lower part of the ROC. The
#                         truly-observed subset (120 pairs) is analysed too.
#
# INTERPRETATION LIMIT (do not overstate):
#   This asks whether sarcosine separates tumour tissue from that patient's own
#   adjacent normal tissue. A within-patient contrast removes between-subject
#   variability, so the AUC is NOT comparable to a between-patient diagnostic
#   AUC such as Lee 2024's. Report it as "discriminates tumour from matched
#   normal tissue", never as diagnostic biomarker performance.
#
# NOTE ON FILE I/O: this project lives under a Korean-named directory. macOS
# stores those names in NFD while R emits NFC, so readxl cannot open the xlsx
# here even via a relative path (base-R I/O and ggsave are unaffected). The
# workbook is therefore copied to tempdir() with base file.copy() and read from
# there. This is a path-encoding workaround only; the bytes are unchanged.
# =============================================================================

suppressPackageStartupMessages({
  library(readxl); library(pROC); library(ggplot2)
})
set.seed(42)
N_BOOT <- 2000

out_dir <- file.path("Outputs", "ROC_stackedbar_26.07.27")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

XLSX_SRC <- "1-s2.0-S1535610815004687-mmc2.xlsx"
XLSX     <- file.path(tempdir(), "hakimi_mmc2.xlsx")
stopifnot(file.exists(XLSX_SRC), file.copy(XLSX_SRC, XLSX, overwrite = TRUE))

COL_LOW <- "#0072B2"; COL_HIGH <- "#E69F00"
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

auc_rank <- function(ctrl, case) {          # Mann-Whitney concordance, tie-aware
  r <- rank(c(ctrl, case)); n0 <- length(ctrl); n1 <- length(case)
  (sum(r[(n0 + 1):(n0 + n1)]) - n1 * (n1 + 1) / 2) / (n0 * n1)
}

# -----------------------------------------------------------------------------
# 1. Extract sarcosine by LABEL matching (not fixed positions)
#    CAUTION: the two sheets label tissue differently --
#      "Median Normalized" uses TUMOR / NORMAL,  "Raw Data" uses T / N.
#    Matching only TUMOR/NORMAL silently selects ZERO columns from Raw Data and
#    reports "0 missing" — a wrong answer that runs cleanly. Both spellings are
#    accepted, and the column count is asserted so it cannot fail silently.
# -----------------------------------------------------------------------------
find_one_row <- function(M, label) {
  h <- which(apply(M, 1, function(r) any(trimws(r) == label, na.rm = TRUE)))
  if (length(h) != 1L) stop(sprintf("Expected exactly 1 '%s' row, found %d", label, length(h)))
  h
}
grab <- function(sheet) {
  M <- suppressMessages(as.matrix(read_excel(XLSX, sheet = sheet, col_names = FALSE,
                                             col_types = "text", .name_repair = "minimal")))
  tt <- find_one_row(M, "TISSUE TYPE"); mt <- find_one_row(M, "MATCHING INDEX")
  sr <- which(apply(M, 1, function(r) any(grepl("sarcosine", r, ignore.case = TRUE))))
  if (length(sr) != 1L) stop("Expected exactly one sarcosine row in ", sheet, "; found ", length(sr))
  tis <- toupper(trimws(M[tt, ])); sc <- which(tis %in% c("TUMOR", "NORMAL", "T", "N"))
  out <- data.frame(matching = trimws(M[mt, sc]),
                    tissue   = ifelse(tis[sc] %in% c("TUMOR", "T"), "TUMOR", "NORMAL"),
                    value    = suppressWarnings(as.numeric(M[sr, sc])),
                    stringsAsFactors = FALSE)
  stopifnot(nrow(out) == 276L,
            sum(out$tissue == "TUMOR") == 138L, sum(out$tissue == "NORMAL") == 138L)
  say("  [%-18s] 276 sample columns (138 TUMOR / 138 NORMAL), sarcosine row %d", sheet, sr)
  out
}
mn  <- grab("Median Normalized")
raw <- grab("Raw Data")

wide <- reshape(mn, idvar = "matching", timevar = "tissue", direction = "wide")
names(wide) <- sub("^value\\.", "", names(wide))
stopifnot(nrow(wide) == 138L, !anyNA(wide$TUMOR), !anyNA(wide$NORMAL))
wide$NORMAL_l2 <- log2(wide$NORMAL); wide$TUMOR_l2 <- log2(wide$TUMOR)

# Integrity: reproduce the PUBLISHED differential-abundance value (Table S3).
fc_pub <- log2(mean(wide$TUMOR) / mean(wide$NORMAL))
stopifnot(abs(fc_pub - 1.4873) < 1e-3)
say("Integrity OK: log2FC(mean T/N) = %.4f reproduces published Table S3 value 1.487", fc_pub)

# Which pairs contain a min-imputed (below-detection) value?
miss_T <- sum(raw$tissue == "TUMOR"  & is.na(raw$value))
miss_N <- sum(raw$tissue == "NORMAL" & is.na(raw$value))
ok_ids <- names(which(tapply(!is.na(raw$value), raw$matching, function(z) length(z) == 2 && all(z))))
stopifnot(miss_T == 4L, miss_N == 16L, length(ok_ids) == 120L)
say("Raw Data below-detection: TUMOR %d / NORMAL %d  ->  %d complete-case pairs",
    miss_T, miss_N, length(ok_ids))
wide$complete_case <- wide$matching %in% ok_ids

# -----------------------------------------------------------------------------
# 2. Analysis, run on all pairs and again on the complete-case subset
# -----------------------------------------------------------------------------
analyse <- function(d, label) {
  vN <- d$NORMAL_l2; vT <- d$TUMOR_l2     # NB: never name these N / T -- in R,
                                          # `T` and `F` abbreviate TRUE / FALSE.
  say("")
  say("--- %s (n = %d pairs) ---", label, nrow(d))
  say("NORMAL median %.3f mean %.3f | TUMOR median %.3f mean %.3f",
      median(vN), mean(vN), median(vT), mean(vT))
  say("within-patient correlation r = %+.3f", cor(vN, vT))
  say("tumour > its own normal: %d/%d (%.1f%%)", sum(vT > vN), nrow(d), 100 * mean(vT > vN))

  ## ROC. pROC semantics verified empirically:
  ##   levels=c(control,case)+direction="<" -> AUC = P(case > control).
  ##   Flipping BOTH levels and direction cancels out; flip only one to invert.
  roc_obj <- roc(response = rep(c("NORMAL", "TUMOR"), each = nrow(d)),
                 predictor = c(vN, vT), levels = c("NORMAL", "TUMOR"),
                 direction = "<", quiet = TRUE)
  auc_pt <- as.numeric(auc(roc_obj))
  brute  <- mean(outer(vT, vN, function(x, y) (x > y) + 0.5 * (x == y)))
  stopifnot(isTRUE(all.equal(auc_pt, brute, tolerance = 1e-10)))
  ci_delong <- ci.auc(roc_obj, method = "delong")          # invalid here; for contrast

  ## Patient-level cluster bootstrap: resample PATIENTS so both tissues of a
  ## patient move together. This is the CI that respects the design.
  bt <- replicate(N_BOOT, { i <- sample.int(nrow(d), replace = TRUE)
                            auc_rank(vN[i], vT[i]) })
  ci_cl <- quantile(bt, c(0.025, 0.975), names = FALSE)

  say("ROC  AUC = %.4f  (brute-force cross-check %.4f OK)", auc_pt, brute)
  say("     95%% CI DeLong (assumes independence -> NOT valid here): %.4f-%.4f",
      ci_delong[1], ci_delong[3])
  say("     95%% CI patient-level cluster bootstrap (VALID, %d reps): %.4f-%.4f",
      N_BOOT, ci_cl[1], ci_cl[2])

  ## Median split
  med <- median(c(vN, vT))
  lvN <- ifelse(vN > med, "High", "Low"); lvT <- ifelse(vT > med, "High", "Low")
  tab <- table(tissue = factor(rep(c("NORMAL", "TUMOR"), each = nrow(d)),
                               levels = c("NORMAL", "TUMOR")),
               level  = factor(c(lvN, lvT), levels = c("Low", "High")))
  ptab <- table(NORMAL = factor(lvN, levels = c("Low", "High")),
                TUMOR  = factor(lvT, levels = c("Low", "High")))
  b_d <- ptab["Low", "High"]; c_d <- ptab["High", "Low"]
  mc  <- mcnemar.test(ptab, correct = TRUE)
  ex  <- binom.test(b_d, b_d + c_d, 0.5)
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))

  say("Median split at pooled median %.4f (%d value(s) exactly on the split point)",
      med, sum(c(vN, vT) == med))
  say("  NORMAL Low %3d / High %3d   |  TUMOR Low %3d / High %3d",
      tab["NORMAL","Low"], tab["NORMAL","High"], tab["TUMOR","Low"], tab["TUMOR","High"])
  say("  patient-level paired table: LowLow %d, Low->High %d, High->Low %d, HighHigh %d",
      ptab["Low","Low"], b_d, c_d, ptab["High","High"])
  say("  McNemar  X2 = %.2f, p = %.3e   <-- CORRECT TEST (%d informative discordant pairs)",
      mc$statistic, mc$p.value, b_d + c_d)
  say("  Exact McNemar (binomial on discordants) p = %.3e", ex$p.value)
  say("  Pearson chi-square p = %.3e  <-- IGNORES PAIRING, anticonservative; not used",
      chi$p.value)

  list(auc = auc_pt, ci_cl = ci_cl, ci_delong = c(ci_delong[1], ci_delong[3]),
       med = med, tab = tab, ptab = ptab, mcnemar_p = mc$p.value, exact_p = ex$p.value,
       chisq_p = chi$p.value, n = nrow(d), roc = roc_obj,
       lvN = lvN, lvT = lvT, label = label)
}

r_all <- analyse(wide, "ALL PAIRS (includes min-imputed values)")
r_cc  <- analyse(wide[wide$complete_case, ], "COMPLETE-CASE (imputed pairs removed)")

say("")
say("Sensitivity: AUC %.4f (all, n=%d) vs %.4f (complete-case, n=%d); McNemar p %.2e vs %.2e",
    r_all$auc, r_all$n, r_cc$auc, r_cc$n, r_all$mcnemar_p, r_cc$mcnemar_p)
say("-> the conclusion does not depend on the min-imputed values.")

# -----------------------------------------------------------------------------
# 3. Figures (all-pairs analysis is the primary display)
# -----------------------------------------------------------------------------

p_roc <- ggplot(data.frame(x = 100 * (1 - r_all$roc$specificities),
                           y = 100 * r_all$roc$sensitivities), aes(x, y)) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed", linewidth = .4) +
  geom_step(direction = "hv", linewidth = .7) +
  annotate("text", x = 98, y = 10, hjust = 1, size = ANNOT_SIZE,
           label = sprintf("AUC = %.4f\n95%% CI %.3f-%.3f",
                           r_all$auc, r_all$ci_cl[1], r_all$ci_cl[2])) +
  coord_equal() + scale_x_continuous(limits = c(0, 100)) + scale_y_continuous(limits = c(0, 100)) +
  labs(title = "Tissue Sarcosine",
       x = "100% - Specificity, %", y = "Sensitivity, %") + theme_publication()

bar_df <- as.data.frame(r_all$tab); names(bar_df) <- c("grp", "level", "n")
bar_df$pct <- 100 * bar_df$n / ave(bar_df$n, bar_df$grp, FUN = sum)
bar_df$level <- factor(bar_df$level, levels = c("High", "Low"))
p_bar <- ggplot(bar_df, aes(grp, pct, fill = level)) +
  geom_col(width = .6) +
  scale_fill_manual(values = c(Low = COL_LOW, High = COL_HIGH)) +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  annotate("text", x = 1.5, y = 106, size = ANNOT_SIZE,
           label = sprintf("McNemar p = %.1e", r_all$mcnemar_p)) +
  labs(title = "Tissue Sarcosine", x = NULL, y = "% within group") +
  theme_publication_legend()

save_png(file.path(out_dir, "Fig1c_ROC_tumour_vs_normal.png"),     p_roc, 6.4, 6.6)
save_png(file.path(out_dir, "Fig1c_stackedbar_median_split.png"), p_bar, 5.8, 6.6)

# -----------------------------------------------------------------------------
# 4. Results, source data, validation log
# -----------------------------------------------------------------------------
res <- data.frame(
  metric = c("n_pairs_all", "AUC_all", "AUC_95CI_low_clusterboot", "AUC_95CI_high_clusterboot",
             "AUC_95CI_low_DeLong_INVALID", "AUC_95CI_high_DeLong_INVALID",
             "within_patient_r", "tumour_gt_own_normal", "pooled_median_all",
             "mcnemar_p_all", "mcnemar_exact_p_all", "chisq_p_all_INVALID",
             "discordant_NlowThigh", "discordant_NhighTlow",
             "n_pairs_complete_case", "AUC_complete_case", "mcnemar_p_complete_case",
             "raw_below_detection_tumour", "raw_below_detection_normal",
             "published_S3_log2FC_reproduced"),
  value = as.character(c(r_all$n, round(r_all$auc, 4), round(r_all$ci_cl[1], 4), round(r_all$ci_cl[2], 4),
            round(r_all$ci_delong[1], 4), round(r_all$ci_delong[2], 4),
            round(cor(wide$NORMAL_l2, wide$TUMOR_l2), 4),
            sprintf("%d/%d", sum(wide$TUMOR_l2 > wide$NORMAL_l2), nrow(wide)),
            round(r_all$med, 4), signif(r_all$mcnemar_p, 4), signif(r_all$exact_p, 4),
            signif(r_all$chisq_p, 4), r_all$ptab["Low","High"], r_all$ptab["High","Low"],
            r_cc$n, round(r_cc$auc, 4), signif(r_cc$mcnemar_p, 4),
            miss_T, miss_N, round(fc_pub, 4))),
  stringsAsFactors = FALSE)
write.csv(res, file.path(out_dir, "Fig1c_ROC_stackedbar_stats.csv"), row.names = FALSE)

src <- data.frame(matching_index = wide$matching,
                  NORMAL_log2 = wide$NORMAL_l2, TUMOR_log2 = wide$TUMOR_l2,
                  NORMAL_level = r_all$lvN, TUMOR_level = r_all$lvT,
                  complete_case = wide$complete_case, stringsAsFactors = FALSE)
write.csv(src, file.path(out_dir, "Fig1c_ROC_stackedbar_source_data.csv"), row.names = FALSE)

writeLines(c("VALIDATION LOG - Fig1c ROC + median-split stacked bar (PAIRED design)",
             format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), strrep("-", 70), log_lines,
             strrep("-", 70),
             "Assertions enforced: 276 sample columns per sheet (138 TUMOR / 138",
             "NORMAL, accepting both TUMOR/NORMAL and T/N spellings); exactly one",
             "sarcosine row per sheet; 138 complete pairs; published Table S3",
             "log2FC = 1.487 reproduced; below-detection counts 4 tumour / 16",
             "normal and 120 complete-case pairs match the project's",
             "REPRODUCIBILITY_sarcosine_tumor_vs_normal.txt; pROC AUC equals a",
             "brute-force concordance count. Upstream results not overwritten.",
             "",
             "Design note: observations are PAIRED (138 patients x 2 tissues), so",
             "the AUC CI uses a patient-level cluster bootstrap and the 2x2 test",
             "is McNemar. The DeLong CI and Pearson chi-square are recorded only",
             "to document that they are not valid for this design."),
           file.path(out_dir, "validation_report.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

cat("\nWrote ->", out_dir, "\n")
print(res, row.names = FALSE)
