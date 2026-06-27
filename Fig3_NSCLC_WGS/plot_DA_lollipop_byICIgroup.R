#!/usr/bin/env Rscript
# ============================================================================
# Top-species differential-abundance LOLLIPOP plots, by ICI-response group.
#   Per-cohort (4) + pooled meta (1). A complement to the volcano figures.
#
# WHY: a volcano shows the whole point cloud; a lollipop RANKS the top species by
# significance and shows each one's effect size (MaAsLin2 coefficient) as a
# stem + dot -- the cleanest way to read "which gut species are most enriched in
# Responders (R) vs Non-responders (NR) to ICI in NSCLC".
#
# INTEGRITY: this script ONLY READS the already-verified differential-abundance
# result CSVs and RE-DRAWS. It does NOT re-run MaAsLin2 / metafor or recompute any
# statistic -- identical philosophy to replot_DA_volcano_byICIgroup.R. The numbers
# plotted are exactly those in:
#   per-cohort : <cohort>/results/DA_species_group_<cohort>.csv  (cols feature,coef,qval,...)
#   pooled     : pooled_analysis/results/pooled_species_meta.csv (cols species,pooled_coef,qval,valid_concordant,...)
#
# ENCODING (per species):
#   stem  : from x=0 to the (pooled) MaAsLin2 coefficient
#   dot   : at the coefficient; COLOR = the ICI group it is enriched in when q<0.05
#           (R / NR; grey if q>=0.05). NR is the model reference => coef>0 higher in
#           R, coef<0 higher in NR (same convention as the volcano).
#   size  : -log10(q)  (bigger dot = more significant)
#   shape : POOLED panel only -- Korean (n=25) sign-concordance of the discovery
#           meta direction (filled = reproduced in Korea; open = not; x = n/a).
#   rows  : the TOP_N species by FDR q (fewer if a cohort has fewer tested).
#
# OUTPUTS (new files):
#   <cohort>/results/DA_lollipop_<cohort>.png                       (4 cohorts)
#   pooled_analysis/results/pooled_species_meta_lollipop.png
#
# Run: Rscript plot_DA_lollipop_byICIgroup.R
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

TOP_N  <- 20                                   # top-N species ranked by FDR q
SIG_Q  <- 0.05                                 # BH-FDR significance threshold
PAL    <- c(NR = "#B2182B", R = "#1B7837")     # group scheme: R=green (favorable), NR=red (unfavorable) (display colour only)
LAB_R  <- "Enriched in Responder (R)"
LAB_NR <- "Enriched in Non-responder (NR)"
LAB_NS <- "n.s. (q >= 0.05)"                   # ASCII only: installed fonts lack >= / arrow glyphs
LOLLI_COLS <- setNames(c(unname(PAL["R"]), unname(PAL["NR"]), "grey75"), c(LAB_R, LAB_NR, LAB_NS))

# Run with the R working directory set to the analysis-folder root; paths below
# are relative to it.
main_dir <- "."

# ragg (ggplot2 4.0 default device) cannot open the non-ASCII (Korean) path on
# macOS; render to an ASCII temp file then copy to the UTF-8 destination.
save_png <- function(p, out_png, width, height) {
  tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = width, height = height, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out_png, overwrite = TRUE)); unlink(tmp)
}

# dt must carry: feature, coef, qval  (+ optional logical 'concord' for the pooled panel)
make_lollipop <- function(dt, title, subtitle, out_png, width = 8.4, height = 6, title_size = 15) {
  d <- copy(as.data.table(dt))
  d <- d[is.finite(coef) & is.finite(qval)]
  setorder(d, qval)
  d <- head(d, TOP_N)                                              # top-N most significant
  n_sig <- sum(d$qval < SIG_Q)
  d[, grp := ifelse(qval >= SIG_Q, LAB_NS, ifelse(coef > 0, LAB_R, LAB_NR))]
  d[, grp := factor(grp, levels = c(LAB_R, LAB_NR, LAB_NS))]
  d[, feature := gsub("_", " ", feature)]
  d[, feature := factor(feature, levels = d[order(-coef)]$feature)] # R-enriched at BOTTOM, NR at top (R-first); display order only

  # ONE legend only (ICI response): the -log10(q) size scale and the Korea-concordance shape
  # scale were dropped per request to declutter. x-axis is REVERSED (scale_x_reverse) so
  # R-enriched (coef>0) point LEFT and NR-enriched point RIGHT -- matching the boxplots'
  # R-left/NR-right convention. scale_x_reverse mirrors the axis only; coef VALUES are unchanged.
  p <- ggplot(d, aes(x = coef, y = feature)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_segment(aes(x = 0, xend = coef, y = feature, yend = feature, color = grp), linewidth = 0.6) +
    geom_point(aes(color = grp), size = 3, shape = 16) +
    scale_color_manual(values = LOLLI_COLS, name = "ICI response", drop = FALSE) +
    scale_x_reverse() +     # flip: R-enriched (coef>0) to the LEFT, NR-enriched to the RIGHT
    labs(x = "MaAsLin2 coefficient   (<-- higher in R    |    higher in NR -->)", y = NULL,
         title = title, subtitle = subtitle) +
    theme_classic(base_size = 16) +
    theme(plot.title = element_text(face = "bold", size = title_size),
          plot.subtitle = element_text(size = 11), legend.position = "right",
          axis.text.y = element_text(face = "italic"))
  save_png(p, out_png, width, height)
  cat(sprintf("  %-44s top %d shown; q<%.2f among them: %d\n", basename(out_png), nrow(d), SIG_Q, n_sig))
}

## ---- per-cohort (4): read the verified DA_species_group_<cohort>.csv ----
cohorts <- c(PRJNA751792 = "NSCLC_PRJNA751792", PRJNA1023797 = "NSCLC_PRJNA1023797",
             PRJEB22863  = "NSCLC_RCC_PRJEB22863", PRJEB26531  = "NSCLC_PRJEB26531")
COVARS  <- list(PRJNA751792 = "age+sex+BMI", PRJNA1023797 = "age+BMI",
                PRJEB22863 = "age+sex+antibiotic", PRJEB26531 = "age")   # mirrors each run_analysis.R
cat("Per-cohort top-species DA lollipops (colored by ICI response group):\n")
for (cn in names(cohorts)) {
  f <- file.path(main_dir, cohorts[cn], "results", paste0("DA_species_group_", cn, ".csv"))
  stopifnot(file.exists(f))
  d <- fread(f)
  stopifnot(all(c("feature", "coef", "qval") %in% names(d)))
  make_lollipop(d,
    title    = sprintf("NSCLC %s: top differential gut species (R vs NR)", cn),
    subtitle = sprintf("MaAsLin2, adj. %s; top %d by FDR q; color = enrichment (q<%.2f BH)",
                       COVARS[[cn]], TOP_N, SIG_Q),
    out_png  = file.path(main_dir, cohorts[cn], "results", paste0("DA_lollipop_", cn, ".png")),
    width = 10.5, height = 6)
}

## ---- pooled meta (1): read the verified pooled_species_meta.csv ----
cat("Pooled meta top-species lollipop:\n")
fp <- file.path(main_dir, "pooled_analysis", "results", "pooled_species_meta.csv")
stopifnot(file.exists(fp))
dp <- fread(fp)
stopifnot(all(c("species", "pooled_coef", "qval") %in% names(dp)))
dp2 <- dp[, .(feature = species, coef = pooled_coef, qval = qval)]
make_lollipop(dp2,
  title    = "Top differential gut species (R vs NR)",
  subtitle = sprintf("Random-effects (DL) inverse-variance meta over 3 discovery cohorts; top %d by meta FDR q", TOP_N),
  out_png  = file.path(main_dir, "pooled_analysis", "results", "pooled_species_meta_lollipop.png"),
  width = 11.5, height = 6, title_size = 22)

writeLines(capture.output(sessionInfo()), file.path(main_dir, "sessionInfo_lollipop.txt"))
cat("DONE.\n")
