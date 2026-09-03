#!/usr/bin/env Rscript
# ============================================================================
# POOLED (combined-sample) SARCOSINE functional figure - the CONVENTIONAL view,
# done COHORT-AWARE. Complements pooled_sarcosine.R (the meta-analytic KO +
# score synthesis, which remains the PRIMARY pooled inference). Here all WGS
# samples are pooled into ONE box per pathway score (R vs NR), as is
# conventional - but because the 4 cohorts differ by country/platform
# (3 France/Ion-Torrent discovery + 1 Korea/Illumina validation), that confound
# is made VISIBLE (jitter colored by cohort) and ADJUSTED FOR in the test. A
# naive R-vs-NR pooled test would confound ICI response with platform; we never
# report that as the result.
#
#   Pathway scores (per sample) = SUM of relative abundance (%) of the KEGG KOs
#     in each sarcosine group: degradation (now incl. K18897 sarcosine->betaine),
#     production, bidirectional - exactly as defined in each cohort's run_sarcosine.R.
#   Read from the verified per-cohort sarcosine_scores_<cohort>.csv (NOT
#     recomputed); medians sanity-checked against sarcosine_score_stats_<cohort>.csv.
#   Test per score (annotated):
#     - naive Wilcoxon (all samples pooled)                  [descriptive]
#     - cohort-ADJUSTED rank ANCOVA: lm(rank(score) ~ cohort + response)
#       (scores are zero-inflated / right-skewed, so a rank-based test is used
#        rather than a raw-value lm).
#   PRIMARY pooled inference stays the discovery Stouffer / SMD meta in
#     pooled_sarcosine.R (degradation = all sarcosine-consuming KOs, incl. K18897);
#     its per-score combined-p is carried into pooled_sarcosine_combined_stats.csv
#     for cross-reference.
#
# INTEGRITY: no statistic in pooled_sarcosine.R is modified; this only reads the
# saved per-sample scores and re-draws. ragg cannot write the non-ASCII (Korean)
# path on macOS -> render to ASCII temp then base R file.copy().
#
# Run: Rscript pooled_analysis/R_scripts/pooled_sarcosine_combined.R
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

COH  <- c(PRJNA751792 = "NSCLC_PRJNA751792", PRJNA1023797 = "NSCLC_PRJNA1023797",
          PRJEB22863  = "NSCLC_RCC_PRJEB22863", PRJEB26531  = "NSCLC_PRJEB26531")
CLAB <- c(PRJNA751792 = "PRJNA751792 (FR, IonTorrent)", PRJNA1023797 = "PRJNA1023797 (FR, IonTorrent)",
          PRJEB22863  = "PRJEB22863 (FR, IonTorrent)",  PRJEB26531  = "PRJEB26531 (KR, Illumina)")
SCORES <- c("degradation", "production", "prod_deg_log2ratio")
SCLAB  <- c(degradation = "Degradation\n(sarcosine consumed)", production = "Production\n(-> sarcosine)",
            prod_deg_log2ratio = "log2(Prod / Deg)\n(balance; both > 0)")
RESP_PAL <- c(NR = "#B2182B", R = "#1B7837")   # group scheme: R=green (favorable), NR=red (unfavorable) (display colour only)
COH_PAL  <- c(PRJNA751792 = "#E69F00", PRJNA1023797 = "#009E73", PRJEB22863 = "#56B4E9", PRJEB26531 = "#CC79A7")

save_png <- function(p, out_png, width, height) {
  tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = width, height = height, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out_png, overwrite = TRUE)); unlink(tmp)
}
# Run with the R working directory set to the analysis-folder root; paths below
# are relative to it.
pooled_dir <- file.path(".", "pooled_analysis"); main_dir <- "."
res <- file.path(pooled_dir, "results")
logcon <- file(file.path(res, "pooled_sarcosine_combined_log.txt"), open = "wt")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logcon) }
say("==== POOLED sarcosine COMBINED (cohort-aware) | %s ====", format(Sys.time(), tz = "UTC", usetz = TRUE))

## ---- load saved per-cohort scores; sanity-check medians vs saved stats ------
rows <- list()
for (cn in names(COH)) {
  cdir <- file.path(main_dir, COH[cn])
  sc <- fread(file.path(cdir, "results", paste0("sarcosine_scores_", cn, ".csv")))
  st <- fread(file.path(cdir, "results", paste0("sarcosine_score_stats_", cn, ".csv")))
  for (s in SCORES) {
    mr  <- median(sc[group == "R"][[s]], na.rm = TRUE); mnr <- median(sc[group == "NR"][[s]], na.rm = TRUE)
    stopifnot(isTRUE(all.equal(mr,  st[score == s]$median_R)),
              isTRUE(all.equal(mnr, st[score == s]$median_NR)))
  }
  sc[, cohort := cn]; rows[[cn]] <- sc
  say("%s: score medians MATCH saved (R=%d, NR=%d)", cn, sum(sc$group == "R"), sum(sc$group == "NR"))
}
A <- rbindlist(rows, use.names = TRUE)
A[, group  := factor(group, levels = c("NR", "R"))]
A[, cohort := factor(cohort, levels = names(COH))]
say("Pooled: %d WGS samples (R=%d, NR=%d) across %d cohorts",
    nrow(A), sum(A$group == "R"), sum(A$group == "NR"), length(COH))

## ---- primary inference (existing Stouffer meta) for cross-reference ---------
stf <- fread(file.path(res, "pooled_sarcosine_score_meta.csv"))   # score, combined_p, combined_q, pooled_dir

## ---- tests: naive Wilcoxon + cohort-adjusted rank ANCOVA --------------------
stat <- rbindlist(lapply(SCORES, function(s) {
  v <- A[[s]]; fin <- is.finite(v)
  vv <- v[fin]; gg <- A$group[fin]; cc <- A$cohort[fin]      # drop NA (log-ratio undefined where a score is 0)
  p_naive <- wilcox.test(vv ~ gg)$p.value
  cf <- summary(lm(rank(vv) ~ cc + gg))$coefficients          # rank ANCOVA, response adj. for cohort
  data.table(score = s, p_naive_wilcox = p_naive,
             p_adj_rankancova = cf["ggR", "Pr(>|t|)"],
             dir = ifelse(cf["ggR", "Estimate"] > 0, "higher in R", "higher in NR"),
             stouffer_disc_p = stf[score == s]$combined_p, stouffer_disc_q = stf[score == s]$combined_q,
             stouffer_disc_dir = stf[score == s]$pooled_dir,
             disc_dirs = stf[score == s]$percoh_dirs, disc_concord = stf[score == s]$concordant_dir)
}))
fwrite(stat, file.path(res, "pooled_sarcosine_combined_stats.csv"))
say("scores: %s", paste(stat$score, sprintf("adjP=%.3g(%s) naiveP=%.3g StoufferDiscP=%.3g",
    stat$p_adj_rankancova, stat$dir, stat$p_naive_wilcox, stat$stouffer_disc_p), collapse = " | "))

## ---- combined box + jitter figure ------------------------------------------
AL <- melt(A, id.vars = c("RunID", "cohort", "group"), measure.vars = SCORES,
           variable.name = "score", value.name = "value")
AL[, score := factor(score, levels = SCORES)]

## ---- y-axis ZOOM (scores have tiny units; a few extreme outliers compress every box) -------
## Draw each box from the FULL-data five-number summary (boxplot.stats -> Tukey whiskers, which
## exclude outliers) via geom_boxplot(stat="identity"), and show jitter only within the panel's
## whisker range. INTEGRITY: medians / quartiles / whiskers and ALL p-values use the FULL data
## (n=849); ONLY the out-of-whisker POINTS are not drawn (count logged + disclosed in caption).
## NR stays the model reference (factor levels unchanged); scale_x_discrete sets DISPLAY order only.
bx <- AL[is.finite(value), as.list(setNames(boxplot.stats(value)$stats,
          c("ymin","lower","middle","upper","ymax"))), by = .(score, group)]
rng <- bx[, .(lo = min(ymin), hi = max(ymax)), by = score]              # per-panel view range (both groups)
jit <- merge(AL[is.finite(value)], rng, by = "score")
jit <- jit[value >= lo & value <= hi]                                  # points shown (within whiskers)
n_clip <- merge(AL[is.finite(value), .(n = .N), by = score], jit[, .(shown = .N), by = score], by = "score")
n_clip[, clipped := n - shown]
say("y-zoom clipped points (beyond whiskers; NOT dropped from any statistic): %s",
    paste(n_clip$score, n_clip$clipped, sep = "=", collapse = ", "))
ann <- merge(copy(stat), rng, by = "score")
ann[, lab := sprintf("p = %.2g", p_adj_rankancova)]
ann[, y := lo + 0.96 * (hi - lo)]
for (D in list(bx, jit, ann, rng)) D[, score := factor(score, levels = SCORES)]

p <- ggplot(bx, aes(group)) +
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax, fill = group),
               stat = "identity", width = 0.6, alpha = 0.22, colour = "grey30") +
  geom_jitter(data = jit, aes(group, value), width = 0.18, height = 0, size = 0.7, alpha = 0.35, color = "grey25") +
  geom_text(data = ann, aes(x = 1.5, y = y, label = lab), inherit.aes = FALSE,
            size = 5.5, vjust = 1, fontface = "bold") +
  facet_wrap(~ score, nrow = 1, scales = "free_y", labeller = labeller(score = SCLAB)) +
  scale_fill_manual(values = RESP_PAL, guide = "none") +
  scale_x_discrete(limits = c("R", "NR")) +                            # DISPLAY order only: R (left) then NR (right)
  labs(x = "ICI response group", y = NULL,
       title = "Gut sarcosine metabolism (R vs NR)",
       caption = "y-axis zoomed to each score's whiskers for legibility; boxes, medians, and all p-values use the full n=849 (outliers beyond whiskers not drawn).") +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 22),
        plot.caption = element_text(size = 10, hjust = 0, color = "grey35"),
        legend.position = "none", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 16))
save_png(p, file.path(res, "pooled_sarcosine_combined.png"), width = 10.5, height = 5.3)
say("wrote pooled_sarcosine_combined.png")

writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_pooled_sarcosine_combined.txt"))
say("DONE. outputs in %s", res)
close(logcon)
