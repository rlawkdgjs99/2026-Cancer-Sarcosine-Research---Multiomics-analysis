#!/usr/bin/env Rscript
# ============================================================================
# SUPPLEMENTARY: Hungatella hathewayi <-> gut sarcosine PRODUCTION-capacity
# association across the NSCLC ICI cohorts -- the robust, defensible positive
# finding (all 4 cohorts positive Spearman rho; discovery random-effects meta
# q < 1e-4; I^2 = 0). Forest of per-cohort rho (Fisher-z 95% CI) + discovery
# meta + Korea validation. READ-ONLY on the verified correlation outputs; no
# statistic recomputed (per-cohort rho/n from sarcosine_species_assoc_percohort.csv;
# pooled rho/CI/q from sarcosine_species_assoc_meta.csv).
# Run: Rscript suppl_hungatella_production.R
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
# Run with the R working directory set to the analysis-folder root; paths below are relative to it.
main_dir <- "."; res <- file.path(main_dir, "pooled_analysis", "results")
save_png <- function(p, out, w, h) { tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = w, height = h, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out, overwrite = TRUE)); unlink(tmp) }

SP <- "Hungatella_hathewayi"
pc <- fread(file.path(res, "sarcosine_species_assoc_percohort.csv"))[species == SP & score == "production", .(cohort, rho, n)]
mt <- fread(file.path(res, "sarcosine_species_assoc_meta.csv"))[species == SP & score == "production"]
stopifnot(nrow(pc) >= 1, nrow(mt) == 1)

## per-cohort Fisher-z 95% CI
pc[, z := atanh(rho)][, se := 1 / sqrt(n - 3)][, lo := tanh(z - 1.96 * se)][, hi := tanh(z + 1.96 * se)]
CLAB <- c(PRJNA751792 = "PRJNA751792 (FR, n=%d)", PRJNA1023797 = "PRJNA1023797 (FR, n=%d)",
          PRJEB22863 = "PRJEB22863 (FR, n=%d)", PRJEB26531 = "PRJEB26531 (KR validation, n=%d)")
ROLE <- c(PRJNA751792 = "discovery", PRJNA1023797 = "discovery", PRJEB22863 = "discovery", PRJEB26531 = "validation")
pc[, label := sprintf(CLAB[cohort], n)][, role := ROLE[cohort]]
pc[, p := 2 * (1 - pnorm(abs(z / se)))]                       # per-cohort Fisher-z p (consistent with the CI shown)
fmtp <- function(pv) ifelse(pv < 1e-3, sprintf("p=%.1e", pv), sprintf("p=%.3f", pv))

fd <- rbind(
  pc[, .(label, rho, lo, hi, role, p)],
  data.table(label = "Pooled (3 discovery, RE)", rho = mt$pooled_rho, lo = mt$ci_lb, hi = mt$ci_ub, role = "pooled", p = mt$pval)
)
lev <- c("Pooled (3 discovery, RE)", pc[cohort == "PRJEB26531"]$label,
         pc[cohort == "PRJEB22863"]$label, pc[cohort == "PRJNA1023797"]$label, pc[cohort == "PRJNA751792"]$label)
fd[, label := factor(label, levels = lev)]
fd[, stat_label := sprintf("%.2f [%.2f, %.2f], %s", rho, lo, hi, fmtp(p))]
fd[role == "pooled", stat_label := sprintf("%.2f [%.2f, %.2f], %s, q=%.1e", rho, lo, hi, fmtp(p), mt$qval)]

COLR  <- c(discovery = "#4C4C4C", validation = "#0072B2", pooled = "#D55E00")
xstat <- max(fd$hi) + 0.03
p <- ggplot(fd, aes(rho, label, color = role)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.22, linewidth = 0.6) +
  geom_point(aes(size = role, shape = role)) +
  geom_text(aes(x = xstat, label = stat_label), hjust = 0, size = 3, color = "grey15") +
  scale_color_manual(values = COLR, guide = "none") +
  scale_shape_manual(values = c(discovery = 16, validation = 17, pooled = 18), guide = "none") +
  scale_size_manual(values = c(discovery = 3, validation = 3.4, pooled = 4.8), guide = "none") +
  scale_x_continuous(breaks = c(0, 0.5, 1.0), limits = c(min(fd$lo) - 0.03, xstat + 0.84)) +
  labs(x = "Spearman rho   (species abundance vs gut sarcosine PRODUCTION score)", y = NULL,
       title = expression(paste(italic("Hungatella hathewayi"), ": gut sarcosine-production association (NSCLC, R vs NR)")),
       subtitle = "Per-cohort Spearman -> random-effects (DL) meta over 3 discovery cohorts; Korea shown as external validation.",
       caption = "All 4 cohorts positive; discovery meta q<1e-4, I2=0 (homogeneous).\nFunctional (production) association; the species' R-vs-NR differential abundance was a non-significant NR trend (see main text).") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(size = 12), axis.text.y = element_text(size = 10),
        plot.subtitle = element_text(size = 9, color = "grey30"),
        plot.caption = element_text(size = 7.6, color = "grey40", hjust = 0))
save_png(p, file.path(res, "suppl_Hungatella_production_forest.png"), 10.8, 4)
fwrite(fd, file.path(res, "suppl_Hungatella_production_forest.csv"))
writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_suppl_hungatella.txt"))
cat("wrote suppl_Hungatella_production_forest.png/.csv\n"); print(fd[, .(label, rho, lo, hi, role)])
