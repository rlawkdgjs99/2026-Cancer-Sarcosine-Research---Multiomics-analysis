#!/usr/bin/env Rscript
# ============================================================================
# FUNCTION-LEVEL survival: gut-microbial sarcosine-metabolism SCORES vs OS/PFS
# on ICI -- the FUNCTION-level analogue of species_survival_km.R (which tested
# the 3 degradation-associated species). Uses the SAME per-sample scores as
# Fig 3c (degradation = summed K00301-06 + K00314 + K18897; etc.) and the SAME
# survival datasets + cohort-stratified Cox + KM framework.
#   PRE-SPECIFIED primary = DEGRADATION score (the thesis indicator); production
#   and prod:deg log2 ratio reported as secondary (honest, all in summary).
#   OS : cohort-stratified (PRJNA751792 + PRJNA1023797). PFS: PRJNA1023797 only.
#   Cox = HR per 1 SD of (log10) score; KM by cohort-internal median (High/Low)
#   + (stratified) log-rank. EXPLORATORY; seed 42; tempfile->copy PNG.
# ============================================================================
set.seed(42)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(survival); library(broom); library(scales) })
# Run with the R working directory set to the analysis-folder root; paths below
# are relative to it.
pooled <- file.path(".", "pooled_analysis"); main <- "."; res <- file.path(pooled, "results")
logf <- file(file.path(res, "score_survival_km_log.txt"), open = "wt")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logf) }
save_png <- function(p, out, w, h) { tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = w, height = h, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out, overwrite = TRUE)); unlink(tmp) }
fmtp <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-3, sprintf("%.1e", p), sprintf("%.3f", p)))
KMCOL <- c("High" = "#D55E00", "Low" = "#0072B2")

km_plot <- function(fit, dd, lvl, xlab, title_expr, cap, stat_lab, out) {
  td <- as.data.table(broom::tidy(fit)); td[, arm := sub("grp=", "", strata)]
  t0 <- unique(td[, .(arm)])[, .(time = 0, estimate = 1, arm)]
  cens <- td[n.censor > 0, .(time, estimate, arm)]
  td2 <- rbind(t0, td[, .(time, estimate, arm)]); setorder(td2, arm, time)
  td2[, arm := factor(arm, levels = lvl)]; cens[, arm := factor(arm, levels = lvl)]
  lbl <- vapply(lvl, function(L) sprintf("%s (n=%d)", L, sum(dd$grp == L)), character(1))
  max_t <- max(td2$time)
  p <- ggplot(td2, aes(time, estimate, color = arm)) + geom_step(linewidth = 0.9) +
    geom_point(data = cens, aes(time, estimate, color = arm), shape = 3, size = 1.7, alpha = 0.6) +
    annotate("text", x = max_t * 0.015, y = 0.06, hjust = 0, vjust = 0, label = stat_lab, size = 4.8, lineheight = 1.1, color = "grey15") +
    scale_color_manual(values = KMCOL, name = NULL, labels = lbl) +
    coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent) +
    labs(x = xlab, y = "Survival probability", title = title_expr, caption = cap) +
    theme_classic(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(face = "bold", size = 20),
          plot.caption = element_text(size = 8, hjust = 0, color = "grey35"))
  save_png(p, out, 7.2, 5.4)
}

## ---- per-sample scores (Fig 3c definition) joined to survival datasets -------
cmap <- c(PRJNA751792 = "NSCLC_PRJNA751792", PRJNA1023797 = "NSCLC_PRJNA1023797")
scores <- rbindlist(lapply(names(cmap), function(cn)
  fread(file.path(main, cmap[cn], "results", paste0("sarcosine_scores_", cn, ".csv")))[, .(RunID, degradation, production, prod_deg_log2ratio)]),
  use.names = TRUE)
OSD  <- merge(fread(file.path(res, "species_survival_OS_dataset_NSCLC.csv"))[, .(RunID, cohort, time, event)], scores, by = "RunID")
PFSD <- merge(fread(file.path(res, "species_survival_PFS_dataset_NSCLC.csv"))[, .(RunID, cohort, time, event)], scores, by = "RunID")
say("joined: OS n=%d (%s) | PFS n=%d", nrow(OSD), paste(names(table(OSD$cohort)), table(OSD$cohort), sep = "=", collapse = ","), nrow(PFSD))

SUMM <- list()
analyse <- function(dd, metric, transform, endpoint, xlab, stratified, mk_km = FALSE, disp = metric) {
  dd <- copy(dd); dd[, val := get(metric)]; dd <- dd[is.finite(val)]
  if (stratified) dd[, cmed := median(val), by = cohort][, grp := factor(ifelse(val > cmed, "High", "Low"), levels = c("High", "Low"))]
  else { md <- median(dd$val); dd[, grp := factor(ifelse(val > md, "High", "Low"), levels = c("High", "Low"))] }
  if (transform == "log10") { pseudo <- min(dd$val[dd$val > 0]) / 2; dd[, z := as.numeric(scale(log10(val + pseudo)))] } else dd[, z := as.numeric(scale(val)) ]
  cx <- if (stratified) coxph(Surv(time, event) ~ z + strata(cohort), data = dd) else coxph(Surv(time, event) ~ z, data = dd)
  s <- summary(cx)$coefficients; hr <- exp(s[1, "coef"]); ci <- exp(confint(cx))[1, ]; pcox <- s[1, "Pr(>|z|)"]
  fit <- survfit(Surv(time, event) ~ grp, data = dd)
  lr <- if (stratified) survdiff(Surv(time, event) ~ grp + strata(cohort), data = dd) else survdiff(Surv(time, event) ~ grp, data = dd)
  plr <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  mt <- summary(fit)$table; arms <- names(table(dd$grp))
  medv <- if (is.matrix(mt)) setNames(mt[, "median"], sub("grp=", "", rownames(mt))) else setNames(mt["median"], arms[1])
  hrtxt <- sprintf("Cox HR/SD = %.2f [%.2f-%.2f], p=%s", hr, ci[1], ci[2], fmtp(pcox))
  lrtxt <- sprintf("%slog-rank p = %s", if (stratified) "stratified " else "", fmtp(plr))
  say("[%s | %s] %s | %s | %s", endpoint, metric, paste(arms, table(dd$grp), sep = "=", collapse = ","), lrtxt, hrtxt)
  SUMM[[paste(endpoint, metric)]] <<- data.table(endpoint = endpoint, metric = metric, arm = arms, n = as.integer(table(dd$grp)),
      median_surv_months = round(as.numeric(medv[arms]), 1), logrank_p = signif(plr, 3),
      cox_HR_perSD = round(hr, 3), cox_lo = round(ci[1], 3), cox_hi = round(ci[2], 3), cox_p = signif(pcox, 3))
  if (mk_km) {
    stat_lab <- paste0(lrtxt, "\n", hrtxt)
    ttl <- sprintf("Gut sarcosine-%s score and %s", disp, endpoint)
    cap <- sprintf("Pre-specified cohort-internal-median split; cohort-stratified Cox on continuous log10-score = primary test (KM descriptive).%s",
                   if (endpoint == "PFS") "\nPFS event reconstructed = (PFS<OS)|death (PRJNA1023797 only)." else "\nOS pooled PRJNA751792 + PRJNA1023797, cohort-stratified.")
    km_plot(fit, dd, c("High", "Low"), xlab, ttl, cap, stat_lab, file.path(res, sprintf("KM_%s_sarcosine_%s_score_NSCLC.png", endpoint, disp)))
  }
}

## PRIMARY = degradation (KM for OS + PFS); SECONDARY = production, ratio (table) -
analyse(OSD,  "degradation",        "log10", "OS",  "Overall survival (months)",          TRUE,  mk_km = TRUE, disp = "degradation")
analyse(PFSD, "degradation",        "log10", "PFS", "Progression-free survival (months)", FALSE, mk_km = TRUE, disp = "degradation")
analyse(OSD,  "production",         "log10", "OS",  "Overall survival (months)",          TRUE,  disp = "production")
analyse(PFSD, "production",         "log10", "PFS", "Progression-free survival (months)", FALSE, disp = "production")
analyse(OSD,  "prod_deg_log2ratio", "raw",   "OS",  "Overall survival (months)",          TRUE,  disp = "ratio")
analyse(PFSD, "prod_deg_log2ratio", "raw",   "PFS", "Progression-free survival (months)", FALSE, disp = "ratio")

out <- rbindlist(SUMM); fwrite(out, file.path(res, "score_survival_summary_NSCLC.csv"))
say("wrote score_survival_summary_NSCLC.csv + degradation-score KM PNGs (OS, PFS)")
print(out[arm == "High"][, .(endpoint, metric, n_total = NA, cox_HR_perSD, cox_lo, cox_hi, cox_p, logrank_p)])
