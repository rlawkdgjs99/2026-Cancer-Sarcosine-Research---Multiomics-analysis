#!/usr/bin/env Rscript
# =============================================================================
# Clostridium sp AF36_4 -- RAW QUARTILE split with the zero-inflation safeguard
# DELIBERATELY DISABLED, so that the identical procedure used for the other two
# species (Lachnospira eligens, Roseburia faecis) is applied to it as well.
#
# WHY THIS SCRIPT EXISTS
#   The main quartile analysis (species_survival_km_multigroup.R) refuses to
#   quantile-cut this species: its median is 0 in both cohorts (51-53% of
#   patients have exactly 0), so it is analysed as present/absent instead.
#   This script bypasses that rule ON REQUEST, to show what the raw-quartile
#   procedure actually produces for this species.
#
# WHAT TO EXPECT -- READ BEFORE USING THE OUTPUT
#   1. It does NOT reproduce the other two species' 3-group structure.
#      With >50% ties at 0, the 0%, 25% and 50% quantile breaks are ALL 0, so
#      only 3 distinct breaks survive unique() and only TWO groups are formed.
#      (L. eligens and R. faecis have 30-39% zeros, so only the 0% and 25%
#      breaks coincide there, leaving 4 distinct breaks and 3 groups.)
#   2. The lowest group MIXES two biologically different states: patients in
#      whom the species is absent (abundance exactly 0) and patients in whom it
#      is present at low abundance. The exact counts are written to the log and
#      printed on every figure.
#   This is precisely the situation the safeguard was written to prevent. The
#   output is provided for comparison only and should NOT replace the
#   present/absent result in results/quartile_survival/.
#
# Everything else is identical to the main analysis: same patients (datasets
#   written by species_survival_km.R), same cohort stratification for OS,
#   same continuous-abundance Cox reference, two arms per figure.
#
# EXPLORATORY; uncorrected. seed 42.
# Outputs -> pooled_analysis/results/quartile_survival_Clostridium_noSafeguard/
# Run: Rscript pooled_analysis/R_scripts/species_survival_km_quartile_Clostridium_noSafeguard.R
# =============================================================================
set.seed(42)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(survival); library(broom) })

SP <- "Clostridium sp AF36 4"
PAIRCOL <- c("#0072B2", "#D55E00")

get_script_dir <- function() { a <- commandArgs(FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f)))); normalizePath(".") }
script_dir <- get_script_dir(); pooled_dir <- dirname(script_dir)
res_in <- file.path(pooled_dir, "results")
res    <- file.path(res_in, "quartile_survival_Clostridium_noSafeguard")
dir.create(res, showWarnings = FALSE, recursive = TRUE)

logcon <- file(file.path(res, "species_survival_quartile_noSafeguard_log.txt"), open = "wt")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logcon) }
save_png <- function(p, out, w, h) { tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = w, height = h, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out, overwrite = TRUE)); unlink(tmp) }
fmtp  <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-3, sprintf("%.1e", p), sprintf("%.3f", p)))
fslug <- function(x) gsub("[ /]+", "_", x)

OSD  <- fread(file.path(res_in, "species_survival_OS_dataset_NSCLC.csv"))
PFSD <- fread(file.path(res_in, "species_survival_PFS_dataset_NSCLC.csv"))
stopifnot(SP %in% names(OSD))

grp_quartile <- function(v) {
  br <- unique(as.numeric(quantile(v, c(0, .25, .5, .75, 1), type = 7)))
  as.integer(cut(v, breaks = br, include.lowest = TRUE, labels = FALSE))
}

km_plot2 <- function(fit, dd, lvl, xlab, title_expr, cap, stat_lab, out) {
  stopifnot(length(lvl) == 2)
  td <- as.data.table(broom::tidy(fit)); td[, arm := sub("grp=", "", strata)]
  t0 <- unique(td[, .(arm)])[, .(time = 0, estimate = 1, arm)]
  cens <- td[n.censor > 0, .(time, estimate, arm)]
  td2 <- rbind(t0, td[, .(time, estimate, arm)]); setorder(td2, arm, time)
  td2[, arm := factor(arm, levels = lvl)]; cens[, arm := factor(arm, levels = lvl)]
  lbl  <- vapply(lvl, function(L) sprintf("%s (n=%d)", L, sum(dd$grp == L)), character(1))
  cols <- setNames(PAIRCOL, lvl)
  max_t <- max(td2$time)
  p <- ggplot(td2, aes(time, estimate, color = arm)) +
    geom_step(linewidth = 0.9) +
    geom_point(data = cens, aes(time, estimate, color = arm), shape = 3, size = 1.7, alpha = 0.6) +
    annotate("text", x = max_t * 0.015, y = 0.06, hjust = 0, vjust = 0,
             label = stat_lab, size = 4.6, lineheight = 1.15, color = "grey15") +
    scale_color_manual(values = cols, name = NULL, labels = lbl) +
    coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent) +
    labs(x = xlab, y = "Survival probability", title = title_expr, caption = cap) +
    theme_classic(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(face = "bold", size = 15),
          plot.caption = element_text(size = 8, hjust = 0, color = "grey35"))
  save_png(p, out, 7.2, 5.8)
}

PAIRS <- list()
one <- function(dd0, endpoint, xlab, stratified) {
  dd <- copy(dd0); dd[, val := get(SP)]
  ## safeguard intentionally NOT applied here
  if (stratified) dd[, gi := grp_quartile(val), by = cohort] else dd[, gi := grp_quartile(val)]
  k_by <- if (stratified) dd[, .(k = uniqueN(gi)), by = cohort] else data.table(cohort = "all", k = uniqueN(dd$gi))
  k <- max(k_by$k); lvl <- paste0("Q", seq_len(k))
  dd[, grp := factor(lvl[gi], levels = lvl)]

  ## quantify the mixing in the lowest group -- this is the point of the script
  lowest <- lvl[1]
  n_zero_low <- dd[grp == lowest & val == 0, .N]
  n_nz_low   <- dd[grp == lowest & val > 0, .N]
  say("[%s | %s] RAW QUARTILE, safeguard OFF", endpoint, SP)
  say("   zeros overall: %d/%d (%.1f%%)", dd[val == 0, .N], nrow(dd), 100 * dd[val == 0, .N] / nrow(dd))
  say("   quantile breaks per cohort:")
  for (ch in unique(dd$cohort)) {
    v <- dd[cohort == ch, val]
    say("      %-14s %s", ch, paste(sprintf("%.5f", as.numeric(quantile(v, c(0,.25,.5,.75,1), type = 7))), collapse = " | "))
  }
  say("   distinct breaks -> %d group(s) obtained (NOT 4, and NOT the 3 groups the other two species give)",
      k)
  say("   groups: %s", paste(lvl, table(dd$grp), sep = "=", collapse = ", "))
  say("   >>> lowest group %s = %d absent (abundance exactly 0) + %d present-at-low-abundance <<<",
      lowest, n_zero_low, n_nz_low)

  for (i in seq_len(k - 1)) for (j in (i + 1):k) {
    a <- lvl[i]; b <- lvl[j]
    s <- dd[grp %in% c(a, b)]; s[, grp := factor(as.character(grp), levels = c(a, b))]
    lr <- if (stratified) survdiff(Surv(time, event) ~ grp + strata(cohort), data = s) else survdiff(Surv(time, event) ~ grp, data = s)
    plr <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
    cx <- if (stratified) coxph(Surv(time, event) ~ grp + strata(cohort), data = s) else coxph(Surv(time, event) ~ grp, data = s)
    sm <- summary(cx)$coefficients; hr <- exp(sm[1, "coef"]); ci <- exp(confint(cx))[1, ]; pp <- sm[1, "Pr(>|z|)"]
    ## continuous-abundance Cox reference (unchanged by any grouping)
    pseudo <- min(dd$val[dd$val > 0]) / 2; dd[, z := as.numeric(scale(log10(val + pseudo)))]
    cxc <- if (stratified) coxph(Surv(time, event) ~ z + strata(cohort), data = dd) else coxph(Surv(time, event) ~ z, data = dd)
    scc <- summary(cxc)$coefficients; hrc <- exp(scc[1, "coef"]); pc <- scc[1, "Pr(>|z|)"]
    fit <- survfit(Surv(time, event) ~ grp, data = s)
    mt <- summary(fit)$table
    medv <- if (is.matrix(mt)) setNames(mt[, "median"], sub("grp=", "", rownames(mt))) else setNames(mt["median"], a)
    hrtxt <- sprintf("Cox HR (%s vs %s) = %.2f [%.2f-%.2f], p=%s", b, a, hr, ci[1], ci[2], fmtp(pp))
    lrtxt <- sprintf("%slog-rank p = %s", if (stratified) "stratified " else "", fmtp(plr))
    say("   %s vs %s: n=%d/%d | median %s/%s mo | %s | %s", a, b,
        sum(s$grp == a), sum(s$grp == b),
        round(as.numeric(medv[a]), 1), round(as.numeric(medv[b]), 1), lrtxt, hrtxt)
    PAIRS[[paste(endpoint, a, b)]] <<- data.table(
      endpoint = endpoint, species = SP, scheme = "raw quartile, safeguard OFF",
      n_groups_obtained = k, comparison = sprintf("%s vs %s", b, a),
      arm_lower = a, arm_upper = b, n_lower = sum(s$grp == a), n_upper = sum(s$grp == b),
      lowest_group_absent = n_zero_low, lowest_group_present_low = n_nz_low,
      median_surv_lower = round(as.numeric(medv[a]), 1), median_surv_upper = round(as.numeric(medv[b]), 1),
      logrank_p = signif(plr, 3), cox_HR_upper_vs_lower = round(hr, 3),
      cox_lo = round(ci[1], 3), cox_hi = round(ci[2], 3), cox_p = signif(pp, 3),
      cont_cox_HR_perSD = round(hrc, 3), cont_cox_p = signif(pc, 3))

    ttl <- bquote(italic(.(SP)) ~ .(sprintf("- %s: %s vs %s (raw quartile)", endpoint, b, a)))
    ## caption states exactly what happened -- no generic wording
    cap <- sprintf(paste0(
      "COMPARISON ONLY -- zero-inflation safeguard DISABLED so the raw-quartile procedure used for the other\n",
      "two species is applied to this species as well. %.0f%% of patients have abundance exactly 0, so the 0%%, 25%%\n",
      "and 50%% quantile breaks all fall at 0 and only %d group(s) are obtained (the other two species give 3).\n",
      "The lower arm %s pools %d ABSENT patients (abundance = 0) with %d present-at-low-abundance patients.\n",
      "The present/absent result in results/quartile_survival/ is the one to use. Continuous-abundance Cox\n",
      "(unaffected by any grouping): HR/SD = %.2f, p = %s. EXPLORATORY; uncorrected.%s"),
      100 * dd[val == 0, .N] / nrow(dd), k, a, n_zero_low, n_nz_low, hrc, fmtp(pc),
      if (endpoint == "PFS") "\nPFS event reconstructed = (PFS<OS)|death." else "")
    km_plot2(fit, s, c(a, b), xlab, ttl, cap, paste0(lrtxt, "\n", hrtxt),
             file.path(res, sprintf("KM_%s_%s_%s_vs_%s_rawQuartile_NSCLC.png", endpoint, fslug(SP), a, b)))
  }
  say("")
}

say("Clostridium sp AF36_4 -- RAW QUARTILE, zero-inflation safeguard DISABLED (comparison only)")
say("OS: %d patients (%s), deaths=%d | PFS: %d patients (PRJNA1023797 only), events=%d",
    nrow(OSD), paste(OSD[, .N, cohort][, paste0(cohort, "=", N)], collapse = ", "),
    sum(OSD$event), nrow(PFSD), sum(PFSD$event))
one(OSD,  "OS",  "Overall survival (months)",          TRUE)
one(PFSD, "PFS", "Progression-free survival (months)", FALSE)

fwrite(rbindlist(PAIRS), file.path(res, "species_survival_quartile_noSafeguard_Clostridium_NSCLC.csv"))
writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_species_survival_quartile_noSafeguard.txt"))
say("DONE. outputs in %s", res)
close(logcon)
