#!/usr/bin/env Rscript
# =============================================================================
# MULTI-GROUP survival of the 3 cross-cohort sarcosine-degradation-associated
# species -- advisor's request (meeting 2026-07-23) to split abundance into
# more than two ordered groups.
#
# EVERY KM PLOT SHOWS EXACTLY TWO ARMS. Ordered groups are compared PAIRWISE,
# one figure per pair (e.g. Absent vs Low, Absent vs High, Low vs High), so no
# panel ever carries three overlapping curves.
#
# Same patients / same cohort structure / same primary Cox as
# species_survival_km.R; ONLY the grouping changes.
#   OS : PRJNA751792 (338) + PRJNA1023797 (499) = 837, COHORT-STRATIFIED.
#   PFS: PRJNA1023797 only (499). No pooling, no stratification.
#   Cut-points are COHORT-INTERNAL for OS and overall for PFS, as in the
#   median/mean versions.
#
# TWO SCHEMES, WRITTEN TO SEPARATE OUTPUT FOLDERS:
#   [A] results/quartile_survival/
#       Q1..Q4 at the cohort-internal quartiles.
#       ZERO-INFLATION SAFEGUARD (inherited from species_survival_km.R): a
#       species whose median is 0 (>50% zeros) is NOT quantile-cut at all; it
#       is analysed as present/absent, because otherwise the low quantiles
#       would mix "species absent" with "species present at low abundance".
#       Clostridium sp AF36_4 (52% zeros) takes this branch.
#       WARNING for the species that do pass the safeguard: abundance is still
#       right-skewed with 30-39% zeros, so the lowest quartile break sits AT 0
#       and cannot separate tied zeros. Duplicate breaks are collapsed with
#       unique(); the ACTUAL number of groups is logged per species/cohort, and
#       the lowest group still contains all zeros plus some low non-zero
#       samples. The degeneration is reported, not hidden.
#   [B] results/absent_low_high_survival/
#       Zeros form their own group ("Absent"); non-zero samples are split at
#       the cohort-internal median of the NON-ZERO values. Standard handling
#       for a floor-inflated exposure; keeps all 837 patients.
#
# PER-PAIR STATISTICS (shown on each 2-arm figure):
#   - (cohort-stratified for OS) log-rank p for that pair
#   - (cohort-stratified for OS) Cox HR of the higher vs the lower arm
#   - median survival per arm
# ACROSS-GROUP TREND (written to *_trend_*.csv, not on the 2-arm figures):
#   cohort-stratified Cox with the ordered group as a NUMERIC score
#   -> HR per one-group increment; plus the global k-group log-rank.
# Continuous cohort-stratified Cox (HR per 1 SD of log10 abundance) is carried
#   through unchanged for reference -- it does not use the grouping at all and
#   is identical to the median/mean runs.
#
# EXPLORATORY; uncorrected. seed 42.
# Run: Rscript pooled_analysis/R_scripts/species_survival_km_multigroup.R
# =============================================================================
set.seed(42)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(survival); library(broom) })

PAIRCOL <- c(lower = "#0072B2", upper = "#D55E00")   # Okabe-Ito, colourblind-safe

get_script_dir <- function() { a <- commandArgs(FALSE); f <- grep("^--file=", a, value = TRUE)
  if (length(f)) return(dirname(normalizePath(sub("^--file=", "", f)))); normalizePath(".") }
script_dir <- get_script_dir(); pooled_dir <- dirname(script_dir)
res_in <- file.path(pooled_dir, "results")

fmtp  <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-3, sprintf("%.1e", p), sprintf("%.3f", p)))
fslug <- function(x) gsub("[ /]+", "_", x)
save_png <- function(p, out, w, h) { tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = w, height = h, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out, overwrite = TRUE)); unlink(tmp) }

OSD  <- fread(file.path(res_in, "species_survival_OS_dataset_NSCLC.csv"))
PFSD <- fread(file.path(res_in, "species_survival_PFS_dataset_NSCLC.csv"))
SP_DISP <- setdiff(names(OSD), c("RunID", "cohort", "time", "event"))
stopifnot(length(SP_DISP) == 3)

## ---- grouping functions (return integer index, 1 = lowest) ------------------
grp_quartile <- function(v) {
  br <- unique(as.numeric(quantile(v, c(0, .25, .5, .75, 1), type = 7)))
  as.integer(cut(v, breaks = br, include.lowest = TRUE, labels = FALSE))
}
grp_absent_low_high <- function(v) {
  g <- integer(length(v)); nz <- v > 0
  g[!nz] <- 1L
  if (any(nz)) { m <- median(v[nz]); g[nz] <- ifelse(v[nz] > m, 3L, 2L) }
  g
}

## ---- two-arm KM plotter -----------------------------------------------------
km_plot2 <- function(fit, dd, lvl, xlab, title_expr, cap, stat_lab, out) {
  stopifnot(length(lvl) == 2)
  td <- as.data.table(broom::tidy(fit)); td[, arm := sub("grp=", "", strata)]
  t0 <- unique(td[, .(arm)])[, .(time = 0, estimate = 1, arm)]
  cens <- td[n.censor > 0, .(time, estimate, arm)]
  td2 <- rbind(t0, td[, .(time, estimate, arm)]); setorder(td2, arm, time)
  td2[, arm := factor(arm, levels = lvl)]; cens[, arm := factor(arm, levels = lvl)]
  lbl  <- vapply(lvl, function(L) sprintf("%s (n=%d)", L, sum(dd$grp == L)), character(1))
  cols <- setNames(unname(PAIRCOL), lvl)          # lvl[1] = lower group
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
    theme(legend.position = "top", plot.title = element_text(face = "bold", size = 16),
          plot.caption = element_text(size = 8, hjust = 0, color = "grey35"))
  save_png(p, out, 7.2, 5.6)
}

## ---- one scheme -------------------------------------------------------------
run_scheme <- function(scheme, gfun, outdir, labels_of) {
  res <- file.path(res_in, outdir); dir.create(res, showWarnings = FALSE, recursive = TRUE)
  logcon <- file(file.path(res, sprintf("species_survival_%s_log.txt", scheme)), open = "wt")
  say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logcon) }
  PAIRS <- list(); TREND <- list()

  say("SCHEME = %s | every KM figure shows exactly TWO arms (pairwise comparisons)", scheme)
  say("OS: %d patients (%s), deaths=%d | PFS: %d patients (PRJNA1023797 only), events=%d",
      nrow(OSD), paste(OSD[, .N, cohort][, paste0(cohort, "=", N)], collapse = ", "),
      sum(OSD$event), nrow(PFSD), sum(PFSD$event))
  say("Zero-inflation (OS set): %s",
      paste(sapply(SP_DISP, function(s) sprintf("%s %.0f%%", s, 100 * mean(OSD[[s]] == 0))), collapse = " | "))

  one <- function(dd0, sp, endpoint, xlab, stratified) {
    dd <- copy(dd0); dd[, val := get(sp)]

    ## ---- ZERO-INFLATION SAFEGUARD (inherited from species_survival_km.R) ----
    ## Tested on the MEDIAN: "are more than half the samples zero?". A species
    ## that fails this test must NOT be cut into quantiles, because every low
    ## quantile would then contain a mixture of "species absent" (value exactly
    ## 0) and "species present at low abundance" patients -- biologically
    ## different states pooled into one arm. Such a species is analysed as
    ## present/absent, exactly as the original script does.
    ## The absentLowHigh scheme already isolates zeros by construction, so the
    ## safeguard only applies to the quantile scheme.
    if (stratified) { meds <- dd[, .(med = median(val)), by = cohort]; zero_infl <- any(meds$med <= 0) }
    else            { zero_infl <- median(dd$val) <= 0 }
    force_pa <- zero_infl && scheme == "quartile"

    if (force_pa) {
      lvl <- c("Absent", "Present")
      dd[, grp := factor(ifelse(val > 0, "Present", "Absent"), levels = lvl)]
      k <- 2L; k_by <- data.table(cohort = "all", k = 2L)
    } else {
      if (stratified) dd[, gi := gfun(val), by = cohort] else dd[, gi := gfun(val)]
      k_by <- if (stratified) dd[, .(k = uniqueN(gi)), by = cohort] else data.table(cohort = "all", k = uniqueN(dd$gi))
      k <- max(k_by$k); lvl <- labels_of(k)
      dd[, grp := factor(lvl[gi], levels = lvl)]
    }
    degen <- if (force_pa)
      sprintf(" [ZERO-INFLATED: median=0, %.0f%% zeros -> quantiles NOT applicable; analysed as present/absent]",
              100 * mean(dd$val == 0))
      else if (scheme == "quartile" && any(k_by$k < 4))
      sprintf(" [DEGENERATE: %s]", paste(sprintf("%s=%d groups", k_by$cohort, k_by$k), collapse = ", ")) else ""

    ## reference continuous Cox (independent of grouping)
    pseudo <- min(dd$val[dd$val > 0]) / 2; dd[, z := as.numeric(scale(log10(val + pseudo)))]
    cxc <- if (stratified) coxph(Surv(time, event) ~ z + strata(cohort), data = dd) else coxph(Surv(time, event) ~ z, data = dd)
    sc <- summary(cxc)$coefficients; hrc <- exp(sc[1, "coef"]); pc <- sc[1, "Pr(>|z|)"]

    ## across-group trend + global log-rank
    dd[, gnum := as.integer(grp)]
    cxt <- if (stratified) coxph(Surv(time, event) ~ gnum + strata(cohort), data = dd) else coxph(Surv(time, event) ~ gnum, data = dd)
    st <- summary(cxt)$coefficients; hrt <- exp(st[1, "coef"]); cit <- exp(confint(cxt))[1, ]; pt <- st[1, "Pr(>|z|)"]
    lrg <- if (stratified) survdiff(Surv(time, event) ~ grp + strata(cohort), data = dd) else survdiff(Surv(time, event) ~ grp, data = dd)
    plrg <- 1 - pchisq(lrg$chisq, length(lrg$n) - 1)
    say("[%s | %s] %d groups: %s%s | global log-rank p=%s | trend HR/group=%.2f [%.2f-%.2f] p=%s",
        endpoint, sp, k, paste(lvl, table(dd$grp), sep = "=", collapse = ","), degen,
        fmtp(plrg), hrt, cit[1], cit[2], fmtp(pt))
    TREND[[paste(endpoint, sp)]] <<- data.table(grouping = scheme, endpoint = endpoint, species = sp,
      n_groups = k, degenerate = nzchar(degen),
      groups = paste(lvl, table(dd$grp), sep = "=", collapse = "; "),
      global_logrank_p = signif(plrg, 3),
      trend_HR_perGroup = round(hrt, 3), trend_lo = round(cit[1], 3), trend_hi = round(cit[2], 3),
      trend_p = signif(pt, 3), cont_cox_HR_perSD = round(hrc, 3), cont_cox_p = signif(pc, 3))

    ## ---- pairwise 2-arm comparisons ----
    for (i in seq_len(k - 1)) for (j in (i + 1):k) {
      a <- lvl[i]; b <- lvl[j]
      s <- dd[grp %in% c(a, b)]; s[, grp := factor(as.character(grp), levels = c(a, b))]
      lr <- if (stratified) survdiff(Surv(time, event) ~ grp + strata(cohort), data = s) else survdiff(Surv(time, event) ~ grp, data = s)
      plr <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
      cx <- if (stratified) coxph(Surv(time, event) ~ grp + strata(cohort), data = s) else coxph(Surv(time, event) ~ grp, data = s)
      sm <- summary(cx)$coefficients; hr <- exp(sm[1, "coef"]); ci <- exp(confint(cx))[1, ]; pp <- sm[1, "Pr(>|z|)"]
      fit <- survfit(Surv(time, event) ~ grp, data = s)
      mt <- summary(fit)$table
      medv <- if (is.matrix(mt)) setNames(mt[, "median"], sub("grp=", "", rownames(mt))) else setNames(mt["median"], a)
      hrtxt <- sprintf("Cox HR (%s vs %s) = %.2f [%.2f-%.2f], p=%s", b, a, hr, ci[1], ci[2], fmtp(pp))
      lrtxt <- sprintf("%slog-rank p = %s", if (stratified) "stratified " else "", fmtp(plr))
      say("   %s vs %s: n=%d/%d | median %s/%s mo | %s | %s",
          a, b, sum(s$grp == a), sum(s$grp == b),
          round(as.numeric(medv[a]), 1), round(as.numeric(medv[b]), 1), lrtxt, hrtxt)
      PAIRS[[paste(endpoint, sp, a, b)]] <<- data.table(grouping = scheme, endpoint = endpoint, species = sp,
        comparison = sprintf("%s vs %s", b, a), arm_lower = a, arm_upper = b,
        n_lower = sum(s$grp == a), n_upper = sum(s$grp == b),
        median_surv_lower = round(as.numeric(medv[a]), 1), median_surv_upper = round(as.numeric(medv[b]), 1),
        logrank_p = signif(plr, 3), cox_HR_upper_vs_lower = round(hr, 3),
        cox_lo = round(ci[1], 3), cox_hi = round(ci[2], 3), cox_p = signif(pp, 3))
      ttl <- bquote(italic(.(sp)) ~ .(sprintf("- %s: %s vs %s", endpoint, b, a)))
      cap <- sprintf("Advisor-requested multi-group split (%s); two arms per figure. Cut-points %s.\nAcross-group trend and the global %d-group test are in species_survival_%s_trend_NSCLC.csv.\nContinuous-abundance Cox (unchanged reference): HR/SD = %.2f, p = %s. EXPLORATORY; uncorrected.%s%s",
        scheme, if (stratified) "cohort-internal" else "overall", k, scheme, hrc, fmtp(pc),
        if (endpoint == "PFS") "\nPFS event reconstructed = (PFS<OS)|death." else "",
        if (nzchar(degen)) "\nNOTE: zero-inflation collapsed the lower quartile break(s); fewer than 4 groups obtained." else "")
      km_plot2(fit, s, c(a, b), xlab, ttl, cap, paste0(lrtxt, "\n", hrtxt),
               file.path(res, sprintf("KM_%s_%s_%s_vs_%s_NSCLC.png", endpoint, fslug(sp), fslug(a), fslug(b))))
    }
  }

  for (sp in SP_DISP) one(OSD,  sp, "OS",  "Overall survival (months)",          TRUE)
  say("")
  for (sp in SP_DISP) one(PFSD, sp, "PFS", "Progression-free survival (months)", FALSE)

  fwrite(rbindlist(PAIRS), file.path(res, sprintf("species_survival_%s_pairwise_NSCLC.csv", scheme)))
  fwrite(rbindlist(TREND), file.path(res, sprintf("species_survival_%s_trend_NSCLC.csv", scheme)))
  writeLines(capture.output(sessionInfo()), file.path(res, sprintf("sessionInfo_species_survival_%s.txt", scheme)))
  say("DONE. %d pairwise figures in %s", length(PAIRS), res)
  close(logcon)
}

run_scheme("quartile", grp_quartile, "quartile_survival",
           function(k) paste0("Q", seq_len(k)))
run_scheme("absentLowHigh", grp_absent_low_high, "absent_low_high_survival",
           function(k) c("Absent", "Low", "High")[seq_len(k)])
