#!/usr/bin/env Rscript
# =============================================================================
# Kaplan-Meier survival by ABUNDANCE (High vs Low) of the 3 cross-cohort
# sarcosine-DEGRADATION-associated species (the 3-set discovery Venn intersection):
#   Lachnospira eligens, Roseburia faecis, Clostridium sp AF36_4.
#
# ENDPOINTS / COHORTS (selected by data availability; verified counts):
#   - OS : pooled PRJNA751792 (OS + OS_evt, 338) + PRJNA1023797 (OS + Death, 499)
#          = 837 patients, 447 deaths. COHORT-STRATIFIED (platforms differ:
#          Ion Torrent Proton vs S5 XL) -- NOT naive pooling.
#   - PFS: PRJNA1023797 ONLY (499) -- the single cohort with continuous PFS.
#          PFS event RECONSTRUCTED (documented assumption; metadata has PFS time +
#          Death but no progression flag): pfs_event = (PFS<OS) | (Death==1),
#          pfs_time = pmin(PFS, OS) [PFS>OS anomalies censored at OS].
#   PRJEB22863 (categorical 3-month PFS only) and Korea PRJEB26531 (no survival)
#   are excluded.
#
# SPLIT: High vs Low at the COHORT-INTERNAL median; present/absent if the metric
#   is zero-inflated (any cohort median <= 0). Pre-specified (no cutpoint search).
# TEST: primary = COHORT-STRATIFIED Cox on the CONTINUOUS metric, HR per 1 SD of
#   log10(abundance + pseudo); KM + (stratified for OS) log-rank shown alongside.
# HONEST: EXPLORATORY (3 species x 2 endpoints = 6 KM); not corrected here -> state.
#   Separate NSCLC analysis from the CRC work. seed 42; ragg-safe PNG write.
# Run: Rscript pooled_analysis/R_scripts/species_survival_km.R
# =============================================================================
set.seed(42)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(survival); library(broom) })

NSCLC <- "Carcinoma, Non-Small-Cell Lung"
KMCOL <- c("High" = "#D55E00", "Low" = "#0072B2", "Present" = "#D55E00", "Absent" = "#0072B2")

# Run with the R working directory set to the analysis-folder root; paths below
# are relative to it.
pooled_dir <- file.path(".", "pooled_analysis"); main_dir <- "."
res <- file.path(pooled_dir, "results")
logcon <- file(file.path(res, "species_survival_km_log.txt"), open = "wt")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logcon) }
save_png <- function(p, out, w, h) { tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = w, height = h, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out, overwrite = TRUE)); unlink(tmp) }
short <- function(x) sub(".*\\|s__", "", x)
fmtp <- function(p) ifelse(is.na(p), "NA", ifelse(p < 1e-3, sprintf("%.1e", p), sprintf("%.3f", p)))

read_meta <- function(cdir) { f <- list.files(cdir, "^selected_project_.*\\.txt$", full.names = TRUE)[1]
  ln <- sub("\r$", "", readLines(f, warn = FALSE)); hdr <- strsplit(ln[2], "\t", fixed = TRUE)[[1]]
  dat <- ln[-(1:2)]; dat <- dat[nzchar(dat)]
  mat <- do.call(rbind, lapply(strsplit(dat, "\t", fixed = TRUE), function(x) x[seq_len(length(hdr))]))
  d <- as.data.table(mat); setnames(d, hdr); d }
getval <- function(s, k) { p <- strsplit(s, ";", fixed = TRUE)[[1]]; p <- p[startsWith(p, paste0(k, ":"))]; if (length(p)) sub(paste0("^", k, ":"), "", p[1]) else NA_character_ }
gv <- function(x, k) vapply(x, getval, character(1), k = k, USE.NAMES = FALSE)
read_species <- function(cdir) { f <- list.files(cdir, "^Bacteria_.*\\.txt$", full.names = TRUE)[1]
  b <- fread(f, sep = "\t", header = TRUE, quote = "", showProgress = FALSE); setnames(b, c("Taxa", "RunID", "Abundance"))
  sp <- b[grepl("|s__", Taxa, fixed = TRUE) & !grepl("|t__", Taxa, fixed = TRUE)]
  w <- dcast(sp, RunID ~ Taxa, value.var = "Abundance", fill = 0, fun.aggregate = sum)
  m <- as.matrix(w[, -1L, with = FALSE]); rownames(m) <- w$RunID; m }

## ---- targets: the 3 intersection species -----------------------------------
inter <- fread(file.path(res, "cross_cohort_deg_assoc_intersection_discovery_NSCLC.csv")); setorder(inter, -mean_disc_rho)
TARGET <- inter$Taxa; SP_DISP <- gsub("_", " ", short(TARGET)); stopifnot(length(TARGET) == 3)
say("3 species: %s", paste(SP_DISP, collapse = ", "))

## helper: build per-cohort sample table (time, event, 3 species abundance) ----
build_cohort <- function(dirn, time_key, evt_key) {
  cdir <- file.path(main_dir, dirn); m <- read_meta(cdir); m <- m[`Assay type` == "WGS" & `Phenotype name` == NSCLC]
  desc <- m$`Sample description`; rid <- m[["Run ID"]]
  spm <- read_species(cdir); stopifnot(all(rid %in% rownames(spm)))
  ab <- sapply(TARGET, function(tx) if (tx %in% colnames(spm)) spm[rid, tx] else rep(0, length(rid)))
  dt <- data.table(RunID = rid, time = suppressWarnings(as.numeric(gv(desc, time_key))),
                   OS = suppressWarnings(as.numeric(gv(desc, "OS"))),
                   PFS = suppressWarnings(as.numeric(gv(desc, "PFS"))),
                   Death = suppressWarnings(as.integer(gv(desc, "Death"))),
                   evt_raw = suppressWarnings(as.integer(gv(desc, evt_key))))
  for (k in seq_along(TARGET)) dt[[SP_DISP[k]]] <- ab[, k]
  dt
}

## ---- KM plotter --------------------------------------------------------------
km_plot <- function(fit, dd, lvl, xlab, title_expr, cap, stat_lab, out) {
  td <- as.data.table(broom::tidy(fit)); td[, arm := sub("grp=", "", strata)]
  t0 <- unique(td[, .(arm)])[, .(time = 0, estimate = 1, arm)]
  cens <- td[n.censor > 0, .(time, estimate, arm)]
  td2 <- rbind(t0, td[, .(time, estimate, arm)]); setorder(td2, arm, time)
  td2[, arm := factor(arm, levels = lvl)]; cens[, arm := factor(arm, levels = lvl)]
  lbl <- vapply(lvl, function(L) sprintf("%s (n=%d)", L, sum(dd$grp == L)), character(1))
  max_t <- max(td2$time)
  p <- ggplot(td2, aes(time, estimate, color = arm)) +
    geom_step(linewidth = 0.9) +
    geom_point(data = cens, aes(time, estimate, color = arm), shape = 3, size = 1.7, alpha = 0.6) +
    # key statistics in the empty lower-left of the plot (large, readable) -- not in the small caption
    annotate("text", x = max_t * 0.015, y = 0.06, hjust = 0, vjust = 0,
             label = stat_lab, size = 4.8, lineheight = 1.1, color = "grey15") +
    scale_color_manual(values = KMCOL, name = NULL, labels = lbl) +
    coord_cartesian(ylim = c(0, 1)) + scale_y_continuous(labels = scales::percent) +
    labs(x = xlab, y = "Survival probability", title = title_expr, caption = cap) +
    theme_classic(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(face = "bold", size = 20),
          plot.caption = element_text(size = 8, hjust = 0, color = "grey35"))
  save_png(p, out, 7.2, 5.4)
}

## helper: split, stratified Cox, KM + log-rank for one species/endpoint --------
SUMM <- list()
analyse <- function(dd, sp, endpoint, xlab, stratified) {
  dd <- copy(dd); dd[, val := get(sp)]
  if (stratified) { meds <- dd[, .(med = median(val)), by = cohort]; zero_infl <- any(meds$med <= 0) }
  else            { zero_infl <- median(dd$val) <= 0 }
  if (zero_infl) { dd[, grp := factor(ifelse(val > 0, "Present", "Absent"), levels = c("Present", "Absent"))]; split_type <- "present/absent (zero-inflated)"; lvl <- c("Present", "Absent") }
  else if (stratified) { dd[, cmed := median(val), by = cohort][, grp := factor(ifelse(val > cmed, "High", "Low"), levels = c("High", "Low"))]; split_type <- "cohort-internal median"; lvl <- c("High", "Low") }
  else { md <- median(dd$val); dd[, grp := factor(ifelse(val > md, "High", "Low"), levels = c("High", "Low"))]; split_type <- "median"; lvl <- c("High", "Low") }
  dd[, grp := droplevels(grp)]
  pseudo <- min(dd$val[dd$val > 0]) / 2; dd[, z := as.numeric(scale(log10(val + pseudo)))]
  cx <- if (stratified) coxph(Surv(time, event) ~ z + strata(cohort), data = dd) else coxph(Surv(time, event) ~ z, data = dd)
  s <- summary(cx)$coefficients; hr <- exp(s[1, "coef"]); ci <- exp(confint(cx))[1, ]; pcox <- s[1, "Pr(>|z|)"]
  fit <- survfit(Surv(time, event) ~ grp, data = dd)
  lr <- if (stratified) survdiff(Surv(time, event) ~ grp + strata(cohort), data = dd) else survdiff(Surv(time, event) ~ grp, data = dd)
  plr <- 1 - pchisq(lr$chisq, length(lr$n) - 1)
  medtab <- summary(fit)$table; arms <- names(table(dd$grp))
  medv <- if (is.matrix(medtab)) setNames(medtab[, "median"], sub("grp=", "", rownames(medtab))) else setNames(medtab["median"], arms[1])
  hrtxt <- sprintf("Cox HR/SD = %.2f [%.2f-%.2f], p=%s", hr, ci[1], ci[2], fmtp(pcox))
  lrtxt <- sprintf("%slog-rank p = %s", if (stratified) "stratified " else "", fmtp(plr))
  say("[%s | %s] split=%s | %s | %s | %s", endpoint, sp, split_type,
      paste(arms, table(dd$grp), sep = "=", collapse = ","), lrtxt, hrtxt)
  SUMM[[paste(endpoint, sp)]] <<- data.table(endpoint = endpoint, species = sp, split = split_type,
      arm = arms, n = as.integer(table(dd$grp)), median_surv_months = round(as.numeric(medv[arms]), 1),
      logrank_p = signif(plr, 3), cox_HR_perSD = round(hr, 3), cox_lo = round(ci[1], 3), cox_hi = round(ci[2], 3), cox_p = signif(pcox, 3))
  ttl <- bquote(italic(.(sp)) ~ .(paste0("abundance and ", endpoint)))
  stat_lab <- paste0(lrtxt, "\n", hrtxt)   # rendered large inside the plot (lower-left), not in the caption
  cap <- sprintf("Pre-specified split; cohort-stratified Cox on continuous log10-abundance is the primary test (KM descriptive).\nEXPLORATORY (3 species x 2 endpoints; uncorrected). Separate NSCLC analysis from the CRC work.%s",
                 if (endpoint == "PFS") "\nPFS event reconstructed = (PFS<OS)|death." else "")
  km_plot(fit, dd, lvl, xlab, ttl, cap, stat_lab,
          file.path(res, sprintf("KM_%s_%s_NSCLC.png", endpoint, short(TARGET)[match(sp, SP_DISP)])))
}

## ============================ OS (pooled, stratified) =========================
os751 <- build_cohort("NSCLC_PRJNA751792", "OS", "OS_evt")[, cohort := "PRJNA751792"]
os102 <- build_cohort("NSCLC_PRJNA1023797", "OS", "Death")[, cohort := "PRJNA1023797"]
OSD <- rbind(os751, os102)[, event := evt_raw][is.finite(time) & time > 0 & event %in% c(0L, 1L)]
say("OS dataset: %d patients (%s) | deaths=%d", nrow(OSD), paste(OSD[, .N, cohort][, paste0(cohort, "=", N)], collapse = ", "), sum(OSD$event))
fwrite(OSD[, c("RunID", "cohort", "time", "event", SP_DISP), with = FALSE], file.path(res, "species_survival_OS_dataset_NSCLC.csv"))
for (sp in SP_DISP) analyse(OSD, sp, "OS", "Overall survival (months)", stratified = TRUE)

## ============================ PFS (1023797 only) ==============================
pf <- build_cohort("NSCLC_PRJNA1023797", "PFS", "Death")[, cohort := "PRJNA1023797"]
pf[, pfs_time := pmin(PFS, OS)]
pf[, event := as.integer((PFS < OS - 1e-9) | (Death == 1L))]
PFSD <- pf[, time := pfs_time][is.finite(time) & time > 0 & event %in% c(0L, 1L)]
say("PFS dataset (PRJNA1023797): %d patients | events(progression-or-death)=%d | %d PFS>OS censored at OS",
    nrow(PFSD), sum(PFSD$event), sum(pf$PFS > pf$OS + 1e-9, na.rm = TRUE))
fwrite(PFSD[, c("RunID", "cohort", "time", "event", SP_DISP), with = FALSE], file.path(res, "species_survival_PFS_dataset_NSCLC.csv"))
for (sp in SP_DISP) analyse(PFSD, sp, "PFS", "Progression-free survival (months)", stratified = FALSE)

fwrite(rbindlist(SUMM), file.path(res, "species_survival_summary_NSCLC.csv"))
writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_species_survival_km.txt"))
say("DONE. outputs in %s", res)
close(logcon)
