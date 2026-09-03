#!/usr/bin/env Rscript
# ============================================================================
# POOLED (cross-cohort) sarcosine meta-analysis — ICI response (R vs NR).
# Strategy (locked): NO naive pooling. Combine per-cohort effects.
#   Discovery  = 3 French Ion-Torrent cohorts (PRJNA751792, PRJNA1023797, PRJEB22863-NSCLC)
#   Validation = Korean Illumina cohort (PRJEB26531), reported separately.
# Inputs = per-cohort sarcosine outputs (already computed & verified).
#   (1) Individual sarcosine KOs: inverse-variance random-effects meta of the
#       per-cohort MaAsLin2 coef +/- stderr (metafor, DerSimonian-Laird) + I^2.
#   (2) Pathway scores: directional p-combination (Stouffer's Z, weighted by
#       sqrt(n)) across discovery cohorts — robust for the sparse score data;
#       effect-size SMD also reported where estimable.
# Run: Rscript pooled_analysis/R_scripts/pooled_sarcosine.R
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(metafor) })

# Run with the R working directory set to the analysis-folder root; paths below
# are relative to it.
pooled_dir <- file.path(".", "pooled_analysis")
main_dir   <- "."
res <- file.path(pooled_dir, "results"); dir.create(res, showWarnings = FALSE, recursive = TRUE)
logcon <- file(file.path(res, "pooled_sarcosine_log.txt"), open = "wt")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logcon) }
# ragg (ggsave default) silently fails to write the non-ASCII (Korean) project path on macOS:
# render to an ASCII tempfile, then file.copy to the real path (bulletproof, matches sibling pooled scripts).
save_png <- function(p, out, w, h) { tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = w, height = h, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out, overwrite = TRUE)); unlink(tmp) }

DISC <- c(PRJNA751792 = "NSCLC_PRJNA751792",
          PRJNA1023797 = "NSCLC_PRJNA1023797",
          PRJEB22863 = "NSCLC_RCC_PRJEB22863")        # discovery (French Ion Torrent)
VALID <- c(PRJEB26531 = "NSCLC_PRJEB26531")            # external validation (Korea Illumina)
say("==== POOLED sarcosine meta-analysis | discovery=%s | validation=%s ====",
    paste(names(DISC), collapse=","), paste(names(VALID), collapse=","))

rd_ko    <- function(cn, dir) { f <- file.path(main_dir, dir, "results", paste0("sarcosine_KO_DA_", cn, ".csv")); if (!file.exists(f)) return(NULL); d <- fread(f); d[, cohort := cn]; d }
rd_score <- function(cn, dir) fread(file.path(main_dir, dir, "results", paste0("sarcosine_scores_", cn, ".csv")))

## =================== (1) individual KO meta (discovery) ===================
ko_disc <- rbindlist(lapply(names(DISC), function(cn) rd_ko(cn, DISC[cn])), fill = TRUE)
ko_disc <- ko_disc[is.finite(coef) & is.finite(stderr) & stderr > 0]
kos <- ko_disc[, .N, by = feature][N >= 2, feature]                # KO testable in >=2 discovery cohorts
say("KOs meta-analysable (>=2 discovery cohorts): %d -> %s", length(kos), paste(kos, collapse=", "))

meta_ko <- rbindlist(lapply(kos, function(k) {
  s <- ko_disc[feature == k]
  m <- tryCatch(rma(yi = coef, sei = stderr, data = s, method = "DL"), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  data.table(KO = k, k_cohorts = nrow(s),
             pooled_coef = as.numeric(m$beta), ci_lb = m$ci.lb, ci_ub = m$ci.ub,
             pval = m$pval, I2 = m$I2, tau2 = m$tau2,
             dir = ifelse(as.numeric(m$beta) > 0, "higher in R", "higher in NR"))
}))
if (nrow(meta_ko)) {
  meta_ko[, qval := p.adjust(pval, "BH")]
  setorder(meta_ko, pval)
  # add validation (Korean) coef for concordance
  kv <- rd_ko(names(VALID), VALID[1])
  if (!is.null(kv)) meta_ko <- merge(meta_ko, kv[, .(KO = feature, valid_coef = coef, valid_q = qval)], by = "KO", all.x = TRUE)
  fwrite(meta_ko, file.path(res, "pooled_sarcosine_KO_meta.csv"))
  say("pooled KO meta: %d KOs | q<0.10: %d | q<0.05: %d", nrow(meta_ko), sum(meta_ko$qval<0.10), sum(meta_ko$qval<0.05))
  # forest (ggplot) of pooled KO effects, annotated with each KO's KEGG sarcosine role
  koset <- fread(file.path(main_dir, "sarcosine_KO_set.csv"))
  meta_ko <- merge(meta_ko, koset[, .(KO, role)], by = "KO", all.x = TRUE)   # role added to the PLOT only (the CSV was already written above, so its md5 is unchanged)
  meta_ko[, role_lab := c(degradation = "Degradation", production = "Production", "degradation;production" = "Bidirectional")[role]]
  meta_ko[, KOf := factor(KO, levels = meta_ko[order(pooled_coef)]$KO)]
  meta_ko[, sig := ifelse(qval < 0.05, "q<0.05", ifelse(qval < 0.10, "q<0.10", "n.s."))]
  xr <- max(meta_ko$ci_ub) + 0.10
  ROLE_COL <- c(Degradation = "#117733", Production = "#E6A000", Bidirectional = "#882255")   # green / amber / purple
  p <- ggplot(meta_ko, aes(pooled_coef, KOf, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_segment(aes(x = ci_lb, xend = ci_ub, y = KOf, yend = KOf), linewidth = 0.5) +
    geom_point(size = 2.6) +
    geom_text(data = meta_ko[role_lab == "Degradation"],   aes(x = xr, y = KOf, label = role_lab), inherit.aes = FALSE, color = ROLE_COL["Degradation"],   hjust = 0, size = 3.1, fontface = "bold") +
    geom_text(data = meta_ko[role_lab == "Production"],    aes(x = xr, y = KOf, label = role_lab), inherit.aes = FALSE, color = ROLE_COL["Production"],    hjust = 0, size = 3.1, fontface = "bold") +
    geom_text(data = meta_ko[role_lab == "Bidirectional"], aes(x = xr, y = KOf, label = role_lab), inherit.aes = FALSE, color = ROLE_COL["Bidirectional"], hjust = 0, size = 3.1, fontface = "bold") +
    scale_color_manual(values = c("q<0.05"="#D55E00","q<0.10"="#E69F00","n.s."="grey55"), name = NULL) +
    scale_x_continuous(limits = c(min(meta_ko$ci_lb) - 0.05, xr + 1.0)) +
    labs(x = "Pooled MaAsLin2 coefficient  (← NR    |    R →)", y = NULL,
         title = "Pooled sarcosine-KO meta-analysis (3 discovery cohorts, RE/DL)",
         subtitle = "Inverse-variance random-effects; bars = 95% CI.  Right label = KEGG sarcosine role (Bidirectional KOs are excluded from the pathway scores)") +
    theme_classic(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 12), plot.subtitle = element_text(size = 8.5))
  save_png(p, file.path(res, "pooled_sarcosine_KO_meta.png"), 9.2, 0.45*nrow(meta_ko)+2)
}

## =================== (2) pathway-score meta (discovery) ===================
score_names <- c("degradation","production","prod_deg_log2ratio")
# per-cohort per-score: n, Wilcoxon p, direction sign, SMD (Hedges g)
percoh <- rbindlist(lapply(names(DISC), function(cn) {
  sc <- rd_score(cn, DISC[cn])
  rbindlist(lapply(score_names, function(s) {
    r <- sc[sc$group=="R"][[s]]; nr <- sc[sc$group=="NR"][[s]]
    r <- r[is.finite(r)]; nr <- nr[is.finite(nr)]     # drop NA (log-ratio undefined where a score is 0)
    wp <- wilcox.test(sc[[s]] ~ sc$group)$p.value      # formula method drops NA automatically
    sgn <- sign(mean(r) - mean(nr))
    sp <- sqrt(((length(r)-1)*var(r) + (length(nr)-1)*var(nr)) / (length(r)+length(nr)-2))
    smd <- if (is.finite(sp) && sp > 0) (mean(r)-mean(nr))/sp else NA_real_
    data.table(cohort=cn, score=s, n=length(r)+length(nr), nR=length(r), nNR=length(nr),
               mean_R=mean(r), mean_NR=mean(nr), wilcox_p=wp, sign=sgn, SMD=smd)
  }))
}))
fwrite(percoh, file.path(res, "pooled_sarcosine_score_percohort.csv"))

# Stouffer's Z (directional), weighted by sqrt(n), across discovery cohorts
stouffer <- rbindlist(lapply(score_names, function(s) {
  d <- percoh[score==s]
  z <- d$sign * qnorm(1 - d$wilcox_p/2)          # directional z per cohort
  w <- sqrt(d$n)
  Z <- sum(w*z)/sqrt(sum(w^2))
  data.table(score=s, n_cohorts=nrow(d), concordant_dir=max(table(d$sign)),
             combined_Z=Z, combined_p=2*(1-pnorm(abs(Z))),
             pooled_dir=ifelse(Z>0,"higher in R","higher in NR"),
             percoh_dirs=paste(ifelse(d$sign>0,"R","NR"), collapse="/"),
             percoh_p=paste(signif(d$wilcox_p,2), collapse="/"))
}))
stouffer[, combined_q := p.adjust(combined_p, "BH")]
fwrite(stouffer, file.path(res, "pooled_sarcosine_score_meta.csv"))
say("SCORE meta (Stouffer, discovery): %s",
    paste(stouffer$score, sprintf("p=%.3g(%s,%d/3 concordant)", stouffer$combined_p, stouffer$pooled_dir, stouffer$concordant_dir), sep=" ", collapse=" | "))

# SMD forest per score (metafor RE over discovery) where estimable
smd_meta <- rbindlist(lapply(score_names, function(s) {
  d <- percoh[score==s & is.finite(SMD)]
  if (nrow(d) < 2) return(NULL)
  # approximate SE of SMD (Hedges): sqrt((nR+nNR)/(nR*nNR) + SMD^2/(2*(nR+nNR)))
  d[, se := sqrt((nR+nNR)/(nR*nNR) + SMD^2/(2*(nR+nNR)))]
  m <- tryCatch(rma(yi = SMD, sei = se, data = d, method = "DL"), error=function(e) NULL)
  if (is.null(m)) return(NULL)
  data.table(score=s, pooled_SMD=as.numeric(m$beta), ci_lb=m$ci.lb, ci_ub=m$ci.ub, p=m$pval, I2=m$I2)
}))
if (nrow(smd_meta)) {
  fwrite(smd_meta, file.path(res, "pooled_sarcosine_score_SMD_meta.csv"))
  smd_meta[, scoref := factor(score, levels=score_names)]
  smd_meta[, sig := ifelse(p<0.05,"p<0.05",ifelse(p<0.10,"p<0.10","n.s."))]
  pp <- ggplot(smd_meta, aes(pooled_SMD, scoref, color=sig)) +
    geom_vline(xintercept=0, linetype="dashed", color="grey60") +
    geom_segment(aes(x=ci_lb, xend=ci_ub, y=scoref, yend=scoref), linewidth=0.6) +
    geom_point(size=3) +
    scale_color_manual(values=c("p<0.05"="#D55E00","p<0.10"="#E69F00","n.s."="grey55"), name=NULL) +
    labs(x="Pooled standardized mean difference  (← NR    |    R →)", y=NULL,
         title="Pooled sarcosine pathway-score meta-analysis (3 discovery cohorts)",
         subtitle="Random-effects (DL) of per-cohort Hedges' g (R vs NR); bars = 95% CI") +
    theme_classic(base_size=12) + theme(plot.title=element_text(face="bold", size=12))
  save_png(pp, file.path(res, "pooled_sarcosine_score_meta.png"), 7.5, 3.4)
}

writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_pooled_sarcosine.txt"))
say("DONE. outputs in %s", res)
close(logcon)
