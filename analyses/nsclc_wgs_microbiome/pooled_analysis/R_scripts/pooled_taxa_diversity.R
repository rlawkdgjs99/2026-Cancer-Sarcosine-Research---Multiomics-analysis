#!/usr/bin/env Rscript
# ============================================================================
# POOLED (cross-cohort) species DA + diversity meta-analysis — R vs NR.
#   Discovery  = 3 French Ion-Torrent cohorts; Validation = Korean (PRJEB26531).
#   (1) Species DA: inverse-variance random-effects meta (metafor DL) of the
#       per-cohort MaAsLin2 coef +/- stderr; species present in >=2 discovery
#       cohorts; BH-FDR; I^2 heterogeneity; Korean sign concordance.
#   (2) Alpha diversity (Shannon): Stouffer directional p-combination (discovery).
#   (3) Beta diversity: tabulate per-cohort PERMANOVA (R^2, p).
# Inputs = per-cohort diversity/DA outputs (already computed & verified).
# Run: Rscript pooled_analysis/R_scripts/pooled_taxa_diversity.R
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel); library(metafor) })

# Run with the R working directory set to the analysis-folder root; paths below
# are relative to it.
pooled_dir <- file.path(".", "pooled_analysis"); main_dir <- "."
res <- file.path(pooled_dir, "results"); dir.create(res, showWarnings = FALSE, recursive = TRUE)
logcon <- file(file.path(res, "pooled_taxa_diversity_log.txt"), open = "wt")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logcon) }

DISC  <- c(PRJNA751792 = "NSCLC_PRJNA751792", PRJNA1023797 = "NSCLC_PRJNA1023797", PRJEB22863 = "NSCLC_RCC_PRJEB22863")
VALID <- c(PRJEB26531 = "NSCLC_PRJEB26531")
NDISC <- c(PRJNA751792 = 338, PRJNA1023797 = 421, PRJEB22863 = 65)   # R/NR-labeled NSCLC-WGS n (verified)
say("==== POOLED species DA + diversity meta | discovery=%s ====", paste(names(DISC), collapse=","))

## ======================= (1) species DA meta (discovery) =======================
sp_disc <- rbindlist(lapply(names(DISC), function(cn) {
  d <- fread(file.path(main_dir, DISC[cn], "results", paste0("DA_species_group_", cn, ".csv")))
  d[, cohort := cn]; d[, .(cohort, feature, coef, stderr, qval)]
}), fill = TRUE)
sp_disc <- sp_disc[is.finite(coef) & is.finite(stderr) & stderr > 0]
sp_keep <- sp_disc[, .N, by = feature][N >= 2, feature]
say("species meta-analysable (>=2 discovery cohorts): %d", length(sp_keep))

meta_sp <- rbindlist(lapply(sp_keep, function(s) {
  x <- sp_disc[feature == s]
  m <- tryCatch(rma(yi = coef, sei = stderr, data = x, method = "DL"), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  data.table(species = s, k_cohorts = nrow(x), pooled_coef = as.numeric(m$beta),
             ci_lb = m$ci.lb, ci_ub = m$ci.ub, pval = m$pval, I2 = m$I2,
             dir = ifelse(as.numeric(m$beta) > 0, "higher in R", "higher in NR"))
}))
meta_sp[, qval := p.adjust(pval, "BH")]
setorder(meta_sp, pval)
# Korean validation: sign concordance
kv <- fread(file.path(main_dir, VALID[1], "results", paste0("DA_species_group_", names(VALID), ".csv")))
meta_sp <- merge(meta_sp, kv[, .(species = feature, valid_coef = coef)], by = "species", all.x = TRUE)
meta_sp[, valid_concordant := ifelse(is.na(valid_coef), NA, sign(valid_coef) == sign(pooled_coef))]
fwrite(meta_sp, file.path(res, "pooled_species_meta.csv"))
say("species meta: %d | q<0.10: %d | q<0.05: %d", nrow(meta_sp), sum(meta_sp$qval<0.10), sum(meta_sp$qval<0.05))
say("top pooled species: %s", paste(head(meta_sp$species, 8), collapse = ", "))

# volcano of pooled species effects — colored by ICI-response GROUP (R vs NR).
# Filter 1 = significance q < 0.05 (BH-FDR). Color by pooled_coef SIGN (coef>0 =>
# higher in R, coef<0 => higher in NR). Effect-size gate COEF_GATE (filter 2) off by default.
SIG_Q <- 0.05; COEF_GATE <- 0
.LAB_R <- "Enriched in Responder (R)"; .LAB_NR <- "Enriched in Non-responder (NR)"; .LAB_NS <- "n.s. (q >= 0.05)"
meta_sp[, neglogq := -log10(qval)]
meta_sp[, hit := qval < SIG_Q & abs(pooled_coef) >= COEF_GATE]
meta_sp[, ici_grp := ifelse(!hit, .LAB_NS, ifelse(pooled_coef > 0, .LAB_R, .LAB_NR))]
meta_sp[, ici_grp := factor(ici_grp, levels = c(.LAB_R, .LAB_NR, .LAB_NS))]
vol_cols <- setNames(c("#1B7837", "#B2182B", "grey80"), c(.LAB_R, .LAB_NR, .LAB_NS))   # group scheme: R=green (favorable), NR=red (unfavorable) (display colour only)
top <- head(meta_sp[hit == TRUE][order(qval)], 15); if (nrow(top)) top[, lab := gsub("_"," ",species)]
pv <- ggplot(meta_sp, aes(pooled_coef, neglogq)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(SIG_Q), linetype = "dashed", color = "grey60") +
  geom_point(aes(color = ici_grp), size = 1.9, alpha = 0.85) +
  scale_color_manual(values = vol_cols, name = "ICI response", drop = FALSE) +
  labs(x = "Pooled MaAsLin2 coefficient   (<-- higher in NR    |    higher in R -->)", y = expression(-log[10]~"(meta q, BH)"),
       title = "Pooled species meta-analysis across 3 discovery cohorts (R vs NR)",
       subtitle = sprintf("Random-effects (DL) inverse-variance; colored by enrichment group (q < %.2f, BH); %d species in >=2 cohorts", SIG_Q, nrow(meta_sp))) +
  theme_classic(base_size = 16) + theme(plot.title = element_text(face = "bold", size = 16), legend.position = "right")
if (nrow(top)) pv <- pv + ggrepel::geom_text_repel(data = top, aes(pooled_coef, neglogq, label = lab),
      size = 3.8, fontface = "italic", max.overlaps = 20, min.segment.length = 0, segment.color = "grey70")
# ragg cannot open the non-ASCII (Korean) path on macOS; render to ASCII temp then copy.
.tmp_png <- tempfile(fileext = ".png")
ggsave(.tmp_png, pv, width = 7.6, height = 6, dpi = 300, bg = "white")
stopifnot(file.copy(.tmp_png, file.path(res, "pooled_species_meta_volcano.png"), overwrite = TRUE)); unlink(.tmp_png)

## ======================= (2) alpha diversity meta (Shannon) =======================
alpha <- rbindlist(lapply(names(DISC), function(cn) {
  a <- fread(file.path(main_dir, DISC[cn], "results", paste0("alpha_diversity_stats_", cn, ".csv")))
  a <- a[index == "Shannon"]; a[, cohort := cn]; a[, n := NDISC[cn]]
  a[, sgn := sign(median_R - median_NR)]; a
}))
z <- alpha$sgn * qnorm(1 - alpha$p_wilcox/2); w <- sqrt(alpha$n); Z <- sum(w*z)/sqrt(sum(w^2))
alpha_meta <- data.table(index = "Shannon", n_cohorts = nrow(alpha),
  per_cohort = paste(alpha$cohort, sprintf("%s p=%.3f", ifelse(alpha$sgn>0,"R","NR"), alpha$p_wilcox), collapse=" | "),
  combined_Z = Z, combined_p = 2*(1-pnorm(abs(Z))), pooled_dir = ifelse(Z>0,"higher in R","higher in NR"))
fwrite(alpha_meta, file.path(res, "pooled_alpha_Shannon_meta.csv"))
say("ALPHA Shannon meta (Stouffer, discovery): combined p=%.4f (%s)", alpha_meta$combined_p, alpha_meta$pooled_dir)

## ======================= (3) beta diversity tabulation =======================
beta <- rbindlist(lapply(c(names(DISC), names(VALID)), function(cn) {
  dir <- if (cn %in% names(DISC)) DISC[cn] else VALID[cn]
  b <- fread(file.path(main_dir, dir, "results", paste0("beta_permanova_", cn, ".csv")))
  g <- b[term == "group"]
  data.table(cohort = cn, role = ifelse(cn %in% names(DISC), "discovery", "validation"),
             R2 = g$R2, F = g$F, p = g$p)
}))
fwrite(beta, file.path(res, "pooled_beta_permanova_table.csv"))
say("BETA PERMANOVA per cohort: %s", paste(beta$cohort, sprintf("R2=%.3f,p=%.3f", beta$R2, beta$p), collapse=" | "))

writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_pooled_taxa_diversity.txt"))
say("DONE. outputs in %s", res)
close(logcon)
