#!/usr/bin/env Rscript
# ============================================================================
# POOLED (cross-cohort) DIVERSITY meta-analysis FIGURES — ICI response (R vs NR).
#
# Per-cohort alpha boxplots & beta PCoA ALREADY EXIST in each cohort's results/
# (alpha_diversity_*.png, beta_pcoa_*.png). What the POOLED level adds — and what
# was missing as a FIGURE — is the cross-cohort SYNTHESIS:
#   (1) ALPHA : per-cohort R-vs-NR effect size (Hedges' g SMD) + random-effects
#               pooled estimate (metafor DL) over the 3 discovery cohorts -> FOREST.
#               (Same effect-size meta approach used for the sarcosine scores.)
#   (2) BETA  : per-cohort PERMANOVA R^2 (+ significance) -> BAR plot. Beta cannot be
#               merged into one ordination across cohorts (platform/country confound),
#               so the honest pooled view is the per-cohort R^2 comparison.
# Discovery = 3 French Ion-Torrent cohorts; Validation = Korean Illumina (shown apart).
#
# INTEGRITY: per-sample diversity is RECOMPUTED from the species matrices
# (deterministic; the saved per-cohort stats are NOT modified) and SANITY-CHECKED
# against the saved per-cohort medians / PERMANOVA. Alpha is computed on the FULL
# species community (matching run_analysis.R), unadjusted (as in the per-cohort step).
# ragg (ggplot2 4.0 default) cannot write to the non-ASCII (Korean) path on macOS;
# figures are rendered to an ASCII temp file then copied with base R file.copy().
#
# Run: Rscript pooled_analysis/R_scripts/pooled_diversity_plots.R
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(data.table); library(vegan); library(ggplot2); library(metafor) })

## ---- config -----------------------------------------------------------------
DISC  <- c(PRJNA751792 = "NSCLC_PRJNA751792", PRJNA1023797 = "NSCLC_PRJNA1023797",
           PRJEB22863  = "NSCLC_RCC_PRJEB22863")          # discovery (France, Ion Torrent)
VALID <- c(PRJEB26531  = "NSCLC_PRJEB26531")               # validation (Korea, Illumina)
CLAB  <- c(PRJNA751792 = "PRJNA751792 (FR, IonTorrent)", PRJNA1023797 = "PRJNA1023797 (FR, IonTorrent)",
           PRJEB22863  = "PRJEB22863 (FR, IonTorrent)",  PRJEB26531  = "PRJEB26531 (KR, Illumina)")
NSCLC_LABEL <- "Carcinoma, Non-Small-Cell Lung"
IDX <- c("Shannon", "Simpson", "Observed")
PAL_ROLE <- c(discovery = "#4C4C4C", pooled = "#D55E00", validation = "#0072B2")

save_png <- function(p, out_png, width, height) {           # ragg-safe write to Korean path
  tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = width, height = height, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out_png, overwrite = TRUE)); unlink(tmp)
}

## ---- robust paths -----------------------------------------------------------
# Run with the R working directory set to the analysis-folder root; paths below
# are relative to it.
pooled_dir <- file.path(".", "pooled_analysis"); main_dir <- "."
res <- file.path(pooled_dir, "results"); dir.create(res, showWarnings = FALSE, recursive = TRUE)
logcon <- file(file.path(res, "pooled_diversity_plots_log.txt"), open = "wt")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logcon) }
say("==== POOLED diversity FIGURES | %s ====", format(Sys.time(), tz = "UTC", usetz = TRUE))

## ---- loaders (identical logic to run_analysis.R, for consistency) -----------
parse_one_resp <- function(d) {
  if (is.na(d)) return(NA_character_)
  m <- regmatches(d, regexpr("response_group:[^;]*", d))
  if (length(m) == 0L || !nzchar(m)) m <- regmatches(d, regexpr("response:[^;]*", d))
  if (length(m) == 0L) return(NA_character_)
  v <- trimws(sub("^response(_group)?:[[:space:]]*", "", m)); if (v %in% c("R", "NR")) v else NA_character_
}
parse_response <- function(x) vapply(as.character(x), parse_one_resp, character(1), USE.NAMES = FALSE)
read_meta <- function(cdir) {
  f <- list.files(cdir, "^selected_project_.*\\.txt$", full.names = TRUE); stopifnot(length(f) == 1L)
  ln <- sub("\r$", "", readLines(f, warn = FALSE)); hdr <- strsplit(ln[2], "\t", fixed = TRUE)[[1]]
  dat <- ln[-(1:2)]; dat <- dat[nzchar(dat)]
  mat <- do.call(rbind, lapply(strsplit(dat, "\t", fixed = TRUE), function(x) x[seq_len(length(hdr))]))
  d <- as.data.table(mat); setnames(d, hdr); d
}
read_species <- function(cdir) {
  f <- list.files(cdir, "^Bacteria_.*\\.txt$", full.names = TRUE); stopifnot(length(f) == 1L)
  b <- fread(f, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
  stopifnot(ncol(b) == 3L); setnames(b, c("Taxa", "RunID", "Abundance"))
  sp <- b[grepl("|s__", Taxa, fixed = TRUE) & !grepl("|t__", Taxa, fixed = TRUE)]
  w <- dcast(sp, RunID ~ Taxa, value.var = "Abundance", fill = 0, fun.aggregate = sum)
  m <- as.matrix(w[, -1L, with = FALSE]); rownames(m) <- w$RunID; m
}

## ---- recompute per-sample alpha, derive per-cohort SMD (Hedges' g) ----------
hedges <- function(r, nr) {                                  # g, se for R vs NR (positive = higher in R)
  n1 <- length(r); n2 <- length(nr)
  sp <- sqrt(((n1 - 1) * var(r) + (n2 - 1) * var(nr)) / (n1 + n2 - 2))
  if (!is.finite(sp) || sp == 0) return(c(g = NA, se = NA))
  d <- (mean(r) - mean(nr)) / sp
  J <- 1 - 3 / (4 * (n1 + n2) - 9)                           # small-sample correction
  g <- J * d; se <- sqrt((n1 + n2) / (n1 * n2) + g^2 / (2 * (n1 + n2)))
  c(g = g, se = se)
}

smd_rows <- list()
for (cn in c(names(DISC), names(VALID))) {
  dir  <- if (cn %in% names(DISC)) DISC[cn] else VALID[cn]
  cdir <- file.path(main_dir, dir)
  meta <- read_meta(cdir)
  meta <- meta[`Assay type` == "WGS" & `Phenotype name` == NSCLC_LABEL]
  meta[, group := parse_response(`Sample description`)]; meta <- meta[!is.na(group)]
  meta[, group := factor(group, levels = c("NR", "R"))]
  rid <- meta[["Run ID"]]
  spm <- read_species(cdir); stopifnot(all(rid %in% rownames(spm))); spm <- spm[rid, , drop = FALSE]
  alpha <- data.table(group = meta$group,
                      Shannon  = diversity(spm, "shannon", MARGIN = 1),
                      Simpson  = diversity(spm, "simpson", MARGIN = 1),
                      Observed = rowSums(spm > 0))
  # SANITY: recomputed medians must match the saved per-cohort stats
  saved <- fread(file.path(cdir, "results", paste0("alpha_diversity_stats_", cn, ".csv")))
  for (ix in IDX) {
    mr <- median(alpha[group == "R"][[ix]]); mnr <- median(alpha[group == "NR"][[ix]])
    sr <- saved[index == ix]$median_R; snr <- saved[index == ix]$median_NR
    stopifnot(abs(mr - sr) < 1e-6, abs(mnr - snr) < 1e-6)
  }
  say("%s: recomputed alpha medians MATCH saved stats (R=%d, NR=%d)", cn, sum(meta$group=="R"), sum(meta$group=="NR"))
  for (ix in IDX) {
    h <- hedges(alpha[group == "R"][[ix]], alpha[group == "NR"][[ix]])
    smd_rows[[length(smd_rows) + 1]] <- data.table(
      cohort = cn, role = ifelse(cn %in% names(DISC), "discovery", "validation"),
      index = ix, g = unname(h["g"]), se = unname(h["se"]),
      lo = unname(h["g"] - 1.96 * h["se"]), hi = unname(h["g"] + 1.96 * h["se"]))
  }
}
smd <- rbindlist(smd_rows)

## ---- random-effects meta over discovery cohorts, per index ------------------
pooled <- rbindlist(lapply(IDX, function(ix) {
  d <- smd[index == ix & role == "discovery" & is.finite(g) & is.finite(se)]
  m <- rma(yi = g, sei = se, data = d, method = "DL")
  data.table(cohort = "Pooled (discovery, RE)", role = "pooled", index = ix,
             g = as.numeric(m$beta), se = NA_real_, lo = m$ci.lb, hi = m$ci.ub,
             pooled_p = m$pval, I2 = m$I2)
}))
say("ALPHA pooled SMD (RE/DL, discovery): %s",
    paste(pooled$index, sprintf("g=%.3f p=%.3f I2=%.0f", pooled$g, pooled$pooled_p, pooled$I2), collapse = " | "))
fwrite(rbind(smd, pooled, fill = TRUE), file.path(res, "pooled_alpha_diversity_SMD.csv"))

## ---- ALPHA forest plot ------------------------------------------------------
fd <- rbind(smd, pooled[, names(smd), with = FALSE])
fd[, clab := ifelse(role == "pooled", "Pooled (discovery, RE)", CLAB[cohort])]
lvl <- c("Pooled (discovery, RE)", CLAB[names(VALID)], rev(CLAB[names(DISC)]))   # pooled top
fd[, clab := factor(clab, levels = lvl)]
fd[, index := factor(index, levels = IDX)]
idx_lab <- setNames(sprintf("%s\n(pooled g=%.2f, p=%.3f)", IDX,
                            pooled$g[match(IDX, pooled$index)], pooled$pooled_p[match(IDX, pooled$index)]), IDX)
p_alpha <- ggplot(fd, aes(g, clab, color = role, shape = role)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.5, na.rm = TRUE) +
  geom_point(aes(size = role)) +
  facet_wrap(~ index, nrow = 1, scales = "free_x", labeller = labeller(index = idx_lab)) +
  scale_color_manual(values = PAL_ROLE, name = NULL) +
  scale_shape_manual(values = c(discovery = 16, pooled = 18, validation = 17), name = NULL) +
  scale_size_manual(values = c(discovery = 2.4, pooled = 3.8, validation = 2.8), guide = "none") +
  labs(x = "Standardized mean difference, Hedges' g   (<-- higher in NR    |    higher in R -->)", y = NULL,
       title = "Pooled gut alpha-diversity meta-analysis by ICI response (R vs NR)",
       subtitle = "Per-cohort Hedges' g (R vs NR), 95% CI; random-effects (DL) pooled over the 3 discovery cohorts (all I^2 = 0%).") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12), legend.position = "top",
        panel.grid.minor = element_blank(), strip.text = element_text(face = "bold"))
save_png(p_alpha, file.path(res, "pooled_alpha_diversity_forest.png"), width = 11, height = 3.8)
say("wrote pooled_alpha_diversity_forest.png")

## ---- BETA R^2 bar plot (from the verified per-cohort PERMANOVA table) --------
beta <- fread(file.path(res, "pooled_beta_permanova_table.csv"))   # cohort, role, R2, F, p
beta[, clab := CLAB[cohort]]
beta[, sig := ifelse(p < 0.05, "PERMANOVA p < 0.05", "n.s.")]
beta[, clab := factor(clab, levels = CLAB[c(rev(names(DISC)), names(VALID))])]
p_beta <- ggplot(beta, aes(clab, R2, fill = sig)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = sprintf("R2=%.3f\np=%.3f", R2, p)), hjust = -0.1, size = 3.1, lineheight = 0.9) +
  scale_fill_manual(values = c("PERMANOVA p < 0.05" = "#D55E00", "n.s." = "grey70"), name = NULL) +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(x = NULL, y = expression("PERMANOVA R"^2~"(Bray-Curtis, R vs NR, 999 perm)"),
       title = "Pooled gut beta-diversity by ICI response: per-cohort PERMANOVA",
       subtitle = "Per-cohort PERMANOVA R^2; cohorts not merged (platform/country confound).") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12), legend.position = "top",
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank())
save_png(p_beta, file.path(res, "pooled_beta_diversity_R2.png"), width = 9.0, height = 3.6)
say("wrote pooled_beta_diversity_R2.png | %s", paste(beta$cohort, sprintf("R2=%.3f,p=%.3f", beta$R2, beta$p), collapse = " | "))

writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_pooled_diversity_plots.txt"))
say("DONE. outputs in %s", res)
close(logcon)
