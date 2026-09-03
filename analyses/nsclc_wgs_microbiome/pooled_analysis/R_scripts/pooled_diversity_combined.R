#!/usr/bin/env Rscript
# ============================================================================
# POOLED (combined-sample) DIVERSITY figures — the CONVENTIONAL view, done
# COHORT-AWARE. Complements pooled_diversity_plots.R (the meta-analytic
# forest / R^2 synthesis); here all WGS samples are pooled into ONE alpha
# boxplot and ONE Bray-Curtis ordination, as is conventional — but because
# the 4 cohorts differ by country/platform (3 France/Ion-Torrent discovery +
# 1 Korea/Illumina validation), that confound is made VISIBLE (colored by
# cohort) and ADJUSTED FOR in every test. A naive R-vs-NR pooled test would
# confound ICI response with platform; we never report that as the result.
#
#   ALPHA : box + jitter (R vs NR) for Shannon/Simpson/Observed over all 849
#           WGS samples; jitter colored by cohort. Test per index =
#           lm(index ~ cohort + response)  -> cohort-ADJUSTED response p
#           (primary), with naive Wilcoxon shown alongside for transparency.
#           NOTE: Observed richness is the most platform-sensitive index.
#   BETA  : combined Bray-Curtis (union of species, 0-filled, re-TSS per sample)
#           -> PCoA (classical MDS). Two panels of the SAME ordination:
#           (A) colored by cohort  -> shows platform/country dominates the axes;
#           (B) colored by R vs NR -> shows whether response separates globally.
#           Test = adonis2(bray ~ cohort + response, by="terms", 999 perm):
#           reports cohort R^2 (batch) AND response-after-cohort R^2/p (the
#           honest adjusted effect); naive adonis2(bray ~ response) shown for contrast.
#
# INTEGRITY: per-sample alpha is RECOMPUTED from the species matrices and
# SANITY-CHECKED against the saved per-cohort medians (saved stats unchanged).
# Species names are identical MetaPhlAn4 SGB strings across cohorts (verified:
# 914 species shared by all 4; union 2,598). ragg cannot write the non-ASCII
# (Korean) path on macOS -> render to ASCII temp then base R file.copy().
#
# Run: Rscript pooled_analysis/R_scripts/pooled_diversity_combined.R
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({
  library(data.table); library(vegan); library(ggplot2); library(cowplot)
})

## ---- config -----------------------------------------------------------------
COH <- c(PRJNA751792 = "NSCLC_PRJNA751792", PRJNA1023797 = "NSCLC_PRJNA1023797",
         PRJEB22863  = "NSCLC_RCC_PRJEB22863", PRJEB26531  = "NSCLC_PRJEB26531")
ROLE <- c(PRJNA751792 = "discovery", PRJNA1023797 = "discovery",
          PRJEB22863  = "discovery", PRJEB26531  = "validation")
CLAB <- c(PRJNA751792 = "PRJNA751792 (FR, IonTorrent)", PRJNA1023797 = "PRJNA1023797 (FR, IonTorrent)",
          PRJEB22863  = "PRJEB22863 (FR, IonTorrent)",  PRJEB26531  = "PRJEB26531 (KR, Illumina)")
NSCLC_LABEL <- "Carcinoma, Non-Small-Cell Lung"
IDX <- c("Shannon", "Simpson", "Observed")
RESP_PAL <- c(NR = "#B2182B", R = "#1B7837")                        # group scheme: R=green (favorable), NR=red (unfavorable) (display colour only)
COH_PAL  <- c(PRJNA751792 = "#E69F00", PRJNA1023797 = "#009E73",    # colorblind-safe (Okabe-Ito)
              PRJEB22863  = "#56B4E9", PRJEB26531  = "#CC79A7")

save_png <- function(p, out_png, width, height) {                  # ragg-safe write to Korean path
  tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = width, height = height, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out_png, overwrite = TRUE)); unlink(tmp)
}

# Run with the R working directory set to the analysis-folder root; paths below
# are relative to it.
pooled_dir <- file.path(".", "pooled_analysis"); main_dir <- "."
res <- file.path(pooled_dir, "results"); dir.create(res, showWarnings = FALSE, recursive = TRUE)
logcon <- file(file.path(res, "pooled_diversity_combined_log.txt"), open = "wt")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logcon) }
say("==== POOLED diversity COMBINED (cohort-aware) | %s ====", format(Sys.time(), tz = "UTC", usetz = TRUE))

## ---- loaders (identical to run_analysis.R / pooled_diversity_plots.R) --------
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

## ---- load all cohorts: per-sample alpha + species blocks --------------------
alpha_rows <- list(); blocks <- list(); meta_rows <- list()
for (cn in names(COH)) {
  cdir <- file.path(main_dir, COH[cn])
  meta <- read_meta(cdir)
  meta <- meta[`Assay type` == "WGS" & `Phenotype name` == NSCLC_LABEL]
  meta[, group := parse_response(`Sample description`)]; meta <- meta[!is.na(group)]
  meta[, group := factor(group, levels = c("NR", "R"))]
  rid <- meta[["Run ID"]]
  spm <- read_species(cdir); stopifnot(all(rid %in% rownames(spm))); spm <- spm[rid, , drop = FALSE]

  a <- data.table(sample = paste0(cn, "::", rid), cohort = cn, group = meta$group,
                  Shannon  = diversity(spm, "shannon", MARGIN = 1),
                  Simpson  = diversity(spm, "simpson", MARGIN = 1),
                  Observed = rowSums(spm > 0))
  # SANITY: recomputed medians must match the saved per-cohort stats
  saved <- fread(file.path(cdir, "results", paste0("alpha_diversity_stats_", cn, ".csv")))
  for (ix in IDX) {
    mr <- median(a[group == "R"][[ix]]); mnr <- median(a[group == "NR"][[ix]])
    stopifnot(abs(mr - saved[index == ix]$median_R) < 1e-6, abs(mnr - saved[index == ix]$median_NR) < 1e-6)
  }
  say("%s: alpha medians MATCH saved (R=%d, NR=%d)", cn, sum(meta$group == "R"), sum(meta$group == "NR"))

  rownames(spm) <- paste0(cn, "::", rid)
  alpha_rows[[cn]] <- a; blocks[[cn]] <- spm
  meta_rows[[cn]] <- data.table(sample = rownames(spm), cohort = cn, group = meta$group)
}
alphaDT <- rbindlist(alpha_rows)
metaA   <- rbindlist(meta_rows)
metaA[, cohort := factor(cohort, levels = names(COH))]
say("Pooled: %d WGS samples (R=%d, NR=%d) across %d cohorts",
    nrow(alphaDT), sum(alphaDT$group == "R"), sum(alphaDT$group == "NR"), length(COH))

## ====================  ALPHA: pooled box + cohort-adjusted test  =============
alpha_stats <- rbindlist(lapply(IDX, function(ix) {
  v <- alphaDT[[ix]]; g <- alphaDT$group; coh <- alphaDT$cohort
  fit <- lm(v ~ coh + g)                                            # response adjusted for cohort
  cf  <- summary(fit)$coefficients
  p_adj <- cf["gR", "Pr(>|t|)"]; beta_adj <- cf["gR", "Estimate"]   # gR>0 => higher in R (NR is ref)
  p_naive <- wilcox.test(v ~ g)$p.value                            # descriptive, unadjusted
  data.table(index = ix, beta_adj_RvsNR = beta_adj, p_adj_cohort = p_adj, p_naive_wilcox = p_naive,
             median_R = median(v[g == "R"]), median_NR = median(v[g == "NR"]),
             dir = ifelse(beta_adj > 0, "higher in R", "higher in NR"))
}))
say("ALPHA cohort-adjusted (lm ~ cohort + response): %s",
    paste(alpha_stats$index, sprintf("adjP=%.3g(%s) naiveP=%.3g", alpha_stats$p_adj_cohort,
          alpha_stats$dir, alpha_stats$p_naive_wilcox), collapse = " | "))

alphaL <- melt(alphaDT, id.vars = c("sample", "cohort", "group"),
               measure.vars = IDX, variable.name = "index", value.name = "value")
alphaL[, index := factor(index, levels = IDX)]
ann <- merge(alpha_stats[, .(index, p_adj_cohort, p_naive_wilcox, dir)],
             alphaL[, .(y = max(value) + 0.09 * (max(value) - min(value))), by = index], by = "index")
ann[, index := factor(index, levels = IDX)]
ann[, lab := sprintf("p = %.3g", p_adj_cohort)]   # p = cohort-adjusted lm(index ~ cohort + response); test stated in the figure legend, not the in-plot label

p_alpha <- ggplot(alphaL, aes(group, value)) +
  geom_boxplot(aes(fill = group), width = 0.6, alpha = 0.30, outlier.shape = NA, colour = "grey30") +
  geom_jitter(width = 0.18, height = 0, size = 0.8, alpha = 0.28, color = "grey40") +
  geom_text(data = ann, aes(x = 1.5, y = y, label = lab), inherit.aes = FALSE,
            size = 5.6, vjust = 1, fontface = "bold") +
  facet_wrap(~ index, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = RESP_PAL, guide = "none") +
  scale_x_discrete(limits = c("R", "NR")) +     # display order only: R left, NR right (factor/reference unchanged)
  labs(x = "ICI response group", y = "Alpha diversity (per sample)",
       title = "NSCLC ICI response: gut alpha-diversity (R vs NR)") +
  theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 23),
        legend.position = "none", panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 20),
        axis.title = element_text(size = 20), axis.text = element_text(size = 17))
save_png(p_alpha, file.path(res, "pooled_alpha_diversity_combined.png"), width = 11, height = 4.8)
say("wrote pooled_alpha_diversity_combined.png")

## ====================  BETA: combined Bray-Curtis -> PCoA  ====================
allsp <- sort(unique(unlist(lapply(blocks, colnames))))
big <- matrix(0, nrow = nrow(metaA), ncol = length(allsp),
              dimnames = list(metaA$sample, allsp))
for (cn in names(blocks)) { b <- blocks[[cn]]; big[rownames(b), colnames(b)] <- b }
stopifnot(all(rownames(big) == metaA$sample))                      # alignment
big <- big / rowSums(big)                                          # re-TSS each sample to 1
stopifnot(max(abs(rowSums(big) - 1)) < 1e-9)
say("BETA combined matrix: %d samples x %d species (union); all rows TSS=1", nrow(big), ncol(big))

bray <- vegdist(big, method = "bray")
pc   <- cmdscale(bray, k = 2, eig = TRUE)
ve   <- pc$eig / sum(pc$eig[pc$eig > 0]) * 100                     # % of positive eigenvalues
pcoa <- data.table(sample = rownames(big), Axis1 = pc$points[, 1], Axis2 = pc$points[, 2])
pcoa <- merge(pcoa, metaA, by = "sample")
pcoa[, cohort := factor(cohort, levels = names(COH))]

set.seed(42)                                                       # reproducible permutations
ad_full  <- adonis2(bray ~ cohort + group, data = metaA, by = "terms", permutations = 999)
set.seed(42)
ad_naive <- adonis2(bray ~ group, data = metaA, by = "terms", permutations = 999)
r2_coh <- ad_full["cohort", "R2"];  p_coh <- ad_full["cohort", "Pr(>F)"]
r2_grp <- ad_full["group",  "R2"];  p_grp <- ad_full["group",  "Pr(>F)"]
r2_nai <- ad_naive["group", "R2"];  p_nai <- ad_naive["group", "Pr(>F)"]
say("BETA adonis2(~cohort+response): cohort R2=%.3f p=%.3f | response|cohort R2=%.4f p=%.3f || naive(~response) R2=%.4f p=%.3f",
    r2_coh, p_coh, r2_grp, p_grp, r2_nai, p_nai)

base_pcoa <- function() list(
  geom_point(size = 1.5, alpha = 0.75),
  labs(x = sprintf("PCoA 1 (%.1f%%)", ve[1]), y = sprintf("PCoA 2 (%.1f%%)", ve[2])),
  theme_bw(base_size = 16),
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 14), plot.subtitle = element_text(size = 11)))
pA <- ggplot(pcoa, aes(Axis1, Axis2, color = cohort)) + base_pcoa() +
  stat_ellipse(level = 0.95, linewidth = 0.4) +
  scale_color_manual(values = COH_PAL, name = "Cohort", labels = CLAB) +
  labs(title = "Colored by cohort (platform/country)",
       subtitle = sprintf("Cohorts overlap heavily; batch effect small: cohort R2 = %.3f (p = %.3f)", r2_coh, p_coh))
pB <- ggplot(pcoa, aes(Axis1, Axis2, color = group)) + base_pcoa() +
  stat_ellipse(level = 0.95, linewidth = 0.5) +
  scale_color_manual(values = RESP_PAL, name = "ICI response") +
  labs(title = "Colored by ICI response (R vs NR)",
       subtitle = sprintf("Response|cohort R2 = %.4f (p = %.3f); naive ~response R2 = %.4f (p = %.3f)",
                          r2_grp, p_grp, r2_nai, p_nai))
title <- ggdraw() + draw_label(
  "Pooled gut beta-diversity (Bray-Curtis PCoA) by ICI response - all 849 WGS samples, union of 2,598 species",
  fontface = "bold", size = 16, x = 0.01, hjust = 0)
cap <- ggdraw() + draw_label(
  "Same ordination, two colorings. PERMANOVA adonis2(Bray ~ cohort + response, 999 perm): response tested AFTER conditioning on cohort. Cohorts not corrected/merged otherwise - this is a descriptive pooled view; the meta-analytic forest / per-cohort PERMANOVA remain the primary inference.",
  size = 10, x = 0.01, hjust = 0, y = 0.7)
body  <- plot_grid(pA, pB, nrow = 1, align = "h", rel_widths = c(1, 1))
p_beta <- plot_grid(title, body, cap, ncol = 1, rel_heights = c(0.08, 1, 0.10))
save_png(p_beta, file.path(res, "pooled_beta_diversity_PCoA_combined.png"), width = 19, height = 5.2)
say("wrote pooled_beta_diversity_PCoA_combined.png")

## ---- standalone single-panel PCoA by ICI response (drop-in for Main Fig 3a) --
## Same ordination as panel B above, as its own figure with a self-contained
## title (no panel-letter prefix). Reuses pcoa/adonis2 already computed.
pB_solo <- ggplot(pcoa, aes(Axis1, Axis2, color = group)) + base_pcoa() +
  stat_ellipse(level = 0.95, linewidth = 0.5) +
  scale_color_manual(values = RESP_PAL, name = "ICI response") +
  labs(title = "Gut composition by ICI response (PCoA)",
       subtitle = sprintf("Response|cohort R2 = %.4f (p = %.3f); n = %d WGS samples",
                          r2_grp, p_grp, nrow(pcoa))) +
  theme(plot.title = element_text(size = 22))
save_png(pB_solo, file.path(res, "pooled_beta_diversity_PCoA_byresponse.png"), width = 7.2, height = 5.6)
say("wrote pooled_beta_diversity_PCoA_byresponse.png")

## ---- save stats -------------------------------------------------------------
fwrite(alpha_stats, file.path(res, "pooled_alpha_combined_stats.csv"))
fwrite(data.table(
  term = c("cohort", "response|cohort", "response (naive)"),
  R2   = c(r2_coh, r2_grp, r2_nai), p = c(p_coh, p_grp, p_nai),
  model = c("adonis2(bray ~ cohort + response)", "adonis2(bray ~ cohort + response)", "adonis2(bray ~ response)")
), file.path(res, "pooled_beta_combined_permanova.csv"))
writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_pooled_diversity_combined.txt"))
say("DONE. outputs in %s", res)
close(logcon)
