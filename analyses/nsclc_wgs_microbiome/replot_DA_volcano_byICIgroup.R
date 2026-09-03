#!/usr/bin/env Rscript
# ============================================================================
# Re-plot species DA volcano figures, COLORED BY ICI-RESPONSE GROUP (R vs NR).
#
# WHY: the original volcanoes colored points by significance TIER (q<0.05 /
# q<0.10 / n.s.), which does NOT surface this study's actual contrast — gut
# species ENRICHED IN RESPONDERS (R) vs NON-RESPONDERS (NR) to immune-checkpoint
# inhibitors in NSCLC. Here, FDR-significant species are colored by the response
# group they are enriched in (the biological purpose of the analysis).
#
# INTEGRITY: this script ONLY READS the already-verified differential-abundance
# result CSVs and RE-DRAWS the figures. It does NOT re-run MaAsLin2 / metafor or
# recompute any statistic — the underlying numbers are untouched. (Verified by
# md5 of the input CSVs before/after running this script.)
#
# FILTERS (a volcano has two):
#   filter 1  significance : q < 0.05 (BH-FDR)                    [user's choice]
#   filter 2  effect size  : |coef| >= COEF_GATE ; default 0 = OFF.
#     A magnitude cutoff on standardized MaAsLin2 coefficients is arbitrary and
#     would hide genuinely q<0.05 species, so it is OFF by default ("rigorous"
#     here = strict FDR). Set COEF_GATE <- 1.0 for a strict dual-threshold volcano.
#   direction: NR is the model reference level => coef>0 higher in R, coef<0 higher in NR.
#
# OUTPUTS (overwrites the canonical volcano PNGs in place; statistics untouched):
#   <cohort>/results/DA_volcano_<cohort>.png                            (4 cohorts)
#   pooled_analysis/results/pooled_species_meta_volcano.png
# (The same figures are also produced by re-running the source run_analysis.R /
#  pooled_taxa_diversity.R; this script just regenerates them WITHOUT re-running
#  MaAsLin2/metafor, reading the verified DA result CSVs only.)
#
# Run: Rscript replot_DA_volcano_byICIgroup.R
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggrepel) })

SIG_Q     <- 0.05          # filter 1: BH-FDR significance threshold
COEF_GATE <- 0             # filter 2: |coef| gate (0 = off; set 1.0 for strict dual-threshold)
PAL       <- c(NR = "#B2182B", R = "#1B7837")            # group scheme: R=green (favorable), NR=red (unfavorable) (display colour only; matches cohort scripts)
LAB_R  <- "Enriched in Responder (R)"
LAB_NR <- "Enriched in Non-responder (NR)"
LAB_NS <- "n.s. (q >= 0.05)"   # ASCII only: installed fonts lack the >= / arrow glyphs
VOLCOLS <- setNames(c(unname(PAL["R"]), unname(PAL["NR"]), "grey80"), c(LAB_R, LAB_NR, LAB_NS))

# Run with the R working directory set to the analysis-folder root; paths below
# are relative to it.
main_dir <- "."

# Group-colored volcano. dt must have columns: feature, coef, qval.
make_volcano <- function(dt, title, subtitle, out_png,
                         point_size = 2, base_size = 16, label_n = 12,
                         width = 8.8, height = 5.6) {   # widened so the long subtitle is not clipped
  d <- copy(dt)
  d <- d[is.finite(coef) & is.finite(qval)]
  d[, neglogq := -log10(qval)]
  d[, hit := qval < SIG_Q & abs(coef) >= COEF_GATE]                  # passes BOTH filters
  d[, ici_grp := ifelse(!hit, LAB_NS, ifelse(coef > 0, LAB_R, LAB_NR))]
  d[, ici_grp := factor(ici_grp, levels = c(LAB_R, LAB_NR, LAB_NS))]
  nR <- sum(d$hit & d$coef > 0); nNR <- sum(d$hit & d$coef < 0)
  top <- head(d[hit == TRUE][order(qval)], label_n)
  if (nrow(top)) top[, lab := gsub("_", " ", feature)]
  p <- ggplot(d, aes(coef, neglogq)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_hline(yintercept = -log10(SIG_Q), linetype = "dashed", color = "grey60") +
    geom_point(aes(color = ici_grp), size = point_size, alpha = 0.85) +
    scale_color_manual(values = VOLCOLS, name = "ICI response", drop = FALSE) +
    labs(x = "MaAsLin2 coefficient   (<-- higher in NR    |    higher in R -->)",
         y = expression(-log[10]~"(q-value, BH)"),
         title = title, subtitle = subtitle) +
    theme_classic(base_size = base_size) +
    theme(plot.title = element_text(face = "bold", size = base_size),
          plot.subtitle = element_text(size = 11),   # smaller so the long subtitle fits (matches lollipop)
          legend.position = "right")
  if (nrow(top)) p <- p + ggrepel::geom_text_repel(data = top, aes(coef, neglogq, label = lab),
        size = 3.9, fontface = "italic", max.overlaps = 20, min.segment.length = 0, segment.color = "grey70")
  # Render with the default ragg device (which renders the unicode glyphs -- arrows,
  # >= -- correctly) to an ASCII TEMP path, because ragg cannot open the non-ASCII
  # (Korean) destination path on macOS; then copy the finished PNG to the final
  # location with base R file.copy(), which handles UTF-8 paths.
  tmp <- tempfile(fileext = ".png")
  ggsave(tmp, p, width = width, height = height, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out_png, overwrite = TRUE)); unlink(tmp)
  cat(sprintf("  %-50s q<%.2f hits: R=%d  NR=%d  (n.s.=%d)\n",
              basename(out_png), SIG_Q, nR, nNR, nrow(d) - nR - nNR))
}

## ---- per-cohort (4): read the verified DA_species_group_<cohort>.csv ----
cohorts <- c(PRJNA751792 = "NSCLC_PRJNA751792", PRJNA1023797 = "NSCLC_PRJNA1023797",
             PRJEB22863  = "NSCLC_RCC_PRJEB22863", PRJEB26531  = "NSCLC_PRJEB26531")
COVARS  <- list(PRJNA751792 = "age+sex+BMI", PRJNA1023797 = "age+BMI",
                PRJEB22863 = "age+sex+antibiotic", PRJEB26531 = "age")   # mirrors each run_analysis.R
cat("Per-cohort species DA volcanoes (colored by ICI response group):\n")
for (cn in names(cohorts)) {
  f <- file.path(main_dir, cohorts[cn], "results", paste0("DA_species_group_", cn, ".csv"))
  stopifnot(file.exists(f))
  d <- fread(f)                                            # columns include feature, coef, qval
  stopifnot(all(c("feature", "coef", "qval") %in% names(d)))
  make_volcano(d,
    title = sprintf("NSCLC %s: differential species by ICI response", cn),
    subtitle = sprintf("MaAsLin2 (TSS+LOG+LM), adjusted for %s; colored by enrichment group (q < %.2f, BH); %d species",
                       COVARS[[cn]], SIG_Q, nrow(d)),
    out_png = file.path(main_dir, cohorts[cn], "results", paste0("DA_volcano_", cn, ".png")))
}

## ---- pooled meta: read the verified pooled_species_meta.csv ----
cat("Pooled meta volcano:\n")
fp <- file.path(main_dir, "pooled_analysis", "results", "pooled_species_meta.csv")
stopifnot(file.exists(fp))
dp <- fread(fp)                                            # columns: species, pooled_coef, qval, ...
stopifnot(all(c("species", "pooled_coef", "qval") %in% names(dp)))
setnames(dp, c("species", "pooled_coef"), c("feature", "coef"))    # unify to feature/coef/qval
make_volcano(dp,
  title = "Pooled species meta-analysis across 3 discovery cohorts (R vs NR)",
  subtitle = sprintf("Random-effects (DL) inverse-variance; colored by enrichment group (q < %.2f, BH); %d species in >=2 cohorts",
                     SIG_Q, nrow(dp)),
  out_png = file.path(main_dir, "pooled_analysis", "results", "pooled_species_meta_volcano.png"),
  point_size = 1.9, base_size = 16, label_n = 15, width = 11.5, height = 6)

writeLines(capture.output(sessionInfo()), file.path(main_dir, "sessionInfo_replot_volcano.txt"))
cat("DONE.\n")
