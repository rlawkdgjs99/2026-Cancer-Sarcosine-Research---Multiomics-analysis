#!/usr/bin/env Rscript
# ============================================================================
# Sup Fig 13 (b)(c) -- SPECIES vs SARCOSINE-SCORE scatter grids (NSCLC, R vs NR)
#   Parallels the CRC Sup Fig 7 layout (per-species small-multiple scatter:
#   species relative abundance (x) vs the sarcosine DEGRADATION / PRODUCTION
#   per-sample score (y)), but in the ICI R-vs-NR context.
#
#   (b) DEGRADATION grid: top R-enriched degraders
#   (c) PRODUCTION  grid: all NR-enriched producers
#
# FAITHFULNESS / INTEGRITY (this script is DISPLAY-ONLY):
#   - It computes NO new statistic. It never calls cor(). The rho / 95% CI / q
#     printed on each panel are READ VERBATIM from the verified meta table
#     pooled_analysis/results/sarcosine_species_assoc_meta.csv, whose rho is a
#     Fisher-z random-effects (DL) META of the per-cohort Spearman correlations
#     over the 3 discovery cohorts (NOT a naive pooled rho). The KO-score
#     definitions are untouched: y = the existing per-sample degradation /
#     production column from sarcosine_scores_<cohort>.csv (= sum of constituent
#     KO TSS relative abundances), read as-is.
#   - The scatter only VISUALISES the raw per-sample (abundance, score) pairs
#     behind that meta rho. Multi-cohort structure is shown (point SHAPE = cohort,
#     per metagenomics multi-cohort rule); colour = ICI response (R / NR).
#   - NO regression line is drawn: an OLS line on zero-inflated abundance is
#     leverage-dominated and can visually contradict the rank-based Spearman rho,
#     so the rank statistic is reported as text and the points shown raw.
#
#   Species SELECTION (matches the panel-a double-dissociation universe = DA-sig):
#     (b) score=="degradation" & grp=="Enriched in Responder (R)"      & corr q<0.05
#     (c) score=="production"  & grp=="Enriched in Non-responder (NR)" & corr q<0.05
#   where grp already encodes DA q<0.05 AND DA direction (verified in the CSV).
#
# Reads (read-only): sarcosine_species_assoc_meta.csv (selection + stats),
#   <cohort>/results/sarcosine_scores_<cohort>.csv (group + per-sample score),
#   <cohort>/Bacteria_*.txt (species relative abundance). Writes only 2 PNGs.
# Run: Rscript sarcosine_species_scatter_grids.R
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

# Run with the R working directory set to the analysis-folder root; paths below
# are relative to it.
main_dir <- "."
res <- file.path(main_dir, "pooled_analysis", "results")

save_png <- function(p, out_png, width, height, dpi = 150) {   # ragg cannot write the Korean path on macOS
  tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = width, height = height, dpi = dpi, bg = "white")
  stopifnot(file.copy(tmp, out_png, overwrite = TRUE)); unlink(tmp)
}

DISC <- c(PRJNA751792 = "NSCLC_PRJNA751792", PRJNA1023797 = "NSCLC_PRJNA1023797", PRJEB22863 = "NSCLC_RCC_PRJEB22863")
PAL  <- c(R = "#1B7837", NR = "#B2182B")          # group scheme: R=green (favorable), NR=red (unfavorable)
SHP  <- c(PRJNA751792 = 16, PRJNA1023797 = 17, PRJEB22863 = 15)   # cohort point shapes
PREV_MIN <- 0.10
TOP_B <- 8                                          # degradation grid: top-N by corr q

read_species <- function(cdir) {   # identical logic to dualfilter_bar.R / run_analysis.R
  f <- list.files(cdir, "^Bacteria_.*\\.txt$", full.names = TRUE); stopifnot(length(f) == 1L)
  b <- fread(f, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
  stopifnot(ncol(b) == 3L); setnames(b, c("Taxa", "RunID", "Abundance"))
  sp <- b[grepl("|s__", Taxa, fixed = TRUE) & !grepl("|t__", Taxa, fixed = TRUE)]
  w  <- dcast(sp, RunID ~ Taxa, value.var = "Abundance", fill = 0, fun.aggregate = sum)
  m  <- as.matrix(w[, -1L, with = FALSE]); rownames(m) <- w$RunID
  colnames(m) <- make.unique(sub(".*\\|s__", "", colnames(m)))
  m
}

## ---- 1. species selection + stats, READ from the verified meta table --------
am <- fread(file.path(res, "sarcosine_species_assoc_meta.csv"))
selB <- head(am[score == "degradation" & grp == "Enriched in Responder (R)"      & qval < 0.05][order(qval)], TOP_B)
selC <-      am[score == "production"  & grp == "Enriched in Non-responder (NR)" & qval < 0.05][order(qval)]
cat(sprintf("panel (b) degradation / R-enriched: %d species (showing top %d)\n", nrow(am[score=="degradation" & grp=="Enriched in Responder (R)" & qval<0.05]), nrow(selB)))
cat(sprintf("panel (c) production  / NR-enriched: %d species (showing all)\n", nrow(selC)))
stopifnot(nrow(selB) == TOP_B, nrow(selC) >= 1)

## ---- 2. assemble raw per-sample (abundance, score) over the 3 discovery cohorts
build_long <- function(score_name, species_vec) {
  out <- list()
  for (cn in names(DISC)) {
    cdir <- file.path(main_dir, DISC[cn])
    sc <- fread(file.path(cdir, "results", paste0("sarcosine_scores_", cn, ".csv")))
    spm <- read_species(cdir)
    stopifnot(all(sc$RunID %in% rownames(spm)))         # ID alignment (data-inspection rule)
    spm <- spm[sc$RunID, , drop = FALSE]                # align species rows to score samples
    for (sp in species_vec) {
      ab <- if (sp %in% colnames(spm)) spm[, sp] else rep(0, nrow(spm))   # absent in cohort -> 0
      out[[paste(cn, sp)]] <- data.table(cohort = cn, species = sp, group = sc$group,
                                         abund = as.numeric(ab), score = as.numeric(sc[[score_name]]))
    }
  }
  rbindlist(out)
}

make_grid <- function(score_name, sel, score_lab, ttl, sub, out_png, ncol, width, height) {
  long <- build_long(score_name, sel$species)
  ord  <- gsub("_", " ", sel$species)                  # facet order = selection order (by corr q)
  long[, species_lab := factor(gsub("_", " ", species), levels = ord)]
  long[, cohort := factor(cohort, levels = names(DISC))]
  ann <- sel[, .(species_lab = factor(gsub("_", " ", species), levels = ord),
                 lab = sprintf("rho[meta]=%.2f [%.2f, %.2f]\nq=%.1e, k=%d cohorts",
                               pooled_rho, ci_lb, ci_ub, qval, k))]
  p <- ggplot(long, aes(abund, score)) +
    geom_point(aes(colour = group, shape = cohort), size = 1.1, alpha = 0.55) +
    geom_label(data = ann, aes(x = -Inf, y = Inf, label = lab), hjust = -0.03, vjust = 1.1,
               size = 3.1, colour = "grey15", fill = "white", alpha = 0.7,
               linewidth = 0, label.padding = unit(0.6, "mm"), lineheight = 0.9, inherit.aes = FALSE) +
    facet_wrap(~ species_lab, scales = "free", ncol = ncol) +
    scale_colour_manual(values = PAL, name = "ICI response", breaks = c("R", "NR")) +
    scale_shape_manual(values = SHP, name = "Discovery cohort") +
    labs(x = "Species relative abundance (MetaPhlAn4)", y = score_lab, title = ttl, subtitle = sub,
         caption = paste0("Each point = one sample (3 discovery cohorts; shape = cohort, colour = ICI response). ",
                          "rho[meta], 95% CI and q are read verbatim from the random-effects (DL) Fisher-z meta of\n",
                          "per-cohort Spearman correlations (sarcosine_species_assoc_meta.csv); no statistic is recomputed here. ",
                          "No fit line is drawn (rank statistic, zero-inflated x).")) +
    guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1)),
           shape  = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
    theme_bw(base_size = 16) +
    theme(strip.text = element_text(face = "italic", size = 11),
          strip.background = element_rect(fill = "grey92", colour = NA),
          legend.position = "top", legend.box = "horizontal",
          plot.title = element_text(face = "bold", size = 17),
          plot.subtitle = element_text(colour = "grey35", size = 13),
          plot.caption = element_text(size = 9, colour = "grey45", hjust = 0),
          panel.grid.minor = element_blank())
  save_png(p, out_png, width = width, height = height, dpi = 150)
  cat(sprintf("wrote %s  (%d species, %d samples)\n", basename(out_png), nrow(sel), nrow(long) / nrow(sel)))
  invisible(long)
}

## ---- 3. panel (b) degradation, (c) production --------------------------------
make_grid("degradation", selB,
          score_lab = "Sarcosine DEGRADATION score (sum of degradation-KO rel. abundance)",
          ttl = "NSCLC ICI: Responder-enriched bacteria track the sarcosine DEGRADATION score",
          sub = sprintf("Top %d R-enriched degraders (DA q<0.05 & higher in R; correlation q<0.05) | discovery cohorts", nrow(selB)),
          out_png = file.path(res, "sarcosine_species_scatter_degradation.png"),
          ncol = 4, width = 16, height = 9)

make_grid("production", selC,
          score_lab = "Sarcosine PRODUCTION score (KO rel. abundance)",
          ttl = "NR-enriched bacteria track the sarcosine PRODUCTION score",
          sub = sprintf("All %d NR-enriched producers (DA q<0.05 & higher in NR; correlation q<0.05) | discovery cohorts", nrow(selC)),
          out_png = file.path(res, "sarcosine_species_scatter_production.png"),
          ncol = 3, width = 11, height = 6)

writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_species_scatter_grids.txt"))
cat("DONE.\n")
