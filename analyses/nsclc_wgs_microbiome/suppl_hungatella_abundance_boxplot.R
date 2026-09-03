#!/usr/bin/env Rscript
# ============================================================================
# Hungatella hathewayi gut relative-abundance box plot by ICI response (R vs NR).
# Saved PER COHORT into each cohort's own results/ folder, AND a pooled-discovery
# version into pooled_analysis/results/ -- matching the project's per-cohort +
# pooled output layout (one figure per analysis folder, NOT a single combined plot).
#   y = log10(relative abundance % + pseudocount) (zero-inflation / right-skew);
#   each plot annotated with that unit's covariate-adjusted MaAsLin2 DA
#   (direction + FDR q; meta for the pooled). Same y-axis across plots (comparable).
# READ-ONLY on verified inputs; recomputes only per-sample abundance for display.
# Run: Rscript suppl_hungatella_abundance_boxplot.R
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
# Run with the R working directory set to the analysis-folder root; paths below are relative to it.
main_dir <- "."; res <- file.path(main_dir, "pooled_analysis", "results")
save_png <- function(p, out, w, h) { tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = w, height = h, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out, overwrite = TRUE)); unlink(tmp) }
read_species <- function(cdir) {
  f <- list.files(cdir, "^Bacteria_.*\\.txt$", full.names = TRUE); stopifnot(length(f) == 1L)
  b <- fread(f, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
  stopifnot(ncol(b) == 3L); setnames(b, c("Taxa", "RunID", "Abundance"))
  sp <- b[grepl("|s__", Taxa, fixed = TRUE) & !grepl("|t__", Taxa, fixed = TRUE)]
  w  <- dcast(sp, RunID ~ Taxa, value.var = "Abundance", fill = 0, fun.aggregate = sum)
  m  <- as.matrix(w[, -1L, with = FALSE]); rownames(m) <- w$RunID
  colnames(m) <- make.unique(sub(".*\\|s__", "", colnames(m))); m
}

SP   <- "Hungatella_hathewayi"
COH  <- c(PRJNA751792 = "NSCLC_PRJNA751792", PRJNA1023797 = "NSCLC_PRJNA1023797",
          PRJEB22863  = "NSCLC_RCC_PRJEB22863", PRJEB26531 = "NSCLC_PRJEB26531")
DISC <- c("PRJNA751792", "PRJNA1023797", "PRJEB22863")
PAL  <- c(R = "#1B7837", NR = "#B2182B")

## ---- per-cohort abundance + group ------------------------------------------
rows <- list()
for (cn in names(COH)) {
  cdir <- file.path(main_dir, COH[cn])
  sc <- fread(file.path(cdir, "results", paste0("sarcosine_scores_", cn, ".csv")))   # RunID, group
  spm <- read_species(cdir); stopifnot(all(sc$RunID %in% rownames(spm)))
  ab <- if (SP %in% colnames(spm)) spm[sc$RunID, SP] else rep(0, nrow(sc))
  rows[[cn]] <- data.table(cohort = cn, group = as.character(sc$group), ab = as.numeric(ab))
}
A <- rbindlist(rows)
A[, group := factor(group, levels = c("NR", "R"))]
pseudo <- min(A$ab[A$ab > 0]) / 2
A[, logab := log10(ab + pseudo)]
YR <- range(A$logab); YPAD <- diff(YR)                          # global y-range -> comparable axes

## ---- DA annotation (per-cohort MaAsLin2; meta for pooled) -------------------
da <- rbindlist(lapply(names(COH), function(cn) {
  d <- fread(file.path(main_dir, COH[cn], "results", paste0("DA_species_group_", cn, ".csv")))[feature == SP]
  data.table(cohort = cn, dir = ifelse(d$coef > 0, "R", "NR"), q = d$qval) }))
mt <- fread(file.path(res, "pooled_species_meta.csv"))[species == SP]
da <- rbind(da, data.table(cohort = "pooled", dir = ifelse(mt$pooled_coef > 0, "R", "NR"), q = mt$qval))
da[, lab := sprintf("MaAsLin2: %s-dir, q=%.2f", dir, q)]

## ---- one box plot per unit -------------------------------------------------
make_box <- function(d, ann_lab, ttl, out_png) {
  p <- ggplot(d, aes(group, logab, fill = group)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.55, linewidth = 0.4) +
    geom_jitter(width = 0.15, height = 0, size = 0.9, alpha = 0.45, color = "grey25") +
    annotate("text", x = 1.5, y = YR[2] + 0.07 * YPAD, label = ann_lab, size = 4.6, fontface = "bold") +
    scale_fill_manual(values = PAL, guide = "none") +
    scale_x_discrete(limits = c("R", "NR")) +
    coord_cartesian(ylim = c(YR[1], YR[2] + 0.16 * YPAD)) +
    labs(x = "ICI response group", y = expression(log[10]~"(relative abundance %, +pseudo)"),
         title = ttl,
         caption = "Box = log10 relative abundance; annotation = MaAsLin2 DA (direction, FDR q).") +
    theme_bw(base_size = 16) +
    theme(plot.title = element_text(size = 13), plot.caption = element_text(size = 9, color = "grey45", hjust = 0))
  save_png(p, out_png, 5.2, 4.4)
}

## (i) per cohort -> each cohort's results/
for (cn in names(COH)) {
  n <- A[cohort == cn, .N]
  make_box(A[cohort == cn], da[cohort == cn]$lab,
           bquote(italic("Hungatella hathewayi")~.(sprintf(" abundance: %s (n=%d)", cn, n))),
           file.path(main_dir, COH[cn], "results", sprintf("Hungatella_abundance_boxplot_%s.png", cn)))
  cat(sprintf("  wrote %s/results/Hungatella_abundance_boxplot_%s.png\n", COH[cn], cn))
}
## (ii) pooled discovery (3 cohorts) -> pooled_analysis/results/
np <- A[cohort %in% DISC, .N]
make_box(A[cohort %in% DISC], da[cohort == "pooled"]$lab,
         bquote(italic("Hungatella hathewayi")~.(sprintf(" (pooled discovery, n=%d)", np))),
         file.path(res, "Hungatella_abundance_boxplot_pooled.png"))
cat("  wrote pooled_analysis/results/Hungatella_abundance_boxplot_pooled.png\n")

writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_suppl_hungatella_abund.txt"))
cat("pseudocount =", signif(pseudo, 3), "\n")
