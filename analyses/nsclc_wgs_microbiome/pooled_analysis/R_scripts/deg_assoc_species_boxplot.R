#!/usr/bin/env Rscript
# =============================================================================
# SUPPLEMENTARY boxplots: relative abundance of the 3 cross-cohort sarcosine-
# DEGRADATION-associated species (the 3-set discovery Venn intersection), R vs NR.
# NSCLC analog of the CRC supp_deg_assoc_5species_boxplot.R (Healthy vs Cancer).
#
# Figures (all -> pooled_analysis/results/):
#   - pooled DISCOVERY (3 cohorts)   deg_assoc_3sp_box_pooled_discovery_NSCLC.png
#   - pooled ALL 4 (incl. Korea)     deg_assoc_3sp_box_pooled_all4_NSCLC.png
#   - per cohort (x4)                deg_assoc_3sp_box_<cohort>_NSCLC.png
#   - stats                          deg_assoc_3sp_box_stats_NSCLC.csv
#
# INTEGRITY: these 3 species were selected by DEGRADATION-SCORE CORRELATION (the
#   Venn), NOT by R-vs-NR differential abundance. In the verified meta DA
#   (pooled_species_meta.csv) all 3 TREND higher in R (coef>0) but NONE is
#   FDR-significant (Lachnospira_eligens q=0.083, Clostridium_sp_AF36_4 q=0.147,
#   Roseburia_faecis q=0.461). Pooled panels are annotated with this meta DA
#   (cohort-aware, the primary R-vs-NR test); per-cohort panels with that cohort's
#   Wilcoxon. Boxplots are DESCRIPTIVE. Separate analysis from the CRC figure.
#
# y = log10 relative abundance (%); zeros drawn at a per-species floor =
#   min(nonzero)/5 (CRC convention). Wilcoxon on RAW abundance. seed 42; ragg-safe.
# Run: Rscript pooled_analysis/R_scripts/deg_assoc_species_boxplot.R
# =============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

PAL  <- c(R = "#1B7837", NR = "#B2182B")                 # group scheme: R=green (favorable), NR=red (unfavorable)
COHORT_PAL <- c(PRJNA751792 = "#332288", PRJNA1023797 = "#CC6677",
                PRJEB22863 = "#44AA99", PRJEB26531 = "#DDCC77")   # matches the Venn
DISC  <- c(PRJNA751792 = "NSCLC_PRJNA751792", PRJNA1023797 = "NSCLC_PRJNA1023797", PRJEB22863 = "NSCLC_RCC_PRJEB22863")
VALID <- c(PRJEB26531 = "NSCLC_PRJEB26531")
ALLC  <- c(DISC, VALID)

# Run with the R working directory set to the analysis-folder root; paths below
# are relative to it.
pooled_dir <- file.path(".", "pooled_analysis"); main_dir <- "."
res <- file.path(pooled_dir, "results")
logcon <- file(file.path(res, "deg_assoc_3sp_box_log.txt"), open = "wt")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logcon) }
save_png <- function(p, out, w, h) { tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = w, height = h, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out, overwrite = TRUE)); unlink(tmp) }
short <- function(x) sub(".*\\|s__", "", x)

read_species <- function(cdir) {   # full taxa names kept (match the intersection CSV)
  f <- list.files(cdir, "^Bacteria_.*\\.txt$", full.names = TRUE); stopifnot(length(f) == 1L)
  b <- fread(f, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
  stopifnot(ncol(b) == 3L); setnames(b, c("Taxa", "RunID", "Abundance"))
  sp <- b[grepl("|s__", Taxa, fixed = TRUE) & !grepl("|t__", Taxa, fixed = TRUE)]
  w  <- dcast(sp, RunID ~ Taxa, value.var = "Abundance", fill = 0, fun.aggregate = sum)
  m  <- as.matrix(w[, -1L, with = FALSE]); rownames(m) <- w$RunID; m }

## ---- targets: the 3 intersection species, display order by mean discovery rho
inter <- fread(file.path(res, "cross_cohort_deg_assoc_intersection_discovery_NSCLC.csv"))
setorder(inter, -mean_disc_rho)
TARGET <- inter$Taxa; SP_DISP <- gsub("_", " ", short(TARGET)); stopifnot(length(TARGET) == 3)
say("targets (3, by mean discovery rho): %s", paste(SP_DISP, collapse = ", "))

## ---- long table: per-sample abundance of the 3 targets + group + cohort ------
build <- rbindlist(lapply(names(ALLC), function(cn) {
  cdir <- file.path(main_dir, ALLC[cn])
  sc <- fread(file.path(cdir, "results", paste0("sarcosine_scores_", cn, ".csv")))   # RunID, group
  spm <- read_species(cdir); stopifnot(all(sc$RunID %in% rownames(spm)))
  spm <- spm[sc$RunID, , drop = FALSE]
  ab <- sapply(TARGET, function(tx) if (tx %in% colnames(spm)) spm[, tx] else rep(0, nrow(spm)))
  dt <- as.data.table(ab); setnames(dt, SP_DISP)
  dt[, `:=`(RunID = sc$RunID, cohort = cn, group = as.character(sc$group))]
  melt(dt, id.vars = c("RunID", "cohort", "group"), variable.name = "Species", value.name = "Abundance")
}))
build[, group := factor(group, levels = c("NR", "R"))]
build[, Species := factor(Species, levels = SP_DISP)]
build[, cohort := factor(cohort, levels = names(ALLC))]

## ---- verified meta DA (primary R-vs-NR test) for the pooled-panel annotation -
da <- fread(file.path(res, "pooled_species_meta.csv"))[, .(species, dir, qval)]
metaDA <- merge(data.table(Species = SP_DISP, sp_short = short(TARGET)), da,
                by.x = "sp_short", by.y = "species", all.x = TRUE)
metaDA[, Species := factor(Species, levels = SP_DISP)]
metaDA[, label := sprintf("meta DA: %s\nq=%.3f%s", sub("higher in ", "", dir), qval, ifelse(qval < 0.05, "", " (n.s.)"))]

## ---- Wilcoxon stats (raw abundance) per scope --------------------------------
wstat <- function(dt, scope) dt[, .(scope = scope,
    n_R = sum(group == "R"), n_NR = sum(group == "NR"),
    prev_R = round(100*mean(Abundance[group == "R"] > 0), 1), prev_NR = round(100*mean(Abundance[group == "NR"] > 0), 1),
    median_R = median(Abundance[group == "R"]), median_NR = median(Abundance[group == "NR"]),
    wilcox_p = tryCatch(wilcox.test(Abundance ~ group)$p.value, error = function(e) NA_real_)), by = Species]
stats <- rbindlist(c(
  lapply(names(ALLC), function(cn) wstat(build[cohort == cn], cn)),
  list(wstat(build[cohort %in% names(DISC)], "pooled_discovery"), wstat(build, "pooled_all4"))))
stats <- merge(stats, metaDA[, .(Species, meta_dir = dir, meta_q = qval)], by = "Species")
fwrite(stats, file.path(res, "deg_assoc_3sp_box_stats_NSCLC.csv"))

## ---- plot helper -------------------------------------------------------------
plot_box <- function(dt, ann, title, subtitle, out, color_by_cohort) {
  d <- copy(dt)
  d[, floor_k := min(Abundance[Abundance > 0]) / 5, by = Species]
  d[, y := ifelse(Abundance <= 0, floor_k, Abundance)]
  a <- merge(ann, d[, .(yp = max(y) * 1.6), by = Species], by = "Species")
  p <- ggplot(d, aes(group, y)) +
    geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.65, width = 0.6, linewidth = 0.4)
  if (color_by_cohort) p <- p +
    geom_jitter(aes(color = cohort), width = 0.15, height = 0, size = 0.7, alpha = 0.5) +
    scale_color_manual(values = COHORT_PAL, name = "Cohort")
  else p <- p + geom_jitter(width = 0.15, height = 0, size = 0.7, alpha = 0.30, color = "grey30")
  p <- p +
    geom_text(data = a, aes(x = 1.5, y = yp, label = label), inherit.aes = FALSE, size = 4.7, fontface = "italic", vjust = 1, lineheight = 0.95) +
    facet_wrap(~ Species, scales = "free_y", nrow = 1) +
    scale_y_log10(expand = expansion(mult = c(0.05, 0.20))) + annotation_logticks(sides = "l") +
    scale_fill_manual(values = PAL, name = "ICI response") +
    scale_x_discrete(limits = c("R", "NR")) +
    labs(title = title, subtitle = subtitle, x = NULL, y = "Relative abundance (log10, %)",
         caption = "Species selected by degradation-score correlation (the Venn), NOT by R/NR DA; all 3 trend higher in R but are n.s. (FDR).") +
    theme_bw(base_size = 16) +
    theme(plot.title = element_text(face = "bold", size = 18),
          plot.subtitle = element_text(size = 12, color = "grey35"),
          plot.caption = element_text(size = 10, color = "grey45", hjust = 0),
          strip.text = element_text(face = "italic", size = 14),
          legend.position = "bottom")
  save_png(p, out, 11, 5.2)
  say("wrote %s", basename(out))
}

## ---- (1) pooled discovery + (2) pooled all 4 (meta-DA annotation) ------------
nd <- nrow(build[cohort %in% names(DISC)]) / 3; na <- nrow(build) / 3
plot_box(build[cohort %in% names(DISC)], metaDA,
  "Sarcosine degradation-associated species: abundance by ICI response",
  sprintf("Pooled 3 discovery cohorts (n=%d) | panel label = verified meta DA (cohort-aware, primary R-vs-NR test) | points colored by cohort", nd),
  file.path(res, "deg_assoc_3sp_box_pooled_discovery_NSCLC.png"), TRUE)
plot_box(build, metaDA,
  "Sarcosine degradation-associated species: abundance by ICI response",
  sprintf("All 4 cohorts incl. Korea (n=%d) | panel label = discovery meta DA | Korea (n=25) shown descriptively | points colored by cohort", na),
  file.path(res, "deg_assoc_3sp_box_pooled_all4_NSCLC.png"), TRUE)

## ---- (3) per-cohort (that cohort's Wilcoxon annotation) ----------------------
for (cn in names(ALLC)) {
  ann_cn <- stats[scope == cn, .(Species, label = sprintf("Wilcoxon\np=%s",
              ifelse(is.na(wilcox_p), "NA", formatC(wilcox_p, format = "g", digits = 2))))]
  ann_cn[, Species := factor(Species, levels = SP_DISP)]
  tag <- if (cn %in% names(VALID)) " [validation, n=25]" else ""
  plot_box(build[cohort == cn], ann_cn,
    sprintf("Degradation-associated species (R vs NR) - %s%s", cn, tag),
    sprintf("%s (n=%d) | per-panel Wilcoxon R vs NR (raw abundance, single cohort)", cn, nrow(build[cohort == cn]) / 3),
    file.path(res, sprintf("deg_assoc_3sp_box_%s_NSCLC.png", cn)), FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_deg_assoc_3sp_box.txt"))
say("DONE. outputs in %s", res)
close(logcon)
