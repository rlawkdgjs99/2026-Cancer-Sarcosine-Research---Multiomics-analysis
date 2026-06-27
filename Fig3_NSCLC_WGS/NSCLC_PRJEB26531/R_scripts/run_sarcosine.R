#!/usr/bin/env Rscript
# ============================================================================
# Cohort PRJEB26531 — SARCOSINE functional analysis (ICI response, R vs NR)
#   External-validation cohort (Korea). Only the 25 WGS NSCLC runs are used
#   (16S runs & Health excluded by the Assay==WGS & Phenotype==NSCLC filter).
#   (1) Pathway scores (TWO): degradation (incl. K18897 sarcosine->betaine) / production.
#       (Bidirectional KOs K21833/K21834 act in BOTH directions -> excluded from the scores,
#        but kept in the individual-KO differential abundance below.)
#   (2) Individual sarcosine-KO differential abundance: MaAsLin2, adjusted.
# Sarcosine KO set & roles: verified from KEGG via KEGGREST (main-folder
#   sarcosine_KO_derivation.R / sarcosine_KO_verify.R).
# KO abundances are TSS % -> MaAsLin2 normalization='NONE' + transform='LOG'.
# Self-contained. Run: Rscript NSCLC_PRJEB26531/R_scripts/run_sarcosine.R
# ============================================================================
set.seed(42)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")   # system2()/grep below needs a UTF-8 locale to read the non-ASCII (Korean) project path; the default C locale errors
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(Maaslin2) })

## ------------------------------------------------------------------ config
COHORT      <- "PRJEB26531"
COVARS      <- c("age")                    # n=25 (small): age only (parsimony)
NSCLC_LABEL <- "Carcinoma, Non-Small-Cell Lung"
KO_TEST_PREV<- 0.10                        # estimability gate for INDIVIDUAL KO tests
PAL         <- c(NR = "#B2182B", R = "#1B7837")   # group scheme: R=green (favorable), NR=red (unfavorable) (DISPLAY colour only; factor levels/reference unchanged)

## sarcosine KO groups (from verified KEGG reactions; documented in derivation)
DEG   <- c("K00301","K00302","K00303","K00304","K00305","K00306","K00314","K18897") # sarcosine CONSUMED (->glycine demethylation; K18897 ->betaine methylation)
PROD  <- c("K08688","K08687","K00315")                                     # -> sarcosine (creatine/carbamoylsarc/DMG)
BIDIR <- c("K21833","K21834")                                              # DMG <-> sarcosine <-> glycine
ALLKO <- unique(c(DEG, PROD, BIDIR))

## ----------------------------------------------------- relative paths
## Run with the R working directory set to the analysis-folder root; paths below
## are relative to it.
cohort_dir  <- file.path(".", "NSCLC_PRJEB26531")
main_dir    <- "."
results_dir <- file.path(cohort_dir, "results"); dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
logcon <- file(file.path(results_dir, paste0("sarcosine_log_", COHORT, ".txt")), open = "wt")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logcon) }
say("==== Sarcosine | %s | %s ====", COHORT, format(Sys.time(), tz = "UTC", usetz = TRUE))

## sanity: groups are a subset of the KEGG-verified KO set
koset <- fread(file.path(main_dir, "sarcosine_KO_set.csv"))
stopifnot(all(ALLKO %in% koset$KO))

## --------------------------------------------------------------- metadata --
parse_one_resp <- function(d) {
  if (is.na(d)) return(NA_character_)
  m <- regmatches(d, regexpr("response_group:[^;]*", d))
  if (length(m) == 0L || !nzchar(m)) m <- regmatches(d, regexpr("response:[^;]*", d))
  if (length(m) == 0L) return(NA_character_)
  v <- trimws(sub("^response(_group)?:[[:space:]]*", "", m)); if (v %in% c("R","NR")) v else NA_character_
}
parse_response <- function(x) vapply(as.character(x), parse_one_resp, character(1), USE.NAMES = FALSE)
read_meta <- function(cdir) {
  f <- list.files(cdir, "^selected_project_.*\\.txt$", full.names = TRUE); stopifnot(length(f) == 1L)
  ln <- sub("\r$", "", readLines(f, warn = FALSE)); hdr <- strsplit(ln[2], "\t", fixed = TRUE)[[1]]
  dat <- ln[-(1:2)]; dat <- dat[nzchar(dat)]
  mat <- do.call(rbind, lapply(strsplit(dat, "\t", fixed = TRUE), function(x) x[seq_len(length(hdr))]))
  d <- as.data.table(mat); setnames(d, hdr); d
}
meta <- read_meta(cohort_dir)
meta <- meta[`Assay type` == "WGS" & `Phenotype name` == NSCLC_LABEL]
meta[, group := parse_response(`Sample description`)]
meta <- meta[!is.na(group)]
meta[, group := factor(group, levels = c("NR","R"))]
meta[, age := suppressWarnings(as.numeric(age))]
meta[, BMI := suppressWarnings(as.numeric(BMI))]
meta[, sex := ifelse(sex %in% c("male","female"), sex, NA_character_)]
runs <- meta[["Run ID"]]
say("NSCLC-WGS R/NR samples: %d (R=%d, NR=%d)", length(runs), sum(meta$group=="R"), sum(meta$group=="NR"))

## ------------------------------- extract sarcosine-KO records from KO JSON --
# KO_relative_abundance.tsv is a single-line JSON array; pull only sarcosine KOs.
kof <- file.path(cohort_dir, "KO_relative_abundance.tsv")
pat <- sprintf('"ko": "(%s)", "run_id": "[^"]+", "abundance": [0-9.eE+-]+', paste(ALLKO, collapse = "|"))
# pass pattern via file (-f) so no shell parses the embedded double-quotes
patfile <- tempfile(); writeLines(pat, patfile)
recs <- system2("grep", c("-oE", "-f", patfile, kof), stdout = TRUE); unlink(patfile)
say("sarcosine-KO records extracted: %d", length(recs))
stopifnot(length(recs) > 0)   # fail loudly rather than emit empty/NaN results
ko  <- sub('"ko": "([^"]+)".*', "\\1", recs)
rid <- sub('.*"run_id": "([^"]+)".*', "\\1", recs)
ab  <- as.numeric(sub('.*"abundance": ([0-9.eE+-]+)$', "\\1", recs))
kodt <- data.table(KO = ko, RunID = rid, ab = ab)[RunID %in% runs]

## build runs x KO matrix (TSS %), fill absent with 0
M <- matrix(0, length(runs), length(ALLKO), dimnames = list(runs, ALLKO))
if (nrow(kodt)) M[cbind(match(kodt$RunID, runs), match(kodt$KO, ALLKO))] <- kodt$ab
prev <- colMeans(M > 0)
say("KO prevalence (%% of samples): %s",
    paste(names(prev), sprintf("%.0f", 100*prev), sep = "=", collapse = ", "))

## ----------------------------------------------------- 1. PATHWAY SCORES ----
grp_score <- function(kos) rowSums(M[, intersect(kos, colnames(M)), drop = FALSE])
sc <- data.table(RunID = runs, group = meta$group,
                 degradation = grp_score(DEG), production = grp_score(PROD))
# Production:Degradation balance = log2(production / degradation); defined ONLY where BOTH > 0
# (NA otherwise -- the ~2-6%% of samples with production == 0 are dropped from this metric).
sc[, prod_deg_log2ratio := ifelse(production > 0 & degradation > 0, log2(production / degradation), NA_real_)]
fwrite(sc, file.path(results_dir, paste0("sarcosine_scores_", COHORT, ".csv")))
score_names <- c("degradation","production","prod_deg_log2ratio")
score_lab   <- c(degradation = "Degradation\n(summed KO %)", production = "Production\n(summed KO %)",
                 prod_deg_log2ratio = "log2(Production / Degradation)\n(both > 0 only)")
sstat <- rbindlist(lapply(score_names, function(s) {
  wt <- wilcox.test(sc[[s]] ~ sc$group)
  data.table(score = s, median_NR = median(sc[group=="NR"][[s]], na.rm=TRUE), median_R = median(sc[group=="R"][[s]], na.rm=TRUE),
             W = unname(wt$statistic), p_wilcox = wt$p.value)
}))
fwrite(sstat, file.path(results_dir, paste0("sarcosine_score_stats_", COHORT, ".csv")))
say("SCORE Wilcoxon p: %s", paste(sstat$score, signif(sstat$p_wilcox,3), sep="=", collapse=", "))

fmt_p <- function(p) ifelse(p < 1e-3, sprintf("P=%.1e", p), sprintf("P=%.3f", p))   # (pre-existing; kept)
sl <- melt(sc, id.vars = c("RunID","group"), measure.vars = score_names,
           variable.name = "score", value.name = "value")
sl[, score := factor(score, levels = score_names)]
## y-axis ZOOM (scores have tiny units; a few extreme outliers compress every box): draw each box
## from the FULL-data five-number summary (boxplot.stats -> Tukey whiskers exclude outliers) via
## geom_boxplot(stat="identity"), and show jitter only within the panel's whisker range. INTEGRITY:
## medians / quartiles / whiskers and the Wilcoxon p use ALL samples; only the out-of-whisker POINTS
## are not drawn (count logged). NR stays the model reference (factor unchanged); scale_x_discrete
## sets DISPLAY order only (R left, NR right).
bx <- sl[is.finite(value), as.list(setNames(boxplot.stats(value)$stats,
          c("ymin","lower","middle","upper","ymax"))), by = .(score, group)]
rng <- bx[, .(lo = min(ymin), hi = max(ymax)), by = score]
jit <- merge(sl[is.finite(value)], rng, by = "score"); jit <- jit[value >= lo & value <= hi]
ncl <- merge(sl[is.finite(value), .(n=.N), by=score], jit[, .(shown=.N), by=score], by="score")[, .(score, clipped = n - shown)]
say("y-zoom clipped points (beyond whiskers; NOT dropped from any test): %s", paste(ncl$score, ncl$clipped, sep="=", collapse=", "))
pos <- merge(sstat[, .(score, p_wilcox)], rng, by = "score")
pos[, lab := sprintf("p = %.2g", p_wilcox)][, y := lo + 0.96*(hi - lo)]
for (D in list(bx, jit, pos, rng)) D[, score := factor(score, levels = score_names)]
p_sc <- ggplot(bx, aes(group, fill = group)) +
  geom_boxplot(aes(ymin = ymin, lower = lower, middle = middle, upper = upper, ymax = ymax),
               stat = "identity", width = 0.55, alpha = 0.65, linewidth = 0.4) +
  geom_jitter(data = jit, aes(group, value), inherit.aes = FALSE, width = 0.16, size = 0.9, alpha = 0.35, color = "grey25") +
  geom_text(data = pos, aes(x = 1.5, y = y, label = lab), inherit.aes = FALSE, size = 4.2, fontface = "bold") +
  facet_wrap(~ score, scales = "free_y", nrow = 1, labeller = labeller(score = score_lab)) +
  scale_fill_manual(values = PAL) +
  scale_x_discrete(limits = c("R","NR")) +
  labs(x = NULL, y = NULL,
       title = sprintf("NSCLC %s: gut sarcosine metabolism (R vs NR, n=%d)", COHORT, length(runs)),
       caption = "y-axis zoomed to each score's boxplot whiskers; boxes/medians and the Wilcoxon p use all samples (outlier points beyond whiskers not drawn).") +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", strip.background = element_blank(),
        plot.caption = element_text(size = 7.5, hjust = 0, color = "grey35"),
        strip.text = element_text(face = "bold", size = 12), plot.title = element_text(face = "bold", size = 16))
.tmp_png <- tempfile(fileext = ".png")    # ragg cannot write the non-ASCII (Korean) path on macOS
ggsave(.tmp_png, p_sc, width = 9.5, height = 4.7, dpi = 300, bg = "white")
stopifnot(file.copy(.tmp_png, file.path(results_dir, paste0("sarcosine_scores_", COHORT, ".png")), overwrite = TRUE)); unlink(.tmp_png)

## ----------------------------------------- 2. INDIVIDUAL KO DA (MaAsLin2) ---
testable <- names(prev)[prev >= KO_TEST_PREV]
say("individual-KO testable (prev>=%.0f%%): %d -> %s", 100*KO_TEST_PREV, length(testable), paste(testable, collapse=", "))
da_done <- FALSE
if (length(testable) >= 1) {
  md <- as.data.frame(meta[, c("Run ID","group", COVARS), with = FALSE])
  rownames(md) <- md[["Run ID"]]; md[["Run ID"]] <- NULL
  cc <- stats::complete.cases(md[, c("group", COVARS)])
  md <- md[cc, , drop = FALSE]; md$group <- droplevels(md$group)
  say("DA complete-case: %d (dropped %d)", sum(cc), sum(!cc))
  dai <- as.data.frame(M[rownames(md), testable, drop = FALSE])
  da_out <- file.path(results_dir, paste0("maaslin2_sarcosine_", COHORT))
  invisible(Maaslin2(input_data = dai, input_metadata = md, output = da_out,
    fixed_effects = c("group", COVARS), min_prevalence = 0, min_abundance = 0,
    normalization = "NONE", transform = "LOG", analysis_method = "LM",
    max_significance = 1, correction = "BH", standardize = TRUE,
    plot_heatmap = FALSE, plot_scatter = FALSE, cores = 1))
  ar <- fread(file.path(da_out, "all_results.tsv"))
  grp <- ar[metadata == "group"]
  setorder(grp, qval)
  grp[, direction := ifelse(coef > 0, "higher in R", "higher in NR")]
  grp <- merge(grp, data.table(feature = names(prev), prevalence = round(100*prev,1)), by = "feature", all.x = TRUE)
  fwrite(grp, file.path(results_dir, paste0("sarcosine_KO_DA_", COHORT, ".csv")))
  say("KO DA: tested %d | q<0.10: %d | q<0.05: %d", nrow(grp), sum(grp$qval<0.10,na.rm=TRUE), sum(grp$qval<0.05,na.rm=TRUE))
  # forest of coefficients (coef +/- 95% CI), colored by q
  grp[, lo := coef - 1.96*stderr][, hi := coef + 1.96*stderr]
  grp[, sig := ifelse(qval < 0.05, "q<0.05", ifelse(qval < 0.10, "q<0.10", "n.s."))]
  grp[, KOlab := feature]
  grp[, KOlab := factor(KOlab, levels = grp[order(coef)]$feature)]
  p_da <- ggplot(grp, aes(coef, KOlab, color = sig)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_segment(aes(x = lo, xend = hi, y = KOlab, yend = KOlab), linewidth = 0.5) +
    geom_point(size = 2.4) +
    scale_color_manual(values = c("q<0.05"="#D55E00","q<0.10"="#E69F00","n.s."="grey55"), name = NULL) +
    labs(x = "MaAsLin2 coefficient  (← NR    |    R →)", y = NULL,
         title = sprintf("NSCLC %s: individual sarcosine-KO differential abundance (R vs NR)", COHORT),
         subtitle = sprintf("MaAsLin2 NONE+LOG+LM, adjusted for %s; bars = 95%% CI", paste(COVARS, collapse="+"))) +
    theme_classic(base_size = 12) + theme(plot.title = element_text(face = "bold", size = 12))
  ggsave(file.path(results_dir, paste0("sarcosine_KO_DA_", COHORT, ".png")), p_da,
         width = 7.5, height = 0.4*nrow(grp)+2, dpi = 300, bg = "white")
  da_done <- TRUE
}
if (!da_done) say("No KO passed the estimability gate; individual-KO DA skipped.")

writeLines(capture.output(sessionInfo()), file.path(results_dir, paste0("sessionInfo_sarcosine_", COHORT, ".txt")))
say("DONE. outputs in %s", results_dir)
close(logcon)
