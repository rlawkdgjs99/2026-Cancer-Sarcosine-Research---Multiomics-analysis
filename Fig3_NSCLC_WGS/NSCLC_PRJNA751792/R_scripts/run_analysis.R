#!/usr/bin/env Rscript
# ============================================================================
# Cohort PRJNA751792 (NSCLC, WGS, stool) — ICI response (R vs NR) microbiome
#   (1) Alpha diversity: Shannon / Simpson / Observed  + Wilcoxon rank-sum
#   (2) Beta  diversity: Bray-Curtis -> PCoA + PERMANOVA (adonis2, 999 perm)
#   (3) Differential abundance: MaAsLin2 (species level, adjusted)
# Data: MetaPhlAn4 species relative abundance (TSS %, sums to 100 per sample).
# Self-contained. Run:  Rscript NSCLC_PRJNA751792/R_scripts/run_analysis.R
# faithful-coder: inspect-before/after, seeds fixed, thresholds justified.
# ============================================================================

set.seed(42)
suppressPackageStartupMessages({
  library(data.table); library(vegan); library(ggplot2); library(ggrepel); library(Maaslin2)
})

## ------------------------------------------------------------------ config
COHORT       <- "PRJNA751792"
NSCLC_LABEL  <- "Carcinoma, Non-Small-Cell Lung"
COVARS       <- c("age", "sex", "BMI")     # confounders available in THIS cohort
PREV_MIN     <- 0.10                        # prevalence filter for DA (paper: 10%)
PAL          <- c(NR = "#B2182B", R = "#1B7837")  # group scheme: R=green (favorable), NR=red (unfavorable) (DISPLAY colour only; factor/reference unchanged)

## ----------------------------------------------------- relative paths
## Run with the R working directory set to the analysis-folder root; paths below
## are relative to it.
cohort_dir  <- file.path(".", "NSCLC_PRJNA751792")          # <cohort>
results_dir <- file.path(cohort_dir, "results")
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
logfile <- file.path(results_dir, paste0("analysis_log_", COHORT, ".txt"))
logcon  <- file(logfile, open = "wt")
say <- function(...) { msg <- sprintf(...); cat(msg, "\n"); cat(msg, "\n", file = logcon) }

say("==== Cohort %s | %s ====", COHORT, format(Sys.time(), tz = "UTC", usetz = TRUE))

## --------------------------------------------------------------- loaders ---
# metadata: line1 preamble, line2 header(23 cols), data rows have an unnamed
# trailing 24th field -> keep first 23 fields aligned to header.
parse_one_resp <- function(d) {
  if (is.na(d)) return(NA_character_)
  m <- regmatches(d, regexpr("response_group:[^;]*", d))         # prefer response_group:
  if (length(m) == 0L || !nzchar(m)) m <- regmatches(d, regexpr("response:[^;]*", d))
  if (length(m) == 0L) return(NA_character_)
  v <- trimws(sub("^response(_group)?:[[:space:]]*", "", m))
  if (v %in% c("R", "NR")) v else NA_character_
}
parse_response <- function(x) vapply(as.character(x), parse_one_resp, character(1), USE.NAMES = FALSE)

read_meta <- function(cdir) {
  f <- list.files(cdir, pattern = "^selected_project_.*\\.txt$", full.names = TRUE)
  stopifnot(length(f) == 1L)
  ln  <- sub("\r$", "", readLines(f, warn = FALSE))
  hdr <- strsplit(ln[2], "\t", fixed = TRUE)[[1]]
  dat <- ln[-(1:2)]; dat <- dat[nzchar(dat)]
  fl  <- strsplit(dat, "\t", fixed = TRUE)
  mat <- do.call(rbind, lapply(fl, function(x) x[seq_len(length(hdr))]))
  d   <- as.data.table(mat); setnames(d, hdr)
  d
}

read_species <- function(cdir) {
  f <- list.files(cdir, pattern = "^Bacteria_.*\\.txt$", full.names = TRUE)
  stopifnot(length(f) == 1L)
  b <- fread(f, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
  stopifnot(ncol(b) == 3L); setnames(b, c("Taxa", "RunID", "Abundance"))
  sp <- b[grepl("|s__", Taxa, fixed = TRUE) & !grepl("|t__", Taxa, fixed = TRUE)]  # species only
  w  <- dcast(sp, RunID ~ Taxa, value.var = "Abundance", fill = 0, fun.aggregate = sum)
  m  <- as.matrix(w[, -1L, with = FALSE]); rownames(m) <- w$RunID
  m
}

## ----------------------------------------------------------- 1. load+filter
meta <- read_meta(cohort_dir)
meta <- meta[`Assay type` == "WGS" & `Phenotype name` == NSCLC_LABEL]
meta[, group := parse_response(`Sample description`)]
meta <- meta[!is.na(group)]
meta[, group := factor(group, levels = c("NR", "R"))]      # NR = reference
meta[, age := suppressWarnings(as.numeric(age))]
meta[, BMI := suppressWarnings(as.numeric(BMI))]
meta[, sex := ifelse(sex %in% c("male", "female"), sex, NA_character_)]
rid <- meta[["Run ID"]]

spm <- read_species(cohort_dir)
stopifnot(all(rid %in% rownames(spm)))
spm <- spm[rid, , drop = FALSE]                            # align order to meta
# clean, unique species labels (drop k__..|s__ prefix)
colnames(spm) <- make.unique(sub(".*\\|s__", "", colnames(spm)))

say("samples (R/NR labeled NSCLC-WGS): %d  [R=%d, NR=%d]",
    nrow(meta), sum(meta$group == "R"), sum(meta$group == "NR"))
say("species detected: %d | per-sample TSS sum range: %.2f-%.2f",
    ncol(spm), min(rowSums(spm)), max(rowSums(spm)))
stopifnot(all(abs(rowSums(spm) - 100) < 1))                # sanity: TSS ~100
stopifnot(identical(rownames(spm), rid))                   # sanity: alignment

## --------------------------------------------------------- 2. ALPHA diversity
# Diversity on the FULL species community (not prevalence-filtered).
alpha <- data.table(
  RunID    = rownames(spm),
  group    = meta$group,
  Shannon  = diversity(spm, index = "shannon",  MARGIN = 1),
  Simpson  = diversity(spm, index = "simpson",  MARGIN = 1),
  Observed = rowSums(spm > 0)
)
idx_names <- c("Shannon", "Simpson", "Observed")
alpha_stats <- rbindlist(lapply(idx_names, function(ix) {
  wt <- wilcox.test(alpha[[ix]] ~ alpha$group)
  data.table(index = ix,
             median_NR = median(alpha[group == "NR"][[ix]]),
             median_R  = median(alpha[group == "R"][[ix]]),
             W = unname(wt$statistic), p_wilcox = wt$p.value)
}))
fwrite(alpha_stats, file.path(results_dir, paste0("alpha_diversity_stats_", COHORT, ".csv")))
say("ALPHA Wilcoxon p: %s",
    paste(alpha_stats$index, signif(alpha_stats$p_wilcox, 3), sep = "=", collapse = ", "))

fmt_p <- function(p) ifelse(p < 1e-3, sprintf("p = %.1e", p), sprintf("p = %.3f", p))   # plain p (test = Wilcoxon rank-sum; stated in figure legend)
al_long <- melt(alpha, id.vars = c("RunID", "group"),
                measure.vars = idx_names, variable.name = "index", value.name = "value")
al_long[, index := factor(index, levels = idx_names)]
pos <- al_long[, .(y = max(value) + 0.08 * (max(value) - min(value))), by = index]
pos <- merge(pos, alpha_stats[, .(index, p_wilcox)], by = "index")
pos[, lab := fmt_p(p_wilcox)]

p_alpha <- ggplot(al_long, aes(group, value, fill = group)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65, linewidth = 0.4) +
  geom_jitter(aes(color = group), width = 0.16, size = 1.2, alpha = 0.6, show.legend = FALSE) +
  geom_text(data = pos, aes(x = 1.5, y = y, label = lab),
            inherit.aes = FALSE, size = 4.0, fontface = "bold") +
  facet_wrap(~ index, scales = "free_y") +
  scale_fill_manual(values = PAL) + scale_color_manual(values = PAL) +
  scale_x_discrete(limits = c("R", "NR")) +     # display order only: R left, NR right (factor/reference unchanged)
  labs(x = NULL, y = "Alpha diversity",
       title = sprintf("NSCLC %s: gut alpha-diversity (R vs NR, n=%d)", COHORT, nrow(meta))) +
  theme_classic(base_size = 14) +
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 13),
        plot.title = element_text(face = "bold", size = 15))
.tmp_png <- tempfile(fileext = ".png")   # ragg cannot write the non-ASCII (Korean) path on macOS
ggsave(.tmp_png, p_alpha, width = 9, height = 3.6, dpi = 300, bg = "white")
stopifnot(file.copy(.tmp_png, file.path(results_dir, paste0("alpha_diversity_", COHORT, ".png")), overwrite = TRUE)); unlink(.tmp_png)

## --------------------------------------------------------- 3. BETA diversity
bc  <- vegdist(spm, method = "bray")
set.seed(42)
per <- adonis2(bc ~ group, data = as.data.frame(meta), permutations = 999, by = "terms")
R2  <- per$R2[1]; Fv <- per$F[1]; pp <- per$`Pr(>F)`[1]
permdt <- data.table(term = rownames(per), Df = per$Df, SumOfSqs = per$SumOfSqs,
                     R2 = per$R2, F = per$F, p = per$`Pr(>F)`)
fwrite(permdt, file.path(results_dir, paste0("beta_permanova_", COHORT, ".csv")))
say("BETA PERMANOVA: R2=%.4f, F=%.3f, p=%.4f", R2, Fv, pp)

pco <- cmdscale(bc, k = 2, eig = TRUE)
ev  <- pco$eig; ve <- 100 * ev[1:2] / sum(ev[ev > 0])
pcd <- data.table(PCo1 = pco$points[,1], PCo2 = pco$points[,2], group = meta$group)
sub_lab <- sprintf("PERMANOVA: R² = %.3f, P = %s  (Bray-Curtis, 999 perm)",
                   R2, ifelse(pp < 1e-3, sprintf("%.1e", pp), sprintf("%.3f", pp)))
p_beta <- ggplot(pcd, aes(PCo1, PCo2, color = group, fill = group)) +
  stat_ellipse(geom = "polygon", alpha = 0.10, color = NA, level = 0.95) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = PAL) + scale_fill_manual(values = PAL) +
  labs(x = sprintf("PCo1 (%.1f%%)", ve[1]), y = sprintf("PCo2 (%.1f%%)", ve[2]),
       color = "ICI response", fill = "ICI response",
       title = sprintf("NSCLC %s: gut beta-diversity (PCoA) by ICI response", COHORT),
       subtitle = sub_lab) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 13),
        legend.position = "right")
.tmp_png <- tempfile(fileext = ".png")   # ragg cannot write the non-ASCII (Korean) path on macOS
ggsave(.tmp_png, p_beta, width = 6.4, height = 5.2, dpi = 300, bg = "white")
stopifnot(file.copy(.tmp_png, file.path(results_dir, paste0("beta_pcoa_", COHORT, ".png")), overwrite = TRUE)); unlink(.tmp_png)

## ----------------------------------------------------- 4. DA via MaAsLin2 ---
# input: species relative abundance (samples x features); MaAsLin2 does TSS+LOG+LM.
# Adjusted for available confounders; complete-case on covariates.
md <- as.data.frame(meta[, c("Run ID", "group", COVARS), with = FALSE])
rownames(md) <- md[["Run ID"]]; md[["Run ID"]] <- NULL
cc <- stats::complete.cases(md[, c("group", COVARS)])
say("DA complete-case: %d / %d samples (dropped %d with missing covariates)",
    sum(cc), nrow(md), sum(!cc))
md  <- md[cc, , drop = FALSE]; md$group <- droplevels(md$group)
dai <- as.data.frame(spm[rownames(md), , drop = FALSE])

da_out <- file.path(results_dir, paste0("maaslin2_", COHORT))
invisible(Maaslin2(
  input_data       = dai,
  input_metadata   = md,
  output           = da_out,
  fixed_effects    = c("group", COVARS),
  min_prevalence   = PREV_MIN,
  min_abundance    = 0,
  normalization    = "TSS",
  transform        = "LOG",
  analysis_method  = "LM",
  max_significance = 0.25,
  correction       = "BH",
  standardize      = TRUE,
  plot_heatmap     = FALSE,
  plot_scatter     = FALSE,
  cores            = 1
))

allres <- fread(file.path(da_out, "all_results.tsv"))
say("MaAsLin2 all_results columns: %s", paste(names(allres), collapse = ", "))
stopifnot(all(c("feature", "metadata", "value", "coef", "pval", "qval") %in% names(allres)))
grp <- allres[metadata == "group"]                       # the R-vs-NR effect (value=R)
setorder(grp, qval)
grp[, direction := ifelse(coef > 0, "enriched in R", "enriched in NR")]
fwrite(grp, file.path(results_dir, paste0("DA_species_group_", COHORT, ".csv")))
n_sig05 <- sum(grp$qval < 0.05, na.rm = TRUE); n_sig10 <- sum(grp$qval < 0.10, na.rm = TRUE)
say("DA species tested: %d | q<0.10: %d | q<0.05: %d", nrow(grp), n_sig10, n_sig05)

# volcano — points colored by ICI-response GROUP (Responder vs Non-responder).
# Filter 1 = significance q < 0.05 (BH-FDR). Color by coef SIGN (NR is the model
# reference: coef>0 => higher in R, coef<0 => higher in NR). The optional effect-
# size gate COEF_GATE (filter 2) is OFF by default (0); set >0 for a dual-threshold.
SIG_Q <- 0.05; COEF_GATE <- 0
.LAB_R <- "Enriched in Responder (R)"; .LAB_NR <- "Enriched in Non-responder (NR)"; .LAB_NS <- "n.s. (q >= 0.05)"
grp[, neglogq := -log10(qval)]
grp[, hit := qval < SIG_Q & abs(coef) >= COEF_GATE]
grp[, ici_grp := ifelse(!hit, .LAB_NS, ifelse(coef > 0, .LAB_R, .LAB_NR))]
grp[, ici_grp := factor(ici_grp, levels = c(.LAB_R, .LAB_NR, .LAB_NS))]
volcano_cols <- setNames(c(unname(PAL["R"]), unname(PAL["NR"]), "grey80"), c(.LAB_R, .LAB_NR, .LAB_NS))
top <- head(grp[hit == TRUE][order(qval)], 12)
if (nrow(top) > 0) top[, lab := gsub("_", " ", feature)]
p_volc <- ggplot(grp, aes(coef, neglogq)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(yintercept = -log10(SIG_Q), linetype = "dashed", color = "grey60") +
  geom_point(aes(color = ici_grp), size = 2, alpha = 0.85) +
  scale_color_manual(values = volcano_cols, name = "ICI response", drop = FALSE) +
  labs(x = "MaAsLin2 coefficient   (<-- higher in NR    |    higher in R -->)",
       y = expression(-log[10]~"(q-value, BH)"),
       title = sprintf("NSCLC %s: differential species by ICI response", COHORT),
       subtitle = sprintf("MaAsLin2 (TSS+LOG+LM), adjusted for %s; colored by enrichment group (q < %.2f, BH); %d species (prev >= %d%%)",
                          paste(COVARS, collapse = "+"), SIG_Q, nrow(grp), round(PREV_MIN*100))) +
  theme_classic(base_size = 13) +
  theme(plot.title = element_text(face = "bold", size = 13), legend.position = "right")
if (nrow(top) > 0) {
  p_volc <- p_volc + ggrepel::geom_text_repel(
    data = top, aes(coef, neglogq, label = lab),
    size = 3, fontface = "italic", max.overlaps = 20,
    min.segment.length = 0, segment.color = "grey70")
}
# ragg (ggplot2 4.0 default PNG device) cannot open the non-ASCII (Korean) folder
# path on macOS; render to an ASCII temp file, then copy to the final destination.
.tmp_png <- tempfile(fileext = ".png")
ggsave(.tmp_png, p_volc, width = 7.2, height = 5.6, dpi = 300, bg = "white")
stopifnot(file.copy(.tmp_png, file.path(results_dir, paste0("DA_volcano_", COHORT, ".png")), overwrite = TRUE)); unlink(.tmp_png)

## --------------------------------------------------------------- provenance
writeLines(capture.output(sessionInfo()),
           file.path(results_dir, paste0("sessionInfo_", COHORT, ".txt")))
say("DONE. outputs in %s", results_dir)
close(logcon)
