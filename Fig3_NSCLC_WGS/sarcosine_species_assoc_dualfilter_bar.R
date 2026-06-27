#!/usr/bin/env Rscript
# ============================================================================
# SARCOSINE FUNCTION <-> SPECIES association, in the ICI R-vs-NR context.
#
# For each gut species: Spearman correlation between its relative abundance and
# the per-sample sarcosine DEGRADATION / PRODUCTION score, computed PER COHORT,
# then combined by a random-effects (DL) meta-analysis of Fisher-z-transformed
# correlations over the 3 DISCOVERY cohorts (Korea PRJEB26531 = external
# validation, reported in the CSV only). Output: one lollipop per score -- the
# top species by meta FDR q, x = pooled Spearman rho, COLOURED by the ICI-response
# group the species is enriched in (from the verified pooled_species_meta.csv DA).
# Single legend.
#
# DESIGN (faithful-coder):
#   - unit = sequencing run; R vs NR are the NSCLC-WGS samples (as elsewhere). The
#     R-vs-NR contrast enters via the COLOUR (each species' DA enrichment group).
#   - Spearman (rank) correlation: the scores are zero-inflated / compositional and
#     species relative abundances are right-skewed -> Spearman, not Pearson.
#   - prevalence gate: a species must be present in >=10% of a cohort's samples
#     (matches the DA filter) to be correlated in that cohort (avoids spurious
#     correlations dominated by shared zeros).
#   - meta: per-cohort rho -> Fisher z = atanh(rho), se = 1/sqrt(n-3) (standard
#     Fisher-z SE); metafor::rma(method="DL"); pooled rho = tanh(pooled z); I^2;
#     BH-FDR across species, within each score. Species testable in >=2 discovery
#     cohorts are meta-analysed.
#
# INTEGRITY / INTERPRETATION: species (MetaPhlAn4 taxonomy) and the KO-derived
#   scores come from the SAME samples; a species that CARRIES sarcosine genes will
#   mechanically track the score. This identifies the taxonomic CONTRIBUTORS /
#   co-occurring taxa of the function -- association, NOT causation or an independent
#   validation. Reads verified inputs only (sarcosine_scores_*.csv, Bacteria_*.txt,
#   pooled_species_meta.csv); no upstream statistic is modified.
#
# Run: Rscript sarcosine_species_assoc_dualfilter_bar.R
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(data.table); library(metafor); library(ggplot2) })

SCORES   <- c("degradation", "production")
SCORELAB <- c(degradation = "degradation score", production = "production score")
PREV_MIN <- 0.10
TOP_N    <- 20
SIG_Q    <- 0.05
PAL      <- c(R = "#1B7837", NR = "#B2182B")     # group scheme: R=green (favorable), NR=red (unfavorable)
LAB_R  <- "Enriched in Responder (R)"
LAB_NR <- "Enriched in Non-responder (NR)"
LAB_NS <- "not DA-significant (q >= 0.05)"
COLS   <- setNames(c(unname(PAL["R"]), unname(PAL["NR"]), "grey75"), c(LAB_R, LAB_NR, LAB_NS))

DISC  <- c(PRJNA751792 = "NSCLC_PRJNA751792", PRJNA1023797 = "NSCLC_PRJNA1023797", PRJEB22863 = "NSCLC_RCC_PRJEB22863")
VALID <- c(PRJEB26531 = "NSCLC_PRJEB26531")

# Run with the R working directory set to the analysis-folder root; paths below
# are relative to it.
main_dir <- "."
res <- file.path(main_dir, "pooled_analysis", "results")
logcon <- file(file.path(res, "sarcosine_species_assoc_log.txt"), open = "wt")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logcon) }

save_png <- function(p, out_png, width, height) {   # ragg cannot write the non-ASCII (Korean) path on macOS
  tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = width, height = height, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out_png, overwrite = TRUE)); unlink(tmp)
}

read_species <- function(cdir) {   # identical logic to run_analysis.R (cleaned, unique species names)
  f <- list.files(cdir, "^Bacteria_.*\\.txt$", full.names = TRUE); stopifnot(length(f) == 1L)
  b <- fread(f, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
  stopifnot(ncol(b) == 3L); setnames(b, c("Taxa", "RunID", "Abundance"))
  sp <- b[grepl("|s__", Taxa, fixed = TRUE) & !grepl("|t__", Taxa, fixed = TRUE)]
  w  <- dcast(sp, RunID ~ Taxa, value.var = "Abundance", fill = 0, fun.aggregate = sum)
  m  <- as.matrix(w[, -1L, with = FALSE]); rownames(m) <- w$RunID
  colnames(m) <- make.unique(sub(".*\\|s__", "", colnames(m)))
  m
}

## ---- per-cohort Spearman rho(species, score) -------------------------------
percoh <- list()
for (cn in c(names(DISC), names(VALID))) {
  cdir <- file.path(main_dir, if (cn %in% names(DISC)) DISC[cn] else VALID[cn])
  sc <- fread(file.path(cdir, "results", paste0("sarcosine_scores_", cn, ".csv")))   # RunID, group, degradation, production, ...
  spm <- read_species(cdir)
  stopifnot(all(sc$RunID %in% rownames(spm)))
  spm <- spm[sc$RunID, , drop = FALSE]                        # align species rows to the score samples
  keep <- colMeans(spm > 0) >= PREV_MIN                       # >=10% prevalence in THIS cohort
  spm <- spm[, keep, drop = FALSE]
  for (s in SCORES) {
    rho <- suppressWarnings(as.numeric(cor(spm, sc[[s]], method = "spearman")))
    d <- data.table(cohort = cn, score = s, species = colnames(spm), rho = rho, n = nrow(spm))
    percoh[[paste(cn, s)]] <- d[is.finite(rho)]
  }
  say("%s: %d species (>=%.0f%% prev) x %d scores correlated (n=%d)", cn, sum(keep), 100*PREV_MIN, length(SCORES), nrow(spm))
}
PC <- rbindlist(percoh)
fwrite(PC, file.path(res, "sarcosine_species_assoc_percohort.csv"))

## ---- random-effects (DL) meta of Fisher-z correlations over discovery cohorts
meta_one <- function(d) {
  d <- d[is.finite(rho) & abs(rho) < 1 & n > 3]
  if (nrow(d) < 2) return(NULL)
  d[, z := atanh(rho)][, se := 1/sqrt(n - 3)]
  m <- tryCatch(rma(yi = z, sei = se, data = d, method = "DL"), error = function(e) NULL)
  if (is.null(m)) return(NULL)
  data.table(k = nrow(d), pooled_rho = tanh(as.numeric(m$beta)),
             ci_lb = tanh(m$ci.lb), ci_ub = tanh(m$ci.ub), pval = m$pval, I2 = m$I2)
}
res_meta <- rbindlist(lapply(SCORES, function(s) {
  dd  <- PC[score == s & cohort %in% names(DISC)]
  spp <- dd[, .N, by = species][N >= 2, species]
  rbindlist(lapply(spp, function(sp) {
    mm <- meta_one(dd[species == sp]); if (is.null(mm)) return(NULL)
    cbind(data.table(score = s, species = sp), mm)
  }))
}))
res_meta[, qval := p.adjust(pval, "BH"), by = score]

## Korea (validation) rho -- concordance, reported in the CSV only
kv <- PC[cohort == names(VALID), .(score, species, valid_rho = rho)]
res_meta <- merge(res_meta, kv, by = c("score", "species"), all.x = TRUE)

## ---- ICI-response enrichment colour, from the verified DA meta -------------
da <- fread(file.path(res, "pooled_species_meta.csv"))[, .(species, da_dir = dir, da_q = qval)]
res_meta <- merge(res_meta, da, by = "species", all.x = TRUE)
res_meta[, grp := fifelse(!is.na(da_q) & da_q < SIG_Q & da_dir == "higher in R",  LAB_R,
                  fifelse(!is.na(da_q) & da_q < SIG_Q & da_dir == "higher in NR", LAB_NR, LAB_NS))]
setorder(res_meta, score, qval)
fwrite(res_meta, file.path(res, "sarcosine_species_assoc_meta.csv"))
say("meta hits (q<0.05): %s", paste(SCORES, vapply(SCORES, function(s) sprintf("%s=%d/%d", s,
    sum(res_meta[score==s]$qval < SIG_Q), nrow(res_meta[score==s])), character(1)), collapse = " | "))

## ---- species SELECTION (by sarcosine-score correlation) + DA BAR -----------
## IMPORTANT: this is NOT a strict dual filter (despite the file name, kept for
## continuity). Selection uses ONLY filter (a) of the CRC/Wirbel-Thomas scheme:
## species significantly CORRELATED with the degradation and/or production score
## (meta Spearman FDR < 0.05). The CRC second dimension -- R-vs-NR DIFFERENTIAL
## ABUNDANCE -- is NOT an inclusion gate here; it is the BAR itself (meta MaAsLin2
## coef) with its BH q shown as a p-label (BLANK when q >= 0.05). Consequently some
## plotted species are correlation-significant but DA-NON-significant (audit
## 2026-06-07: of the top-20 by |rho|, 8 are DA q<0.05 and 12 are not). An
## exact-CRC strict dual filter (corr FDR<0.05 AND DA FDR<0.05) can be made on request.
## ============================================================================
## CRC-STYLE plot (matches .../sarcosine_bacteria_diff_abundance_pooled.png):
##   - SELECT the top species by |rho| from the sarcosine-KO correlation
##   - BAR = differential abundance (meta MaAsLin2 coef, R vs NR), coloured by
##     DIRECTION (R / NR enriched)
##   - right-hand SARCOSINE-PATHWAY annotation (Deg / Prod / Both = which sarcosine
##     function the species significantly correlates with)
##   - BH-adjusted DA p-value labels; two legends.
## Bars show DIFFERENTIAL ABUNDANCE; selection is by correlation. (KO indicators unchanged.)
## ============================================================================
TOPN <- 20
w <- dcast(res_meta, species ~ score, value.var = c("pooled_rho", "qval"))
w[, deg_sig  := !is.na(qval_degradation) & qval_degradation < SIG_Q]
w[, prod_sig := !is.na(qval_production)  & qval_production  < SIG_Q]
w[, pathway := fifelse(deg_sig & prod_sig, "Deg & Prod", fifelse(deg_sig, "Deg", fifelse(prod_sig, "Prod", NA_character_)))]
w[, max_rho := pmax(abs(pooled_rho_degradation), abs(pooled_rho_production), na.rm = TRUE)]
da_full <- fread(file.path(res, "pooled_species_meta.csv"))[, .(species, da_coef = pooled_coef, da_q = qval)]
assoc <- merge(w[!is.na(pathway)], da_full, by = "species")[order(-max_rho)]
assoc <- head(assoc, TOPN)
say("CRC-style DA plot: %d sarcosine-associated species (sig corr Deg/Prod) -> top %d by |rho|", nrow(w[!is.na(pathway)]), nrow(assoc))

assoc[, direction := factor(fifelse(da_coef > 0, "Enriched in R", "Enriched in NR"), levels = c("Enriched in R", "Enriched in NR"))]
assoc[, pathway := factor(pathway, levels = c("Deg", "Deg & Prod", "Prod"))]
assoc <- assoc[order(-da_coef)]   # display order: NR-enriched at TOP, R-enriched at BOTTOM (match the species lollipop)
assoc[, sp2 := factor(gsub("_", " ", species), levels = gsub("_", " ", species))]
assoc[, plab := fifelse(da_q < 1e-4, "p<0.0001", fifelse(da_q < 0.05, sprintf("p=%.3g", da_q), ""))]
rng <- max(abs(assoc$da_coef)); xr <- -rng * 1.32; xt <- -rng * 1.42   # negative side -> displays on the RIGHT after scale_x_reverse
assoc[, px := da_coef + ifelse(da_coef >= 0, 1, -1) * rng * 0.04]

DIR_COLS  <- c("Enriched in R" = "#1B7837", "Enriched in NR" = "#B2182B")           # group scheme: R=green (favorable), NR=red (unfavorable)
PATH_COLS <- c("Deg" = "#117733", "Deg & Prod" = "#882255", "Prod" = "#E6A000")     # distinct from direction

p <- ggplot(assoc, aes(da_coef, sp2)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55") +
  geom_col(aes(fill = direction), width = 0.68) +
  geom_text(aes(x = px, label = plab), hjust = ifelse(assoc$da_coef >= 0, 1, 0), size = 3.8, fontface = "italic", colour = "grey25") +
  geom_point(aes(x = xr, colour = pathway), shape = 15, size = 4.6) +
  geom_text(aes(x = xt, label = pathway, colour = pathway), hjust = 0, size = 3.9, fontface = "bold", show.legend = FALSE) +
  scale_fill_manual(values = DIR_COLS, name = "Direction", drop = FALSE) +
  scale_colour_manual(values = PATH_COLS, name = "Sarcosine Pathway", drop = FALSE) +
  scale_x_reverse(limits = c(max(assoc$da_coef) + rng * 0.6, -rng * 2.7), expand = c(0, 0)) +
  labs(x = "MaAsLin2 coefficient, R vs NR   (← enriched in R    |    enriched in NR →)", y = NULL,
       title = "Sarcosine-Associated Bacteria: Differential Abundance (R vs NR)",
       subtitle = sprintf("Pooled discovery (3 cohorts) | top %d species by |rho| from sarcosine-KO correlation", nrow(assoc)),
       caption = "Bars: meta MaAsLin2 coef (R vs NR)  |  p: BH-adjusted DA q  |  right column: sarcosine pathway the species correlates with (Deg / Prod / both)") +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 18),
        plot.subtitle = element_text(colour = "grey35", size = 13),
        plot.caption = element_text(colour = "grey45", size = 10, hjust = 0.5),
        axis.text.y = element_text(face = "italic", size = 12),
        panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
        legend.position = "bottom", legend.box = "vertical")
save_png(p, file.path(res, "sarcosine_associated_bacteria_DA_pooled.png"), width = 13, height = 8)
say("wrote sarcosine_associated_bacteria_DA_pooled.png (%d species)", nrow(assoc))
writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_sarcosine_species_assoc.txt"))
say("DONE. outputs in %s", res)
close(logcon)
