#!/usr/bin/env Rscript
# =============================================================================
# Cross-cohort consistency of sarcosine-DEGRADATION-associated species
#   NSCLC ICI Responder-vs-Non-responder analysis  (R 4.6.0 environment)
# =============================================================================
# NSCLC-adapted analogue of the CRC cross_cohort_degradation_consistency.R Venn.
# This is a SEPARATE analysis from the CRC one -- different disease, contrast,
# cohorts, KO panel, and pooling. Results are NOT merged with the CRC Venn.
#
# DESIGN DECISIONS (chosen by the user, 2026-06-09):
#   (1) PRIMARY Venn = the 3 DISCOVERY cohorts (PRJNA751792 n=338, PRJNA1023797
#       n=421, PRJEB22863-NSCLC n=65); Korea PRJEB26531 (n=25) reported as a
#       VALIDATION annotation. A SUPPLEMENTARY 4-set Venn incl. Korea as a full
#       set is ALSO produced (cross_cohort_deg_assoc_venn4_NSCLC.png) -- Korea is
#       underpowered (n=25; 0 species under per-cohort BH), so its regions carry
#       a caution caveat in the caption.
#   (2) PER-COHORT CRITERION = Spearman rho(species abundance, degradation score)
#         > 0.3  AND  nominal p < 0.05,  on species with >=10% prevalence.
#       RATIONALE: CRC used per-cohort BH-FDR<0.05 & rho>0.3, but NSCLC's small
#       cohorts are underpowered for per-cohort FDR (PRJEB22863 n=65, PRJEB26531
#       n=25 -> 0 species after BH). The CROSS-COHORT OVERLAP is itself the
#       consistency / replication control. The FORMAL NSCLC test remains the
#       Fisher-z random-effects meta (pooled_analysis/results/
#       sarcosine_species_assoc_meta.csv). BH p_adj is still REPORTED per cohort
#       in the correlations CSV for full transparency.
#   (3) DEGRADATION SCORE = the verified NSCLC degradation-KO sum
#         (K00301-K00306, K00314, K18897) read from sarcosine_scores_<cohort>.csv
#       -- the NSCLC panel, NOT the CRC 6-KO panel.
#
# Reuses verified inputs only (sarcosine_scores_*.csv + Bacteria_*.txt); the
# per-cohort Spearman (rho + p) is recomputed here. Species matched on the FULL
# taxonomy string. seed 42; sessionInfo + log written; ragg-safe PNG write.
#
# Run: Rscript pooled_analysis/R_scripts/cross_cohort_degradation_venn.R
# =============================================================================
set.seed(42)
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(ggvenn) })

# ---- config -----------------------------------------------------------------
PREV_MIN <- 0.10
RHO_CUT  <- 0.30
P_CUT    <- 0.05            # nominal per-cohort p (see header rationale)
SCORE    <- "degradation"
DISC  <- c(PRJNA751792 = "NSCLC_PRJNA751792",
           PRJNA1023797 = "NSCLC_PRJNA1023797",
           PRJEB22863 = "NSCLC_RCC_PRJEB22863")   # discovery -> Venn sets
VALID <- c(PRJEB26531 = "NSCLC_PRJEB26531")       # Korea -> validation annotation
VENN_COLS <- c("#332288", "#CC6677", "#44AA99")   # Tol muted (indigo/rose/teal) -- intentionally DIFFERENT from CRC's Okabe-Ito orange/skyblue/green

# Run with the R working directory set to the analysis-folder root; paths below
# are relative to it.
pooled_dir <- file.path(".", "pooled_analysis"); main_dir <- "."
res <- file.path(pooled_dir, "results"); dir.create(res, showWarnings = FALSE, recursive = TRUE)
logcon <- file(file.path(res, "cross_cohort_deg_venn_log.txt"), open = "wt")
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logcon) }
save_png <- function(p, out, w, h) {   # ragg cannot write the non-ASCII (Korean) path on macOS
  tmp <- tempfile(fileext = ".png"); ggsave(tmp, p, width = w, height = h, dpi = 300, bg = "white")
  stopifnot(file.copy(tmp, out, overwrite = TRUE)); unlink(tmp)
}
short <- function(x) sub(".*\\|s__", "", x)

read_species <- function(cdir) {   # identical filters to run_analysis.R; FULL taxa names kept
  f <- list.files(cdir, "^Bacteria_.*\\.txt$", full.names = TRUE); stopifnot(length(f) == 1L)
  b <- fread(f, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
  stopifnot(ncol(b) == 3L); setnames(b, c("Taxa", "RunID", "Abundance"))
  sp <- b[grepl("|s__", Taxa, fixed = TRUE) & !grepl("|t__", Taxa, fixed = TRUE)]
  w  <- dcast(sp, RunID ~ Taxa, value.var = "Abundance", fill = 0, fun.aggregate = sum)
  m  <- as.matrix(w[, -1L, with = FALSE]); rownames(m) <- w$RunID; m
}

corr_cohort <- function(cn, dirname_) {   # per-cohort Spearman rho + p of each species vs degradation score
  cdir <- file.path(main_dir, dirname_)
  sc <- fread(file.path(cdir, "results", paste0("sarcosine_scores_", cn, ".csv")))
  stopifnot(SCORE %in% names(sc))
  spm <- read_species(cdir); stopifnot(all(sc$RunID %in% rownames(spm)))
  spm <- spm[sc$RunID, , drop = FALSE]
  spm <- spm[, colMeans(spm > 0) >= PREV_MIN, drop = FALSE]
  deg <- sc[[SCORE]]
  out <- rbindlist(lapply(colnames(spm), function(s) {
    x <- spm[, s]
    if (sd(x) == 0) return(data.table(Taxa = s, rho = NA_real_, p = NA_real_))
    ct <- suppressWarnings(cor.test(x, deg, method = "spearman", exact = FALSE))
    data.table(Taxa = s, rho = unname(ct$estimate), p = ct$p.value)
  }))
  out[, padj := p.adjust(p, "BH")]
  out[, `:=`(cohort = cn, n = nrow(spm), species = short(Taxa))]
  out[, .(cohort, Taxa, species, rho, p, padj, n)]
}

say("=== Cross-cohort sarcosine-DEGRADATION-associated species (NSCLC, ICI R vs NR) ===")
say("criterion: rho > %.2f AND nominal p < %.2f | prevalence >= %.0f%% | score = %s", RHO_CUT, P_CUT, 100*PREV_MIN, SCORE)
say("Venn sets = discovery (%s); Korea %s = validation annotation\n", paste(names(DISC), collapse=", "), names(VALID))

# ---- per-cohort correlations (3 discovery + Korea) --------------------------
all_corr <- rbindlist(lapply(names(c(DISC, VALID)), function(cn) {
  d <- corr_cohort(cn, c(DISC, VALID)[cn])
  say("  %-12s n=%d | species tested=%d | rho>%.1f & p<%.2f = %d | (BH p_adj<0.05 = %d)",
      cn, d$n[1], nrow(d), RHO_CUT, P_CUT, sum(d$rho > RHO_CUT & d$p < P_CUT, na.rm=TRUE),
      sum(d$rho > RHO_CUT & d$padj < 0.05, na.rm=TRUE))
  d
}))
fwrite(all_corr, file.path(res, "cross_cohort_deg_assoc_correlations_NSCLC.csv"))

is_assoc <- function(d) d$Taxa[!is.na(d$rho) & !is.na(d$p) & d$rho > RHO_CUT & d$p < P_CUT]
sets <- lapply(names(DISC), function(cn) is_assoc(all_corr[cohort == cn])); names(sets) <- names(DISC)
korea_set <- is_assoc(all_corr[cohort == names(VALID)])

# ---- regions / intersection -------------------------------------------------
A <- sets[[1]]; B <- sets[[2]]; C <- sets[[3]]
union3 <- Reduce(union, sets); inter3 <- Reduce(intersect, sets)
reg <- list(
  A_only = setdiff(A, union(B, C)), B_only = setdiff(B, union(A, C)), C_only = setdiff(C, union(A, B)),
  AB = setdiff(intersect(A, B), C), AC = setdiff(intersect(A, C), B), BC = setdiff(intersect(B, C), A),
  ABC = inter3)
cnt <- sapply(reg, length)
say("\nset sizes: %s = %d / %d / %d", paste(names(DISC), collapse=" / "), length(A), length(B), length(C))
say("regions: A_only=%d B_only=%d C_only=%d | AB=%d AC=%d BC=%d | ABC=%d", cnt["A_only"],cnt["B_only"],cnt["C_only"],cnt["AB"],cnt["AC"],cnt["BC"],cnt["ABC"])
say("total unique (discovery union) = %d | 3-way intersection = %d", length(union3), length(inter3))

# ---- per-cohort significant lists + summary CSVs ----------------------------
per_cohort_sig <- rbindlist(lapply(names(DISC), function(cn) {
  d <- all_corr[cohort == cn & Taxa %in% sets[[cn]]][order(-rho)]
  d[, .(cohort, Taxa, species, rho, p, padj, n)]
}))
fwrite(per_cohort_sig, file.path(res, "cross_cohort_deg_assoc_per_cohort_lists_NSCLC.csv"))

summary_df <- rbindlist(lapply(names(DISC), function(cn) {
  d <- all_corr[cohort == cn]
  data.table(cohort = cn, n_samples = d$n[1], n_species_tested = nrow(d),
             n_rho_gt_cut_and_p = length(sets[[cn]]),
             n_BH_padj_lt0.05_and_rho = sum(d$rho > RHO_CUT & d$padj < 0.05, na.rm = TRUE))
}))
summary_df <- rbind(summary_df,
  data.table(cohort = "Korea_PRJEB26531(validation)", n_samples = all_corr[cohort==names(VALID)]$n[1],
             n_species_tested = nrow(all_corr[cohort==names(VALID)]),
             n_rho_gt_cut_and_p = length(korea_set), n_BH_padj_lt0.05_and_rho = NA_integer_),
  data.table(cohort = "Discovery_3way_intersection", n_samples = NA_integer_, n_species_tested = NA_integer_,
             n_rho_gt_cut_and_p = length(inter3), n_BH_padj_lt0.05_and_rho = NA_integer_))
fwrite(summary_df, file.path(res, "cross_cohort_deg_assoc_summary_NSCLC.csv"))

# ---- 3-way intersection table + KOREA VALIDATION ----------------------------
inter_tab <- NULL
if (length(inter3) > 0) {
  inter_tab <- rbindlist(lapply(inter3, function(tx) {
    row <- data.table(Taxa = tx, species = short(tx))
    for (cn in names(DISC)) { d <- all_corr[cohort == cn & Taxa == tx]
      row[[paste0(cn, "_rho")]] <- d$rho; row[[paste0(cn, "_p")]] <- d$p }
    kv <- all_corr[cohort == names(VALID) & Taxa == tx]
    row[["Korea_rho"]] <- if (nrow(kv)) kv$rho else NA_real_
    row[["Korea_p"]]   <- if (nrow(kv)) kv$p   else NA_real_
    row[["Korea_concordant_pos"]] <- if (nrow(kv)) (is.finite(kv$rho) && kv$rho > 0) else NA
    row
  }))
  disc_rho_cols <- paste0(names(DISC), "_rho")
  inter_tab[, mean_disc_rho := rowMeans(.SD), .SDcols = disc_rho_cols]
  setorder(inter_tab, -mean_disc_rho)
  fwrite(inter_tab, file.path(res, "cross_cohort_deg_assoc_intersection_discovery_NSCLC.csv"))
  say("\n3-way intersection species (Korea validation rho shown):")
  for (i in seq_len(nrow(inter_tab)))
    say("  %-34s mean discovery rho=%.2f | Korea rho=%.2f (p=%.3g)%s",
        inter_tab$species[i], inter_tab$mean_disc_rho[i], inter_tab$Korea_rho[i], inter_tab$Korea_p[i],
        ifelse(isTRUE(inter_tab$Korea_concordant_pos[i]), " [+concordant]", ""))
}

# ---- draw the Venns with ggvenn (same package as CRC; distinct Tol-muted palette)
# Two figures, both restyled vs the CRC Venn (palette indigo/rose/teal[/sand],
# bigger title + labels, NO subtitle, concise caption): (1) PRIMARY 3-set
# discovery Venn (Korea = validation); (2) SUPPLEMENTARY 4-set Venn incl. Korea
# as a full set (Korea n=25 underpowered -> caption caveat).
inter4  <- Reduce(intersect, c(sets, list(korea_set)))     # species sig in ALL 4 cohorts
sp_join <- function(tx) if (length(tx)) paste(gsub("_", " ", short(tx)), collapse = ", ") else "(none)"
style_venn <- function(g, title, caption) {
  g + labs(title = title, caption = caption) +
    theme(plot.title   = element_text(size = 22, face = "bold", hjust = 0.5, margin = margin(b = 12)),
          plot.caption = element_text(size = 10, hjust = 0, color = "grey30", margin = margin(t = 8)),
          plot.margin  = margin(16, 18, 14, 18))
}

# (1) PRIMARY 3-set (discovery) -------------------------------------------------
korea_allpos <- length(inter3) > 0 && all(inter_tab$Korea_concordant_pos, na.rm = TRUE)
cap3 <- paste0(
  "Per cohort: Spearman rho > 0.3 & nominal p < 0.05 (prevalence >= 10%).\n",
  "Shared by all 3", if (korea_allpos) " (Korea n=25 validation: all rho > 0)" else "", ": ", sp_join(inter3), ".\n",
  "Overlap = cross-cohort consistency; formal test = Fisher-z meta.  Separate analysis from the CRC Venn.")
p3 <- style_venn(
  ggvenn(sets, fill_color = VENN_COLS, fill_alpha = 0.55, stroke_color = "grey25",
         stroke_size = 0.5, set_name_size = 8.5, text_size = 9, show_percentage = FALSE) +
    annotate("text", x = 0, y = 0.2, label = as.character(length(inter3)),     # highlight the 3-cohort central intersection in red (matches the CRC Venns)
             colour = "red", fontface = "bold", size = 9),
  "Sarcosine degradation-associated species\nacross 3 NSCLC discovery cohorts", cap3)
save_png(p3, file.path(res, "cross_cohort_deg_assoc_venn_NSCLC.png"), 10, 10)
say("Saved: cross_cohort_deg_assoc_venn_NSCLC.png (3-set discovery, ggvenn)")

# (2) SUPPLEMENTARY 4-set (incl. Korea as a full set) ---------------------------
sets4 <- c(sets, setNames(list(korea_set), names(VALID)))
cap4 <- paste0(
  "Per cohort: Spearman rho > 0.3 & nominal p < 0.05 (prevalence >= 10%).\n",
  "Korea PRJEB26531 (n=25) underpowered -> interpret Korea regions with caution.  Shared by all 4: ", sp_join(inter4), ".\n",
  "Formal test = Fisher-z meta.  Separate analysis from the CRC Venn.")
p4 <- style_venn(
  ggvenn(sets4, fill_color = c(VENN_COLS, "#DDCC77"), fill_alpha = 0.55, stroke_color = "grey25",
         stroke_size = 0.5, set_name_size = 7.5, text_size = 7, show_percentage = FALSE) +
    annotate("text", x = 0, y = -0.7, label = as.character(length(inter4)),    # highlight the 4-cohort central intersection in red (matches the CRC Venns)
             colour = "red", fontface = "bold", size = 7),
  "Sarcosine degradation-associated species\nacross 4 NSCLC cohorts", cap4)
save_png(p4, file.path(res, "cross_cohort_deg_assoc_venn4_NSCLC.png"), 11, 11)
say("Saved: cross_cohort_deg_assoc_venn4_NSCLC.png (4-set incl. Korea, ggvenn)")

# 4-way (all 4 cohorts) intersection table --------------------------------------
if (length(inter4)) {
  it4 <- rbindlist(lapply(inter4, function(tx) {
    row <- data.table(Taxa = tx, species = short(tx))
    for (cn in c(names(DISC), names(VALID))) { d <- all_corr[cohort == cn & Taxa == tx]
      row[[paste0(cn, "_rho")]] <- d$rho; row[[paste0(cn, "_p")]] <- d$p }
    row }))
  fwrite(it4, file.path(res, "cross_cohort_deg_assoc_intersection_4cohorts_NSCLC.csv"))
}
say("4-way (all 4 cohorts) intersection: %d species -> %s", length(inter4), sp_join(inter4))

writeLines(capture.output(sessionInfo()), file.path(res, "sessionInfo_cross_cohort_deg_venn.txt"))
say("DONE. outputs in %s", res)
close(logcon)
