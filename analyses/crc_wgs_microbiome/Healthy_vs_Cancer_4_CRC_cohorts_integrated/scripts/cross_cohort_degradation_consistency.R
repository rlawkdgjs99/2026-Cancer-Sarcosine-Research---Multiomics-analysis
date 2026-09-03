#!/usr/bin/env Rscript
# =============================================================================
# Cross-cohort consistency of sarcosine-DEGRADATION-associated species
# =============================================================================
# Purpose:
#   Methodological TWIN of GOAL B in cross_cohort_consistency.R, but for the
#   DEGRADATION arm of sarcosine metabolism instead of production. For each
#   cohort it finds the species whose relative abundance correlates with the
#   community's sarcosine-degradation gene capacity (Degradation_sum), then
#   visualizes the cross-cohort overlap as a 4-set Venn diagram so it is
#   directly comparable to cross_cohort_prod_assoc_venn.png.
#
# What is IDENTICAL to Goal B (so the two Venns are comparable):
#   - species: species-level (|s__ & not |t__), >=10% prevalence on the
#     per-cohort correlation samples
#   - correlation: Spearman(species abundance, sarcosine pathway sum)
#   - multiple testing: Benjamini-Hochberg, computed WITHIN each cohort
#   - significance: p_adj < 0.05 AND rho > 0.3 (see DIR for the direction choice)
#   - per-cohort then intersect across the 4 cohorts
#
# What CHANGES vs Goal B:
#   (1) DEGRADATION KOs (Role=="Degradation") instead of Production KOs
#   (2) output filenames: deg_assoc instead of prod_assoc
#
# Efficiency vs. faithfulness (KO side):
#   The KO side is taken from the precomputed per-sample cache
#     results_integrated/sarcosine/sarcosine_KO_per_sample_pooled.csv
#   which is the SAME 1,647-sample analysis set as Goal B (128+260+1066+193).
#   This avoids re-reading the ~700 MB KO JSONs. Faithfulness is enforced by
#   hard stopifnot() checks that the per-cohort common-sample counts equal
#   Goal B's (128/260/1066/193) and the per-cohort species counts equal Goal
#   B's (557/484/378/396). The species side is re-read from the same
#   Bacteria_*.txt files Goal B used, with the identical filters.
#
# Direction of association (DESIGN CHOICE):
#   DIR = "positive" (default) -> rho >  0.3 : species that INCREASE with
#         degradation capacity = the exact mirror of Goal B's production rule.
#   DIR = "negative"           -> rho < -0.3 : species that DECREASE with it.
#   DIR = "both"               -> |rho| > 0.3.
#
# Outputs (results_integrated/cross_cohort/):
#   cross_cohort_deg_assoc_venn.png
#   cross_cohort_deg_assoc_correlations.csv
#   cross_cohort_deg_assoc_per_cohort_lists.csv
#   cross_cohort_deg_assoc_membership_matrix.csv
#   cross_cohort_deg_assoc_intersection_4cohorts.csv
#   cross_cohort_deg_assoc_summary.csv
#
# Author: added 2026-06-09 as an extension of cross_cohort_consistency.R
# =============================================================================

options(bitmapType = "quartz")   # cairo/X11 unavailable on this Mac; quartz works
set.seed(42)                     # reproducibility (Spearman is deterministic; fixed anyway)

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggvenn)
})

# ---- CONFIG (edit BASE_DIR if the folder moves) ------------------------------
BASE_DIR <- "Healthy_vs_Cancer_4_CRC_cohorts_integrated"   # relative to the analysis-folder root (working dir)

DIR      <- "positive"   # "positive" (mirror of Goal B) | "negative" | "both"
RHO_CUT  <- 0.3
PADJ_CUT <- 0.05
PREV_CUT <- 0.10

# Degradation KO panel (verbatim Role=="Degradation" from the existing scripts)
deg_kos_all <- c("K00301", "K00302", "K00303", "K00304", "K00305", "K00306")

# Cohorts live INSIDE BASE_DIR (same layout as integrated_analysis_pooled.R)
COHORT_DIRS <- list(
  PRJEB10878  = file.path(BASE_DIR, "PRJEB10878_CRC"),
  PRJEB27928  = file.path(BASE_DIR, "PRJEB27928_CRC"),
  PRJEB6070   = file.path(BASE_DIR, "PRJEB6070_CRC_AdenomatousPolyps"),
  PRJNA429097 = file.path(BASE_DIR, "PRJNA429097_CRC")
)
# Goal B's verified per-cohort numbers (integrity anchors)
EXPECTED_COMMON  <- c(PRJEB10878 = 128, PRJEB27928 = 260, PRJEB6070 = 1066, PRJNA429097 = 193)
EXPECTED_SPECIES <- c(PRJEB10878 = 557, PRJEB27928 = 484, PRJEB6070 = 378, PRJNA429097 = 396)

OUTDIR <- file.path(BASE_DIR, "results_integrated/cross_cohort")
CACHE  <- file.path(BASE_DIR, "results_integrated/sarcosine/sarcosine_KO_per_sample_pooled.csv")
stopifnot(dir.exists(OUTDIR), file.exists(CACHE))

venn_colors <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")  # same palette as Goal B

# quartz PNG wrapper (ggsave would try the unavailable cairo device)
.ggsave <- function(filename, plot, width, height, dpi = 300) {
  grDevices::png(filename = filename, width = width, height = height,
                 units = "in", res = dpi, type = "quartz", bg = "white")
  print(plot); grDevices::dev.off()
}

short_species <- function(x) {
  parts <- strsplit(x, "\\|")[[1]]
  sp <- parts[grep("^s__", parts)]
  if (length(sp) > 0) gsub("^s__", "", sp[1]) else x
}
dir_keep <- function(rho) {
  if (DIR == "positive") rho >  RHO_CUT
  else if (DIR == "negative") rho < -RHO_CUT
  else abs(rho) > RHO_CUT
}

cat("=== Cross-cohort sarcosine-DEGRADATION-associated species ===\n")
cat("Direction:", DIR, "| rho:", RHO_CUT, "| p_adj:", PADJ_CUT, "| prevalence:", PREV_CUT, "\n\n")

# ---- 1. KO side from cache: Degradation_sum per sample -----------------------
cache <- read.csv(CACHE, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("Run.ID", "Group", "Cohort") %in% colnames(cache)))

deg_kos_present <- intersect(deg_kos_all, colnames(cache))
cat("Degradation KOs in cache:", paste(deg_kos_present, collapse = ", "),
    " | missing:", paste(setdiff(deg_kos_all, deg_kos_present), collapse = ", "), "\n")
stopifnot(length(deg_kos_present) >= 1)

cache$Degradation_sum <- rowSums(cache[, deg_kos_present, drop = FALSE])
deg_sum_all <- setNames(cache$Degradation_sum, cache$Run.ID)

cohort_counts <- as.integer(table(cache$Cohort)[names(EXPECTED_COMMON)])
cat("\nCache per-cohort counts:", paste(names(EXPECTED_COMMON), cohort_counts, sep = "=", collapse = ", "), "\n")
stopifnot(identical(cohort_counts, as.integer(EXPECTED_COMMON)))
cat("[VERIFIED] cache sample set == Goal B correlation set (128/260/1066/193).\n\n")

# degradation KOs actually non-zero per cohort (mirror Goal B's prod_kos_used reporting)
deg_kos_used_by_cohort <- sapply(names(COHORT_DIRS), function(co) {
  sub <- cache[cache$Cohort == co, deg_kos_present, drop = FALSE]
  paste(deg_kos_present[colSums(sub, na.rm = TRUE) > 0], collapse = "+")
})

# ---- 2. Per-cohort species ~ Degradation_sum correlation ---------------------
cohort_corr_results <- list()
cohort_deg_sets     <- list()
cohort_summary      <- list()

for (cohort in names(COHORT_DIRS)) {
  cat("--- ", cohort, " ---\n", sep = "")
  cpath <- COHORT_DIRS[[cohort]]
  cohort_samples <- cache$Run.ID[cache$Cohort == cohort]

  bact_file <- list.files(cpath, pattern = "^Bacteria_.*\\.txt$", full.names = TRUE)[1]
  stopifnot(!is.na(bact_file))
  bact_raw <- read.delim(bact_file, stringsAsFactors = FALSE)
  colnames(bact_raw) <- trimws(colnames(bact_raw))
  stopifnot(all(c("Taxa", "Run.ID", "Abundance") %in% colnames(bact_raw)))

  bact_species <- bact_raw %>%
    filter(grepl("\\|s__", Taxa) & !grepl("\\|t__", Taxa)) %>%
    filter(Run.ID %in% cohort_samples)
  bact_wide <- bact_species %>%
    select(Taxa, Run.ID, Abundance) %>%
    pivot_wider(names_from = Taxa, values_from = Abundance, values_fill = 0)
  bact_mat_c <- as.data.frame(bact_wide)
  rownames(bact_mat_c) <- bact_mat_c$Run.ID
  bact_mat_c$Run.ID <- NULL
  bact_mat_c <- as.matrix(bact_mat_c)
  rm(bact_raw, bact_species, bact_wide); gc(verbose = FALSE)

  # common samples = bacteria ∩ cache  (== Goal B's common_c)
  common_c <- intersect(rownames(bact_mat_c), cohort_samples)
  cat("  common (bacteria ∩ KO):", length(common_c),
      " (Goal B:", EXPECTED_COMMON[cohort], ")\n")
  stopifnot(length(common_c) == EXPECTED_COMMON[cohort])

  bact_mat_c <- bact_mat_c[common_c, , drop = FALSE]
  deg_sum_c  <- deg_sum_all[common_c]
  stopifnot(!any(is.na(deg_sum_c)))

  prev <- colSums(bact_mat_c > 0) / nrow(bact_mat_c)
  bact_filt_c <- bact_mat_c[, prev >= PREV_CUT, drop = FALSE]
  cat("  species after >=", PREV_CUT * 100, "% prevalence: ", ncol(bact_filt_c),
      " (Goal B:", EXPECTED_SPECIES[cohort], ")\n", sep = "")
  stopifnot(ncol(bact_filt_c) == EXPECTED_SPECIES[cohort])

  corr_df_c <- bind_rows(lapply(seq_len(ncol(bact_filt_c)), function(i) {
    sp_abund <- bact_filt_c[, i]
    if (sd(sp_abund) == 0) {
      return(data.frame(Species_full = colnames(bact_filt_c)[i],
                        rho = NA_real_, p_value = NA_real_, stringsAsFactors = FALSE))
    }
    ct <- suppressWarnings(cor.test(sp_abund, deg_sum_c, method = "spearman", exact = FALSE))
    data.frame(Species_full = colnames(bact_filt_c)[i],
               rho = unname(ct$estimate), p_value = ct$p.value, stringsAsFactors = FALSE)
  }))
  corr_df_c$p_adj   <- p.adjust(corr_df_c$p_value, method = "BH")
  corr_df_c$Species <- sapply(corr_df_c$Species_full, short_species)
  corr_df_c$Cohort  <- cohort
  corr_df_c$deg_kos_summed <- deg_kos_used_by_cohort[cohort]
  corr_df_c$n_samples <- length(common_c)
  corr_df_c <- corr_df_c %>%
    select(Cohort, Species_full, Species, rho, p_value, p_adj, deg_kos_summed, n_samples)

  sig_c <- corr_df_c %>% filter(p_adj < PADJ_CUT & dir_keep(rho))
  cat("  degradation-associated species (", DIR, " rho, p_adj<", PADJ_CUT, "): ",
      nrow(sig_c), "\n\n", sep = "")

  cohort_corr_results[[cohort]] <- corr_df_c
  cohort_deg_sets[[cohort]]     <- sig_c$Species_full
  cohort_summary[[cohort]] <- data.frame(
    Cohort = cohort,
    n_samples_in_correlation = length(common_c),
    n_species_tested = nrow(corr_df_c),
    n_significant_padj_only = sum(corr_df_c$p_adj < PADJ_CUT, na.rm = TRUE),
    n_degradation_associated_strict = nrow(sig_c),
    deg_kos_used = deg_kos_used_by_cohort[cohort],
    direction = DIR,
    stringsAsFactors = FALSE
  )
}

# ---- 3. Cross-cohort union / intersection ------------------------------------
union_deg      <- Reduce(union, cohort_deg_sets)
deg_intersect4 <- Reduce(intersect, cohort_deg_sets)
overlap_counts <- sapply(union_deg, function(sp) sum(sapply(cohort_deg_sets, function(v) sp %in% v)))
cat("Total unique degradation-associated species:", length(union_deg), "\n")
cat("Distribution (in how many cohorts each is significant):\n"); print(table(overlap_counts))
cat("4-way intersection size:", length(deg_intersect4), "\n\n")

# ---- 4. Outputs --------------------------------------------------------------
# 4.1 Venn
venn_plot <- ggvenn(
  cohort_deg_sets, fill_color = venn_colors,
  stroke_size = 0.6, set_name_size = 8, text_size = 7, show_percentage = FALSE
) +
  labs(
    title = paste0("Sarcosine-degradation-associated species\nacross 4 WGS cohorts\n",
                   "(Spearman rho ",
                   ifelse(DIR == "negative", "< -0.3", ifelse(DIR == "both", "|rho| > 0.3", "> 0.3")),
                   " & p_adj < 0.05)"),
    subtitle = paste0("Degradation_sum = sum of available Degradation KOs per cohort\n",
                      "Total unique species: ", length(union_deg),
                      "  |  4-way intersection: ", length(deg_intersect4))
  ) +
  theme(plot.title    = element_text(size = 26, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 18, hjust = 0.5, margin = margin(b = 14)),
        plot.margin   = margin(20, 20, 20, 20)) +
  # highlight the 4-cohort central intersection count in bright red
  annotate("text", x = 0, y = -0.7,
           label = as.character(length(deg_intersect4)),
           colour = "red", fontface = "bold", size = 7)
.ggsave(file.path(OUTDIR, "cross_cohort_deg_assoc_venn.png"), venn_plot, width = 11, height = 11, dpi = 300)
cat("Saved: cross_cohort_deg_assoc_venn.png\n")

# 4.2 All correlations (audit trail)
all_corr_df <- bind_rows(cohort_corr_results)
write.csv(all_corr_df, file.path(OUTDIR, "cross_cohort_deg_assoc_correlations.csv"), row.names = FALSE)

# 4.3 Per-cohort significant lists
per_cohort_sig <- bind_rows(lapply(names(cohort_corr_results), function(co) {
  cohort_corr_results[[co]] %>% filter(p_adj < PADJ_CUT & dir_keep(rho)) %>%
    select(Cohort, Species_full, Species, rho, p_adj) %>% arrange(desc(rho))
}))
write.csv(per_cohort_sig, file.path(OUTDIR, "cross_cohort_deg_assoc_per_cohort_lists.csv"), row.names = FALSE)

# 4.4 Membership matrix (wide)
membership_rows <- lapply(union_deg, function(sp) {
  out <- list(Species_full = sp)
  short <- NA
  for (co in names(cohort_corr_results)) {
    df <- cohort_corr_results[[co]]
    if (sp %in% df$Species_full) { short <- df$Species[df$Species_full == sp][1]; break }
  }
  out$Species <- short
  for (co in names(cohort_corr_results)) {
    df <- cohort_corr_results[[co]]; idx <- which(df$Species_full == sp)
    if (length(idx) == 1) {
      out[[paste0(co, "_sig")]]   <- (df$p_adj[idx] < PADJ_CUT & dir_keep(df$rho[idx]))
      out[[paste0(co, "_rho")]]   <- df$rho[idx]
      out[[paste0(co, "_p_adj")]] <- df$p_adj[idx]
    } else {
      out[[paste0(co, "_sig")]] <- FALSE
      out[[paste0(co, "_rho")]] <- NA_real_
      out[[paste0(co, "_p_adj")]] <- NA_real_
    }
  }
  out$n_cohorts_significant <- sum(unlist(lapply(names(cohort_corr_results),
                                                 function(co) out[[paste0(co, "_sig")]])))
  out
})
membership_df <- bind_rows(membership_rows)
membership_df$mean_rho <- apply(membership_df[, paste0(names(cohort_corr_results), "_rho")], 1,
                                function(x) mean(x, na.rm = TRUE))
membership_df <- membership_df %>% arrange(desc(n_cohorts_significant), desc(mean_rho)) %>% select(-mean_rho)
write.csv(membership_df, file.path(OUTDIR, "cross_cohort_deg_assoc_membership_matrix.csv"), row.names = FALSE)

# 4.5 4-way intersection
if (length(deg_intersect4) == 0) {
  empty <- data.frame(Species_full = character(), Species = character(),
                      PRJEB10878_rho = numeric(), PRJEB10878_p_adj = numeric(),
                      PRJEB27928_rho = numeric(), PRJEB27928_p_adj = numeric(),
                      PRJEB6070_rho = numeric(),  PRJEB6070_p_adj = numeric(),
                      PRJNA429097_rho = numeric(), PRJNA429097_p_adj = numeric(),
                      mean_rho = numeric(), max_p_adj = numeric(), stringsAsFactors = FALSE)
  write.csv(empty, file.path(OUTDIR, "cross_cohort_deg_assoc_intersection_4cohorts.csv"), row.names = FALSE)
  cat("4-way intersection EMPTY -> header-only file written.\n")
} else {
  rows <- lapply(deg_intersect4, function(sp) {
    out <- list(Species_full = sp)
    short <- NA
    for (co in names(cohort_corr_results)) {
      df <- cohort_corr_results[[co]]
      if (sp %in% df$Species_full) { short <- df$Species[df$Species_full == sp][1]; break }
    }
    out$Species <- short; rho_vec <- c(); padj_vec <- c()
    for (co in names(cohort_corr_results)) {
      df <- cohort_corr_results[[co]]; idx <- which(df$Species_full == sp)
      out[[paste0(co, "_rho")]] <- df$rho[idx]; out[[paste0(co, "_p_adj")]] <- df$p_adj[idx]
      rho_vec <- c(rho_vec, df$rho[idx]); padj_vec <- c(padj_vec, df$p_adj[idx])
    }
    out$mean_rho <- mean(rho_vec); out$max_p_adj <- max(padj_vec); out
  })
  inter_df <- bind_rows(rows) %>% arrange(desc(mean_rho))
  write.csv(inter_df, file.path(OUTDIR, "cross_cohort_deg_assoc_intersection_4cohorts.csv"), row.names = FALSE)
  cat("4-way intersection species:", nrow(inter_df), "\n")
}

# 4.6 Summary
summary_df <- bind_rows(cohort_summary)
summary_df <- bind_rows(summary_df, data.frame(
  Cohort = "Intersection_4cohorts", n_samples_in_correlation = NA_integer_,
  n_species_tested = NA_integer_, n_significant_padj_only = NA_integer_,
  n_degradation_associated_strict = length(deg_intersect4),
  deg_kos_used = NA_character_, direction = DIR, stringsAsFactors = FALSE))
write.csv(summary_df, file.path(OUTDIR, "cross_cohort_deg_assoc_summary.csv"), row.names = FALSE)

cat("\n=== DONE ===\n")
print(summary_df)
cat("\nR version:", R.version.string, "| ggvenn", as.character(packageVersion("ggvenn")), "\n")
