###############################################################################
# [PARTIALLY SUPERSEDED 2026-05-21]  Section 15 only.
#
# Section 15 (pathway balance -> sarcosine_prod_vs_deg_pooled.png) is superseded
# by Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/regenerate_plots.R,
# which re-renders that figure at dissertation font sizes (6.5-in canvas,
# 300 dpi). If you re-run THIS script, skip Section 15 or it will overwrite the
# pooled pathway figure with the older smaller-font version. All other sections
# of this script (alpha/beta diversity, volcano, heatmap, enzyme plots) remain
# current.
###############################################################################
###############################################################################
# Regenerate Diversity Plots (Integrated 4 CRC Cohorts)
#
# Purpose: Re-render the alpha + beta diversity figures with:
#   1. Corrected group order: Healthy (left / first), Cancer (right / second)
#   2. New colorblind-friendly palette (blue / orange) shared across all plots
#   3. Additional alpha-diversity bar plot (mean +/- SEM)
#   4. Larger, more readable fonts across all plots
#
# IMPORTANT:
#   - Uses the SAME data-loading logic as integrated_analysis_pooled.R so
#     results are bit-for-bit identical. The only changes are visualization.
#   - Statistical test results (Wilcoxon, PERMANOVA) are symmetric in Group;
#     p-values are unchanged by flipping factor levels. set.seed(42) is used
#     before adonis2 so PERMANOVA p-values are reproducible with the
#     original run.
#   - Data files are read from the current (Sarcosine_분석_모음) path.
###############################################################################

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(jsonlite)  # for KO_relative_abundance.tsv (JSON format)
  library(pheatmap)  # for Section 16 differential species heatmap
  library(ggrepel)   # for Section 17 volcano plot labels
})

# Helper: parse species name out of a full taxonomy string (from main script)
get_species_name <- function(taxa_str) {
  parts <- strsplit(taxa_str, "\\|")[[1]]
  sp <- parts[grep("^s__", parts)]
  if (length(sp) > 0) return(gsub("^s__", "", sp[1]))
  return(taxa_str)
}

# --- Sarcosine KO definitions (verbatim from integrated_analysis_pooled.R L52-73) ---
sarcosine_kos <- data.frame(
  KO = c("K00301", "K00302", "K00303", "K00304", "K00305", "K00306",
         "K00315", "K00552", "K08688"),
  Enzyme = c("Sarcosine oxidase (monomeric)",
             "Sarcosine oxidase subunit alpha (soxA)",
             "Sarcosine oxidase subunit beta (soxB)",
             "Sarcosine oxidase subunit delta (soxD)",
             "Sarcosine oxidase subunit gamma (soxG)",
             "Sarcosine oxidase / L-pipecolate oxidase (PIPOX)",
             "Dimethylglycine dehydrogenase (DMGDH)",
             "Glycine N-methyltransferase (GNMT)",
             "Creatinase"),
  EC = c("EC:1.5.3.1", "EC:1.5.3.1/1.5.3.24", "EC:1.5.3.1/1.5.3.24",
         "EC:1.5.3.1/1.5.3.24", "EC:1.5.3.1/1.5.3.24", "EC:1.5.3.1/1.5.3.7",
         "EC:1.5.8.4", "EC:2.1.1.20", "EC:3.5.3.3"),
  Role = c("Degradation", "Degradation", "Degradation", "Degradation",
           "Degradation", "Degradation",
           "Production", "Production", "Production"),
  Short = c("SOX (mono)", "soxA", "soxB", "soxD", "soxG", "PIPOX",
            "DMGDH", "GNMT", "Creatinase"),
  stringsAsFactors = FALSE
)

PREV_THRESHOLD <- 10  # percent (same as main script)

# Use quartz backend for PNG on macOS (cairo not available; agg fails)
options(bitmapType = "quartz")

# Helper: save a ggplot using grDevices::png with quartz backend,
# avoiding ggsave's default agg_png device which fails on this system.
save_png_quartz <- function(plot_obj, filename, width_in, height_in, dpi = 200) {
  grDevices::png(filename = filename,
                 width    = width_in,
                 height   = height_in,
                 units    = "in",
                 res      = dpi,
                 type     = "quartz")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot_obj)
  invisible(NULL)
}

# Current (moved) location of the cohort files
BASE_DIR <- "Healthy_vs_Cancer_4_CRC_cohorts_integrated"   # relative to the analysis-folder root (working dir)
setwd(BASE_DIR)

COHORT_DIRS <- list(
  PRJEB6070   = file.path(BASE_DIR, "PRJEB6070_CRC_AdenomatousPolyps"),
  PRJNA429097 = file.path(BASE_DIR, "PRJNA429097_CRC"),
  PRJEB10878  = file.path(BASE_DIR, "PRJEB10878_CRC"),
  PRJEB27928  = file.path(BASE_DIR, "PRJEB27928_CRC")
)

format_pval <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.0001) return("p<0.0001")
  if (p < 0.001)  return(sprintf("p=%.4f", p))
  return(sprintf("p=%.3f", p))
}

###############################################################################
# Section 1: Load & merge metadata (same logic as main script)
###############################################################################
cat("=== Regenerating alpha diversity plots ===\n\n")
cat("--- Loading metadata ---\n")

meta_pooled <- data.frame()
for (cohort_id in names(COHORT_DIRS)) {
  cohort_path <- COHORT_DIRS[[cohort_id]]
  meta_file <- list.files(cohort_path, pattern = "^selected_project_.*\\.txt$",
                          full.names = TRUE)
  meta_file <- sort(meta_file, decreasing = TRUE)[1]
  cat("  Loading:", cohort_id, "from", basename(meta_file), "\n")

  meta_lines <- readLines(meta_file)
  header <- strsplit(meta_lines[2], "\t")[[1]]
  data_lines <- meta_lines[3:length(meta_lines)]
  data_lines <- data_lines[data_lines != ""]
  data_list <- strsplit(data_lines, "\t")
  data_list <- lapply(data_list, function(x) {
    if (length(x) >= length(header)) x[1:length(header)]
    else c(x, rep(NA, length(header) - length(x)))
  })
  meta <- as.data.frame(do.call(rbind, data_list), stringsAsFactors = FALSE)
  colnames(meta) <- gsub(" ", ".", header)

  all_pheno <- unique(meta$Phenotype.name)
  cancer_pheno <- all_pheno[!all_pheno %in% c("Health", "Adenomatous Polyps")]
  if (length(cancer_pheno) == 0) next
  cancer_label <- cancer_pheno[1]

  meta_hc <- meta %>%
    filter(Phenotype.name %in% c("Health", cancer_label)) %>%
    mutate(Group = ifelse(Phenotype.name == "Health", "Healthy", "Cancer"))

  if ("Assay.type" %in% colnames(meta_hc)) {
    meta_hc <- meta_hc %>% filter(Assay.type == "WGS")
  }

  meta_hc$Cohort <- cohort_id
  meta_pooled <- rbind(meta_pooled, meta_hc)
  cat("    Healthy:", sum(meta_hc$Group == "Healthy"),
      ", Cancer:", sum(meta_hc$Group == "Cancer"), "\n")
}

cat("Total pooled samples (before matching to bacteria):", nrow(meta_pooled), "\n\n")

# Verification invariants
stopifnot(all(meta_pooled$Assay.type == "WGS" | is.na(meta_pooled$Assay.type)))
stopifnot(!any(meta_pooled$Phenotype.name == "Adenomatous Polyps"))
cat("[VERIFIED] No 16S, no Adenomatous Polyps.\n\n")

###############################################################################
# Section 2: Load bacteria data (same logic as main script)
###############################################################################
cat("--- Loading bacteria data ---\n")
bact_all <- data.frame()

for (cohort_id in names(COHORT_DIRS)) {
  cohort_path <- COHORT_DIRS[[cohort_id]]
  bact_file <- list.files(cohort_path, pattern = "^Bacteria_.*\\.txt$",
                          full.names = TRUE)[1]
  cat("  Loading:", cohort_id, "from", basename(bact_file), "\n")

  bact <- read.delim(bact_file, stringsAsFactors = FALSE)
  colnames(bact) <- trimws(colnames(bact))

  cohort_runs <- meta_pooled$Run.ID[meta_pooled$Cohort == cohort_id]
  bact_species <- bact %>%
    filter(grepl("\\|s__", Taxa) & !grepl("\\|t__", Taxa)) %>%
    filter(Run.ID %in% cohort_runs)

  bact_all <- rbind(bact_all, bact_species)
  cat("    Rows:", nrow(bact_species), "\n")
  rm(bact, bact_species); gc(verbose = FALSE)
}

bact_wide <- bact_all %>%
  select(Taxa, Run.ID, Abundance) %>%
  pivot_wider(names_from = Taxa, values_from = Abundance, values_fill = 0)

bact_mat <- as.data.frame(bact_wide)
rownames(bact_mat) <- bact_mat$Run.ID
bact_mat$Run.ID <- NULL
bact_mat <- as.matrix(bact_mat)

sample_order <- intersect(rownames(bact_mat), meta_pooled$Run.ID)
bact_mat <- bact_mat[sample_order, , drop = FALSE]
meta_matched <- meta_pooled %>% filter(Run.ID %in% sample_order)
meta_matched <- meta_matched[match(sample_order, meta_matched$Run.ID), ]

cat("  Bacteria matrix:", nrow(bact_mat), "samples x",
    ncol(bact_mat), "species\n\n")
rm(bact_all, bact_wide); gc(verbose = FALSE)

###############################################################################
# Section 3: Compute alpha diversity (identical to main script)
###############################################################################
cat("--- Computing alpha diversity ---\n")
alpha_div <- data.frame(
  Run.ID   = rownames(bact_mat),
  Shannon  = diversity(bact_mat, index = "shannon"),
  Simpson  = diversity(bact_mat, index = "simpson"),
  Richness = specnumber(bact_mat),
  stringsAsFactors = FALSE
)
alpha_div <- merge(alpha_div, meta_matched[, c("Run.ID", "Group")], by = "Run.ID")

# FIX 1: Force group order so Healthy appears first (left)
alpha_div$Group <- factor(alpha_div$Group, levels = c("Healthy", "Cancer"))

# Wilcoxon tests (symmetric: unaffected by factor level order)
shannon_test  <- wilcox.test(Shannon  ~ Group, data = alpha_div)
simpson_test  <- wilcox.test(Simpson  ~ Group, data = alpha_div)
richness_test <- wilcox.test(Richness ~ Group, data = alpha_div)

cat("Shannon  p:", format_pval(shannon_test$p.value),  "\n")
cat("Simpson  p:", format_pval(simpson_test$p.value),  "\n")
cat("Richness p:", format_pval(richness_test$p.value), "\n")

# Group size sanity check
cat("\nGroup sizes:\n"); print(table(alpha_div$Group))

###############################################################################
# Section 4: Build plot-ready long dataframe + p-value labels
###############################################################################
# Explicit facet order: Richness -> Shannon -> Simpson (matches original layout)
metric_order <- c("Richness", "Shannon", "Simpson")

alpha_long <- alpha_div %>%
  pivot_longer(cols = c(Shannon, Simpson, Richness),
               names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = metric_order))

pval_labels <- data.frame(
  Metric = factor(metric_order, levels = metric_order),
  label  = c(format_pval(richness_test$p.value),
             format_pval(shannon_test$p.value),
             format_pval(simpson_test$p.value)),
  stringsAsFactors = FALSE
)
ypos <- alpha_long %>%
  group_by(Metric) %>%
  summarise(ymax = max(Value, na.rm = TRUE), .groups = "drop")
pval_labels <- merge(pval_labels, ypos, by = "Metric")

###############################################################################
# Section 5: Colors  (unified group palette: green Healthy / red Cancer)
###############################################################################
group_colors <- c("Healthy" = "#1B7837",  # green (favorable)
                  "Cancer"  = "#B2182B")  # red (unfavorable)

###############################################################################
# Shared theme — larger fonts for readability (user request)
###############################################################################
bigger_theme <- theme_bw() +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(size = 24, face = "bold"),
    plot.subtitle    = element_text(size = 16),
    strip.text       = element_text(size = 17, face = "bold"),
    axis.title.x     = element_text(size = 16),
    axis.title.y     = element_text(size = 16),
    axis.text.x      = element_text(size = 16),
    axis.text.y      = element_text(size = 14),
    legend.title     = element_text(size = 16),
    legend.text      = element_text(size = 14)
  )

PVAL_LABEL_SIZE <- 5.9  # ggplot geom_text size; ~12pt

###############################################################################
# Section 6: FIX 3a — Boxplot (same design as before, new order + colors)
###############################################################################
p_box <- ggplot(alpha_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.15) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  geom_text(data = pval_labels,
            aes(x = 1.5, y = ymax * 1.1, label = label),
            inherit.aes = FALSE, size = PVAL_LABEL_SIZE) +
  labs(title = paste0("Alpha diversity (4 CRC cohorts, n=",
                       nrow(bact_mat), ")"),
       subtitle = "Healthy vs Colorectal Neoplasms (Pooled)",
       x = "", y = "Value") +
  bigger_theme

save_png_quartz(p_box,
                "results_integrated/bacteria/alpha_diversity_pooled.png",
                width_in = 10, height_in = 5, dpi = 200)
cat("\n[SAVED] results_integrated/bacteria/alpha_diversity_pooled.png (boxplot)\n")

###############################################################################
# Section 7: FIX 3b — Bar plot (mean +/- SEM)
###############################################################################
bar_summary <- alpha_long %>%
  group_by(Metric, Group) %>%
  summarise(
    mean_val = mean(Value, na.rm = TRUE),
    sd_val   = sd(Value,   na.rm = TRUE),
    n_val    = sum(!is.na(Value)),
    sem      = sd_val / sqrt(n_val),
    .groups  = "drop"
  )
cat("\nBar plot summary (mean +/- SEM):\n")
print(bar_summary)

# p-value y-position for bars: top of error bar + small padding
bar_pval_pos <- bar_summary %>%
  group_by(Metric) %>%
  summarise(ymax = max(mean_val + sem, na.rm = TRUE), .groups = "drop")
bar_pval_labels <- merge(
  data.frame(Metric = factor(metric_order, levels = metric_order),
             label  = c(format_pval(richness_test$p.value),
                        format_pval(shannon_test$p.value),
                        format_pval(simpson_test$p.value)),
             stringsAsFactors = FALSE),
  bar_pval_pos, by = "Metric")

p_bar <- ggplot(bar_summary,
                aes(x = Group, y = mean_val, fill = Group)) +
  geom_col(width = 0.6, alpha = 0.85, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = mean_val - sem, ymax = mean_val + sem),
                width = 0.15, linewidth = 0.4) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(values = group_colors) +
  geom_text(data = bar_pval_labels,
            aes(x = 1.5, y = ymax * 1.1, label = label),
            inherit.aes = FALSE, size = PVAL_LABEL_SIZE) +
  labs(title = paste0("Alpha diversity (4 CRC cohorts, n=",
                       nrow(bact_mat), ")"),
       subtitle = "Healthy vs Colorectal Neoplasms (Pooled) - mean +/- SEM",
       x = "", y = "Mean value") +
  bigger_theme

save_png_quartz(p_bar,
                "results_integrated/bacteria/alpha_diversity_pooled_barplot.png",
                width_in = 10, height_in = 5, dpi = 200)
cat("[SAVED] results_integrated/bacteria/alpha_diversity_pooled_barplot.png (barplot)\n")

###############################################################################
# Section 8: Beta Diversity PCoA (Bray-Curtis)
#   - Identical computation to integrated_analysis_pooled.R lines 309-352
#   - Changes: group color palette -> blue/orange; factor levels ->
#     Healthy first (legend consistency); larger fonts.
###############################################################################
cat("\n--- Beta Diversity (PCoA, Bray-Curtis) ---\n")

# Force Healthy-first factor levels on meta_matched BEFORE building pcoa_df
meta_matched$Group <- factor(meta_matched$Group, levels = c("Healthy", "Cancer"))

bc_dist  <- vegan::vegdist(bact_mat, method = "bray")
pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)
# eigenvalue percent — match main script formula exactly
eig_pct  <- round(pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]) * 100, 1)

pcoa_df <- data.frame(
  PC1    = pcoa_res$points[, 1],
  PC2    = pcoa_res$points[, 2],
  Group  = meta_matched$Group,    # already factor with Healthy first
  Cohort = meta_matched$Cohort,
  stringsAsFactors = FALSE
)

# Reproducibility: same seed, same permutation count as original run
set.seed(42)
perm_simple <- vegan::adonis2(bc_dist ~ Group,          data = meta_matched, permutations = 999)
perm_cohort <- vegan::adonis2(bc_dist ~ Group + Cohort, data = meta_matched, permutations = 999)

perm_p1   <- perm_simple$`Pr(>F)`[1]
perm_r2_1 <- round(perm_simple$R2[1], 4)
perm_p2   <- perm_cohort$`Pr(>F)`[1]
perm_r2_2 <- round(perm_cohort$R2[1], 4)

cat("PERMANOVA (Group only)     : R2=", perm_r2_1, ", p=", format_pval(perm_p1), "\n")
cat("PERMANOVA (Group + Cohort) : Group R2=", perm_r2_2, ", p=", format_pval(perm_p2), "\n")

# Same cohort shape mapping as main script
cohort_shapes <- c("PRJEB6070" = 16, "PRJNA429097" = 17,
                   "PRJEB10878" = 15, "PRJEB27928" = 18)

p_beta <- ggplot(pcoa_df, aes(x = PC1, y = PC2,
                               color = Group, shape = Cohort)) +
  geom_point(size = 1.5, alpha = 0.5) +
  stat_ellipse(aes(group = Group), level = 0.95,
               linetype = 2, linewidth = 0.7) +
  scale_color_manual(values = group_colors) +
  scale_shape_manual(values = cohort_shapes) +
  labs(title = paste0("Bray-Curtis PCoA (4 CRC cohorts, n=",
                       nrow(bact_mat), ")"),
       x = paste0("PCoA1 (", eig_pct[1], "%)"),
       y = paste0("PCoA2 (", eig_pct[2], "%)"),
       subtitle = paste0("PERMANOVA (Group): R2=", perm_r2_1, ", ",
                         format_pval(perm_p1),
                         "  |  PERMANOVA (Group+Cohort): Group R2=",
                         perm_r2_2, ", ", format_pval(perm_p2))) +
  bigger_theme +
  theme(plot.subtitle = element_text(size = 13))  # slightly smaller to fit long text
# Override default legend dot/ellipse visibility: make color legend dots bigger
# so the color is actually readable in the legend (points in plot are size 1.5)
p_beta <- p_beta +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)),
         shape = guide_legend(override.aes = list(size = 3, alpha = 1)))

save_png_quartz(p_beta,
                "results_integrated/bacteria/beta_diversity_pcoa_pooled.png",
                width_in = 11.5, height_in = 7, dpi = 200)
cat("[SAVED] results_integrated/bacteria/beta_diversity_pcoa_pooled.png\n")

###############################################################################
# Section 9: Load KO data (required for sarcosine plots)
#   Verbatim logic from integrated_analysis_pooled.R lines 209-249.
###############################################################################
cat("\n--- Loading KO data for sarcosine plots ---\n")
ko_all <- data.frame()
for (cohort_id in names(COHORT_DIRS)) {
  cohort_path <- COHORT_DIRS[[cohort_id]]
  ko_file <- file.path(cohort_path, "KO_relative_abundance.tsv")
  cat("  Loading:", cohort_id, "...\n")

  ko_raw <- jsonlite::fromJSON(ko_file)
  cohort_runs <- meta_pooled$Run.ID[meta_pooled$Cohort == cohort_id]
  ko_filtered <- ko_raw %>%
    filter(run_id %in% cohort_runs) %>%
    select(ko, run_id, abundance)

  ko_all <- rbind(ko_all, ko_filtered)
  cat("    Rows:", nrow(ko_filtered), "\n")
  rm(ko_raw, ko_filtered); gc(verbose = FALSE)
}

cat("  Total KO rows:", nrow(ko_all), "\n")
ko_wide <- ko_all %>%
  pivot_wider(names_from = ko, values_from = abundance, values_fill = 0)

ko_mat <- as.data.frame(ko_wide)
rownames(ko_mat) <- ko_mat$run_id
ko_mat$run_id <- NULL
ko_mat <- as.matrix(ko_mat)

ko_samples <- intersect(rownames(ko_mat), meta_pooled$Run.ID)
ko_mat <- ko_mat[ko_samples, , drop = FALSE]
meta_ko <- meta_pooled %>% filter(Run.ID %in% ko_samples)
meta_ko <- meta_ko[match(ko_samples, meta_ko$Run.ID), ]

cat("  KO matrix:", nrow(ko_mat), "samples x", ncol(ko_mat), "KOs\n\n")
rm(ko_all, ko_wide); gc(verbose = FALSE)

# FIX: force Group factor levels on meta_ko as well (Healthy first)
meta_ko$Group <- factor(meta_ko$Group, levels = c("Healthy", "Cancer"))

h_idx <- which(meta_ko$Group == "Healthy")
c_idx <- which(meta_ko$Group == "Cancer")
cat("KO sample counts: Healthy =", length(h_idx),
    ", Cancer =", length(c_idx), "\n")

###############################################################################
# Section 10: Sarcosine KO computation (verbatim from main script)
###############################################################################
cat("\n--- Computing sarcosine KO statistics ---\n")

sarc_available <- intersect(sarcosine_kos$KO, colnames(ko_mat))
sarc_info <- sarcosine_kos %>% filter(KO %in% sarc_available)
sarc_mat <- ko_mat[, sarc_available, drop = FALSE]
cat("Sarcosine KOs available:", length(sarc_available), "/",
    nrow(sarcosine_kos), "\n")

sarc_results <- data.frame()
for (i in 1:nrow(sarc_info)) {
  ko_id <- sarc_info$KO[i]
  h_vals <- sarc_mat[h_idx, ko_id]
  c_vals <- sarc_mat[c_idx, ko_id]
  test_res <- tryCatch(wilcox.test(h_vals, c_vals),
                       error = function(e) list(p.value = NA))
  sarc_results <- rbind(sarc_results, data.frame(
    KO = ko_id, Enzyme = sarc_info$Enzyme[i],
    EC = sarc_info$EC[i], Role = sarc_info$Role[i],
    Mean_Healthy = mean(h_vals), Mean_Cancer = mean(c_vals),
    Median_Healthy = median(h_vals), Median_Cancer = median(c_vals),
    p_value = test_res$p.value, stringsAsFactors = FALSE
  ))
}
sarc_results$p_adj <- p.adjust(sarc_results$p_value, method = "BH")

# Prevalence
sarc_results$Prev_Overall <- NA
sarc_results$Prev_Healthy <- NA
sarc_results$Prev_Cancer  <- NA
for (i in 1:nrow(sarc_results)) {
  ko_id <- sarc_results$KO[i]
  vals <- sarc_mat[, ko_id]
  sarc_results$Prev_Overall[i] <- round(sum(vals > 0) / length(vals) * 100, 1)
  sarc_results$Prev_Healthy[i] <- round(sum(vals[h_idx] > 0) / length(h_idx) * 100, 1)
  sarc_results$Prev_Cancer[i]  <- round(sum(vals[c_idx] > 0) / length(c_idx) * 100, 1)
}
sarc_results$Pass_Prevalence <- sarc_results$Prev_Overall >= PREV_THRESHOLD

cat("\nSarcosine KO results:\n")
print(sarc_results[, c("KO", "Enzyme", "Role", "Mean_Healthy",
                       "Mean_Cancer", "p_value", "Prev_Overall",
                       "Pass_Prevalence")])

sarc_sig       <- sarc_results[!is.na(sarc_results$p_value) &
                                  sarc_results$p_value < 0.05, ]
sarc_filtered  <- sarc_results[sarc_results$Pass_Prevalence, ]
sarc_filt_sig  <- sarc_filtered[!is.na(sarc_filtered$p_value) &
                                   sarc_filtered$p_value < 0.05, ]

cat("\nSignificant sarcosine KOs (p<0.05):",
    nrow(sarc_sig), "(all),",
    nrow(sarc_filt_sig), "(prevalence filtered)\n")

###############################################################################
# Helper: build a tweaked theme for sarcosine plots
#   Different from bigger_theme because:
#   - axis.text.x for bar plots needs angle=45, size ~10
#   - strip.text for box plots needs a bit smaller size (~10) due to
#     long multi-line facet labels
###############################################################################
sarc_bar_theme <- theme_bw() +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(size = 24, face = "bold"),
    plot.subtitle    = element_text(size = 14),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y      = element_text(size = 14),
    axis.title.x     = element_text(size = 16),
    axis.title.y     = element_text(size = 16),
    legend.title     = element_text(size = 16),
    legend.text      = element_text(size = 14),
    plot.margin      = margin(t = 10, r = 10, b = 20, l = 20)
  )

sarc_box_theme <- theme_bw() +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(size = 24, face = "bold"),
    plot.subtitle    = element_text(size = 14),
    strip.text       = element_text(size = 13, face = "bold"),
    axis.text.x      = element_text(size = 16),
    axis.text.y      = element_text(size = 14),
    axis.title.x     = element_text(size = 16),
    axis.title.y     = element_text(size = 16),
    legend.title     = element_text(size = 16),
    legend.text      = element_text(size = 14)
  )

sarc_pathway_theme <- theme_bw() +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(size = 24, face = "bold"),
    plot.subtitle    = element_text(size = 13),
    strip.text       = element_text(size = 16, face = "bold"),
    axis.text.x      = element_text(size = 16),
    axis.text.y      = element_text(size = 14),
    axis.title.x     = element_text(size = 16),
    axis.title.y     = element_text(size = 16),
    legend.title     = element_text(size = 16),
    legend.text      = element_text(size = 14)
  )

# P-value label sizes (increased proportionally from originals)
SARC_PVAL_BAR_SIZE  <- 5.2   # was 3
SARC_PVAL_BOX_SIZE  <- 5.2   # was 3
SARC_PVAL_PATH_SIZE <- 5.9   # was 3.5

###############################################################################
# Section 11: Sarcosine plot #1 — enzyme bar plot (unfiltered, significant KOs)
#   Main-script reference: lines 549-568
#   Layout preserved: reorder(Label, -Abundance), angled x labels,
#     width = min(max(8, nrow(sig)*1.5+2), 16), height = 7
###############################################################################
if (nrow(sarc_sig) > 0) {
  cat("\n--- Plot 1: enzyme bar plot (unfiltered) ---\n")

  sarc_plot_data <- sarc_sig %>%
    select(KO, Enzyme, EC, Role, Mean_Healthy, Mean_Cancer, p_value) %>%
    pivot_longer(cols = c(Mean_Healthy, Mean_Cancer),
                 names_to = "Group", values_to = "Abundance") %>%
    mutate(Group = gsub("Mean_", "", Group),
           Group = factor(Group, levels = c("Healthy", "Cancer")),
           Label = paste0(Enzyme, "\n(", Role, ") [", EC, "]"))

  sarc_pval_labels <- sarc_sig %>%
    mutate(Label = paste0(Enzyme, "\n(", Role, ") [", EC, "]"),
           pval_text = sapply(p_value, format_pval))
  ypos_data <- sarc_plot_data %>%
    group_by(Label) %>%
    summarise(ymax = max(Abundance), .groups = "drop")
  sarc_pval_labels <- merge(sarc_pval_labels, ypos_data, by = "Label")

  p_sarc_bar <- ggplot(sarc_plot_data,
                       aes(x = reorder(Label, -Abundance),
                           y = Abundance, fill = Group)) +
    geom_bar(stat = "identity",
             position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(data = sarc_pval_labels,
              aes(x = Label, y = ymax * 1.15, label = pval_text),
              inherit.aes = FALSE, size = SARC_PVAL_BAR_SIZE,
              fontface = "italic") +
    scale_fill_manual(values = group_colors) +
    labs(title = paste0("4 CRC cohorts: sarcosine enzyme gene abundance (n=",
                        nrow(ko_mat), ")"),
         subtitle = "Healthy vs cancer (significant only, p<0.05)",
         x = "", y = "Mean Relative Abundance") +
    sarc_bar_theme

  fig_w <- max(8, nrow(sarc_sig) * 1.5 + 2)
  save_png_quartz(p_sarc_bar,
                  "results_integrated/sarcosine/sarcosine_enzyme_barplot_pooled.png",
                  width_in = min(fig_w, 16), height_in = 7, dpi = 200)
  cat("[SAVED] sarcosine_enzyme_barplot_pooled.png (", nrow(sarc_sig), "KOs)\n")
}

###############################################################################
# Section 12: Sarcosine plot #2 — enzyme box plot (unfiltered, significant KOs)
#   Main-script reference: lines 589-605
#   Layout preserved: facet_wrap(~Label, scales="free_y", ncol=3),
#     width = 12, height = max(4, ceiling(nsig/3) * 3.5)
###############################################################################
if (nrow(sarc_sig) > 0) {
  cat("\n--- Plot 2: enzyme box plot (unfiltered) ---\n")

  sarc_box_data <- data.frame()
  for (i in 1:nrow(sarc_sig)) {
    ko_id <- sarc_sig$KO[i]
    tmp <- data.frame(
      Abundance = sarc_mat[, ko_id],
      Group     = meta_ko$Group,     # already factor with Healthy first
      Label     = paste0(sarc_sig$KO[i], ": ", sarc_sig$Enzyme[i],
                         "\n(", sarc_sig$Role[i], ") [", sarc_sig$EC[i], "]"),
      stringsAsFactors = FALSE)
    sarc_box_data <- rbind(sarc_box_data, tmp)
  }

  box_pval <- sarc_sig %>%
    mutate(Label = paste0(KO, ": ", Enzyme, "\n(", Role, ") [", EC, "]"),
           pval_text = sapply(p_value, format_pval))
  box_ypos <- sarc_box_data %>%
    group_by(Label) %>%
    summarise(ymax = max(Abundance, na.rm = TRUE), .groups = "drop")
  box_pval <- merge(box_pval, box_ypos, by = "Label")

  p_sarc_box <- ggplot(sarc_box_data,
                       aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.75) +
    geom_jitter(width = 0.15, size = 0.3, alpha = 0.15) +
    facet_wrap(~Label, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = group_colors) +
    geom_text(data = box_pval,
              aes(x = 1.5, y = ymax * 1.15, label = pval_text),
              inherit.aes = FALSE, size = SARC_PVAL_BOX_SIZE) +
    labs(title = paste0("Sarcosine Enzyme KOs (n=",
                        nrow(ko_mat), ")"),
         subtitle = "Healthy vs Colorectal Neoplasms (p<0.05)",
         x = "", y = "Relative Abundance") +
    sarc_box_theme

  box_nrow <- ceiling(nrow(sarc_sig) / 3)
  save_png_quartz(p_sarc_box,
                  "results_integrated/sarcosine/sarcosine_enzyme_boxplot_pooled.png",
                  width_in = 17, height_in = max(4, box_nrow * 3.5),
                  dpi = 200)
  cat("[SAVED] sarcosine_enzyme_boxplot_pooled.png (", nrow(sarc_sig), "KOs,",
      box_nrow, "facet rows)\n")
}

###############################################################################
# Section 13: Sarcosine plot #3 — enzyme bar plot (prevalence-filtered)
#   Main-script reference: lines 654-673
###############################################################################
if (nrow(sarc_filt_sig) > 0) {
  cat("\n--- Plot 3: enzyme bar plot (prevalence filtered) ---\n")

  filt_plot_data <- sarc_filt_sig %>%
    select(KO, Enzyme, EC, Role, Mean_Healthy, Mean_Cancer, p_value) %>%
    pivot_longer(cols = c(Mean_Healthy, Mean_Cancer),
                 names_to = "Group", values_to = "Abundance") %>%
    mutate(Group = gsub("Mean_", "", Group),
           Group = factor(Group, levels = c("Healthy", "Cancer")),
           Label = paste0(Enzyme, "\n(", Role, ") [", EC, "]"))

  filt_pval_labels <- sarc_filt_sig %>%
    mutate(Label = paste0(Enzyme, "\n(", Role, ") [", EC, "]"),
           pval_text = sapply(p_value, format_pval))
  filt_ypos <- filt_plot_data %>%
    group_by(Label) %>%
    summarise(ymax = max(Abundance), .groups = "drop")
  filt_pval_labels <- merge(filt_pval_labels, filt_ypos, by = "Label")

  p_filt_bar <- ggplot(filt_plot_data,
                       aes(x = reorder(Label, -Abundance),
                           y = Abundance, fill = Group)) +
    geom_bar(stat = "identity",
             position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(data = filt_pval_labels,
              aes(x = Label, y = ymax * 1.15, label = pval_text),
              inherit.aes = FALSE, size = SARC_PVAL_BAR_SIZE,
              fontface = "italic") +
    scale_fill_manual(values = group_colors) +
    labs(title = paste0("4 CRC cohorts: sarcosine enzymes (filtered, n=",
                        nrow(ko_mat), ")"),
         subtitle = paste0("Prevalence >=", PREV_THRESHOLD,
                           "% | Significant only (p<0.05)"),
         x = "", y = "Mean Relative Abundance") +
    sarc_bar_theme

  filt_fig_w <- max(10, nrow(sarc_filt_sig) * 1.5 + 2)
  save_png_quartz(p_filt_bar,
                  "results_integrated/sarcosine/sarcosine_enzyme_barplot_filtered_pooled.png",
                  width_in = min(filt_fig_w, 16), height_in = 7, dpi = 200)
  cat("[SAVED] sarcosine_enzyme_barplot_filtered_pooled.png (",
      nrow(sarc_filt_sig), "KOs)\n")
}

###############################################################################
# Section 14: Sarcosine plot #4 — enzyme box plot (prevalence-filtered)
#   Main-script reference: lines 697-714
###############################################################################
if (nrow(sarc_filt_sig) > 0) {
  cat("\n--- Plot 4: enzyme box plot (prevalence filtered) ---\n")

  filt_box_data <- data.frame()
  for (i in 1:nrow(sarc_filt_sig)) {
    ko_id <- sarc_filt_sig$KO[i]
    tmp <- data.frame(
      Abundance = sarc_mat[, ko_id],
      Group     = meta_ko$Group,
      Label     = paste0(sarc_filt_sig$KO[i], ": ", sarc_filt_sig$Enzyme[i],
                         "\n(", sarc_filt_sig$Role[i], ") [",
                         sarc_filt_sig$EC[i], "]",
                         "\nprev=", sarc_filt_sig$Prev_Overall[i], "%"),
      stringsAsFactors = FALSE)
    filt_box_data <- rbind(filt_box_data, tmp)
  }

  filt_box_pval <- sarc_filt_sig %>%
    mutate(Label = paste0(KO, ": ", Enzyme,
                          "\n(", Role, ") [", EC, "]",
                          "\nprev=", Prev_Overall, "%"),
           pval_text = sapply(p_value, format_pval))
  filt_box_ypos <- filt_box_data %>%
    group_by(Label) %>%
    summarise(ymax = max(Abundance, na.rm = TRUE), .groups = "drop")
  filt_box_pval <- merge(filt_box_pval, filt_box_ypos, by = "Label")

  p_filt_box <- ggplot(filt_box_data,
                       aes(x = Group, y = Abundance, fill = Group)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.75) +
    geom_jitter(width = 0.15, size = 0.3, alpha = 0.15) +
    facet_wrap(~Label, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = group_colors) +
    geom_text(data = filt_box_pval,
              aes(x = 1.5, y = ymax * 1.15, label = pval_text),
              inherit.aes = FALSE, size = SARC_PVAL_BOX_SIZE) +
    labs(title = paste0("Sarcosine Enzyme KOs (Prevalence Filtered, n=",
                        nrow(ko_mat), ")"),
         subtitle = paste0("Prevalence >=", PREV_THRESHOLD, "% | p<0.05"),
         x = "", y = "Relative Abundance") +
    sarc_box_theme

  filt_box_nrow <- ceiling(nrow(sarc_filt_sig) / 3)
  save_png_quartz(p_filt_box,
                  "results_integrated/sarcosine/sarcosine_enzyme_boxplot_filtered_pooled.png",
                  width_in = 12,
                  height_in = max(4, filt_box_nrow * 3.5),
                  dpi = 200)
  cat("[SAVED] sarcosine_enzyme_boxplot_filtered_pooled.png (",
      nrow(sarc_filt_sig), "KOs)\n")
}

###############################################################################
# Section 15: Sarcosine plot #5 — Pathway Balance (Production vs Degradation)
#   Main-script reference: lines 1081-1166 (Part 5)
#   Layout preserved: facet_wrap(~Metric, scales="free_y", nrow=1),
#     width=12, height=5.5
###############################################################################
cat("\n--- Plot 5: Pathway balance (Prod vs Deg) ---\n")

deg_kos  <- intersect(sarcosine_kos$KO[sarcosine_kos$Role == "Degradation"],
                      colnames(ko_mat))
prod_kos <- intersect(sarcosine_kos$KO[sarcosine_kos$Role == "Production"],
                      colnames(ko_mat))

cat("Degradation KOs available:", length(deg_kos), "(",
    paste(deg_kos, collapse = ","), ")\n")
cat("Production  KOs available:", length(prod_kos), "(",
    paste(prod_kos, collapse = ","), ")\n")

if (length(deg_kos) > 0 && length(prod_kos) > 0) {
  deg_sum  <- rowSums(ko_mat[, deg_kos,  drop = FALSE])
  prod_sum <- rowSums(ko_mat[, prod_kos, drop = FALSE])
  pseudo   <- 1e-8
  log2_ratio <- log2((prod_sum + pseudo) / (deg_sum + pseudo))

  ratio_df <- data.frame(
    Run.ID              = rownames(ko_mat),
    Group               = meta_ko$Group,   # factor
    Cohort              = meta_ko$Cohort,
    Degradation_sum     = deg_sum,
    Production_sum      = prod_sum,
    Log2_Prod_Deg_Ratio = log2_ratio,
    stringsAsFactors = FALSE
  )

  ratio_test <- tryCatch(
    wilcox.test(Log2_Prod_Deg_Ratio ~ Group, data = ratio_df),
    error = function(e) list(p.value = NA))
  deg_test <- tryCatch(
    wilcox.test(Degradation_sum ~ Group, data = ratio_df),
    error = function(e) list(p.value = NA))
  prod_test <- tryCatch(
    wilcox.test(Production_sum ~ Group, data = ratio_df),
    error = function(e) list(p.value = NA))

  cat("Degradation sum p:", format_pval(deg_test$p.value),  "\n")
  cat("Production  sum p:", format_pval(prod_test$p.value), "\n")
  cat("log2(Prod/Deg)  p:", format_pval(ratio_test$p.value),"\n")

  metric_levels_pd <- c(paste0("Degradation\n(sum of ", length(deg_kos), " KOs)"),
                        paste0("Production\n(sum of ", length(prod_kos), " KOs)"),
                        "log2(Prod/Deg)\nRatio")

  plot_data_pd <- ratio_df %>%
    select(Group, Degradation_sum, Production_sum, Log2_Prod_Deg_Ratio) %>%
    pivot_longer(cols = c(Degradation_sum, Production_sum, Log2_Prod_Deg_Ratio),
                 names_to = "Metric", values_to = "Value") %>%
    mutate(Metric = factor(Metric,
                           levels = c("Degradation_sum", "Production_sum",
                                      "Log2_Prod_Deg_Ratio"),
                           labels = metric_levels_pd))

  pval_labels_pd <- data.frame(
    Metric = factor(metric_levels_pd, levels = metric_levels_pd),
    pval_text = c(format_pval(deg_test$p.value),
                  format_pval(prod_test$p.value),
                  format_pval(ratio_test$p.value)),
    stringsAsFactors = FALSE
  )
  ypos_pd <- plot_data_pd %>%
    group_by(Metric) %>%
    summarise(ymax = max(Value, na.rm = TRUE), .groups = "drop")
  pval_labels_pd <- merge(pval_labels_pd, ypos_pd, by = "Metric")

  p_ratio <- ggplot(plot_data_pd,
                    aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.75) +
    geom_jitter(width = 0.15, size = 0.3, alpha = 0.15) +
    facet_wrap(~Metric, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = group_colors) +
    geom_text(data = pval_labels_pd,
              aes(x = 1.5, y = ymax * 1.1, label = pval_text),
              inherit.aes = FALSE, size = SARC_PVAL_PATH_SIZE,
              fontface = "italic") +
    labs(title = paste0("Sarcosine Metabolism Pathway Balance (n=",
                        nrow(ko_mat), ")"),
         subtitle = paste0("Deg KOs: ", paste(deg_kos, collapse = ","),
                           "  |  Prod KOs: ", paste(prod_kos, collapse = ",")),
         x = "", y = "Abundance / Ratio") +
    sarc_pathway_theme

  save_png_quartz(p_ratio,
                  "results_integrated/sarcosine/sarcosine_prod_vs_deg_pooled.png",
                  width_in = 12, height_in = 5.5, dpi = 200)
  cat("[SAVED] sarcosine_prod_vs_deg_pooled.png\n")
} else {
  cat("Cannot compute Prod/Deg ratio (need both pathway KOs available).\n")
}

###############################################################################
# Section 16: Differential species heatmap
#   Main-script reference: lines 354-438 (diff abundance + heatmap)
#
#   Changes vs original:
#     (a) Group annotation colors: green/red -> blue/orange (group_colors)
#     (b) Samples reordered: within-Healthy hclust, then within-Cancer hclust,
#         concatenated [Healthy, Cancer] with a visual gap (gaps_col).
#         cluster_cols disabled (the prior column dendrogram is removed).
#     (c) Larger pheatmap fonts (fontsize 10 -> 13, fontsize_row 7 -> 10)
#     (d) Width 12 -> 14, height nudged up to accommodate bigger fonts
#
#   Kept unchanged:
#     - Prevalence >=10% filter
#     - Top 30 significant species by p_adj ascending
#     - Z-score scaling with cap at +/-3
#     - Blue-white-red cell color gradient (z-score convention)
#     - Row (species) hierarchical clustering
###############################################################################
cat("\n--- Differential abundance (species) for heatmap ---\n")
prevalence <- colSums(bact_mat > 0) / nrow(bact_mat)
bact_filt  <- bact_mat[, prevalence >= 0.1, drop = FALSE]
cat("Species after prevalence filter:", ncol(bact_filt), "\n")

healthy_idx_bact <- which(meta_matched$Group == "Healthy")
cancer_idx_bact  <- which(meta_matched$Group == "Cancer")
cat("Bacteria sample counts: Healthy =", length(healthy_idx_bact),
    ", Cancer =", length(cancer_idx_bact), "\n")

# Per-species Wilcoxon H vs C  (same loop as main script)
diff_results <- data.frame(
  Species_full = colnames(bact_filt),
  Species      = sapply(colnames(bact_filt), get_species_name),
  Mean_Healthy = colMeans(bact_filt[healthy_idx_bact, , drop = FALSE]),
  Mean_Cancer  = colMeans(bact_filt[cancer_idx_bact,  , drop = FALSE]),
  stringsAsFactors = FALSE
)
diff_results$log2FC <- log2((diff_results$Mean_Cancer + 1e-6) /
                            (diff_results$Mean_Healthy + 1e-6))

pvals <- sapply(1:ncol(bact_filt), function(i) {
  tryCatch(wilcox.test(bact_filt[healthy_idx_bact, i],
                       bact_filt[cancer_idx_bact,  i])$p.value,
           error = function(e) NA)
})
diff_results$p_value <- pvals
diff_results$p_adj   <- p.adjust(pvals, method = "BH")
diff_results <- diff_results %>% arrange(p_adj)

sig_species <- diff_results %>% filter(p_adj < 0.05)
cat("Significant species (FDR < 0.05):", nrow(sig_species), "\n")

if (nrow(sig_species) > 1) {
  top_n  <- min(30, nrow(sig_species))
  top_sp <- sig_species$Species_full[1:top_n]

  cat("\nTop", top_n, "species by p_adj (should match original heatmap rows):\n")
  print(sig_species[1:top_n, c("Species", "log2FC", "p_adj")])

  hm_mat    <- bact_filt[, top_sp, drop = FALSE]
  hm_scaled <- scale(hm_mat)
  hm_scaled[hm_scaled > 3]  <-  3
  hm_scaled[hm_scaled < -3] <- -3
  colnames(hm_scaled) <- sapply(colnames(hm_scaled), get_species_name)

  # Sanity check: no NaN in z-scored matrix (would break hclust/dist)
  n_nan <- sum(is.nan(hm_scaled))
  if (n_nan > 0) {
    stop("FATAL: ", n_nan, " NaN values in z-scored heatmap matrix — ",
         "a top species has zero variance across all samples. ",
         "Investigate before proceeding.")
  }

  # --- FIX: Within-group sample clustering ---
  cat("\nClustering samples within each group...\n")
  # dist() computes pairwise distances between rows (samples)
  hm_healthy_block <- hm_scaled[healthy_idx_bact, , drop = FALSE]
  hm_cancer_block  <- hm_scaled[cancer_idx_bact,  , drop = FALSE]
  hc_h_order <- hclust(dist(hm_healthy_block), method = "complete")$order
  hc_c_order <- hclust(dist(hm_cancer_block),  method = "complete")$order

  new_col_order <- c(healthy_idx_bact[hc_h_order],
                     cancer_idx_bact[hc_c_order])
  stopifnot(length(new_col_order) == nrow(hm_scaled))
  stopifnot(!anyDuplicated(new_col_order))

  hm_reordered <- hm_scaled[new_col_order, , drop = FALSE]

  # Column annotation — Healthy/Cancer factor in display order
  annotation_col_df <- data.frame(
    Group = factor(meta_matched$Group[new_col_order],
                   levels = c("Healthy", "Cancer"))
  )
  rownames(annotation_col_df) <- rownames(hm_reordered)

  ann_colors_new <- list(Group = group_colors)  # blue/orange from Section 5

  # Plot dimensions — preserve original aspect but nudge up for larger fonts
  hm_width  <- 14
  hm_height <- max(10, top_n * 0.3 + 3)

  cat("Saving heatmap: width=", hm_width, ", height=", hm_height, "\n")
  grDevices::png(
    filename = "results_integrated/bacteria/heatmap_top_species_pooled.png",
    width    = hm_width,
    height   = hm_height,
    units    = "in",
    res      = 200,
    type     = "quartz"
  )
  tryCatch({
    pheatmap(
      t(hm_reordered),
      annotation_col    = annotation_col_df,
      annotation_colors = ann_colors_new,
      color             = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
      cluster_rows      = TRUE,     # species still clustered
      cluster_cols      = FALSE,    # samples pre-ordered within groups
      gaps_col          = length(healthy_idx_bact),  # visual split
      show_colnames     = FALSE,
      fontsize          = 17,       # up from pheatmap default 10
      fontsize_row      = 13,       # up from main-script 7
      treeheight_row    = 40,
      main = paste0("4 CRC cohorts: top differential species (Z-score, n=",
                    nrow(bact_mat), ")")
    )
  }, finally = grDevices::dev.off())
  cat("[SAVED] results_integrated/bacteria/heatmap_top_species_pooled.png\n")
} else {
  cat("[SKIP] Not enough significant species for heatmap.\n")
}

###############################################################################
# Section 17: Volcano plot (differential species)
#   Main-script reference: lines 386-411
#
#   Changes vs original:
#     (a) Color palette: unified group scheme
#           "Enriched in Cancer"  -> #B2182B (red, unfavorable)
#           "Enriched in Healthy" -> #1B7837 (green, favorable)
#           "Not Significant"     grey60 (unchanged)
#     (b) Larger fonts (title, axis titles, axis text, legend, labels)
#     (c) Point size 1.5 -> 2 for readability
#     (d) geom_text_repel size 2.5 -> 3.5
#
#   Kept unchanged:
#     - Significance thresholds: p_adj < 0.05 AND |log2FC| > 1
#     - Top 15 labels by p_adj ascending
#     - Horizontal/vertical dashed cutoff lines
#     - Axes, plot dimensions (10 x 7 in, 200 dpi)
#
#   This section REUSES diff_results already computed in Section 16, so there
#   is no additional computation beyond the Significance column + labels.
###############################################################################
cat("\n--- Volcano plot (differential species) ---\n")

if (!exists("diff_results")) {
  stop("FATAL: diff_results not found — Section 16 must run before 17.")
}

diff_results$Significance <- "Not Significant"
diff_results$Significance[diff_results$p_adj < 0.05 &
                          diff_results$log2FC > 1]  <- "Enriched in Cancer"
diff_results$Significance[diff_results$p_adj < 0.05 &
                          diff_results$log2FC < -1] <- "Enriched in Healthy"

cat("Volcano significance counts:\n")
print(table(diff_results$Significance))

# Several top species have a BH-adjusted p below double precision (p_adj == 0 ->
# -log10 = Inf), which would plot them off-scale, clamped to the panel top as an
# indistinguishable pile. Cap the plotted y just above the largest finite value
# so those points are visible on-scale; the cap is stated in the caption.
neglog_all <- -log10(diff_results$p_adj)
finite_max <- max(neglog_all[is.finite(neglog_all)], na.rm = TRUE)
y_cap      <- ceiling(finite_max) + 5
diff_results$neglog_plot <- pmin(neglog_all, y_cap)
n_capped   <- sum(!is.finite(neglog_all))
cat(sprintf("Volcano: %d species with p_adj below precision capped at -log10 = %d\n",
            n_capped, y_cap))

# Label only the genuinely significant species (p_adj<0.05 AND |log2FC|>1); the
# grey 'Not Significant' points are NOT labelled. Top 15 by p_adj.
top_to_label <- diff_results %>%
  filter(Significance != "Not Significant") %>%
  arrange(p_adj) %>%
  head(15)

# Colors — Cancer -> red (same as group_colors["Cancer"]),
# Healthy -> green (same as group_colors["Healthy"])
volcano_colors <- c(
  "Enriched in Cancer"  = unname(group_colors["Cancer"]),   # #B2182B
  "Enriched in Healthy" = unname(group_colors["Healthy"]),  # #1B7837
  "Not Significant"     = "grey60"
)

p_volcano <- ggplot(diff_results,
                    aes(x = log2FC, y = neglog_plot, color = Significance)) +
  geom_point(alpha = 0.65, size = 2) +
  scale_color_manual(values = volcano_colors) +
  scale_y_continuous(limits = c(0, y_cap * 1.22),
                     expand = expansion(mult = c(0.02, 0.02))) +  # headroom above the cap for the top labels
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1),
             linetype = "dashed", color = "grey40") +
  geom_text_repel(data = top_to_label, aes(label = Species),
                  size = 4.6, max.overlaps = Inf,
                  min.segment.length = 0.2,
                  box.padding = 0.5, point.padding = 0.3,
                  force = 3, seed = 7,
                  ylim = c(NA, y_cap * 1.20),
                  show.legend = FALSE) +
  labs(title = paste0("Differential species (4 CRC cohorts, n=",
                      nrow(bact_mat), ")"),
       subtitle = "Healthy vs Colorectal Neoplasms (Pooled)",
       x = "log2FC (Cancer/Healthy)",
       y = "-log10(adj. p-value)",
       caption = sprintf("p_adj below machine precision capped at -log10 = %d (n = %d species)",
                         y_cap, n_capped)) +
  theme_bw() +
  theme(
    legend.position  = "bottom",
    plot.title       = element_text(size = 24, face = "bold"),
    plot.subtitle    = element_text(size = 16),
    axis.title.x     = element_text(size = 17),
    axis.title.y     = element_text(size = 17),
    axis.text.x      = element_text(size = 14),
    axis.text.y      = element_text(size = 14),
    legend.title     = element_text(size = 16),
    legend.text      = element_text(size = 14),
    plot.caption     = element_text(size = 12, color = "grey35")
  ) +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

save_png_quartz(p_volcano,
                "results_integrated/bacteria/volcano_species_pooled.png",
                width_in = 11, height_in = 8, dpi = 200)
cat("[SAVED] results_integrated/bacteria/volcano_species_pooled.png\n")

cat("\n=== Done. ===\n")
