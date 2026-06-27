###############################################################################
# Regenerate Species x Sarcosine-Indicator Scatter Plots with Larger Fonts
#
# Purpose: re-render the 5 species-vs-sarcosine-indicator scatter figures at
#   larger, more readable font sizes. The summary bubble plot is rendered at its
#   ORIGINAL styling (enlarging its 40 species labels worsened their overlap).
#
# IMPORTANT - this is a RE-RENDER:
#   - It reproduces the data loading, indicator computation, Spearman
#     correlation, dual filter and species selection of
#     scatter_sarcosine_species.R VERBATIM, so the species shown, the points,
#     regression lines, rho and FDR values are identical to the original
#     figures. The ONLY changes are font sizes and a proportionally larger
#     panel/canvas so the larger text does not clip.
#   - The original scatter_sarcosine_species.R is NOT modified. This script does
#     NOT re-write the two CSV outputs (species_sarcosine_correlation_all.csv,
#     selected_species_for_scatter.csv) - it produces figures only.
#   - Output PNGs overwrite the existing figures at the same paths.
#   - Font sizes follow regenerate_diversity_plots.R for a consistent figure set.
#
# Heavy step: re-loads all 4 cohorts' KO_relative_abundance.tsv (JSON, ~0.8 GB
#   total) - expect a few minutes of run time.
###############################################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(jsonlite)
library(scales)

options(bitmapType = "quartz")

# quartz PNG saver (cairo unavailable on this machine; same approach as
# regenerate_diversity_plots.R - avoids ggsave's default device).
.ggsave <- function(filename, plot, width, height, dpi = 200) {
  grDevices::png(filename = filename, width = width, height = height,
                 units = "in", res = dpi, type = "quartz", bg = "white")
  print(plot)
  grDevices::dev.off()
  invisible(filename)
}

# --- Larger font sizes (consistent with regenerate_diversity_plots.R) --------
# Original sizes shown in comments for reference.
FS_TITLE   <- 26    # plot title                       (was 10; raised to 20 per user request)
FS_SUBT    <- 14    # subtitle                         (was 7)
FS_STRIP   <- 16    # facet strip = species headers    (was 7)
FS_AXIS_T  <- 17    # axis titles                      (was ~11 default)
FS_AXIS_X  <- 14    # axis tick text                   (was ~9 default)
FS_LEG_T   <- 16    # legend title                     (was ~11 default)
FS_LEG_X   <- 14    # legend text                      (was 7)
FS_ANNOT   <- 5.2   # in-panel rho / p annotation      (was 2.5)

# --- Current project root ----------------------------------------------------
BASE_DIR <- "Healthy_vs_Cancer_4_CRC_cohorts_integrated"   # relative to the analysis-folder root (working dir)
setwd(BASE_DIR)

COHORT_DIRS <- list(
  PRJEB6070   = file.path(BASE_DIR, "PRJEB6070_CRC_AdenomatousPolyps"),
  PRJNA429097 = file.path(BASE_DIR, "PRJNA429097_CRC"),
  PRJEB10878  = file.path(BASE_DIR, "PRJEB10878_CRC"),
  PRJEB27928  = file.path(BASE_DIR, "PRJEB27928_CRC")
)

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
  Role = c("Degradation", "Degradation", "Degradation", "Degradation",
           "Degradation", "Degradation",
           "Production", "Production", "Production"),
  stringsAsFactors = FALSE
)

PREV_THRESHOLD <- 10

format_pval <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.0001) return("p<0.0001")
  if (p < 0.001)  return(sprintf("p=%.4f", p))
  return(sprintf("p=%.3f", p))
}

get_species_name <- function(taxa_str) {
  parts <- strsplit(taxa_str, "\\|")[[1]]
  sp <- parts[grep("^s__", parts)]
  if (length(sp) > 0) return(gsub("^s__", "", sp[1]))
  return(taxa_str)
}

dir.create("results_integrated/sarcosine/scatter", recursive = TRUE, showWarnings = FALSE)

cat("=== Re-render: Species x Sarcosine Scatter Plots (larger fonts) ===\n\n")

###############################################################################
# SECTION 1: DATA LOADING  (verbatim from scatter_sarcosine_species.R)
###############################################################################
cat("========== SECTION 1: DATA LOADING ==========\n\n")

# --- 1.1 Metadata ---
cat("--- 1.1 Loading metadata ---\n")
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
    n_before <- nrow(meta_hc)
    meta_hc <- meta_hc %>% filter(Assay.type == "WGS")
    cat("    WGS filter: kept", nrow(meta_hc), "out of", n_before, "\n")
  }

  meta_hc$Cohort <- cohort_id
  meta_pooled <- rbind(meta_pooled, meta_hc)
  cat("    Healthy:", sum(meta_hc$Group == "Healthy"),
      ", Cancer:", sum(meta_hc$Group == "Cancer"), "\n")
}

cat("\nTotal samples:", nrow(meta_pooled), "\n")
stopifnot(all(meta_pooled$Assay.type == "WGS" | is.na(meta_pooled$Assay.type)))
stopifnot(!any(meta_pooled$Phenotype.name == "Adenomatous Polyps"))
cat("[VERIFIED] No 16S, no Adenomatous Polyps.\n\n")

# --- 1.2 Bacteria data ---
cat("--- 1.2 Loading bacteria data ---\n")
bact_all <- data.frame()
for (cohort_id in names(COHORT_DIRS)) {
  cohort_path <- COHORT_DIRS[[cohort_id]]
  bact_file <- list.files(cohort_path, pattern = "^Bacteria_.*\\.txt$",
                          full.names = TRUE)[1]
  cat("  Loading:", cohort_id, "\n")
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

cat("  Pivoting to wide format...\n")
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
cat("  Bacteria matrix:", nrow(bact_mat), "x", ncol(bact_mat), "\n\n")
rm(bact_all, bact_wide); gc(verbose = FALSE)

# --- 1.3 KO data ---
cat("--- 1.3 Loading KO data ---\n")
ko_all <- data.frame()
for (cohort_id in names(COHORT_DIRS)) {
  cohort_path <- COHORT_DIRS[[cohort_id]]
  ko_file <- file.path(cohort_path, "KO_relative_abundance.tsv")
  cat("  Loading:", cohort_id, "\n")
  ko_raw <- fromJSON(ko_file)
  cohort_runs <- meta_pooled$Run.ID[meta_pooled$Cohort == cohort_id]
  ko_filtered <- ko_raw %>%
    filter(run_id %in% cohort_runs) %>%
    select(ko, run_id, abundance)
  ko_all <- rbind(ko_all, ko_filtered)
  cat("    Rows:", nrow(ko_filtered), "\n")
  rm(ko_raw, ko_filtered); gc(verbose = FALSE)
}

cat("  Pivoting to wide format...\n")
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
cat("  KO matrix:", nrow(ko_mat), "x", ncol(ko_mat), "\n\n")
rm(ko_all, ko_wide); gc(verbose = FALSE)

###############################################################################
# SECTION 2: COMPUTE 5 SARCOSINE INDICATORS  (verbatim)
###############################################################################
cat("========== SECTION 2: SARCOSINE INDICATORS ==========\n\n")

common_samples <- intersect(rownames(bact_mat), rownames(ko_mat))
bact_common <- bact_mat[common_samples, , drop = FALSE]
ko_common <- ko_mat[common_samples, , drop = FALSE]
meta_common <- meta_ko %>% filter(Run.ID %in% common_samples)
meta_common <- meta_common[match(common_samples, meta_common$Run.ID), ]
cat("Common samples:", length(common_samples), "\n")

bact_prev <- colSums(bact_common > 0) / nrow(bact_common)
bact_filt <- bact_common[, bact_prev >= (PREV_THRESHOLD / 100), drop = FALSE]
cat("Species after prevalence filter (>=", PREV_THRESHOLD, "%):", ncol(bact_filt), "\n")

ko_K00303 <- ko_common[, "K00303"]
ko_K08688 <- ko_common[, "K08688"]

deg_kos <- intersect(sarcosine_kos$KO[sarcosine_kos$Role == "Degradation"], colnames(ko_common))
prod_kos <- intersect(sarcosine_kos$KO[sarcosine_kos$Role == "Production"], colnames(ko_common))
cat("Degradation KOs:", length(deg_kos), "(", paste(deg_kos, collapse = ", "), ")\n")
cat("Production KOs:", length(prod_kos), "(", paste(prod_kos, collapse = ", "), ")\n")

deg_sum <- rowSums(ko_common[, deg_kos, drop = FALSE])
prod_sum <- rowSums(ko_common[, prod_kos, drop = FALSE])
pseudo <- 1e-8
log2_ratio <- log2((prod_sum + pseudo) / (deg_sum + pseudo))

indicators <- list(
  "K00303_soxB"       = ko_K00303,
  "K08688_Creatinase" = ko_K08688,
  "Degradation_sum"   = deg_sum,
  "Production_sum"    = prod_sum,
  "ProdDeg_ratio"     = log2_ratio
)

indicator_labels <- c(
  "K00303_soxB"       = "K00303 soxB (Degradation)",
  "K08688_Creatinase" = "K08688 Creatinase (Production)",
  "Degradation_sum"   = "Degradation Sum (5 KOs)",
  "Production_sum"    = "Production Sum (2 KOs)",
  "ProdDeg_ratio"     = "log2(Prod/Deg) Ratio"
)

indicator_ylabels <- c(
  "K00303_soxB"       = "K00303 Relative Abundance",
  "K08688_Creatinase" = "K08688 Relative Abundance",
  "Degradation_sum"   = "Degradation Sum (Relative Abundance)",
  "Production_sum"    = "Production Sum (Relative Abundance)",
  "ProdDeg_ratio"     = "log2(Production / Degradation)"
)

###############################################################################
# SECTION 3: SPEARMAN CORRELATION  (verbatim; CSV output intentionally omitted)
###############################################################################
cat("\n========== SECTION 3: SPEARMAN CORRELATION ==========\n\n")

all_corr <- data.frame()

for (ind_name in names(indicators)) {
  ind_vals <- indicators[[ind_name]]
  cat("Computing correlations for:", ind_name, "...\n")

  corr_list <- list()
  for (j in seq_len(ncol(bact_filt))) {
    sp_vals <- bact_filt[, j]
    test_res <- tryCatch(
      cor.test(sp_vals, ind_vals, method = "spearman", exact = FALSE),
      error = function(e) list(estimate = NA, p.value = NA)
    )
    corr_list[[j]] <- data.frame(
      Indicator = ind_name,
      Species_full = colnames(bact_filt)[j],
      Species = get_species_name(colnames(bact_filt)[j]),
      rho = as.numeric(test_res$estimate),
      p_value = test_res$p.value,
      stringsAsFactors = FALSE
    )
  }
  corr_df <- do.call(rbind, corr_list)
  corr_df$p_adj <- p.adjust(corr_df$p_value, method = "BH")
  all_corr <- rbind(all_corr, corr_df)
}

cat("Total correlation tests:", nrow(all_corr), "\n")

###############################################################################
# SECTION 4: DUAL FILTER - SPECIES SELECTION  (verbatim; CSV output omitted)
###############################################################################
cat("\n========== SECTION 4: DUAL FILTER ==========\n\n")

diff_abund <- read.csv("results_integrated/bacteria/diff_abundance_species_pooled.csv",
                        stringsAsFactors = FALSE)
sig_diff_species <- diff_abund$Species_full[diff_abund$p_adj < 0.05]
cat("Species with significant diff abundance:", length(sig_diff_species), "\n")

sig_corr <- all_corr %>%
  filter(p_adj < 0.05, abs(rho) > 0.3, Species_full %in% sig_diff_species)
cat("Correlations passing dual filter (FDR<0.05, |rho|>0.3, diff_abund sig):",
    nrow(sig_corr), "\n")

MAX_SPECIES_PER_INDICATOR <- 8

selected_species <- sig_corr %>%
  group_by(Indicator) %>%
  arrange(desc(abs(rho))) %>%
  slice_head(n = MAX_SPECIES_PER_INDICATOR) %>%
  ungroup()

cat("\nSelected species per indicator:\n")
for (ind_name in names(indicators)) {
  sel <- selected_species %>% filter(Indicator == ind_name)
  cat("  ", ind_name, ":", nrow(sel), "species\n")
}

selected_species <- selected_species %>%
  left_join(diff_abund %>% select(Species_full, log2FC, p_adj_diff = p_adj),
            by = "Species_full") %>%
  mutate(Direction = ifelse(log2FC > 0, "Cancer-enriched", "Healthy-enriched"))

###############################################################################
# SECTION 5: SCATTER PLOTS  (verbatim layout; larger fonts + larger panels)
###############################################################################
cat("\n========== SECTION 5: SCATTER PLOTS ==========\n\n")

make_scatter <- function(ind_name, sel_species_df, bact_data, ind_vals,
                         meta_df, ylabel, title_suffix, filename) {

  if (nrow(sel_species_df) == 0) {
    cat("  No species for", ind_name, "- skipping.\n")
    return(invisible(NULL))
  }

  plot_data <- data.frame()
  annotation_data <- data.frame()

  for (i in seq_len(nrow(sel_species_df))) {
    sp_full  <- sel_species_df$Species_full[i]
    sp_short <- sel_species_df$Species[i]
    sp_rho   <- sel_species_df$rho[i]
    sp_padj  <- sel_species_df$p_adj[i]
    sp_dir   <- sel_species_df$Direction[i]

    sp_abund <- bact_data[, sp_full]

    sp_label <- paste0(gsub("_", " ", sp_short), "\n(",
                       ifelse(is.na(sp_dir), "", sp_dir), ")")

    tmp <- data.frame(
      Species = sp_label,
      Species_short = sp_short,
      Species_abundance = sp_abund,
      Indicator_value = ind_vals,
      Group = meta_df$Group,
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(plot_data, tmp)

    annotation_data <- rbind(annotation_data, data.frame(
      Species = sp_label,
      label = paste0("rho=", round(sp_rho, 3), ", ", format_pval(sp_padj)),
      stringsAsFactors = FALSE
    ))
  }

  sp_order <- unique(plot_data$Species)
  plot_data$Species <- factor(plot_data$Species, levels = sp_order)
  annotation_data$Species <- factor(annotation_data$Species, levels = sp_order)

  ann_pos <- plot_data %>%
    group_by(Species) %>%
    summarise(
      x_pos = min(Species_abundance, na.rm = TRUE),
      y_pos = max(Indicator_value, na.rm = TRUE) * 0.95,
      .groups = "drop"
    )
  annotation_data <- merge(annotation_data, ann_pos, by = "Species")

  n_sp <- nrow(sel_species_df)
  ncol_facet <- min(4, n_sp)
  nrow_facet <- ceiling(n_sp / ncol_facet)

  p <- ggplot(plot_data, aes(x = Species_abundance, y = Indicator_value, color = Group)) +
    geom_point(size = 0.8, alpha = 0.3) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, alpha = 0.15,
                aes(group = 1), color = "grey30") +
    facet_wrap(~Species, scales = "free", ncol = ncol_facet) +
    scale_color_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
    geom_text(data = annotation_data,
              aes(x = x_pos, y = y_pos, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = 1,
              size = FS_ANNOT, color = "grey20") +
    labs(
      title = paste0("Integrated 4 CRC Cohorts: Species vs ", title_suffix,
                     " (n=", length(ind_vals), ")"),
      subtitle = "Dual filter: Spearman FDR<0.05 + |rho|>0.3 + Diff. Abundance FDR<0.05",
      x = "Species Relative Abundance (%)",
      y = ylabel
    ) +
    theme_bw() +
    theme(
      plot.title    = element_text(size = FS_TITLE, face = "bold"),
      plot.subtitle = element_text(size = FS_SUBT, color = "grey40"),
      strip.text    = element_text(size = FS_STRIP, face = "bold"),
      axis.title    = element_text(size = FS_AXIS_T),
      axis.text     = element_text(size = FS_AXIS_X),
      legend.title  = element_text(size = FS_LEG_T),
      legend.text   = element_text(size = FS_LEG_X),
      legend.position = "bottom"
    )

  # Larger panels so the bigger fonts have room (was 3.5 in / facet).
  plot_width  <- max(12, ncol_facet * 4.0)
  plot_height <- max(6,  nrow_facet * 4.0)

  .ggsave(filename, p, width = plot_width, height = plot_height, dpi = 200)
  cat("  Saved:", basename(filename), "\n")
}

# --- 5a. K00303 (soxB) scatter ---
sel_K00303 <- selected_species %>% filter(Indicator == "K00303_soxB")
make_scatter("K00303_soxB", sel_K00303, bact_filt, ko_K00303, meta_common,
             indicator_ylabels["K00303_soxB"],
             "K00303 soxB (Degradation)",
             "results_integrated/sarcosine/scatter/scatter_K00303_soxB_species.png")

# --- 5a. K08688 (Creatinase) scatter ---
sel_K08688 <- selected_species %>% filter(Indicator == "K08688_Creatinase")
make_scatter("K08688_Creatinase", sel_K08688, bact_filt, ko_K08688, meta_common,
             indicator_ylabels["K08688_Creatinase"],
             "K08688 Creatinase (Production)",
             "results_integrated/sarcosine/scatter/scatter_K08688_creatinase_species.png")

# --- 5b. Degradation sum scatter ---
sel_deg <- selected_species %>% filter(Indicator == "Degradation_sum")
make_scatter("Degradation_sum", sel_deg, bact_filt, deg_sum, meta_common,
             indicator_ylabels["Degradation_sum"],
             paste0("Degradation Sum (", paste(deg_kos, collapse=","), ")"),
             "results_integrated/sarcosine/scatter/scatter_degradation_sum_species.png")

# --- 5b. Production sum scatter ---
sel_prod <- selected_species %>% filter(Indicator == "Production_sum")
make_scatter("Production_sum", sel_prod, bact_filt, prod_sum, meta_common,
             indicator_ylabels["Production_sum"],
             paste0("Production Sum (", paste(prod_kos, collapse=","), ")"),
             "results_integrated/sarcosine/scatter/scatter_production_sum_species.png")

# --- 5b. Prod/Deg ratio scatter ---
sel_ratio <- selected_species %>% filter(Indicator == "ProdDeg_ratio")
make_scatter("ProdDeg_ratio", sel_ratio, bact_filt, log2_ratio, meta_common,
             indicator_ylabels["ProdDeg_ratio"],
             "Prod/Deg Ratio",
             "results_integrated/sarcosine/scatter/scatter_ProdDeg_ratio_species.png")

# --- 5c. Summary bubble plot ---
cat("\n--- Summary plot ---\n")
if (nrow(selected_species) > 0) {
  summary_df <- selected_species %>%
    mutate(
      neg_log10_fdr = -log10(p_adj),
      abs_log2FC = abs(log2FC),
      Indicator_label = indicator_labels[Indicator]
    )

  p_summary <- ggplot(summary_df,
                       aes(x = rho, y = neg_log10_fdr,
                           size = abs_log2FC, color = Indicator_label)) +
    geom_point(alpha = 0.7) +
    geom_text(aes(label = gsub("_", " ", Species)),
              size = 2, vjust = -1.2, show.legend = FALSE) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    scale_size_continuous(range = c(2, 8), name = "|log2FC|\n(Cancer/Healthy)") +
    scale_color_brewer(palette = "Set1", name = "Sarcosine Indicator") +
    labs(
      title = paste0("Species-Sarcosine Correlation Summary (n=", length(common_samples), ")"),
      subtitle = "Dual filter: Spearman FDR<0.05, |rho|>0.3, Diff. Abundance FDR<0.05",
      x = "Spearman rho",
      y = "-log10(FDR)"
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 7, color = "grey40"),
      legend.position = "right",
      legend.text = element_text(size = 7)
    )

  .ggsave("results_integrated/sarcosine/scatter/summary_species_sarcosine_correlation.png",
          p_summary, width = 12, height = 8, dpi = 200)
  cat("  Saved: summary_species_sarcosine_correlation.png\n")
}

cat("\n=== Done. Scatter figures re-rendered with larger fonts. ===\n")
gc(verbose = FALSE)
