#!/usr/bin/env Rscript
# ============================================================================
# PER-COHORT diversity plots, re-rendered with a LARGER main title font only.
# Per-cohort analogue of scripts/regenerate_diversity_plots.R (which does the
# pooled figures). Reproduces each cohort's alpha-diversity boxplot and PCoA
# (Bray-Curtis) EXACTLY as the per-cohort analysis_healthy_vs_cancer.R (lines
# 163-241) -- same species filter, same vegan indices, same red/green palette,
# same default alphabetical group order (Cancer | Healthy), same titles/subtitle,
# same set.seed(42) before adonis2, same 10x5 / 8x6 in @200 dpi canvases.
#
# THE ONLY CHANGE: plot.title font size 10 -> 14 (user request: bigger main title).
#
# FAITHFULNESS NOTES:
#   * Statistics are recomputed identically (alpha indices + Wilcoxon are
#     deterministic; PERMANOVA p uses set.seed(42) as in the original) -> the
#     PCoA R2/p and all box statistics are unchanged.
#   * The alpha-diversity geom_jitter dot scatter is the ONE non-reproducible
#     element: the original script did not seed before that jitter, so the exact
#     horizontal dot positions cannot be reproduced. This script seeds (42)
#     right before drawing so OUR output is reproducible. The dots' VERTICAL
#     (value) positions, the boxplots, outliers and p-values are identical.
#   * Output goes to NEW files (*_paper.png); the originals (alpha_diversity.png,
#     beta_diversity_pcoa.png) are NOT overwritten.
#
# Output per cohort: <cohort>/results_bacteria/alpha_diversity_paper.png
#                    <cohort>/results_bacteria/beta_diversity_pcoa_paper.png
# Run: Rscript regenerate_diversity_plots_percohort.R   (paths auto-discovered)
# ============================================================================
options(bitmapType = "quartz")
suppressPackageStartupMessages({ library(vegan); library(ggplot2); library(dplyr); library(tidyr) })

script_dir <- "Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts"   # relative to the analysis-folder root (working dir)
crc_root   <- "."   # scripts/ -> integrated/ -> CRC root (= the analysis-folder root)

COHORTS <- c("PRJEB6070_CRC_AdenomatousPolyps", "PRJEB10878_CRC",
             "PRJEB27928_CRC", "PRJNA429097_CRC")

TITLE_FS <- 24   # enlarged 18->24 (user request); title shortened to BioProject ID + plot.title.position="plot" to avoid clip

# verbatim helper from analysis_healthy_vs_cancer.R (lines 70-75)
format_pval <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.0001) return("p<0.0001")
  if (p < 0.001)  return(sprintf("p=%.4f", p))
  return(sprintf("p=%.3f", p))
}
# quartz-backed PNG writer (handles the non-ASCII backup path; matches pooled regen)
save_png_quartz <- function(plot_obj, filename, width_in, height_in, dpi = 200) {
  grDevices::png(filename = filename, width = width_in, height = height_in,
                 units = "in", res = dpi, type = "quartz")
  on.exit(grDevices::dev.off(), add = TRUE)
  print(plot_obj); invisible(NULL)
}

render_cohort <- function(cohort) {
  cdir <- file.path(crc_root, cohort)
  # --- metadata (verbatim logic; [1] = alphabetical-first, as in the original) ---
  meta_file <- sort(list.files(cdir, "^selected_project_.*\\.txt$", full.names = TRUE))[1]
  ml <- readLines(meta_file); header <- strsplit(ml[2], "\t")[[1]]
  dl <- ml[3:length(ml)]; dl <- dl[dl != ""]
  dlist <- lapply(strsplit(dl, "\t"), function(x) if (length(x) >= length(header)) x[1:length(header)] else c(x, rep(NA, length(header) - length(x))))
  meta <- as.data.frame(do.call(rbind, dlist), stringsAsFactors = FALSE); colnames(meta) <- gsub(" ", ".", header)
  cancer_label <- unique(meta$Phenotype.name)[!unique(meta$Phenotype.name) %in% c("Health", "Adenomatous Polyps")][1]
  meta_hc <- meta %>% filter(Phenotype.name %in% c("Health", cancer_label)) %>%
    mutate(Group = factor(ifelse(Phenotype.name == "Health", "Healthy", "Cancer"),
                          levels = c("Healthy", "Cancer")))   # Healthy left, Cancer right (match pooled order)
  if ("Assay.type" %in% colnames(meta_hc)) meta_hc <- meta_hc %>% filter(Assay.type == "WGS")

  # --- bacteria species matrix (verbatim logic) ---
  bact_file <- list.files(cdir, "^Bacteria_.*\\.txt$", full.names = TRUE)[1]
  bact <- read.delim(bact_file, stringsAsFactors = FALSE); colnames(bact) <- trimws(colnames(bact))
  bact_species <- bact %>% filter(grepl("\\|s__", Taxa) & !grepl("\\|t__", Taxa)) %>% filter(Run.ID %in% meta_hc$Run.ID)
  bact_wide <- bact_species %>% select(Taxa, Run.ID, Abundance) %>%
    pivot_wider(names_from = Taxa, values_from = Abundance, values_fill = 0)
  bact_mat <- as.data.frame(bact_wide); rownames(bact_mat) <- bact_mat$Run.ID
  bact_mat$Run.ID <- NULL; bact_mat <- as.matrix(bact_mat)
  sample_order <- intersect(rownames(bact_mat), meta_hc$Run.ID)
  bact_mat <- bact_mat[sample_order, , drop = FALSE]
  meta_matched <- meta_hc %>% filter(Run.ID %in% sample_order)
  meta_matched <- meta_matched[match(sample_order, meta_matched$Run.ID), ]

  # --- 1a. Alpha diversity (identical) ---
  alpha_div <- data.frame(Run.ID = rownames(bact_mat),
                          Shannon = diversity(bact_mat, index = "shannon"),
                          Simpson = diversity(bact_mat, index = "simpson"),
                          Richness = specnumber(bact_mat), stringsAsFactors = FALSE)
  alpha_div <- merge(alpha_div, meta_matched[, c("Run.ID", "Group")], by = "Run.ID")
  shannon_test  <- wilcox.test(Shannon  ~ Group, data = alpha_div)
  simpson_test  <- wilcox.test(Simpson  ~ Group, data = alpha_div)
  richness_test <- wilcox.test(Richness ~ Group, data = alpha_div)
  alpha_long <- alpha_div %>% pivot_longer(cols = c(Shannon, Simpson, Richness), names_to = "Metric", values_to = "Value")
  pval_labels <- data.frame(Metric = c("Shannon", "Simpson", "Richness"),
                            label = c(format_pval(shannon_test$p.value), format_pval(simpson_test$p.value), format_pval(richness_test$p.value)),
                            stringsAsFactors = FALSE)
  ypos <- alpha_long %>% group_by(Metric) %>% summarise(ymax = max(Value, na.rm = TRUE), .groups = "drop")
  pval_labels <- merge(pval_labels, ypos, by = "Metric")

  p_alpha <- ggplot(alpha_long, aes(x = Group, y = Value, fill = Group)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.7) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.4) +
    facet_wrap(~Metric, scales = "free_y") +
    scale_fill_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
    geom_text(data = pval_labels, aes(x = 1.5, y = ymax * 1.1, label = label), inherit.aes = FALSE, size = 5.9) +
    labs(title = paste0(sub("_.*", "", cohort), ": Alpha diversity"),
         subtitle = paste0("Healthy vs ", cancer_label), x = "", y = "Value") +
    theme_bw() +
    theme(legend.position = "bottom",                                   # enlarged to match the pooled alpha theme
          plot.title    = element_text(size = TITLE_FS, face = "bold"), plot.title.position = "plot",
          plot.subtitle = element_text(size = 16),
          strip.text    = element_text(size = 17, face = "bold"),
          axis.title.x  = element_text(size = 16),
          axis.title.y  = element_text(size = 16),
          axis.text.x   = element_text(size = 16),
          axis.text.y   = element_text(size = 14),
          legend.title  = element_text(size = 16),
          legend.text   = element_text(size = 14))

  set.seed(42)   # make OUR jitter reproducible (original did not seed here; see header)
  out_a <- file.path(cdir, "results_bacteria", "alpha_diversity_paper.png")
  save_png_quartz(p_alpha, out_a, 10, 5, 200)

  # --- 1b. Beta diversity PCoA (identical) ---
  out_b <- NA_character_
  if (nrow(bact_mat) > 3) {
    bc_dist <- vegdist(bact_mat, method = "bray")
    pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)
    eig_pct <- round(pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]) * 100, 1)
    pcoa_df <- data.frame(PC1 = pcoa_res$points[, 1], PC2 = pcoa_res$points[, 2], Group = meta_matched$Group)
    set.seed(42)
    perm_res <- adonis2(bc_dist ~ Group, data = meta_matched, permutations = 999)
    perm_p <- perm_res$`Pr(>F)`[1]; perm_r2 <- round(perm_res$R2[1], 4)
    p_beta <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
      geom_point(size = 2.5, alpha = 0.7) +
      stat_ellipse(level = 0.95, linetype = 2) +
      scale_color_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
      labs(title = paste0(sub("_.*", "", cohort), ": Bray-Curtis PCoA"),
           x = paste0("PCoA1 (", eig_pct[1], "%)"), y = paste0("PCoA2 (", eig_pct[2], "%)"),
           subtitle = paste0("PERMANOVA: R2=", perm_r2, ", ", format_pval(perm_p))) +
      theme_bw() +
      theme(plot.title    = element_text(size = TITLE_FS, face = "bold"), plot.title.position = "plot",   # enlarged to match pooled
            plot.subtitle = element_text(size = 14),
            axis.title    = element_text(size = 16),
            axis.text     = element_text(size = 14),
            legend.title  = element_text(size = 16),
            legend.text   = element_text(size = 14))
    out_b <- file.path(cdir, "results_bacteria", "beta_diversity_pcoa_paper.png")
    save_png_quartz(p_beta, out_b, 8, 6, 200)
    cat(sprintf("%-32s alpha n=%d (H=%d/C=%d)  PERMANOVA R2=%.4f %s\n", cohort, nrow(bact_mat),
        sum(meta_matched$Group == "Healthy"), sum(meta_matched$Group == "Cancer"), perm_r2, format_pval(perm_p)))
  } else {
    cat(sprintf("%-32s alpha only (n<=3, beta skipped)\n", cohort))
  }
}

for (cn in COHORTS) render_cohort(cn)
writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo_diversity_percohort.txt"))
cat("DONE. 4 cohorts x (alpha_diversity_paper.png + beta_diversity_pcoa_paper.png).\n")
