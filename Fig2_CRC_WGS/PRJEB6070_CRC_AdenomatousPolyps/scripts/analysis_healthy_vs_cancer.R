###############################################################################
# PRJEB6070_CRC_AdenomatousPolyps: Healthy vs Cancer Comprehensive Analysis
# Cancer type: CRC / AdenomatousPolyps
#
# Analysis pipeline:
#   Part 1: Bacteria differential abundance (alpha/beta diversity, volcano, heatmap)
#   Part 2: KEGG KO differential abundance (volcano)
#   Part 3: Sarcosine metabolism enzyme analysis (barplot, boxplot)
#   Part 3b: Sarcosine KO analysis with prevalence filter (>=10%)
#   Part 4: Sarcosine-associated species (correlation heatmap with group annotation)
#   Part 4b: Sarcosine-associated species with prevalence filter (>=10%)
#   Part 5: Sarcosine Production vs Degradation pathway balance
#   Part 6: Sarcosine-associated bacteria differential abundance
###############################################################################

library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(jsonlite)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(ggrepel)

# Fix: ragg backend fails with Korean paths — force base png device for ggsave
options(bitmapType = "cairo")
.ggsave <- function(...) ggplot2::ggsave(..., device = grDevices::png)

# Set working directory
project_dir <- "PRJEB6070_CRC_AdenomatousPolyps"   # relative to the analysis-folder root (working dir)
setwd(project_dir)

# Create output directories
dir.create("results_bacteria", showWarnings = FALSE)
dir.create("results_kegg", showWarnings = FALSE)
dir.create("results_sarcosine", showWarnings = FALSE)

cat("=== PRJEB6070_CRC_AdenomatousPolyps: Healthy vs Cancer Analysis ===\n\n")

###############################################################################
# Sarcosine-related KOs (verified from KEGG)
###############################################################################
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
  EC = c("EC:1.5.3.1", "EC:1.5.3.1/1.5.3.24", "EC:1.5.3.1/1.5.3.24", "EC:1.5.3.1/1.5.3.24",
         "EC:1.5.3.1/1.5.3.24", "EC:1.5.3.1/1.5.3.7",
         "EC:1.5.8.4", "EC:2.1.1.20", "EC:3.5.3.3"),
  Role = c("Degradation", "Degradation", "Degradation", "Degradation",
           "Degradation", "Degradation",
           "Production", "Production", "Production"),
  Short = c("SOX (mono)", "soxA", "soxB", "soxD", "soxG", "PIPOX",
            "DMGDH", "GNMT", "Creatinase"),
  stringsAsFactors = FALSE
)

###############################################################################
# Helper functions
###############################################################################
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

###############################################################################
# 1. Load metadata and identify groups
###############################################################################
cat("--- Loading metadata ---\n")
meta_file <- list.files(pattern = "^selected_project_.*\\.txt$")[1]
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

cat("Phenotypes found:\n")
print(table(meta$Phenotype.name))

# Auto-detect cancer group (non-Health, non-Adenomatous Polyps)
all_pheno <- unique(meta$Phenotype.name)
cancer_pheno <- all_pheno[!all_pheno %in% c("Health", "Adenomatous Polyps")]

if (length(cancer_pheno) == 0) {
  stop("No cancer phenotype found in metadata.")
}

cancer_label <- cancer_pheno[1]
cat("Cancer phenotype:", cancer_label, "\n")

meta_hc <- meta %>%
  filter(Phenotype.name %in% c("Health", cancer_label)) %>%
  mutate(Group = ifelse(Phenotype.name == "Health", "Healthy", "Cancer"))

# Filter to WGS samples only (exclude 16S)
if ("Assay.type" %in% colnames(meta_hc)) {
  n_before <- nrow(meta_hc)
  meta_hc <- meta_hc %>% filter(Assay.type == "WGS")
  cat("Assay type filter: kept", nrow(meta_hc), "WGS samples out of", n_before, "total\n")
}

cat("Sample counts: Healthy=", sum(meta_hc$Group == "Healthy"),
    ", Cancer=", sum(meta_hc$Group == "Cancer"), "\n")

if (sum(meta_hc$Group == "Healthy") < 5 | sum(meta_hc$Group == "Cancer") < 5) {
  stop("Too few samples in one group (< 5). Cannot proceed.")
}

###############################################################################
# PART 1: BACTERIA ANALYSIS
###############################################################################
cat("\n\n========== PART 1: BACTERIA ANALYSIS ==========\n")

bact_file <- list.files(pattern = "^Bacteria_.*\\.txt$")[1]
bact <- read.delim(bact_file, stringsAsFactors = FALSE)
colnames(bact) <- trimws(colnames(bact))

# Species level only
bact_species <- bact %>%
  filter(grepl("\\|s__", Taxa) & !grepl("\\|t__", Taxa)) %>%
  filter(Run.ID %in% meta_hc$Run.ID)

bact_wide <- bact_species %>%
  select(Taxa, Run.ID, Abundance) %>%
  pivot_wider(names_from = Taxa, values_from = Abundance, values_fill = 0)

bact_mat <- as.data.frame(bact_wide)
rownames(bact_mat) <- bact_mat$Run.ID
bact_mat$Run.ID <- NULL
bact_mat <- as.matrix(bact_mat)

sample_order <- intersect(rownames(bact_mat), meta_hc$Run.ID)
bact_mat <- bact_mat[sample_order, , drop = FALSE]
meta_matched <- meta_hc %>% filter(Run.ID %in% sample_order)
meta_matched <- meta_matched[match(sample_order, meta_matched$Run.ID), ]

cat("Bacteria matrix:", nrow(bact_mat), "samples x", ncol(bact_mat), "species\n")

# 1a. Alpha Diversity
cat("\n--- Alpha Diversity ---\n")
alpha_div <- data.frame(
  Run.ID = rownames(bact_mat),
  Shannon = diversity(bact_mat, index = "shannon"),
  Simpson = diversity(bact_mat, index = "simpson"),
  Richness = specnumber(bact_mat),
  stringsAsFactors = FALSE
)
alpha_div <- merge(alpha_div, meta_matched[, c("Run.ID", "Group")], by = "Run.ID")

shannon_test <- wilcox.test(Shannon ~ Group, data = alpha_div)
simpson_test <- wilcox.test(Simpson ~ Group, data = alpha_div)
richness_test <- wilcox.test(Richness ~ Group, data = alpha_div)

alpha_long <- alpha_div %>%
  pivot_longer(cols = c(Shannon, Simpson, Richness), names_to = "Metric", values_to = "Value")

pval_labels <- data.frame(
  Metric = c("Shannon", "Simpson", "Richness"),
  label = c(format_pval(shannon_test$p.value),
            format_pval(simpson_test$p.value),
            format_pval(richness_test$p.value)),
  stringsAsFactors = FALSE
)
ypos <- alpha_long %>%
  group_by(Metric) %>%
  summarise(ymax = max(Value, na.rm = TRUE), .groups = "drop")
pval_labels <- merge(pval_labels, ypos, by = "Metric")

p_alpha <- ggplot(alpha_long, aes(x = Group, y = Value, fill = Group)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.4) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
  geom_text(data = pval_labels, aes(x = 1.5, y = ymax * 1.1, label = label),
            inherit.aes = FALSE, size = 3.5) +
  labs(title = paste0("PRJEB6070_CRC_AdenomatousPolyps: Alpha Diversity (Healthy vs ", cancer_label, ")"),
       x = "", y = "Value") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 10))

.ggsave("results_bacteria/alpha_diversity.png", p_alpha, width = 10, height = 5, dpi = 200)
cat("Alpha diversity plot saved.\n")

# 1b. Beta Diversity
cat("\n--- Beta Diversity ---\n")
if (nrow(bact_mat) > 3) {
  bc_dist <- vegdist(bact_mat, method = "bray")
  pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)
  eig_pct <- round(pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]) * 100, 1)

  pcoa_df <- data.frame(
    PC1 = pcoa_res$points[, 1],
    PC2 = pcoa_res$points[, 2],
    Group = meta_matched$Group
  )

  set.seed(42)
  perm_res <- adonis2(bc_dist ~ Group, data = meta_matched, permutations = 999)
  perm_p <- perm_res$`Pr(>F)`[1]
  perm_r2 <- round(perm_res$R2[1], 4)

  p_beta <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group)) +
    geom_point(size = 2.5, alpha = 0.7) +
    stat_ellipse(level = 0.95, linetype = 2) +
    scale_color_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
    labs(title = paste0("PRJEB6070_CRC_AdenomatousPolyps: PCoA (Bray-Curtis)"),
         x = paste0("PCoA1 (", eig_pct[1], "%)"),
         y = paste0("PCoA2 (", eig_pct[2], "%)"),
         subtitle = paste0("PERMANOVA: R2=", perm_r2, ", ", format_pval(perm_p))) +
    theme_bw() +
    theme(plot.title = element_text(size = 10, face = "bold"))

  .ggsave("results_bacteria/beta_diversity_pcoa.png", p_beta, width = 8, height = 6, dpi = 200)
  cat("Beta diversity plot saved.\n")
}

# 1c. Differential Abundance
cat("\n--- Differential Abundance (Species) ---\n")
prevalence <- colSums(bact_mat > 0) / nrow(bact_mat)
bact_filt <- bact_mat[, prevalence >= 0.1, drop = FALSE]
cat("Species after prevalence filter:", ncol(bact_filt), "\n")

healthy_idx <- which(meta_matched$Group == "Healthy")
cancer_idx <- which(meta_matched$Group == "Cancer")

diff_results <- data.frame(
  Species_full = colnames(bact_filt),
  Species = sapply(colnames(bact_filt), get_species_name),
  Mean_Healthy = colMeans(bact_filt[healthy_idx, , drop = FALSE]),
  Mean_Cancer = colMeans(bact_filt[cancer_idx, , drop = FALSE]),
  stringsAsFactors = FALSE
)
diff_results$log2FC <- log2((diff_results$Mean_Cancer + 1e-6) / (diff_results$Mean_Healthy + 1e-6))

pvals <- sapply(1:ncol(bact_filt), function(i) {
  tryCatch(wilcox.test(bact_filt[healthy_idx, i], bact_filt[cancer_idx, i])$p.value,
           error = function(e) NA)
})
diff_results$p_value <- pvals
diff_results$p_adj <- p.adjust(pvals, method = "BH")
diff_results <- diff_results %>% arrange(p_adj)

write.csv(diff_results, "results_bacteria/diff_abundance_species.csv", row.names = FALSE)

sig_species <- diff_results %>% filter(p_adj < 0.05)
cat("Significant species (FDR < 0.05):", nrow(sig_species), "\n")

# Volcano plot
diff_results$Significance <- "Not Significant"
diff_results$Significance[diff_results$p_adj < 0.05 & diff_results$log2FC > 1] <- "Enriched in Cancer"
diff_results$Significance[diff_results$p_adj < 0.05 & diff_results$log2FC < -1] <- "Enriched in Healthy"

top_to_label <- diff_results %>% filter(p_adj < 0.05) %>% arrange(p_adj) %>% head(15)

p_volcano <- ggplot(diff_results, aes(x = log2FC, y = -log10(p_adj), color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Enriched in Cancer" = "#B2182B",
                                 "Enriched in Healthy" = "#1B7837",
                                 "Not Significant" = "grey60")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = top_to_label, aes(label = Species),
                  size = 2.5, max.overlaps = 20) +
  labs(title = paste0("PRJEB6070_CRC_AdenomatousPolyps: Differential Species"),
       subtitle = paste0("Healthy vs ", cancer_label),
       x = "log2FC (Cancer/Healthy)", y = "-log10(adj. p-value)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = "bottom")

.ggsave("results_bacteria/volcano_species.png", p_volcano, width = 10, height = 7, dpi = 200)
cat("Volcano plot saved.\n")

# Heatmap of top significant species
if (nrow(sig_species) > 1) {
  top_n <- min(30, nrow(sig_species))
  top_sp <- sig_species$Species_full[1:top_n]
  hm_mat <- bact_filt[, top_sp, drop = FALSE]
  hm_scaled <- scale(hm_mat)
  hm_scaled[hm_scaled > 3] <- 3
  hm_scaled[hm_scaled < -3] <- -3
  colnames(hm_scaled) <- sapply(colnames(hm_scaled), get_species_name)

  annotation_row <- data.frame(Group = meta_matched$Group)
  rownames(annotation_row) <- rownames(hm_scaled)
  ann_colors <- list(Group = c("Healthy" = "#1B7837", "Cancer" = "#B2182B"))

  png("results_bacteria/heatmap_top_species.png",
      width = 12, height = max(8, top_n * 0.25 + 2), units = "in", res = 200)
  pheatmap(t(hm_scaled),
           annotation_col = annotation_row,
           annotation_colors = ann_colors,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           cluster_rows = TRUE, cluster_cols = TRUE,
           show_colnames = FALSE, fontsize_row = 7,
           main = paste0("PRJEB6070_CRC_AdenomatousPolyps: Top Differential Species (Z-score)"))
  dev.off()
  cat("Heatmap saved.\n")
}

###############################################################################
# PART 2: KEGG KO ANALYSIS
###############################################################################
cat("\n\n========== PART 2: KEGG KO ANALYSIS ==========\n")
cat("Loading KO data...\n")
ko_raw <- fromJSON("KO_relative_abundance.tsv")
ko_raw <- ko_raw %>% filter(run_id %in% meta_hc$Run.ID)

ko_wide <- ko_raw %>%
  select(ko, run_id, abundance) %>%
  pivot_wider(names_from = ko, values_from = abundance, values_fill = 0)

ko_mat <- as.data.frame(ko_wide)
rownames(ko_mat) <- ko_mat$run_id
ko_mat$run_id <- NULL
ko_mat <- as.matrix(ko_mat)

ko_samples <- intersect(rownames(ko_mat), meta_hc$Run.ID)
ko_mat <- ko_mat[ko_samples, , drop = FALSE]
meta_ko <- meta_hc %>% filter(Run.ID %in% ko_samples)
meta_ko <- meta_ko[match(ko_samples, meta_ko$Run.ID), ]

cat("KO matrix:", nrow(ko_mat), "samples x", ncol(ko_mat), "KOs\n")
rm(ko_raw); gc(verbose = FALSE)

# Differential KO
ko_prev <- colSums(ko_mat > 0) / nrow(ko_mat)
ko_filt <- ko_mat[, ko_prev >= 0.1, drop = FALSE]

h_idx <- which(meta_ko$Group == "Healthy")
c_idx <- which(meta_ko$Group == "Cancer")

ko_diff <- data.frame(
  KO = colnames(ko_filt),
  Mean_Healthy = colMeans(ko_filt[h_idx, , drop = FALSE]),
  Mean_Cancer = colMeans(ko_filt[c_idx, , drop = FALSE]),
  stringsAsFactors = FALSE
)
ko_diff$log2FC <- log2((ko_diff$Mean_Cancer + 1e-8) / (ko_diff$Mean_Healthy + 1e-8))

ko_pvals <- sapply(1:ncol(ko_filt), function(i) {
  tryCatch(wilcox.test(ko_filt[h_idx, i], ko_filt[c_idx, i])$p.value,
           error = function(e) NA)
})
ko_diff$p_value <- ko_pvals
ko_diff$p_adj <- p.adjust(ko_pvals, method = "BH")
ko_diff <- ko_diff %>% arrange(p_adj)

write.csv(ko_diff, "results_kegg/diff_KO_abundance.csv", row.names = FALSE)

sig_kos <- ko_diff %>% filter(p_adj < 0.05)
cat("Significant KOs (FDR < 0.05):", nrow(sig_kos), "\n")

# KO Volcano
ko_diff$Significance <- "Not Significant"
ko_diff$Significance[ko_diff$p_adj < 0.05 & ko_diff$log2FC > 0.5] <- "Enriched in Cancer"
ko_diff$Significance[ko_diff$p_adj < 0.05 & ko_diff$log2FC < -0.5] <- "Enriched in Healthy"

top_kos <- ko_diff %>% filter(p_adj < 0.05) %>% head(20)

p_ko_volcano <- ggplot(ko_diff, aes(x = log2FC, y = -log10(p_adj), color = Significance)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Enriched in Cancer" = "#B2182B",
                                 "Enriched in Healthy" = "#1B7837",
                                 "Not Significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(data = top_kos, aes(label = KO), size = 2, max.overlaps = 20) +
  labs(title = paste0("PRJEB6070_CRC_AdenomatousPolyps: Differential KO Abundance"),
       subtitle = paste0("Healthy vs ", cancer_label),
       x = "log2FC (Cancer/Healthy)", y = "-log10(adj. p-value)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"), legend.position = "bottom")

.ggsave("results_kegg/volcano_KO.png", p_ko_volcano, width = 10, height = 7, dpi = 200)
cat("KO volcano plot saved.\n")

###############################################################################
# PART 3: SARCOSINE METABOLISM ANALYSIS
###############################################################################
cat("\n\n========== PART 3: SARCOSINE METABOLISM ==========\n")

available_kos <- intersect(sarcosine_kos$KO, colnames(ko_mat))
cat("Sarcosine KOs found:", length(available_kos), "-", paste(available_kos, collapse = ", "), "\n")

if (length(available_kos) > 0) {
  sarc_info <- sarcosine_kos %>% filter(KO %in% available_kos)
  sarc_mat <- ko_mat[, available_kos, drop = FALSE]

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
  write.csv(sarc_results, "results_sarcosine/sarcosine_KO_comparison.csv", row.names = FALSE)

  cat("Sarcosine KO results:\n")
  print(sarc_results[, c("KO", "Enzyme", "Role", "Mean_Healthy", "Mean_Cancer", "p_value")])

  # Bar plot and Boxplot for significant sarcosine KOs
  sarc_sig <- sarc_results %>% filter(p_value < 0.05)
  cat("Significant sarcosine KOs (p<0.05):", nrow(sarc_sig), "\n")

  if (nrow(sarc_sig) > 0) {
    sarc_plot_data <- sarc_sig %>%
      select(KO, Enzyme, EC, Role, Mean_Healthy, Mean_Cancer, p_value) %>%
      pivot_longer(cols = c(Mean_Healthy, Mean_Cancer),
                   names_to = "Group", values_to = "Abundance") %>%
      mutate(Group = gsub("Mean_", "", Group),
             Label = paste0(Enzyme, "\n(", Role, ") [", EC, "]"))

    sarc_pval_labels <- sarc_sig %>%
      mutate(Label = paste0(Enzyme, "\n(", Role, ") [", EC, "]"),
             pval_text = sapply(p_value, format_pval))
    ypos_data <- sarc_plot_data %>%
      group_by(Label) %>%
      summarise(ymax = max(Abundance), .groups = "drop")
    sarc_pval_labels <- merge(sarc_pval_labels, ypos_data, by = "Label")

    p_sarc <- ggplot(sarc_plot_data,
                     aes(x = reorder(Label, -Abundance), y = Abundance, fill = Group)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(data = sarc_pval_labels,
                aes(x = Label, y = ymax * 1.15, label = pval_text),
                inherit.aes = FALSE, size = 3, fontface = "italic") +
      scale_fill_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
      labs(title = paste0("PRJEB6070_CRC_AdenomatousPolyps: Sarcosine Enzyme Gene Abundance"),
           subtitle = paste0("Healthy vs ", cancer_label, " (significant only, p<0.05)"),
           x = "", y = "Mean Relative Abundance") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            plot.title = element_text(size = 10, face = "bold"),
            legend.position = "bottom",
            plot.margin = margin(t = 10, r = 10, b = 20, l = 20))

    fig_w <- max(8, nrow(sarc_sig) * 1.5 + 2)
    .ggsave("results_sarcosine/sarcosine_enzyme_barplot.png", p_sarc,
           width = min(fig_w, 16), height = 7, dpi = 200)

    # Boxplot per significant KO
    sarc_box_data <- data.frame()
    for (i in 1:nrow(sarc_sig)) {
      ko_id <- sarc_sig$KO[i]
      tmp <- data.frame(
        Abundance = sarc_mat[, ko_id], Group = meta_ko$Group,
        Label = paste0(sarc_sig$KO[i], ": ", sarc_sig$Enzyme[i],
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

    p_sarc_box <- ggplot(sarc_box_data, aes(x = Group, y = Abundance, fill = Group)) +
      geom_boxplot(outlier.shape = 21, alpha = 0.7) +
      geom_jitter(width = 0.15, size = 0.8, alpha = 0.3) +
      facet_wrap(~Label, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
      geom_text(data = box_pval, aes(x = 1.5, y = ymax * 1.15, label = pval_text),
                inherit.aes = FALSE, size = 3) +
      labs(title = paste0("PRJEB6070_CRC_AdenomatousPolyps: Sarcosine Enzyme KOs"),
           subtitle = paste0("Healthy vs ", cancer_label, " (p<0.05)"),
           x = "", y = "Relative Abundance") +
      theme_bw() +
      theme(plot.title = element_text(size = 10, face = "bold"),
            legend.position = "bottom", strip.text = element_text(size = 7))

    box_nrow <- ceiling(nrow(sarc_sig) / 3)
    .ggsave("results_sarcosine/sarcosine_enzyme_boxplot.png", p_sarc_box,
           width = 12, height = max(4, box_nrow * 3.5), dpi = 200)
    cat("Sarcosine enzyme plots saved.\n")
  }

  #############################################################################
  # PART 3b: SARCOSINE KO ANALYSIS WITH PREVALENCE FILTER (>=10%)
  #############################################################################
  cat("\n\n========== PART 3b: SARCOSINE KOs (PREVALENCE FILTERED) ==========\n")

  PREV_THRESHOLD <- 10  # percent

  # Calculate prevalence for each KO
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

  write.csv(sarc_results,
            "results_sarcosine/sarcosine_KO_comparison_filtered.csv", row.names = FALSE)

  sarc_filtered <- sarc_results[sarc_results$Pass_Prevalence, ]
  sarc_filt_sig <- sarc_filtered[!is.na(sarc_filtered$p_value) & sarc_filtered$p_value < 0.05, ]

  cat("KOs passing prevalence (>=", PREV_THRESHOLD, "%):", nrow(sarc_filtered), "/", nrow(sarc_results), "\n")
  cat("Significant after filter:", nrow(sarc_filt_sig), "\n")

  if (nrow(sarc_filt_sig) > 0) {
    # Filtered Barplot
    filt_plot_data <- sarc_filt_sig %>%
      select(KO, Enzyme, EC, Role, Mean_Healthy, Mean_Cancer, p_value) %>%
      pivot_longer(cols = c(Mean_Healthy, Mean_Cancer),
                   names_to = "Group", values_to = "Abundance") %>%
      mutate(Group = gsub("Mean_", "", Group),
             Label = paste0(Enzyme, "\n(", Role, ") [", EC, "]"))

    filt_pval_labels <- sarc_filt_sig %>%
      mutate(Label = paste0(Enzyme, "\n(", Role, ") [", EC, "]"),
             pval_text = sapply(p_value, format_pval))
    filt_ypos <- filt_plot_data %>%
      group_by(Label) %>%
      summarise(ymax = max(Abundance), .groups = "drop")
    filt_pval_labels <- merge(filt_pval_labels, filt_ypos, by = "Label")

    p_filt_bar <- ggplot(filt_plot_data,
                    aes(x = reorder(Label, -Abundance), y = Abundance, fill = Group)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      geom_text(data = filt_pval_labels,
                aes(x = Label, y = ymax * 1.15, label = pval_text),
                inherit.aes = FALSE, size = 3, fontface = "italic") +
      scale_fill_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
      labs(title = paste0("PRJEB6070_CRC_AdenomatousPolyps: Sarcosine Enzyme Gene Abundance (Prevalence Filtered)"),
           subtitle = paste0("Healthy vs ", cancer_label,
                            " | Prevalence >=", PREV_THRESHOLD, "% | Significant only (p<0.05)"),
           x = "", y = "Mean Relative Abundance") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            plot.title = element_text(size = 10, face = "bold"),
            legend.position = "bottom",
            plot.margin = margin(t = 10, r = 10, b = 20, l = 20))

    filt_fig_w <- max(8, nrow(sarc_filt_sig) * 1.5 + 2)
    .ggsave("results_sarcosine/sarcosine_enzyme_barplot_filtered.png", p_filt_bar,
           width = min(filt_fig_w, 16), height = 7, dpi = 200)

    # Filtered Boxplot
    filt_box_data <- data.frame()
    for (i in 1:nrow(sarc_filt_sig)) {
      ko_id <- sarc_filt_sig$KO[i]
      tmp <- data.frame(
        Abundance = sarc_mat[, ko_id], Group = meta_ko$Group,
        Label = paste0(sarc_filt_sig$KO[i], ": ", sarc_filt_sig$Enzyme[i],
                       "\n(", sarc_filt_sig$Role[i], ") [", sarc_filt_sig$EC[i], "]",
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

    p_filt_box <- ggplot(filt_box_data, aes(x = Group, y = Abundance, fill = Group)) +
      geom_boxplot(outlier.shape = 21, alpha = 0.7) +
      geom_jitter(width = 0.15, size = 0.8, alpha = 0.3) +
      facet_wrap(~Label, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
      geom_text(data = filt_box_pval, aes(x = 1.5, y = ymax * 1.15, label = pval_text),
                inherit.aes = FALSE, size = 3) +
      labs(title = paste0("PRJEB6070_CRC_AdenomatousPolyps: Sarcosine Enzyme KOs (Prevalence Filtered)"),
           subtitle = paste0("Healthy vs ", cancer_label,
                            " | Prevalence >=", PREV_THRESHOLD, "% | p<0.05"),
           x = "", y = "Relative Abundance") +
      theme_bw() +
      theme(plot.title = element_text(size = 10, face = "bold"),
            legend.position = "bottom", strip.text = element_text(size = 7))

    filt_box_nrow <- ceiling(nrow(sarc_filt_sig) / 3)
    .ggsave("results_sarcosine/sarcosine_enzyme_boxplot_filtered.png", p_filt_box,
           width = 12, height = max(4, filt_box_nrow * 3.5), dpi = 200)
    cat("Filtered sarcosine enzyme plots saved.\n")
  } else {
    cat("No significant KOs after prevalence filter.\n")
  }

  #############################################################################
  # PART 4: SARCOSINE-ASSOCIATED SPECIES (with heatmap group annotation)
  #############################################################################
  cat("\n\n========== PART 4: SARCOSINE-ASSOCIATED SPECIES ==========\n")

  common_samples <- intersect(rownames(bact_filt), rownames(sarc_mat))
  if (length(common_samples) > 10) {
    bact_common <- bact_filt[common_samples, , drop = FALSE]
    sarc_common <- sarc_mat[common_samples, , drop = FALSE]

    cor_results <- data.frame()
    for (ko_id in available_kos) {
      ko_vals <- sarc_common[, ko_id]
      enzyme_info <- sarc_info %>% filter(KO == ko_id)
      if (sd(ko_vals) == 0) next
      for (sp in colnames(bact_common)) {
        sp_vals <- bact_common[, sp]
        if (sd(sp_vals) == 0) next
        cor_test <- suppressWarnings(cor.test(ko_vals, sp_vals, method = "spearman"))
        cor_results <- rbind(cor_results, data.frame(
          KO = ko_id, Enzyme = enzyme_info$Enzyme, Role = enzyme_info$Role,
          Species_full = sp, Species = get_species_name(sp),
          rho = cor_test$estimate, p_value = cor_test$p.value,
          stringsAsFactors = FALSE))
      }
    }

    if (nrow(cor_results) > 0) {
      cor_results$p_adj <- p.adjust(cor_results$p_value, method = "BH")
      cor_results <- cor_results %>% arrange(p_adj)
      write.csv(cor_results, "results_sarcosine/sarcosine_species_correlation.csv",
                row.names = FALSE)

      sig_cors <- cor_results %>% filter(p_adj < 0.05 & abs(rho) > 0.3)
      cat("Significant species-KO correlations:", nrow(sig_cors), "\n")

      if (nrow(sig_cors) > 0) {
        sig_species_list <- unique(sig_cors$Species_full)
        sig_ko_list <- unique(sig_cors$KO)

        if (length(sig_species_list) > 1 & length(sig_ko_list) >= 1) {
          # --- Heatmap with group annotation ---

          # Compute enrichment (species-level aggregation + zero-fill, consistent with Part 6)
          bact_sp_agg <- bact_species %>%
            filter(Run.ID %in% meta_hc$Run.ID) %>%
            mutate(Species_taxa = sub("\\|t__.*$", "", Taxa)) %>%
            group_by(Species_taxa, Run.ID) %>%
            summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
          bact_sp_wide <- bact_sp_agg %>%
            pivot_wider(names_from = Run.ID, values_from = Abundance, values_fill = 0)
          healthy_runs_enr <- meta_hc$Run.ID[meta_hc$Group == "Healthy"]
          cancer_runs_enr  <- meta_hc$Run.ID[meta_hc$Group == "Cancer"]
          h_cols <- intersect(healthy_runs_enr, colnames(bact_sp_wide))
          c_cols <- intersect(cancer_runs_enr, colnames(bact_sp_wide))
          sp_enrichment_df <- data.frame(
            Species_taxa = bact_sp_wide$Species_taxa,
            Healthy = if (length(h_cols) > 0) rowMeans(as.matrix(bact_sp_wide[, h_cols, drop = FALSE])) else 0,
            Cancer  = if (length(c_cols) > 0) rowMeans(as.matrix(bact_sp_wide[, c_cols, drop = FALSE])) else 0,
            stringsAsFactors = FALSE
          ) %>% mutate(Enriched = ifelse(Cancer >= Healthy, "Cancer", "Healthy"))
          sp_enrichment_map <- setNames(sp_enrichment_df$Enriched, sp_enrichment_df$Species_taxa)
          # Map back to full Taxa strings
          taxa_to_sp <- bact_species %>%
            mutate(Species_taxa = sub("\\|t__.*$", "", Taxa)) %>%
            select(Taxa, Species_taxa) %>% distinct()
          enrichment_map <- setNames(
            sp_enrichment_map[taxa_to_sp$Species_taxa],
            taxa_to_sp$Taxa
          )

          # Top species by max |rho|
          sp_max_rho <- sig_cors %>%
            group_by(Species_full) %>%
            summarise(max_rho = max(abs(rho)), .groups = "drop") %>%
            arrange(desc(max_rho)) %>% head(30)

          top_sp_full <- sp_max_rho$Species_full
          top_sp_names <- sapply(top_sp_full, get_species_name)

          # Build correlation matrix
          cor_matrix <- matrix(0, nrow = length(top_sp_full), ncol = length(sig_ko_list))
          rownames(cor_matrix) <- top_sp_names
          colnames(cor_matrix) <- sig_ko_list

          pval_matrix <- matrix(NA, nrow = length(top_sp_full), ncol = length(sig_ko_list))
          rownames(pval_matrix) <- top_sp_names
          colnames(pval_matrix) <- sig_ko_list

          for (r in 1:nrow(sig_cors)) {
            sp_name <- get_species_name(sig_cors$Species_full[r])
            ko_name <- sig_cors$KO[r]
            if (sp_name %in% rownames(cor_matrix) & ko_name %in% colnames(cor_matrix)) {
              cor_matrix[sp_name, ko_name] <- sig_cors$rho[r]
              pval_matrix[sp_name, ko_name] <- sig_cors$p_adj[r]
            }
          }

          # Significance labels
          sig_labels <- matrix("", nrow = nrow(cor_matrix), ncol = ncol(cor_matrix))
          rownames(sig_labels) <- rownames(cor_matrix)
          colnames(sig_labels) <- colnames(cor_matrix)
          for (i in 1:nrow(pval_matrix)) {
            for (j in 1:ncol(pval_matrix)) {
              p <- pval_matrix[i, j]
              if (!is.na(p)) {
                if (p < 0.001) sig_labels[i, j] <- "***"
                else if (p < 0.01) sig_labels[i, j] <- "**"
                else if (p < 0.05) sig_labels[i, j] <- "*"
              }
            }
          }

          # Map enrichment
          sp_enrichment <- sapply(top_sp_full, function(sp) {
            if (sp %in% names(enrichment_map)) return(enrichment_map[sp])
            return("Unknown")
          })
          names(sp_enrichment) <- top_sp_names

          row_annotation <- data.frame(
            `Enriched in` = sp_enrichment[rownames(cor_matrix)],
            check.names = FALSE, stringsAsFactors = FALSE
          )
          rownames(row_annotation) <- rownames(cor_matrix)

          # Order species by group
          healthy_species <- rownames(row_annotation)[row_annotation$`Enriched in` == "Healthy"]
          cancer_species <- rownames(row_annotation)[row_annotation$`Enriched in` == "Cancer"]

          order_within_group <- function(species_names, mat) {
            if (length(species_names) <= 1) return(species_names)
            sub_mat <- mat[species_names, , drop = FALSE]
            d <- dist(sub_mat, method = "euclidean")
            hc <- hclust(d, method = "complete")
            return(species_names[hc$order])
          }

          healthy_ordered <- order_within_group(healthy_species, cor_matrix)
          cancer_ordered <- order_within_group(cancer_species, cor_matrix)
          final_order <- c(healthy_ordered, cancer_ordered)

          cor_matrix_ordered <- cor_matrix[final_order, , drop = FALSE]
          sig_labels_ordered <- sig_labels[final_order, , drop = FALSE]
          row_annotation_ordered <- row_annotation[final_order, , drop = FALSE]

          # KO column labels: KO / Enzyme (EC:x.x.x.x) / Role
          ko_labels <- sapply(colnames(cor_matrix_ordered), function(k) {
            info <- sarc_info %>% filter(KO == k)
            if (nrow(info) > 0) {
              paste0(k, "\n", info$Short[1], " (", info$EC[1], ")\n", info$Role[1])
            } else k
          })
          colnames(cor_matrix_ordered) <- ko_labels
          colnames(sig_labels_ordered) <- ko_labels

          ann_colors <- list(`Enriched in` = c("Healthy" = "#1B7837", "Cancer" = "#B2182B"))

          n_healthy_sp <- length(healthy_ordered)
          gaps_row <- if (n_healthy_sp > 0 & n_healthy_sp < nrow(cor_matrix_ordered)) n_healthy_sp else NULL

          n_sp <- nrow(cor_matrix_ordered)
          n_ko <- ncol(cor_matrix_ordered)
          cell_w <- max(90, 130 - n_ko * 5)
          cell_h <- max(14, 20 - n_sp * 0.15)

          pheatmap(cor_matrix_ordered,
                   filename = "results_sarcosine/sarcosine_species_correlation_heatmap.png",
                   color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   gaps_row = gaps_row,
                   annotation_row = row_annotation_ordered,
                   annotation_colors = ann_colors,
                   display_numbers = sig_labels_ordered,
                   number_color = "black", fontsize_number = 10,
                   fontsize_row = 7, fontsize_col = 8,
                   cellwidth = cell_w, cellheight = cell_h,
                   width = max(10, n_ko * 1.5 + 5),
                   height = max(7, n_sp * 0.35 + 4),
                   main = paste0("PRJEB6070_CRC_AdenomatousPolyps: Species-Sarcosine KO Correlation\n",
                                 "(Healthy vs ", cancer_label,
                                 ")  [* p<0.05, ** p<0.01, *** p<0.001]"),
                   angle_col = 0)
          cat("Sarcosine species correlation heatmap saved (with group annotation).\n")
        }

        # Species barplot
        species_summary <- sig_cors %>%
          group_by(Species, Role) %>%
          summarise(mean_rho = mean(rho), max_rho = max(abs(rho)),
                    n_kos = n(), .groups = "drop") %>%
          arrange(desc(max_rho)) %>% head(20)

        p_sp <- ggplot(species_summary,
                       aes(x = reorder(Species, mean_rho), y = mean_rho, fill = Role)) +
          geom_bar(stat = "identity") +
          coord_flip() +
          scale_fill_manual(values = c("Degradation" = "#377EB8", "Production" = "#FF7F00")) +
          labs(title = paste0("PRJEB6070_CRC_AdenomatousPolyps: Species Associated with Sarcosine Metabolism"),
               subtitle = "Spearman correlation (FDR<0.05, |rho|>0.3)",
               x = "", y = "Mean Spearman rho") +
          theme_bw() +
          theme(plot.title = element_text(size = 10, face = "bold"),
                axis.text.y = element_text(size = 7), legend.position = "bottom")

        .ggsave("results_sarcosine/sarcosine_associated_species.png", p_sp,
               width = 10, height = max(5, nrow(species_summary) * 0.3 + 2), dpi = 200)
        cat("Sarcosine associated species barplot saved.\n")
      }
    }
  }

  #############################################################################
  # PART 4b: SARCOSINE-ASSOCIATED SPECIES WITH PREVALENCE FILTER (>=10%)
  #############################################################################
  cat("\n\n========== PART 4b: SARCOSINE SPECIES CORRELATION (FILTERED) ==========\n")

  if (exists("cor_results") && nrow(cor_results) > 0 && exists("sig_cors") && nrow(sig_cors) > 0) {
    # Read filtered CSV to get KOs matching Part 3b (prevalence + Wilcoxon significance)
    filt_csv <- "results_sarcosine/sarcosine_KO_comparison_filtered.csv"
    if (file.exists(filt_csv)) {
      filt_data <- read.csv(filt_csv, stringsAsFactors = FALSE)
      pass_kos <- filt_data$KO[filt_data$Pass_Prevalence == TRUE &
                                !is.na(filt_data$p_value) &
                                filt_data$p_value < 0.05]
    } else {
      # Fallback: compute prevalence + Wilcoxon inline
      n_total_samples <- length(common_samples)
      ko_prevalence <- sapply(available_kos, function(k) {
        if (k %in% colnames(sarc_common)) sum(sarc_common[, k] > 0) / n_total_samples * 100 else 0
      })
      prev_kos <- names(ko_prevalence[ko_prevalence >= PREV_THRESHOLD])
      pass_kos <- c()
      for (k in prev_kos) {
        if (k %in% colnames(sarc_mat)) {
          h_vals <- sarc_mat[rownames(sarc_mat) %in% meta_hc$Run.ID[meta_hc$Group == "Healthy"], k]
          c_vals <- sarc_mat[rownames(sarc_mat) %in% meta_hc$Run.ID[meta_hc$Group == "Cancer"], k]
          wt <- suppressWarnings(wilcox.test(h_vals, c_vals))
          if (!is.na(wt$p.value) && wt$p.value < 0.05) pass_kos <- c(pass_kos, k)
        }
      }
    }
    sig_cors_filt <- sig_cors %>% filter(KO %in% pass_kos)

    cat("KOs passing prevalence + significance:", length(pass_kos), "\n")
    cat("Significant correlations after filter:", nrow(sig_cors_filt), "\n")

    if (nrow(sig_cors_filt) > 0) {
      filt_species_list <- unique(sig_cors_filt$Species_full)
      filt_ko_list <- unique(sig_cors_filt$KO)

      if (length(filt_species_list) > 1 & length(filt_ko_list) >= 1) {
        # Top species by max |rho|
        filt_sp_max_rho <- sig_cors_filt %>%
          group_by(Species_full) %>%
          summarise(max_rho = max(abs(rho)), .groups = "drop") %>%
          arrange(desc(max_rho)) %>% head(30)

        filt_top_sp_full <- filt_sp_max_rho$Species_full
        filt_top_sp_names <- sapply(filt_top_sp_full, get_species_name)

        # Build correlation matrix
        filt_cor_matrix <- matrix(0, nrow = length(filt_top_sp_full), ncol = length(filt_ko_list))
        rownames(filt_cor_matrix) <- filt_top_sp_names
        colnames(filt_cor_matrix) <- filt_ko_list

        filt_pval_matrix <- matrix(NA, nrow = length(filt_top_sp_full), ncol = length(filt_ko_list))
        rownames(filt_pval_matrix) <- filt_top_sp_names
        colnames(filt_pval_matrix) <- filt_ko_list

        for (r in seq_len(nrow(sig_cors_filt))) {
          sp_name <- get_species_name(sig_cors_filt$Species_full[r])
          ko_name <- sig_cors_filt$KO[r]
          if (sp_name %in% rownames(filt_cor_matrix) & ko_name %in% colnames(filt_cor_matrix)) {
            filt_cor_matrix[sp_name, ko_name] <- sig_cors_filt$rho[r]
            filt_pval_matrix[sp_name, ko_name] <- sig_cors_filt$p_adj[r]
          }
        }

        # Significance labels
        filt_sig_labels <- matrix("", nrow = nrow(filt_cor_matrix), ncol = ncol(filt_cor_matrix))
        rownames(filt_sig_labels) <- rownames(filt_cor_matrix)
        colnames(filt_sig_labels) <- colnames(filt_cor_matrix)
        for (i in seq_len(nrow(filt_pval_matrix))) {
          for (j in seq_len(ncol(filt_pval_matrix))) {
            p <- filt_pval_matrix[i, j]
            if (!is.na(p)) {
              if (p < 0.001)      filt_sig_labels[i, j] <- "***"
              else if (p < 0.01)  filt_sig_labels[i, j] <- "**"
              else if (p < 0.05)  filt_sig_labels[i, j] <- "*"
            }
          }
        }

        # Enrichment mapping
        filt_sp_enrichment <- sapply(filt_top_sp_full, function(sp) {
          if (sp %in% names(enrichment_map)) return(enrichment_map[sp])
          return("Unknown")
        })
        names(filt_sp_enrichment) <- filt_top_sp_names

        filt_row_ann <- data.frame(
          `Enriched in` = filt_sp_enrichment[rownames(filt_cor_matrix)],
          check.names = FALSE, stringsAsFactors = FALSE
        )
        rownames(filt_row_ann) <- rownames(filt_cor_matrix)

        filt_healthy_sp <- rownames(filt_row_ann)[filt_row_ann$`Enriched in` == "Healthy"]
        filt_cancer_sp  <- rownames(filt_row_ann)[filt_row_ann$`Enriched in` == "Cancer"]

        filt_healthy_ordered <- order_within_group(filt_healthy_sp, filt_cor_matrix)
        filt_cancer_ordered  <- order_within_group(filt_cancer_sp, filt_cor_matrix)
        filt_final_order <- c(filt_healthy_ordered, filt_cancer_ordered)

        filt_cor_ordered  <- filt_cor_matrix[filt_final_order, , drop = FALSE]
        filt_sig_ordered  <- filt_sig_labels[filt_final_order, , drop = FALSE]
        filt_row_ordered  <- filt_row_ann[filt_final_order, , drop = FALSE]

        filt_ko_labels <- sapply(colnames(filt_cor_ordered), function(k) {
          info <- sarc_info %>% filter(KO == k)
          if (nrow(info) > 0) {
            paste0(k, "\n", info$Short[1], " (", info$EC[1], ")\n", info$Role[1])
          } else k
        })
        colnames(filt_cor_ordered) <- filt_ko_labels
        colnames(filt_sig_ordered) <- filt_ko_labels

        filt_ann_colors <- list(`Enriched in` = c("Healthy" = "#1B7837", "Cancer" = "#B2182B"))
        filt_n_healthy_sp <- length(filt_healthy_ordered)
        filt_gaps_row <- if (filt_n_healthy_sp > 0 & filt_n_healthy_sp < nrow(filt_cor_ordered)) filt_n_healthy_sp else NULL

        filt_n_sp <- nrow(filt_cor_ordered)
        filt_n_ko <- ncol(filt_cor_ordered)
        filt_cell_w <- max(90, 130 - filt_n_ko * 5)
        filt_cell_h <- max(14, 20 - filt_n_sp * 0.15)

        pheatmap(filt_cor_ordered,
                 filename = "results_sarcosine/sarcosine_species_correlation_heatmap_filtered.png",
                 color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                 cluster_rows = FALSE, cluster_cols = FALSE,
                 gaps_row = filt_gaps_row,
                 annotation_row = filt_row_ordered,
                 annotation_colors = filt_ann_colors,
                 display_numbers = filt_sig_ordered,
                 number_color = "black", fontsize_number = 10,
                 fontsize_row = 7, fontsize_col = 8,
                 cellwidth = filt_cell_w, cellheight = filt_cell_h,
                 width = max(10, filt_n_ko * 1.5 + 5),
                 height = max(7, filt_n_sp * 0.35 + 4),
                 main = paste0("PRJEB6070_CRC_AdenomatousPolyps: Species-Sarcosine KO Correlation\n",
                               "(Healthy vs ", cancer_label,
                               ")  [* p<0.05, ** p<0.01, *** p<0.001]\n",
                               "(Prevalence >= ", PREV_THRESHOLD, "% + Wilcoxon p<0.05)"),
                 angle_col = 0)
        cat("Filtered species correlation heatmap saved.\n")
      } else {
        cat("Not enough species/KOs for filtered heatmap.\n")
      }
    } else {
      cat("No significant correlations after prevalence filter.\n")
    }
  } else {
    cat("No correlation data available for filtering.\n")
  }

  #############################################################################
  # PART 5: PRODUCTION vs DEGRADATION PATHWAY BALANCE
  #############################################################################
  cat("\n\n========== PART 5: PRODUCTION vs DEGRADATION PATHWAY ==========\n")

  deg_kos <- intersect(sarcosine_kos$KO[sarcosine_kos$Role == "Degradation"], colnames(ko_mat))
  prod_kos <- intersect(sarcosine_kos$KO[sarcosine_kos$Role == "Production"], colnames(ko_mat))

  cat("Degradation KOs available:", length(deg_kos), "(", paste(deg_kos, collapse=", "), ")\n")
  cat("Production KOs available:", length(prod_kos), "(", paste(prod_kos, collapse=", "), ")\n")

  if (length(deg_kos) > 0 & length(prod_kos) > 0) {
    deg_sum <- rowSums(ko_mat[, deg_kos, drop = FALSE])
    prod_sum <- rowSums(ko_mat[, prod_kos, drop = FALSE])
    pseudo <- 1e-8
    log2_ratio <- log2((prod_sum + pseudo) / (deg_sum + pseudo))

    ratio_df <- data.frame(
      Run.ID = rownames(ko_mat),
      Group = meta_ko$Group,
      Degradation_sum = deg_sum,
      Production_sum = prod_sum,
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

    cat("Degradation sum:", format_pval(deg_test$p.value), "\n")
    cat("Production sum:", format_pval(prod_test$p.value), "\n")
    cat("log2(Prod/Deg) ratio:", format_pval(ratio_test$p.value), "\n")

    # Boxplot figure
    plot_data <- ratio_df %>%
      select(Group, Degradation_sum, Production_sum, Log2_Prod_Deg_Ratio) %>%
      pivot_longer(cols = c(Degradation_sum, Production_sum, Log2_Prod_Deg_Ratio),
                   names_to = "Metric", values_to = "Value") %>%
      mutate(Metric = factor(Metric,
                             levels = c("Degradation_sum","Production_sum","Log2_Prod_Deg_Ratio"),
                             labels = c(paste0("Degradation\n(sum of ", length(deg_kos), " KOs)"),
                                       paste0("Production\n(sum of ", length(prod_kos), " KOs)"),
                                       "log2(Prod/Deg)\nRatio")))

    pval_labels_pd <- data.frame(
      Metric = factor(c(paste0("Degradation\n(sum of ", length(deg_kos), " KOs)"),
                       paste0("Production\n(sum of ", length(prod_kos), " KOs)"),
                       "log2(Prod/Deg)\nRatio"),
                     levels = c(paste0("Degradation\n(sum of ", length(deg_kos), " KOs)"),
                               paste0("Production\n(sum of ", length(prod_kos), " KOs)"),
                               "log2(Prod/Deg)\nRatio")),
      pval_text = c(
        format_pval(deg_test$p.value),
        format_pval(prod_test$p.value),
        format_pval(ratio_test$p.value)
      ),
      stringsAsFactors = FALSE
    )

    ypos_pd <- plot_data %>%
      group_by(Metric) %>%
      summarise(ymax = max(Value, na.rm = TRUE), .groups = "drop")
    pval_labels_pd <- merge(pval_labels_pd, ypos_pd, by = "Metric")

    p_ratio <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot(outlier.shape = 21, alpha = 0.7) +
      geom_jitter(width = 0.15, size = 0.6, alpha = 0.3) +
      facet_wrap(~Metric, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
      geom_text(data = pval_labels_pd, aes(x = 1.5, y = ymax * 1.1, label = pval_text),
                inherit.aes = FALSE, size = 3.5, fontface = "italic") +
      labs(title = paste0("PRJEB6070_CRC_AdenomatousPolyps: Sarcosine Metabolism Pathway Balance"),
           subtitle = paste0("Healthy vs ", cancer_label,
                            "  |  Deg KOs: ", paste(deg_kos, collapse=","),
                            "  |  Prod KOs: ", paste(prod_kos, collapse=",")),
           x = "", y = "Abundance / Ratio") +
      theme_bw() +
      theme(plot.title = element_text(size = 10, face = "bold"),
            plot.subtitle = element_text(size = 7),
            strip.text = element_text(size = 9),
            legend.position = "bottom")

    .ggsave("results_sarcosine/sarcosine_prod_vs_deg.png", p_ratio,
           width = 12, height = 5.5, dpi = 200)
    cat("Production vs Degradation figure saved.\n")
  } else {
    cat("Cannot compute Prod/Deg ratio (need both pathway KOs available).\n")
  }
}

################################################################################
# Part 6: Sarcosine-Associated Bacteria Differential Abundance
################################################################################
cat("\n\n========================================\n")
cat("Part 6: Sarcosine-Associated Bacteria Differential Abundance\n")
cat("========================================\n\n")

sp_corr_file <- file.path("results_sarcosine", "sarcosine_species_correlation.csv")
if (file.exists(sp_corr_file)) {
  sp_corr <- read.csv(sp_corr_file, stringsAsFactors = FALSE)
  sig_sp <- sp_corr[sp_corr$p_adj < 0.05 & abs(sp_corr$rho) > 0.3, ]

  if (nrow(sig_sp) > 0) {
    # Get top species by max |rho|
    sp_max_rho <- sig_sp %>%
      group_by(Species, Species_full) %>%
      summarise(
        max_abs_rho = max(abs(rho)),
        roles = paste(sort(unique(Role)), collapse = " & "),
        n_sig_kos = n_distinct(KO),
        .groups = "drop"
      ) %>%
      arrange(desc(max_abs_rho))

    n_top <- min(20, nrow(sp_max_rho))
    top_species_df <- sp_max_rho[1:n_top, ]

    # Read bacteria abundance
    bact_file <- list.files(".", pattern = "^Bacteria_.*\\.txt$", full.names = TRUE)[1]
    bact <- read.delim(bact_file, stringsAsFactors = FALSE)
    colnames(bact) <- c("Taxa", "Run_ID", "Abundance")

    # Species-level aggregation
    bact_species <- bact %>%
      filter(grepl("\\|s__", Taxa)) %>%
      mutate(Species_taxa = sub("\\|t__.*$", "", Taxa)) %>%
      group_by(Species_taxa, Run_ID) %>%
      summarise(Abundance = sum(Abundance), .groups = "drop")

    bact_target <- bact_species %>%
      filter(Species_taxa %in% top_species_df$Species_full)

    if (nrow(bact_target) > 0) {
      runs_hc <- meta_hc$Run.ID
      diff_results <- data.frame()

      for (i in seq_len(nrow(top_species_df))) {
        sp_taxa <- top_species_df$Species_full[i]
        sp_name <- top_species_df$Species[i]
        sp_roles <- top_species_df$roles[i]

        sp_abund <- bact_target %>%
          filter(Species_taxa == sp_taxa, Run_ID %in% runs_hc)

        abund_vec <- setNames(rep(0, length(runs_hc)), runs_hc)
        if (nrow(sp_abund) > 0) {
          for (j in seq_len(nrow(sp_abund))) {
            abund_vec[sp_abund$Run_ID[j]] <- sp_abund$Abundance[j]
          }
        }

        healthy_runs <- meta_hc$Run.ID[meta_hc$Group == "Healthy"]
        cancer_runs  <- meta_hc$Run.ID[meta_hc$Group == "Cancer"]
        healthy_vals <- abund_vec[healthy_runs]
        cancer_vals  <- abund_vec[cancer_runs]

        mean_healthy <- mean(healthy_vals, na.rm = TRUE)
        mean_cancer  <- mean(cancer_vals, na.rm = TRUE)
        pseudo_fc <- 1e-6
        log2fc <- log2((mean_cancer + pseudo_fc) / (mean_healthy + pseudo_fc))

        wt_p <- tryCatch({
          wilcox.test(cancer_vals, healthy_vals)$p.value
        }, error = function(e) NA)

        diff_results <- rbind(diff_results, data.frame(
          Species = sp_name, Species_full = sp_taxa, Role = sp_roles,
          Mean_Healthy = mean_healthy, Mean_Cancer = mean_cancer,
          log2FC = log2fc, p_value = wt_p,
          Direction = ifelse(log2fc > 0, "Cancer-enriched", "Healthy-enriched"),
          max_abs_rho = top_species_df$max_abs_rho[i],
          stringsAsFactors = FALSE))
      }

      diff_results$p_adj <- p.adjust(diff_results$p_value, method = "BH")
      write.csv(diff_results, "results_sarcosine/sarcosine_bacteria_diff_abundance.csv",
                row.names = FALSE)

      cat("Significant species (p_adj<0.05):", sum(diff_results$p_adj < 0.05, na.rm=TRUE), "\n")

      # Plot
      diff_results$Species_display <- gsub("_", " ", diff_results$Species)
      diff_results <- diff_results %>% arrange(log2FC)
      diff_results$Species_display <- factor(diff_results$Species_display,
                                              levels = diff_results$Species_display)

      diff_results$p_text <- sapply(seq_len(nrow(diff_results)), function(i) {
        p <- diff_results$p_adj[i]
        if (is.na(p) || p >= 0.05) return("")
        format_pval(p)
      })
      diff_results$Role_label <- ifelse(
        diff_results$Role == "Degradation", "Deg",
        ifelse(diff_results$Role == "Production", "Prod", "Deg & Prod"))

      x_max <- max(abs(diff_results$log2FC), na.rm = TRUE) * 1.35
      n_healthy <- sum(meta_hc$Group == "Healthy")
      n_cancer  <- sum(meta_hc$Group == "Cancer")
      role_colors <- c("Deg" = "#2196F3", "Prod" = "#FF9800", "Deg & Prod" = "#9C27B0")

      p_diffab <- ggplot(diff_results, aes(x = log2FC, y = Species_display)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
        geom_col(aes(fill = Direction), width = 0.7, alpha = 0.85) +
        geom_text(data = diff_results[diff_results$p_text != "", ],
          aes(x = ifelse(log2FC > 0, log2FC + x_max * 0.03, log2FC - x_max * 0.03),
              label = p_text),
          hjust = ifelse(diff_results$log2FC[diff_results$p_text != ""] > 0, 0, 1),
          size = 2.8, color = "gray20", fontface = "italic") +
        geom_point(aes(x = x_max * 0.92, color = Role_label), size = 3, shape = 15) +
        geom_text(aes(x = x_max * 0.97, label = Role_label, color = Role_label),
                  size = 2.5, hjust = 0, fontface = "bold") +
        scale_fill_manual(values = c("Cancer-enriched" = "#B2182B", "Healthy-enriched" = "#1B7837"),
                          name = "Direction") +
        scale_color_manual(values = role_colors, name = "Sarcosine\nPathway") +
        scale_x_continuous(limits = c(-x_max, x_max * 1.15), expand = c(0, 0)) +
        labs(title = "Sarcosine-Associated Bacteria: Differential Abundance",
             subtitle = paste0("PRJEB6070_CRC_AdenomatousPolyps | Healthy (n=", n_healthy, ") vs Cancer (n=", n_cancer, ")"),
             x = expression(log[2]~"Fold Change (Cancer / Healthy)"), y = NULL,
             caption = "Top species by |rho| from sarcosine KO correlation analysis\nBars: log2FC of mean abundance | p-values: BH-adjusted") +
        theme_minimal(base_size = 11) +
        theme(plot.title = element_text(face = "bold", size = 13, hjust = 0),
              plot.subtitle = element_text(size = 10, color = "gray40"),
              plot.caption = element_text(size = 8, color = "gray50", hjust = 0),
              axis.text.y = element_text(face = "italic", size = 9),
              panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
              legend.position = "bottom", legend.box = "horizontal",
              legend.title = element_text(size = 9, face = "bold"),
              plot.margin = margin(10, 15, 10, 10)) +
        guides(fill = guide_legend(order = 1), color = guide_legend(order = 2))

      fig_h <- max(6, n_top * 0.35 + 2)
      .ggsave("results_sarcosine/sarcosine_bacteria_diff_abundance.png",
             p_diffab, width = 10, height = fig_h, dpi = 200, bg = "white")
      cat("Sarcosine bacteria diff. abundance figure saved.\n")
    } else {
      cat("No matching bacteria in abundance data.\n")
    }
  } else {
    cat("No significant species-KO correlations for diff. abundance analysis.\n")
  }
} else {
  cat("Species correlation file not found. Run Part 4 first.\n")
}

cat("\n\n=== Analysis Complete for PRJEB6070_CRC_AdenomatousPolyps ===\n")
cat("Results saved in: results_bacteria/, results_kegg/, results_sarcosine/\n")

