###############################################################################
# Integrated Pooled Analysis: 4 CRC Cohorts (Healthy vs Cancer)
# HGMT Database, WGS-based
#
# Cohorts:
#   PRJEB6070  (France+Germany) - has 16S + Adenomatous Polyps to exclude
#   PRJNA429097 (China)
#   PRJEB10878  (China)
#   PRJEB27928  (Germany)
#
# Analysis pipeline:
#   Section 0: Setup
#   Section 1: Data loading & merging (metadata, bacteria, KO)
#   Part 1: Bacteria analysis (alpha/beta diversity, differential abundance)
#   Part 2: KEGG KO differential abundance
#   Part 3: Sarcosine KO analysis (+ Part 3b: prevalence filtered)
#   Part 4: Sarcosine-associated species correlation (+ Part 4b: filtered)
#   Part 5: Production vs Degradation pathway balance
#   Part 6: Sarcosine-associated bacteria differential abundance
#   Supplementary: Per-cohort consistency verification
###############################################################################

###############################################################################
# Section 0: SETUP
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

# Fix: cairo DLL not available on this system - use quartz backend for PNG
options(bitmapType = "quartz")

BASE_DIR <- "Healthy_vs_Cancer_4_CRC_cohorts_integrated"   # relative to the analysis-folder root (working dir)
setwd(BASE_DIR)

# Cohort directory mapping (each has uniquely timestamped files inside)
COHORT_DIRS <- list(
  PRJEB6070  = file.path(BASE_DIR, "PRJEB6070_CRC_AdenomatousPolyps"),
  PRJNA429097 = file.path(BASE_DIR, "PRJNA429097_CRC"),
  PRJEB10878  = file.path(BASE_DIR, "PRJEB10878_CRC"),
  PRJEB27928  = file.path(BASE_DIR, "PRJEB27928_CRC")
)

# Sarcosine-related KOs (verified from KEGG)
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

PREV_THRESHOLD <- 10  # percent

# Helper functions
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

# Create output directories
dir.create("results_integrated/bacteria", recursive = TRUE, showWarnings = FALSE)
dir.create("results_integrated/kegg", recursive = TRUE, showWarnings = FALSE)
dir.create("results_integrated/sarcosine", recursive = TRUE, showWarnings = FALSE)
dir.create("results_integrated/supplementary", recursive = TRUE, showWarnings = FALSE)

cat("=== Integrated Pooled Analysis: 4 CRC Cohorts ===\n\n")

###############################################################################
# Section 1: DATA LOADING & MERGING
###############################################################################
cat("========== SECTION 1: DATA LOADING ==========\n\n")

# --- 1.1 Metadata ---
cat("--- 1.1 Loading metadata ---\n")
meta_pooled <- data.frame()

for (cohort_id in names(COHORT_DIRS)) {
  cohort_path <- COHORT_DIRS[[cohort_id]]
  meta_file <- list.files(cohort_path, pattern = "^selected_project_.*\\.txt$",
                          full.names = TRUE)
  # Use the last (most recent) file if multiple exist
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

  # Auto-detect cancer phenotype
  all_pheno <- unique(meta$Phenotype.name)
  cancer_pheno <- all_pheno[!all_pheno %in% c("Health", "Adenomatous Polyps")]
  if (length(cancer_pheno) == 0) {
    cat("    WARNING: No cancer phenotype found. Skipping.\n")
    next
  }
  cancer_label <- cancer_pheno[1]

  # Filter: Health + Cancer only (exclude Adenomatous Polyps)
  meta_hc <- meta %>%
    filter(Phenotype.name %in% c("Health", cancer_label)) %>%
    mutate(Group = ifelse(Phenotype.name == "Health", "Healthy", "Cancer"))

  # Filter: WGS only (exclude 16S)
  if ("Assay.type" %in% colnames(meta_hc)) {
    n_before <- nrow(meta_hc)
    meta_hc <- meta_hc %>% filter(Assay.type == "WGS")
    cat("    Assay filter: kept", nrow(meta_hc), "WGS out of", n_before, "total\n")
  }

  meta_hc$Cohort <- cohort_id
  meta_pooled <- rbind(meta_pooled, meta_hc)
  cat("    Healthy:", sum(meta_hc$Group == "Healthy"),
      ", Cancer:", sum(meta_hc$Group == "Cancer"), "\n")
}

cat("\n--- Pooled metadata summary ---\n")
print(table(meta_pooled$Cohort, meta_pooled$Group))
cat("Total samples:", nrow(meta_pooled), "\n\n")

# Verification: no 16S, no Adenomatous Polyps
stopifnot(all(meta_pooled$Assay.type == "WGS" | is.na(meta_pooled$Assay.type)))
stopifnot(!any(meta_pooled$Phenotype.name == "Adenomatous Polyps"))
cat("[VERIFIED] No 16S samples, no Adenomatous Polyps.\n\n")

# --- 1.2 Bacteria data ---
cat("--- 1.2 Loading bacteria data ---\n")
bact_all <- data.frame()

for (cohort_id in names(COHORT_DIRS)) {
  cohort_path <- COHORT_DIRS[[cohort_id]]
  bact_file <- list.files(cohort_path, pattern = "^Bacteria_.*\\.txt$",
                          full.names = TRUE)[1]
  cat("  Loading:", cohort_id, "from", basename(bact_file), "\n")

  bact <- read.delim(bact_file, stringsAsFactors = FALSE)
  colnames(bact) <- trimws(colnames(bact))

  # Species level only, filter to pooled samples
  cohort_runs <- meta_pooled$Run.ID[meta_pooled$Cohort == cohort_id]
  bact_species <- bact %>%
    filter(grepl("\\|s__", Taxa) & !grepl("\\|t__", Taxa)) %>%
    filter(Run.ID %in% cohort_runs)

  bact_all <- rbind(bact_all, bact_species)
  cat("    Rows:", nrow(bact_species), "(species-level, filtered to WGS Healthy/Cancer)\n")
  rm(bact, bact_species); gc(verbose = FALSE)
}

cat("  Total bacteria rows:", nrow(bact_all), "\n")
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

cat("  Bacteria matrix:", nrow(bact_mat), "samples x", ncol(bact_mat), "species\n\n")
rm(bact_all, bact_wide); gc(verbose = FALSE)

# --- 1.3 KO data (sequential loading for memory) ---
cat("--- 1.3 Loading KO data ---\n")
ko_all <- data.frame()

for (cohort_id in names(COHORT_DIRS)) {
  cohort_path <- COHORT_DIRS[[cohort_id]]
  ko_file <- file.path(cohort_path, "KO_relative_abundance.tsv")
  cat("  Loading:", cohort_id, "...\n")

  ko_raw <- fromJSON(ko_file)
  cohort_runs <- meta_pooled$Run.ID[meta_pooled$Cohort == cohort_id]
  ko_filtered <- ko_raw %>%
    filter(run_id %in% cohort_runs) %>%
    select(ko, run_id, abundance)

  ko_all <- rbind(ko_all, ko_filtered)
  cat("    Rows:", nrow(ko_filtered), "\n")
  rm(ko_raw, ko_filtered); gc(verbose = FALSE)
}

cat("  Total KO rows:", nrow(ko_all), "\n")
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

cat("  KO matrix:", nrow(ko_mat), "samples x", ncol(ko_mat), "KOs\n\n")
rm(ko_all, ko_wide); gc(verbose = FALSE)

# Group indices for KO analysis
h_idx <- which(meta_ko$Group == "Healthy")
c_idx <- which(meta_ko$Group == "Cancer")

###############################################################################
# PART 1: BACTERIA ANALYSIS
###############################################################################
cat("\n\n========== PART 1: BACTERIA ANALYSIS ==========\n")

# --- 1a. Alpha Diversity ---
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

cat("Shannon p:", format_pval(shannon_test$p.value), "\n")
cat("Simpson p:", format_pval(simpson_test$p.value), "\n")
cat("Richness p:", format_pval(richness_test$p.value), "\n")

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
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.15) +
  facet_wrap(~Metric, scales = "free_y") +
  scale_fill_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
  geom_text(data = pval_labels, aes(x = 1.5, y = ymax * 1.1, label = label),
            inherit.aes = FALSE, size = 3.5) +
  labs(title = paste0("Integrated 4 CRC Cohorts: Alpha Diversity (n=", nrow(bact_mat), ")"),
       subtitle = "Healthy vs Colorectal Neoplasms (Pooled)",
       x = "", y = "Value") +
  theme_bw() +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 10, face = "bold"),
        strip.text = element_text(size = 10))

ggsave("results_integrated/bacteria/alpha_diversity_pooled.png", p_alpha,
        width = 10, height = 5, dpi = 200)
cat("Alpha diversity plot saved.\n")

# --- 1b. Beta Diversity ---
cat("\n--- Beta Diversity ---\n")
bc_dist <- vegdist(bact_mat, method = "bray")
pcoa_res <- cmdscale(bc_dist, k = 2, eig = TRUE)
eig_pct <- round(pcoa_res$eig / sum(pcoa_res$eig[pcoa_res$eig > 0]) * 100, 1)

pcoa_df <- data.frame(
  PC1 = pcoa_res$points[, 1],
  PC2 = pcoa_res$points[, 2],
  Group = meta_matched$Group,
  Cohort = meta_matched$Cohort
)

set.seed(42)
perm_simple <- adonis2(bc_dist ~ Group, data = meta_matched, permutations = 999)
perm_cohort <- adonis2(bc_dist ~ Group + Cohort, data = meta_matched, permutations = 999)

perm_p1 <- perm_simple$`Pr(>F)`[1]
perm_r2_1 <- round(perm_simple$R2[1], 4)
perm_p2 <- perm_cohort$`Pr(>F)`[1]
perm_r2_2 <- round(perm_cohort$R2[1], 4)

cat("PERMANOVA (Group only): R2=", perm_r2_1, ", p=", format_pval(perm_p1), "\n")
cat("PERMANOVA (Group+Cohort): Group R2=", perm_r2_2, ", p=", format_pval(perm_p2), "\n")

p_beta <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = Group, shape = Cohort)) +
  geom_point(size = 1.5, alpha = 0.5) +
  stat_ellipse(aes(group = Group), level = 0.95, linetype = 2) +
  scale_color_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
  scale_shape_manual(values = c("PRJEB6070" = 16, "PRJNA429097" = 17,
                                "PRJEB10878" = 15, "PRJEB27928" = 18)) +
  labs(title = paste0("Integrated 4 CRC Cohorts: PCoA (Bray-Curtis, n=", nrow(bact_mat), ")"),
       x = paste0("PCoA1 (", eig_pct[1], "%)"),
       y = paste0("PCoA2 (", eig_pct[2], "%)"),
       subtitle = paste0("PERMANOVA (Group): R2=", perm_r2_1, ", ", format_pval(perm_p1),
                         "  |  PERMANOVA (Group+Cohort): Group R2=", perm_r2_2,
                         ", ", format_pval(perm_p2))) +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7))

ggsave("results_integrated/bacteria/beta_diversity_pcoa_pooled.png", p_beta,
        width = 10, height = 7, dpi = 200)
cat("Beta diversity plot saved.\n")

# --- 1c. Differential Abundance (Species) ---
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

write.csv(diff_results, "results_integrated/bacteria/diff_abundance_species_pooled.csv",
          row.names = FALSE)

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
  labs(title = paste0("Integrated 4 CRC Cohorts: Differential Species (n=", nrow(bact_mat), ")"),
       subtitle = "Healthy vs Colorectal Neoplasms (Pooled)",
       x = "log2FC (Cancer/Healthy)", y = "-log10(adj. p-value)") +
  theme_bw() +
  theme(plot.title = element_text(size = 10, face = "bold"),
        legend.position = "bottom")

ggsave("results_integrated/bacteria/volcano_species_pooled.png", p_volcano,
        width = 10, height = 7, dpi = 200)
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

  png("results_integrated/bacteria/heatmap_top_species_pooled.png",
      width = 12, height = max(8, top_n * 0.25 + 2), units = "in", res = 200)
  pheatmap(t(hm_scaled),
           annotation_col = annotation_row,
           annotation_colors = ann_colors,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           cluster_rows = TRUE, cluster_cols = TRUE,
           show_colnames = FALSE, fontsize_row = 7,
           main = paste0("Integrated 4 CRC Cohorts: Top Differential Species (Z-score, n=",
                         nrow(bact_mat), ")"))
  dev.off()
  cat("Heatmap saved.\n")
}

###############################################################################
# PART 2: KEGG KO ANALYSIS
###############################################################################
cat("\n\n========== PART 2: KEGG KO ANALYSIS ==========\n")

ko_prev <- colSums(ko_mat > 0) / nrow(ko_mat)
ko_filt <- ko_mat[, ko_prev >= 0.1, drop = FALSE]
cat("KOs after prevalence filter:", ncol(ko_filt), "\n")

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

write.csv(ko_diff, "results_integrated/kegg/diff_KO_abundance_pooled.csv", row.names = FALSE)

sig_kos <- ko_diff %>% filter(p_adj < 0.05)
cat("Significant KOs (FDR < 0.05):", nrow(sig_kos), "\n")

# KO Volcano — paper style: larger fonts; label ONLY significant sarcosine KOs.
# (kept in sync with scripts/regenerate_volcano_KO_paper.R)
ko_diff$Significance <- "Not Significant"
ko_diff$Significance[ko_diff$p_adj < 0.05 & ko_diff$log2FC > 0.5] <- "Enriched in Cancer"
ko_diff$Significance[ko_diff$p_adj < 0.05 & ko_diff$log2FC < -0.5] <- "Enriched in Healthy"

# Significant sarcosine KOs present on the volcano (prevalence>=10% & p_adj<0.05)
sarc_lab <- ko_diff %>%
  inner_join(sarcosine_kos[, c("KO", "Short", "Role")], by = "KO") %>%
  filter(p_adj < 0.05) %>%
  mutate(fillcol = ifelse(log2FC > 0, "#B2182B", "#1B7837"),
         lab = paste0(KO, " — ", Short, "\n(", Role, ")"),
         nx  = ifelse(log2FC > 0, 1.0, -1.0))

p_ko_volcano <- ggplot(ko_diff, aes(x = log2FC, y = -log10(p_adj))) +
  geom_point(aes(color = Significance), alpha = 0.6, size = 1.6) +
  scale_color_manual(values = c("Enriched in Cancer" = "#B2182B",
                                 "Enriched in Healthy" = "#1B7837",
                                 "Not Significant" = "grey70")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  # label only the significant sarcosine KOs — segment + text, no marker
  geom_text_repel(data = sarc_lab, aes(label = lab),
                  color = "black", fontface = "bold", size = 5,
                  bg.color = "white", bg.r = 0.12, lineheight = 0.9,
                  box.padding = 1.0, point.padding = 0.6,
                  segment.color = "grey20", segment.size = 0.6,
                  min.segment.length = 0, max.overlaps = Inf,
                  nudge_x = sarc_lab$nx, nudge_y = 4, show.legend = FALSE) +
  labs(title = paste0("Integrated 4 CRC Cohorts: Differential KO Abundance (n=", nrow(ko_mat), ")"),
       subtitle = "Healthy vs Colorectal Neoplasms (Pooled) — significant sarcosine KOs highlighted",
       x = "log2FC (Cancer/Healthy)", y = "-log10(adj. p-value)") +
  theme_bw() +
  theme(plot.title    = element_text(size = 20, face = "bold"),
        plot.subtitle = element_text(size = 15),
        axis.title    = element_text(size = 17),
        axis.text     = element_text(size = 14),
        legend.title  = element_text(size = 15),
        legend.text   = element_text(size = 14),
        legend.position = "bottom") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

ggsave("results_integrated/kegg/volcano_KO_pooled.png", p_ko_volcano,
        width = 12, height = 8.5, dpi = 300, bg = "white")
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
  write.csv(sarc_results, "results_integrated/sarcosine/sarcosine_KO_comparison_pooled.csv",
            row.names = FALSE)

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
      labs(title = paste0("Integrated 4 CRC Cohorts: Sarcosine Enzyme Gene Abundance (n=",
                          nrow(ko_mat), ")"),
           subtitle = "Healthy vs Colorectal Neoplasms (significant only, p<0.05)",
           x = "", y = "Mean Relative Abundance") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            plot.title = element_text(size = 10, face = "bold"),
            legend.position = "bottom",
            plot.margin = margin(t = 10, r = 10, b = 20, l = 20))

    fig_w <- max(8, nrow(sarc_sig) * 1.5 + 2)
    ggsave("results_integrated/sarcosine/sarcosine_enzyme_barplot_pooled.png", p_sarc,
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
      geom_jitter(width = 0.15, size = 0.3, alpha = 0.15) +
      facet_wrap(~Label, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
      geom_text(data = box_pval, aes(x = 1.5, y = ymax * 1.15, label = pval_text),
                inherit.aes = FALSE, size = 3) +
      labs(title = paste0("Integrated 4 CRC Cohorts: Sarcosine Enzyme KOs (n=", nrow(ko_mat), ")"),
           subtitle = "Healthy vs Colorectal Neoplasms (p<0.05)",
           x = "", y = "Relative Abundance") +
      theme_bw() +
      theme(plot.title = element_text(size = 10, face = "bold"),
            legend.position = "bottom", strip.text = element_text(size = 7))

    box_nrow <- ceiling(nrow(sarc_sig) / 3)
    ggsave("results_integrated/sarcosine/sarcosine_enzyme_boxplot_pooled.png", p_sarc_box,
           width = 12, height = max(4, box_nrow * 3.5), dpi = 200)
    cat("Sarcosine enzyme plots saved.\n")
  }

  #############################################################################
  # PART 3b: SARCOSINE KO ANALYSIS WITH PREVALENCE FILTER (>=10%)
  #############################################################################
  cat("\n\n========== PART 3b: SARCOSINE KOs (PREVALENCE FILTERED) ==========\n")

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
            "results_integrated/sarcosine/sarcosine_KO_comparison_filtered_pooled.csv",
            row.names = FALSE)

  sarc_filtered <- sarc_results[sarc_results$Pass_Prevalence, ]
  sarc_filt_sig <- sarc_filtered[!is.na(sarc_filtered$p_value) & sarc_filtered$p_value < 0.05, ]

  cat("KOs passing prevalence (>=", PREV_THRESHOLD, "%):", nrow(sarc_filtered),
      "/", nrow(sarc_results), "\n")
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
      labs(title = paste0("Integrated 4 CRC Cohorts: Sarcosine Enzymes (Prevalence Filtered, n=",
                          nrow(ko_mat), ")"),
           subtitle = paste0("Prevalence >=", PREV_THRESHOLD, "% | Significant only (p<0.05)"),
           x = "", y = "Mean Relative Abundance") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            plot.title = element_text(size = 10, face = "bold"),
            legend.position = "bottom",
            plot.margin = margin(t = 10, r = 10, b = 20, l = 20))

    filt_fig_w <- max(8, nrow(sarc_filt_sig) * 1.5 + 2)
    ggsave("results_integrated/sarcosine/sarcosine_enzyme_barplot_filtered_pooled.png",
           p_filt_bar, width = min(filt_fig_w, 16), height = 7, dpi = 200)

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
      geom_jitter(width = 0.15, size = 0.3, alpha = 0.15) +
      facet_wrap(~Label, scales = "free_y", ncol = 3) +
      scale_fill_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
      geom_text(data = filt_box_pval, aes(x = 1.5, y = ymax * 1.15, label = pval_text),
                inherit.aes = FALSE, size = 3) +
      labs(title = paste0("Integrated 4 CRC Cohorts: Sarcosine Enzyme KOs (Prevalence Filtered, n=",
                          nrow(ko_mat), ")"),
           subtitle = paste0("Prevalence >=", PREV_THRESHOLD, "% | p<0.05"),
           x = "", y = "Relative Abundance") +
      theme_bw() +
      theme(plot.title = element_text(size = 10, face = "bold"),
            legend.position = "bottom", strip.text = element_text(size = 7))

    filt_box_nrow <- ceiling(nrow(sarc_filt_sig) / 3)
    ggsave("results_integrated/sarcosine/sarcosine_enzyme_boxplot_filtered_pooled.png",
           p_filt_box, width = 12, height = max(4, filt_box_nrow * 3.5), dpi = 200)
    cat("Filtered sarcosine enzyme plots saved.\n")
  } else {
    cat("No significant KOs after prevalence filter.\n")
  }

  #############################################################################
  # PART 4: SARCOSINE-ASSOCIATED SPECIES
  #############################################################################
  cat("\n\n========== PART 4: SARCOSINE-ASSOCIATED SPECIES ==========\n")

  common_samples <- intersect(rownames(bact_filt), rownames(sarc_mat))
  cat("Common samples for correlation:", length(common_samples), "\n")

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
      write.csv(cor_results,
                "results_integrated/sarcosine/sarcosine_species_correlation_pooled.csv",
                row.names = FALSE)

      sig_cors <- cor_results %>% filter(p_adj < 0.05 & abs(rho) > 0.3)
      cat("Significant species-KO correlations:", nrow(sig_cors), "\n")

      if (nrow(sig_cors) > 0) {
        sig_species_list <- unique(sig_cors$Species_full)
        sig_ko_list <- unique(sig_cors$KO)

        # Define order_within_group function (used by Part 4 and 4b)
        order_within_group <- function(species_names, mat) {
          if (length(species_names) <= 1) return(species_names)
          sub_mat <- mat[species_names, , drop = FALSE]
          d <- dist(sub_mat, method = "euclidean")
          hc <- hclust(d, method = "complete")
          return(species_names[hc$order])
        }

        # Compute enrichment for heatmap annotation (used by Part 4 and 4b)
        bact_long_species <- data.frame()
        for (cohort_id in names(COHORT_DIRS)) {
          cohort_path <- COHORT_DIRS[[cohort_id]]
          bact_file <- list.files(cohort_path, pattern = "^Bacteria_.*\\.txt$",
                                  full.names = TRUE)[1]
          bact_tmp <- read.delim(bact_file, stringsAsFactors = FALSE)
          colnames(bact_tmp) <- trimws(colnames(bact_tmp))
          cohort_runs <- meta_pooled$Run.ID[meta_pooled$Cohort == cohort_id]
          bact_tmp <- bact_tmp %>%
            filter(grepl("\\|s__", Taxa)) %>%
            filter(Run.ID %in% cohort_runs) %>%
            mutate(Species_taxa = sub("\\|t__.*$", "", Taxa)) %>%
            group_by(Species_taxa, Run.ID) %>%
            summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
          bact_long_species <- rbind(bact_long_species, bact_tmp)
          rm(bact_tmp); gc(verbose = FALSE)
        }

        bact_sp_wide <- bact_long_species %>%
          pivot_wider(names_from = Run.ID, values_from = Abundance, values_fill = 0)
        healthy_runs_enr <- meta_pooled$Run.ID[meta_pooled$Group == "Healthy"]
        cancer_runs_enr  <- meta_pooled$Run.ID[meta_pooled$Group == "Cancer"]
        h_cols <- intersect(healthy_runs_enr, colnames(bact_sp_wide))
        c_cols <- intersect(cancer_runs_enr, colnames(bact_sp_wide))
        sp_enrichment_df <- data.frame(
          Species_taxa = bact_sp_wide$Species_taxa,
          Healthy = if (length(h_cols) > 0)
            rowMeans(as.matrix(bact_sp_wide[, h_cols, drop = FALSE])) else 0,
          Cancer = if (length(c_cols) > 0)
            rowMeans(as.matrix(bact_sp_wide[, c_cols, drop = FALSE])) else 0,
          stringsAsFactors = FALSE
        ) %>% mutate(Enriched = ifelse(Cancer >= Healthy, "Cancer", "Healthy"))
        enrichment_map <- setNames(sp_enrichment_df$Enriched, sp_enrichment_df$Species_taxa)

        rm(bact_sp_wide, bact_long_species); gc(verbose = FALSE)

        if (length(sig_species_list) > 1 & length(sig_ko_list) >= 1) {
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
            sp_taxa <- sub("\\|t__.*$", "", sp)
            if (sp_taxa %in% names(enrichment_map)) return(enrichment_map[sp_taxa])
            if (sp %in% names(enrichment_map)) return(enrichment_map[sp])
            return("Unknown")
          })
          names(sp_enrichment) <- top_sp_names

          row_annotation <- data.frame(
            `Enriched in` = sp_enrichment[rownames(cor_matrix)],
            check.names = FALSE, stringsAsFactors = FALSE
          )
          rownames(row_annotation) <- rownames(cor_matrix)

          # Order species by enrichment group
          healthy_species <- rownames(row_annotation)[row_annotation$`Enriched in` == "Healthy"]
          cancer_species <- rownames(row_annotation)[row_annotation$`Enriched in` == "Cancer"]

          healthy_ordered <- order_within_group(healthy_species, cor_matrix)
          cancer_ordered <- order_within_group(cancer_species, cor_matrix)
          final_order <- c(healthy_ordered, cancer_ordered)

          cor_matrix_ordered <- cor_matrix[final_order, , drop = FALSE]
          sig_labels_ordered <- sig_labels[final_order, , drop = FALSE]
          row_annotation_ordered <- row_annotation[final_order, , drop = FALSE]

          # KO column labels
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
          gaps_row <- if (n_healthy_sp > 0 & n_healthy_sp < nrow(cor_matrix_ordered)) {
            n_healthy_sp
          } else NULL

          n_sp <- nrow(cor_matrix_ordered)
          n_ko <- ncol(cor_matrix_ordered)
          cell_w <- max(90, 130 - n_ko * 5)
          cell_h <- max(14, 20 - n_sp * 0.15)

          png("results_integrated/sarcosine/sarcosine_species_correlation_heatmap_pooled.png",
              width = max(10, n_ko * 1.5 + 5), height = max(7, n_sp * 0.35 + 4),
              units = "in", res = 200)
          pheatmap(cor_matrix_ordered,
                   color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                   cluster_rows = FALSE, cluster_cols = FALSE,
                   gaps_row = gaps_row,
                   annotation_row = row_annotation_ordered,
                   annotation_colors = ann_colors,
                   display_numbers = sig_labels_ordered,
                   number_color = "black", fontsize_number = 10,
                   fontsize_row = 7, fontsize_col = 8,
                   cellwidth = cell_w, cellheight = cell_h,
                   main = paste0("Integrated 4 CRC Cohorts: Species-Sarcosine KO Correlation (n=",
                                 length(common_samples), ")\n",
                                 "[* p<0.05, ** p<0.01, *** p<0.001]"),
                   angle_col = 0)
          dev.off()
          cat("Sarcosine species correlation heatmap saved.\n")
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
          labs(title = paste0("Integrated 4 CRC Cohorts: Species Associated with Sarcosine Metabolism"),
               subtitle = "Spearman correlation (FDR<0.05, |rho|>0.3)",
               x = "", y = "Mean Spearman rho") +
          theme_bw() +
          theme(plot.title = element_text(size = 10, face = "bold"),
                axis.text.y = element_text(size = 7), legend.position = "bottom")

        ggsave("results_integrated/sarcosine/sarcosine_associated_species_pooled.png", p_sp,
               width = 10, height = max(5, nrow(species_summary) * 0.3 + 2), dpi = 200)
        cat("Sarcosine associated species barplot saved.\n")
      }
    }
  }

  #############################################################################
  # PART 4b: SARCOSINE-ASSOCIATED SPECIES WITH PREVALENCE FILTER
  #############################################################################
  cat("\n\n========== PART 4b: SARCOSINE SPECIES CORRELATION (FILTERED) ==========\n")

  if (exists("cor_results") && nrow(cor_results) > 0 &&
      exists("sig_cors") && nrow(sig_cors) > 0) {

    pass_kos <- sarc_results$KO[sarc_results$Pass_Prevalence == TRUE &
                                 !is.na(sarc_results$p_value) &
                                 sarc_results$p_value < 0.05]
    sig_cors_filt <- sig_cors %>% filter(KO %in% pass_kos)

    cat("KOs passing prevalence + significance:", length(pass_kos), "\n")
    cat("Significant correlations after filter:", nrow(sig_cors_filt), "\n")

    if (nrow(sig_cors_filt) > 0) {
      filt_species_list <- unique(sig_cors_filt$Species_full)
      filt_ko_list <- unique(sig_cors_filt$KO)

      if (length(filt_species_list) > 1 & length(filt_ko_list) >= 1) {
        filt_sp_max_rho <- sig_cors_filt %>%
          group_by(Species_full) %>%
          summarise(max_rho = max(abs(rho)), .groups = "drop") %>%
          arrange(desc(max_rho)) %>% head(30)

        filt_top_sp_full <- filt_sp_max_rho$Species_full
        filt_top_sp_names <- sapply(filt_top_sp_full, get_species_name)

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

        filt_sig_labels <- matrix("", nrow = nrow(filt_cor_matrix), ncol = ncol(filt_cor_matrix))
        rownames(filt_sig_labels) <- rownames(filt_cor_matrix)
        colnames(filt_sig_labels) <- colnames(filt_cor_matrix)
        for (i in seq_len(nrow(filt_pval_matrix))) {
          for (j in seq_len(ncol(filt_pval_matrix))) {
            p <- filt_pval_matrix[i, j]
            if (!is.na(p)) {
              if (p < 0.001) filt_sig_labels[i, j] <- "***"
              else if (p < 0.01) filt_sig_labels[i, j] <- "**"
              else if (p < 0.05) filt_sig_labels[i, j] <- "*"
            }
          }
        }

        # Enrichment mapping
        filt_sp_enrichment <- sapply(filt_top_sp_full, function(sp) {
          sp_taxa <- sub("\\|t__.*$", "", sp)
          if (sp_taxa %in% names(enrichment_map)) return(enrichment_map[sp_taxa])
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
        filt_gaps_row <- if (filt_n_healthy_sp > 0 & filt_n_healthy_sp < nrow(filt_cor_ordered)) {
          filt_n_healthy_sp
        } else NULL

        # Paper-style rendering: larger fonts + horizontally-centered title.
        # (kept in sync with scripts/regenerate_correlation_heatmap_paper.R)
        ph_filt <- pheatmap(filt_cor_ordered,
                 color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                 cluster_rows = FALSE, cluster_cols = FALSE,
                 gaps_row = filt_gaps_row,
                 annotation_row = filt_row_ordered,
                 annotation_colors = filt_ann_colors,
                 display_numbers = filt_sig_ordered,
                 number_color = "black",
                 fontsize = 15, fontsize_number = 16,
                 fontsize_row = 14, fontsize_col = 14,
                 cellwidth = 205, cellheight = 24,
                 angle_col = 0, main = "", silent = TRUE)

        filt_title1 <- "Integrated 4 CRC Cohorts: Species-Sarcosine KO Correlation"
        filt_title2 <- paste0("(Prevalence >= ", PREV_THRESHOLD,
                              "% + Wilcoxon p<0.05)   [* p<0.05, ** p<0.01, *** p<0.001]")

        png("results_integrated/sarcosine/sarcosine_species_correlation_heatmap_filtered_pooled.png",
            width = 12.5, height = 12.4, units = "in", res = 300, type = "quartz", bg = "white")
        grid::grid.newpage()
        # outer margins so the bottom-left "Enriched in" label isn't clipped
        grid::pushViewport(grid::viewport(x = grid::unit(0.55, "in"),
            width = grid::unit(1, "npc") - grid::unit(0.65, "in"), just = "left"))
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(
            2, 1, heights = grid::unit.c(grid::unit(1.05, "in"), grid::unit(1, "null")))))
        # title band (centered over the full plot region)
        grid::pushViewport(grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
        grid::grid.text(filt_title1, x = 0.5, y = 0.62, just = c("center", "center"),
                        gp = grid::gpar(fontsize = 19, fontface = "bold"))
        grid::grid.text(filt_title2, x = 0.5, y = 0.24, just = c("center", "center"),
                        gp = grid::gpar(fontsize = 13, fontface = "bold"))
        grid::popViewport()
        # heatmap band
        grid::pushViewport(grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
        grid::grid.draw(ph_filt$gtable)
        grid::popViewport()
        grid::popViewport()
        grid::popViewport()
        dev.off()
        cat("Filtered species correlation heatmap saved.\n")
      }
    } else {
      cat("No significant correlations after prevalence filter.\n")
    }
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
      Cohort = meta_ko$Cohort,
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

    plot_data <- ratio_df %>%
      select(Group, Degradation_sum, Production_sum, Log2_Prod_Deg_Ratio) %>%
      pivot_longer(cols = c(Degradation_sum, Production_sum, Log2_Prod_Deg_Ratio),
                   names_to = "Metric", values_to = "Value") %>%
      mutate(Metric = factor(Metric,
                             levels = c("Degradation_sum", "Production_sum", "Log2_Prod_Deg_Ratio"),
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
      pval_text = c(format_pval(deg_test$p.value),
                    format_pval(prod_test$p.value),
                    format_pval(ratio_test$p.value)),
      stringsAsFactors = FALSE
    )
    ypos_pd <- plot_data %>%
      group_by(Metric) %>%
      summarise(ymax = max(Value, na.rm = TRUE), .groups = "drop")
    pval_labels_pd <- merge(pval_labels_pd, ypos_pd, by = "Metric")

    p_ratio <- ggplot(plot_data, aes(x = Group, y = Value, fill = Group)) +
      geom_boxplot(outlier.shape = 21, alpha = 0.7) +
      geom_jitter(width = 0.15, size = 0.3, alpha = 0.15) +
      facet_wrap(~Metric, scales = "free_y", nrow = 1) +
      scale_fill_manual(values = c("Healthy" = "#1B7837", "Cancer" = "#B2182B")) +
      geom_text(data = pval_labels_pd, aes(x = 1.5, y = ymax * 1.1, label = pval_text),
                inherit.aes = FALSE, size = 3.5, fontface = "italic") +
      labs(title = paste0("Integrated 4 CRC Cohorts: Sarcosine Metabolism Pathway Balance (n=",
                          nrow(ko_mat), ")"),
           subtitle = paste0("Deg KOs: ", paste(deg_kos, collapse=","),
                            "  |  Prod KOs: ", paste(prod_kos, collapse=",")),
           x = "", y = "Abundance / Ratio") +
      theme_bw() +
      theme(plot.title = element_text(size = 10, face = "bold"),
            plot.subtitle = element_text(size = 7),
            strip.text = element_text(size = 9),
            legend.position = "bottom")

    ggsave("results_integrated/sarcosine/sarcosine_prod_vs_deg_pooled.png", p_ratio,
           width = 12, height = 5.5, dpi = 200)
    cat("Production vs Degradation figure saved.\n")
  } else {
    cat("Cannot compute Prod/Deg ratio (need both pathway KOs available).\n")
  }
}

###############################################################################
# PART 6: SARCOSINE-ASSOCIATED BACTERIA DIFFERENTIAL ABUNDANCE
###############################################################################
cat("\n\n========== PART 6: SARCOSINE BACTERIA DIFF ABUNDANCE ==========\n")

sp_corr_file <- "results_integrated/sarcosine/sarcosine_species_correlation_pooled.csv"
if (file.exists(sp_corr_file)) {
  sp_corr <- read.csv(sp_corr_file, stringsAsFactors = FALSE)
  sig_sp_corr <- sp_corr[sp_corr$p_adj < 0.05 & abs(sp_corr$rho) > 0.3, ]

  if (nrow(sig_sp_corr) > 0) {
    sp_max_rho <- sig_sp_corr %>%
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

    # Re-read bacteria data for species-level aggregation
    bact_long_sp <- data.frame()
    for (cohort_id in names(COHORT_DIRS)) {
      cohort_path <- COHORT_DIRS[[cohort_id]]
      bact_file <- list.files(cohort_path, pattern = "^Bacteria_.*\\.txt$",
                              full.names = TRUE)[1]
      bact <- read.delim(bact_file, stringsAsFactors = FALSE)
      colnames(bact) <- trimws(colnames(bact))
      # Rename to consistent column names
      colnames(bact) <- c("Taxa", "Run_ID", "Abundance")

      cohort_runs <- meta_pooled$Run.ID[meta_pooled$Cohort == cohort_id]
      bact_sp <- bact %>%
        filter(grepl("\\|s__", Taxa)) %>%
        filter(Run_ID %in% cohort_runs) %>%
        mutate(Species_taxa = sub("\\|t__.*$", "", Taxa)) %>%
        group_by(Species_taxa, Run_ID) %>%
        summarise(Abundance = sum(Abundance), .groups = "drop")
      bact_long_sp <- rbind(bact_long_sp, bact_sp)
      rm(bact, bact_sp); gc(verbose = FALSE)
    }

    bact_target <- bact_long_sp %>%
      filter(Species_taxa %in% top_species_df$Species_full)

    if (nrow(bact_target) > 0) {
      runs_hc <- meta_pooled$Run.ID
      diff_results_p6 <- data.frame()

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

        healthy_runs <- meta_pooled$Run.ID[meta_pooled$Group == "Healthy"]
        cancer_runs  <- meta_pooled$Run.ID[meta_pooled$Group == "Cancer"]
        healthy_vals <- abund_vec[healthy_runs]
        cancer_vals  <- abund_vec[cancer_runs]

        mean_healthy <- mean(healthy_vals, na.rm = TRUE)
        mean_cancer  <- mean(cancer_vals, na.rm = TRUE)
        pseudo_fc <- 1e-6
        log2fc <- log2((mean_cancer + pseudo_fc) / (mean_healthy + pseudo_fc))

        wt_p <- tryCatch({
          wilcox.test(cancer_vals, healthy_vals)$p.value
        }, error = function(e) NA)

        diff_results_p6 <- rbind(diff_results_p6, data.frame(
          Species = sp_name, Species_full = sp_taxa, Role = sp_roles,
          Mean_Healthy = mean_healthy, Mean_Cancer = mean_cancer,
          log2FC = log2fc, p_value = wt_p,
          Direction = ifelse(log2fc > 0, "Cancer-enriched", "Healthy-enriched"),
          max_abs_rho = top_species_df$max_abs_rho[i],
          stringsAsFactors = FALSE))
      }

      diff_results_p6$p_adj <- p.adjust(diff_results_p6$p_value, method = "BH")
      write.csv(diff_results_p6,
                "results_integrated/sarcosine/sarcosine_bacteria_diff_abundance_pooled.csv",
                row.names = FALSE)

      cat("Significant species (p_adj<0.05):",
          sum(diff_results_p6$p_adj < 0.05, na.rm = TRUE), "\n")

      # Plot
      diff_results_p6$Species_display <- gsub("_", " ", diff_results_p6$Species)
      diff_results_p6 <- diff_results_p6 %>% arrange(log2FC)
      diff_results_p6$Species_display <- factor(diff_results_p6$Species_display,
                                                 levels = diff_results_p6$Species_display)

      diff_results_p6$p_text <- sapply(seq_len(nrow(diff_results_p6)), function(i) {
        p <- diff_results_p6$p_adj[i]
        if (is.na(p) || p >= 0.05) return("")
        format_pval(p)
      })
      diff_results_p6$Role_label <- ifelse(
        diff_results_p6$Role == "Degradation", "Deg",
        ifelse(diff_results_p6$Role == "Production", "Prod", "Deg & Prod"))

      x_max <- max(abs(diff_results_p6$log2FC), na.rm = TRUE) * 1.35
      n_healthy_total <- sum(meta_pooled$Group == "Healthy")
      n_cancer_total  <- sum(meta_pooled$Group == "Cancer")
      role_colors <- c("Deg" = "#2196F3", "Prod" = "#FF9800", "Deg & Prod" = "#9C27B0")

      p_diffab <- ggplot(diff_results_p6, aes(x = log2FC, y = Species_display)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
        geom_col(aes(fill = Direction), width = 0.7, alpha = 0.85) +
        geom_text(data = diff_results_p6[diff_results_p6$p_text != "", ],
          aes(x = ifelse(log2FC > 0, log2FC + x_max * 0.03, log2FC - x_max * 0.03),
              label = p_text),
          hjust = ifelse(diff_results_p6$log2FC[diff_results_p6$p_text != ""] > 0, 0, 1),
          size = 2.8, color = "gray20", fontface = "italic") +
        geom_point(aes(x = x_max * 0.92, color = Role_label), size = 3, shape = 15) +
        geom_text(aes(x = x_max * 0.97, label = Role_label, color = Role_label),
                  size = 2.5, hjust = 0, fontface = "bold") +
        scale_fill_manual(values = c("Cancer-enriched" = "#B2182B",
                                      "Healthy-enriched" = "#1B7837"), name = "Direction") +
        scale_color_manual(values = role_colors, name = "Sarcosine\nPathway") +
        scale_x_continuous(limits = c(-x_max, x_max * 1.15), expand = c(0, 0)) +
        labs(title = "Sarcosine-Associated Bacteria: Differential Abundance (Pooled)",
             subtitle = paste0("Integrated 4 CRC Cohorts | Healthy (n=", n_healthy_total,
                              ") vs Cancer (n=", n_cancer_total, ")"),
             x = expression(log[2]~"Fold Change (Cancer / Healthy)"), y = NULL,
             caption = "Top species by |rho| from sarcosine KO correlation\nBars: log2FC | p-values: BH-adjusted") +
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
      ggsave("results_integrated/sarcosine/sarcosine_bacteria_diff_abundance_pooled.png",
             p_diffab, width = 10, height = fig_h, dpi = 200, bg = "white")
      cat("Sarcosine bacteria diff. abundance figure saved.\n")
    }
    rm(bact_long_sp); gc(verbose = FALSE)
  } else {
    cat("No significant species-KO correlations for diff. abundance analysis.\n")
  }
} else {
  cat("Species correlation file not found. Run Part 4 first.\n")
}

###############################################################################
# SUPPLEMENTARY: PER-COHORT CONSISTENCY VERIFICATION
###############################################################################
cat("\n\n========== SUPPLEMENTARY: PER-COHORT CONSISTENCY ==========\n")

# --- S1. Cohort sample summary ---
cat("\n--- S1. Cohort sample summary ---\n")
cohort_summary <- meta_pooled %>%
  group_by(Cohort, Group) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = n, values_fill = 0)

# Add country info
country_info <- meta_pooled %>%
  group_by(Cohort) %>%
  summarise(Country = paste(unique(country), collapse = ", "), .groups = "drop")
cohort_summary <- merge(cohort_summary, country_info, by = "Cohort")
cohort_summary$Total <- cohort_summary$Healthy + cohort_summary$Cancer

write.csv(cohort_summary, "results_integrated/supplementary/cohort_sample_summary.csv",
          row.names = FALSE)
print(cohort_summary)

# --- S2. Per-cohort sarcosine KO consistency ---
cat("\n--- S2. Per-cohort sarcosine KO consistency ---\n")

if (exists("available_kos") && length(available_kos) > 0) {
  consistency <- data.frame()

  for (cohort_id in unique(meta_pooled$Cohort)) {
    cohort_runs <- meta_ko$Run.ID[meta_ko$Cohort == cohort_id]
    if (length(cohort_runs) < 5) next

    cohort_ko <- ko_mat[cohort_runs, , drop = FALSE]
    cohort_meta <- meta_ko[meta_ko$Cohort == cohort_id, ]
    ch_idx <- which(cohort_meta$Group == "Healthy")
    cc_idx <- which(cohort_meta$Group == "Cancer")

    if (length(ch_idx) < 3 | length(cc_idx) < 3) next

    for (ko_id in available_kos) {
      h_vals <- cohort_ko[ch_idx, ko_id]
      c_vals <- cohort_ko[cc_idx, ko_id]
      wt <- tryCatch(wilcox.test(h_vals, c_vals),
                     error = function(e) list(p.value = NA))
      mean_h <- mean(h_vals)
      mean_c <- mean(c_vals)
      direction <- ifelse(mean_c > mean_h, "Cancer UP", "Healthy UP")

      consistency <- rbind(consistency, data.frame(
        Cohort = cohort_id, KO = ko_id,
        Enzyme = sarcosine_kos$Enzyme[sarcosine_kos$KO == ko_id],
        Role = sarcosine_kos$Role[sarcosine_kos$KO == ko_id],
        Mean_Healthy = mean_h, Mean_Cancer = mean_c,
        p_value = wt$p.value, Direction = direction,
        stringsAsFactors = FALSE))
    }
  }

  # Add pooled results for comparison
  pooled_directions <- sarc_results %>%
    mutate(Pooled_Direction = ifelse(Mean_Cancer > Mean_Healthy, "Cancer UP", "Healthy UP"),
           Pooled_p = p_value) %>%
    select(KO, Pooled_Direction, Pooled_p)

  consistency <- merge(consistency, pooled_directions, by = "KO", all.x = TRUE)

  # Count direction consistency
  dir_consistency <- consistency %>%
    group_by(KO) %>%
    summarise(
      n_same_as_pooled = sum(Direction == Pooled_Direction[1]),
      n_total_cohorts = n(),
      .groups = "drop"
    )
  consistency <- merge(consistency, dir_consistency, by = "KO")

  write.csv(consistency,
            "results_integrated/supplementary/per_cohort_sarcosine_KO_consistency.csv",
            row.names = FALSE)
  cat("Per-cohort sarcosine KO consistency saved.\n")
  print(consistency[, c("Cohort", "KO", "Role", "Direction", "p_value",
                        "Pooled_Direction", "n_same_as_pooled")])
}

# --- S3. Per-cohort alpha diversity ---
cat("\n--- S3. Per-cohort alpha diversity ---\n")
alpha_cohort <- data.frame()

for (cohort_id in unique(meta_matched$Cohort)) {
  cohort_runs <- meta_matched$Run.ID[meta_matched$Cohort == cohort_id]
  if (length(cohort_runs) < 5) next

  cohort_bact <- bact_mat[cohort_runs, , drop = FALSE]
  cohort_meta <- meta_matched[meta_matched$Cohort == cohort_id, ]

  cohort_alpha <- data.frame(
    Run.ID = rownames(cohort_bact),
    Shannon = diversity(cohort_bact, index = "shannon"),
    Richness = specnumber(cohort_bact),
    stringsAsFactors = FALSE
  )
  cohort_alpha <- merge(cohort_alpha, cohort_meta[, c("Run.ID", "Group")], by = "Run.ID")

  shan_p <- tryCatch(wilcox.test(Shannon ~ Group, data = cohort_alpha)$p.value,
                     error = function(e) NA)
  rich_p <- tryCatch(wilcox.test(Richness ~ Group, data = cohort_alpha)$p.value,
                     error = function(e) NA)

  alpha_cohort <- rbind(alpha_cohort, data.frame(
    Cohort = cohort_id,
    n_Healthy = sum(cohort_meta$Group == "Healthy"),
    n_Cancer = sum(cohort_meta$Group == "Cancer"),
    Shannon_p = shan_p,
    Richness_p = rich_p,
    stringsAsFactors = FALSE
  ))
}

write.csv(alpha_cohort, "results_integrated/supplementary/per_cohort_alpha_diversity.csv",
          row.names = FALSE)
print(alpha_cohort)

# --- S4. Per-cohort PERMANOVA ---
cat("\n--- S4. Per-cohort PERMANOVA ---\n")
perm_cohort_results <- data.frame()

for (cohort_id in unique(meta_matched$Cohort)) {
  cohort_runs <- meta_matched$Run.ID[meta_matched$Cohort == cohort_id]
  if (length(cohort_runs) < 5) next

  cohort_bact <- bact_mat[cohort_runs, , drop = FALSE]
  cohort_meta <- meta_matched[meta_matched$Cohort == cohort_id, ]

  cohort_dist <- vegdist(cohort_bact, method = "bray")
  set.seed(42)
  perm_res <- tryCatch({
    adonis2(cohort_dist ~ Group, data = cohort_meta, permutations = 999)
  }, error = function(e) NULL)

  if (!is.null(perm_res)) {
    perm_cohort_results <- rbind(perm_cohort_results, data.frame(
      Cohort = cohort_id,
      R2 = round(perm_res$R2[1], 4),
      p_value = perm_res$`Pr(>F)`[1],
      stringsAsFactors = FALSE
    ))
  }
}

write.csv(perm_cohort_results,
          "results_integrated/supplementary/per_cohort_permanova.csv", row.names = FALSE)
print(perm_cohort_results)

# --- S5. Per-cohort Prod/Deg ratio ---
cat("\n--- S5. Per-cohort Prod/Deg ratio ---\n")

if (exists("deg_kos") && exists("prod_kos") &&
    length(deg_kos) > 0 && length(prod_kos) > 0) {

  proddeg_cohort <- data.frame()

  for (cohort_id in unique(meta_ko$Cohort)) {
    cohort_runs <- meta_ko$Run.ID[meta_ko$Cohort == cohort_id]
    if (length(cohort_runs) < 5) next

    cohort_ko_sub <- ko_mat[cohort_runs, , drop = FALSE]
    cohort_meta_sub <- meta_ko[meta_ko$Cohort == cohort_id, ]

    ch_idx <- which(cohort_meta_sub$Group == "Healthy")
    cc_idx <- which(cohort_meta_sub$Group == "Cancer")
    if (length(ch_idx) < 3 | length(cc_idx) < 3) next

    d_sum <- rowSums(cohort_ko_sub[, deg_kos, drop = FALSE])
    p_sum <- rowSums(cohort_ko_sub[, prod_kos, drop = FALSE])
    l2r <- log2((p_sum + 1e-8) / (d_sum + 1e-8))

    ratio_test_c <- tryCatch(
      wilcox.test(l2r[ch_idx], l2r[cc_idx]),
      error = function(e) list(p.value = NA))

    mean_h <- mean(l2r[ch_idx])
    mean_c <- mean(l2r[cc_idx])
    direction <- ifelse(mean_c > mean_h, "Cancer higher Prod/Deg",
                        "Healthy higher Prod/Deg")

    proddeg_cohort <- rbind(proddeg_cohort, data.frame(
      Cohort = cohort_id,
      n_Healthy = length(ch_idx),
      n_Cancer = length(cc_idx),
      mean_ratio_Healthy = round(mean_h, 4),
      mean_ratio_Cancer = round(mean_c, 4),
      p_value = ratio_test_c$p.value,
      Direction = direction,
      stringsAsFactors = FALSE
    ))
  }

  write.csv(proddeg_cohort,
            "results_integrated/supplementary/per_cohort_prod_deg_ratio.csv",
            row.names = FALSE)
  print(proddeg_cohort)
}

cat("\n\n=== INTEGRATED POOLED ANALYSIS COMPLETE ===\n")
cat("Results saved in: results_integrated/\n")
cat("Total samples analyzed:", nrow(meta_pooled), "\n")
cat("Cohorts:", paste(unique(meta_pooled$Cohort), collapse = ", "), "\n")
