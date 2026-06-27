###############################################################################
# Regenerate the FILTERED Species-Sarcosine KO correlation heatmap
# for paper use: larger fonts + centered title.
#
# Faithfully reproduces Part 4 / 4b of integrated_analysis_pooled.R, but instead
# of re-loading the giant KO JSON files it reads the already-computed
# correlation table and re-reads only the (small) Bacteria_*.txt files to
# recover the Healthy/Cancer enrichment annotation.
#
# Output: results_integrated/sarcosine/
#           sarcosine_species_correlation_heatmap_filtered_pooled.png
# (paper-style: larger fonts + centered title; same content as Part 4b)
###############################################################################

options(bitmapType = "quartz")

library(dplyr)
library(tidyr)
library(pheatmap)
library(grid)

# ----------------------------------------------------------------------------
# CONFIG — tweak these to taste
# ----------------------------------------------------------------------------
BASE_DIR <- "Healthy_vs_Cancer_4_CRC_cohorts_integrated"   # relative to the analysis-folder root (working dir)

FS_BASE   <- 20   # base font (legend numbers, annotation legend)
FS_ROW    <- 18   # species names (rows)
FS_COL    <- 18   # KO labels (columns, bottom)
FS_NUMBER <- 21   # the *** significance marks inside cells
FS_TITLE1 <- 34   # title line 1 (enlarged for readability)
FS_TITLE2 <- 21   # title line 2 (thresholds)

CELL_W <- 205     # cell width  (points) — wide enough for the EC labels below
CELL_H <- 24      # cell height (points)

TITLE_X    <- 0.5   # horizontal position of title, 0=left .. 1=right (device npc)
TITLE_H_IN <- 1.25  # vertical band reserved for the title (inches; widened for the larger title font)
MARGIN_L   <- 1.9  # left margin (inches) so "Enriched in" label isn't clipped
MARGIN_R   <- 1.1   # right margin (inches) so the "Enriched in" legend isn't clipped

DEV_W <- 14.9     # device width  (inches)
DEV_H <- 13.6     # device height (inches; raised to give the 3-line column labels bottom room)
DEV_RES <- 300
BOTTOM_PAD_IN <- 0.7  # reserved band at the device bottom so the 3-line col labels ("…/Degradation"/"…/Production") are not clipped

PREV_THRESHOLD <- 10

COHORT_DIRS <- list(
  PRJEB6070   = file.path(BASE_DIR, "PRJEB6070_CRC_AdenomatousPolyps"),
  PRJNA429097 = file.path(BASE_DIR, "PRJNA429097_CRC"),
  PRJEB10878  = file.path(BASE_DIR, "PRJEB10878_CRC"),
  PRJEB27928  = file.path(BASE_DIR, "PRJEB27928_CRC")
)

# Sarcosine KO panel (verbatim from integrated_analysis_pooled.R)
sarcosine_kos <- data.frame(
  KO = c("K00301","K00302","K00303","K00304","K00305","K00306",
         "K00315","K00552","K08688"),
  Enzyme = c("Sarcosine oxidase (monomeric)",
             "Sarcosine oxidase subunit alpha (soxA)",
             "Sarcosine oxidase subunit beta (soxB)",
             "Sarcosine oxidase subunit delta (soxD)",
             "Sarcosine oxidase subunit gamma (soxG)",
             "Sarcosine oxidase / L-pipecolate oxidase (PIPOX)",
             "Dimethylglycine dehydrogenase (DMGDH)",
             "Glycine N-methyltransferase (GNMT)",
             "Creatinase"),
  EC = c("EC:1.5.3.1","EC:1.5.3.1/1.5.3.24","EC:1.5.3.1/1.5.3.24",
         "EC:1.5.3.1/1.5.3.24","EC:1.5.3.1/1.5.3.24","EC:1.5.3.1/1.5.3.7",
         "EC:1.5.8.4","EC:2.1.1.20","EC:3.5.3.3"),
  Role = c("Degradation","Degradation","Degradation","Degradation",
           "Degradation","Degradation","Production","Production","Production"),
  Short = c("SOX (mono)","soxA","soxB","soxD","soxG","PIPOX",
            "DMGDH","GNMT","Creatinase"),
  stringsAsFactors = FALSE
)

get_species_name <- function(taxa_str) {
  parts <- strsplit(taxa_str, "\\|")[[1]]
  sp <- parts[grep("^s__", parts)]
  if (length(sp) > 0) return(gsub("^s__", "", sp[1]))
  return(taxa_str)
}

setwd(BASE_DIR)

# ----------------------------------------------------------------------------
# 1. Metadata (pooled) — replicate Section 1.1
# ----------------------------------------------------------------------------
cat("--- Loading metadata ---\n")
meta_pooled <- data.frame()
for (cohort_id in names(COHORT_DIRS)) {
  cohort_path <- COHORT_DIRS[[cohort_id]]
  meta_file <- list.files(cohort_path, pattern = "^selected_project_.*\\.txt$",
                          full.names = TRUE)
  meta_file <- sort(meta_file, decreasing = TRUE)[1]
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
  cancer_label <- cancer_pheno[1]
  meta_hc <- meta %>%
    filter(Phenotype.name %in% c("Health", cancer_label)) %>%
    mutate(Group = ifelse(Phenotype.name == "Health", "Healthy", "Cancer"))
  if ("Assay.type" %in% colnames(meta_hc)) {
    meta_hc <- meta_hc %>% filter(Assay.type == "WGS")
  }
  meta_hc$Cohort <- cohort_id
  meta_pooled <- rbind(meta_pooled, meta_hc)
}
cat("  Pooled samples:", nrow(meta_pooled),
    "| Healthy:", sum(meta_pooled$Group == "Healthy"),
    "| Cancer:", sum(meta_pooled$Group == "Cancer"), "\n")

# ----------------------------------------------------------------------------
# 2. Correlation table + filtered KO set
# ----------------------------------------------------------------------------
cat("--- Loading correlation + KO tables ---\n")
cor_results <- read.csv("results_integrated/sarcosine/sarcosine_species_correlation_pooled.csv",
                        stringsAsFactors = FALSE)
sarc_filtered_csv <- read.csv("results_integrated/sarcosine/sarcosine_KO_comparison_filtered_pooled.csv",
                              stringsAsFactors = FALSE)

# Same selection as Part 4b:
pass_kos <- sarc_filtered_csv$KO[sarc_filtered_csv$Pass_Prevalence == TRUE &
                                   !is.na(sarc_filtered_csv$p_value) &
                                   sarc_filtered_csv$p_value < 0.05]
cat("  KOs passing prevalence>=", PREV_THRESHOLD, "% & Wilcoxon p<0.05:",
    paste(pass_kos, collapse = ", "), "\n")

sig_cors <- cor_results %>% filter(p_adj < 0.05 & abs(rho) > 0.3)
sig_cors_filt <- sig_cors %>% filter(KO %in% pass_kos)
cat("  Significant correlations after filter:", nrow(sig_cors_filt), "\n")

filt_ko_list <- unique(sig_cors_filt$KO)
filt_sp_max_rho <- sig_cors_filt %>%
  group_by(Species_full) %>%
  summarise(max_rho = max(abs(rho)), .groups = "drop") %>%
  arrange(desc(max_rho)) %>% head(30)
filt_top_sp_full  <- filt_sp_max_rho$Species_full
filt_top_sp_names <- sapply(filt_top_sp_full, get_species_name)
cat("  Species shown:", length(filt_top_sp_full), "\n")

# ----------------------------------------------------------------------------
# 3. Healthy/Cancer enrichment annotation — replicate Part 4 (re-read Bacteria)
# ----------------------------------------------------------------------------
cat("--- Computing Healthy/Cancer enrichment (re-reading Bacteria files) ---\n")
bact_long_species <- data.frame()
for (cohort_id in names(COHORT_DIRS)) {
  cohort_path <- COHORT_DIRS[[cohort_id]]
  bact_file <- list.files(cohort_path, pattern = "^Bacteria_.*\\.txt$", full.names = TRUE)[1]
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
  Healthy = rowMeans(as.matrix(bact_sp_wide[, h_cols, drop = FALSE])),
  Cancer  = rowMeans(as.matrix(bact_sp_wide[, c_cols, drop = FALSE])),
  stringsAsFactors = FALSE
) %>% mutate(Enriched = ifelse(Cancer >= Healthy, "Cancer", "Healthy"))
enrichment_map <- setNames(sp_enrichment_df$Enriched, sp_enrichment_df$Species_taxa)
rm(bact_sp_wide, bact_long_species); gc(verbose = FALSE)

# ----------------------------------------------------------------------------
# 4. Build matrices + ordering (replicate Part 4b)
# ----------------------------------------------------------------------------
order_within_group <- function(species_names, mat) {
  if (length(species_names) <= 1) return(species_names)
  sub_mat <- mat[species_names, , drop = FALSE]
  d <- dist(sub_mat, method = "euclidean")
  hc <- hclust(d, method = "complete")
  return(species_names[hc$order])
}

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

filt_sp_enrichment <- sapply(filt_top_sp_full, function(sp) {
  sp_taxa <- sub("\\|t__.*$", "", sp)
  if (sp_taxa %in% names(enrichment_map)) return(enrichment_map[sp_taxa])
  if (sp %in% names(enrichment_map)) return(enrichment_map[sp])
  return("Unknown")
})
names(filt_sp_enrichment) <- filt_top_sp_names

filt_row_ann <- data.frame(`Enriched in` = filt_sp_enrichment[rownames(filt_cor_matrix)],
                           check.names = FALSE, stringsAsFactors = FALSE)
rownames(filt_row_ann) <- rownames(filt_cor_matrix)

filt_healthy_sp <- rownames(filt_row_ann)[filt_row_ann$`Enriched in` == "Healthy"]
filt_cancer_sp  <- rownames(filt_row_ann)[filt_row_ann$`Enriched in` == "Cancer"]
filt_healthy_ordered <- order_within_group(filt_healthy_sp, filt_cor_matrix)
filt_cancer_ordered  <- order_within_group(filt_cancer_sp, filt_cor_matrix)
filt_final_order <- c(filt_healthy_ordered, filt_cancer_ordered)

filt_cor_ordered <- filt_cor_matrix[filt_final_order, , drop = FALSE]
filt_sig_ordered <- filt_sig_labels[filt_final_order, , drop = FALSE]
filt_row_ordered <- filt_row_ann[filt_final_order, , drop = FALSE]

filt_ko_labels <- sapply(colnames(filt_cor_ordered), function(k) {
  info <- sarcosine_kos %>% filter(KO == k)
  if (nrow(info) > 0) paste0(k, "\n", info$Short[1], " (", info$EC[1], ")\n", info$Role[1]) else k
})
colnames(filt_cor_ordered) <- filt_ko_labels
colnames(filt_sig_ordered) <- filt_ko_labels

ann_colors <- list(`Enriched in` = c("Healthy" = "#1B7837", "Cancer" = "#B2182B"))
n_healthy_sp <- length(filt_healthy_ordered)
gaps_row <- if (n_healthy_sp > 0 & n_healthy_sp < nrow(filt_cor_ordered)) n_healthy_sp else NULL

cat("  Healthy-enriched species:", n_healthy_sp,
    "| Cancer-enriched species:", length(filt_cancer_ordered), "\n")

# ----------------------------------------------------------------------------
# 5. Render: pheatmap grob (no built-in title) + manually centered title
# ----------------------------------------------------------------------------
ph <- pheatmap(filt_cor_ordered,
               color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
               cluster_rows = FALSE, cluster_cols = FALSE,
               gaps_row = gaps_row,
               annotation_row = filt_row_ordered,
               annotation_colors = ann_colors,
               display_numbers = filt_sig_ordered,
               number_color = "black",
               fontsize = FS_BASE,
               fontsize_number = FS_NUMBER,
               fontsize_row = FS_ROW,
               fontsize_col = FS_COL,
               cellwidth = CELL_W, cellheight = CELL_H,
               angle_col = 0,
               main = "",
               silent = TRUE)

title1 <- "Species-Sarcosine KO Correlation (4 CRC cohorts)"
title2 <- paste0("(Prevalence ≥ ", PREV_THRESHOLD,
                 "% + Wilcoxon p<0.05)   [* p<0.05, ** p<0.01, *** p<0.001]")

out_file <- "results_integrated/sarcosine/sarcosine_species_correlation_heatmap_filtered_pooled.png"
png(out_file, width = DEV_W, height = DEV_H, units = "in", res = DEV_RES,
    type = "quartz", bg = "white")
grid.newpage()
# Outer viewport with left/right margins so the row-annotation name
# ("Enriched in") at the bottom-left is not clipped by the device edge.
pushViewport(viewport(x = unit(MARGIN_L, "in"),
                      width = unit(1, "npc") - unit(MARGIN_L + MARGIN_R, "in"),
                      just = "left"))
pushViewport(viewport(layout = grid.layout(
  3, 1, heights = unit.c(unit(TITLE_H_IN, "in"), unit(1, "null"), unit(BOTTOM_PAD_IN, "in")))))

# Title band
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.text(title1, x = TITLE_X, y = 0.62, just = c("center", "center"),
          gp = gpar(fontsize = FS_TITLE1, fontface = "bold"))
grid.text(title2, x = TITLE_X, y = 0.24, just = c("center", "center"),
          gp = gpar(fontsize = FS_TITLE2, fontface = "bold"))
popViewport()

# Heatmap band
pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1))
grid.draw(ph$gtable)
popViewport()
popViewport()   # layout
popViewport()   # outer margin
dev.off()

cat("\nSaved:", out_file, "\n")
