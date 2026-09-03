#!/usr/bin/env Rscript
# ============================================================================
# PER-COHORT paper-style Species-Sarcosine KO correlation heatmaps (FILTERED).
# Per-cohort analogue of scripts/regenerate_correlation_heatmap_paper.R (pooled).
# For each cohort, re-renders the Part-4b filtered heatmap (species x the
# sarcosine KOs that pass prevalence >=10% AND Wilcoxon p<0.05) in paper style
# (concise title + larger fonts), reading the already-computed per-cohort inputs.
#
# READ-ONLY: reads <cohort>/results_sarcosine/sarcosine_species_correlation.csv,
#   sarcosine_KO_comparison_filtered.csv, and the cohort Bacteria_*.txt (for the
#   Healthy/Cancer enrichment annotation only). NO statistic is recomputed.
# Selection (verbatim from Part 4b): species-KO Spearman p_adj<0.05 & |rho|>0.3,
#   restricted to prevalence+Wilcoxon-significant KOs; top 30 species by max|rho|.
#
# Output: <cohort>/results_sarcosine/sarcosine_species_correlation_heatmap_filtered_paper.png
# Run: Rscript regenerate_correlation_heatmap_paper_percohort.R   (paths auto-discovered)
# ============================================================================
options(bitmapType = "quartz")
suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(pheatmap); library(grid) })

script_dir <- "Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts"   # relative to the analysis-folder root (working dir)
crc_root   <- "."   # scripts/ -> integrated/ -> CRC root (= the analysis-folder root)

COHORTS <- c(PRJEB6070   = "PRJEB6070_CRC_AdenomatousPolyps",
             PRJEB10878  = "PRJEB10878_CRC",
             PRJEB27928  = "PRJEB27928_CRC",
             PRJNA429097 = "PRJNA429097_CRC")

PREV_THRESHOLD <- 10
sarcosine_kos <- data.frame(
  KO    = c("K00301","K00302","K00303","K00304","K00305","K00306","K00315","K00552","K08688"),
  EC    = c("EC:1.5.3.1","EC:1.5.3.1/1.5.3.24","EC:1.5.3.1/1.5.3.24","EC:1.5.3.1/1.5.3.24",
            "EC:1.5.3.1/1.5.3.24","EC:1.5.3.1/1.5.3.7","EC:1.5.8.4","EC:2.1.1.20","EC:3.5.3.3"),
  Role  = c("Degradation","Degradation","Degradation","Degradation","Degradation","Degradation",
            "Production","Production","Production"),
  Short = c("SOX (mono)","soxA","soxB","soxD","soxG","PIPOX","DMGDH","GNMT","Creatinase"),
  stringsAsFactors = FALSE)

# paper style
FS_BASE <- 20; FS_ROW <- 16; FS_COL <- 18; FS_NUMBER <- 21; FS_TITLE1 <- 31; FS_TITLE2 <- 17
CELL_W <- 265; CELL_H <- 22; TITLE_H_IN <- 1.2; MARGIN_L <- 1.9; MARGIN_R <- 1.1; DEV_RES <- 300  # margins match the verified pooled script so the "Enriched in" row-title (left) and legend-title (right) are not clipped
BOTTOM_PAD_IN <- 0.7  # reserved band at the device bottom so the 3-line col labels (".../Degradation" / ".../Production") are not clipped (mirrors the verified pooled fix)

get_species_name <- function(taxa_str) {
  parts <- strsplit(taxa_str, "\\|")[[1]]; sp <- parts[grep("^s__", parts)]
  if (length(sp) > 0) gsub("^s__", "", sp[1]) else taxa_str
}
order_within_group <- function(species_names, mat) {
  if (length(species_names) <= 1) return(species_names)
  hc <- hclust(dist(mat[species_names, , drop = FALSE], method = "euclidean"), method = "complete")
  species_names[hc$order]
}
read_meta_hc <- function(cohort_dir) {
  mf <- sort(list.files(cohort_dir, "^selected_project_.*\\.txt$", full.names = TRUE), decreasing = TRUE)[1]
  ml <- readLines(mf); hdr <- strsplit(ml[2], "\t")[[1]]
  dl <- ml[3:length(ml)]; dl <- dl[dl != ""]
  rows <- lapply(strsplit(dl, "\t"), function(x) if (length(x) >= length(hdr)) x[1:length(hdr)] else c(x, rep(NA, length(hdr) - length(x))))
  m <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE); colnames(m) <- gsub(" ", ".", hdr)
  cancer <- unique(m$Phenotype.name)[!unique(m$Phenotype.name) %in% c("Health", "Adenomatous Polyps")][1]
  m <- m[m$Phenotype.name %in% c("Health", cancer), ]
  if ("Assay.type" %in% colnames(m)) m <- m[m$Assay.type == "WGS", ]
  data.frame(Run.ID = m$Run.ID, Group = ifelse(m$Phenotype.name == "Health", "Healthy", "Cancer"), stringsAsFactors = FALSE)
}

render_cohort <- function(cn, folder) {
  cohort_dir <- file.path(crc_root, folder)
  res <- file.path(cohort_dir, "results_sarcosine")
  cor_results <- read.csv(file.path(res, "sarcosine_species_correlation.csv"), stringsAsFactors = FALSE)
  filt <- read.csv(file.path(res, "sarcosine_KO_comparison_filtered.csv"), stringsAsFactors = FALSE)
  pass_kos <- filt$KO[filt$Pass_Prevalence == TRUE & !is.na(filt$p_value) & filt$p_value < 0.05]
  sig <- cor_results %>% filter(p_adj < 0.05 & abs(rho) > 0.3 & KO %in% pass_kos)
  if (nrow(sig) == 0 || length(unique(sig$Species_full)) < 2) { cat(cn, ": no filtered sig corr; skipped\n"); return(invisible()) }
  ko_list <- unique(sig$KO)
  top <- sig %>% group_by(Species_full) %>% summarise(max_rho = max(abs(rho)), .groups = "drop") %>%
    arrange(desc(max_rho)) %>% head(30)
  sp_full <- top$Species_full; sp_names <- sapply(sp_full, get_species_name)

  # --- Healthy/Cancer enrichment from the cohort Bacteria file ---
  meta_hc <- read_meta_hc(cohort_dir)
  bf <- list.files(cohort_dir, "^Bacteria_.*\\.txt$", full.names = TRUE)[1]
  bact <- read.delim(bf, stringsAsFactors = FALSE); colnames(bact) <- trimws(colnames(bact))
  bsp <- bact %>% filter(grepl("\\|s__", Taxa) & Run.ID %in% meta_hc$Run.ID) %>%
    mutate(Species_taxa = sub("\\|t__.*$", "", Taxa)) %>%
    group_by(Species_taxa, Run.ID) %>% summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop")
  bw <- bsp %>% pivot_wider(names_from = Run.ID, values_from = Abundance, values_fill = 0)
  hcol <- intersect(meta_hc$Run.ID[meta_hc$Group == "Healthy"], colnames(bw))
  ccol <- intersect(meta_hc$Run.ID[meta_hc$Group == "Cancer"], colnames(bw))
  enr <- data.frame(Species_taxa = bw$Species_taxa,
                    H = rowMeans(as.matrix(bw[, hcol, drop = FALSE])),
                    C = rowMeans(as.matrix(bw[, ccol, drop = FALSE])), stringsAsFactors = FALSE)
  emap <- setNames(ifelse(enr$C >= enr$H, "Cancer", "Healthy"), enr$Species_taxa)

  # --- matrices ---
  M <- matrix(0, length(sp_full), length(ko_list), dimnames = list(sp_names, ko_list))
  P <- matrix(NA_real_, length(sp_full), length(ko_list), dimnames = list(sp_names, ko_list))
  for (r in seq_len(nrow(sig))) {
    sn <- get_species_name(sig$Species_full[r]); kn <- sig$KO[r]
    if (sn %in% rownames(M) && kn %in% colnames(M)) { M[sn, kn] <- sig$rho[r]; P[sn, kn] <- sig$p_adj[r] }
  }
  L <- matrix("", nrow(M), ncol(M), dimnames = dimnames(M))
  L[!is.na(P) & P < 0.05]  <- "*"; L[!is.na(P) & P < 0.01] <- "**"; L[!is.na(P) & P < 0.001] <- "***"

  spenr <- sapply(sp_full, function(s) { st <- sub("\\|t__.*$", "", s)
    if (st %in% names(emap)) emap[st] else if (s %in% names(emap)) emap[s] else "Unknown" })
  names(spenr) <- sp_names
  rowann <- data.frame(`Enriched in` = spenr[rownames(M)], check.names = FALSE, stringsAsFactors = FALSE)
  rownames(rowann) <- rownames(M)
  hsp <- rownames(rowann)[rowann$`Enriched in` == "Healthy"]; csp <- rownames(rowann)[rowann$`Enriched in` == "Cancer"]
  ord <- c(order_within_group(hsp, M), order_within_group(csp, M))
  Mo <- M[ord, , drop = FALSE]; Lo <- L[ord, , drop = FALSE]; rowanno <- rowann[ord, , drop = FALSE]

  ko_lab <- sapply(colnames(Mo), function(k) { i <- sarcosine_kos[sarcosine_kos$KO == k, ]
    if (nrow(i)) paste0(k, "\n", i$Short[1], " (", i$EC[1], ")\n", i$Role[1]) else k })
  colnames(Mo) <- ko_lab; colnames(Lo) <- ko_lab
  anncol <- list(`Enriched in` = c("Healthy" = "#1B7837", "Cancer" = "#B2182B"))
  nH <- length(order_within_group(hsp, M)); gaps <- if (nH > 0 && nH < nrow(Mo)) nH else NULL
  nsp <- nrow(Mo); nko <- ncol(Mo)

  ph <- pheatmap(Mo, color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
                 cluster_rows = FALSE, cluster_cols = FALSE, gaps_row = gaps,
                 annotation_row = rowanno, annotation_colors = anncol,
                 display_numbers = Lo, number_color = "black",
                 fontsize = FS_BASE, fontsize_number = FS_NUMBER, fontsize_row = FS_ROW, fontsize_col = FS_COL,
                 cellwidth = CELL_W, cellheight = CELL_H, angle_col = 0, main = "", silent = TRUE)

  title1 <- paste0("Species-Sarcosine KO Correlation (", cn, ")")
  title2 <- paste0("(Prevalence >= ", PREV_THRESHOLD, "% + Wilcoxon p<0.05)   [* p<0.05, ** p<0.01, *** p<0.001]")
  DEV_W <- nko * (CELL_W / 72) + 9.85   # +2.35 vs the old +7.5 to absorb the wider L/R margins (1.9+1.1 vs 0.55+0.10), preserving the inner heatmap/species/legend width
  DEV_H <- nsp * (CELL_H / 72) + 3.2 + BOTTOM_PAD_IN
  out <- file.path(res, "sarcosine_species_correlation_heatmap_filtered_paper.png")
  png(out, width = DEV_W, height = DEV_H, units = "in", res = DEV_RES, type = "quartz", bg = "white")
  grid.newpage()
  pushViewport(viewport(x = unit(MARGIN_L, "in"), width = unit(1, "npc") - unit(MARGIN_L + MARGIN_R, "in"), just = "left"))
  pushViewport(viewport(layout = grid.layout(3, 1, heights = unit.c(unit(TITLE_H_IN, "in"), unit(1, "null"), unit(BOTTOM_PAD_IN, "in")))))
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
  grid.text(title1, x = 0.5, y = 0.62, just = c("center", "center"), gp = gpar(fontsize = FS_TITLE1, fontface = "bold"))
  grid.text(title2, x = 0.5, y = 0.24, just = c("center", "center"), gp = gpar(fontsize = FS_TITLE2, fontface = "bold"))
  popViewport()
  pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1)); grid.draw(ph$gtable); popViewport()
  popViewport(); popViewport(); dev.off()
  cat(sprintf("%-12s: %d species x %d KO (%s) -> %s (%.1f x %.1f in)\n",
      cn, nsp, nko, paste(ko_list, collapse = ","), basename(out), DEV_W, DEV_H))
}

for (cn in names(COHORTS)) render_cohort(cn, COHORTS[cn])
writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo_correlation_heatmap_paper_percohort.txt"))
cat("DONE.\n")
