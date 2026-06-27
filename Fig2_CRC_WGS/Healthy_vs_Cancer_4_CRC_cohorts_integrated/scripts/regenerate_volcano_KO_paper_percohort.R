#!/usr/bin/env Rscript
# ============================================================================
# PER-COHORT paper-style KEGG-KO differential-abundance volcano plots.
# Re-renders each of the 4 CRC cohorts' KO differential abundance
# (<cohort>/results_kegg/diff_KO_abundance.csv) as a paper-style volcano with
# the SIGNIFICANT sarcosine-metabolism KOs highlighted -- the per-cohort analogue
# of the pooled figure (scripts/regenerate_volcano_KO_paper.R).
#
# READ-ONLY: reads the already-computed per-cohort DA CSVs and RE-DRAWS.
#   NO statistic is recomputed (CSVs from analysis_healthy_vs_cancer.R Part 2:
#   >=10% prevalence, Wilcoxon rank-sum, BH). Significance (verbatim, same as the
#   pooled volcano and the point colouring): p_adj < 0.05 & |log2FC| > 0.5.
#
# LABELLING (faithful): labels EVERY sarcosine-panel KO that is SIGNIFICANT in
#   that cohort by the above rule (degradation: K00301-K00306; production:
#   K00315, K00552, K08688) -- so each labelled enzyme is a coloured (green/red)
#   significant point. A cohort with no significant sarcosine KO is annotated
#   as such. This reflects the actually-significant sarcosine genes per cohort
#   (e.g., PRJEB27928's significant sarcosine enzyme is soxA, not soxB).
#
# Output (one per cohort): <cohort>/results_kegg/volcano_KO_paper.png
# Run: Rscript regenerate_volcano_KO_paper_percohort.R   (paths auto-discovered)
# ============================================================================
options(bitmapType = "quartz")   # cairo unavailable; quartz writes non-ASCII paths (ragg cannot)
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(ggrepel) })

script_dir <- "Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts"   # relative to the analysis-folder root (working dir)
crc_root   <- "."   # scripts/ -> integrated/ -> CRC root (cohort folders live here) (= the analysis-folder root)

COHORTS <- c(PRJEB6070   = "PRJEB6070_CRC_AdenomatousPolyps",
             PRJEB10878  = "PRJEB10878_CRC",
             PRJEB27928  = "PRJEB27928_CRC",
             PRJNA429097 = "PRJNA429097_CRC")

# Full KEGG sarcosine-metabolism KO panel (verbatim roles from the project scripts).
SARC <- data.frame(
  KO    = c("K00301","K00302","K00303","K00304","K00305","K00306","K00315","K00552","K08688"),
  Short = c("SOX(mono)","soxA","soxB","soxD","soxG","PIPOX","DMGDH","GNMT","Creatinase"),
  Role  = c("Degradation","Degradation","Degradation","Degradation","Degradation","Degradation",
            "Production","Production","Production"),
  stringsAsFactors = FALSE)

FS_TITLE <- 26; FS_SUB <- 18; FS_AXIS_T <- 22; FS_AXIS_X <- 18; FS_LEG_T <- 20; FS_LEG_X <- 18
LABEL_SIZE <- 6.5; POINT_SIZE <- 1.6
COLS <- c("Enriched in Cancer" = "#B2182B", "Enriched in Healthy" = "#1B7837", "Not Significant" = "grey70")

save_png <- function(p, out, w, h) {   # quartz device handles both the non-ASCII path and unicode glyphs
  grDevices::png(out, width = w, height = h, units = "in", res = 300, type = "quartz", bg = "white")
  print(p); grDevices::dev.off()
}

for (cn in names(COHORTS)) {
  f <- file.path(crc_root, COHORTS[cn], "results_kegg", "diff_KO_abundance.csv")
  stopifnot(file.exists(f))
  d <- read.csv(f, stringsAsFactors = FALSE)
  stopifnot(all(c("KO", "log2FC", "p_adj") %in% names(d)))

  d$Significance <- "Not Significant"
  d$Significance[d$p_adj < 0.05 & d$log2FC >  0.5] <- "Enriched in Cancer"
  d$Significance[d$p_adj < 0.05 & d$log2FC < -0.5] <- "Enriched in Healthy"
  d$Significance <- factor(d$Significance, levels = names(COLS))

  # label only the SIGNIFICANT sarcosine-panel KOs in this cohort (label = coloured point)
  lab <- d %>% inner_join(SARC, by = "KO") %>%
    filter(p_adj < 0.05 & abs(log2FC) > 0.5) %>%
    mutate(txt = paste0(KO, " — ", Short, "\n(", Role, ")"),
           nx  = ifelse(log2FC > 0, 1.0, -1.0))

  ncan <- sum(d$Significance == "Enriched in Cancer")
  nhea <- sum(d$Significance == "Enriched in Healthy")

  p <- ggplot(d, aes(log2FC, -log10(p_adj))) +
    geom_point(aes(color = Significance), alpha = 0.6, size = POINT_SIZE) +
    scale_color_manual(values = COLS, drop = FALSE) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(title = paste0(cn, ": Differential KO Abundance"),
         subtitle = "Healthy vs CRC — significant sarcosine KOs highlighted (q < 0.05 & |log2FC| > 0.5)",
         x = "log2FC (Cancer/Healthy)", y = "-log10(adj. p-value)", color = "Significance") +
    theme_bw() +
    theme(plot.title = element_text(size = FS_TITLE, face = "bold"),
          plot.subtitle = element_text(size = FS_SUB),
          axis.title = element_text(size = FS_AXIS_T), axis.text = element_text(size = FS_AXIS_X),
          legend.title = element_text(size = FS_LEG_T), legend.text = element_text(size = FS_LEG_X),
          legend.position = "bottom") +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1)))

  if (nrow(lab) > 0) {
    p <- p + ggrepel::geom_text_repel(data = lab, aes(label = txt),
        color = "black", fontface = "bold", size = LABEL_SIZE,
        bg.color = "white", bg.r = 0.12, lineheight = 0.9,
        box.padding = 1.0, point.padding = 0.6, segment.color = "grey20", segment.size = 0.6,
        min.segment.length = 0, max.overlaps = Inf, nudge_x = lab$nx, nudge_y = 4, show.legend = FALSE)
  } else {
    p <- p + annotate("text", x = -Inf, y = Inf, hjust = -0.06, vjust = 1.6, size = 6.0,
        fontface = "italic", color = "grey25",
        label = "No sarcosine KO reached significance\n(q < 0.05 & |log2FC| > 0.5)")
  }

  out <- file.path(crc_root, COHORTS[cn], "results_kegg", "volcano_KO_paper.png")
  save_png(p, out, 12, 8.5)
  labstr <- if (nrow(lab)) paste(sprintf("%s/%s(%s)", lab$KO, lab$Short, ifelse(lab$log2FC>0,"Cancer","Healthy")), collapse=", ") else "(none)"
  cat(sprintf("%-12s KOs=%d | sig Cancer/Healthy=%d/%d | labelled sarcosine KOs: %s -> %s\n",
      cn, nrow(d), ncan, nhea, labstr, basename(out)))
}

writeLines(capture.output(sessionInfo()), file.path(script_dir, "sessionInfo_volcano_KO_paper_percohort.txt"))
cat("DONE. 4 per-cohort paper-style KO volcanoes (significant sarcosine KOs highlighted) written.\n")
