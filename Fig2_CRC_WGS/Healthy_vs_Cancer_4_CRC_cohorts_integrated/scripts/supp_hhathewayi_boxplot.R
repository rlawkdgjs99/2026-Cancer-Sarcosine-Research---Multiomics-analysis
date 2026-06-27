#!/usr/bin/env Rscript
# =============================================================================
# Healthy vs Cancer abundance of Hungatella hathewayi
# -----------------------------------------------------------------------------
# H. hathewayi is the single species that was sarcosine-PRODUCTION-associated in
# 3 of the 4 cohorts (cross_cohort_prod_assoc Venn: significant in PRJEB27928,
# PRJEB6070, PRJNA429097; rho=0.325 but not significant in PRJEB10878). It is
# also a CRC-enriched species (Goal A 4-way intersection), so its cancer
# enrichment is established independently of the production selection.
#
# Same display conventions as supp_deg_assoc_5species_boxplot.R:
#   log10 y-axis, zeros at floor = min(nonzero)/5, box+jitter, Wilcoxon on RAW.
# Single species -> single two-sided Wilcoxon test, raw P (no multiple testing).
#
# Inputs:  cross_cohort_prod_assoc_membership_matrix.csv (to fetch exact lineage),
#          sarcosine_KO_per_sample_pooled.csv (Run.ID, Group; 1647 = Venn set),
#          <cohort>/Bacteria_*.txt
# Outputs: results_integrated/cross_cohort/hhathewayi_boxplot_log10.png
#          results_integrated/cross_cohort/hhathewayi_stats.csv
# =============================================================================

options(bitmapType = "quartz")
set.seed(42)
suppressMessages({ library(dplyr); library(tidyr); library(ggplot2) })

BASE_DIR <- "Healthy_vs_Cancer_4_CRC_cohorts_integrated"   # relative to the analysis-folder root (working dir)
setwd(BASE_DIR)
COHORT_DIRS <- c("PRJEB10878_CRC", "PRJEB27928_CRC",
                 "PRJEB6070_CRC_AdenomatousPolyps", "PRJNA429097_CRC")
OUTDIR <- "results_integrated/cross_cohort"
MEMB   <- file.path(OUTDIR, "cross_cohort_prod_assoc_membership_matrix.csv")
CACHE  <- "results_integrated/sarcosine/sarcosine_KO_per_sample_pooled.csv"
stopifnot(file.exists(MEMB), file.exists(CACHE))

COL_H <- "#1B7837"; COL_C <- "#B2182B"; PSEUDO_DIV <- 5

# ---- exact lineage of H. hathewayi (verify it is the 3-cohort species) -------
memb <- read.csv(MEMB, stringsAsFactors = FALSE, check.names = FALSE)
row <- memb[memb$Species == "Hungatella_hathewayi", ]
stopifnot(nrow(row) == 1, row$n_cohorts_significant == 3)
target <- row$Species_full
cat("Target:", target, "\n  production-significant in", row$n_cohorts_significant, "cohorts\n")

# ---- group labels (1647-sample Venn set) -------------------------------------
cache <- read.csv(CACHE, stringsAsFactors = FALSE, check.names = FALSE)
grp <- setNames(cache$Group, cache$Run.ID); samples <- cache$Run.ID

# ---- species-level abundance of H. hathewayi ---------------------------------
ab <- setNames(rep(0, length(samples)), samples)
for (cd in COHORT_DIRS) {
  bf <- list.files(cd, pattern = "^Bacteria_.*\\.txt$", full.names = TRUE)[1]
  b <- read.delim(bf, stringsAsFactors = FALSE); colnames(b) <- trimws(colnames(b))
  keep <- grepl("\\|s__", b$Taxa) & !grepl("\\|t__", b$Taxa) &
          b$Taxa == target & b$Run.ID %in% samples
  b <- b[keep, , drop = FALSE]
  if (nrow(b) > 0) for (i in seq_len(nrow(b))) ab[b$Run.ID[i]] <- b$Abundance[i]
}
stopifnot(max(ab) > 0)
g <- factor(grp[names(ab)], levels = c("Healthy", "Cancer"))
h <- ab[g == "Healthy"]; c <- ab[g == "Cancer"]

# ---- stats (raw abundance) ---------------------------------------------------
p <- suppressWarnings(wilcox.test(h, c)$p.value)
stats <- data.frame(
  Species = "Hungatella_hathewayi", Species_full = target,
  n_Healthy = length(h), n_Cancer = length(c),
  prev_Healthy = round(mean(h > 0) * 100, 1), prev_Cancer = round(mean(c > 0) * 100, 1),
  median_Healthy = median(h), median_Cancer = median(c),
  mean_Healthy = mean(h), mean_Cancer = mean(c),
  log2FC = log2((mean(c) + 1e-6) / (mean(h) + 1e-6)),
  p_value = p, stringsAsFactors = FALSE)
write.csv(stats, file.path(OUTDIR, "hhathewayi_stats.csv"), row.names = FALSE)
cat("\nHealthy vs Cancer:\n"); print(stats[, c("prev_Healthy","prev_Cancer","median_Healthy",
      "median_Cancer","log2FC","p_value")], digits = 3)

p_txt <- if (p < 0.0001) "p<0.0001" else if (p < 0.001) sprintf("p=%.4f", p) else sprintf("p=%.3f", p)

# ---- plot (log10 box + jitter; zeros at floor = min(nonzero)/5) --------------
floor_v <- min(ab[ab > 0]) / PSEUDO_DIV
df <- data.frame(Group = g, Abundance = ab) %>%
  mutate(y = ifelse(Abundance <= 0, floor_v, Abundance))

pl <- ggplot(df, aes(x = Group, y = y, fill = Group)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75, width = 0.55) +
  geom_jitter(position = position_jitter(width = 0.15, seed = 42), size = 0.35, alpha = 0.12) +
  scale_y_log10() + annotation_logticks(sides = "l") +
  scale_fill_manual(values = c("Healthy" = COL_H, "Cancer" = COL_C)) +
  annotate("text", x = 1.5, y = max(df$y) * 2.2, label = p_txt, fontface = "italic", size = 6.5) +
  labs(title = expression(italic("Hungatella hathewayi")*": Healthy vs Cancer"),
       subtitle = paste0("Prod-assoc. in 3/4 cohorts; CRC-enriched in all 4  |  ",
                         "H n=", length(h), ", C n=", length(c), "  |  Wilcoxon"),
       x = "", y = "Relative Abundance (log10)") +
  theme_bw() +
  theme(plot.title = element_text(size = 21, face = "bold"),
        plot.subtitle = element_text(size = 12, color = "gray35"),
        axis.title = element_text(size = 18), axis.text = element_text(size = 16),
        legend.position = "bottom")

outfig <- file.path(OUTDIR, "hhathewayi_boxplot_log10.png")
ggsave(outfig, pl, width = 7.5, height = 6, dpi = 300, bg = "white",
       device = grDevices::png, type = "quartz")
cat("\nSaved:", outfig, "\n")
