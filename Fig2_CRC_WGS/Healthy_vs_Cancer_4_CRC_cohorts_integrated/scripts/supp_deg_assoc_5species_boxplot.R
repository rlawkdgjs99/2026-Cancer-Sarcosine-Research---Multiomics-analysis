#!/usr/bin/env Rscript
# =============================================================================
# SUPPLEMENTARY figure for the degradation Venn (cross_cohort_deg_assoc_venn.png)
# -----------------------------------------------------------------------------
# Relative abundance of the 5 cross-cohort sarcosine-DEGRADATION-associated
# species (the 4-way Venn intersection), Healthy vs Cancer, as a log10 box plot
# with jittered points and per-species Wilcoxon p (BH-adjusted across the 5).
#
# WHY box(log)+points, not a mean bar:
#   These are zero-inflated, right-skewed relative abundances (35-76% detection;
#   2 of the 5 have median 0), so a mean +/- error bar would hide the
#   distribution / mislead. Box+points on log10 shows medians, spread and zeros.
#
# Conventions mirrored from regenerate_sarcosine_boxplot_log10.R (paper style):
#   - zeros displayed at a per-species floor = min(nonzero)/5
#   - scale_y_log10() + annotation_logticks(); Wilcoxon on RAW abundance
#   - paper-size fonts
# Difference: this panel uses BH-adjusted p across the 5 species (focused panel),
#   and the unified group palette (#1B7837 Healthy / #B2182B Cancer).
#   The sister log10 KO figure now uses the same unified palette (#1B7837/#B2182B).
#
# Inputs (all already on disk):
#   results_integrated/cross_cohort/cross_cohort_deg_assoc_intersection_4cohorts.csv  (the 5 species)
#   results_integrated/sarcosine/sarcosine_KO_per_sample_pooled.csv                   (Run.ID, Group; 1647 = Venn set)
#   <cohort>/Bacteria_*.txt                                                           (species relative abundance)
# Outputs:
#   results_integrated/cross_cohort/cross_cohort_deg_assoc_5species_boxplot_log10.png
#   results_integrated/cross_cohort/cross_cohort_deg_assoc_5species_stats.csv
#
# Author: added 2026-06-09 (supplementary to cross_cohort_degradation_consistency.R)
# =============================================================================

options(bitmapType = "quartz")
set.seed(42)
suppressMessages({ library(dplyr); library(tidyr); library(ggplot2) })

BASE_DIR <- "Healthy_vs_Cancer_4_CRC_cohorts_integrated"   # relative to the analysis-folder root (working dir)
setwd(BASE_DIR)

COHORT_DIRS <- c("PRJEB10878_CRC", "PRJEB27928_CRC",
                 "PRJEB6070_CRC_AdenomatousPolyps", "PRJNA429097_CRC")
OUTDIR   <- "results_integrated/cross_cohort"
INTERSECT<- file.path(OUTDIR, "cross_cohort_deg_assoc_intersection_4cohorts.csv")
CACHE    <- "results_integrated/sarcosine/sarcosine_KO_per_sample_pooled.csv"
stopifnot(file.exists(INTERSECT), file.exists(CACHE))

# ---- style ----
COL_H <- "#1B7837"; COL_C <- "#B2182B"            # unified group green/red
PSEUDO_DIV <- 5                                    # zero-floor = min(nonzero)/5
FS_TITLE <- 30; FS_SUB <- 14; FS_STRIP <- 14; FS_AXIS_T <- 20; FS_AXIS <- 16; FS_LEG <- 16

q_label <- function(q) {
  if (is.na(q)) return("")
  if (q < 0.0001) return("q<0.0001")
  if (q < 0.001)  return(sprintf("q=%.4f", q))
  sprintf("q=%.3f", q)
}

# ---- 1. targets (5 species) + group labels (1647-sample Venn set) ------------
inter  <- read.csv(INTERSECT, stringsAsFactors = FALSE, check.names = FALSE)
target <- inter$Species_full
stopifnot(length(target) == 5)
sp_short <- gsub("_", " ", inter$Species)            # display names, intersection (mean_rho) order

cache <- read.csv(CACHE, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(all(c("Run.ID", "Group") %in% colnames(cache)))
grp <- setNames(cache$Group, cache$Run.ID)
samples <- cache$Run.ID
cat("Samples:", length(samples), " (Healthy =", sum(cache$Group == "Healthy"),
    ", Cancer =", sum(cache$Group == "Cancer"), ")\n")

# ---- 2. species-level relative abundance for the 5 targets -------------------
ab <- matrix(0, nrow = length(samples), ncol = length(target),
             dimnames = list(samples, target))
for (cd in COHORT_DIRS) {
  bf <- list.files(cd, pattern = "^Bacteria_.*\\.txt$", full.names = TRUE)[1]
  stopifnot(!is.na(bf))
  b <- read.delim(bf, stringsAsFactors = FALSE); colnames(b) <- trimws(colnames(b))
  keep <- grepl("\\|s__", b$Taxa) & !grepl("\\|t__", b$Taxa) &
          b$Taxa %in% target & b$Run.ID %in% samples
  b <- b[keep, , drop = FALSE]
  if (nrow(b) > 0) for (i in seq_len(nrow(b))) ab[b$Run.ID[i], b$Taxa[i]] <- b$Abundance[i]
}
stopifnot(all(apply(ab, 2, max) > 0))               # every target detected somewhere

# ---- 3. per-species Wilcoxon (raw abundance) + BH across the 5 ---------------
g <- factor(grp[rownames(ab)], levels = c("Healthy", "Cancer"))
stats <- lapply(seq_along(target), function(j) {
  h <- ab[g == "Healthy", j]; c <- ab[g == "Cancer", j]
  data.frame(
    Species = inter$Species[j], Species_full = target[j],
    n_Healthy = length(h), n_Cancer = length(c),
    prev_Healthy = round(mean(h > 0) * 100, 1), prev_Cancer = round(mean(c > 0) * 100, 1),
    median_Healthy = median(h), median_Cancer = median(c),
    mean_Healthy = mean(h), mean_Cancer = mean(c),
    log2FC = log2((mean(c) + 1e-6) / (mean(h) + 1e-6)),
    p_value = suppressWarnings(wilcox.test(h, c)$p.value),
    stringsAsFactors = FALSE)
})
stats <- bind_rows(stats)
stats$p_adj <- p.adjust(stats$p_value, method = "BH")
write.csv(stats, file.path(OUTDIR, "cross_cohort_deg_assoc_5species_stats.csv"), row.names = FALSE)
cat("\nPer-species Healthy vs Cancer (Wilcoxon, BH across 5):\n")
print(stats[, c("Species", "prev_Healthy", "prev_Cancer", "log2FC", "p_value", "p_adj")], digits = 3)

# ---- 4. long format; log10 zero-floor = min(nonzero)/5 per species -----------
long <- as.data.frame(ab) %>%
  mutate(Run.ID = rownames(ab), Group = g) %>%
  pivot_longer(all_of(target), names_to = "Species_full", values_to = "Abundance") %>%
  left_join(data.frame(Species_full = target, Species = sp_short, stringsAsFactors = FALSE),
            by = "Species_full") %>%
  mutate(Species = factor(Species, levels = sp_short))    # intersection (mean_rho) order
long <- long %>% group_by(Species) %>%
  mutate(floor_k = min(Abundance[Abundance > 0], na.rm = TRUE) / PSEUDO_DIV,
         y = ifelse(Abundance <= 0, floor_k, Abundance)) %>% ungroup()

# p-value annotation positions (top of each panel)
ann <- stats %>%
  mutate(Species = factor(gsub("_", " ", Species), levels = sp_short),
         qtext = sapply(p_adj, q_label)) %>%
  left_join(long %>% group_by(Species) %>% summarise(ymax = max(y, na.rm = TRUE), .groups = "drop"),
            by = "Species")

# ---- 5. plot -----------------------------------------------------------------
p <- ggplot(long, aes(x = Group, y = y, fill = Group)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75, width = 0.6) +
  geom_jitter(position = position_jitter(width = 0.15, seed = 42), size = 0.3, alpha = 0.12) +
  facet_wrap(~Species, scales = "free_y", nrow = 1) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  scale_fill_manual(values = c("Healthy" = COL_H, "Cancer" = COL_C)) +
  geom_text(data = ann, aes(x = 1.5, y = ymax * 2.2, label = qtext),
            inherit.aes = FALSE, size = 5.2, fontface = "italic") +
  labs(title = "Sarcosine-degradation-associated species (4-cohort Venn intersection): Healthy vs Cancer",
       subtitle = paste0("Relative abundance (log10; zeros at per-species floor = min(nonzero)/5)  |  ",
                         "Healthy n=", sum(g == "Healthy"), ", Cancer n=", sum(g == "Cancer"),
                         "  |  Wilcoxon, BH-adjusted across 5 species"),
       x = "", y = "Relative Abundance (log10)") +
  theme_bw() +
  theme(plot.title = element_text(size = FS_TITLE, face = "bold"),
        plot.subtitle = element_text(size = FS_SUB, color = "gray35"),
        strip.text = element_text(size = FS_STRIP, face = "italic"),
        axis.title = element_text(size = FS_AXIS_T),
        axis.text = element_text(size = FS_AXIS),
        legend.title = element_text(size = FS_LEG), legend.text = element_text(size = FS_LEG),
        legend.position = "bottom")

outfig <- file.path(OUTDIR, "cross_cohort_deg_assoc_5species_boxplot_log10.png")
ggsave(outfig, p, width = 24, height = 5.5, dpi = 300, bg = "white",
       device = grDevices::png, type = "quartz")
cat("\nSaved figure:", outfig, "\n")
cat("Saved stats :", file.path(OUTDIR, "cross_cohort_deg_assoc_5species_stats.csv"), "\n")
cat("R:", R.version.string, "| ggplot2", as.character(packageVersion("ggplot2")), "\n")
