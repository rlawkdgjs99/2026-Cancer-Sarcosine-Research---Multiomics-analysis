#!/usr/bin/env Rscript
# =============================================================================
# SUPPLEMENTARY figure for the CRC-enriched Venn (cross_cohort_CRC_enriched_venn.png)
# -----------------------------------------------------------------------------
# Relative abundance of the 3 cross-cohort CRC-ENRICHED species (the 4-way Venn
# intersection of panel d: per-cohort p_adj<0.05 & log2FC>1 in ALL 4 cohorts),
# Healthy vs Cancer, as a log10 box plot with jittered points and per-species
# Wilcoxon p (BH-adjusted across the 3).
#
# SISTER of supp_deg_assoc_5species_boxplot.R (which shows the 5 DEGRADATION-
# associated/Venn-intersection species, the depleted/healthy-enriched set).
# This one shows the OPPOSITE pole: the cancer-enriched intersection species
# (Parvimonas micra, Hungatella hathewayi, Clostridium symbiosum). Unlike the
# degradation panel (species selected by score correlation, not by DA), here the
# 3 species ARE the DA-significant 4-way intersection itself, so the box plot
# shows the selected species directly. The pooled Wilcoxon below is DESCRIPTIVE;
# the formal R-vs-NR... (here Healthy-vs-Cancer) inference is the per-cohort DA +
# 4-way intersection already shown in panels c/d.
#
# WHY box(log)+points, not a mean bar: zero-inflated, right-skewed relative
#   abundances -> a mean +/- error bar would mislead. Box+points on log10 shows
#   medians, spread and zeros. Conventions mirrored from the sister script
#   (zeros at per-species floor = min(nonzero)/5; Wilcoxon on RAW abundance;
#   unified palette #1B7837 Healthy / #B2182B Cancer; paper-size fonts).
#
# Inputs (all already on disk):
#   results_integrated/cross_cohort/cross_cohort_CRC_enriched_intersection_4cohorts.csv  (the 3 species)
#   results_integrated/sarcosine/sarcosine_KO_per_sample_pooled.csv                      (Run.ID, Group)
#   <cohort>/Bacteria_*.txt                                                              (species relative abundance)
# Outputs:
#   results_integrated/cross_cohort/cross_cohort_CRC_enriched_3species_boxplot_log10.png
#   results_integrated/cross_cohort/cross_cohort_CRC_enriched_3species_stats.csv
#
# Author: added 2026-06-26 (supplementary companion to cross_cohort_consistency.R, panel d)
# =============================================================================

options(bitmapType = "quartz")
set.seed(42)
suppressMessages({ library(dplyr); library(tidyr); library(ggplot2) })

BASE_DIR <- "Healthy_vs_Cancer_4_CRC_cohorts_integrated"   # relative to the analysis-folder root (working dir)
setwd(BASE_DIR)

COHORT_DIRS <- c("PRJEB10878_CRC", "PRJEB27928_CRC",
                 "PRJEB6070_CRC_AdenomatousPolyps", "PRJNA429097_CRC")
OUTDIR   <- "results_integrated/cross_cohort"
INTERSECT<- file.path(OUTDIR, "cross_cohort_CRC_enriched_intersection_4cohorts.csv")
stopifnot(file.exists(INTERSECT))

# ---- style (identical to the sister figure) ----
COL_H <- "#1B7837"; COL_C <- "#B2182B"            # unified group green/red
PSEUDO_DIV <- 5                                    # zero-floor = min(nonzero)/5
FS_TITLE <- 30; FS_SUB <- 14; FS_STRIP <- 15; FS_AXIS_T <- 20; FS_AXIS <- 16; FS_LEG <- 16
# Title kept LARGE (30) for readability; shortened to a concise core so it fits width 15.

q_label <- function(q) {
  if (is.na(q)) return("")
  if (q < 0.0001) return("q<0.0001")
  if (q < 0.001)  return(sprintf("q=%.4f", q))
  sprintf("q=%.3f", q)
}

# ---- 1. targets (3 species, by descending mean cross-cohort log2FC) + group labels
inter  <- read.csv(INTERSECT, stringsAsFactors = FALSE, check.names = FALSE)
inter  <- inter[order(-inter$mean_log2FC), ]        # display order: strongest cancer-enrichment first
target <- inter$Species_full
stopifnot(length(target) == 3)
sp_short <- gsub("_", " ", inter$Species)           # display names

# group labels from cohort metadata = the SAME n=1,649 taxonomic set as panels a-d
# (selected_project_*.txt: keep Health + the cohort's cancer phenotype, WGS only) — NOT
# the smaller KO-functional subset — so panel e's n matches the rest of Supplementary Figure 2.
read_groups <- function(cd) {
  mf <- sort(list.files(cd, "^selected_project_.*\\.txt$", full.names = TRUE), decreasing = TRUE)[1]
  ml <- readLines(mf); hdr <- strsplit(ml[2], "\t")[[1]]
  dl <- ml[3:length(ml)]; dl <- dl[dl != ""]
  rows <- lapply(strsplit(dl, "\t"), function(x)
    if (length(x) >= length(hdr)) x[seq_along(hdr)] else c(x, rep(NA, length(hdr) - length(x))))
  m <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE); colnames(m) <- gsub(" ", ".", hdr)
  cl <- setdiff(unique(m$Phenotype.name), c("Health", "Adenomatous Polyps"))[1]
  m <- m[m$Phenotype.name %in% c("Health", cl), , drop = FALSE]
  if ("Assay.type" %in% colnames(m)) m <- m[m$Assay.type == "WGS", , drop = FALSE]
  data.frame(RunID = m$Run.ID, Group = ifelse(m$Phenotype.name == "Health", "Healthy", "Cancer"),
             stringsAsFactors = FALSE)
}
meta <- do.call(rbind, lapply(COHORT_DIRS, read_groups))
meta <- meta[!is.na(meta$RunID) & meta$RunID != "", , drop = FALSE]
grp_all <- setNames(meta$Group, meta$RunID)

# ---- 2. species-level relative abundance for the 3 targets; roster = metadata ∩ bacteria
long_rows <- list(); present <- character(0)
for (cd in COHORT_DIRS) {
  bf <- list.files(cd, pattern = "^Bacteria_.*\\.txt$", full.names = TRUE)[1]
  stopifnot(!is.na(bf))
  b <- read.delim(bf, stringsAsFactors = FALSE); colnames(b) <- trimws(colnames(b))
  sp <- b[grepl("\\|s__", b$Taxa) & !grepl("\\|t__", b$Taxa), , drop = FALSE]
  present <- union(present, unique(sp$Run.ID))
  long_rows[[cd]] <- sp[sp$Taxa %in% target, c("Run.ID", "Taxa", "Abundance"), drop = FALSE]
}
samples <- intersect(names(grp_all), present)       # = the n=1,649 taxonomic set (matches a-d)
grp <- grp_all[samples]
cat("Samples:", length(samples), " (Healthy =", sum(grp == "Healthy"),
    ", Cancer =", sum(grp == "Cancer"), ")\n")

ab <- matrix(0, nrow = length(samples), ncol = length(target), dimnames = list(samples, target))
allrows <- do.call(rbind, long_rows); allrows <- allrows[allrows$Run.ID %in% samples, , drop = FALSE]
if (nrow(allrows) > 0) for (i in seq_len(nrow(allrows))) ab[allrows$Run.ID[i], allrows$Taxa[i]] <- allrows$Abundance[i]
stopifnot(all(apply(ab, 2, max) > 0))               # every target detected somewhere

# ---- 3. per-species Wilcoxon (raw abundance) + BH across the 3 ---------------
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
write.csv(stats, file.path(OUTDIR, "cross_cohort_CRC_enriched_3species_stats.csv"), row.names = FALSE)
cat("\nPer-species Healthy vs Cancer (Wilcoxon, BH across 3):\n")
print(stats[, c("Species", "prev_Healthy", "prev_Cancer", "log2FC", "p_value", "p_adj")], digits = 3)

# ---- 4. long format; log10 zero-floor = min(nonzero)/5 per species -----------
long <- as.data.frame(ab) %>%
  mutate(Run.ID = rownames(ab), Group = g) %>%
  pivot_longer(all_of(target), names_to = "Species_full", values_to = "Abundance") %>%
  left_join(data.frame(Species_full = target, Species = sp_short, stringsAsFactors = FALSE),
            by = "Species_full") %>%
  mutate(Species = factor(Species, levels = sp_short))
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
  labs(title = "CRC-enriched species (4-cohort intersection)",
       subtitle = paste0("Relative abundance (log10)  |  Healthy n=", sum(g == "Healthy"),
                         " vs Cancer n=", sum(g == "Cancer")),
       x = "", y = "Relative Abundance (log10)") +
  theme_bw() +
  theme(plot.title = element_text(size = FS_TITLE, face = "bold"),
        plot.subtitle = element_text(size = FS_SUB, color = "gray35"),
        strip.text = element_text(size = FS_STRIP, face = "italic"),
        axis.title = element_text(size = FS_AXIS_T),
        axis.text = element_text(size = FS_AXIS),
        legend.title = element_text(size = FS_LEG), legend.text = element_text(size = FS_LEG),
        legend.position = "bottom",
        plot.margin = margin(5.5, 16, 5.5, 5.5))   # extra right padding so title/subtitle never touch the edge

outfig <- file.path(OUTDIR, "cross_cohort_CRC_enriched_3species_boxplot_log10.png")
ggsave(outfig, p, width = 15, height = 5.5, dpi = 300, bg = "white",
       device = grDevices::png, type = "quartz")
cat("\nSaved figure:", outfig, "\n")
cat("Saved stats :", file.path(OUTDIR, "cross_cohort_CRC_enriched_3species_stats.csv"), "\n")
cat("R:", R.version.string, "| ggplot2", as.character(packageVersion("ggplot2")), "\n")
