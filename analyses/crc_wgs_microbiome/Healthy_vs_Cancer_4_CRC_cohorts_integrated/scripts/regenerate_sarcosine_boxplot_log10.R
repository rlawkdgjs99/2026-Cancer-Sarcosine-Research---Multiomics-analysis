###############################################################################
# Sarcosine enzyme KO boxplot — log10 y-axis, well-detected KOs (soxB, creatinase)
# Paper style: larger fonts. Reads a small cached per-sample CSV (built once from
# jq-extracted sarcosine KO records) so re-styling is fast.
#
# Build inputs (created once via jq, see chat):
#   results_integrated/sarcosine/.sarc_records.tsv   (KO  run_id  abundance)
#   results_integrated/sarcosine/.ko_samples.txt     (all KO-matrix run_ids)
# Cache:
#   results_integrated/sarcosine/sarcosine_KO_per_sample_pooled.csv
# Output (replaces the filtered boxplot; log10 + paper fonts, soxB & creatinase):
#   results_integrated/sarcosine/sarcosine_enzyme_boxplot_filtered_pooled.png
# NOTE: the heavy pipeline (Part 3b / regenerate_diversity_plots.R) was NOT changed,
#       so re-running it will overwrite this file with the old linear version.
#       Re-run THIS script to restore the log10 paper figure (uses the cached CSV).
###############################################################################

options(bitmapType = "quartz")
library(dplyr); library(tidyr); library(ggplot2)

BASE_DIR <- "Healthy_vs_Cancer_4_CRC_cohorts_integrated"   # relative to the analysis-folder root (working dir)
setwd(BASE_DIR)
SARC_DIR <- "results_integrated/sarcosine"
CACHE <- file.path(SARC_DIR, "sarcosine_KO_per_sample_pooled.csv")

COHORT_DIRS <- list(
  PRJEB6070   = "PRJEB6070_CRC_AdenomatousPolyps",
  PRJNA429097 = "PRJNA429097_CRC",
  PRJEB10878  = "PRJEB10878_CRC",
  PRJEB27928  = "PRJEB27928_CRC")

# ---- CONFIG (style) --------------------------------------------------------
COL_H <- "#1B7837"; COL_C <- "#B2182B"      # green Healthy / red Cancer
FS_TITLE <- 23; FS_SUB <- 17; FS_STRIP <- 16; FS_AXIS_T <- 21; FS_AXIS <- 17; FS_LEG <- 17
PLOT_KOS <- c("K00303", "K08688")           # well-detected: soxB, creatinase
PSEUDO_DIV <- 5                              # per-KO zero-floor = min(nonzero)/PSEUDO_DIV

format_pval <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.0001) return("p<0.0001")
  if (p < 0.001)  return(sprintf("p=%.4f", p))
  return(sprintf("p=%.3f", p))
}

# ---- Metadata (pooled Healthy/Cancer, WGS) ---------------------------------
meta_pooled <- data.frame()
for (cid in names(COHORT_DIRS)) {
  mf <- list.files(COHORT_DIRS[[cid]], pattern = "^selected_project_.*\\.txt$", full.names = TRUE)
  mf <- sort(mf, decreasing = TRUE)[1]
  ml <- readLines(mf); hdr <- strsplit(ml[2], "\t")[[1]]
  dl <- ml[3:length(ml)]; dl <- dl[dl != ""]
  dl <- lapply(strsplit(dl, "\t"), function(x)
    if (length(x) >= length(hdr)) x[1:length(hdr)] else c(x, rep(NA, length(hdr) - length(x))))
  m <- as.data.frame(do.call(rbind, dl), stringsAsFactors = FALSE); colnames(m) <- gsub(" ", ".", hdr)
  cl <- unique(m$Phenotype.name); cl <- cl[!cl %in% c("Health", "Adenomatous Polyps")][1]
  m <- m %>% filter(Phenotype.name %in% c("Health", cl)) %>%
    mutate(Group = ifelse(Phenotype.name == "Health", "Healthy", "Cancer"))
  if ("Assay.type" %in% colnames(m)) m <- m %>% filter(Assay.type == "WGS")
  m$Cohort <- cid
  meta_pooled <- rbind(meta_pooled, m[, c("Run.ID", "Group", "Cohort")])
}

# ---- Build per-sample table (cache) ----------------------------------------
if (file.exists(CACHE)) {
  per_sample <- read.csv(CACHE, stringsAsFactors = FALSE, check.names = FALSE)
  cat("Loaded cache:", CACHE, "(", nrow(per_sample), "samples )\n")
} else {
  rec <- read.delim(file.path(SARC_DIR, ".sarc_records.tsv"), header = FALSE,
                    col.names = c("KO", "Run.ID", "Abundance"), stringsAsFactors = FALSE)
  ko_samples <- readLines(file.path(SARC_DIR, ".ko_samples.txt"))
  pooled <- intersect(ko_samples, meta_pooled$Run.ID)
  cat("KO-samples total:", length(ko_samples), "| pooled (Healthy/Cancer):", length(pooled), "\n")

  rec <- rec[rec$Run.ID %in% pooled, ]
  kos <- sort(unique(rec$KO))
  mat <- matrix(0, nrow = length(pooled), ncol = length(kos), dimnames = list(pooled, kos))
  mat[cbind(rec$Run.ID, rec$KO)] <- rec$Abundance
  per_sample <- data.frame(Run.ID = pooled, mat, check.names = FALSE, stringsAsFactors = FALSE)
  per_sample <- merge(per_sample, meta_pooled, by = "Run.ID")
  write.csv(per_sample, CACHE, row.names = FALSE)
  cat("Cached:", CACHE, "\n")
}

# ---- Long format for the two plotted KOs -----------------------------------
panel <- data.frame(
  KO = c("K00303", "K08688"),
  Enzyme = c("Sarcosine oxidase subunit beta (soxB)", "Creatinase"),
  EC = c("EC:1.5.3.1/1.5.3.24", "EC:3.5.3.3"),
  Role = c("Degradation", "Production"), stringsAsFactors = FALSE)

long <- per_sample %>%
  select(Run.ID, Group, all_of(PLOT_KOS)) %>%
  pivot_longer(all_of(PLOT_KOS), names_to = "KO", values_to = "Abundance") %>%
  left_join(panel, by = "KO") %>%
  mutate(Group = factor(Group, levels = c("Healthy", "Cancer")),
         Label = paste0(KO, ": ", Enzyme, "\n(", Role, ") [", EC, "]"))

# Wilcoxon on RAW abundance (per KO)
pv <- long %>% group_by(KO, Label) %>%
  summarise(p = tryCatch(wilcox.test(Abundance ~ Group)$p.value, error = function(e) NA),
            .groups = "drop") %>%
  mutate(ptext = sapply(p, format_pval))
cat("p-values:\n"); for (i in 1:nrow(pv)) cat(" ", pv$KO[i], pv$ptext[i], "\n")

# log10 display: zeros -> per-KO floor (that KO's min nonzero / PSEUDO_DIV)
long <- long %>% group_by(KO) %>%
  mutate(pseudo_k = min(Abundance[Abundance > 0], na.rm = TRUE) / PSEUDO_DIV,
         y = ifelse(Abundance <= 0, pseudo_k, Abundance)) %>% ungroup()
cat("per-KO zero-floors:\n"); print(as.data.frame(distinct(long, KO, pseudo_k)))

ypos <- long %>% group_by(Label) %>% summarise(ymax = max(y, na.rm = TRUE), .groups = "drop")
pv <- merge(pv, ypos, by = "Label")

# ---- Plot ------------------------------------------------------------------
p <- ggplot(long, aes(x = Group, y = y, fill = Group)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.75, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.35, alpha = 0.12) +
  facet_wrap(~Label, scales = "free_y", nrow = 1) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  scale_fill_manual(values = c("Healthy" = COL_H, "Cancer" = COL_C)) +
  geom_text(data = pv, aes(x = 1.5, y = ymax * 1.6, label = ptext),
            inherit.aes = FALSE, size = 5.9, fontface = "italic") +
  labs(title = "Integrated 4 CRC Cohorts: Sarcosine Enzyme KOs (log10 scale)",
       subtitle = "Healthy vs Cancer | well-detected KOs (prevalence: soxB 92.8%, creatinase 55.5%) | zeros shown at floor",
       x = "", y = "Relative Abundance (log10)") +
  theme_bw() +
  theme(plot.title = element_text(size = FS_TITLE, face = "bold"),
        plot.subtitle = element_text(size = FS_SUB, color = "gray35"),
        strip.text = element_text(size = FS_STRIP, face = "bold"),
        axis.title = element_text(size = FS_AXIS_T),
        axis.text = element_text(size = FS_AXIS),
        legend.title = element_text(size = FS_LEG), legend.text = element_text(size = FS_LEG),
        legend.position = "bottom")

ggsave(file.path(SARC_DIR, "sarcosine_enzyme_boxplot_filtered_pooled.png"), p,
       width = 11, height = 6.5, dpi = 300, bg = "white")
cat("Saved:", file.path(SARC_DIR, "sarcosine_enzyme_boxplot_filtered_pooled.png"), "\n")
