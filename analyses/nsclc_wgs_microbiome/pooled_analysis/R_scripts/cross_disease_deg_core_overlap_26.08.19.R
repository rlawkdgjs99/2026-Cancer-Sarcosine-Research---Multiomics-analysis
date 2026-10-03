#!/usr/bin/env Rscript
# =============================================================================
# cross_disease_deg_core_overlap_26.08.19.R
# Purpose: Visualize the overlap between the CRC degradation core (Figure 1,
#   4 CRC cohorts) and the NSCLC ICI degradation core (Figure 2, 3 NSCLC
#   discovery cohorts). Two species — L. eligens and R. faecis — are shared.
#   This fact is already stated in the manuscript text (para 30) but not yet
#   shown as a single panel. This script produces:
#     (1) a 2-set Venn diagram (CRC core ∩ NSCLC core)
#     (2) a summary table (species x direction x disease)
#
# Data sources (read-only; the CSVs already exist):
#   CRC  : HGMT_CRC_WGS.../results_integrated/cross_cohort/
#            cross_cohort_deg_assoc_intersection_4cohorts.csv   (5 species)
#   NSCLC: HGMT_NSCLC_ICI.../pooled_analysis/results/
#            cross_cohort_deg_assoc_intersection_discovery_NSCLC.csv (3 species)
#
# Definition note (transparently stated in the Methods, paras 63 & 77):
#   CRC core  : Spearman rho > 0.3 AND BH-FDR p_adj < 0.05, 4 cohorts, all 4
#   NSCLC core: Spearman rho > 0.3 AND nominal p < 0.05,    3 discovery cohorts
#   The NSCLC criterion is deliberately nominal because its small cohorts are
#   underpowered for per-cohort BH (n=65, n=25); the cross-cohort overlap itself
#   serves as the consistency control. This difference is acknowledged in the
#   figure caption so the reader sees it.
#
# Output (new folder; originals untouched):
#   HGMT_NSCLC_ICI.../pooled_analysis/results/cross_disease_deg_overlap_26.08.19/
#     cross_disease_deg_core_overlap_venn.png
#     cross_disease_deg_core_overlap_table.png
#     cross_disease_deg_core_overlap_stats.csv
#     sessionInfo.txt
#
# Run from: HGMT_NSCLC_ICI_RvsNR_WGS/pooled_analysis/
# (chosen because the NSCLC core lives here; the CRC CSV is read by absolute path)
# =============================================================================

set.seed(42)
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggvenn)
  library(ragg)
})

# ---- Paths -----------------------------------------------------------------
# NSCLC pooled dir is the cwd
CRC_CSV <- file.path("..", "crc_wgs_microbiome", "Healthy_vs_Cancer_4_CRC_cohorts_integrated", "results_integrated", "cross_cohort", "cross_cohort_deg_assoc_intersection_4cohorts.csv")
NSCLC_CSV <- "results/cross_cohort_deg_assoc_intersection_discovery_NSCLC.csv"
stopifnot(file.exists(CRC_CSV), file.exists(NSCLC_CSV))

out_dir <- "results/cross_disease_deg_overlap_26.08.19"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Shared NC theme (same source as 99h/99i) -----------------------------
# Absolute path: from pooled_analysis, ../..  == 공공_Metabolomics&Metagenomics_분석모음/
.SH <- normalizePath(file.path("..", "_shared", "theme_nc_26.08.18.R"))   # run from analyses/<module>/
stopifnot(file.exists(.SH)); source(.SH)

# ---- Read the two cores ----------------------------------------------------
crc_core <- read.csv(CRC_CSV, stringsAsFactors = FALSE)
nsc_core <- read.csv(NSCLC_CSV, stringsAsFactors = FALSE)
stopifnot("Species" %in% names(crc_core), "species" %in% names(nsc_core))

# Normalize species names (both use underscored binomial like "Roseburia_faecis")
crc_species  <- unique(crc_core$Species)
nsc_species  <- unique(nsc_core$species)
cat("CRC core (5 expected):", length(crc_species), "->", paste(crc_species, collapse=", "), "\n")
cat("NSCLC core (3 expected):", length(nsc_species), "->", paste(nsc_species, collapse=", "), "\n")
stopifnot(length(crc_species) == 5L, length(nsc_species) == 3L)

# Overlap (match on the underscored binomial)
overlap <- intersect(crc_species, nsc_species)
cat("Overlap:", length(overlap), "->", paste(overlap, collapse=", "), "\n")
stopifnot(length(overlap) == 2L)
stopifnot(all(c("Lachnospira_eligens","Roseburia_faecis") %in% overlap))

# Pretty species labels (italic binomial)
pretty_sp <- function(x) {
  x <- gsub("_", " ", x)
  ifelse(x == "Clostridiales bacterium KLE1615" | x == "Lachnospiraceae bacterium AM48 27BH",
         x, paste0("italic('", x, "')"))
}
# For plain-text labels (CSV), keep underscores
# For plot labels, use italics for binomials via parse-friendly strings

# ---- Build summary table for the CSV + the table figure --------------------
# Direction in CRC: all healthy-enriched (from the CRC analysis)
# Direction in NSCLC: all responder-enriched (from the NSCLC analysis)
# Mean rho values from each CSV
crc_mean_rho <- setNames(crc_core$mean_rho, crc_core$Species)
nsc_mean_rho <- setNames(nsc_core$mean_disc_rho, nsc_core$species)

all_species <- union(crc_species, nsc_species)
tbl <- data.frame(
  Species = all_species,
  in_CRC_core = all_species %in% crc_species,
  in_NSCLC_core = all_species %in% nsc_species,
  overlap = all_species %in% overlap,
  CRC_direction = ifelse(all_species %in% crc_species, "Healthy-enriched", "\u2014"),
  NSCLC_direction = ifelse(all_species %in% nsc_species, "Responder-enriched", "\u2014"),
  CRC_mean_rho = round(crc_mean_rho[all_species], 2),
  NSCLC_mean_rho = round(nsc_mean_rho[all_species], 2),
  stringsAsFactors = FALSE,
  row.names = NULL
)
# Order: overlap first, then CRC-only, then NSCLC-only
tbl$order_key <- ifelse(tbl$overlap, 0, ifelse(tbl$in_CRC_core, 1, 2))
tbl <- tbl[order(tbl$order_key, -rowSums(cbind(!is.na(tbl$CRC_mean_rho), !is.na(tbl$NSCLC_mean_rho)))), ]
tbl$order_key <- NULL
write.csv(tbl, file.path(out_dir, "cross_disease_deg_core_overlap_stats.csv"), row.names = FALSE)
cat("\nSummary table:\n"); print(tbl, row.names = FALSE)

# ---- (1) 2-set Venn --------------------------------------------------------
# ggvenn expects a named list of character vectors
venn_sets <- list(
  "CRC 4-cohort core (n=5)" = crc_species,
  "NSCLC 3-cohort core (n=3)" = nsc_species
)
# Colors: use CRC green (healthy) and NSCLC blue (responder) from the shared theme
VENN_COLS <- c(COL_HEALTHY, COL_R)

# Render Venn at NC footprint — wider to give set names room
p_venn <- ggvenn(
  venn_sets,
  fill_color = VENN_COLS,
  fill_alpha = 0.45,
  stroke_color = "grey25",
  stroke_size = pt_lw(1.0),
  set_name_size = mm_text(NC_TITLE_PT - 0.5),
  text_size = mm_text(NC_ANNOT_PT + 1),
  show_percentage = FALSE,
  label_sep = "\n"
) +
  # Remove all axes, ticks, gridlines — Venn is a set diagram, not a plot
  theme_void(base_family = "Arial", base_size = NC_TICK_PT) +
  theme(plot.margin = margin(8, 14, 8, 14))

save_nc(p_venn, file.path(out_dir, "cross_disease_deg_core_overlap_venn.png"), 3.4, 2.6)
cat("written cross_disease_deg_core_overlap_venn.png\n")

# ---- (2) Summary table figure ----------------------------------------------
# Build a ggplot table: rows = species, columns = CRC dir / NSCLC dir / overlap
# Use italic for binomials
tbl$Species_label <- ifelse(
  tbl$Species %in% c("Clostridiales_bacterium_KLE1615", "Lachnospiraceae_bacterium_AM48_27BH",
                    "Clostridium_sp_AF36_4"),
  gsub("_", " ", tbl$Species),
  paste0("italic('", gsub("_", " ", tbl$Species), "')")
)
# Mark overlap with a filled dot
tbl$overlap_mark <- ifelse(tbl$overlap, "\u25CF", "")

# Order rows: overlap first at top
tbl$Species_label <- factor(tbl$Species_label, levels = rev(tbl$Species_label))

# Build a 4-column table: Species | CRC core | NSCLC core | Shared
# Use geom_tile + geom_text for a clean table-like figure
tbl$CRC_mark <- ifelse(tbl$in_CRC_core, "\u25CF", "")
tbl$NSCLC_mark <- ifelse(tbl$in_NSCLC_core, "\u25CF", "")

# Long format for presence dots
# Use single-word column names to prevent any overlap; full meaning in caption
long <- data.frame(
  Species = factor(rep(as.character(tbl$Species_label), each = 3),
                   levels = rev(as.character(tbl$Species_label))),
  column = factor(c("CRC", "NSCLC", "Shared"),
                  levels = c("CRC", "NSCLC", "Shared")),
  present = c(rbind(tbl$in_CRC_core, tbl$in_NSCLC_core, tbl$overlap)),
  stringsAsFactors = FALSE
)

# Color the dots by column: CRC=green, NSCLC=blue, Shared=violet
dot_colors <- c("CRC" = COL_HEALTHY, "NSCLC" = COL_R, "Shared" = COL_HIGH)

p_tbl <- ggplot(long, aes(column, Species)) +
  geom_tile(fill = "grey98", colour = "grey85", linewidth = pt_lw(0.5)) +
  geom_point(data = long[long$present, , drop = FALSE],
             aes(colour = column), shape = 16, size = 2.6) +
  scale_colour_manual(values = dot_colors, guide = "none") +
  scale_y_discrete(labels = function(x) {
    sapply(x, function(lbl) {
      if (grepl("^italic\\(", lbl)) parse(text = lbl) else lbl
    }, USE.NAMES = FALSE)
  }) +
  labs(x = NULL, y = NULL) +
  coord_fixed(ratio = 2.8) +
  theme_nc() +
  theme(axis.text.x = element_text(size = NC_TICK_PT + 0.5, colour = "black",
                                    margin = margin(t = 4)),
        axis.text.y = element_text(size = NC_TICK_PT, colour = "black",
                                    hjust = 1, margin = margin(r = 5)),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        plot.margin = margin(3, 6, 3, 6))

save_nc(p_tbl, file.path(out_dir, "cross_disease_deg_core_overlap_table.png"), 3.4, 2.6)
cat("written cross_disease_deg_core_overlap_table.png\n")

# ---- sessionInfo -----------------------------------------------------------
sink(file.path(out_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()
cat("\nDone. Outputs in:", out_dir, "\n")
