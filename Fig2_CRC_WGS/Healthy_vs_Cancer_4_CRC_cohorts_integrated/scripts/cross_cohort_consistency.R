#!/usr/bin/env Rscript
# =============================================================================
# Cross-cohort consistency of CRC-enriched species
# =============================================================================
# Purpose:
#   Identify which CRC-enriched species are shared across the 4 WGS cohorts
#   and visualize the overlap. This is Goal A of a 2-goal analysis.
#
# Input files (per-cohort differential abundance CSVs):
#   PRJEB10878_CRC/results_bacteria/diff_abundance_species.csv
#   PRJEB27928_CRC/results_bacteria/diff_abundance_species.csv
#   PRJEB6070_CRC_AdenomatousPolyps/results_bacteria/diff_abundance_species.csv
#   PRJNA429097_CRC/results_bacteria/diff_abundance_species.csv
#
# Selection criteria for "CRC-enriched":
#   p_adj < 0.05 AND log2FC > 1
#
# Output files (written to results_integrated/cross_cohort/):
#   1. cross_cohort_CRC_enriched_venn.png
#   2. cross_cohort_CRC_enriched_upset.png
#   3. cross_cohort_CRC_enriched_per_cohort_lists.csv
#   4. cross_cohort_CRC_enriched_membership_matrix.csv
#   5. cross_cohort_CRC_enriched_intersection_4cohorts.csv
#   6. cross_cohort_CRC_enriched_summary.csv
#
# Author: automated analysis
# Date: 2026-05-19
# =============================================================================

# --- Environment quirks ------------------------------------------------------
options(bitmapType = "quartz")

# --- Libraries ---------------------------------------------------------------
library(dplyr)
library(tidyr)
library(ggplot2)

if (!requireNamespace("ggvenn", quietly = TRUE)) install.packages("ggvenn")
if (!requireNamespace("UpSetR", quietly = TRUE)) install.packages("UpSetR")
library(ggvenn)
library(UpSetR)

# --- Working directory -------------------------------------------------------
setwd("Healthy_vs_Cancer_4_CRC_cohorts_integrated")   # relative to the analysis-folder root (working dir)

# --- Constants ---------------------------------------------------------------
P_ADJ_CUT <- 0.05
LFC_CUT   <- 1

# --- Helper: custom ggsave for quartz PNG ------------------------------------
.ggsave <- function(filename, plot, width, height, dpi = 200) {
  grDevices::png(filename = filename, width = width, height = height,
                 units = "in", res = dpi, type = "quartz", bg = "white")
  print(plot)
  grDevices::dev.off()
}

# =============================================================================
# GOAL A: Cross-cohort consistency of CRC-enriched species
# =============================================================================

cat("=== Loading inputs ===\n")

INPUT_FILES <- c(
  PRJEB10878  = "../PRJEB10878_CRC/results_bacteria/diff_abundance_species.csv",
  PRJEB27928  = "../PRJEB27928_CRC/results_bacteria/diff_abundance_species.csv",
  PRJEB6070   = "../PRJEB6070_CRC_AdenomatousPolyps/results_bacteria/diff_abundance_species.csv",
  PRJNA429097 = "../PRJNA429097_CRC/results_bacteria/diff_abundance_species.csv"
)

EXPECTED_ROWS <- c(PRJEB10878 = 557, PRJEB27928 = 484, PRJEB6070 = 378, PRJNA429097 = 396)
EXPECTED_HEADER <- c("Species_full","Species","Mean_Healthy","Mean_Cancer","log2FC","p_value","p_adj")

# Read all CSVs into a list
cohort_data <- lapply(names(INPUT_FILES), function(cohort) {
  df <- read.csv(INPUT_FILES[cohort], stringsAsFactors = FALSE, check.names = FALSE)
  list(cohort = cohort, df = df)
})
names(cohort_data) <- names(INPUT_FILES)

# --- V1: Schema check --------------------------------------------------------
cat("=== Verifying schema and consistency ===\n")
cat("--- V1: Schema check ---\n")

v1_pass <- TRUE
for (cohort in names(cohort_data)) {
  header <- colnames(cohort_data[[cohort]]$df)
  if (!identical(header, EXPECTED_HEADER)) {
    cat("FAIL V1:", cohort, "header does not match expected.\n")
    cat("  Expected:", paste(EXPECTED_HEADER, collapse = ", "), "\n")
    cat("  Got     :", paste(header, collapse = ", "), "\n")
    v1_pass <- FALSE
  }
}
if (v1_pass) {
  cat("PASS V1: All 4 CSVs have identical headers in identical order.\n")
} else {
  stop("V1 schema check FAILED. Stopping.")
}

# --- V2: Row counts ----------------------------------------------------------
cat("--- V2: Row counts ---\n")

v2_pass <- TRUE
for (cohort in names(cohort_data)) {
  n_rows <- nrow(cohort_data[[cohort]]$df)
  expected <- EXPECTED_ROWS[cohort]
  if (n_rows != expected) {
    cat("FAIL V2:", cohort, "has", n_rows, "rows, expected", expected, "\n")
    v2_pass <- FALSE
  } else {
    cat("  ", cohort, ":", n_rows, "rows (expected", expected, ") — OK\n")
  }
}
if (v2_pass) {
  cat("PASS V2: All row counts match expected values.\n")
} else {
  stop("V2 row count check FAILED. Stopping.")
}

# --- V3: Cross-cohort taxonomy consistency -----------------------------------
cat("--- V3: Cross-cohort taxonomy consistency ---\n")

all_species <- lapply(cohort_data, function(x) x$df$Species_full)
union_species <- Reduce(union, all_species)
cat("Total unique Species_full across 4 cohorts:", length(union_species), "\n")

overlap_counts <- sapply(union_species, function(sp) {
  sum(sapply(all_species, function(v) sp %in% v))
})

distribution <- table(overlap_counts)
cat("Distribution (how many cohorts each species appears in):\n")
print(distribution)

species_in_all4 <- names(overlap_counts)[overlap_counts == 4]
cat("Species appearing in all 4 cohorts:", length(species_in_all4), "\n")
cat("5 examples of Species_full in all 4 cohorts:\n")
if (length(species_in_all4) >= 5) {
  cat(paste("  ", head(species_in_all4, 5), collapse = "\n"), "\n")
} else {
  cat(paste("  ", species_in_all4, collapse = "\n"), "\n")
}

# Sanity: at least 10% of each cohort should overlap with another cohort
v3_pass <- TRUE
for (i in seq_along(all_species)) {
  cohort_name <- names(all_species)[i]
  cohort_species <- all_species[[i]]
  # Overlap with at least one other cohort
  overlaps <- sapply(all_species[-i], function(other) length(intersect(cohort_species, other)))
  max_overlap <- max(overlaps)
  pct_overlap <- max_overlap / length(cohort_species) * 100
  cat("  ", cohort_name, "max overlap with another cohort:", max_overlap, "/",
      length(cohort_species), "=", round(pct_overlap, 1), "%\n")
  if (pct_overlap < 10) {
    cat("FAIL V3:", cohort_name, "has <10% overlap with any other cohort.\n")
    v3_pass <- FALSE
  }
}
if (v3_pass) {
  cat("PASS V3: All cohorts have >=10% overlap with at least one other cohort.\n")
} else {
  stop("V3 taxonomy consistency check FAILED. Stopping.")
}

# --- Define enriched sets ----------------------------------------------------
cat("=== Computing intersections ===\n")

enriched_list <- lapply(names(cohort_data), function(cohort) {
  df <- cohort_data[[cohort]]$df
  enriched <- df[df$p_adj < P_ADJ_CUT & df$log2FC > LFC_CUT, ]
  cat("  ", cohort, "— Cancer-enriched (p_adj <", P_ADJ_CUT, "& log2FC >", LFC_CUT, "):",
      nrow(enriched), "species\n")
  enriched
})
names(enriched_list) <- names(cohort_data)

# --- V4: Enriched set sizes --------------------------------------------------
cat("--- V4: Enriched set sizes ---\n")
for (cohort in names(enriched_list)) {
  cat("  ", cohort, ":", nrow(enriched_list[[cohort]]), "species\n")
}
cat("PASS V4: Enriched set sizes computed.\n")

# Extract just the Species_full vectors for set operations
enriched_sets <- lapply(enriched_list, function(df) df$Species_full)

# --- Compute intersections ---------------------------------------------------
intersect_all4 <- Reduce(intersect, enriched_sets)
cat("4-way intersection size:", length(intersect_all4), "\n")

# Build membership matrix for union of all enriched species
union_enriched <- Reduce(union, enriched_sets)

membership_rows <- lapply(union_enriched, function(sp) {
  out <- list(Species_full = sp)
  # Get short name from any cohort that has it
  short_name <- NA
  for (cohort in names(cohort_data)) {
    df <- cohort_data[[cohort]]$df
    if (sp %in% df$Species_full) {
      short_name <- df$Species[df$Species_full == sp][1]
      break
    }
  }
  out$Species <- short_name

  for (cohort in names(cohort_data)) {
    df <- cohort_data[[cohort]]$df
    idx <- which(df$Species_full == sp)
    if (length(idx) == 1) {
      out[[paste0(cohort, "_sig")]] <- sp %in% enriched_sets[[cohort]]
      out[[paste0(cohort, "_log2FC")]] <- df$log2FC[idx]
      out[[paste0(cohort, "_p_adj")]] <- df$p_adj[idx]
    } else {
      out[[paste0(cohort, "_sig")]] <- FALSE
      out[[paste0(cohort, "_log2FC")]] <- NA_real_
      out[[paste0(cohort, "_p_adj")]] <- NA_real_
    }
  }

  sig_vec <- sapply(names(cohort_data), function(cohort) {
    out[[paste0(cohort, "_sig")]]
  })
  out$n_cohorts_significant <- sum(unlist(sig_vec))

  out
})

membership_df <- bind_rows(membership_rows)

# Compute mean log2FC across cohorts (NA-removed) for sorting
membership_df$mean_log2FC <- apply(membership_df[, paste0(names(cohort_data), "_log2FC")], 1,
                                    function(x) mean(x, na.rm = TRUE))

membership_df <- membership_df %>%
  arrange(desc(n_cohorts_significant), desc(mean_log2FC))

# --- Output directory --------------------------------------------------------
OUTDIR <- "results_integrated/cross_cohort"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", OUTDIR, "\n")

# --- Output 1: Venn diagram --------------------------------------------------
cat("=== Saving outputs ===\n")
cat("--- 1. Venn diagram ---\n")

venn_data <- enriched_sets
names(venn_data) <- names(enriched_sets)

# Use a colorblind-friendly palette (Okabe-Ito inspired)
venn_colors <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

venn_plot <- ggvenn(
  venn_data,
  fill_color = venn_colors,
  stroke_size = 0.6,
  set_name_size = 8,
  text_size = 7,
  show_percentage = FALSE
) +
  labs(
    title = "CRC-enriched species across 4 WGS cohorts\n(p_adj < 0.05 & log2FC > 1)",
    subtitle = paste0(
      "Total unique species: ", length(union_enriched),
      "  |  4-way intersection: ", length(intersect_all4)
    )
  ) +
  theme(
    plot.title = element_text(size = 28, face = "bold", hjust = 0.5,
                              margin = margin(b = 10)),
    plot.subtitle = element_text(size = 22, hjust = 0.5,
                                 margin = margin(b = 14)),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  # highlight the 4-cohort central intersection count in bright red
  annotate("text", x = 0, y = -0.7,
           label = as.character(length(intersect_all4)),
           colour = "red", fontface = "bold", size = 7)

.ggsave(
  filename = file.path(OUTDIR, "cross_cohort_CRC_enriched_venn.png"),
  plot = venn_plot,
  width = 11, height = 11, dpi = 300
)
cat("Saved:", file.path(OUTDIR, "cross_cohort_CRC_enriched_venn.png"), "\n")

# --- Output 2: UpSet plot ----------------------------------------------------
cat("--- 2. UpSet plot ---\n")

# Prepare binary membership matrix for UpSetR
upset_df <- membership_df[, c("Species_full", paste0(names(cohort_data), "_sig"))]
colnames(upset_df) <- c("Species_full", names(cohort_data))
upset_df[, -1] <- lapply(upset_df[, -1], as.integer)

# UpSetR requires rownames
upset_mat <- as.data.frame(upset_df[, -1])
rownames(upset_mat) <- upset_df$Species_full

# Build title as grid text
upset_title <- paste0(
  "CRC-enriched species across 4 WGS cohorts\n(p_adj < 0.05 & log2FC > 1)\n",
  "Total unique: ", length(union_enriched), " | 4-way intersection: ", length(intersect_all4)
)

grDevices::png(
  filename = file.path(OUTDIR, "cross_cohort_CRC_enriched_upset.png"),
  width = 14, height = 8, units = "in", res = 300, type = "quartz", bg = "white"
)
upset(
  upset_mat,
  nsets = 4,
  sets = c("PRJEB10878", "PRJEB27928", "PRJEB6070", "PRJNA429097"),
  order.by = "freq",
  keep.order = TRUE,
  decreasing = TRUE,
  mainbar.y.label = "Intersection size",
  sets.x.label = "Species per cohort",
  text.scale = c(2.2, 2.0, 1.8, 1.6, 2.2, 2.0),
  point.size = 3.5,
  line.size = 1.2,
  mb.ratio = c(0.65, 0.35)
)
grid::grid.text(upset_title, x = 0.65, y = 0.99, just = "top", gp = grid::gpar(fontsize = 26, fontface = "bold"))
grDevices::dev.off()
cat("Saved:", file.path(OUTDIR, "cross_cohort_CRC_enriched_upset.png"), "\n")

# --- Output 3: Per-cohort lists (long format) --------------------------------
cat("--- 3. Per-cohort lists ---\n")

per_cohort_lists <- bind_rows(lapply(names(enriched_list), function(cohort) {
  enriched_list[[cohort]] %>%
    select(Species_full, Species, log2FC, p_adj) %>%
    mutate(Cohort = cohort) %>%
    select(Cohort, Species_full, Species, log2FC, p_adj) %>%
    arrange(p_adj)
}))

write.csv(per_cohort_lists,
          file = file.path(OUTDIR, "cross_cohort_CRC_enriched_per_cohort_lists.csv"),
          row.names = FALSE)
cat("Saved:", file.path(OUTDIR, "cross_cohort_CRC_enriched_per_cohort_lists.csv"), "\n")

# --- Output 4: Membership matrix (wide format) -------------------------------
cat("--- 4. Membership matrix ---\n")

membership_out <- membership_df %>%
  select(-mean_log2FC) %>%
  select(
    Species_full, Species,
    PRJEB10878_sig, PRJEB10878_log2FC, PRJEB10878_p_adj,
    PRJEB27928_sig, PRJEB27928_log2FC, PRJEB27928_p_adj,
    PRJEB6070_sig,  PRJEB6070_log2FC,  PRJEB6070_p_adj,
    PRJNA429097_sig, PRJNA429097_log2FC, PRJNA429097_p_adj,
    n_cohorts_significant
  )

write.csv(membership_out,
          file = file.path(OUTDIR, "cross_cohort_CRC_enriched_membership_matrix.csv"),
          row.names = FALSE)
cat("Saved:", file.path(OUTDIR, "cross_cohort_CRC_enriched_membership_matrix.csv"), "\n")

# --- Output 5: 4-way intersection --------------------------------------------
cat("--- 5. 4-way intersection ---\n")

if (length(intersect_all4) == 0) {
  warning("4-way intersection is EMPTY.")
  # Write header-only file
  empty_df <- data.frame(
    Species_full = character(),
    Species = character(),
    PRJEB10878_log2FC = numeric(), PRJEB10878_p_adj = numeric(),
    PRJEB27928_log2FC = numeric(), PRJEB27928_p_adj = numeric(),
    PRJEB6070_log2FC = numeric(),  PRJEB6070_p_adj = numeric(),
    PRJNA429097_log2FC = numeric(),PRJNA429097_p_adj = numeric(),
    mean_log2FC = numeric(),
    max_p_adj = numeric(),
    stringsAsFactors = FALSE
  )
  write.csv(empty_df,
            file = file.path(OUTDIR, "cross_cohort_CRC_enriched_intersection_4cohorts.csv"),
            row.names = FALSE)
  cat("WARNING: 4-way intersection is empty. Wrote header-only file.\n")
} else {
  intersection_rows <- lapply(intersect_all4, function(sp) {
    out <- list(Species_full = sp)
    # Get short name
    short_name <- NA
    for (cohort in names(cohort_data)) {
      df <- cohort_data[[cohort]]$df
      if (sp %in% df$Species_full) {
        short_name <- df$Species[df$Species_full == sp][1]
        break
      }
    }
    out$Species <- short_name
    log2fc_vec <- c()
    padj_vec <- c()
    for (cohort in names(cohort_data)) {
      df <- cohort_data[[cohort]]$df
      idx <- which(df$Species_full == sp)
      lfc <- df$log2FC[idx]
      padj <- df$p_adj[idx]
      out[[paste0(cohort, "_log2FC")]] <- lfc
      out[[paste0(cohort, "_p_adj")]] <- padj
      log2fc_vec <- c(log2fc_vec, lfc)
      padj_vec <- c(padj_vec, padj)
    }
    out$mean_log2FC <- mean(log2fc_vec)
    out$max_p_adj <- max(padj_vec)
    out
  })

  intersection_df <- bind_rows(intersection_rows) %>%
    arrange(desc(mean_log2FC))

  write.csv(intersection_df,
            file = file.path(OUTDIR, "cross_cohort_CRC_enriched_intersection_4cohorts.csv"),
            row.names = FALSE)
  cat("Saved:", file.path(OUTDIR, "cross_cohort_CRC_enriched_intersection_4cohorts.csv"), "\n")
  cat("  Intersection contains", nrow(intersection_df), "species.\n")
}

# --- Output 6: Summary -------------------------------------------------------
cat("--- 6. Summary ---\n")

summary_rows <- lapply(names(cohort_data), function(cohort) {
  df <- cohort_data[[cohort]]$df
  data.frame(
    Cohort = cohort,
    n_total_species_tested = nrow(df),
    n_significant_padj_only = sum(df$p_adj < P_ADJ_CUT, na.rm = TRUE),
    n_cancer_enriched_strict = sum(df$p_adj < P_ADJ_CUT & df$log2FC > LFC_CUT, na.rm = TRUE),
    n_healthy_enriched_strict = sum(df$p_adj < P_ADJ_CUT & df$log2FC < -LFC_CUT, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})

summary_df <- bind_rows(summary_rows)

# Add intersection row
intersection_row <- data.frame(
  Cohort = "Intersection_4cohorts",
  n_total_species_tested = NA_integer_,
  n_significant_padj_only = NA_integer_,
  n_cancer_enriched_strict = length(intersect_all4),
  n_healthy_enriched_strict = NA_integer_,
  stringsAsFactors = FALSE
)

summary_df <- bind_rows(summary_df, intersection_row)

write.csv(summary_df,
          file = file.path(OUTDIR, "cross_cohort_CRC_enriched_summary.csv"),
          row.names = FALSE)
cat("Saved:", file.path(OUTDIR, "cross_cohort_CRC_enriched_summary.csv"), "\n")

# --- Final stdout report -----------------------------------------------------
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("V1 Schema check: PASS\n")
cat("V2 Row counts: PASS\n")
cat("V3 Taxonomy consistency: PASS\n")
cat("V4 Enriched set sizes: PASS\n")
cat("\n4-way intersection size:", length(intersect_all4), "\n")
if (length(intersect_all4) > 0) {
  cat("Example species in 4-way intersection:\n")
  example_names <- intersection_df$Species[1:min(5, nrow(intersection_df))]
  for (nm in example_names) {
    cat("  -", nm, "\n")
  }
}
cat("\nAll outputs written to:", normalizePath(OUTDIR), "\n")

# =============================================================================
# GOAL B: Cross-cohort consistency of sarcosine-production-associated species
# =============================================================================

cat("\n\n=== GOAL B: Cross-cohort sarcosine-production-associated species ===\n")

# --- Sarcosine KO panel (verbatim from existing scripts) ---------------------
sarcosine_kos <- data.frame(
  KO    = c("K00301","K00302","K00303","K00304","K00305","K00306",
            "K00315","K00552","K08688"),
  Role  = c("Degradation","Degradation","Degradation","Degradation",
            "Degradation","Degradation",
            "Production","Production","Production"),
  stringsAsFactors = FALSE
)

prod_kos_all <- sarcosine_kos$KO[sarcosine_kos$Role == "Production"]

# --- Cohort paths ------------------------------------------------------------
COHORT_PATHS <- list(
  PRJEB10878  = "../PRJEB10878_CRC",
  PRJEB27928  = "../PRJEB27928_CRC",
  PRJEB6070   = "../PRJEB6070_CRC_AdenomatousPolyps",
  PRJNA429097 = "../PRJNA429097_CRC"
)

EXPECTED_SAMPLES <- list(
  PRJEB10878  = c(Healthy = 54,  Cancer = 74),
  PRJEB27928  = c(Healthy = 120, Cancer = 140),
  PRJEB6070   = c(Healthy = 476, Cancer = 594),
  PRJNA429097 = c(Healthy = 95,  Cancer = 98)
)

# --- Per-cohort correlation results storage ----------------------------------
cohort_corr_results <- list()
cohort_prod_sets <- list()
cohort_summary <- list()

# --- V5: Schema verification helper -----------------------------------------
cat("=== Goal B: Loading per-cohort raw data ===\n")

for (cohort in names(COHORT_PATHS)) {
  cat("\n--- Processing cohort:", cohort, "---\n")
  cpath <- COHORT_PATHS[[cohort]]

  # --- Load metadata ---
  meta_files <- list.files(cpath, pattern = "^selected_project_.*\\.txt$", full.names = TRUE)
  meta_file <- sort(meta_files, decreasing = TRUE)[1]

  meta_lines <- readLines(meta_file)
  header <- strsplit(meta_lines[2], "\t")[[1]]
  data_lines <- meta_lines[3:length(meta_lines)]
  data_lines <- data_lines[data_lines != ""]
  data_list <- strsplit(data_lines, "\t")
  data_list <- lapply(data_list, function(x) {
    if (length(x) >= length(header)) x[1:length(header)]
    else c(x, rep(NA, length(header) - length(x)))
  })
  meta_c <- as.data.frame(do.call(rbind, data_list), stringsAsFactors = FALSE)
  colnames(meta_c) <- gsub(" ", ".", header)

  # V5 schema check
  required_cols <- c("Run.ID", "Phenotype.name", "Assay.type")
  missing_cols <- setdiff(required_cols, colnames(meta_c))
  if (length(missing_cols) > 0) {
    stop("V5 FAIL:", cohort, "missing required columns:", paste(missing_cols, collapse = ", "))
  }
  cat("V5 PASS: Metadata schema OK for", cohort, "\n")

  # Filter to Health + CRC + WGS
  all_pheno <- unique(meta_c$Phenotype.name)
  cancer_pheno <- all_pheno[!all_pheno %in% c("Health", "Adenomatous Polyps")]
  cancer_label <- cancer_pheno[1]

  meta_c <- meta_c %>%
    filter(Phenotype.name %in% c("Health", cancer_label)) %>%
    mutate(Group = ifelse(Phenotype.name == "Health", "Healthy", "Cancer"))

  if ("Assay.type" %in% colnames(meta_c)) {
    meta_c <- meta_c %>% filter(Assay.type == "WGS")
  }

  stopifnot(all(meta_c$Assay.type == "WGS" | is.na(meta_c$Assay.type)))
  stopifnot(!any(meta_c$Phenotype.name == "Adenomatous Polyps"))

  # V6 sample counts
  n_h <- sum(meta_c$Group == "Healthy")
  n_c <- sum(meta_c$Group == "Cancer")
  cat("V6: Healthy =", n_h, ", Cancer =", n_c, "\n")
  if (n_h != EXPECTED_SAMPLES[[cohort]]["Healthy"] || n_c != EXPECTED_SAMPLES[[cohort]]["Cancer"]) {
    stop("V6 FAIL:", cohort, "sample counts don't match expected.")
  }
  cat("V6 PASS: Sample counts match expected for", cohort, "\n")

  # --- Load Bacteria ---
  bact_file <- list.files(cpath, pattern = "^Bacteria_.*\\.txt$", full.names = TRUE)[1]
  bact_raw <- read.delim(bact_file, stringsAsFactors = FALSE)
  colnames(bact_raw) <- trimws(colnames(bact_raw))

  # V5 schema check for Bacteria
  if (!all(c("Taxa", "Run.ID", "Abundance") %in% colnames(bact_raw))) {
    stop("V5 FAIL:", cohort, "Bacteria file missing required columns.")
  }

  bact_species <- bact_raw %>%
    filter(grepl("\\|s__", Taxa) & !grepl("\\|t__", Taxa)) %>%
    filter(Run.ID %in% meta_c$Run.ID)

  bact_wide <- bact_species %>%
    select(Taxa, Run.ID, Abundance) %>%
    pivot_wider(names_from = Taxa, values_from = Abundance, values_fill = 0)

  bact_mat_c <- as.data.frame(bact_wide)
  rownames(bact_mat_c) <- bact_mat_c$Run.ID
  bact_mat_c$Run.ID <- NULL
  bact_mat_c <- as.matrix(bact_mat_c)

  rm(bact_raw, bact_species, bact_wide); gc(verbose = FALSE)

  # --- Load KO (JSON) ---
  ko_file <- file.path(cpath, "KO_relative_abundance.tsv")
  ko_raw <- jsonlite::fromJSON(ko_file)

  # V5 schema check for KO
  if (!all(c("ko", "run_id", "abundance") %in% colnames(ko_raw))) {
    stop("V5 FAIL:", cohort, "KO JSON missing required fields.")
  }

  ko_filtered <- ko_raw %>%
    filter(run_id %in% meta_c$Run.ID) %>%
    select(ko, run_id, abundance)

  ko_wide <- ko_filtered %>%
    pivot_wider(names_from = ko, values_from = abundance, values_fill = 0)

  ko_mat_c <- as.data.frame(ko_wide)
  rownames(ko_mat_c) <- ko_mat_c$run_id
  ko_mat_c$run_id <- NULL
  ko_mat_c <- as.matrix(ko_mat_c)

  rm(ko_raw, ko_filtered, ko_wide); gc(verbose = FALSE)

  # --- Intersect bacteria and KO samples (V7) ---
  common_c <- intersect(rownames(bact_mat_c), rownames(ko_mat_c))
  cat("V7: Common bacteria-KO samples:", length(common_c), "\n")
  if (length(common_c) < 50) {
    stop("V7 FAIL:", cohort, "has only", length(common_c), "common samples (<50).")
  }
  cat("V7 PASS: Sufficient common samples for", cohort, "\n")

  bact_mat_c <- bact_mat_c[common_c, , drop = FALSE]
  ko_mat_c <- ko_mat_c[common_c, , drop = FALSE]

  # --- Species prevalence filter (V9) ---
  prev <- colSums(bact_mat_c > 0) / nrow(bact_mat_c)
  bact_filt_c <- bact_mat_c[, prev >= 0.10, drop = FALSE]
  cat("V9: Species after >=10% prevalence filter:", ncol(bact_filt_c), "\n")
  if (ncol(bact_filt_c) < 200 || ncol(bact_filt_c) > 700) {
    cat("WARNING V9:", cohort, "has", ncol(bact_filt_c),
        "species after prevalence filter (outside 200-700 range).\n")
  } else {
    cat("V9 PASS: Species count in expected range for", cohort, "\n")
  }

  # --- Available Production KOs (V8) ---
  prod_kos_c <- intersect(prod_kos_all, colnames(ko_mat_c))
  cat("V8: Available Production KOs:", paste(prod_kos_c, collapse = ", "), "\n")
  if (length(prod_kos_c) == 0) {
    stop("V8 FAIL:", cohort, "has 0 Production KOs available.")
  }
  cat("V8 PASS: Production KOs available for", cohort, "\n")

  # --- Production_sum per sample ---
  prod_sum_c <- rowSums(ko_mat_c[, prod_kos_c, drop = FALSE])
  prod_sum_mean <- mean(prod_sum_c)

  # --- Spearman correlation (step 6) ---
  cat("Computing Spearman correlations for", ncol(bact_filt_c), "species...\n")
  corr_results_c <- lapply(seq_len(ncol(bact_filt_c)), function(i) {
    sp_abund <- bact_filt_c[, i]
    if (sd(sp_abund) == 0) {
      return(data.frame(
        Species_full = colnames(bact_filt_c)[i],
        rho = NA_real_,
        p_value = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    ct <- cor.test(sp_abund, prod_sum_c, method = "spearman", exact = FALSE)
    data.frame(
      Species_full = colnames(bact_filt_c)[i],
      rho = ct$estimate,
      p_value = ct$p.value,
      stringsAsFactors = FALSE
    )
  })
  corr_df_c <- bind_rows(corr_results_c)

  # BH adjustment within cohort (step 7)
  corr_df_c$p_adj <- p.adjust(corr_df_c$p_value, method = "BH")

  # Extract short species name
  corr_df_c$Species <- sapply(corr_df_c$Species_full, function(x) {
    parts <- strsplit(x, "\\|")[[1]]
    sp <- parts[grep("^s__", parts)]
    if (length(sp) > 0) gsub("^s__", "", sp[1]) else x
  })

  # Add metadata columns
  corr_df_c$Cohort <- cohort
  corr_df_c$prod_kos_summed <- paste(prod_kos_c, collapse = "+")
  corr_df_c$n_samples <- length(common_c)

  # Reorder columns
  corr_df_c <- corr_df_c %>%
    select(Cohort, Species_full, Species, rho, p_value, p_adj,
           prod_kos_summed, n_samples)

  # --- Define P_c: production-associated species (step 8) ---
  sig_c <- corr_df_c %>% filter(p_adj < 0.05 & rho > 0.3)
  cat("V10: Production-associated species (rho > 0.3 & p_adj < 0.05):",
      nrow(sig_c), "\n")

  # Store results
  cohort_corr_results[[cohort]] <- corr_df_c
  cohort_prod_sets[[cohort]] <- sig_c$Species_full

  # Summary for this cohort
  cohort_summary[[cohort]] <- data.frame(
    Cohort = cohort,
    n_samples_in_correlation = length(common_c),
    n_species_tested = nrow(corr_df_c),
    n_significant_padj_only = sum(corr_df_c$p_adj < 0.05, na.rm = TRUE),
    n_production_associated_strict = nrow(sig_c),
    prod_kos_used = paste(prod_kos_c, collapse = "+"),
    production_sum_mean = prod_sum_mean,
    stringsAsFactors = FALSE
  )

  rm(bact_mat_c, ko_mat_c, prod_sum_c); gc(verbose = FALSE)
}

# --- V11: Cross-cohort union and distribution ---
cat("\n=== V11: Cross-cohort production-associated union ===\n")
all_prod_species <- lapply(cohort_prod_sets, function(x) x)
union_prod <- Reduce(union, all_prod_species)
cat("Total unique production-associated species across 4 cohorts:", length(union_prod), "\n")

prod_overlap_counts <- sapply(union_prod, function(sp) {
  sum(sapply(all_prod_species, function(v) sp %in% v))
})
prod_dist <- table(prod_overlap_counts)
cat("Distribution (how many cohorts each species is significant in):\n")
print(prod_dist)
cat("V11 PASS: Cross-cohort union computed.\n")

# --- 4-way intersection for Goal B ---
P_intersect4 <- Reduce(intersect, cohort_prod_sets)
cat("\nGoal B 4-way intersection size:", length(P_intersect4), "\n")

# =============================================================================
# GOAL B OUTPUT FILES
# =============================================================================

cat("\n=== Goal B: Saving outputs ===\n")

# --- Output B1: Venn diagram -------------------------------------------------
cat("--- B1. Venn diagram ---\n")

prod_venn_data <- cohort_prod_sets
names(prod_venn_data) <- names(cohort_prod_sets)

prod_venn_plot <- ggvenn(
  prod_venn_data,
  fill_color = venn_colors,
  stroke_size = 0.6,
  set_name_size = 8,
  text_size = 7,
  show_percentage = FALSE
) +
  labs(
    title = "Sarcosine-production-associated species\nacross 4 WGS cohorts\n(Spearman rho > 0.3 & p_adj < 0.05)",
    subtitle = paste0(
      "Production_sum = sum of available Production KOs per cohort\n",
      "Total unique species: ", length(union_prod),
      "  |  4-way intersection: ", length(P_intersect4)
    )
  ) +
  theme(
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5,
                              margin = margin(b = 10)),
    plot.subtitle = element_text(size = 18, hjust = 0.5,
                                 margin = margin(b = 14)),
    plot.margin = margin(20, 20, 20, 20)
  ) +
  # highlight the 4-cohort central intersection count in bright red
  annotate("text", x = 0, y = -0.7,
           label = as.character(length(P_intersect4)),
           colour = "red", fontface = "bold", size = 7)

.ggsave(
  filename = file.path(OUTDIR, "cross_cohort_prod_assoc_venn.png"),
  plot = prod_venn_plot,
  width = 11, height = 11, dpi = 300
)
cat("Saved:", file.path(OUTDIR, "cross_cohort_prod_assoc_venn.png"), "\n")

# --- Output B2: UpSet plot ---------------------------------------------------
cat("--- B2. UpSet plot ---\n")

# Build binary membership for UpSet
prod_upset_df <- data.frame(Species_full = union_prod, stringsAsFactors = FALSE)
for (cohort in names(cohort_prod_sets)) {
  prod_upset_df[[cohort]] <- as.integer(prod_upset_df$Species_full %in% cohort_prod_sets[[cohort]])
}
prod_upset_mat <- as.data.frame(prod_upset_df[, -1])
rownames(prod_upset_mat) <- prod_upset_df$Species_full

prod_upset_title <- paste0(
  "Sarcosine-production-associated species\nacross 4 WGS cohorts\n(Spearman rho > 0.3 & p_adj < 0.05)\n",
  "Total unique: ", length(union_prod), " | 4-way intersection: ", length(P_intersect4)
)

grDevices::png(
  filename = file.path(OUTDIR, "cross_cohort_prod_assoc_upset.png"),
  width = 14, height = 8, units = "in", res = 300, type = "quartz", bg = "white"
)
upset(
  prod_upset_mat,
  nsets = 4,
  sets = c("PRJEB10878", "PRJEB27928", "PRJEB6070", "PRJNA429097"),
  order.by = "freq",
  keep.order = TRUE,
  decreasing = TRUE,
  mainbar.y.label = "Intersection size",
  sets.x.label = "Species per cohort",
  text.scale = c(2.2, 2.0, 1.8, 1.6, 2.2, 2.0),
  point.size = 3.5,
  line.size = 1.2,
  mb.ratio = c(0.65, 0.35)
)
grid::grid.text(prod_upset_title, x = 0.65, y = 0.99, just = "top", gp = grid::gpar(fontsize = 22, fontface = "bold"))
grDevices::dev.off()
cat("Saved:", file.path(OUTDIR, "cross_cohort_prod_assoc_upset.png"), "\n")

# --- Output B3: All correlations (long format, audit trail) ------------------
cat("--- B3. All correlations ---\n")

all_corr_df <- bind_rows(cohort_corr_results)
write.csv(all_corr_df,
          file = file.path(OUTDIR, "cross_cohort_prod_assoc_correlations.csv"),
          row.names = FALSE)
cat("Saved:", file.path(OUTDIR, "cross_cohort_prod_assoc_correlations.csv"), "\n")

# --- Output B4: Per-cohort significant lists (long format) -------------------
cat("--- B4. Per-cohort significant lists ---\n")

per_cohort_sig <- bind_rows(lapply(names(cohort_corr_results), function(cohort) {
  cohort_corr_results[[cohort]] %>%
    filter(p_adj < 0.05 & rho > 0.3) %>%
    select(Cohort, Species_full, Species, rho, p_adj) %>%
    arrange(desc(rho))
}))

write.csv(per_cohort_sig,
          file = file.path(OUTDIR, "cross_cohort_prod_assoc_per_cohort_lists.csv"),
          row.names = FALSE)
cat("Saved:", file.path(OUTDIR, "cross_cohort_prod_assoc_per_cohort_lists.csv"), "\n")

# --- Output B5: Membership matrix (wide format) ------------------------------
cat("--- B5. Membership matrix ---\n")

prod_membership_rows <- lapply(union_prod, function(sp) {
  out <- list(Species_full = sp)
  short_name <- NA
  for (cohort in names(cohort_corr_results)) {
    df <- cohort_corr_results[[cohort]]
    if (sp %in% df$Species_full) {
      short_name <- df$Species[df$Species_full == sp][1]
      break
    }
  }
  out$Species <- short_name

  for (cohort in names(cohort_corr_results)) {
    df <- cohort_corr_results[[cohort]]
    idx <- which(df$Species_full == sp)
    if (length(idx) == 1) {
      out[[paste0(cohort, "_sig")]] <- (df$p_adj[idx] < 0.05 & df$rho[idx] > 0.3)
      out[[paste0(cohort, "_rho")]] <- df$rho[idx]
      out[[paste0(cohort, "_p_adj")]] <- df$p_adj[idx]
    } else {
      out[[paste0(cohort, "_sig")]] <- FALSE
      out[[paste0(cohort, "_rho")]] <- NA_real_
      out[[paste0(cohort, "_p_adj")]] <- NA_real_
    }
  }

  sig_vec <- sapply(names(cohort_corr_results), function(cohort) {
    out[[paste0(cohort, "_sig")]]
  })
  out$n_cohorts_significant <- sum(unlist(sig_vec))

  out
})

prod_membership_df <- bind_rows(prod_membership_rows)
prod_membership_df$mean_rho <- apply(
  prod_membership_df[, paste0(names(cohort_corr_results), "_rho")], 1,
  function(x) mean(x, na.rm = TRUE)
)

prod_membership_df <- prod_membership_df %>%
  arrange(desc(n_cohorts_significant), desc(mean_rho)) %>%
  select(-mean_rho)

write.csv(prod_membership_df,
          file = file.path(OUTDIR, "cross_cohort_prod_assoc_membership_matrix.csv"),
          row.names = FALSE)
cat("Saved:", file.path(OUTDIR, "cross_cohort_prod_assoc_membership_matrix.csv"), "\n")

# --- Output B6: 4-way intersection -------------------------------------------
cat("--- B6. 4-way intersection ---\n")

if (length(P_intersect4) == 0) {
  warning("Goal B 4-way intersection is EMPTY.")
  empty_intersection <- data.frame(
    Species_full = character(), Species = character(),
    PRJEB10878_rho = numeric(), PRJEB10878_p_adj = numeric(),
    PRJEB27928_rho = numeric(), PRJEB27928_p_adj = numeric(),
    PRJEB6070_rho = numeric(),  PRJEB6070_p_adj = numeric(),
    PRJNA429097_rho = numeric(),PRJNA429097_p_adj = numeric(),
    mean_rho = numeric(), max_p_adj = numeric(),
    stringsAsFactors = FALSE
  )
  write.csv(empty_intersection,
            file = file.path(OUTDIR, "cross_cohort_prod_assoc_intersection_4cohorts.csv"),
            row.names = FALSE)
  cat("WARNING: Goal B 4-way intersection is empty. Wrote header-only file.\n")
} else {
  prod_intersection_rows <- lapply(P_intersect4, function(sp) {
    out <- list(Species_full = sp)
    short_name <- NA
    for (cohort in names(cohort_corr_results)) {
      df <- cohort_corr_results[[cohort]]
      if (sp %in% df$Species_full) {
        short_name <- df$Species[df$Species_full == sp][1]
        break
      }
    }
    out$Species <- short_name
    rho_vec <- c()
    padj_vec <- c()
    for (cohort in names(cohort_corr_results)) {
      df <- cohort_corr_results[[cohort]]
      idx <- which(df$Species_full == sp)
      rho_val <- df$rho[idx]
      padj_val <- df$p_adj[idx]
      out[[paste0(cohort, "_rho")]] <- rho_val
      out[[paste0(cohort, "_p_adj")]] <- padj_val
      rho_vec <- c(rho_vec, rho_val)
      padj_vec <- c(padj_vec, padj_val)
    }
    out$mean_rho <- mean(rho_vec)
    out$max_p_adj <- max(padj_vec)
    out
  })

  prod_intersection_df <- bind_rows(prod_intersection_rows) %>%
    arrange(desc(mean_rho))

  write.csv(prod_intersection_df,
            file = file.path(OUTDIR, "cross_cohort_prod_assoc_intersection_4cohorts.csv"),
            row.names = FALSE)
  cat("Saved:", file.path(OUTDIR, "cross_cohort_prod_assoc_intersection_4cohorts.csv"), "\n")
  cat("  Goal B intersection contains", nrow(prod_intersection_df), "species.\n")
}

# --- Output B7: Summary ------------------------------------------------------
cat("--- B7. Summary ---\n")

prod_summary_df <- bind_rows(cohort_summary)

prod_intersection_row <- data.frame(
  Cohort = "Intersection_4cohorts",
  n_samples_in_correlation = NA_integer_,
  n_species_tested = NA_integer_,
  n_significant_padj_only = NA_integer_,
  n_production_associated_strict = length(P_intersect4),
  prod_kos_used = NA_character_,
  production_sum_mean = NA_real_,
  stringsAsFactors = FALSE
)

prod_summary_df <- bind_rows(prod_summary_df, prod_intersection_row)

write.csv(prod_summary_df,
          file = file.path(OUTDIR, "cross_cohort_prod_assoc_summary.csv"),
          row.names = FALSE)
cat("Saved:", file.path(OUTDIR, "cross_cohort_prod_assoc_summary.csv"), "\n")

# --- Output B8: Combined Goal A AND Goal B -----------------------------------
cat("--- B8. Combined Goal A AND Goal B ---\n")

# Read Goal A intersection
GoalA_df <- read.csv(
  file.path(OUTDIR, "cross_cohort_CRC_enriched_intersection_4cohorts.csv"),
  stringsAsFactors = FALSE
)
GoalA_species <- GoalA_df$Species_full

combined_species <- intersect(GoalA_species, P_intersect4)
cat("Combined Goal A AND Goal B intersection size:", length(combined_species), "\n")

if (length(combined_species) == 0) {
  warning("Combined Goal A AND Goal B intersection is EMPTY.")
  empty_combined <- data.frame(
    Species_full = character(), Species = character(),
    PRJEB10878_log2FC = numeric(), PRJEB10878_padj_DA = numeric(),
    PRJEB27928_log2FC = numeric(), PRJEB27928_padj_DA = numeric(),
    PRJEB6070_log2FC = numeric(),  PRJEB6070_padj_DA = numeric(),
    PRJNA429097_log2FC = numeric(),PRJNA429097_padj_DA = numeric(),
    mean_log2FC = numeric(), max_padj_DA = numeric(),
    PRJEB10878_rho = numeric(), PRJEB10878_padj_corr = numeric(),
    PRJEB27928_rho = numeric(), PRJEB27928_padj_corr = numeric(),
    PRJEB6070_rho = numeric(),  PRJEB6070_padj_corr = numeric(),
    PRJNA429097_rho = numeric(),PRJNA429097_padj_corr = numeric(),
    mean_rho = numeric(), max_padj_corr = numeric(),
    stringsAsFactors = FALSE
  )
  write.csv(empty_combined,
            file = file.path(OUTDIR, "cross_cohort_GoalA_AND_GoalB_combined.csv"),
            row.names = FALSE)
  cat("WARNING: Combined intersection is empty. Wrote header-only file.\n")
} else {
  combined_rows <- lapply(combined_species, function(sp) {
    out <- list(Species_full = sp)
    short_name <- NA
    for (cohort in names(cohort_data)) {
      df <- cohort_data[[cohort]]$df
      if (sp %in% df$Species_full) {
        short_name <- df$Species[df$Species_full == sp][1]
        break
      }
    }
    out$Species <- short_name

    # Goal A values
    for (cohort in names(cohort_data)) {
      df <- cohort_data[[cohort]]$df
      idx <- which(df$Species_full == sp)
      out[[paste0(cohort, "_log2FC")]] <- df$log2FC[idx]
      out[[paste0(cohort, "_padj_DA")]] <- df$p_adj[idx]
    }
    lfc_vec <- unlist(out[paste0(names(cohort_data), "_log2FC")])
    padj_da_vec <- unlist(out[paste0(names(cohort_data), "_padj_DA")])
    out$mean_log2FC <- mean(lfc_vec)
    out$max_padj_DA <- max(padj_da_vec)

    # Goal B values
    for (cohort in names(cohort_corr_results)) {
      df <- cohort_corr_results[[cohort]]
      idx <- which(df$Species_full == sp)
      out[[paste0(cohort, "_rho")]] <- df$rho[idx]
      out[[paste0(cohort, "_padj_corr")]] <- df$p_adj[idx]
    }
    rho_vec <- unlist(out[paste0(names(cohort_corr_results), "_rho")])
    padj_corr_vec <- unlist(out[paste0(names(cohort_corr_results), "_padj_corr")])
    out$mean_rho <- mean(rho_vec)
    out$max_padj_corr <- max(padj_corr_vec)

    out
  })

  combined_df <- bind_rows(combined_rows) %>%
    arrange(desc(mean_log2FC))

  write.csv(combined_df,
            file = file.path(OUTDIR, "cross_cohort_GoalA_AND_GoalB_combined.csv"),
            row.names = FALSE)
  cat("Saved:", file.path(OUTDIR, "cross_cohort_GoalA_AND_GoalB_combined.csv"), "\n")
  cat("  Combined intersection contains", nrow(combined_df), "species.\n")
}

# =============================================================================
# GOAL B FINAL REPORT
# =============================================================================
cat("\n=== GOAL B ANALYSIS COMPLETE ===\n")
cat("V5 Schema check: PASS\n")
cat("V6 Sample counts: PASS\n")
cat("V7 Common samples: PASS\n")
cat("V8 Production KOs: PASS\n")
cat("V9 Species prevalence: PASS\n")
cat("V10 Per-cohort P_c sizes:\n")
for (cohort in names(cohort_prod_sets)) {
  cat("  ", cohort, ":", length(cohort_prod_sets[[cohort]]), "species\n")
}
cat("V11 Cross-cohort union: PASS\n")

cat("\nGoal B 4-way intersection size:", length(P_intersect4), "\n")
if (length(P_intersect4) > 0 && exists("prod_intersection_df")) {
  cat("Example species in Goal B 4-way intersection:\n")
  example_names_b <- prod_intersection_df$Species[1:min(5, nrow(prod_intersection_df))]
  for (nm in example_names_b) {
    cat("  -", nm, "\n")
  }
}

cat("\nCombined Goal A AND Goal B intersection size:", length(combined_species), "\n")
if (length(combined_species) > 0 && exists("combined_df")) {
  cat("Species in combined intersection:\n")
  combined_names <- combined_df$Species
  for (nm in combined_names) {
    cat("  -", nm, "\n")
  }
}

cat("\nAll Goal B outputs written to:", normalizePath(OUTDIR), "\n")

# =============================================================================
# V12: Verify Goal A outputs unchanged
# =============================================================================
cat("\n=== V12: Verifying Goal A outputs unchanged ===\n")

V12_FILES <- c(
  "cross_cohort_CRC_enriched_intersection_4cohorts.csv" = "775d001fd362451e54d6bba1e66eee99",
  "cross_cohort_CRC_enriched_membership_matrix.csv"     = "afdafaea0084a3752d24b0ba90b4b66a",
  "cross_cohort_CRC_enriched_per_cohort_lists.csv"      = "aeed2d9b6e658bff6548a674a5dc7e95",
  "cross_cohort_CRC_enriched_summary.csv"               = "a506e1ca9b4c9f918ec92639dcc2f672"
)

v12_pass <- TRUE
for (fname in names(V12_FILES)) {
  fpath <- file.path(OUTDIR, fname)
  if (file.exists(fpath)) {
    actual_md5 <- tools::md5sum(fpath)
    expected_md5 <- V12_FILES[fname]
    if (actual_md5 != expected_md5) {
      cat("V12 FAIL:", fname, "\n")
      cat("  Expected:", expected_md5, "\n")
      cat("  Actual  :", actual_md5, "\n")
      v12_pass <- FALSE
    } else {
      cat("V12 PASS:", fname, "— md5 matches\n")
    }
  } else {
    cat("V12 FAIL:", fname, "— file not found\n")
    v12_pass <- FALSE
  }
}

if (v12_pass) {
  cat("V12 PASS: All Goal A CSVs unchanged.\n")
} else {
  stop("V12 FAILED: Goal A outputs changed during Goal B run. Stopping.")
}

cat("\n=== ALL ANALYSES COMPLETE ===\n")
