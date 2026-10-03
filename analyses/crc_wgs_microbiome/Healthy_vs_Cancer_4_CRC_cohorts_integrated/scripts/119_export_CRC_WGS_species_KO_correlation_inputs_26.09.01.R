#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dplyr))

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", file_arg[grep("^--file=", file_arg)][1])
script_file <- normalizePath(file_arg, mustWork = TRUE)
integrated_root <- dirname(dirname(script_file))
analysis_root <- dirname(integrated_root)
pooled_candidates <- list.dirs(analysis_root, recursive = TRUE, full.names = TRUE)
pooled_dir <- pooled_candidates[basename(pooled_candidates) == "CRC_WGS_4cohort_pooled_26.09.01"]
stopifnot(length(pooled_dir) == 1)
use_data_root <- dirname(pooled_dir)
out_dir <- file.path(use_data_root, "CRC_WGS_species_KO_correlation_26.09.01")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

species_file <- file.path(pooled_dir, "CRC_WGS_4cohort_pooled_species_relative_abundance_matrix.csv")
ko_file <- file.path(pooled_dir, "CRC_WGS_4cohort_pooled_KEGG_KO_relative_abundance_matrix.csv")
cor_file <- file.path(integrated_root, "results_integrated", "sarcosine", "sarcosine_species_correlation_pooled.csv")
ko_stats_file <- file.path(integrated_root, "results_integrated", "sarcosine", "sarcosine_KO_comparison_filtered_pooled.csv")
stopifnot(file.exists(species_file), file.exists(ko_file), file.exists(cor_file), file.exists(ko_stats_file))

message("Reading pooled species matrix ...")
species_all <- read.csv(species_file, check.names = FALSE, stringsAsFactors = FALSE)
message("Reading pooled KO matrix ...")
ko_all <- read.csv(ko_file, check.names = FALSE, stringsAsFactors = FALSE)
cor_frozen <- read.csv(cor_file, check.names = FALSE, stringsAsFactors = FALSE)
ko_stats <- read.csv(ko_stats_file, check.names = FALSE, stringsAsFactors = FALSE)

id_cols <- c("Run_ID", "Cohort", "Analysis_Group")
stopifnot(all(id_cols %in% names(species_all)), all(id_cols %in% names(ko_all)))
stopifnot(identical(species_all$Run_ID, ko_all$Run_ID))
stopifnot(identical(species_all$Cohort, ko_all$Cohort))
stopifnot(identical(species_all$Analysis_Group, ko_all$Analysis_Group))

# Preserve the exact tested feature sets from the frozen correlation grid.
tested_kos <- unique(cor_frozen$KO)
tested_species <- unique(cor_frozen$Species_full)
stopifnot(length(tested_kos) * length(tested_species) == nrow(cor_frozen))
stopifnot(all(tested_kos %in% names(ko_all)))
stopifnot(all(tested_species %in% names(species_all)))

sample_metadata <- species_all[, id_cols, drop = FALSE]
ko_input <- cbind(sample_metadata, ko_all[, tested_kos, drop = FALSE])
species_input <- cbind(sample_metadata, species_all[, tested_species, drop = FALSE])
combined_input <- cbind(sample_metadata, ko_all[, tested_kos, drop = FALSE], species_all[, tested_species, drop = FALSE])

write.csv(combined_input, file.path(out_dir, "CRC_WGS_species_KO_correlation_exact_sample_input.csv"), row.names = FALSE, na = "")
write.csv(ko_input, file.path(out_dir, "CRC_WGS_species_KO_correlation_KO_input.csv"), row.names = FALSE, na = "")
write.csv(species_input, file.path(out_dir, "CRC_WGS_species_KO_correlation_species_input.csv"), row.names = FALSE, na = "")
write.csv(cor_frozen, file.path(out_dir, "CRC_WGS_species_KO_correlation_all_results_frozen.csv"), row.names = FALSE, na = "")

ko_dictionary <- ko_stats %>%
  filter(KO %in% tested_kos) %>%
  mutate(
    Used_in_all_species_KO_tests = TRUE,
    Passed_heatmap_KO_filter = Pass_Prevalence & !is.na(p_value) & p_value < 0.05
  ) %>%
  arrange(match(KO, tested_kos))
write.csv(ko_dictionary, file.path(out_dir, "CRC_WGS_species_KO_definitions_and_heatmap_filter.csv"), row.names = FALSE, na = "")

# Reconstruct the exact filtered publication heatmap selection from frozen data.
pass_kos <- ko_dictionary$KO[ko_dictionary$Passed_heatmap_KO_filter]
sig_frozen <- cor_frozen %>% filter(p_adj < 0.05, abs(rho) > 0.3, KO %in% pass_kos)
species_rank <- sig_frozen %>%
  group_by(Species_full, Species) %>%
  summarise(max_abs_rho = max(abs(rho)), .groups = "drop") %>%
  arrange(desc(max_abs_rho)) %>%
  slice_head(n = 30) %>%
  mutate(Species_selection_rank = row_number())
heatmap_cells <- cor_frozen %>%
  filter(KO %in% pass_kos, Species_full %in% species_rank$Species_full) %>%
  left_join(species_rank[, c("Species_full", "Species_selection_rank", "max_abs_rho")], by = "Species_full") %>%
  mutate(
    Passed_association_filter = p_adj < 0.05 & abs(rho) > 0.3,
    Heatmap_value = ifelse(Passed_association_filter, rho, 0),
    Significance_symbol = case_when(
      Passed_association_filter & p_adj < 0.001 ~ "***",
      Passed_association_filter & p_adj < 0.01 ~ "**",
      Passed_association_filter & p_adj < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  arrange(Species_selection_rank, match(KO, pass_kos))
write.csv(heatmap_cells, file.path(out_dir, "CRC_WGS_species_KO_correlation_filtered_heatmap_cells.csv"), row.names = FALSE, na = "")

parameters <- data.frame(
  Parameter = c(
    "Cohorts", "Analysis unit", "Common sample count", "Species prevalence threshold",
    "Tested species count", "Tested sarcosine KO count", "Association test",
    "Multiple-testing correction", "Association display threshold",
    "KO filter for filtered heatmap", "Species selection for filtered heatmap"
  ),
  Value = c(
    paste(unique(sample_metadata$Cohort), collapse = "; "), "WGS Run/sample", nrow(sample_metadata),
    ">=10% in the original pooled species matrix", length(tested_species), length(tested_kos),
    "Spearman correlation of individual KO relative abundance vs species relative abundance",
    "Benjamini-Hochberg across all tested KO-species pairs", "BH p_adj <0.05 and absolute rho >0.3",
    "KO prevalence >=10% and Healthy-vs-Cancer Wilcoxon p<0.05",
    "Top 30 species by maximum absolute rho among associations passing both filters"
  ), stringsAsFactors = FALSE
)
write.csv(parameters, file.path(out_dir, "CRC_WGS_species_KO_correlation_analysis_parameters.csv"), row.names = FALSE, na = "")

readme <- c(
  "CRC WGS species-KO correlation data package (2026-09-01)", "",
  "Exact sample-level inputs and frozen results for the pooled 4-cohort CRC WGS analysis",
  "of individual sarcosine-related KO abundance versus species abundance.", "",
  sprintf("Samples: %d", nrow(sample_metadata)),
  sprintf("Tested species: %d", length(tested_species)),
  sprintf("Tested KOs: %d (%s)", length(tested_kos), paste(tested_kos, collapse = ", ")), "",
  "This analysis did NOT correlate species with Production/Degradation sums.",
  "It tested individual KOs using Spearman correlation and BH correction over the full grid.", "",
  "Primary source tables:", paste0("- ", species_file), paste0("- ", ko_file),
  paste0("- ", cor_file), paste0("- ", ko_stats_file), "",
  "CSV guide:",
  "- exact_sample_input: identifiers + 7 KO columns + 411 species columns",
  "- KO_input/species_input: the same input split for inspection",
  "- all_results_frozen: all 2,877 stored KO-species correlations",
  "- filtered_heatmap_cells: exact filtered heatmap cell data",
  "- definitions_and_heatmap_filter: KO definitions and filter status",
  "- analysis_parameters: concise analysis specification"
)
writeLines(readme, file.path(out_dir, "README.txt"), useBytes = TRUE)

message("Export complete: ", out_dir)
message("  samples = ", nrow(sample_metadata))
message("  tested KOs = ", length(tested_kos))
message("  tested species = ", length(tested_species))
message("  frozen correlations = ", nrow(cor_frozen))
message("  filtered heatmap cells = ", nrow(heatmap_cells))
