#!/usr/bin/env Rscript

file_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", file_arg[grep("^--file=", file_arg)][1])
script_file <- normalizePath(file_arg, mustWork = TRUE)
integrated_root <- dirname(dirname(script_file))
analysis_root <- dirname(integrated_root)
pooled_candidates <- list.dirs(analysis_root, recursive = TRUE, full.names = TRUE)
source_dir <- pooled_candidates[basename(pooled_candidates) == "CRC_WGS_4cohort_pooled_26.09.01"]
stopifnot(length(source_dir) == 1)
out_dir <- file.path(dirname(source_dir), "CRC_WGS_species_KO_correlation_26.09.01")

input_file <- file.path(out_dir, "CRC_WGS_species_KO_correlation_exact_sample_input.csv")
ko_input_file <- file.path(out_dir, "CRC_WGS_species_KO_correlation_KO_input.csv")
species_input_file <- file.path(out_dir, "CRC_WGS_species_KO_correlation_species_input.csv")
result_file <- file.path(out_dir, "CRC_WGS_species_KO_correlation_all_results_frozen.csv")
heatmap_file <- file.path(out_dir, "CRC_WGS_species_KO_correlation_filtered_heatmap_cells.csv")
dict_file <- file.path(out_dir, "CRC_WGS_species_KO_definitions_and_heatmap_filter.csv")

checks <- data.frame(Check = character(), Pass = logical(), Detail = character(), stringsAsFactors = FALSE)
record <- function(check, pass, detail = "") {
  checks <<- rbind(checks, data.frame(Check = check, Pass = isTRUE(pass), Detail = as.character(detail), stringsAsFactors = FALSE))
}

required <- c(input_file, ko_input_file, species_input_file, result_file, heatmap_file, dict_file)
record("All required exported CSV files exist", all(file.exists(required)), paste(basename(required[!file.exists(required)]), collapse = "; "))
stopifnot(all(file.exists(required)))

x <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
kx <- read.csv(ko_input_file, check.names = FALSE, stringsAsFactors = FALSE)
sx <- read.csv(species_input_file, check.names = FALSE, stringsAsFactors = FALSE)
frozen <- read.csv(result_file, check.names = FALSE, stringsAsFactors = FALSE)
heat <- read.csv(heatmap_file, check.names = FALSE, stringsAsFactors = FALSE)
dict <- read.csv(dict_file, check.names = FALSE, stringsAsFactors = FALSE)

id_cols <- c("Run_ID", "Cohort", "Analysis_Group")
tested_kos <- unique(frozen$KO)
tested_species <- unique(frozen$Species_full)

record("Exact input has 1,647 samples", nrow(x) == 1647, nrow(x))
record("Run IDs are unique", !anyDuplicated(x$Run_ID), length(unique(x$Run_ID)))
record("KO input has 7 tested KO columns", ncol(kx) - length(id_cols) == 7, paste(setdiff(names(kx), id_cols), collapse = ", "))
record("Species input has 411 tested species columns", ncol(sx) - length(id_cols) == 411, ncol(sx) - length(id_cols))
record("Combined input has expected columns", ncol(x) == 3 + 7 + 411, ncol(x))
record("Split inputs have identical Run order", identical(x$Run_ID, kx$Run_ID) && identical(x$Run_ID, sx$Run_ID))
record("Frozen results form complete KO x species grid", nrow(frozen) == length(tested_kos) * length(tested_species), sprintf("%d = %d x %d", nrow(frozen), length(tested_kos), length(tested_species)))
record("All tested KO columns are present", all(tested_kos %in% names(x)))
record("All tested species columns are present", all(tested_species %in% names(x)))

source_species <- read.csv(file.path(source_dir, "CRC_WGS_4cohort_pooled_species_relative_abundance_matrix.csv"), check.names = FALSE, stringsAsFactors = FALSE)
source_ko <- read.csv(file.path(source_dir, "CRC_WGS_4cohort_pooled_KEGG_KO_relative_abundance_matrix.csv"), check.names = FALSE, stringsAsFactors = FALSE)
record("Run order matches pooled source species matrix", identical(x$Run_ID, source_species$Run_ID))
record("Run order matches pooled source KO matrix", identical(x$Run_ID, source_ko$Run_ID))
ko_equal <- isTRUE(all.equal(as.matrix(x[, tested_kos, drop = FALSE]), as.matrix(source_ko[, tested_kos, drop = FALSE]), tolerance = 0, check.attributes = TRUE))
sp_equal <- isTRUE(all.equal(as.matrix(x[, tested_species, drop = FALSE]), as.matrix(source_species[, tested_species, drop = FALSE]), tolerance = 0, check.attributes = TRUE))
record("All exported KO values exactly match source", ko_equal)
record("All exported species values exactly match source", sp_equal)
record("No missing numeric input values", !anyNA(x[, c(tested_kos, tested_species), drop = FALSE]))

cohort_counts <- table(x$Cohort)
expected_counts <- c(PRJEB10878 = 128L, PRJEB27928 = 260L, PRJEB6070 = 1066L, PRJNA429097 = 193L)
record("Per-cohort sample counts match canonical roster", identical(as.integer(cohort_counts[names(expected_counts)]), as.integer(expected_counts)), paste(names(expected_counts), cohort_counts[names(expected_counts)], sep = "=", collapse = "; "))

# Independent full recomputation of every Spearman test and global BH family.
recalc <- frozen
recalc$rho_recalc <- NA_real_
recalc$p_value_recalc <- NA_real_
for (i in seq_len(nrow(recalc))) {
  ct <- suppressWarnings(cor.test(x[[recalc$KO[i]]], x[[recalc$Species_full[i]]], method = "spearman"))
  recalc$rho_recalc[i] <- unname(ct$estimate)
  recalc$p_value_recalc[i] <- ct$p.value
}
recalc$p_adj_recalc <- p.adjust(recalc$p_value_recalc, method = "BH")
rho_diff <- max(abs(recalc$rho - recalc$rho_recalc), na.rm = TRUE)
p_diff <- max(abs(recalc$p_value - recalc$p_value_recalc), na.rm = TRUE)
padj_diff <- max(abs(recalc$p_adj - recalc$p_adj_recalc), na.rm = TRUE)
record("All 2,877 Spearman rho values reproduce", rho_diff < 1e-12, format(rho_diff, scientific = TRUE))
record("All 2,877 raw p values reproduce", p_diff < 1e-12, format(p_diff, scientific = TRUE))
record("All 2,877 BH-adjusted p values reproduce", padj_diff < 1e-12, format(padj_diff, scientific = TRUE))

pass_kos_expected <- dict$KO[dict$Pass_Prevalence & !is.na(dict$p_value) & dict$p_value < 0.05]
record("Filtered heatmap KO set is K00303 and K08688", identical(pass_kos_expected, c("K00303", "K08688")), paste(pass_kos_expected, collapse = ", "))
record("Filtered heatmap has 30 species x 2 KOs", nrow(heat) == 60, sprintf("%d rows; %d species; %d KOs", nrow(heat), length(unique(heat$Species_full)), length(unique(heat$KO))))
record("Heatmap plotted values follow the stored threshold", all(heat$Heatmap_value == ifelse(heat$Passed_association_filter, heat$rho, 0)))

validation_file <- file.path(out_dir, "CRC_WGS_species_KO_correlation_validation.csv")
write.csv(checks, validation_file, row.names = FALSE, na = "")
cat("\nCRC WGS species-KO correlation export validation\n")
print(checks, row.names = FALSE)
cat("\nPassed:", sum(checks$Pass), "/", nrow(checks), "\n")
if (any(!checks$Pass)) stop("Validation failed: see ", validation_file)
cat("VALIDATION PASSED\n")
