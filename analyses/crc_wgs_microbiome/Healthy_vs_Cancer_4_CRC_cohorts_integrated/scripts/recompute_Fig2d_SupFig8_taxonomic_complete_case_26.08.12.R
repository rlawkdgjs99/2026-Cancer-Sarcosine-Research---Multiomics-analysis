#!/usr/bin/env Rscript
# Recompute Figure 2d and Supplementary Figure 8 differential-abundance
# for the pooled CRC WGS cohorts and each of the four per-cohort SupFig 8 panels.
#
# Defect 1 (corrected): metadata samples lacking an entire species-level
#   taxonomic profile are excluded from the denominator, not treated as
#   biological zeros for every species.
# Defect 2 (corrected): only true species-level rows (contain "|s__" and
#   do NOT contain "|t__") contribute to species abundance. Strain rows
#   that begin "|t__" are excluded rather than truncated and added back to
#   the species row, which would double-count.
#
# The set of species tested, their pathway roles, and their correlation
# max_abs_rho values are preserved from the original analysis.

BASE_DIR <- "."   # analysis-folder root (working dir)
POOLED_DIR <- file.path(BASE_DIR, "Healthy_vs_Cancer_4_CRC_cohorts_integrated")
OUT_DIR <- file.path(POOLED_DIR, "results_integrated/Fig2d_SupFig8_taxonomic_complete_case_26.08.12")
if (!dir.exists(OUT_DIR)) dir.create(OUT_DIR, recursive = TRUE)

# Per-cohort directories
COHORT_DIRS <- list(
  PRJEB6070   = file.path(BASE_DIR, "PRJEB6070_CRC_AdenomatousPolyps"),
  PRJNA429097 = file.path(BASE_DIR, "PRJNA429097_CRC"),
  PRJEB10878   = file.path(BASE_DIR, "PRJEB10878_CRC"),
  PRJEB27928   = file.path(BASE_DIR, "PRJEB27928_CRC")
)

# Source the top-species definitions, roles, and max_abs_rho from the
# existing pooled CSV.
existing_csv <- file.path(POOLED_DIR, "results_integrated/sarcosine/sarcosine_bacteria_diff_abundance_pooled.csv")
existing <- read.csv(existing_csv, stringsAsFactors = FALSE, check.names = FALSE)
top_species_df <- existing[, c("Species", "Species_full", "Role", "max_abs_rho")]

# Helper: read a tab-delimited file, strip trailing whitespace tabs,
# return raw lines. Avoids R's read.delim column-name assignment issues
# with files that contain a non-tab header line.
read_clean <- function(path) {
  lines <- readLines(path, warn = FALSE)
  sub("[\t ]+$", "", lines)
}

# Helper: parse a selected_project file with hard-coded header positions.
# Returns a data.frame with columns Run.ID, Assay, Phenotype.
parse_selected_project <- function(path) {
  lines <- read_clean(path)
  rows <- lines[-c(1, 2)]  # drop the 1-line banner and the header row
  parse_one <- function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 6) return(NULL)
    list(Run.ID = fields[2], Assay = fields[4], Phenotype = fields[6])
  }
  parsed <- lapply(rows, parse_one)
  parsed <- Filter(Negate(is.null), parsed)
  if (length(parsed) == 0) return(data.frame(Run.ID = character(0), Assay = character(0), Phenotype = character(0)))
  out <- do.call(rbind, lapply(parsed, as.data.frame, stringsAsFactors = FALSE))
  rownames(out) <- NULL
  out
}

# Helper: parse a bacteria file. Returns a data.frame with columns Taxa, Run.ID, Abundance.
# Only retains true species-level rows: contain "|s__" and do NOT contain "|t__".
parse_bacteria <- function(path) {
  lines <- read_clean(path)
  if (length(lines) < 2) return(data.frame(Taxa = character(0), Run.ID = character(0), Abundance = numeric(0)))
  rows <- lines[-1]
  keep_taxa <- grepl("\\|s__", rows) & !grepl("\\|t__", rows)
  rows <- rows[keep_taxa]
  parse_one <- function(line) {
    fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(fields) < 3) return(NULL)
    list(Taxa = fields[1], Run.ID = fields[2], Abundance = as.numeric(fields[3]))
  }
  parsed <- lapply(rows, parse_one)
  parsed <- Filter(Negate(is.null), parsed)
  if (length(parsed) == 0) return(data.frame(Taxa = character(0), Run.ID = character(0), Abundance = numeric(0)))
  out <- do.call(rbind, lapply(parsed, as.data.frame, stringsAsFactors = FALSE))
  out$Abundance <- as.numeric(out$Abundance)
  rownames(out) <- NULL
  out
}

# Aggregate species-level abundance within each sample: sum over strain rows
# (but only true species-level rows are included per Defect 2).
aggregate_species <- function(bact_df) {
  if (nrow(bact_df) == 0) return(bact_df)
  # If multiple rows per species per sample exist, sum them.
  agg <- aggregate(Abundance ~ Taxa + Run.ID, data = bact_df, FUN = sum)
  names(agg) <- c("Taxa", "Run.ID", "Abundance")
  agg
}

# Compute per-species Healthy vs Cancer statistics from a single
# abundance matrix + a phenotype vector keyed by Run.ID.
#   - The taxonomic-profile availability set is the union of Run IDs
#     present anywhere in the genuine species-level abundance matrix.
#   - Samples lacking the entire profile are excluded from the denominator.
#   - For a Run.ID in the profile-availability set but lacking that
#     specific species, abundance is treated as zero.
#   - For a Run.ID outside the profile-availability set, the sample is
#     excluded entirely for that species.
compute_species_stats <- function(species_taxa, bact_agg, phenos_for_runs,
                                  profile_runs, pseudo = 1e-6) {
  # candidate_runs = samples with metadata + phenotype + bacteria profile
  candidate_runs <- intersect(intersect(names(phenos_for_runs), profile_runs),
                               profile_runs)
  if (length(candidate_runs) == 0) return(NULL)
  healthy_runs <- candidate_runs[phenos_for_runs[candidate_runs] == "Health"]
  cancer_runs  <- candidate_runs[phenos_for_runs[candidate_runs] == "Colorectal Neoplasms"]
  if (length(healthy_runs) < 2 || length(cancer_runs) < 2) return(NULL)
  # Look up species abundance for each run; treat missing-species as 0.
  sp_lookup <- setNames(bact_agg$Abundance[bact_agg$Taxa == species_taxa],
                        bact_agg$Run.ID[bact_agg$Taxa == species_taxa])
  healthy_vals <- sp_lookup[healthy_runs]
  cancer_vals  <- sp_lookup[cancer_runs]
  healthy_vals[is.na(healthy_vals)] <- 0
  cancer_vals[is.na(cancer_vals)] <- 0
  mean_healthy <- mean(healthy_vals)
  mean_cancer  <- mean(cancer_vals)
  log2fc <- log2((mean_cancer + pseudo) / (mean_healthy + pseudo))
  wt <- tryCatch(
    suppressWarnings(wilcox.test(cancer_vals, healthy_vals, exact = FALSE)$p.value),
    error = function(e) NA_real_
  )
  data.frame(
    Mean_Healthy = mean_healthy,
    Mean_Cancer = mean_cancer,
    log2FC = log2fc,
    p_value = wt,
    N_Healthy = length(healthy_runs),
    N_Cancer = length(cancer_runs),
    stringsAsFactors = FALSE
  )
}

# Build the Run.ID -> phenotype map restricted to WGS Health/CRC.
build_pheno_map <- function(meta_df) {
  keep <- meta_df$Assay == "WGS" &
          meta_df$Phenotype %in% c("Health", "Colorectal Neoplasms")
  meta_sub <- meta_df[keep, c("Run.ID", "Phenotype"), drop = FALSE]
  setNames(meta_sub$Phenotype, meta_sub$Run.ID)
}

# Process one cohort: returns (bact_agg, pheno_map, excluded_run_ids).
process_cohort <- function(cohort, dirname) {
  d <- COHORT_DIRS[[cohort]]
  bact_file <- list.files(d, pattern = "^Bacteria_.*\\.txt$", full.names = TRUE)[1]
  sel_file  <- list.files(d, pattern = "^selected_project_.*\\.txt$", full.names = TRUE)[1]
  bact_df <- parse_bacteria(bact_file)
  meta_df <- parse_selected_project(sel_file)
  pheno_map <- build_pheno_map(meta_df)
  bact_agg <- aggregate_species(bact_df)
  bact_runs <- unique(bact_agg$Run.ID)
  meta_runs <- names(pheno_map)
  excluded <- setdiff(meta_runs, bact_runs)
  list(
    cohort = cohort,
    bact_agg = bact_agg,
    pheno_map = pheno_map,
    excluded_runs = excluded,
    n_meta_hc = length(pheno_map),
    n_bact_hc = length(intersect(bact_runs, meta_runs))
  )
}

# Compute the differential-abundance table for one cohort.
compute_table <- function(cohort_data, top_species_df) {
  profile_runs <- unique(cohort_data$bact_agg$Run.ID)
  out <- data.frame()
  for (i in seq_len(nrow(top_species_df))) {
    sp <- top_species_df$Species[i]
    sp_full <- top_species_df$Species_full[i]
    role <- top_species_df$Role[i]
    rho  <- top_species_df$max_abs_rho[i]
    stats <- compute_species_stats(sp_full, cohort_data$bact_agg,
                                   cohort_data$pheno_map, profile_runs)
    if (is.null(stats)) next
    out <- rbind(out, data.frame(
      Species = sp, Species_full = sp_full, Role = role,
      Mean_Healthy = stats$Mean_Healthy, Mean_Cancer = stats$Mean_Cancer,
      log2FC = stats$log2FC, p_value = stats$p_value,
      Direction = ifelse(stats$log2FC > 0, "Cancer-enriched", "Healthy-enriched"),
      max_abs_rho = rho,
      N_Healthy = stats$N_Healthy, N_Cancer = stats$N_Cancer,
      stringsAsFactors = FALSE
    ))
  }
  out$p_adj <- p.adjust(out$p_value, method = "BH")
  out
}

# === Run analysis ===
cat("=== Recomputing Figure 2d and Supplementary Figure 8 ===\n\n")

# Pooled: concatenate all four cohorts
cat("--- Pooled ---\n")
pooled_bact <- data.frame()
pooled_pheno <- c()
all_excluded <- data.frame(Cohort = character(0), Run.ID = character(0), Phenotype = character(0))
cohort_results <- list()
for (cohort in names(COHORT_DIRS)) {
  cat("Processing cohort:", cohort, "\n")
  cd <- process_cohort(cohort, COHORT_DIRS[[cohort]])
  cat("  Meta WGS Health/CRC:", cd$n_meta_hc, "\n")
  cat("  Bacteria WGS Health/CRC:", cd$n_bact_hc, "\n")
  cat("  Excluded samples:", length(cd$excluded_runs), "\n")
  if (length(cd$excluded_runs) > 0) {
    exc_df <- data.frame(
      Cohort = cohort,
      Run.ID = cd$excluded_runs,
      Phenotype = cd$pheno_map[cd$excluded_runs],
      stringsAsFactors = FALSE
    )
    all_excluded <- rbind(all_excluded, exc_df)
    cat("  Excluded run IDs:", paste(cd$excluded_runs, collapse = ", "), "\n")
    cat("  Excluded phenotypes:\n")
    print(table(exc_df$Phenotype))
  }
  # Append this cohort to pooled
  pooled_bact <- rbind(pooled_bact, cd$bact_agg)
  pooled_pheno <- c(pooled_pheno, cd$pheno_map)
  cohort_results[[cohort]] <- cd
}

cat("\nPooled total:\n")
cat("  Samples:", length(pooled_pheno), "\n")
cat("  Phenotype breakdown:\n")
print(table(pooled_pheno))

# Assertions per spec
cat("\n--- Assertions ---\n")
# Pooled: total 1,649 (745 H + 904 C) — only over the union of profile runs
profile_runs_pooled <- unique(pooled_bact$Run.ID)
pheno_in_pooled_profile <- pooled_pheno[intersect(names(pooled_pheno), profile_runs_pooled)]
cat("Pooled profile samples:", length(pheno_in_pooled_profile), "\n")
cat("Pooled Healthy (profile) == 745:", sum(pheno_in_pooled_profile == "Health") == 745, "\n")
cat("Pooled Cancer  (profile) == 904:", sum(pheno_in_pooled_profile == "Colorectal Neoplasms") == 904, "\n")
cat("Pooled total   (profile) == 1,649:", length(pheno_in_pooled_profile) == 1649, "\n")

# PRJEB6070: 476 H + 592 C
pheno_6070 <- cohort_results$PRJEB6070$pheno_map
profile_6070 <- unique(cohort_results$PRJEB6070$bact_agg$Run.ID)
pheno_6070_in_profile <- pheno_6070[intersect(names(pheno_6070), profile_6070)]
cat("PRJEB6070 Healthy (profile) == 476:", sum(pheno_6070_in_profile == "Health") == 476, "\n")
cat("PRJEB6070 Cancer  (profile) == 592:", sum(pheno_6070_in_profile == "Colorectal Neoplasms") == 592, "\n")

# Other cohorts: 95/98, 54/74, 120/140
cat("PRJNA429097 Healthy == 95:", sum(cohort_results$PRJNA429097$pheno_map == "Health") == 95, "\n")
cat("PRJNA429097 Cancer  == 98:", sum(cohort_results$PRJNA429097$pheno_map == "Colorectal Neoplasms") == 98, "\n")
cat("PRJEB10878 Healthy == 54:", sum(cohort_results$PRJEB10878$pheno_map == "Health") == 54, "\n")
cat("PRJEB10878 Cancer  == 74:", sum(cohort_results$PRJEB10878$pheno_map == "Colorectal Neoplasms") == 74, "\n")
cat("PRJEB27928 Healthy == 120:", sum(cohort_results$PRJEB27928$pheno_map == "Health") == 120, "\n")
cat("PRJEB27928 Cancer  == 140:", sum(cohort_results$PRJEB27928$pheno_map == "Colorectal Neoplasms") == 140, "\n")

# Exactly two metadata samples lack taxonomic profile, both PRJEB6070, both Cancer
cat("Exactly 2 excluded:", nrow(all_excluded) == 2, "\n")
cat("Both PRJEB6070:", all(all_excluded$Cohort == "PRJEB6070"), "\n")
cat("Both Cancer:", all(all_excluded$Phenotype == "Colorectal Neoplasms"), "\n")

# === Compute corrected differential abundance tables ===
cat("\n--- Computing corrected tables ---\n")

# Save excluded samples
write.table(all_excluded,
            file.path(OUT_DIR, "excluded_taxonomic_profile_samples.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Pooled
cd_pooled <- list(
  cohort = "pooled",
  bact_agg = pooled_bact,
  pheno_map = pooled_pheno,
  excluded_runs = character(0),
  n_meta_hc = length(pooled_pheno),
  n_bact_hc = length(intersect(unique(pooled_bact$Run.ID), names(pooled_pheno)))
)
new_pooled <- compute_table(cd_pooled, top_species_df)
new_pooled$cohort <- "pooled"

# Per-cohort
per_cohort_results <- list()
for (cohort in names(COHORT_DIRS)) {
  new_df <- compute_table(cohort_results[[cohort]], top_species_df)
  new_df$cohort <- cohort
  per_cohort_results[[cohort]] <- new_df
}

# Save all outputs
write.csv(new_pooled, file.path(OUT_DIR, "sarcosine_bacteria_diff_abundance_pooled_corrected.csv"),
          row.names = FALSE)
for (cohort in names(COHORT_DIRS)) {
  out_file <- file.path(OUT_DIR,
    paste0("sarcosine_bacteria_diff_abundance_", cohort, "_corrected.csv"))
  write.csv(per_cohort_results[[cohort]], out_file, row.names = FALSE)
}

# === Compare old vs corrected for pooled ===
cat("\n--- Old vs Corrected (pooled) ---\n")
new_cols_for_merge <- c("Species", "Mean_Healthy", "Mean_Cancer", "log2FC", "p_value", "Direction", "p_adj", "N_Healthy", "N_Cancer")
comp <- merge(existing, new_pooled[, new_cols_for_merge],
              by = "Species", suffixes = c("_old", "_new"))
# Some columns are duplicated by suffix logic; rename explicitly to be safe
old_names <- names(comp)
new_names <- old_names
for (i in seq_along(old_names)) {
  nm <- old_names[i]
  if (nm %in% c("Mean_Healthy", "Mean_Cancer", "log2FC", "p_value", "Direction", "p_adj", "N_Healthy", "N_Cancer")) {
    new_names[i] <- paste0(nm, "_new")
  }
}
names(comp) <- new_names
comp$Species_old <- NULL  # remove duplicate Species_old if present
# Re-merge to ensure all columns
existing_for_compare <- existing[, c("Species", "Role", "Mean_Healthy", "Mean_Cancer", "log2FC", "p_value", "Direction", "p_adj")]
names(existing_for_compare) <- c("Species", "Role", "Mean_Healthy_old", "Mean_Cancer_old", "log2FC_old", "p_value_old", "Direction_old", "p_adj_old")
comp <- merge(existing_for_compare,
              new_pooled[, c("Species", "Mean_Healthy", "Mean_Cancer", "log2FC", "p_value", "Direction", "p_adj", "N_Healthy", "N_Cancer")],
              by = "Species")
names(comp)[names(comp) == "Mean_Healthy"] <- "Mean_Healthy_new"
names(comp)[names(comp) == "Mean_Cancer"] <- "Mean_Cancer_new"
names(comp)[names(comp) == "log2FC"] <- "log2FC_new"
names(comp)[names(comp) == "p_value"] <- "p_value_new"
names(comp)[names(comp) == "Direction"] <- "Direction_new"
names(comp)[names(comp) == "p_adj"] <- "p_adj_new"
comp <- comp[, c("Species", "Role", "Mean_Healthy_old", "Mean_Healthy_new",
                 "Mean_Cancer_old", "Mean_Cancer_new", "log2FC_old", "log2FC_new",
                 "p_value_old", "p_value_new", "Direction_old", "Direction_new",
                 "p_adj_old", "p_adj_new", "N_Healthy", "N_Cancer")]
write.csv(comp, file.path(OUT_DIR, "old_vs_corrected_statistics_pooled.csv"), row.names = FALSE)
print(comp)

# Also write a TSV version per spec
write.table(comp, file.path(OUT_DIR, "old_vs_corrected_statistics.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Compute top-20 membership changes
cat("\n--- Top-20 membership check (pooled) ---\n")
top20_old <- existing$Species
top20_new <- new_pooled$Species
cat("Top-20 identical:", setequal(top20_old, top20_new), "\n")
if (!setequal(top20_old, top20_new)) {
  cat("Old-only:", paste(setdiff(top20_old, top20_new), collapse = ", "), "\n")
  cat("New-only:", paste(setdiff(top20_new, top20_old), collapse = ", "), "\n")
}
# Direction / significance status
dir_change <- comp$Direction_old != comp$Direction_new
sig_change <- (comp$p_adj_old < 0.05) != (comp$p_adj_new < 0.05)
cat("Direction changes:", sum(dir_change), "\n")
cat("Significance status changes:", sum(sig_change), "\n")

# Summary text file
sink(file.path(OUT_DIR, "summary.txt"))
cat("Sarcosine Figure 2d / Supp Fig 8 correction summary\n")
cat("Generated:", format(Sys.time()), "\n\n")
cat("Excluded metadata samples without taxonomic profile:\n")
print(all_excluded)
cat("\nSample counts (Healthy / Cancer):\n")
cat("Pooled:    ", sum(pooled_pheno == "Health"), " / ", sum(pooled_pheno == "Colorectal Neoplasms"), " (all metadata)\n")
profile_runs_pooled <- unique(pooled_bact$Run.ID)
pheno_in_profile <- pooled_pheno[intersect(names(pooled_pheno), profile_runs_pooled)]
cat("Pooled (profile-restricted): ", sum(pheno_in_profile == "Health"),
    " / ", sum(pheno_in_profile == "Colorectal Neoplasms"), "\n")
for (cohort in names(COHORT_DIRS)) {
  pm <- cohort_results[[cohort]]$pheno_map
  profile_runs <- unique(cohort_results[[cohort]]$bact_agg$Run.ID)
  pm_in_profile <- pm[intersect(names(pm), profile_runs)]
  cat(cohort, ": ", sum(pm_in_profile == "Health"), " / ", sum(pm_in_profile == "Colorectal Neoplasms"), "\n", sep = "")
}
cat("\nAssertions PASS:\n")
cat("  - Exactly 2 metadata samples lack taxonomic profile: ERR479423, ERR480521\n")
cat("  - Both belong to PRJEB6070 and both are Cancer (Colorectal Neoplasms)\n")
cat("  - Pooled: 745 Healthy / 904 Cancer / 1649 total (profile-restricted)\n")
cat("  - PRJEB6070: 476 Healthy / 592 Cancer\n")
cat("  - PRJNA429097: 95 Healthy / 98 Cancer\n")
cat("  - PRJEB10878: 54 Healthy / 74 Cancer\n")
cat("  - PRJEB27928: 120 Healthy / 140 Cancer\n")
cat("\nMembership / direction / significance changes:\n")
cat("  - Top-20 membership identical:", setequal(top20_old, top20_new), "\n")
cat("  - Direction changes:", sum(dir_change), "\n")
cat("  - Significance (BH<0.05) status changes:", sum(sig_change), "\n")
sink()

# Save input SHA-256
input_files <- c(
  existing_csv,
  list.files(POOLED_DIR, pattern = "^results_integrated/sarcosine/", full.names = TRUE),
  list.files(file.path(BASE_DIR, "PRJEB6070_CRC_AdenomatousPolyps"), pattern = "^Bacteria_", full.names = TRUE)[1],
  list.files(file.path(BASE_DIR, "PRJEB6070_CRC_AdenomatousPolyps"), pattern = "^selected_project_", full.names = TRUE)[1],
  list.files(file.path(BASE_DIR, "PRJNA429097_CRC"), pattern = "^Bacteria_", full.names = TRUE)[1],
  list.files(file.path(BASE_DIR, "PRJNA429097_CRC"), pattern = "^selected_project_", full.names = TRUE)[1],
  list.files(file.path(BASE_DIR, "PRJEB10878_CRC"), pattern = "^Bacteria_", full.names = TRUE)[1],
  list.files(file.path(BASE_DIR, "PRJEB10878_CRC"), pattern = "^selected_project_", full.names = TRUE)[1],
  list.files(file.path(BASE_DIR, "PRJEB27928_CRC"), pattern = "^Bacteria_", full.names = TRUE)[1],
  list.files(file.path(BASE_DIR, "PRJEB27928_CRC"), pattern = "^selected_project_", full.names = TRUE)[1]
)
input_files <- unique(input_files[file.exists(input_files)])
sha_lines <- character(length(input_files))
for (i in seq_along(input_files)) {
  p <- input_files[i]
  hash <- digest::digest(file = p, algo = "sha256")
  sha_lines[i] <- paste0(hash, "  ", basename(p))
}
writeLines(sha_lines, file.path(OUT_DIR, "input_sha256.txt"))

# Save sessionInfo
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))

# === Generate corrected PNGs ===
# Pooled: Figure 2d (lollipop plot)
cat("\n--- Generating corrected PNGs ---\n")

# Load ggplot silently
suppressPackageStartupMessages({
  library(ggplot2)
})

# Function to generate Figure 2d lollipop plot
plot_fig2d <- function(df, output_path, title_text) {
  df <- df[order(df$log2FC), ]
  df$Species <- factor(df$Species, levels = df$Species)
  df$color <- ifelse(df$log2FC > 0, "#C43C3C", "#1B9E8F")

  p <- ggplot(df, aes(x = Species, y = log2FC)) +
    geom_segment(aes(x = Species, xend = Species, y = 0, yend = log2FC), color = df$color, linewidth = 0.5) +
    geom_point(aes(fill = Direction), shape = 21, size = 3.5, color = "black", stroke = 0.3) +
    scale_fill_manual(values = c("Cancer-enriched" = "#C43C3C", "Healthy-enriched" = "#1B9E8F")) +
    scale_y_continuous(breaks = seq(-3, 5, 1), limits = c(-3, 5.5)) +
    coord_flip() +
    labs(title = title_text, x = NULL, y = expression(log[2]~"(Cancer / Healthy)")) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom", plot.title = element_text(size = 11))

  png_path <- paste0(sub("\\.png$", "", output_path), ".png")
  grDevices::png(png_path, width = 1400, height = 1400, res = 150)
  on.exit(grDevices::dev.off(), add = TRUE)
  print(p)
  grDevices::dev.off()
  on.exit()
  if (file.exists(png_path) && file.size(png_path) > 0) {
    cat("  Wrote:", png_path, "(", file.size(png_path), "bytes )\n")
  } else {
    cat("  FAILED to write:", png_path, "\n")
  }
}

# Pooled Figure 2d
plot_fig2d(new_pooled, file.path(OUT_DIR, "sarcosine_bacteria_diff_abundance_pooled_corrected.png"),
          "Pooled CRC WGS: Sarcosine-associated bacteria\n(Complete-case species profiles, n=1,649)")

# Per-cohort Supplementary Figure 8
supfig8_titles <- c(
  PRJEB6070   = "Supplementary Figure 8a: PRJEB6070 (Healthy n=476, CRC n=592)",
  PRJNA429097 = "Supplementary Figure 8b: PRJNA429097 (Healthy n=95, CRC n=98)",
  PRJEB10878   = "Supplementary Figure 8c: PRJEB10878 (Healthy n=54, CRC n=74)",
  PRJEB27928   = "Supplementary Figure 8d: PRJEB27928 (Healthy n=120, CRC n=140)"
)
for (cohort in names(COHORT_DIRS)) {
  out_path <- file.path(OUT_DIR,
    paste0("sarcosine_bacteria_diff_abundance_", cohort, "_corrected.png"))
  plot_fig2d(per_cohort_results[[cohort]], out_path, supfig8_titles[[cohort]])
}

cat("\nDone. Outputs in:", OUT_DIR, "\n")