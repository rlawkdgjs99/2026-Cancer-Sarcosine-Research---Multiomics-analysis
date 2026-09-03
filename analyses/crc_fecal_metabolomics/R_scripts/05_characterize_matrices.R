# Step 1: Load & characterize MTBLS10232 intensity matrices
# No FTP calls — all files already local

sessionInfo_path <- file.path("R_scripts", paste0("sessionInfo_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

log_file <- file.path("analysis", sprintf("01_characterize_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create("analysis", showWarnings = FALSE)
sink(log_file, split = TRUE)

cat("=== Loading intensity matrices ===\n\n")

pos_file <- "MTBLS10232_inventory/m_MTBLS10232_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv"
neg_file <- "MTBLS10232_inventory/m_MTBLS10232_LC-MS_negative_reverse-phase_metabolite_profiling_v2_maf.tsv"

pos_raw <- read.table(pos_file, header = TRUE, sep = "\t",
                       check.names = FALSE, quote = "", comment.char = "",
                       stringsAsFactors = FALSE)
neg_raw <- read.table(neg_file, header = TRUE, sep = "\t",
                       check.names = FALSE, quote = "", comment.char = "",
                       stringsAsFactors = FALSE)

cat("POS raw:", nrow(pos_raw), "features x", ncol(pos_raw), "columns\n")
cat("NEG raw:", nrow(neg_raw), "features x", ncol(neg_raw), "columns\n\n")

# Split metadata vs sample columns
meta_cols <- colnames(pos_raw)[1:23]
sample_cols_pos <- setdiff(colnames(pos_raw), meta_cols)
sample_cols_neg <- setdiff(colnames(neg_raw), meta_cols)

cat("Metadata columns (23):\n")
cat(paste("  ", meta_cols, collapse = "\n"), "\n\n")
cat("Sample columns POS:", length(sample_cols_pos), "\n")
cat("Sample columns NEG:", length(sample_cols_neg), "\n")

# --- Characterize metadata columns ---
cat("\n=== Metabolite identification ===\n")
cat("POS unique metabolite_identification:", length(unique(pos_raw$metabolite_identification)), "\n")
cat("NEG unique metabolite_identification:", length(unique(neg_raw$metabolite_identification)), "\n")

cat("\n=== Database identifier ===\n")
tbl_db_pos <- sort(table(pos_raw$database_identifier), decreasing = TRUE)
tbl_db_neg <- sort(table(neg_raw$database_identifier), decreasing = TRUE)
cat("POS top 10:\n")
for (i in seq_len(min(10, length(tbl_db_pos)))) {
    cat(sprintf("  %s: %d\n", names(tbl_db_pos)[i], tbl_db_pos[i]))
}
cat("NEG top 10:\n")
for (i in seq_len(min(10, length(tbl_db_neg)))) {
    cat(sprintf("  %s: %d\n", names(tbl_db_neg)[i], tbl_db_neg[i]))
}

cat("\n=== Reliability (MetaboLights convention: 1=identified, 2=putatively annotated, 3=putatively characterized, 4=unknown) ===\n")
cat("POS reliability distribution:\n")
print(table(pos_raw$reliability, useNA = "always"))
cat("NEG reliability distribution:\n")
print(table(neg_raw$reliability, useNA = "always"))

# --- Extract intensity matrices ---
cat("\n=== Extracting numeric intensity submatrices ===\n")

intensity_pos <- as.matrix(pos_raw[, sample_cols_pos])
intensity_neg <- as.matrix(neg_raw[, sample_cols_neg])

mode(intensity_pos) <- "numeric"
mode(intensity_neg) <- "numeric"

cat("POS matrix:", nrow(intensity_pos), "x", ncol(intensity_pos), "\n")
cat("NEG matrix:", nrow(intensity_neg), "x", ncol(intensity_neg), "\n")

# --- NA patterns ---
cat("\n=== Missingness per feature ===\n")
na_per_feature_pos <- rowMeans(is.na(intensity_pos))
na_per_feature_neg <- rowMeans(is.na(intensity_neg))
cat(sprintf("POS: %d features (%.1f%%) with >50%% NA\n",
    sum(na_per_feature_pos > 0.5), mean(na_per_feature_pos > 0.5) * 100))
cat(sprintf("NEG: %d features (%.1f%%) with >50%% NA\n",
    sum(na_per_feature_neg > 0.5), mean(na_per_feature_neg > 0.5) * 100))

cat(sprintf("POS: %d features (%.1f%%) with >80%% NA\n",
    sum(na_per_feature_pos > 0.8), mean(na_per_feature_pos > 0.8) * 100))
cat(sprintf("NEG: %d features (%.1f%%) with >80%% NA\n",
    sum(na_per_feature_neg > 0.8), mean(na_per_feature_neg > 0.8) * 100))

cat("\n=== Missingness per sample ===\n")
na_per_sample_pos <- colMeans(is.na(intensity_pos)) * 100
na_per_sample_neg <- colMeans(is.na(intensity_neg)) * 100
cat(sprintf("POS sample NA%%: median=%.1f%%, IQR=%.1f-%.1f%%, max=%.1f%%\n",
    median(na_per_sample_pos), quantile(na_per_sample_pos, 0.25),
    quantile(na_per_sample_pos, 0.75), max(na_per_sample_pos)))
cat(sprintf("NEG sample NA%%: median=%.1f%%, IQR=%.1f-%.1f%%, max=%.1f%%\n",
    median(na_per_sample_neg), quantile(na_per_sample_neg, 0.25),
    quantile(na_per_sample_neg, 0.75), max(na_per_sample_neg)))

# --- Zero values ---
cat("\n=== Zero values ===\n")
zero_fraction_pos <- mean(intensity_pos == 0, na.rm = TRUE)
zero_fraction_neg <- mean(intensity_neg == 0, na.rm = TRUE)
cat(sprintf("POS zero fraction: %.4f%%\n", zero_fraction_pos * 100))
cat(sprintf("NEG zero fraction: %.4f%%\n", zero_fraction_neg * 100))

# --- Intensity distribution stats ---
cat("\n=== Intensity distribution (raw, before transformation) ===\n")
cat("POS:\n")
pos_vals <- intensity_pos[!is.na(intensity_pos) & intensity_pos > 0]
summary_pos <- summary(pos_vals)
cat(sprintf("  Min=%.3f, Q1=%.3f, Median=%.3f, Mean=%.3f, Q3=%.3f, Max=%.3f\n",
    summary_pos["Min."], summary_pos["1st Qu."], summary_pos["Median"],
    summary_pos["Mean"], summary_pos["3rd Qu."], summary_pos["Max."]))
cat(sprintf("  Range span: %.1f orders of magnitude\n",
    log10(summary_pos["Max."]) - log10(summary_pos["Min."])))

cat("NEG:\n")
neg_vals <- intensity_neg[!is.na(intensity_neg) & intensity_neg > 0]
summary_neg <- summary(neg_vals)
cat(sprintf("  Min=%.3f, Q1=%.3f, Median=%.3f, Mean=%.3f, Q3=%.3f, Max=%.3f\n",
    summary_neg["Min."], summary_neg["1st Qu."], summary_neg["Median"],
    summary_neg["Mean"], summary_neg["3rd Qu."], summary_neg["Max."]))
cat(sprintf("  Range span: %.1f orders of magnitude\n",
    log10(summary_neg["Max."]) - log10(summary_neg["Min."])))

# --- Decision: log-transform needed? ---
# Scan the distribution of a few representative samples
cat("\n=== Sample-level intensity histograms (first 6 samples) ===\n")
for (i in 1:min(6, ncol(intensity_pos))) {
    vals <- intensity_pos[, i]
    vals <- vals[!is.na(vals) & vals > 0]
    if (length(vals) > 0) {
        skew <- mean((vals - mean(vals))^3) / (sd(vals)^3)
        cat(sprintf("  %s: median=%.1f, skew=%.2f\n",
            colnames(intensity_pos)[i], median(vals), skew))
    }
}
cat("  → High positive skew in raw intensities → log-transform likely needed.\n")

# --- Save characterization RDS for downstream steps ---
saveRDS(list(
    pos_raw = pos_raw, neg_raw = neg_raw,
    meta_cols = meta_cols,
    sample_cols_pos = sample_cols_pos,
    sample_cols_neg = sample_cols_neg,
    intensity_pos = intensity_pos,
    intensity_neg = intensity_neg,
    feature_meta_pos = pos_raw[, meta_cols, drop = FALSE],
    feature_meta_neg = neg_raw[, meta_cols, drop = FALSE]
), file = "analysis/intensity_matrices.rds")
cat("\nSaved intensity_matrices.rds\n")

sink()
writeLines(capture.output(sessionInfo()), sessionInfo_path)
cat("Done.\n")
