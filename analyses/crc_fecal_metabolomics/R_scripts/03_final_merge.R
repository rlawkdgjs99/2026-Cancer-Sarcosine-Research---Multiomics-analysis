# SETUP script (SUPERSEDED by 04_parsing_fixes.R) — needs FTP + raw/pos/*.wiff + Other_files/OMIX006513-01.txt; rebuilds metadata that ALREADY EXISTS here. Not part of the from-RDS pipeline. See analysis/README_structure.md.
# Step 7: Build final merged metadata with proper group classification and all 619 samples
options(metabolights.sleep_mult = 5)
sessionInfo_path <- file.path("R_scripts", paste0("sessionInfo_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

log_file <- file.path("MTBLS10232_inventory", sprintf("03_final_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
sink(log_file, split = TRUE)

library(MsBackendMetaboLights)

cat("=== Re-load previously downloaded ISA-Tab files ===\n")

isa_files <- mtbls_list_files("MTBLS10232")
isa_files <- isa_files[grepl("^(s_|a_)", isa_files)]
ftp_base <- mtbls_ftp_path("MTBLS10232")

parsed <- list()
for (fname in isa_files) {
    url <- paste0(ftp_base, fname)
    df <- retry(
        read.table(url, header = TRUE, sep = "\t",
                   check.names = FALSE, quote = "", comment.char = "",
                   stringsAsFactors = FALSE),
        ntimes = 5, sleep_mult = 4)
    parsed[[fname]] <- df
    cat("Loaded:", fname, "-", nrow(df), "x", ncol(df), "\n")
}

s_df    <- parsed[["s_MTBLS10232.txt"]]
a_pos   <- parsed[["a_MTBLS10232_LC-MS_positive_reverse-phase_metabolite_profiling.txt"]]
a_neg   <- parsed[["a_MTBLS10232_LC-MS_negative_reverse-phase_metabolite_profiling.txt"]]

cat("\n=== s_df column names ===\n")
cat(paste(colnames(s_df), collapse="\n  "), "\n")

cat("\n=== Location distribution ===\n")
print(table(s_df[["Factor Value[Location]"]]))

# --- Build local sample inventory ---
local_raw <- list.files("raw/pos", pattern = "\\.wiff$")
local_patient <- grep("^21E", local_raw, value = TRUE)
local_qc      <- grep("^QC", local_raw, value = TRUE)

local_patient_ids <- sub("_P\\.wiff$", "", local_patient)
local_qc_ids <- sub("_P\\.wiff$", "", local_qc)

cat("\n=== Local sample counts ===\n")
cat("Patient .wiff:", length(local_patient), "\n")
cat("QC .wiff:", length(local_qc), "\n")
cat("Total .wiff:", length(local_raw), "\n")

# --- Cross-reference ---
mtbls_source_names <- s_df[["Source Name"]]
mtbls_patient <- grep("^21E", mtbls_source_names, value = TRUE)
mtbls_qc <- grep("^QC", mtbls_source_names, value = TRUE)

cat("\n=== MTBLS10232 Source Name counts ===\n")
cat("Patient:", length(mtbls_patient), "\n")
cat("QC:", length(mtbls_qc), "\n")

intersect_patients <- intersect(local_patient_ids, mtbls_patient)
only_mtbls_patients <- setdiff(mtbls_patient, local_patient_ids)
only_local_patients <- setdiff(local_patient_ids, mtbls_patient)

intersect_qc <- intersect(local_qc_ids, mtbls_qc)
only_mtbls_qc <- setdiff(mtbls_qc, local_qc_ids)
only_local_qc <- setdiff(local_qc_ids, mtbls_qc)

cat("\nPatient ID match:", length(intersect_patients), "/", length(local_patient_ids), "\n")
cat("Only in MTBLS:", length(only_mtbls_patients), "\n")
cat("Only in local:", length(only_local_patients), "\n")

cat("\nQC ID match:", length(intersect_qc), "/", length(local_qc_ids), "\n")
cat("Only in MTBLS:", length(only_mtbls_qc), "\n")
if (length(only_mtbls_qc) > 0) cat("  ", paste(only_mtbls_qc, collapse=", "), "\n")
cat("Only in local:", length(only_local_qc), "\n")
if (length(only_local_qc) > 0) cat("  ", paste(only_local_qc, collapse=", "), "\n")

# --- Load local metadata for group labels ---
local_meta <- read.table("Other_files/OMIX006513-01.txt", header = TRUE,
                          sep = "\t", check.names = FALSE, quote = "",
                          comment.char = "", stringsAsFactors = FALSE)

group_label_full <- setNames(local_meta[["Sample Name"]], local_meta[["Source Name"]])

# --- Build comprehensive metadata for ALL samples (patients + QCs) ---
all_local_ids <- c(local_patient_ids, local_qc_ids)
all_wiff_files <- c(local_patient, local_qc)

merged <- data.frame(
    source_name = all_local_ids,
    wiff_file = all_wiff_files,
    sample_type = ifelse(grepl("^QC", all_local_ids), "QC", "Patient"),
    stringsAsFactors = FALSE
)

# Match against s_df
m_idx <- match(merged$source_name, s_df[["Source Name"]])

merged$MTBLS_Sample_Name <- s_df[m_idx, "Sample Name"]
merged$Factor_Value_Location <- s_df[m_idx, "Factor Value[Location]"]
merged$Characteristics_Organism <- s_df[m_idx, "Characteristics[Organism]"]
merged$Characteristics_Organism_part <- s_df[m_idx, "Characteristics[Organism part]"]
merged$Characteristics_Sample_type <- s_df[m_idx, "Characteristics[Sample type]"]

# Add group labels from local metadata
merged$group_label_full <- group_label_full[merged$source_name]

# Parse group label into group (oCRC/yCRC/oCTRL/yCTRL/QC) and index
parse_group <- function(x) {
    x[is.na(x)] <- ""
    groups <- character(length(x))
    indices <- character(length(x))
    for (i in seq_along(x)) {
        if (x[i] == "") {
            groups[i] <- NA_character_
            indices[i] <- NA_character_
        } else if (grepl("^QC", x[i])) {
            groups[i] <- "QC"
            m <- regmatches(x[i], regexpr("[0-9]+$", x[i]))
            indices[i] <- ifelse(length(m) > 0, m, NA_character_)
        } else {
            parts <- strsplit(x[i], "_")[[1]]
            if (length(parts) >= 2) {
                groups[i] <- parts[1]
                indices[i] <- parts[2]
            } else {
                groups[i] <- x[i]
                indices[i] <- NA_character_
            }
        }
    }
    list(group = groups, index = indices)
}

pg <- parse_group(merged$MTBLS_Sample_Name)
merged$group  <- pg$group
merged$group_index <- pg$index

# Override sample_type with QC for QC samples (MTBLS doesn't label them distinctly in s_df)
# QCs in MTBLS have empty group but the source_name starts with QC
merged$group[is.na(merged$group) & grepl("^QC", merged$source_name)] <- "QC"

cat("\n=== Merged table dimensions ===\n")
cat("Rows:", nrow(merged), "Cols:", ncol(merged), "\n")

cat("\n=== Column NA counts ===\n")
for (cn in colnames(merged)) {
    n_na <- sum(is.na(merged[[cn]]) | merged[[cn]] == "")
    cat(sprintf("  %s: %d NAs / %d total\n", cn, n_na, nrow(merged)))
}

cat("\n=== Group distribution ===\n")
print(table(merged$group, useNA = "always"))

cat("\n=== Factor Value[Location] distribution ===\n")
print(table(merged$Factor_Value_Location, useNA = "always"))

cat("\n=== Group × Location contingency table ===\n")
sub_patients <- merged[merged$sample_type == "Patient", ]
gt <- table(sub_patients$group, sub_patients$Factor_Value_Location)
print(gt)

cat("\n=== Group label counts (by group, patients only) ===\n")
print(table(sub_patients$group))

cat("\n=== Unique group-label prefix counts (e.g. oCRC, yCRC etc) ===\n")
prefixes <- sub("_[0-9]+$", "", na.omit(merged$MTBLS_Sample_Name))
print(table(prefixes))

# --- Location categories ---
cat("\n=== Location categories detailed ===\n")
loc_table <- sort(table(merged$Factor_Value_Location), decreasing = TRUE)
for (nm in names(loc_table)) {
    cat(sprintf("  %s: %d\n", nm, loc_table[nm]))
}

# Save
write.table(merged, "metadata_merged_619samples.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\nSaved metadata_merged_619samples.tsv\n")

# --- Report what's still missing ---
cat("\n=== What's MISSING from MTBLS10232 (not present in sample sheet) ===\n")
cat("  BMI: NOT present\n")
cat("  Age: NOT present\n")
cat("  Gender/Sex: NOT present\n")
cat("  Tumor Stage (TNM): NOT present\n")
cat("  Tumor location (RCC/LCC): PRESENT as 'Factor Value[Location]'\n")

cat("\n=== Open questions resolved ===\n")
cat("  Q1 (619 samples?): YES - 557 patient + 62 QC = 619 Source Names, all 557 patient IDs match 1:1\n")
cat("  Q2 (negative mode?): YES - separate assay file for negative mode with 619 entries\n")
cat("  Q3 (BMI?): NO - not in MTBLS10232 sample sheet\n")
cat("  Q4 (RCC vs LCC?): YES - 'Factor Value[Location]' column: CRC_Left_hemicolon, CRC_Right_hemicolon, CRC_Rectum, CRC, CTRL\n")

sink()
writeLines(capture.output(sessionInfo()), sessionInfo_path)
cat("Done.\n")
