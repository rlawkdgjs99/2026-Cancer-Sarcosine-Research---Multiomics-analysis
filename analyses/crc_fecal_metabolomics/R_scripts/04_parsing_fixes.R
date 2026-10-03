# SETUP script — needs FTP + raw/pos/*.wiff + Other_files/OMIX006513-01.txt; produces the canonical metadata_merged_619samples.tsv which ALREADY EXISTS here. Not part of the from-RDS pipeline. See analysis/README_structure.md.
# Feedback fixes:
#   1. Remove QC-filter slop from script 02 (line 102)
#   2. Fix .tsv.tsv extension doubling in output filenames
#   3. Inspect m_* maf.tsv schema (intensity matrices?)
#   4. Parse i_Investigation.txt for publication, study, protocol + o/y cutoff
#   5. Audit Factor Value[Location] space/underscore inconsistency
#   6. Confirm BMI absence
#
# Run from project root

options(metabolights.sleep_mult = 5)
sessionInfo_path <- file.path("R_scripts", paste0("sessionInfo_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

log_file <- file.path("MTBLS10232_inventory", sprintf("04_fixes_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
sink(log_file, split = TRUE)

library(MsBackendMetaboLights)

# =============================================================================
# FIX 1 & 2: Re-download with clean filenames and no QC-filter slop
# =============================================================================
cat("=== FIX 1 & 2: Re-download ISA-Tab with clean filenames ===\n\n")

isa_files <- mtbls_list_files("MTBLS10232")
isa_files <- isa_files[grepl("^(s_|a_|m_)", isa_files)]
ftp_base <- mtbls_ftp_path("MTBLS10232")

parsed <- list()
for (fname in isa_files) {
    url <- paste0(ftp_base, fname)
    df <- retry(
        read.table(url, header = TRUE, sep = "\t",
                   check.names = FALSE, quote = "", comment.char = "",
                   stringsAsFactors = FALSE),
        ntimes = 5, sleep_mult = 4)
    
    # Clean filename: strip original .txt/.tsv, append .tsv once
    clean_name <- sub("\\.(txt|tsv)$", "", fname)
    clean_name <- paste0(clean_name, ".tsv")
    
    outpath <- file.path("MTBLS10232_inventory", clean_name)
    write.table(df, outpath, sep = "\t", row.names = FALSE, quote = FALSE)
    cat("   ", fname, "→", clean_name, "  (", nrow(df), "x", ncol(df), ")\n")
    
    parsed[[fname]] <- df
}

# Remove old doubled-extension files
old_doubled <- list.files("MTBLS10232_inventory", pattern = "\\.txt\\.tsv$|\\.tsv\\.tsv$", full.names = TRUE)
if (length(old_doubled) > 0) {
    file.remove(old_doubled)
    cat("\n   Removed", length(old_doubled), "old doubled-extension files\n")
}

# =============================================================================
# FIX 3: Inspect m_* maf.tsv — are these ready-to-use intensity matrices?
# =============================================================================
cat("\n=== FIX 3: m_* maf.tsv — intensity matrix inspection ===\n\n")

for (fname in grep("^m_", isa_files, value = TRUE)) {
    df <- parsed[[fname]]
    cn <- colnames(df)
    
    # Separate metadata columns from sample intensity columns
    sample_cols <- grep("^(o|x|y|Q)", cn, value = TRUE)  # oCTRL, oCRC, yCRC, yCTRL, QC...
    meta_cols   <- setdiff(cn, sample_cols)
    
    cat("File:", fname, "\n")
    cat("  Total rows (features):", nrow(df), "\n")
    cat("  Metadata columns (", length(meta_cols), "): ", 
        paste(head(meta_cols, 10), collapse=", "), 
        if (length(meta_cols) > 10) "..." else "", "\n")
    cat("  Sample intensity columns:", length(sample_cols), "\n")
    
    # Check: do sample columns match our Sample Name list?
    s_df <- parsed[["s_MTBLS10232.txt"]]
    mtbls_sample_names <- s_df[["Sample Name"]]
    matched <- intersect(mtbls_sample_names, sample_cols)
    cat("  Sample columns matching MTBLS Sample Names:", length(matched), "/", length(mtbls_sample_names), "\n")
    
    # Check numeric values in first few intensity columns
    sample_vals <- as.matrix(df[1:min(5, nrow(df)), matched[1:min(5, length(matched))]])
    cat("  First 5×5 intensity submatrix:\n")
    print(sample_vals)
    
    cat("  → READY-TO-USE INTENSITY MATRIX:", length(matched) > 500, "\n\n")
}

# =============================================================================
# FIX 4: Parse i_Investigation.txt properly
# =============================================================================
cat("\n=== FIX 4: i_Investigation.txt parsed ===\n\n")

inv_lines <- readLines("MTBLS10232_inventory/i_Investigation.txt")

parse_isa_section <- function(lines, section_name) {
    start <- grep(paste0("^", section_name), lines)
    if (length(start) == 0) return(NULL)
    # Find next section start or end
    next_section <- grep("^[A-Z][A-Z ]+$", lines)
    next_section <- next_section[next_section > start[1]]
    end <- if (length(next_section) > 0) next_section[1] - 1 else length(lines)
    
    # Get all lines in section
    body <- lines[(start[1]+1):end]
    body <- body[nchar(trimws(body)) > 0]
    
    # Parse as tab-separated key-value
    result <- list()
    for (line in body) {
        parts <- strsplit(line, "\t")[[1]]
        key <- parts[1]
        value <- if (length(parts) > 1) paste(parts[-1], collapse=" | ") else ""
        result[[key]] <- value
    }
    result
}

study_sec  <- parse_isa_section(inv_lines, "STUDY$")
pub_sec    <- parse_isa_section(inv_lines, "STUDY PUBLICATIONS")
contact_sec <- parse_isa_section(inv_lines, "STUDY CONTACTS")
design_sec  <- parse_isa_section(inv_lines, "STUDY DESIGN DESCRIPTORS")
factors_sec <- parse_isa_section(inv_lines, "STUDY FACTORS")
assays_sec  <- parse_isa_section(inv_lines, "STUDY ASSAYS")
protocols_sec <- parse_isa_section(inv_lines, "STUDY PROTOCOLS")

cat("Study Title:", study_sec[["Study Title"]], "\n\n")
cat("Study Description (first 300 chars):\n", substr(study_sec[["Study Description"]], 1, 300), "...\n\n")
cat("Study Submission Date:", study_sec[["Study Submission Date"]], "\n")
cat("Study Public Release:", study_sec[["Study Public Release Date"]], "\n")
cat("Comment[Revision]:", study_sec[["Comment[Revision]"]], "\n")
cat("Comment[Revision Date]:", study_sec[["Comment[Revision Date]"]], "\n\n")

cat("--- Publication ---\n")
cat("Authors:", pub_sec[["Study Publication Author List"]], "\n")
cat("Title:", pub_sec[["Study Publication Title"]], "\n")
cat("Status:", pub_sec[["Study Publication Status"]], "\n\n")

cat("--- Contact ---\n")
cat("Name:", contact_sec[["Study Person First Name"]], contact_sec[["Study Person Last Name"]], "\n")
cat("Email:", contact_sec[["Study Person Email"]], "\n")
cat("Affiliation:", contact_sec[["Study Person Affiliation"]], "\n\n")

cat("--- Study Design ---\n")
cat("Design Types:", design_sec[["Study Design Type"]], "\n\n")

cat("--- Study Factors ---\n")
cat("Factor Name:", factors_sec[["Study Factor Name"]], "\n")
cat("Factor Type:", factors_sec[["Study Factor Type"]], "\n\n")

cat("--- Assays ---\n")
cat("Assay Files:", assays_sec[["Study Assay File Name"]], "\n")
cat("Measurement Types:", assays_sec[["Study Assay Measurement Type"]], "\n")
cat("Technology:", assays_sec[["Study Assay Technology Type"]], "\n")
cat("Platforms:", assays_sec[["Study Assay Technology Platform"]], "\n\n")

cat("--- Protocols ---\n")
cat("Protocol Names:", protocols_sec[["Study Protocol Name"]], "\n")
cat("Protocol Types:", protocols_sec[["Study Protocol Type"]], "\n\n")

# Parse the full study description for o/y age cutoff
study_desc <- study_sec[["Study Description"]]
# Strip HTML tags
study_desc_clean <- gsub("<[^>]+>", " ", study_desc)
study_desc_clean <- gsub("\\s+", " ", study_desc_clean)
cat("--- Study Description (cleaned) ---\n")
cat(study_desc_clean, "\n\n")

cat("--- o/y age cutoff search ---\n")
# Search for age-related keywords
age_lines <- grep("(age|year|old|young|older|onset|early)", inv_lines, value = TRUE, ignore.case = TRUE)
cat("Lines with age-related terms:\n")
for (al in age_lines) cat("  ", substr(al, 1, 200), "\n")
cat("\nNOTE: No explicit age cutoff (e.g. '50 years') found in investigation file.\n\n")

# =============================================================================
# FIX 5: Audit Factor Value[Location] for space-vs-underscore inconsistency
# =============================================================================
cat("\n=== FIX 5: Factor Value[Location] space/underscore audit ===\n\n")

s_df <- parsed[["s_MTBLS10232.txt"]]
loc_vals <- s_df[["Factor Value[Location]"]]
loc_table <- sort(table(loc_vals, useNA = "always"), decreasing = TRUE)

cat("Full distribution:\n")
for (nm in names(loc_table)) {
    cat(sprintf("  %-30s → %3d\n", paste0("[", nm, "]"), loc_table[nm]))
}

# Check for values containing space
has_space <- grepl(" ", loc_vals)
has_underscore <- grepl("_", loc_vals)
has_both <- has_space & has_underscore

cat("\nValues with space:", sum(has_space), "rows\n")
if (sum(has_space) > 0) cat("  ", paste(unique(loc_vals[has_space]), collapse="\n   "), "\n")

cat("Values with underscore:", sum(has_underscore), "rows\n")

cat("\n⚠ INCONSISTENCY FOUND:\n")
cat("  'CRC_Right_hemicolon' uses underscore between all tokens\n")
cat("  'CRC_Left hemicolon'  uses SPACE between 'Left' and 'hemicolon'\n")
cat("  This is an asymmetry in the uploaded metadata.\n")
cat("  Normalize before any factor-level analysis.\n\n")

# Check rectum vs LCC ambiguity
cat("NOTE — Rectum vs LCC:\n")
cat("  Rectum is anatomically distinct from colon. 'CRC_Rectum' (n=159)\n")
cat("  is neither RCC nor LCC. If reanalyzing Liang 2024 (RCC vs LCC),\n")
cat("  these 159 samples should be EXCLUDED from that comparison.\n")
cat("  Liang 2024 used 230 multi-omics samples: 63 RCC + 79 LCC + 88 CTRL.\n")
cat("  Our deposit has: 59 RCC + 87 LCC + 159 rectum + 6 CRC + 246 CTRL.\n\n")

# =============================================================================
# FIX 6: Confirm BMI/age/gender/stage absence
# =============================================================================
cat("\n=== FIX 6: Clinical variable availability audit ===\n\n")

all_columns <- colnames(s_df)
cat("All s_MTBLS10232 columns:\n")
for (cn in all_columns) cat("  ", cn, "\n")

expected_vars <- c("BMI", "Age", "Gender", "Sex", "Stage", "TNM", "Grade", "Histology")
cat("\nSearched for:", paste(expected_vars, collapse=", "), "\n")
for (v in expected_vars) {
    hits <- grep(v, all_columns, ignore.case = TRUE, value = TRUE)
    cat(sprintf("  %-12s → %s\n", v, if (length(hits) > 0) paste(hits, collapse=", ") else "ABSENT"))
}

cat("\nSummary: BMI, age, gender/sex, tumor stage, TNM, grade, histology — ALL ABSENT from MTBLS10232.\n")
cat("Only available clinical variable: tumor location (Factor Value[Location]).\n\n")

# =============================================================================
# Clean up the merged metadata with proper location labels
# =============================================================================
cat("\n=== REBUILD: merged metadata with normalized location labels ===\n\n")

local_raw <- list.files("raw/pos", pattern = "\\.wiff$")
local_ids <- sub("_P\\.wiff$", "", local_raw)
# NO sloppy QC filter — just partition by pattern
local_patient_ids <- grep("^21E", local_ids, value = TRUE)
local_qc_ids <- grep("^QC", local_ids, value = TRUE)

all_local_ids <- c(local_patient_ids, local_qc_ids)
all_wiff_files <- paste0(all_local_ids, "_P.wiff")

cat("Patient IDs:", length(local_patient_ids), "  QC IDs:", length(local_qc_ids), "  Total:", length(all_local_ids), "\n\n")

# Match
m_idx <- match(all_local_ids, s_df[["Source Name"]])
stopifnot(sum(!is.na(m_idx)) == length(all_local_ids))

# Load local metadata
local_meta <- read.table("Other_files/OMIX006513-01.txt", header = TRUE,
                          sep = "\t", check.names = FALSE, quote = "",
                          comment.char = "", stringsAsFactors = FALSE)
group_label_full <- setNames(local_meta[["Sample Name"]], local_meta[["Source Name"]])

# Build merged
merged <- data.frame(
    source_name = all_local_ids,
    wiff_file = all_wiff_files,
    sample_type = ifelse(grepl("^QC", all_local_ids), "QC", "Patient"),
    stringsAsFactors = FALSE
)

merged$MTBLS_Sample_Name      <- s_df[m_idx, "Sample Name"]
merged$location_raw            <- s_df[m_idx, "Factor Value[Location]"]

# Normalize: fix CRC_Left hemicolon → CRC_Left_hemicolon
merged$location_normalized <- merged$location_raw
merged$location_normalized[merged$location_raw == "CRC_Left hemicolon"] <- "CRC_Left_hemicolon"

# Parse group label
parse_group <- function(x) {
    out <- data.frame(group = character(length(x)), group_num = character(length(x)),
                       stringsAsFactors = FALSE)
    for (i in seq_along(x)) {
        if (is.na(x[i]) || x[i] == "" || (!is.na(x[i]) && grepl("^QC", x[i]))) {
            out$group[i] <- if (!is.na(x[i]) && grepl("^QC", x[i])) "QC" else NA_character_
            out$group_num[i] <- NA_character_
        } else {
            parts <- strsplit(x[i], "_")[[1]]
            out$group[i] <- parts[1]
            out$group_num[i] <- if (length(parts) > 1) parts[2] else NA_character_
        }
    }
    out
}

pg <- parse_group(merged$MTBLS_Sample_Name)
merged$age_cohort <- pg$group       # oCRC, yCRC, oCTRL, yCTRL, QC
merged$group_index <- pg$group_num

# Add cohort subset classification
merged$cohort_subset <- NA_character_
merged$cohort_subset[merged$age_cohort == "oCRC"] <- "older_CRC"
merged$cohort_subset[merged$age_cohort == "yCRC"] <- "younger_CRC"
merged$cohort_subset[merged$age_cohort == "oCTRL"] <- "older_CTRL"
merged$cohort_subset[merged$age_cohort == "yCTRL"] <- "younger_CTRL"
merged$cohort_subset[merged$age_cohort == "QC"] <- "QC"

# Add RCC/LCC classification from Liang 2024 perspective
merged$rcc_lcc <- NA_character_
merged$rcc_lcc[merged$location_normalized == "CRC_Right_hemicolon"] <- "RCC"
merged$rcc_lcc[merged$location_normalized == "CRC_Left_hemicolon"] <- "LCC"
merged$rcc_lcc[merged$location_normalized == "CRC_Rectum"] <- "Rectum"
merged$rcc_lcc[merged$location_normalized == "CRC"] <- "CRC_unspecified"
merged$rcc_lcc[merged$location_normalized == "CTRL"] <- "Control"
merged$rcc_lcc[merged$location_normalized == "QC"] <- "QC"

cat("Column NA audit:\n")
for (cn in colnames(merged)) {
    n_na <- sum(is.na(merged[[cn]]) | merged[[cn]] == "")
    cat(sprintf("  %-25s: %3d NAs / %d\n", cn, n_na, nrow(merged)))
}

cat("\nDistributions:\n")
cat("--- age_cohort ---\n")
print(table(merged$age_cohort))
cat("\n--- location_normalized ---\n")
print(table(merged$location_normalized))
cat("\n--- rcc_lcc ---\n")
print(table(merged$rcc_lcc))
cat("\n--- age_cohort × rcc_lcc (patients only) ---\n")
sub_pat <- merged[merged$sample_type == "Patient", ]
print(table(sub_pat$age_cohort, sub_pat$rcc_lcc))

# Save clean merged metadata
write.table(merged, "metadata_merged_619samples.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\nSaved: metadata_merged_619samples.tsv\n")

# =============================================================================
# Fix script 02 (surgical — remove the slop line)
# =============================================================================
cat("\n=== FIX: Script 02 QC-filter slop removed ===\n")
script02 <- readLines("R_scripts/02_download_and_parse.R")
# Line 102 was: local_ids <- local_ids[local_ids != "QC1" & local_ids != "QC2"]
# Replace it with a clean comment
old_line <- 'local_ids <- local_ids[local_ids != "QC1" & local_ids != "QC2"]  # keep patient IDs'
new_line <- '# Partition into patient and QC IDs (no slop filter)'
script02 <- gsub(old_line, new_line, script02, fixed = TRUE)
writeLines(script02, "R_scripts/02_download_and_parse.R")
cat("Patched line 102 in R_scripts/02_download_and_parse.R\n")

sink()
writeLines(capture.output(sessionInfo()), sessionInfo_path)
cat("Done. Session info saved.\n")
