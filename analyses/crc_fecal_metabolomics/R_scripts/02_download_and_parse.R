# SETUP script — needs MetaboLights FTP + raw/pos/*.wiff + Other_files/OMIX006513-01.txt; rebuilds metadata_merged_619samples.tsv which ALREADY EXISTS here. Not part of the from-RDS pipeline. See analysis/README_structure.md.
# Steps 3-6: Download ISA-Tab metadata, examine, cross-reference, and build merged table
options(metabolights.sleep_mult = 5)
sessionInfo_path <- file.path("R_scripts", paste0("sessionInfo_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

log_file <- file.path("MTBLS10232_inventory", sprintf("02_log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
sink(log_file, split = TRUE)

library(MsBackendMetaboLights)

cat("=== Step 3: Download ISA-Tab metadata files ===\n\n")

isa_files <- mtbls_list_files("MTBLS10232")
isa_files <- isa_files[grepl("^(i_|s_|a_|m_)", isa_files)]
ftp_base <- mtbls_ftp_path("MTBLS10232")

parsed <- list()

for (fname in isa_files) {
    cat("Downloading:", fname, "\n")
    url <- paste0(ftp_base, fname)
    
    if (grepl("^i_", fname)) {
        raw_lines <- retry(readLines(url), ntimes = 5, sleep_mult = 4)
        cat("  Lines:", length(raw_lines), "\n")
        cat("  First 10 lines:\n")
        for (li in seq_len(min(10, length(raw_lines)))) {
            cat("    ", raw_lines[li], "\n")
        }
        outpath <- file.path("MTBLS10232_inventory", paste0(fname, ".txt"))
        writeLines(raw_lines, outpath)
        cat("  Raw file saved to:", outpath, "\n\n")
        parsed[[fname]] <- raw_lines
    } else {
        df <- retry(
            read.table(url, header = TRUE, sep = "\t",
                       check.names = FALSE, quote = "", comment.char = "",
                       stringsAsFactors = FALSE),
            ntimes = 5, sleep_mult = 4)
        cat("  Dimensions:", nrow(df), "x", ncol(df), "\n")
        
        outpath <- file.path("MTBLS10232_inventory", paste0(fname, ".tsv"))
        write.table(df, outpath, sep = "\t", row.names = FALSE, quote = FALSE)
        cat("  Saved to:", outpath, "\n\n")
        
        parsed[[fname]] <- df
    }
}

cat("=== Step 4: Examine sample file (s_MTBLS10232.txt) ===\n\n")

s_name <- grep("^s_", isa_files, value = TRUE)
s_df <- parsed[[s_name]]

cat("Columns (", ncol(s_df), "):\n")
for (j in seq_len(ncol(s_df))) {
    cat(sprintf("  [%d] %s\n", j, colnames(s_df)[j]))
}
cat("\n")

cat("head(5):\n")
print(head(s_df, 5))
cat("\n")

cat("summary:\n")
print(summary(s_df))
cat("\n")

# Identify sample ID column and factor/characteristics columns
cat("--- Column content analysis ---\n")
for (j in seq_len(ncol(s_df))) {
    cn <- colnames(s_df)[j]
    vals <- s_df[[j]]
    n_na <- sum(is.na(vals) | vals == "")
    n_unique <- length(unique(vals))
    cat(sprintf("[%d] %s: %d NAs, %d unique values\n", j, cn, n_na, n_unique))
}
cat("\n")

# Determine which column is the sample identifier
mtbls_ids <- NULL
id_col <- NULL
for (j in seq_len(ncol(s_df))) {
    vals <- s_df[[j]]
    if (length(unique(vals)) == nrow(s_df) && any(grepl("^21E", vals))) {
        mtbls_ids <- vals
        id_col <- colnames(s_df)[j]
        cat("Sample ID column identified:", id_col, "\n")
        break
    }
}
if (is.null(mtbls_ids)) {
    cat("No column with 21E-prefix IDs found. Trying all columns...\n")
    for (j in seq_len(ncol(s_df))) {
        cat("  Column", colnames(s_df)[j], "sample values:", head(unique(s_df[[j]]), 10), "\n")
    }
}

cat("\n=== Cross-reference with local IDs ===\n\n")

local_raw <- list.files("raw/pos", pattern = "\\.wiff$")
local_ids <- sub("_P\\.wiff$", "", local_raw)
# Partition into patient and QC IDs (no slop filter)
local_ids_patients <- grep("^21E", local_ids, value = TRUE)
local_ids_qc <- grep("^QC", local_ids, value = TRUE)

cat("Local .wiff files:", length(local_raw), "\n")
cat("Patient IDs (21E):", length(local_ids_patients), "\n")
cat("QC IDs:", length(local_ids_qc), "\n")

if (!is.null(mtbls_ids)) {
    intersect_n <- length(intersect(mtbls_ids, local_ids_patients))
    only_mtbls <- setdiff(mtbls_ids, local_ids_patients)
    only_local <- setdiff(local_ids_patients, mtbls_ids)
    
    cat("Intersect (matching IDs):", intersect_n, "\n")
    cat("Only in MTBLS10232:", length(only_mtbls), "\n")
    if (length(only_mtbls) <= 20) cat("  ", paste(only_mtbls, collapse=", "), "\n")
    cat("Only in local:", length(only_local), "\n")
    if (length(only_local) <= 20) cat("  ", paste(only_local, collapse=", "), "\n")
}

cat("\n=== Step 5: Examine assay files ===\n\n")

for (a_name in grep("^a_", isa_files, value = TRUE)) {
    cat("---", a_name, "---\n")
    a_df <- parsed[[a_name]]
    cat("Columns (", ncol(a_df), "):\n")
    for (j in seq_len(ncol(a_df))) {
        cn <- colnames(a_df)[j]
        vals <- a_df[[j]]
        n_na <- sum(is.na(vals) | vals == "")
        n_unique <- length(unique(vals))
        cat(sprintf("  [%d] %s: %d NAs, %d unique vals\n", j, cn, n_na, n_unique))
        if (n_unique <= 10) {
            cat("      Values:", paste(unique(vals), collapse=" | "), "\n")
        }
    }
    cat("head(3):\n")
    print(head(a_df, 3))
    cat("\n")
}

cat("=== Step 6: Build merged per-sample metadata table ===\n\n")

# Load local metadata (Other_files/OMIX006513-01.txt) for group labels
local_meta <- read.table("Other_files/OMIX006513-01.txt", header = TRUE,
                          sep = "\t", check.names = FALSE, quote = "",
                          comment.char = "", stringsAsFactors = FALSE)
cat("Local metadata columns:", paste(colnames(local_meta), collapse=", "), "\n")
cat("Local metadata rows:", nrow(local_meta), "\n")

# Build mapping from Source Name to Sample Name (group label)
group_map <- setNames(local_meta[["Sample Name"]], local_meta[["Source Name"]])

# Identify factor/characteristics columns in s_df
factor_cols <- grep("^Factor Value|^Characteristics\\[", colnames(s_df), value = TRUE)

merged <- data.frame(
    source_name = local_ids_patients,
    wiff_file = paste0(local_ids_patients, "_P.wiff"),
    stringsAsFactors = FALSE
)

# Match against s_df (using Source Name column)
if (!is.null(mtbls_ids) && id_col == "Source Name") {
    m_idx <- match(merged$source_name, mtbls_ids)
    for (fc in factor_cols) {
        merged[[fc]] <- s_df[m_idx, fc]
    }
    # Also add Sample Name if present
    if ("Sample Name" %in% colnames(s_df)) {
        merged[["MTBLS_Sample_Name"]] <- s_df[m_idx, "Sample Name"]
    }
}

# Add group label from local metadata
merged$group_label <- group_map[merged$source_name]

cat("Merged table dimensions:", nrow(merged), "x", ncol(merged), "\n")
cat("Columns:", paste(colnames(merged), collapse=", "), "\n\n")

cat("Per-column NA counts:\n")
for (cn in colnames(merged)) {
    cat(sprintf("  %s: %d NAs\n", cn, sum(is.na(merged[[cn]]))))
}
cat("\n")

cat("Unique-value distributions:\n")
for (cn in colnames(merged)) {
    if (is.character(merged[[cn]]) || is.factor(merged[[cn]])) {
        uv <- unique(merged[[cn]])
        if (length(uv) <= 30) {
            cat(sprintf("  %s: %d unique values\n", cn, length(uv)))
            for (v in sort(uv)) {
                cat(sprintf("    %s: %d\n", v, sum(merged[[cn]] == v, na.rm = TRUE)))
            }
        } else {
            cat(sprintf("  %s: %d unique values (too many to list)\n", cn, length(uv)))
        }
    }
}
cat("\n")

# Save merged table
write.table(merged, "metadata_merged_619samples.tsv", sep = "\t", 
            row.names = FALSE, quote = FALSE)
cat("Saved metadata_merged_619samples.tsv\n")

cat("\n=== Investigation file summary ===\n")
i_name <- grep("^i_", isa_files, value = TRUE)
i_content <- parsed[[i_name]]
cat("Lines:", length(i_content), "\n")
for (li in seq_along(i_content)) {
    if (nchar(i_content[li]) > 0) cat("  ", i_content[li], "\n")
}

sink()
writeLines(capture.output(sessionInfo()), sessionInfo_path)
cat("Done. Log saved to", log_file, "\n")
