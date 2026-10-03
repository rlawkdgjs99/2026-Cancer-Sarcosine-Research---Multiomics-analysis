# Step 1 & 2: Install MsBackendMetaboLights, verify, and inventory MTBLS10232
options(metabolights.sleep_mult = 5)
sessionInfo_path <- file.path("R_scripts", paste0("sessionInfo_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".txt"))

log_file <- file.path("MTBLS10232_inventory", sprintf("01_log_%s.txt", format(Sys.time(), "%Y%m%d_%H%M%S")))
sink(log_file, split = TRUE)

cat("=== Step 1: Install and verify ===\n\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
BiocManager::install("MsBackendMetaboLights", ask = FALSE, update = FALSE)

library(MsBackendMetaboLights)

cat("Package version: ")
cat(as.character(packageVersion("MsBackendMetaboLights")), "\n")

stopifnot(packageVersion("MsBackendMetaboLights") >= "1.4.2")

ftp_url <- mtbls_ftp_path("MTBLS10232")
cat("MTBLS10232 FTP URL: ", ftp_url, "\n\n")

cat("=== Step 2: Inventory MTBLS10232 ===\n\n")

all_files <- mtbls_list_files("MTBLS10232")
cat("Total files in MTBLS10232:", length(all_files), "\n\n")

writeLines(all_files, file.path("MTBLS10232_inventory", "all_files.txt"))
cat("Full listing saved to MTBLS10232_inventory/all_files.txt\n\n")

isa_patterns <- c(i = "^i_", s = "^s_", a = "^a_", m = "^m_")
for (tag in names(isa_patterns)) {
    pat <- isa_patterns[tag]
    matches <- all_files[grepl(pat, all_files)]
    cat(sprintf("ISA-Tab %s-files (%s): %d\n", tag, pat, length(matches)))
    if (length(matches) > 0) {
        for (f in matches) cat("  ", f, "\n")
    }
    cat("\n")
}

cat("=== All non-ISA files ===\n")
non_isa <- all_files[!grepl("^(i_|s_|a_|m_)", all_files)]
for (f in non_isa) cat("  ", f, "\n")

sink()
writeLines(capture.output(sessionInfo()), sessionInfo_path)
cat("Done. Log saved to", log_file, "\n")
cat("SessionInfo saved to", sessionInfo_path, "\n")
