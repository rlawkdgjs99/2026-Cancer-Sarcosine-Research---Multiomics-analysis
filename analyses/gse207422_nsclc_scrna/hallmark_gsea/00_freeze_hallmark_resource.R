#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
.libPaths(c(file.path(analysis_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(msigdbr)
})

sha256 <- function(path) {
  sub(" .*", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
}

plan_file <- file.path(analysis_dir, "ANALYSIS_PLAN_FROZEN.md")
expected_plan_sha <- "5ca4942bb26e7b45f7a0ab4c7f46f8718a34e705ca3a1a581012bae0d40f0cb4"
if (!identical(sha256(plan_file), expected_plan_sha)) stop("Frozen plan SHA-256 mismatch")
if (as.character(packageVersion("msigdbr")) != "26.1.1") stop("Unexpected msigdbr version")

resource_dir <- file.path(analysis_dir, "resources")
log_dir <- file.path(analysis_dir, "logs")
dir.create(resource_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

message("Retrieving MSigDB Human Hallmark collection from msigdbr")
raw <- as.data.table(msigdbr(
  db_species = "HS", species = "Homo sapiens", collection = "H"
))
if (nrow(raw) != 7331L || uniqueN(raw$gs_name) != 50L) {
  stop("Unexpected Hallmark resource dimensions")
}
if (!identical(unique(raw$db_version), "2026.1.Hs")) stop("Unexpected MSigDB release")
if (!identical(unique(raw$gs_collection), "H")) stop("Unexpected MSigDB collection")

resource <- unique(raw[, .(
  gs_name, gene_symbol, gs_description, gs_source_species,
  gs_collection, gs_collection_name, db_version
)])
setorder(resource, gs_name, gene_symbol)
if (nrow(resource) != 7322L || uniqueN(resource$gs_name) != 50L) {
  stop("Unexpected number of unique Hallmark pathway-gene memberships")
}
if (resource[, anyDuplicated(paste(gs_name, gene_symbol, sep = "\r"))]) {
  stop("Duplicate pathway-gene membership")
}

csv_file <- file.path(resource_dir, "MSigDB_Hallmark_2026.1.Hs_gene_symbols.csv")
gmt_file <- file.path(resource_dir, "MSigDB_Hallmark_2026.1.Hs_gene_symbols.gmt")
fwrite(resource, csv_file, quote = TRUE)

descriptions <- resource[, .(description = unique(gs_description)[1]), by = gs_name]
genes <- resource[, .(genes = list(sort(unique(gene_symbol)))), by = gs_name]
gmt <- merge(descriptions, genes, by = "gs_name", sort = TRUE)
gmt_lines <- vapply(seq_len(nrow(gmt)), function(i) {
  paste(c(gmt$gs_name[i], gmt$description[i], gmt$genes[[i]]), collapse = "\t")
}, character(1))
writeLines(gmt_lines, gmt_file, useBytes = TRUE)

manifest <- data.table(
  resource = c("sorted_csv", "gmt"),
  path = c(csv_file, gmt_file),
  sha256 = c(sha256(csv_file), sha256(gmt_file)),
  db_version = "2026.1.Hs",
  msigdbr_version = as.character(packageVersion("msigdbr")),
  pathways = 50L,
  raw_rows = 7331L,
  unique_memberships = 7322L
)
fwrite(manifest, file.path(resource_dir, "00_hallmark_resource_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "00_resource_sessionInfo.txt"))

message("Frozen Hallmark resource complete")
print(manifest)
