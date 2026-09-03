#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
input_dir <- normalizePath(file.path(analysis_dir, "..", ".."), mustWork = TRUE)
local_lib <- file.path(analysis_dir, "R_libs")
.libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(Seurat)
})

options(stringsAsFactors = FALSE, warn = 1)
set.seed(260824)

table_dir <- file.path(analysis_dir, "results", "tables")
log_dir <- file.path(analysis_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) {
  out <- system2("shasum", c("-a", "256", path), stdout = TRUE)
  sub(" .*", "", out)
}

sc_file <- file.path(input_dir, "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz")
meta_file <- file.path(input_dir, "GSE207422_NSCLC_scRNAseq_metadata.xlsx")
plan_file <- file.path(analysis_dir, "ANALYSIS_PLAN_FROZEN.md")
ref_cluster <- file.path(analysis_dir, "references", "original_02_All_cell_clustering.R")
ref_qc <- file.path(analysis_dir, "references", "original_01.2_QC_each_matrix.R")
ref_scrublet <- file.path(analysis_dir, "references", "original_01.3_Scrublet.py")
required <- c(sc_file, meta_file, plan_file, ref_cluster, ref_qc, ref_scrublet)
if (!all(file.exists(required))) stop("Missing required input/reference: ", paste(required[!file.exists(required)], collapse = ", "))

expected <- c(
  scRNA_count_matrix = "aba15960fc7ee6a2443511bce5177e4d71b964131b6e98597e7a85d0a213ba36",
  scRNA_metadata = "d098a750c7ebc595994c929b666177f25c4da6fe3d3e38bb795269a1fb21053e"
)
observed <- c(scRNA_count_matrix = sha256(sc_file), scRNA_metadata = sha256(meta_file))
if (!identical(observed, expected)) stop("Protected input checksum mismatch")

con <- gzfile(sc_file, open = "rt")
on.exit(close(con), add = TRUE)
header <- readLines(con, n = 1L, warn = FALSE)
if (length(header) != 1L) stop("Could not read count-matrix header")
fields <- strsplit(header, "\t", fixed = TRUE)[[1]]
if (fields[1] != "Gene") stop("Unexpected first count-matrix field")
cell_id <- fields[-1]
if (length(cell_id) != 92330L || anyDuplicated(cell_id)) stop("Unexpected or duplicated cell identifiers")
sample <- sub("_[^_]+$", "", cell_id)

meta <- as.data.table(read_excel(meta_file))
required_meta <- c("Sample", "Patient", "Resource", "Pathologic Response")
if (!all(required_meta %in% names(meta))) stop("Missing metadata columns")
idx <- match(sample, meta$Sample)
if (anyNA(idx)) stop("Count-matrix samples not found in metadata")
cell_meta <- meta[idx]
cell_meta[, `:=`(cell_id = cell_id, sample = sample)]
if (!identical(cell_meta$sample, sample)) stop("Cell metadata alignment failure")

input_manifest <- data.table(
  role = c(names(observed), "analysis_plan", "source_code_all_cell", "source_code_qc", "source_code_scrublet"),
  path = c(sc_file, meta_file, plan_file, ref_cluster, ref_qc, ref_scrublet),
  sha256 = vapply(c(sc_file, meta_file, plan_file, ref_cluster, ref_qc, ref_scrublet), sha256, character(1)),
  bytes = file.info(c(sc_file, meta_file, plan_file, ref_cluster, ref_qc, ref_scrublet))$size
)
fwrite(input_manifest, file.path(table_dir, "00_input_manifest.csv"))
fwrite(cell_meta, file.path(table_dir, "00_cell_metadata_aligned.csv"))

summary_table <- data.table(
  metric = c("matrix_cells", "matrix_samples", "metadata_rows", "matched_metadata_cells"),
  value = c(length(cell_id), uniqueN(sample), nrow(meta), nrow(cell_meta))
)
fwrite(summary_table, file.path(table_dir, "00_preflight_summary.csv"))

pkg <- c("R", "Seurat", "SeuratObject", "Matrix", "data.table", "readxl", "future")
ver <- c(R.version.string, vapply(pkg[-1], function(x) as.character(packageVersion(x)), character(1)))
fwrite(data.table(package = pkg, version = ver), file.path(table_dir, "00_software_versions.csv"))
capture.output(sessionInfo(), file = file.path(log_dir, "00_sessionInfo.txt"))

cat("Preflight PASS\n")
cat("Cells:", length(cell_id), "\n")
cat("Samples:", uniqueN(sample), "\n")
cat("Seurat:", as.character(packageVersion("Seurat")), "\n")
