#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
input_dir <- normalizePath(file.path(analysis_dir, "..", ".."), mustWork = TRUE)
.libPaths(c(file.path(analysis_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Rcpp)
})

options(stringsAsFactors = FALSE, warn = 1)
set.seed(260824)

table_dir <- file.path(analysis_dir, "results", "tables")
intermediate_dir <- file.path(analysis_dir, "intermediate")
log_dir <- file.path(analysis_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(intermediate_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

plan_sha <- sub(" .*", "", system2("shasum", c("-a", "256", file.path(analysis_dir, "ANALYSIS_PLAN_FROZEN.md")), stdout = TRUE))
if (!identical(plan_sha, "c26ed1600d4134ab897e9078efd900e70deef786981520453b0c7459c4dff2a4")) stop("Frozen plan hash mismatch")

sc_file <- file.path(input_dir, "GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz")
meta_file <- file.path(table_dir, "00_cell_metadata_aligned.csv")
if (!all(file.exists(c(sc_file, meta_file)))) stop("Missing input/preflight metadata")
sc_sha <- sub(" .*", "", system2("shasum", c("-a", "256", sc_file), stdout = TRUE))
if (!identical(sc_sha, "aba15960fc7ee6a2443511bce5177e4d71b964131b6e98597e7a85d0a213ba36")) stop("Count-matrix checksum mismatch")

Rcpp::sourceCpp(code = '
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector parse_qc_accumulate(
    std::string s, NumericVector library_size, IntegerVector detected_features,
    NumericVector mt_count, NumericVector ribo_count, NumericVector housekeeping_count,
    bool is_mt, bool is_ribo, bool is_housekeeping) {
  const size_t n = s.size();
  long long value = 0;
  bool has_digit = false;
  int field = 0;
  double gene_total = 0.0;
  int gene_detected = 0;
  for (size_t pos = 0; pos <= n; ++pos) {
    const char ch = (pos < n ? s[pos] : "\\t"[0]);
    if (ch >= "0"[0] && ch <= "9"[0]) {
      has_digit = true;
      value = value * 10 + (ch - "0"[0]);
      if (value > 2147483647LL) stop("Count exceeds signed 32-bit range");
    } else if (ch == "\\t"[0] || ch == "\\r"[0]) {
      if (!has_digit) stop("Empty or malformed count field");
      if (field >= library_size.size()) stop("Too many count fields");
      const double v = static_cast<double>(value);
      library_size[field] += v;
      gene_total += v;
      if (value > 0) {
        detected_features[field] += 1;
        gene_detected += 1;
      }
      if (is_mt) mt_count[field] += v;
      if (is_ribo) ribo_count[field] += v;
      if (is_housekeeping) housekeeping_count[field] += v;
      ++field;
      value = 0;
      has_digit = false;
      if (ch == "\\r"[0]) break;
    } else {
      stop("Non-integer character in count matrix");
    }
  }
  return NumericVector::create(gene_total, gene_detected, field);
}

// [[Rcpp::export]]
int parse_fill_sparse(
    std::string s, IntegerVector old_to_new, IntegerVector next_pos,
    IntegerVector i_out, NumericVector x_out, int gene_index) {
  const size_t n = s.size();
  long long value = 0;
  bool has_digit = false;
  int field = 0;
  int inserted = 0;
  for (size_t pos = 0; pos <= n; ++pos) {
    const char ch = (pos < n ? s[pos] : "\\t"[0]);
    if (ch >= "0"[0] && ch <= "9"[0]) {
      has_digit = true;
      value = value * 10 + (ch - "0"[0]);
      if (value > 2147483647LL) stop("Count exceeds signed 32-bit range");
    } else if (ch == "\\t"[0] || ch == "\\r"[0]) {
      if (!has_digit) stop("Empty or malformed count field");
      if (field >= old_to_new.size()) stop("Too many count fields");
      const int new_col = old_to_new[field];
      if (value > 0 && new_col >= 0) {
        const int out_pos = next_pos[new_col];
        if (out_pos < 0 || out_pos >= i_out.size()) stop("Sparse output position out of range");
        i_out[out_pos] = gene_index;
        x_out[out_pos] = static_cast<double>(value);
        next_pos[new_col] += 1;
        inserted += 1;
      }
      ++field;
      value = 0;
      has_digit = false;
      if (ch == "\\r"[0]) break;
    } else {
      stop("Non-integer character in count matrix");
    }
  }
  if (field != old_to_new.size()) stop("Sparse row width mismatch");
  return inserted;
}
')

read_header <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  line <- readLines(con, n = 1L, warn = FALSE)
  fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
  if (fields[1] != "Gene") stop("Matrix first field is not Gene")
  fields[-1]
}

cell_id <- read_header(sc_file)
n_cells <- length(cell_id)
if (n_cells != 92330L || anyDuplicated(cell_id)) stop("Unexpected or duplicated cell IDs")
cell_meta <- fread(meta_file)
if (!identical(cell_meta$cell_id, cell_id)) stop("Preflight metadata/cell order mismatch")

n_expected_genes <- 24292L
gene <- character(n_expected_genes)
gene_total <- numeric(n_expected_genes)
gene_detected <- integer(n_expected_genes)
library_size <- numeric(n_cells)
detected_features <- integer(n_cells)
mt_count <- numeric(n_cells)
ribo_count <- numeric(n_cells)
housekeeping_count <- numeric(n_cells)

message("Pass 1/2: full-matrix QC scan")
con <- gzfile(sc_file, open = "rt")
invisible(readLines(con, n = 1L, warn = FALSE))
g <- 0L
repeat {
  line <- readLines(con, n = 1L, warn = FALSE)
  if (!length(line)) break
  g <- g + 1L
  if (g > n_expected_genes) stop("More gene rows than expected")
  first_tab <- regexpr("\t", line, fixed = TRUE)[1]
  if (first_tab < 1L) stop("Malformed row ", g + 1L)
  this_gene <- substr(line, 1L, first_tab - 1L)
  if (!nzchar(this_gene)) stop("Blank gene identifier")
  gene[g] <- this_gene
  stats <- parse_qc_accumulate(
    substr(line, first_tab + 1L, nchar(line)), library_size, detected_features,
    mt_count, ribo_count, housekeeping_count,
    grepl("^MT-", this_gene), grepl("^RP[SL]", this_gene),
    this_gene %in% c("ACTB", "GAPDH", "MALAT1")
  )
  if (as.integer(stats[3]) != n_cells) stop("Row-width mismatch for ", this_gene)
  gene_total[g] <- stats[1]
  gene_detected[g] <- as.integer(stats[2])
  if (g %% 2000L == 0L) message("QC scan: ", g, " / ", n_expected_genes, " genes")
}
close(con)
if (g != n_expected_genes || anyDuplicated(gene)) stop("Gene-row audit failed")
if (any(library_size <= 0) || any(!is.finite(library_size))) stop("Invalid library sizes")
if (!identical(sum(library_size), sum(gene_total))) stop("Row/column UMI totals disagree")

qc <- copy(cell_meta)
qc[, `:=`(
  library_size = library_size,
  detected_features = detected_features,
  percent_mt = 100 * mt_count / library_size,
  percent_ribo = 100 * ribo_count / library_size,
  housekeeping_umi = housekeeping_count
)]
qc[, source_study_qc_audit := detected_features >= 500 & percent_mt <= 20 & percent_ribo <= 50 & housekeeping_umi >= 1]
qc[, keep_drug_resistance_2024 := percent_mt <= 20]

fwrite(qc, file.path(table_dir, "01_cell_qc_full_matrix.csv"))
fwrite(data.table(gene = gene, total_umi = gene_total, detected_cells = gene_detected), file.path(table_dir, "01_gene_qc_full_matrix.csv"))

qc_summary <- data.table(
  metric = c(
    "input_cells", "cells_mt_le_20", "cells_mt_gt_20", "cells_meeting_source_study_full_qc_audit",
    "published_2024_cell_count", "observed_minus_published_cells",
    "input_genes", "total_umi", "zero_library_cells", "duplicate_cell_ids", "duplicate_gene_ids"
  ),
  value = c(
    n_cells, sum(qc$keep_drug_resistance_2024), sum(!qc$keep_drug_resistance_2024),
    sum(qc$source_study_qc_audit), 92031L, sum(qc$keep_drug_resistance_2024) - 92031L,
    length(gene), sum(library_size), sum(library_size == 0),
    anyDuplicated(cell_id), anyDuplicated(gene)
  )
)
fwrite(qc_summary, file.path(table_dir, "01_qc_summary.csv"))
message("Cells after mt<=20%: ", sum(qc$keep_drug_resistance_2024), " (2024 paper reports 92,031)")
if (!file.exists(file.path(log_dir, "01_cell_count_discrepancy_decision.md"))) stop("Cell-count discrepancy decision is not documented")
if (sum(qc$keep_drug_resistance_2024) != 92053L) stop("Verified mt<=20% cell count changed unexpectedly")
warning("Published total 92,031 is not reproducible from the deposited matrix: retaining the 92,053 cells that exactly meet percent_mt <= 20; see logs/01_cell_count_discrepancy_decision.md")

keep <- qc$keep_drug_resistance_2024
kept_cell_id <- cell_id[keep]
kept_feature_count <- detected_features[keep]
total_nnz <- sum(as.double(kept_feature_count))
if (total_nnz <= 0 || total_nnz > .Machine$integer.max) stop("Sparse nonzero count outside dgCMatrix limit")
message("Pass 2/2: building full sparse matrix with ", format(total_nnz, big.mark = ","), " nonzero entries")

p_double <- c(0, cumsum(as.double(kept_feature_count)))
if (tail(p_double, 1) != total_nnz) stop("Column pointer total mismatch")
p <- as.integer(p_double)
old_to_new <- rep.int(-1L, n_cells)
old_to_new[keep] <- seq_len(sum(keep)) - 1L
next_pos <- p[seq_len(sum(keep))]
i_out <- integer(as.integer(total_nnz))
x_out <- numeric(as.integer(total_nnz))

con <- gzfile(sc_file, open = "rt")
invisible(readLines(con, n = 1L, warn = FALSE))
g <- 0L
inserted_total <- 0
repeat {
  line <- readLines(con, n = 1L, warn = FALSE)
  if (!length(line)) break
  g <- g + 1L
  first_tab <- regexpr("\t", line, fixed = TRUE)[1]
  this_gene <- substr(line, 1L, first_tab - 1L)
  if (!identical(this_gene, gene[g])) stop("Gene order changed between passes")
  inserted_total <- inserted_total + parse_fill_sparse(
    substr(line, first_tab + 1L, nchar(line)), old_to_new, next_pos,
    i_out, x_out, g - 1L
  )
  if (g %% 2000L == 0L) message("Sparse build: ", g, " / ", n_expected_genes, " genes")
}
close(con)
if (g != n_expected_genes || inserted_total != total_nnz) stop("Sparse build total mismatch")
if (!identical(next_pos, p[2:(length(p))])) stop("Sparse column pointers were not filled exactly")

counts <- new(
  "dgCMatrix", i = i_out, p = p, x = x_out,
  Dim = as.integer(c(n_expected_genes, sum(keep))),
  Dimnames = list(gene, kept_cell_id), factors = list()
)
if (!identical(as.numeric(Matrix::colSums(counts)), library_size[keep])) stop("Sparse matrix library-size verification failed")
if (!identical(as.integer(Matrix::colSums(counts > 0)), detected_features[keep])) stop("Sparse matrix feature-count verification failed")
if (any(counts@x < 0) || any(counts@x != floor(counts@x))) stop("Sparse matrix contains negative/non-integer counts")

saveRDS(counts, file.path(intermediate_dir, "01_full_counts_mt20_qc.rds"), compress = FALSE)
fwrite(qc[keep], file.path(table_dir, "01_cell_metadata_mt20_qc.csv"))
capture.output(sessionInfo(), file = file.path(log_dir, "01_sessionInfo.txt"))

cat("QC and sparse build PASS\n")
cat("Cells:", ncol(counts), "\n")
cat("Genes:", nrow(counts), "\n")
cat("Nonzeros:", length(counts@x), "\n")
