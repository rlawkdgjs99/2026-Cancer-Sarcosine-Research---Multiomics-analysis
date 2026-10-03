#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(digest))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", script_arg))
base_dir <- normalizePath(file.path(dirname(script_path), ".."))
analysis_root <- normalizePath(file.path(base_dir, ".."))
input_path <- file.path(base_dir, "results_integrated", "sarcosine", "sarcosine_KO_per_sample_pooled.csv")
dest_dir <- file.path(analysis_root, "사용데이터_모음", "CRC_WGS_7KO_4cohort_forest_source_data_26.09.01")

input <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
manifest <- read.csv(file.path(dest_dir, "CRC_WGS_7KO_4cohort_forest_SOURCE_PROVENANCE.csv"), check.names = FALSE)
stopifnot(nrow(manifest) == 7L, nrow(input) == 1647L, !anyDuplicated(input$Run.ID))

rank_biserial_direct <- function(healthy, cancer) {
  mean(outer(healthy, cancer, FUN = function(x, y) sign(x - y)))
}

verified_rows <- 0L
for (i in seq_len(nrow(manifest))) {
  ko <- manifest$KO[i]
  plot_path <- file.path(dest_dir, manifest$Plot_data_file[i])
  raw_path <- file.path(dest_dir, manifest$Raw_values_file[i])
  plot <- read.csv(plot_path, check.names = FALSE, stringsAsFactors = FALSE)
  raw <- read.csv(raw_path, check.names = FALSE, stringsAsFactors = FALSE)

  stopifnot(nrow(plot) == 4L, nrow(raw) == 1647L)
  stopifnot(identical(as.character(raw$Run.ID), as.character(input$Run.ID)))
  stopifnot(identical(as.character(raw$Cohort), as.character(input$Cohort)))
  stopifnot(identical(as.character(raw$Group), as.character(input$Group)))
  stopifnot(all(raw$KO == ko))
  stopifnot(isTRUE(all.equal(raw$KO_relative_abundance, input[[ko]], tolerance = 0)))
  stopifnot(identical(
    digest(plot_path, algo = "sha256", file = TRUE),
    digest(manifest$Plot_data_source[i], algo = "sha256", file = TRUE)
  ))

  for (j in seq_len(nrow(plot))) {
    z <- raw[raw$Cohort == plot$Cohort[j], , drop = FALSE]
    healthy <- z$KO_relative_abundance[z$Group == "Healthy"]
    cancer <- z$KO_relative_abundance[z$Group == "Cancer"]
    stopifnot(length(healthy) == plot$n_Healthy[j], length(cancer) == plot$n_CRC[j])
    stopifnot(sum(healthy == 0) == plot$zero_Healthy[j], sum(cancer == 0) == plot$zero_CRC[j])
    if (all(c(healthy, cancer) == 0)) {
      stopifnot(abs(plot$rank_biserial_Healthy_vs_CRC[j]) < 1e-15)
      stopifnot(is.na(plot$wilcox_p[j]))
    } else {
      effect <- rank_biserial_direct(healthy, cancer)
      p_value <- suppressWarnings(wilcox.test(healthy, cancer, exact = FALSE)$p.value)
      stopifnot(abs(effect - plot$rank_biserial_Healthy_vs_CRC[j]) < 1e-10)
      stopifnot(abs(p_value - plot$wilcox_p[j]) < 1e-12)
    }
    verified_rows <- verified_rows + 1L
  }

  q_col <- grep("^wilcox_BH_q", names(plot), value = TRUE)
  stopifnot(length(q_col) == 1L)
  estimable <- !is.na(plot$wilcox_p)
  stopifnot(max(abs(p.adjust(plot$wilcox_p[estimable], method = "BH") - plot[[q_col]][estimable])) < 1e-12)
  if (any(!estimable)) stopifnot(all(is.na(plot[[q_col]][!estimable])))
}

stopifnot(verified_rows == 28L)
cat("PASS: 7 KOs; 28 cohort estimates; 11,529 raw observations; exact source-copy and r/p/q checks\n")
