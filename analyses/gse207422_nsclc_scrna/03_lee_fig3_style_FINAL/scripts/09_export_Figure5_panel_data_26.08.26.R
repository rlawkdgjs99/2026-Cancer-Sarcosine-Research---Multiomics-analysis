#!/usr/bin/env Rscript

# Export the exact cell-level or plotted-term data used by manuscript Figure 5a-g.
# This script does not recompute module scores, differential expression, enrichment,
# thresholds, or statistics. It only selects and labels values from frozen outputs.

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript 09_export_Figure5_panel_data_26.08.26.R <Figure_Panel_Data_26.08.26>")
}
output_root <- normalizePath(args[[1L]], mustWork = TRUE)

suppressPackageStartupMessages(library(data.table))

sha256 <- function(path) {
  if (!file.exists(path)) stop("Missing required file: ", path)
  out <- system2("shasum", c("-a", "256", path), stdout = TRUE, stderr = TRUE)
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) stop("SHA-256 command failed: ", path)
  sub(" .*", "", out[[1L]])
}

assert_hash <- function(path, expected) {
  observed <- sha256(path)
  if (!identical(observed, expected)) {
    stop("Frozen source SHA-256 mismatch for ", path, ": ", observed)
  }
  observed
}

table_dir <- file.path(analysis_dir, "results", "tables")
quartile_dir <- file.path(
  analysis_dir, "lineage_quartile_GOBP_26.08.26", "results", "tables"
)

source_files <- c(
  scores = file.path(analysis_dir, "intermediate", "01_cell_paper_style_scores.rds"),
  umap_caps = file.path(table_dir, "plotdata_01_UMAP_display_caps.csv"),
  ratio_cap = file.path(table_dir, "plotdata_01_ratio_violin_display_cap.csv"),
  violin = file.path(table_dir, "plotdata_01_violin_scores.csv.gz"),
  global_go = file.path(table_dir, "plotdata_03_top_GO_BP_terms.csv"),
  lineage_go = file.path(quartile_dir, "plotdata_top_GO_BP_terms.csv")
)
expected_hashes <- c(
  scores = "44a1c56f947910c2930ec0f1ce1e27a3b36cf059facd1835a7ae29ea05068262",
  umap_caps = "e1555187e6508786341acbfc2e36911e41ce03dce69d0ab24be80c38bb7c79bd",
  ratio_cap = "1055c93fdea3f6c0627c4947906d41bba8c01092162afd252db9ac2bf268d1e5",
  violin = "dabb7a22923ee13c883d7cba3a82685ed8994eb68beb3e73f7a6519bdc8961ea",
  global_go = "12632a3167816690ff41bde1065d7edccbc4c1cb5b82cee1768446af3f6f00e2",
  lineage_go = "3cb9e159efba6856ccca8823f46a4b0995544fdb41599bafe149133a8ffb027f"
)
observed_hashes <- mapply(assert_hash, source_files, expected_hashes, USE.NAMES = TRUE)

scores <- as.data.table(readRDS(source_files[["scores"]]))
required_score_cols <- c(
  "cell_id", "Sample", "Patient", "final_lineage", "analysis_eligible",
  "umap_1", "umap_2", "production_module_shifted",
  "degradation_module_shifted", "ratio_defined",
  "production_degradation_ratio", "ratio_group"
)
if (!all(required_score_cols %in% names(scores))) stop("Frozen score table schema is incomplete")
if (nrow(scores) != 92053L || anyDuplicated(scores$cell_id)) {
  stop("Frozen score table has unexpected dimensions or duplicate cell IDs")
}
if (scores[, sum(analysis_eligible)] != 90512L || scores[, sum(ratio_defined)] != 92052L) {
  stop("Frozen score eligibility/ratio counts do not match the verified analysis")
}
if (any(!is.finite(scores$umap_1)) || any(!is.finite(scores$umap_2))) {
  stop("Non-finite UMAP coordinate in frozen score table")
}

umap_caps <- fread(source_files[["umap_caps"]])
if (nrow(umap_caps) != 3L || !all(c("axis", "display_cap_99_5") %in% names(umap_caps))) {
  stop("Unexpected UMAP-cap table")
}
cap_for <- function(axis_name) {
  z <- umap_caps[axis == axis_name, display_cap_99_5]
  if (length(z) != 1L || !is.finite(z)) stop("Missing/invalid UMAP cap: ", axis_name)
  z
}
prod_umap_cap <- cap_for("Production")
deg_umap_cap <- cap_for("Degradation")

ratio_cap <- fread(source_files[["ratio_cap"]])
if (nrow(ratio_cap) != 1L || ratio_cap$quantile != 0.99 || ratio_cap$cells != 90511L) {
  stop("Unexpected ratio-violin cap table")
}
ratio_median <- median(scores[ratio_defined == TRUE, production_degradation_ratio])

# Figure 5a: four UMAP views share these coordinates and score columns. The
# lineage view excludes non-final/unresolved lineages; the other three views
# draw all 92,053 cells, exactly as the installed plotting script does.
fig5a <- scores[, .(
  cell_id, Sample, Patient, final_lineage, analysis_eligible,
  included_in_lineage_umap = analysis_eligible,
  included_in_ratio_group_umap = TRUE,
  included_in_production_umap = TRUE,
  included_in_degradation_umap = TRUE,
  umap_1, umap_2,
  production_module_shifted,
  production_display_cap_99_5 = prod_umap_cap,
  production_value_display = pmin(production_module_shifted, prod_umap_cap),
  degradation_module_shifted,
  degradation_display_cap_99_5 = deg_umap_cap,
  degradation_value_display = pmin(degradation_module_shifted, deg_umap_cap),
  ratio_defined, production_degradation_ratio, ratio_group,
  ratio_group_cutoff_median = ratio_median
)]
setorder(fig5a, cell_id)
if (nrow(fig5a) != 92053L || fig5a[, sum(included_in_lineage_umap)] != 90512L) {
  stop("Figure 5a export dimensions are wrong")
}

# Figure 5b: exact cells and 99th-percentile display winsorization used by the
# dedicated production/degradation-ratio violin.
fig5b <- scores[analysis_eligible == TRUE & ratio_defined == TRUE, .(
  cell_id, Sample, Patient, final_lineage, ratio_defined,
  production_degradation_ratio, ratio_group,
  ratio_display_cap_99 = ratio_cap$display_cap,
  ratio_value_display = pmin(production_degradation_ratio, ratio_cap$display_cap)
)]
setorder(fig5b, final_lineage, cell_id)
if (nrow(fig5b) != 90511L || anyDuplicated(fig5b$cell_id)) {
  stop("Figure 5b export dimensions are wrong")
}

# Figures 5c-d: select the two exact long-form plot-data blocks generated by
# the verified violin pipeline (including the plotted capped values).
violin <- fread(source_files[["violin"]])
if (nrow(violin) != 271535L || !all(c(
  "axis", "cell_id", "final_lineage", "value", "display_cap_99_5",
  "maximum", "cells", "value_display"
) %in% names(violin))) stop("Unexpected verified violin plot-data table")

make_module_export <- function(axis_name, value_name, display_name) {
  z <- copy(violin[axis == axis_name])
  if (nrow(z) != 90512L || anyDuplicated(z$cell_id)) {
    stop("Unexpected cell count for ", axis_name, " violin")
  }
  idx <- match(z$cell_id, scores$cell_id)
  if (anyNA(idx)) stop("Violin cell IDs do not map to frozen score table")
  z[, `:=`(Sample = scores$Sample[idx], Patient = scores$Patient[idx])]
  setnames(z, c("value", "maximum", "cells", "value_display"),
           c(value_name, "raw_maximum", "n_cells_in_panel", display_name))
  z <- z[, c(
    "cell_id", "Sample", "Patient", "final_lineage", value_name,
    "display_cap_99_5", "raw_maximum", "n_cells_in_panel", display_name
  ), with = FALSE]
  setorder(z, final_lineage, cell_id)
  z
}
fig5c <- make_module_export(
  "Production", "production_module_shifted", "production_value_display"
)
fig5d <- make_module_export(
  "Degradation", "degradation_module_shifted", "degradation_value_display"
)

# Figure 5e: exactly the 10 plotted GO:BP terms per global ratio direction.
fig5e <- fread(source_files[["global_go"]])
if (nrow(fig5e) != 20L || !setequal(
  fig5e$direction, c("High-ratio cells", "Low-ratio cells")
)) stop("Unexpected Figure 5e GO:BP plot-data table")
fig5e[, neg_log10_BH_q := -log10(q)]
setorder(fig5e, direction, q, p, -gene_fraction)

# Figures 5f-g: exactly the 10 Q4 and 10 Q1 terms plotted for each lineage.
lineage_go <- fread(source_files[["lineage_go"]])
if (nrow(lineage_go) != 40L || !setequal(lineage_go$lineage, c("Epithelial", "CAF")) ||
    !setequal(lineage_go$direction, c("Top 25%", "Bottom 25%"))) {
  stop("Unexpected lineage-restricted GO:BP plot-data table")
}
if (any(abs(lineage_go$gene_fraction - lineage_go$DE / lineage_go$N) > 1e-15)) {
  stop("Frozen lineage GO gene_fraction is not exactly DE/N")
}
lineage_go[, `:=`(
  quartile = fifelse(direction == "Top 25%", "Q4", "Q1"),
  neg_log10_BH_q = -log10(q)
)]
fig5f <- copy(lineage_go[lineage == "Epithelial"])
fig5g <- copy(lineage_go[lineage == "CAF"])
setorder(fig5f, direction, q, p, -gene_fraction)
setorder(fig5g, direction, q, p, -gene_fraction)
if (nrow(fig5f) != 20L || nrow(fig5g) != 20L) {
  stop("Figure 5f/g term counts are wrong")
}

exports <- list(
  Fig5a = fig5a, Fig5b = fig5b, Fig5c = fig5c, Fig5d = fig5d,
  Fig5e = fig5e, Fig5f = fig5f, Fig5g = fig5g
)
filenames <- setNames(
  paste0(names(exports), "_작도용_및_통계용.csv"), names(exports)
)

output_paths <- character(length(exports))
for (panel in names(exports)) {
  panel_dir <- file.path(output_root, panel)
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
  output_paths[[panel]] <- file.path(panel_dir, filenames[[panel]])
  fwrite(exports[[panel]], output_paths[[panel]], na = "NA", bom = FALSE)
}

# Independent CSV round-trip checks before declaring success.
for (panel in names(exports)) {
  back <- fread(output_paths[[panel]], na.strings = "NA")
  original <- exports[[panel]]
  if (nrow(back) != nrow(original) || ncol(back) != ncol(original) ||
      !identical(names(back), names(original))) {
    stop("CSV round-trip dimension/header failure: ", panel)
  }
  if ("cell_id" %in% names(back) && anyDuplicated(back$cell_id)) {
    stop("Duplicate cell ID after CSV round trip: ", panel)
  }
  numeric_cols <- names(original)[vapply(original, is.numeric, logical(1))]
  for (col in numeric_cols) {
    x <- original[[col]]
    y <- back[[col]]
    keep <- !is.na(x)
    tolerance <- 1e-12 * pmax(1, abs(x[keep]))
    if (!identical(is.na(x), is.na(y)) ||
        any(abs(x[keep] - y[keep]) > tolerance)) {
      stop("Numeric round-trip mismatch in ", panel, ": ", col)
    }
  }
}

source_label <- c(
  Fig5a = "frozen score RDS + frozen UMAP display caps",
  Fig5b = "frozen score RDS + frozen ratio-violin display cap",
  Fig5c = "verified violin plot-data (Production rows)",
  Fig5d = "verified violin plot-data (Degradation rows)",
  Fig5e = "verified global High/Low GO:BP plot-data",
  Fig5f = "verified lineage-quartile GO:BP plot-data (Epithelial rows)",
  Fig5g = "verified lineage-quartile GO:BP plot-data (CAF rows)"
)
source_hash_label <- c(
  Fig5a = paste(observed_hashes[c("scores", "umap_caps")], collapse = ";"),
  Fig5b = paste(observed_hashes[c("scores", "ratio_cap")], collapse = ";"),
  Fig5c = observed_hashes[["violin"]],
  Fig5d = observed_hashes[["violin"]],
  Fig5e = observed_hashes[["global_go"]],
  Fig5f = observed_hashes[["lineage_go"]],
  Fig5g = observed_hashes[["lineage_go"]]
)
filter_label <- c(
  Fig5a = "All 92,053 cells; analysis_eligible flags the 90,512-cell lineage UMAP",
  Fig5b = "analysis_eligible == TRUE and ratio_defined == TRUE",
  Fig5c = "analysis_eligible == TRUE; Production axis",
  Fig5d = "analysis_eligible == TRUE; Degradation axis",
  Fig5e = "Frozen displayed top 10 GO:BP terms per High/Low direction",
  Fig5f = "Frozen displayed top 10 GO:BP terms per Epithelial Q4/Q1 direction",
  Fig5g = "Frozen displayed top 10 GO:BP terms per CAF Q4/Q1 direction"
)

manifest <- rbindlist(lapply(names(exports), function(panel) {
  data.table(
    panel = panel,
    output_file = file.path(panel, filenames[[panel]]),
    rows = nrow(exports[[panel]]),
    columns = ncol(exports[[panel]]),
    column_names = paste(names(exports[[panel]]), collapse = "; "),
    source = source_label[[panel]],
    source_sha256 = source_hash_label[[panel]],
    inclusion_or_selection_rule = filter_label[[panel]],
    output_sha256 = sha256(output_paths[[panel]]),
    verification = "PASS: frozen-source identity, schema/count assertions, CSV round trip"
  )
}))
manifest_file <- file.path(output_root, "FIGURE_5_PROVENANCE.csv")
fwrite(manifest, manifest_file, bom = FALSE)

message("Figure 5 panel-data export complete")
print(manifest[, .(panel, rows, columns, output_sha256)])
