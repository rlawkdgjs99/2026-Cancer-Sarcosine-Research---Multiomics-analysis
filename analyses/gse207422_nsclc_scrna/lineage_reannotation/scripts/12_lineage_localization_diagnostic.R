#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
.libPaths(c(file.path(analysis_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
})

options(stringsAsFactors = FALSE, warn = 1)

table_dir <- file.path(analysis_dir, "results", "tables")
figure_dir <- file.path(analysis_dir, "results", "figures_diagnostic")
intermediate_dir <- file.path(analysis_dir, "intermediate")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

counts_file <- file.path(intermediate_dir, "01_full_counts_mt20_qc.rds")
lineage_file <- file.path(table_dir, "09_final_cell_lineages_FROZEN.csv")
if (!all(file.exists(c(counts_file, lineage_file)))) stop("Required inputs are missing")

counts <- readRDS(counts_file)
lineages <- fread(lineage_file)
if (!inherits(counts, "dgCMatrix")) stop("Counts are not a sparse dgCMatrix")
if (ncol(counts) != 92053L || nrow(lineages) != 92053L) stop("Unexpected input dimensions")

ord <- match(colnames(counts), lineages$cell_id)
if (anyNA(ord)) stop("Count/lineage cell-ID mismatch")
lineages <- lineages[ord]
if (!identical(colnames(counts), lineages$cell_id)) stop("Count/lineage ordering failed")
if (!identical(as.numeric(Matrix::colSums(counts)), as.numeric(lineages$library_size))) {
  stop("Lineage-table library sizes disagree with counts")
}

targets <- c("SARDH", "PIPOX")
if (!all(targets %in% rownames(counts))) stop("A target gene is absent from the matrix")

# Descriptive patient-level localization. Patient-lineage combinations with at
# least 10 cells are summarized; cells are never treated as replicates.
eligible <- which(lineages$analysis_eligible == TRUE)
if (length(unique(lineages$final_lineage[eligible])) != 14L) {
  stop("Expected 14 analysis-eligible lineages")
}

cell_rows <- vector("list", length(targets))
for (i in seq_along(targets)) {
  target <- targets[i]
  raw <- as.numeric(counts[target, eligible])
  lib <- as.numeric(lineages$library_size[eligible])
  if (any(lib <= 0)) stop("Non-positive cell library size encountered")
  cell_rows[[i]] <- data.table(
    Patient = lineages$Patient[eligible],
    Sample = lineages$Sample[eligible],
    final_lineage = lineages$final_lineage[eligible],
    target = target,
    detected = as.integer(raw > 0),
    log1p_CP10k = log1p(raw / lib * 10000)
  )
}
cell_values <- rbindlist(cell_rows)

patient_lineage <- cell_values[, .(
  cells = .N,
  detected_cells = as.integer(sum(detected)),
  detection_pct = as.numeric(mean(detected) * 100),
  mean_log1p_CP10k = as.numeric(mean(log1p_CP10k))
), by = .(Patient, Sample, final_lineage, target)]
patient_lineage[, eligible_ge10_cells := cells >= 10L]

if (any(patient_lineage$detected_cells > patient_lineage$cells)) {
  stop("Detected-cell count exceeds total cells")
}
if (any(patient_lineage$detection_pct < 0 | patient_lineage$detection_pct > 100)) {
  stop("Detection percentage outside [0, 100]")
}

# Within-patient ranks quantify whether localization is consistent across
# patients. Rank 1 is the highest-valued lineage in that patient.
ranked <- patient_lineage[eligible_ge10_cells == TRUE]
ranked[, detection_rank := frank(-detection_pct, ties.method = "average"),
       by = .(Patient, target)]
ranked[, expression_rank := frank(-mean_log1p_CP10k, ties.method = "average"),
       by = .(Patient, target)]
ranked[, in_detection_top3 := detection_rank <= 3]
ranked[, in_expression_top3 := expression_rank <= 3]

summary_table <- ranked[, .(
  patients = .N,
  median_cells = as.numeric(median(cells)),
  median_detection_pct = as.numeric(median(detection_pct)),
  Q1_detection_pct = as.numeric(quantile(detection_pct, 0.25, names = FALSE)),
  Q3_detection_pct = as.numeric(quantile(detection_pct, 0.75, names = FALSE)),
  median_mean_log1p_CP10k = as.numeric(median(mean_log1p_CP10k)),
  Q1_mean_log1p_CP10k = as.numeric(quantile(mean_log1p_CP10k, 0.25, names = FALSE)),
  Q3_mean_log1p_CP10k = as.numeric(quantile(mean_log1p_CP10k, 0.75, names = FALSE)),
  median_detection_rank = as.numeric(median(detection_rank)),
  detection_top3_fraction = as.numeric(mean(in_detection_top3)),
  median_expression_rank = as.numeric(median(expression_rank)),
  expression_top3_fraction = as.numeric(mean(in_expression_top3))
), by = .(final_lineage, target)]

setorder(patient_lineage, target, Patient, final_lineage)
setorder(ranked, target, Patient, detection_rank, final_lineage)
setorder(summary_table, target, median_detection_rank, -median_detection_pct, final_lineage)
fwrite(patient_lineage, file.path(table_dir, "12_patient_lineage_target_celllevel_summary.csv"))
fwrite(ranked, file.path(table_dir, "12_patient_within_lineage_target_ranks.csv"))
fwrite(summary_table, file.path(table_dir, "12_lineage_localization_summary.csv"))

# Diagnostic only: point area = median cell detection percentage; fill = median
# patient-level mean log1p(CP10k). This is not an inferential main-figure panel.
lineage_order <- summary_table[, .(
  order_score = mean(median_detection_rank, na.rm = TRUE)
), by = final_lineage][order(order_score, final_lineage)]$final_lineage
plot_dt <- copy(summary_table)
plot_dt[, final_lineage := factor(final_lineage, levels = rev(lineage_order))]
plot_dt[, target := factor(target, levels = targets)]

p <- ggplot(plot_dt, aes(x = target, y = final_lineage)) +
  geom_point(aes(size = median_detection_pct, fill = median_mean_log1p_CP10k),
             shape = 21, colour = "#202020", stroke = 0.35) +
  scale_size_continuous(range = c(1.5, 10), limits = c(0, NA),
                        name = "Median detected cells (%)") +
  scale_fill_gradient(low = "#F2F2F2", high = "#2E5F8A",
                      name = "Median mean\nlog1p(CP10k)") +
  labs(x = NULL, y = NULL,
       caption = "Descriptive patient-level summary; only patient-lineages with >=10 cells") +
  theme_classic(base_family = "Arial", base_size = 8) +
  theme(
    axis.text = element_text(colour = "black", size = 7),
    axis.text.x = element_text(face = "italic"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6.5),
    plot.caption = element_text(size = 6, colour = "#404040", hjust = 0),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

png_file <- file.path(figure_dir, "12_lineage_localization_dotplot.png")
pdf_file <- file.path(figure_dir, "12_lineage_localization_dotplot.pdf")
ggsave(png_file, p, width = 4.6, height = 4.4, units = "in", dpi = 600,
       device = "png", bg = "white")
ggsave(pdf_file, p, width = 4.6, height = 4.4, units = "in", device = cairo_pdf,
       bg = "white")
if (!all(file.exists(c(png_file, pdf_file)))) stop("Diagnostic figure export failed")

audit <- data.table(
  check = c(
    "input_cells", "analysis_eligible_cells", "eligible_lineages",
    "targets", "patients", "patient_lineage_rows_ge10", "output_png_exists",
    "output_pdf_exists"
  ),
  value = c(
    ncol(counts), length(eligible), length(unique(lineages$final_lineage[eligible])),
    length(targets), length(unique(lineages$Patient)), nrow(ranked),
    file.exists(png_file), file.exists(pdf_file)
  )
)
fwrite(audit, file.path(table_dir, "12_lineage_localization_audit.csv"))

message("PASS: lineage-localization diagnostic completed")
message("  Patient-level eligible rows: ", nrow(ranked))
message("  Summary rows: ", nrow(summary_table))
