#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
analysis_root <- normalizePath(file.path(analysis_dir, ".."), mustWork = TRUE)
prior_dir <- file.path(analysis_root, "02_lineage_reannotation")
.libPaths(c(file.path(analysis_dir, "R_libs"), file.path(prior_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

table_dir <- file.path(analysis_dir, "results", "tables")
pub_dir <- file.path(analysis_dir, "results", "figures_publication")
log_dir <- file.path(analysis_dir, "logs")
for (d in c(table_dir, pub_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

score_file <- file.path(table_dir, "01_cell_paper_style_scores.csv.gz")
lineage_file <- file.path(prior_dir, "results", "tables", "09_final_cell_lineages_FROZEN.csv")
if (!all(file.exists(c(score_file, lineage_file)))) stop("Required verified input is missing")

sha256 <- function(path) sub(" .*", "", system2("shasum", c("-a", "256", path), stdout = TRUE))
expected_sha <- c(
  scores = "17724b2b1fcbad436e6ce4cee91a34bda851361aac7b14c26a3557e5d006fb61",
  lineages = "97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1"
)
observed_sha <- c(scores = sha256(score_file), lineages = sha256(lineage_file))
if (!identical(observed_sha, expected_sha)) stop("Verified input SHA-256 mismatch")

COL_TEXT <- "#202020"
lineage_order <- c(
  "Epithelial", "CAF", "B cell", "Plasma cell", "CD4 T cell",
  "CD8 T cell", "Cycling T cell", "NK cell", "Mast cell",
  "Neutrophil", "Monocyte", "Macrophage", "Conventional DC", "pDC"
)
lineage_cols <- setNames(hcl.colors(length(lineage_order), palette = "Dark 3"), lineage_order)

theme_nc <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = COL_TEXT),
      axis.text = element_text(size = 6, colour = COL_TEXT),
      axis.title = element_text(size = 7, colour = COL_TEXT),
      axis.line = element_line(linewidth = 0.35, colour = COL_TEXT),
      axis.ticks = element_line(linewidth = 0.35, colour = COL_TEXT),
      axis.ticks.length = unit(1.3, "mm"),
      plot.title = element_text(size = 7.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6, hjust = 0, colour = "#555555"),
      plot.margin = margin(3, 3, 3, 3, unit = "pt")
    )
}

save_plot <- function(stem, plot, width = 7.20, height = 3.25) {
  png_file <- file.path(pub_dir, paste0(stem, ".png"))
  pdf_file <- file.path(pub_dir, paste0(stem, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, units = "in", dpi = 600,
         device = "png", bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, units = "in",
         device = cairo_pdf, bg = "white")
  if (!file.exists(png_file) || !file.exists(pdf_file)) stop("Figure write failed: ", stem)
}

scores <- fread(
  score_file,
  select = c(
    "cell_id", "final_lineage", "analysis_eligible",
    "production_module_shifted", "degradation_module_shifted"
  )
)
if (nrow(scores) != 92053L || anyDuplicated(scores$cell_id)) {
  stop("Unexpected score-table dimensions or duplicate cell IDs")
}
plot_cells <- scores[analysis_eligible == TRUE]
if (nrow(plot_cells) != 90512L) stop("Unexpected number of analysis-eligible cells")
if (!setequal(unique(plot_cells$final_lineage), lineage_order)) stop("Unexpected eligible lineage set")
plot_cells[, final_lineage := factor(final_lineage, levels = lineage_order)]

long <- rbindlist(list(
  plot_cells[, .(
    cell_id, final_lineage, module = "Production",
    genes = "GNMT + DMGDH", raw_value = production_module_shifted
  )],
  plot_cells[, .(
    cell_id, final_lineage, module = "Degradation",
    genes = "SARDH + PIPOX", raw_value = degradation_module_shifted
  )]
))
if (any(!is.finite(long$raw_value)) || any(long$raw_value < 0)) {
  stop("Module-score display input contains non-finite or negative shifted values")
}

caps <- long[, .(
  display_rule = "Upper 0.5% winsorized for display only",
  quantile = 0.995,
  display_cap = as.numeric(quantile(raw_value, 0.995, na.rm = TRUE)),
  raw_maximum = max(raw_value),
  cells = .N
), by = .(module, genes)]
long <- merge(long, caps[, .(module, display_cap)], by = "module", all.x = TRUE, sort = FALSE)
long[, display_value := pmin(raw_value, display_cap)]
setorder(long, module, final_lineage, cell_id)
fwrite(caps, file.path(table_dir, "plotdata_06_dedicated_module_violin_display_caps.csv"))

make_violin <- function(module_name, stem, y_title) {
  d <- long[module == module_name]
  gene_label <- unique(d$genes)
  if (length(gene_label) != 1L) stop("Unexpected module gene-label count")
  p <- ggplot(d, aes(final_lineage, display_value, fill = final_lineage)) +
    geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, colour = COL_TEXT) +
    geom_boxplot(
      width = 0.08, outlier.shape = NA, fill = "white", alpha = 0.80,
      linewidth = 0.28, colour = COL_TEXT
    ) +
    scale_fill_manual(values = lineage_cols, drop = FALSE) +
    labs(
      title = paste0(module_name, " module score across cell lineages"),
      subtitle = paste0(
        gene_label,
        "; cell-level shifted AddModuleScore; upper 0.5% winsorized for display only"
      ),
      x = NULL,
      y = y_title
    ) +
    theme_nc() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.5)
    )
  save_plot(stem, p)
}

make_violin(
  "Production",
  "Fig3_style_D_production_violin",
  "Production module score"
)
make_violin(
  "Degradation",
  "Fig3_style_D_degradation_violin",
  "Degradation module score"
)

# Read-only clinical-metadata audit for deciding whether patient-group analyses
# are supportable. This does not perform or imply a response-association test.
lin <- fread(
  lineage_file,
  select = c(
    "cell_id", "Patient", "Sample", "Resource", "Sex", "Age",
    "Clinical.Stage", "Pathology", "PD1.Antibody", "Chemotherapy",
    "Pathologic.Response", "Residual.Tumor", "RECIST", "final_lineage"
  )
)
if (nrow(lin) != 92053L || anyDuplicated(lin$cell_id)) {
  stop("Unexpected frozen-lineage dimensions or duplicate cell IDs")
}

patient_fields <- c(
  "Sample", "Resource", "Sex", "Age", "Clinical.Stage", "Pathology",
  "PD1.Antibody", "Chemotherapy", "Pathologic.Response", "Residual.Tumor", "RECIST"
)
uniqueness <- lin[, lapply(.SD, uniqueN), by = Patient, .SDcols = patient_fields]
if (any(as.matrix(uniqueness[, ..patient_fields]) != 1L)) {
  stop("A patient maps to multiple values for a supposedly patient-level metadata field")
}
patient_meta <- unique(lin[, c("Patient", patient_fields), with = FALSE])
cell_counts <- lin[, .(
  all_cells = .N,
  CD8_T_cells = sum(final_lineage == "CD8 T cell"),
  NK_cells = sum(final_lineage == "NK cell")
), by = Patient]
patient_audit <- merge(patient_meta, cell_counts, by = "Patient", all = TRUE, sort = TRUE)
patient_audit[, pathologic_group_candidate := fifelse(
  Pathologic.Response %chin% c("MPR", "pCR"), "MPR/pCR",
  fifelse(Pathologic.Response == "NMPR", "NMPR", "Not evaluable")
)]
patient_audit[, RECIST_group_candidate := fifelse(
  RECIST == "PR", "PR", fifelse(RECIST == "SD", "SD", "Not evaluable")
)]
setcolorder(patient_audit, c(
  "Patient", "Sample", "Resource", "Pathologic.Response", "pathologic_group_candidate",
  "RECIST", "RECIST_group_candidate", "Residual.Tumor", "PD1.Antibody", "Chemotherapy",
  "Sex", "Age", "Clinical.Stage", "Pathology", "all_cells", "CD8_T_cells", "NK_cells"
))
fwrite(patient_audit, file.path(table_dir, "06_patient_group_metadata_audit.csv"))

summarize_candidate <- function(group_col, definition) {
  patient_audit[, .(
    patients = .N,
    pre_treatment_biopsy_patients = sum(Resource == "Pre-treatment biopsy"),
    post_treatment_surgery_patients = sum(Resource == "Post-treatment surgery"),
    CD8_T_cells = sum(CD8_T_cells),
    NK_cells = sum(NK_cells)
  ), by = .(group = get(group_col))][, definition := definition]
}
candidate_summary <- rbindlist(list(
  summarize_candidate(
    "pathologic_group_candidate",
    "Pathologic response after neoadjuvant anti-PD-1 plus chemotherapy"
  ),
  summarize_candidate(
    "RECIST_group_candidate",
    "Radiographic RECIST response after neoadjuvant anti-PD-1 plus chemotherapy"
  )
), use.names = TRUE)
setcolorder(candidate_summary, c(
  "definition", "group", "patients", "pre_treatment_biopsy_patients",
  "post_treatment_surgery_patients", "CD8_T_cells", "NK_cells"
))
setorder(candidate_summary, definition, group)
fwrite(candidate_summary, file.path(table_dir, "06_candidate_patient_group_summary.csv"))

writeLines(capture.output(sessionInfo()), file.path(log_dir, "06_dedicated_module_violins_sessionInfo.txt"))
message("Dedicated module violins and read-only patient metadata audit completed")
