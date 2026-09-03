#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
base_dir <- normalizePath(file.path(analysis_dir, ".."), mustWork = TRUE)
analysis_root <- normalizePath(file.path(base_dir, ".."), mustWork = TRUE)
prior_dir <- file.path(analysis_root, "02_lineage_reannotation")
output_root <- Sys.getenv("LINEAGE_QUARTILE_OUTPUT_DIR", unset = analysis_dir)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(file.path(base_dir, "R_libs"), file.path(prior_dir, "R_libs"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

options(stringsAsFactors = FALSE, warn = 1, future.globals.maxSize = 16 * 1024^3)
set.seed(260826)

table_dir <- file.path(output_root, "results", "tables")
pub_dir <- file.path(output_root, "results", "figures_publication")
diag_dir <- file.path(output_root, "results", "figures_diagnostic")
log_dir <- file.path(output_root, "logs")
for (d in c(table_dir, pub_dir, diag_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

sha256 <- function(path) {
  sub(" .*", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
}

plan_file <- file.path(analysis_dir, "ANALYSIS_PLAN_FROZEN.md")
expected_plan_sha <- "35d9220054c68d493f7513bc3d275fe163e81aaa3279082b75785fbd5b0d4fd2"
if (!identical(sha256(plan_file), expected_plan_sha)) stop("Frozen plan SHA-256 mismatch")

counts_file <- file.path(prior_dir, "intermediate", "01_full_counts_mt20_qc.rds")
lineage_file <- file.path(prior_dir, "results", "tables", "09_final_cell_lineages_FROZEN.csv")
score_file <- file.path(base_dir, "intermediate", "01_cell_paper_style_scores.rds")
input_files <- c(counts = counts_file, lineages = lineage_file, scores = score_file)
expected_input_sha <- c(
  counts = "7552ba882c306fa68380390a8172e0fc96af2c7be7cb3f5778345204e9fd61b2",
  lineages = "97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1",
  scores = "44a1c56f947910c2930ec0f1ce1e27a3b36cf059facd1835a7ae29ea05068262"
)
if (!all(file.exists(input_files))) stop("One or more frozen inputs are missing")
observed_input_sha <- vapply(input_files, sha256, character(1))
if (!identical(observed_input_sha, expected_input_sha)) stop("Frozen input SHA-256 mismatch")
fwrite(data.table(
  input = names(input_files), path = unname(input_files),
  sha256 = unname(observed_input_sha), verified = observed_input_sha == expected_input_sha
), file.path(table_dir, "00_input_manifest.csv"))

message("Loading frozen score and lineage records")
score <- as.data.table(readRDS(score_file))
lin <- fread(lineage_file)
if (nrow(score) != 92053L || uniqueN(score$cell_id) != nrow(score)) {
  stop("Unexpected score dimensions or duplicate cell IDs")
}
required_score <- c(
  "cell_id", "Patient", "Pathologic.Response", "final_lineage",
  "production_degradation_ratio", "ratio_defined"
)
if (!all(required_score %in% names(score))) stop("Score schema is incomplete")
if (!setequal(score$cell_id, lin$cell_id)) stop("Score/lineage cell-ID sets differ")

lineages <- c("Epithelial", "CAF")
target <- copy(score[final_lineage %chin% lineages & is.finite(production_degradation_ratio)])
if (target[, uniqueN(final_lineage)] != 2L) stop("Target lineages are incomplete")
target[, `:=`(
  q1 = as.numeric(quantile(production_degradation_ratio, 0.25, type = 7)),
  q3 = as.numeric(quantile(production_degradation_ratio, 0.75, type = 7))
), by = final_lineage]
target[, quartile_group := fifelse(
  production_degradation_ratio <= q1, "Bottom 25%",
  fifelse(production_degradation_ratio >= q3, "Top 25%", "Middle 50%")
)]
target[, quartile_group := factor(quartile_group, levels = c("Bottom 25%", "Middle 50%", "Top 25%"))]

quartile_summary <- target[, .(
  finite_ratio_cells = .N,
  q1 = unique(q1), q3 = unique(q3),
  bottom_cells = sum(quartile_group == "Bottom 25%"),
  middle_cells = sum(quartile_group == "Middle 50%"),
  top_cells = sum(quartile_group == "Top 25%"),
  q1_ties = sum(production_degradation_ratio == unique(q1)),
  q3_ties = sum(production_degradation_ratio == unique(q3)),
  total_patients = uniqueN(Patient),
  bottom_patients = uniqueN(Patient[quartile_group == "Bottom 25%"]),
  top_patients = uniqueN(Patient[quartile_group == "Top 25%"])
), by = final_lineage]
fwrite(quartile_summary, file.path(table_dir, "01_quartile_thresholds_and_counts.csv"))

patient_contribution <- target[quartile_group != "Middle 50%", .(
  cells = .N,
  min_ratio = min(production_degradation_ratio),
  median_ratio = median(production_degradation_ratio),
  max_ratio = max(production_degradation_ratio),
  response = paste(sort(unique(na.omit(Pathologic.Response))), collapse = ";")
), by = .(final_lineage, quartile_group, Patient)]
fwrite(patient_contribution, file.path(table_dir, "01_patient_contribution_by_quartile.csv"))

selected <- target[quartile_group != "Middle 50%"]
if (selected[, uniqueN(cell_id)] != nrow(selected)) stop("Selected-cell IDs are not unique")
fwrite(selected[, .(
  cell_id, Patient, Pathologic.Response, final_lineage,
  production_degradation_ratio, q1, q3, quartile_group
)], file.path(table_dir, "01_selected_cell_quartile_membership.csv.gz"))

message("Loading frozen sparse counts and selecting lineage quartiles")
counts <- readRDS(counts_file)
if (!inherits(counts, "dgCMatrix") || nrow(counts) != 24292L || ncol(counts) != 92053L) {
  stop("Frozen count object has unexpected class or dimensions")
}
if (anyDuplicated(colnames(counts)) || !setequal(colnames(counts), score$cell_id)) {
  stop("Count/score cell-ID integrity failure")
}
selected_ids <- selected$cell_id
counts_sub <- counts[, match(selected_ids, colnames(counts)), drop = FALSE]
if (!identical(colnames(counts_sub), selected_ids)) stop("Selected count alignment failed")
rm(counts)
invisible(gc())

meta <- as.data.frame(selected[, .(cell_id, Patient, final_lineage, quartile_group)])
rownames(meta) <- meta$cell_id
meta$cell_id <- NULL
obj <- CreateSeuratObject(counts = counts_sub, meta.data = meta, min.cells = 0, min.features = 0)
rm(counts_sub, meta)
invisible(gc())
obj <- NormalizeData(
  obj, assay = "RNA", normalization.method = "LogNormalize",
  scale.factor = 10000, margin = 1, verbose = TRUE
)
DefaultAssay(obj) <- "RNA"

run_one_lineage <- function(lineage_name) {
  message("Running Top-versus-Bottom analysis for ", lineage_name)
  cells_lineage <- selected[final_lineage == lineage_name, cell_id]
  obj_lineage <- subset(obj, cells = cells_lineage)
  Idents(obj_lineage) <- "quartile_group"
  observed_idents <- sort(as.character(unique(Idents(obj_lineage))))
  if (!identical(observed_idents, sort(c("Bottom 25%", "Top 25%")))) {
    stop("Unexpected quartile identities for ", lineage_name)
  }

  deg <- FindMarkers(
    obj_lineage, ident.1 = "Top 25%", ident.2 = "Bottom 25%",
    assay = "RNA", slot = "data", test.use = "wilcox",
    logfc.threshold = 0.1, min.pct = 0.01, min.diff.pct = -Inf,
    only.pos = FALSE, max.cells.per.ident = Inf,
    random.seed = 260826, densify = FALSE, verbose = TRUE
  )
  deg <- as.data.table(deg, keep.rownames = "gene")
  if (!all(c("gene", "avg_log2FC", "p_val_adj") %in% names(deg))) {
    stop("FindMarkers output schema is incomplete for ", lineage_name)
  }
  deg[, direction := fifelse(
    p_val_adj < 0.05 & avg_log2FC >= 0.25, "Top 25%",
    fifelse(p_val_adj < 0.05 & avg_log2FC <= -0.25, "Bottom 25%", "Not selected")
  )]

  entrez <- AnnotationDbi::mapIds(
    org.Hs.eg.db, keys = deg$gene, keytype = "SYMBOL", column = "ENTREZID",
    multiVals = "first"
  )
  deg[, ENTREZID := unname(entrez[gene])]
  universe <- unique(na.omit(deg$ENTREZID))
  top_genes <- unique(na.omit(deg[direction == "Top 25%", ENTREZID]))
  bottom_genes <- unique(na.omit(deg[direction == "Bottom 25%", ENTREZID]))
  if (!length(top_genes) || !length(bottom_genes)) {
    stop("No selected Entrez genes in one or both directions for ", lineage_name)
  }

  go <- as.data.table(
    limma::goana(list(Top = top_genes, Bottom = bottom_genes), universe = universe, species = "Hs"),
    keep.rownames = "GO_ID"
  )
  required_go <- c("GO_ID", "Term", "Ont", "N", "Top", "P.Top", "Bottom", "P.Bottom")
  if (!all(required_go %in% names(go))) stop("Unexpected goana schema for ", lineage_name)
  go <- go[Ont == "BP"]
  go[, `:=`(
    BH_Top = p.adjust(P.Top, method = "BH"),
    BH_Bottom = p.adjust(P.Bottom, method = "BH"),
    lineage = lineage_name
  )]

  go_long <- rbindlist(list(
    go[, .(
      lineage = lineage_name, GO_ID, Term, N, direction = "Top 25%",
      DE = Top, p = P.Top, q = BH_Top,
      selected_entrez = length(top_genes), universe_entrez = length(universe)
    )],
    go[, .(
      lineage = lineage_name, GO_ID, Term, N, direction = "Bottom 25%",
      DE = Bottom, p = P.Bottom, q = BH_Bottom,
      selected_entrez = length(bottom_genes), universe_entrez = length(universe)
    )]
  ))
  go_long[, `:=`(
    gene_fraction = DE / N,
    background_fraction = selected_entrez / universe_entrez,
    fold_enrichment = (DE / N) / (selected_entrez / universe_entrez)
  )]

  audit <- data.table(
    lineage = lineage_name,
    metric = c(
      "bottom_cells", "top_cells", "tested_genes", "mapped_universe_entrez",
      "selected_top_symbols", "selected_bottom_symbols",
      "selected_top_entrez", "selected_bottom_entrez"
    ),
    value = c(
      sum(selected$final_lineage == lineage_name & selected$quartile_group == "Bottom 25%"),
      sum(selected$final_lineage == lineage_name & selected$quartile_group == "Top 25%"),
      nrow(deg), length(universe),
      sum(deg$direction == "Top 25%"), sum(deg$direction == "Bottom 25%"),
      length(top_genes), length(bottom_genes)
    )
  )

  list(deg = deg[, lineage := lineage_name], go = go, go_long = go_long, audit = audit)
}

analyses <- lapply(lineages, run_one_lineage)
names(analyses) <- lineages
deg_all <- rbindlist(lapply(analyses, `[[`, "deg"), use.names = TRUE, fill = TRUE)
go_all <- rbindlist(lapply(analyses, `[[`, "go"), use.names = TRUE, fill = TRUE)
go_long_all <- rbindlist(lapply(analyses, `[[`, "go_long"), use.names = TRUE, fill = TRUE)
deg_audit <- rbindlist(lapply(analyses, `[[`, "audit"))

fwrite(deg_all, file.path(table_dir, "02_lineage_quartile_FindMarkers_all.csv"))
fwrite(go_all, file.path(table_dir, "03_lineage_quartile_GO_BP_all_terms.csv"))
fwrite(go_long_all, file.path(table_dir, "03_lineage_quartile_GO_BP_directional_all.csv"))
fwrite(deg_audit, file.path(table_dir, "02_DEG_enrichment_input_audit.csv"))

plot_data <- go_long_all[
  N >= 10 & N <= 500 & DE > 0 & is.finite(q) & q < 0.05 & is.finite(fold_enrichment)
]
setorder(plot_data, lineage, direction, q, p, -fold_enrichment)
plot_data <- plot_data[, head(.SD, 10L), by = .(lineage, direction)]
fwrite(plot_data, file.path(table_dir, "plotdata_top_GO_BP_terms.csv"))
if (!nrow(plot_data)) stop("No significant GO BP terms passed the publication filter")

neglog_range <- range(-log10(pmax(plot_data$q, .Machine$double.xmin)), finite = TRUE)
size_range <- range(plot_data$DE, finite = TRUE)
COL_LOW_Q <- "#D5E5F2"
COL_HIGH_Q <- "#C43C3C"
COL_TEXT <- "#202020"

theme_nc <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(colour = COL_TEXT),
      axis.text = element_text(size = 6, colour = COL_TEXT),
      axis.text.y = element_text(size = 5.8),
      axis.title = element_text(size = 7, colour = COL_TEXT),
      axis.line = element_line(linewidth = 0.35, colour = COL_TEXT),
      axis.ticks = element_line(linewidth = 0.35, colour = COL_TEXT),
      plot.title = element_text(size = 7.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6, colour = "#555555"),
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6),
      plot.margin = margin(4, 5, 4, 4, "pt")
    )
}

make_go_panel <- function(lineage_name, direction_name) {
  z <- copy(plot_data[lineage == lineage_name & direction == direction_name])
  direction_short <- if (identical(direction_name, "Top 25%")) "Q4" else "Q1"
  if (!nrow(z)) {
    return(ggplot() + theme_void() +
      annotate("text", x = 0.5, y = 0.5, label = "No BH q < 0.05 GO BP terms", size = 2.4) +
      ggtitle(paste0(lineage_name, " — ", direction_short)))
  }
  z[, Term_plot := factor(Term, levels = rev(Term))]
  ggplot(z, aes(fold_enrichment, Term_plot)) +
    geom_vline(xintercept = 1, colour = "#A0A0A0", linewidth = 0.3, linetype = 2) +
    geom_point(aes(size = DE, colour = -log10(pmax(q, .Machine$double.xmin))), alpha = 0.95) +
    scale_colour_gradient(
      low = COL_LOW_Q, high = COL_HIGH_Q, limits = neglog_range,
      name = expression(-log[10]("BH q"))
    ) +
    scale_size_continuous(limits = size_range, range = c(1.8, 5.4), name = "DE genes") +
    labs(
      title = paste0(lineage_name, " — ", direction_short),
      x = "Fold enrichment", y = NULL
    ) +
    theme_nc()
}

panels <- list(
  epithelial_top = make_go_panel("Epithelial", "Top 25%"),
  epithelial_bottom = make_go_panel("Epithelial", "Bottom 25%"),
  caf_top = make_go_panel("CAF", "Top 25%"),
  caf_bottom = make_go_panel("CAF", "Bottom 25%")
)

save_plot <- function(stem, plot, width, height, directory = pub_dir) {
  png_file <- file.path(directory, paste0(stem, ".png"))
  pdf_file <- file.path(directory, paste0(stem, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, units = "in", dpi = 600,
         device = "png", bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, units = "in",
         device = cairo_pdf, bg = "white")
  if (!file.exists(png_file) || !file.exists(pdf_file)) stop("Figure write failed: ", stem)

  srgb_profile <- "/System/Library/ColorSync/Profiles/sRGB Profile.icc"
  sips_bin <- Sys.which("sips")
  if (nzchar(sips_bin) && file.exists(srgb_profile)) {
    tmp_png <- paste0(png_file, ".srgb.png")
    out <- system2(
      sips_bin,
      c("-m", shQuote(srgb_profile), shQuote(png_file), "--out", shQuote(tmp_png)),
      stdout = TRUE, stderr = TRUE
    )
    status <- attr(out, "status")
    if (!is.null(status) && status != 0L) stop("sRGB embedding failed for ", png_file)
    if (!file.exists(tmp_png)) stop("sRGB derivative was not created for ", png_file)
    unlink(png_file)
    if (!file.rename(tmp_png, png_file)) stop("Could not install sRGB PNG: ", png_file)
  }
  xattr_bin <- Sys.which("xattr")
  if (nzchar(xattr_bin)) system2(xattr_bin, c("-c", shQuote(png_file)), stdout = FALSE, stderr = FALSE)
  chflags_bin <- Sys.which("chflags")
  if (nzchar(chflags_bin)) system2(chflags_bin, c("nohidden", shQuote(png_file)), stdout = FALSE, stderr = FALSE)
}

epithelial_plot <- (panels$epithelial_top | panels$epithelial_bottom) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Epithelial-cell GO Biological Process enrichment",
    subtitle = "Within-lineage top versus bottom quartile of production/degradation ratio",
    theme = theme(
      plot.title = element_text(family = "Arial", size = 8, face = "bold", colour = COL_TEXT),
      plot.subtitle = element_text(family = "Arial", size = 6.2, colour = "#555555")
    )
  ) & theme(legend.position = "bottom")

caf_plot <- (panels$caf_top | panels$caf_bottom) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "CAF GO Biological Process enrichment",
    subtitle = "Within-lineage top versus bottom quartile of production/degradation ratio",
    theme = theme(
      plot.title = element_text(family = "Arial", size = 8, face = "bold", colour = COL_TEXT),
      plot.subtitle = element_text(family = "Arial", size = 6.2, colour = "#555555")
    )
  ) & theme(legend.position = "bottom")

combined_plot <- ((panels$epithelial_top | panels$epithelial_bottom) /
                    (panels$caf_top | panels$caf_bottom)) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Lineage-restricted GO Biological Process enrichment",
    subtitle = "Top and bottom quartiles of the production/degradation ratio within each lineage",
    theme = theme(
      plot.title = element_text(family = "Arial", size = 8.2, face = "bold", colour = COL_TEXT),
      plot.subtitle = element_text(family = "Arial", size = 6.2, colour = "#555555")
    )
  ) & theme(legend.position = "bottom")

save_plot("Fig_Epithelial_ratio_Q4_vs_Q1_GO_BP", epithelial_plot, 7.20, 4.65)
save_plot("Fig_CAF_ratio_Q4_vs_Q1_GO_BP", caf_plot, 7.20, 4.65)
save_plot("Fig_Epithelial_CAF_ratio_Q4_vs_Q1_GO_BP_combined", combined_plot, 7.20, 8.40)

selection_plot <- ggplot(target, aes(1, production_degradation_ratio, fill = final_lineage)) +
  geom_violin(width = 0.82, colour = COL_TEXT, linewidth = 0.35, trim = TRUE) +
  geom_hline(
    data = unique(target[, .(final_lineage, q1, q3)]),
    aes(yintercept = q1), linewidth = 0.35, linetype = 2, colour = "#2E5F8A"
  ) +
  geom_hline(
    data = unique(target[, .(final_lineage, q1, q3)]),
    aes(yintercept = q3), linewidth = 0.35, linetype = 2, colour = "#C43C3C"
  ) +
  facet_wrap(~final_lineage, scales = "free_y") +
  scale_fill_manual(values = c(Epithelial = "#E56B8A", CAF = "#D97850")) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 0.05)) +
  labs(
    title = "Within-lineage ratio-quartile definition",
    subtitle = "Blue dashed line: Q1; red dashed line: Q3",
    x = NULL, y = "Production/degradation ratio"
  ) +
  theme_nc() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 7, face = "bold")
  )
save_plot("Diagnostic_ratio_quartile_definition", selection_plot, 5.20, 3.20, directory = diag_dir)

package_versions <- data.table(
  package = c("R", "Seurat", "data.table", "Matrix", "limma", "AnnotationDbi", "org.Hs.eg.db", "ggplot2", "patchwork"),
  version = c(
    paste(R.version$major, R.version$minor, sep = "."),
    vapply(c("Seurat", "data.table", "Matrix", "limma", "AnnotationDbi", "org.Hs.eg.db", "ggplot2", "patchwork"),
           function(x) as.character(packageVersion(x)), character(1))
  )
)
fwrite(package_versions, file.path(table_dir, "00_package_versions.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "01_analysis_sessionInfo.txt"))

core_files <- sort(list.files(file.path(output_root, "results"), recursive = TRUE, full.names = TRUE))
core_files <- core_files[file.info(core_files)$isdir %in% FALSE]
hash_table <- data.table(
  relative_path = sub(paste0("^", normalizePath(output_root), "/?"), "", normalizePath(core_files)),
  bytes = file.info(core_files)$size,
  sha256 = vapply(core_files, sha256, character(1))
)
fwrite(hash_table, file.path(log_dir, "output_sha256.csv"))

message("Lineage-restricted quartile GO BP analysis complete")
