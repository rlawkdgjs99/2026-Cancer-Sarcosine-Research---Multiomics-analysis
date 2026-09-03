#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
base_dir <- normalizePath(file.path(analysis_dir, ".."), mustWork = TRUE)
analysis_root <- normalizePath(file.path(base_dir, ".."), mustWork = TRUE)
prior_dir <- file.path(analysis_root, "02_lineage_reannotation")
output_root <- Sys.getenv("HALLMARK_GSEA_OUTPUT_DIR", unset = analysis_dir)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(
  file.path(analysis_dir, "R_libs"), file.path(base_dir, "R_libs"),
  file.path(prior_dir, "R_libs"), .libPaths()
))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(fgsea)
  library(BiocParallel)
  library(limma)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE, warn = 1, future.globals.maxSize = 16 * 1024^3)
RNGkind("L'Ecuyer-CMRG")
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
assert_equal <- function(observed, expected, label, tolerance = NULL) {
  ok <- if (is.null(tolerance)) identical(observed, expected) else {
    length(observed) == length(expected) && all(abs(observed - expected) <= tolerance)
  }
  if (!isTRUE(ok)) {
    stop(label, " mismatch. Observed: ", paste(observed, collapse = ", "),
         "; expected: ", paste(expected, collapse = ", "))
  }
}

plan_file <- file.path(analysis_dir, "ANALYSIS_PLAN_FROZEN.md")
expected_plan_sha <- "5ca4942bb26e7b45f7a0ab4c7f46f8718a34e705ca3a1a581012bae0d40f0cb4"
assert_equal(sha256(plan_file), expected_plan_sha, "Frozen plan SHA-256")

counts_file <- file.path(prior_dir, "intermediate", "01_full_counts_mt20_qc.rds")
lineage_file <- file.path(prior_dir, "results", "tables", "09_final_cell_lineages_FROZEN.csv")
score_file <- file.path(base_dir, "intermediate", "01_cell_paper_style_scores.rds")
resource_csv <- file.path(analysis_dir, "resources", "MSigDB_Hallmark_2026.1.Hs_gene_symbols.csv")
resource_gmt <- file.path(analysis_dir, "resources", "MSigDB_Hallmark_2026.1.Hs_gene_symbols.gmt")
input_files <- c(
  counts = counts_file, lineages = lineage_file, scores = score_file,
  hallmark_csv = resource_csv, hallmark_gmt = resource_gmt
)
expected_sha <- c(
  counts = "7552ba882c306fa68380390a8172e0fc96af2c7be7cb3f5778345204e9fd61b2",
  lineages = "97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1",
  scores = "44a1c56f947910c2930ec0f1ce1e27a3b36cf059facd1835a7ae29ea05068262",
  hallmark_csv = "476a3ad64f6fa2d8cc6a0a418191e288224b68b7d589177e63c897919d0fe8da",
  hallmark_gmt = "e3452a2952123f8f44fcf29cc0d050846e493d5a19a454d00e72d043723d9a9b"
)
if (!all(file.exists(input_files))) stop("One or more frozen inputs are missing")
observed_sha <- vapply(input_files, sha256, character(1))
assert_equal(observed_sha, expected_sha, "Frozen input/resource SHA-256")
fwrite(data.table(
  input = names(input_files), path = unname(input_files),
  sha256 = unname(observed_sha), verified = unname(observed_sha == expected_sha)
), file.path(table_dir, "00_input_manifest.csv"))

packages <- c(
  "R", "data.table", "Matrix", "Seurat", "fgsea", "BiocParallel",
  "limma", "ggplot2", "patchwork"
)
versions <- c(
  paste(R.version$major, R.version$minor, sep = "."),
  vapply(packages[-1], function(x) as.character(packageVersion(x)), character(1))
)
fwrite(data.table(package = packages, version = versions),
       file.path(table_dir, "00_package_versions.csv"))
if (as.character(packageVersion("fgsea")) != "1.38.0") stop("Unexpected fgsea version")

message("Loading frozen scores, lineages and Hallmark resource")
score <- as.data.table(readRDS(score_file))
lin <- fread(lineage_file)
resource <- fread(resource_csv)
if (nrow(score) != 92053L || uniqueN(score$cell_id) != 92053L) {
  stop("Unexpected score dimensions or duplicate cell IDs")
}
if (nrow(lin) != 92053L || uniqueN(lin$cell_id) != 92053L ||
    !setequal(score$cell_id, lin$cell_id)) {
  stop("Frozen lineage table is inconsistent with the score table")
}
required_score <- c(
  "cell_id", "Patient", "final_lineage", "production_degradation_ratio", "ratio_defined"
)
if (!all(required_score %in% names(score))) stop("Score schema is incomplete")
if (nrow(resource) != 7322L || uniqueN(resource$gs_name) != 50L ||
    resource[, anyDuplicated(paste(gs_name, gene_symbol, sep = "\r"))]) {
  stop("Frozen Hallmark resource integrity failure")
}

lineages <- c("Epithelial", "CAF")
target <- copy(score[
  final_lineage %chin% lineages & ratio_defined == TRUE &
    is.finite(production_degradation_ratio)
])
target[, `:=`(
  q1 = as.numeric(quantile(production_degradation_ratio, 0.25, type = 7)),
  q3 = as.numeric(quantile(production_degradation_ratio, 0.75, type = 7))
), by = final_lineage]
target[, quartile := fifelse(
  production_degradation_ratio <= q1, "Q1",
  fifelse(production_degradation_ratio >= q3, "Q4", "Middle")
)]

quartile_summary <- target[, .(
  finite_ratio_cells = .N,
  q1 = unique(q1), q3 = unique(q3),
  q1_cells = sum(quartile == "Q1"),
  q4_cells = sum(quartile == "Q4"),
  middle_cells = sum(quartile == "Middle"),
  q1_boundary_ties = sum(production_degradation_ratio == unique(q1)),
  q3_boundary_ties = sum(production_degradation_ratio == unique(q3)),
  all_patients = uniqueN(Patient),
  q1_patients = uniqueN(Patient[quartile == "Q1"]),
  q4_patients = uniqueN(Patient[quartile == "Q4"])
), by = final_lineage]
quartile_summary[, lineage_order__ := match(final_lineage, lineages)]
setorder(quartile_summary, lineage_order__)
quartile_summary[, lineage_order__ := NULL]
expected_counts <- data.table(
  final_lineage = lineages,
  finite_ratio_cells = c(11912L, 1102L),
  q1 = c(0.376174483868376, 0.304701212119105),
  q3 = c(0.536637221540880, 0.501039946395598),
  q1_cells = c(2978L, 276L), q4_cells = c(2978L, 276L)
)
for (nm in c("finite_ratio_cells", "q1_cells", "q4_cells")) {
  assert_equal(as.integer(quartile_summary[[nm]]), as.integer(expected_counts[[nm]]),
               paste("Quartile", nm))
}
for (nm in c("q1", "q3")) {
  assert_equal(quartile_summary[[nm]], expected_counts[[nm]], paste("Quartile", nm), 1e-15)
}
if (any(quartile_summary$q1_boundary_ties != 0L | quartile_summary$q3_boundary_ties != 0L)) {
  stop("Unexpected cells exactly tied to a quartile boundary")
}
fwrite(quartile_summary, file.path(table_dir, "01_quartile_thresholds_and_counts.csv"))

selected <- target[quartile %chin% c("Q1", "Q4")]
setorder(selected, final_lineage, quartile, cell_id)
if (nrow(selected) != 6508L || uniqueN(selected$cell_id) != 6508L ||
    length(intersect(selected[quartile == "Q1", cell_id], selected[quartile == "Q4", cell_id]))) {
  stop("Selected-cell membership integrity failure")
}
membership <- selected[, .(
  cell_id, Patient, final_lineage, production_degradation_ratio, q1, q3, quartile
)]
fwrite(membership, file.path(table_dir, "01_selected_cell_quartile_membership.csv"))
patient_contribution <- membership[, .(
  cells = .N,
  min_ratio = min(production_degradation_ratio),
  median_ratio = median(production_degradation_ratio),
  max_ratio = max(production_degradation_ratio)
), by = .(final_lineage, quartile, Patient)]
fwrite(patient_contribution, file.path(table_dir, "01_patient_contribution_audit.csv"))

method_contract <- data.table(
  field = c(
    "biological_question", "contrast", "grouping", "response_labels_used",
    "normalization", "gene_ranking", "primary_rank_exclusions", "gsea_engine",
    "hallmark_release", "multiplicity", "replication_boundary"
  ),
  value = c(
    "Within-lineage Hallmark programs associated with sarcosine ratio Q4 versus Q1",
    "Q4 minus Q1, separately in Epithelial and CAF",
    "Within-lineage type-7 quartiles of the finite shifted AddModuleScore Production/Degradation ratio",
    "No",
    "Seurat LogNormalize, scale factor 10000",
    "Seurat FindMarkers Wilcoxon avg_log2FC; no DEG threshold before GSEA",
    "GNMT;DMGDH;SARDH;PIPOX",
    "fgseaMultilevel; minSize 15; maxSize 500; eps 0; scoreType std; nproc 1; seed 260826",
    "MSigDB Human Hallmark 2026.1.Hs, 50 sets",
    "BH within each lineage (primary); pooled BH across both lineages (sensitivity)",
    "Exploratory cell-level association; cells are nested within 15 patients"
  )
)
fwrite(method_contract, file.path(table_dir, "00_method_contract.csv"))

message("Loading frozen sparse counts and selecting 6,508 Q1/Q4 cells")
counts <- readRDS(counts_file)
if (!inherits(counts, "dgCMatrix") || nrow(counts) != 24292L || ncol(counts) != 92053L ||
    anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)) ||
    !setequal(colnames(counts), score$cell_id)) {
  stop("Frozen count object integrity failure")
}
idx <- match(membership$cell_id, colnames(counts))
if (anyNA(idx)) stop("Selected cell missing from count matrix")
counts_sub <- counts[, idx, drop = FALSE]
if (!identical(colnames(counts_sub), membership$cell_id)) stop("Selected count alignment failure")
rm(counts)
invisible(gc())

meta <- as.data.frame(membership[, .(cell_id, Patient, final_lineage, quartile)])
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

pathways <- split(resource$gene_symbol, resource$gs_name)
pathways <- lapply(pathways, unique)
if (length(pathways) != 50L) stop("Expected exactly 50 frozen Hallmark pathways")
score_defining_genes <- c("GNMT", "DMGDH", "SARDH", "PIPOX")

make_unique_rank <- function(marker_dt, lineage_name, exclude_score_genes) {
  z <- marker_dt[is.finite(avg_log2FC), .(gene, avg_log2FC)]
  if (exclude_score_genes) z <- z[!gene %chin% score_defining_genes]
  if (anyDuplicated(z$gene)) stop("Duplicate genes in rank for ", lineage_name)
  vals <- sort(unique(z$avg_log2FC))
  gaps <- diff(vals)
  min_gap <- if (any(gaps > 0)) min(gaps[gaps > 0]) else Inf
  rank_scale <- max(1, max(abs(z$avg_log2FC)))
  epsilon <- min(1e-12 * rank_scale, min_gap / 4)
  if (!is.finite(epsilon) || epsilon <= 0) epsilon <- 1e-12 * rank_scale
  setorder(z, avg_log2FC, gene)
  z[, tie_n := .N, by = avg_log2FC]
  z[, tie_offset := if (.N == 1L) 0 else seq(-epsilon, epsilon, length.out = .N),
    by = avg_log2FC]
  z[, rank_final := avg_log2FC + tie_offset]
  if (anyDuplicated(z$rank_final)) stop("Tie perturbation did not produce a unique rank for ", lineage_name)
  if (max(abs(z$tie_offset)) > 1e-12 * rank_scale * (1 + 1e-12)) {
    stop("Tie perturbation exceeded the frozen bound")
  }
  setorder(z, -rank_final, gene)
  rank_vec <- z$rank_final
  names(rank_vec) <- z$gene
  list(
    table = z,
    vector = rank_vec,
    audit = data.table(
      lineage = lineage_name,
      score_genes_excluded = exclude_score_genes,
      genes = nrow(z),
      tied_genes = sum(z$tie_n > 1L),
      tie_groups = uniqueN(z[tie_n > 1L, avg_log2FC]),
      max_tie_group = if (any(z$tie_n > 1L)) max(z$tie_n) else 1L,
      perturbation_epsilon = epsilon,
      maximum_absolute_perturbation = max(abs(z$tie_offset)),
      rank_scale = rank_scale,
      unique_final_rank = !anyDuplicated(z$rank_final)
    )
  )
}

run_fgsea <- function(rank_vec, lineage_name, rank_type) {
  set.seed(260826)
  ans <- as.data.table(fgseaMultilevel(
    pathways = pathways, stats = rank_vec,
    minSize = 15, maxSize = 500, eps = 0,
    scoreType = "std", nproc = 1,
    BPPARAM = BiocParallel::SerialParam(progressbar = FALSE)
  ))
  if (nrow(ans) != 50L || uniqueN(ans$pathway) != 50L) {
    stop("fgsea did not return all 50 Hallmark sets for ", lineage_name)
  }
  ans[, `:=`(
    lineage = lineage_name,
    rank_type = rank_type,
    leadingEdge_genes = vapply(leadingEdge, paste, collapse = ";", character(1)),
    leading_edge_count = lengths(leadingEdge),
    q_within_lineage = p.adjust(pval, method = "BH")
  )]
  ans[, leadingEdge := NULL]
  setcolorder(ans, c(
    "lineage", "rank_type", "pathway", "size", "ES", "NES", "pval",
    "q_within_lineage", "log2err", "leading_edge_count", "leadingEdge_genes"
  ))
  ans
}

run_camera <- function(rank_vec, lineage_name) {
  indices <- lapply(pathways, function(g) which(names(rank_vec) %chin% g))
  indices <- indices[lengths(indices) >= 15L & lengths(indices) <= 500L]
  cam <- as.data.table(cameraPR(rank_vec, index = indices, use.ranks = TRUE, sort = FALSE),
                       keep.rownames = "pathway")
  if (nrow(cam) != 50L || uniqueN(cam$pathway) != 50L) {
    stop("cameraPR did not return all 50 Hallmark sets for ", lineage_name)
  }
  cam[, `:=`(
    lineage = lineage_name,
    q_within_lineage = p.adjust(PValue, method = "BH")
  )]
  setcolorder(cam, c("lineage", "pathway", "NGenes", "Direction", "PValue", "q_within_lineage"))
  cam
}

run_one_lineage <- function(lineage_name) {
  message("Differential ranking and Hallmark GSEA: ", lineage_name)
  ids <- membership[final_lineage == lineage_name, cell_id]
  obj_l <- subset(obj, cells = ids)
  Idents(obj_l) <- "quartile"
  if (!setequal(levels(droplevels(Idents(obj_l))), c("Q1", "Q4"))) {
    stop("Unexpected identities for ", lineage_name)
  }
  markers <- FindMarkers(
    obj_l, ident.1 = "Q4", ident.2 = "Q1", assay = "RNA", slot = "data",
    test.use = "wilcox", features = rownames(obj_l),
    logfc.threshold = 0, min.pct = 0, min.diff.pct = -Inf,
    only.pos = FALSE, max.cells.per.ident = Inf,
    random.seed = 260826, densify = FALSE, verbose = TRUE
  )
  markers <- as.data.table(markers, keep.rownames = "gene")
  if (!all(c("gene", "avg_log2FC", "p_val", "p_val_adj") %in% names(markers))) {
    stop("FindMarkers output schema is incomplete for ", lineage_name)
  }
  markers[, lineage := lineage_name]
  primary <- make_unique_rank(markers, lineage_name, TRUE)
  included <- make_unique_rank(markers, lineage_name, FALSE)
  primary$table[, `:=`(lineage = lineage_name, rank_type = "primary_score_genes_excluded")]
  included$table[, `:=`(lineage = lineage_name, rank_type = "sensitivity_score_genes_included")]

  pathway_map <- data.table(
    lineage = lineage_name,
    pathway = names(pathways),
    resource_genes = lengths(pathways),
    primary_rank_genes = vapply(pathways, function(g) sum(g %chin% names(primary$vector)), integer(1)),
    sensitivity_rank_genes = vapply(pathways, function(g) sum(g %chin% names(included$vector)), integer(1))
  )

  gsea_primary <- run_fgsea(primary$vector, lineage_name, "primary_score_genes_excluded")
  gsea_included <- run_fgsea(included$vector, lineage_name, "sensitivity_score_genes_included")
  gsea_reversed <- run_fgsea(-primary$vector, lineage_name, "sensitivity_primary_rank_reversed")
  camera <- run_camera(primary$vector, lineage_name)

  list(
    markers = markers,
    primary_rank = primary$table,
    included_rank = included$table,
    tie_audit = rbind(primary$audit, included$audit),
    pathway_map = pathway_map,
    gsea_primary = gsea_primary,
    gsea_included = gsea_included,
    gsea_reversed = gsea_reversed,
    camera = camera
  )
}

analyses <- lapply(lineages, run_one_lineage)
names(analyses) <- lineages
rm(obj)
invisible(gc())

markers_all <- rbindlist(lapply(analyses, `[[`, "markers"), use.names = TRUE, fill = TRUE)
primary_rank_all <- rbindlist(lapply(analyses, `[[`, "primary_rank"), use.names = TRUE, fill = TRUE)
included_rank_all <- rbindlist(lapply(analyses, `[[`, "included_rank"), use.names = TRUE, fill = TRUE)
tie_audit <- rbindlist(lapply(analyses, `[[`, "tie_audit"), use.names = TRUE, fill = TRUE)
pathway_map <- rbindlist(lapply(analyses, `[[`, "pathway_map"), use.names = TRUE, fill = TRUE)
gsea_primary <- rbindlist(lapply(analyses, `[[`, "gsea_primary"), use.names = TRUE, fill = TRUE)
gsea_included <- rbindlist(lapply(analyses, `[[`, "gsea_included"), use.names = TRUE, fill = TRUE)
gsea_reversed <- rbindlist(lapply(analyses, `[[`, "gsea_reversed"), use.names = TRUE, fill = TRUE)
camera_all <- rbindlist(lapply(analyses, `[[`, "camera"), use.names = TRUE, fill = TRUE)

gsea_primary[, q_global_100_tests := p.adjust(pval, method = "BH")]
gsea_included[, q_global_100_tests := p.adjust(pval, method = "BH")]
gsea_reversed[, q_global_100_tests := p.adjust(pval, method = "BH")]
camera_all[, q_global_100_tests := p.adjust(PValue, method = "BH")]

fwrite(markers_all, file.path(table_dir, "02_FindMarkers_all_genes.csv"))
fwrite(primary_rank_all, file.path(table_dir, "03_primary_rank_score_genes_excluded.csv"))
fwrite(included_rank_all, file.path(table_dir, "03_sensitivity_rank_score_genes_included.csv"))
fwrite(tie_audit, file.path(table_dir, "03_rank_tie_audit.csv"))
fwrite(pathway_map, file.path(table_dir, "04_hallmark_rank_intersection_counts.csv"))
fwrite(gsea_primary, file.path(table_dir, "05_fgsea_primary_all_50_sets.csv"))
fwrite(gsea_included, file.path(table_dir, "05_fgsea_sensitivity_score_genes_included.csv"))
fwrite(gsea_reversed, file.path(table_dir, "05_fgsea_sensitivity_reversed_rank.csv"))
fwrite(camera_all, file.path(table_dir, "06_cameraPR_sensitivity_all_50_sets.csv"))

direction_agreement <- merge(
  gsea_primary[, .(lineage, pathway, NES_primary = NES, q_primary = q_within_lineage)],
  gsea_included[, .(lineage, pathway, NES_score_genes_included = NES,
                    q_score_genes_included = q_within_lineage)],
  by = c("lineage", "pathway")
)
direction_agreement <- merge(
  direction_agreement,
  camera_all[, .(lineage, pathway, camera_direction = Direction,
                 camera_q = q_within_lineage)],
  by = c("lineage", "pathway")
)
direction_agreement[, `:=`(
  fgsea_direction_stable = sign(NES_primary) == sign(NES_score_genes_included),
  camera_direction_expected = fifelse(NES_primary > 0, "Up", "Down"),
  camera_direction_agrees = camera_direction == fifelse(NES_primary > 0, "Up", "Down")
)]
fwrite(direction_agreement, file.path(table_dir, "06_method_and_rank_sensitivity_agreement.csv"))

reverse_check <- merge(
  gsea_primary[, .(lineage, pathway, NES_primary = NES)],
  gsea_reversed[, .(lineage, pathway, NES_reversed = NES)],
  by = c("lineage", "pathway")
)
reverse_check[, `:=`(
  expected_reversal = -NES_primary,
  absolute_error = abs(NES_reversed + NES_primary),
  sign_reversed = sign(NES_reversed) == -sign(NES_primary)
)]
fwrite(reverse_check, file.path(table_dir, "06_reversed_rank_NES_check.csv"))

plot_data <- copy(gsea_primary[is.finite(q_within_lineage) & q_within_lineage < 0.05])
plot_data[, direction := fifelse(NES > 0, "Q4 (high ratio)", "Q1 (low ratio)")]
plot_data[, abs_NES__ := abs(NES)]
setorder(plot_data, lineage, direction, q_within_lineage, -abs_NES__, pathway)
plot_data <- plot_data[, head(.SD, 10L), by = .(lineage, direction)]
plot_data[, abs_NES__ := NULL]
plot_data[, pathway_label := tools::toTitleCase(tolower(gsub("_", " ", sub("^HALLMARK_", "", pathway))))]
fwrite(plot_data, file.path(table_dir, "plotdata_07_significant_hallmark_terms.csv"))

display_audit <- data.table(
  lineage = rep(lineages, each = 2L),
  direction = rep(c("Q4 (high ratio)", "Q1 (low ratio)"), times = 2L)
)
display_audit[, significant_pathways := mapply(
  function(l, d) nrow(gsea_primary[
    lineage == l & q_within_lineage < 0.05 &
      fifelse(NES > 0, "Q4 (high ratio)", "Q1 (low ratio)") == d
  ]), lineage, direction
)]
display_audit[, displayed_pathways := mapply(
  function(l, d) nrow(plot_data[lineage == l & direction == d]), lineage, direction
)]
fwrite(display_audit, file.path(table_dir, "plotdata_07_display_audit.csv"))

COL_Q1 <- "#9868B1"
COL_Q4 <- "#642F7F"
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
      strip.background = element_blank(),
      strip.text = element_text(size = 6.5, face = "bold", colour = COL_TEXT),
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6),
      plot.margin = margin(4, 6, 4, 4, "pt")
    )
}

make_lineage_panel <- function(lineage_name) {
  z <- copy(plot_data[lineage == lineage_name])
  missing_directions <- display_audit[lineage == lineage_name & significant_pathways == 0L, direction]
  if (!nrow(z)) {
    return(
      ggplot() + theme_void(base_family = "Arial") +
        annotate("text", x = 0.5, y = 0.55,
                 label = "No Hallmark pathway passed within-lineage BH q < 0.05",
                 size = 2.5, colour = COL_TEXT) +
        labs(title = lineage_name) +
        theme(plot.title = element_text(size = 6.5, face = "bold", hjust = 0))
    )
  }
  z[, pathway_axis := factor(pathway_label, levels = rev(unique(pathway_label)))]
  q_strength <- -log10(pmax(z$q_within_lineage, .Machine$double.xmin))
  z[, q_strength := q_strength]
  subtitle <- if (length(missing_directions)) {
    paste("No q < 0.05:", paste(missing_directions, collapse = "; "))
  } else NULL
  ggplot(z, aes(x = NES, y = pathway_axis)) +
    geom_vline(xintercept = 0, linewidth = 0.35, colour = "#777777") +
    geom_segment(aes(x = 0, xend = NES, yend = pathway_axis),
                 linewidth = 0.5, colour = "#A0A0A0") +
    geom_point(aes(size = leading_edge_count, fill = direction, alpha = q_strength),
               shape = 21, colour = COL_TEXT, stroke = 0.35) +
    scale_fill_manual(values = c(
      "Q4 (high ratio)" = COL_Q4,
      "Q1 (low ratio)" = COL_Q1
    ), name = "Enriched in") +
    scale_alpha_continuous(range = c(0.60, 1.00), name = expression(-log[10](q))) +
    scale_size_continuous(range = c(2.0, 5.0), name = "Leading-edge\ngenes") +
    labs(title = lineage_name, subtitle = subtitle, x = "Normalized enrichment score (NES)", y = NULL) +
    theme_nc() +
    theme(
      plot.title = element_text(size = 6.5, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 5.8, colour = "#555555")
    )
}

p_epi <- make_lineage_panel("Epithelial")
p_caf <- make_lineage_panel("CAF")
combined_plot <- (p_epi / p_caf) +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(legend.position = "bottom")

save_plot <- function(stem, plot, width, height) {
  png_file <- file.path(pub_dir, paste0(stem, ".png"))
  pdf_file <- file.path(pub_dir, paste0(stem, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, units = "in", dpi = 600,
         device = "png", bg = "white", limitsize = FALSE)
  ggsave(pdf_file, plot, width = width, height = height, units = "in",
         device = cairo_pdf, bg = "white", limitsize = FALSE)
  if (!file.exists(png_file) || !file.exists(pdf_file)) stop("Figure write failed: ", stem)
  srgb_profile <- "/System/Library/ColorSync/Profiles/sRGB Profile.icc"
  if (nzchar(Sys.which("sips")) && file.exists(srgb_profile)) {
    tmp_png <- paste0(png_file, ".srgb.png")
    status <- system2("sips", c("-m", shQuote(srgb_profile), shQuote(png_file),
                                "--out", shQuote(tmp_png)), stdout = TRUE, stderr = TRUE)
    if (!file.exists(tmp_png)) stop("sRGB conversion failed: ", paste(status, collapse = "\n"))
    if (!file.rename(tmp_png, png_file)) stop("Could not install sRGB PNG")
  }
  invisible(system2("xattr", c("-c", shQuote(png_file)), stdout = TRUE, stderr = TRUE))
  invisible(system2("chflags", c("nohidden", shQuote(png_file)), stdout = TRUE, stderr = TRUE))
  c(png = png_file, pdf = pdf_file)
}

figures <- save_plot(
  "Fig_Hallmark_GSEA_Epithelial_CAF_ratio_Q4_vs_Q1",
  combined_plot, width = 7.20, height = 5.20
)

summary_table <- gsea_primary[, .(
  tested_pathways = .N,
  q4_enriched_q_lt_0_05 = sum(NES > 0 & q_within_lineage < 0.05),
  q1_enriched_q_lt_0_05 = sum(NES < 0 & q_within_lineage < 0.05),
  smallest_q = min(q_within_lineage, na.rm = TRUE),
  pathway_at_smallest_q = pathway[which.min(q_within_lineage)],
  NES_at_smallest_q = NES[which.min(q_within_lineage)]
), by = lineage]
fwrite(summary_table, file.path(table_dir, "07_primary_result_summary.csv"))

figure_manifest <- data.table(
  file = unname(figures),
  type = names(figures),
  sha256 = vapply(unname(figures), sha256, character(1)),
  intended_width_in = 7.20,
  intended_height_in = 5.20,
  dpi = c(600, NA_integer_)
)
fwrite(figure_manifest, file.path(table_dir, "08_figure_manifest.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "01_sessionInfo.txt"))

machine_files <- sort(c(
  list.files(table_dir, full.names = TRUE, pattern = "\\.csv$"),
  list.files(pub_dir, full.names = TRUE, pattern = "\\.png$")
))
hash_manifest <- data.table(
  file = machine_files,
  sha256 = vapply(machine_files, sha256, character(1))
)
fwrite(hash_manifest, file.path(log_dir, "01_machine_readable_sha256.csv"))

message("Hallmark GSEA workflow complete")
print(summary_table)
