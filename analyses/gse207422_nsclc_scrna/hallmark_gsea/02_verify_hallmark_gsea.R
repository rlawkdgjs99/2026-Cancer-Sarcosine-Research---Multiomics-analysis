#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
base_dir <- normalizePath(file.path(analysis_dir, ".."), mustWork = TRUE)
analysis_root <- normalizePath(file.path(base_dir, ".."), mustWork = TRUE)
prior_dir <- file.path(analysis_root, "02_lineage_reannotation")
.libPaths(c(file.path(analysis_dir, "R_libs"), file.path(base_dir, "R_libs"),
            file.path(prior_dir, "R_libs"), .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(png)
})

table_dir <- file.path(analysis_dir, "results", "tables")
pub_dir <- file.path(analysis_dir, "results", "figures_publication")
log_dir <- file.path(analysis_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) {
  sub(" .*", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
}
checks <- data.table(check = character(), passed = logical(), detail = character())
record <- function(name, pass, detail) {
  checks <<- rbind(checks, data.table(check = name, passed = isTRUE(pass), detail = as.character(detail)))
}

plan_file <- file.path(analysis_dir, "ANALYSIS_PLAN_FROZEN.md")
record("frozen_plan_hash", sha256(plan_file) == "5ca4942bb26e7b45f7a0ab4c7f46f8718a34e705ca3a1a581012bae0d40f0cb4",
       sha256(plan_file))

input_files <- c(
  counts = file.path(prior_dir, "intermediate", "01_full_counts_mt20_qc.rds"),
  lineages = file.path(prior_dir, "results", "tables", "09_final_cell_lineages_FROZEN.csv"),
  scores = file.path(base_dir, "intermediate", "01_cell_paper_style_scores.rds"),
  hallmark_csv = file.path(analysis_dir, "resources", "MSigDB_Hallmark_2026.1.Hs_gene_symbols.csv"),
  hallmark_gmt = file.path(analysis_dir, "resources", "MSigDB_Hallmark_2026.1.Hs_gene_symbols.gmt")
)
expected_sha <- c(
  counts = "7552ba882c306fa68380390a8172e0fc96af2c7be7cb3f5778345204e9fd61b2",
  lineages = "97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1",
  scores = "44a1c56f947910c2930ec0f1ce1e27a3b36cf059facd1835a7ae29ea05068262",
  hallmark_csv = "476a3ad64f6fa2d8cc6a0a418191e288224b68b7d589177e63c897919d0fe8da",
  hallmark_gmt = "e3452a2952123f8f44fcf29cc0d050846e493d5a19a454d00e72d043723d9a9b"
)
observed_sha <- vapply(input_files, sha256, character(1))
record("all_input_and_resource_hashes", identical(observed_sha, expected_sha),
       paste(names(observed_sha), substr(observed_sha, 1, 12), collapse = "; "))

score <- as.data.table(readRDS(input_files[["scores"]]))
lin <- fread(input_files[["lineages"]])
record("score_dimensions_and_unique_ids", nrow(score) == 92053L && uniqueN(score$cell_id) == 92053L,
       paste(nrow(score), uniqueN(score$cell_id), sep = "/"))
record("score_lineage_cell_id_identity", nrow(lin) == 92053L && uniqueN(lin$cell_id) == 92053L &&
         setequal(score$cell_id, lin$cell_id), "92,053 IDs in both tables")

lineages <- c("Epithelial", "CAF")
target <- copy(score[
  final_lineage %chin% lineages & ratio_defined == TRUE & is.finite(production_degradation_ratio)
])
target[, `:=`(
  q1 = as.numeric(quantile(production_degradation_ratio, 0.25, type = 7)),
  q3 = as.numeric(quantile(production_degradation_ratio, 0.75, type = 7))
), by = final_lineage]
target[, quartile := fifelse(
  production_degradation_ratio <= q1, "Q1",
  fifelse(production_degradation_ratio >= q3, "Q4", "Middle")
)]
recomputed_summary <- target[, .(
  finite_ratio_cells = .N, q1 = unique(q1), q3 = unique(q3),
  q1_cells = sum(quartile == "Q1"), q4_cells = sum(quartile == "Q4"),
  q1_ties = sum(production_degradation_ratio == unique(q1)),
  q3_ties = sum(production_degradation_ratio == unique(q3))
), by = final_lineage]
setorder(recomputed_summary, final_lineage)
expected_summary <- data.table(
  final_lineage = c("CAF", "Epithelial"),
  finite_ratio_cells = c(1102L, 11912L),
  q1 = c(0.304701212119105, 0.376174483868376),
  q3 = c(0.501039946395598, 0.536637221540880),
  q1_cells = c(276L, 2978L), q4_cells = c(276L, 2978L),
  q1_ties = c(0L, 0L), q3_ties = c(0L, 0L)
)
record("quartile_thresholds_counts_and_boundary_ties",
       identical(recomputed_summary[, c("final_lineage", "finite_ratio_cells", "q1_cells", "q4_cells", "q1_ties", "q3_ties")],
                 expected_summary[, c("final_lineage", "finite_ratio_cells", "q1_cells", "q4_cells", "q1_ties", "q3_ties")]) &&
         max(abs(recomputed_summary$q1 - expected_summary$q1)) <= 1e-15 &&
         max(abs(recomputed_summary$q3 - expected_summary$q3)) <= 1e-15,
       "Exact locked Q1/Q3 selection reproduced")

membership <- fread(file.path(table_dir, "01_selected_cell_quartile_membership.csv"))
recomputed_membership <- target[quartile %chin% c("Q1", "Q4"), .(
  cell_id, Patient, final_lineage, production_degradation_ratio, q1, q3, quartile
)]
setorder(membership, final_lineage, quartile, cell_id)
setorder(recomputed_membership, final_lineage, quartile, cell_id)
record("exact_selected_cell_membership",
       identical(membership$cell_id, recomputed_membership$cell_id) &&
         identical(membership$Patient, recomputed_membership$Patient) &&
         identical(membership$final_lineage, recomputed_membership$final_lineage) &&
         identical(membership$quartile, recomputed_membership$quartile) &&
         max(abs(membership$production_degradation_ratio -
                   recomputed_membership$production_degradation_ratio)) <= 1e-12 &&
         max(abs(membership$q1 - recomputed_membership$q1)) <= 1e-12 &&
         max(abs(membership$q3 - recomputed_membership$q3)) <= 1e-12 &&
         nrow(membership) == 6508L,
       "6,508 exact IDs/labels; numeric CSV round-trip error <= 1e-12")
record("q1_q4_cell_sets_disjoint",
       length(intersect(membership[quartile == "Q1", cell_id], membership[quartile == "Q4", cell_id])) == 0L,
       "zero overlap")
record("response_labels_absent_from_design_outputs",
       !"Pathologic.Response" %in% names(membership) &&
         !any(grepl("Pathologic.Response",
                    readLines(file.path(analysis_dir, "scripts", "01_run_hallmark_GSEA.R"), warn = FALSE),
                    fixed = TRUE)),
       "Response field absent from membership and analysis script")

resource <- fread(input_files[["hallmark_csv"]])
record("hallmark_resource_dimensions_and_uniqueness",
       nrow(resource) == 7322L && uniqueN(resource$gs_name) == 50L &&
         !resource[, anyDuplicated(paste(gs_name, gene_symbol, sep = "\r"))],
       "50 pathways; 7,322 unique memberships")
pathways <- split(resource$gene_symbol, resource$gs_name)
pathways <- lapply(pathways, unique)

markers <- fread(file.path(table_dir, "02_FindMarkers_all_genes.csv"))
rank_primary <- fread(file.path(table_dir, "03_primary_rank_score_genes_excluded.csv"))
rank_included <- fread(file.path(table_dir, "03_sensitivity_rank_score_genes_included.csv"))
tie_audit <- fread(file.path(table_dir, "03_rank_tie_audit.csv"))
score_genes <- c("GNMT", "DMGDH", "SARDH", "PIPOX")
record("primary_rank_excludes_score_genes",
       !any(rank_primary$gene %chin% score_genes), "GNMT/DMGDH/SARDH/PIPOX absent")
record("sensitivity_rank_includes_score_genes",
       all(score_genes %chin% rank_included$gene), "all four score genes present")
rank_integrity <- rank_primary[, .(
  gene_ids_unique = !anyDuplicated(gene),
  final_ranks_finite = all(is.finite(rank_final)),
  final_ranks_unique = !anyDuplicated(rank_final),
  perturbation_bounded = max(abs(tie_offset)) <=
    1e-12 * max(1, max(abs(avg_log2FC))) * (1 + 1e-12),
  unperturbed_matches_markers = {
    current_lineage <- unique(lineage)
    mm <- markers[lineage == current_lineage, .(gene, avg_log2FC)]
    joined <- merge(.SD[, .(gene, avg_log2FC)], mm, by = "gene", suffixes = c("_rank", "_marker"))
    nrow(joined) == .N && max(abs(joined$avg_log2FC_rank - joined$avg_log2FC_marker)) == 0
  }
), by = lineage]
rank_check_cols <- setdiff(names(rank_integrity), "lineage")
record("primary_rank_integrity", all(unlist(rank_integrity[, ..rank_check_cols])),
       paste(rank_integrity$lineage,
             apply(rank_integrity[, ..rank_check_cols], 1, paste, collapse = "/"),
             collapse = "; "))
record("tie_audit_reports_unique_bounded_ranks",
       all(tie_audit$unique_final_rank) &&
         all(tie_audit$maximum_absolute_perturbation <= 1e-12 * tie_audit$rank_scale * (1 + 1e-12)),
       paste("rows", nrow(tie_audit)))

mapping <- fread(file.path(table_dir, "04_hallmark_rank_intersection_counts.csv"))
mapping_recomputed <- rbindlist(lapply(lineages, function(l) {
  rp <- rank_primary[lineage == l, gene]
  ri <- rank_included[lineage == l, gene]
  data.table(
    lineage = l, pathway = names(pathways), resource_genes = lengths(pathways),
    primary_rank_genes = vapply(pathways, function(g) sum(g %chin% rp), integer(1)),
    sensitivity_rank_genes = vapply(pathways, function(g) sum(g %chin% ri), integer(1))
  )
}))
setorder(mapping, lineage, pathway)
setorder(mapping_recomputed, lineage, pathway)
record("pathway_rank_intersections", identical(mapping, mapping_recomputed),
       "All 100 lineage-pathway intersections reproduced")

gsea <- fread(file.path(table_dir, "05_fgsea_primary_all_50_sets.csv"))
gsea_included <- fread(file.path(table_dir, "05_fgsea_sensitivity_score_genes_included.csv"))
gsea_reversed <- fread(file.path(table_dir, "05_fgsea_sensitivity_reversed_rank.csv"))
record("fgsea_complete_50_sets_per_lineage",
       nrow(gsea) == 100L && gsea[, all(.N == 50L), by = lineage]$V1 |> all(),
       paste(gsea[, .N, by = lineage]$lineage, gsea[, .N, by = lineage]$N, collapse = "; "))

bh_tolerance <- 5e-15
bh_within_max <- gsea[, max(abs(q_within_lineage - p.adjust(pval, "BH")), na.rm = TRUE), by = lineage]
bh_global_max <- max(abs(gsea$q_global_100_tests - p.adjust(gsea$pval, "BH")), na.rm = TRUE)
record("BH_within_lineage_recomputed", all(bh_within_max$V1 <= bh_tolerance),
       paste(bh_within_max$lineage, format(bh_within_max$V1, scientific = TRUE), collapse = "; "))
record("BH_global_100_tests_recomputed", bh_global_max <= bh_tolerance,
       format(bh_global_max, scientific = TRUE))

leading_ok <- TRUE
for (i in seq_len(nrow(gsea))) {
  le <- if (is.na(gsea$leadingEdge_genes[i]) || !nzchar(gsea$leadingEdge_genes[i])) character() else {
    strsplit(gsea$leadingEdge_genes[i], ";", fixed = TRUE)[[1]]
  }
  rgenes <- rank_primary[lineage == gsea$lineage[i], gene]
  leading_ok <- leading_ok && length(le) == gsea$leading_edge_count[i] &&
    all(le %chin% pathways[[gsea$pathway[i]]]) && all(le %chin% rgenes) && !anyDuplicated(le)
}
record("leading_edge_membership_and_counts", leading_ok,
       "Every leading-edge gene belongs to its pathway and lineage rank")

reverse <- merge(
  gsea[, .(lineage, pathway, NES_primary = NES)],
  gsea_reversed[, .(lineage, pathway, NES_reversed = NES)],
  by = c("lineage", "pathway")
)
reverse_summary <- reverse[, .(
  all_signs_reversed = all(sign(NES_reversed) == -sign(NES_primary)),
  correlation = cor(NES_primary, NES_reversed),
  maximum_absolute_error = max(abs(NES_reversed + NES_primary))
), by = lineage]
record("reversed_rank_NES_direction",
       all(reverse_summary$all_signs_reversed) && all(reverse_summary$correlation < -0.999) &&
         all(reverse_summary$maximum_absolute_error < 0.05),
       paste(reverse_summary$lineage, sprintf("r=%.6f,maxerr=%.4f", reverse_summary$correlation,
                                              reverse_summary$maximum_absolute_error), collapse = "; "))

sensitivity_direction <- merge(
  gsea[, .(lineage, pathway, NES_primary = NES)],
  gsea_included[, .(lineage, pathway, NES_included = NES)],
  by = c("lineage", "pathway")
)
record("score_gene_inclusion_direction_sensitivity",
       all(sign(sensitivity_direction$NES_primary) == sign(sensitivity_direction$NES_included)),
       paste(sum(sign(sensitivity_direction$NES_primary) == sign(sensitivity_direction$NES_included)),
             "of", nrow(sensitivity_direction), "directions agree"))

expected_plot <- copy(gsea[is.finite(q_within_lineage) & q_within_lineage < 0.05])
expected_plot[, direction := fifelse(NES > 0, "Q4 (high ratio)", "Q1 (low ratio)")]
expected_plot[, abs_NES__ := abs(NES)]
setorder(expected_plot, lineage, direction, q_within_lineage, -abs_NES__, pathway)
expected_plot <- expected_plot[, head(.SD, 10L), by = .(lineage, direction)]
expected_plot[, abs_NES__ := NULL]
expected_plot[, pathway_label := tools::toTitleCase(tolower(gsub("_", " ", sub("^HALLMARK_", "", pathway))))]
observed_plot <- fread(file.path(table_dir, "plotdata_07_significant_hallmark_terms.csv"))
setorder(expected_plot, lineage, direction, q_within_lineage, pathway)
setorder(observed_plot, lineage, direction, q_within_lineage, pathway)
record("display_filter_and_top10_selection", identical(expected_plot, observed_plot),
       paste("display rows", nrow(observed_plot)))

display_audit <- fread(file.path(table_dir, "plotdata_07_display_audit.csv"))
record("CAF_null_direction_disclosed",
       display_audit[lineage == "CAF", all(significant_pathways == 0L & displayed_pathways == 0L)],
       "No primary CAF Hallmark passed q < 0.05 in either direction")

png_file <- file.path(pub_dir, "Fig_Hallmark_GSEA_Epithelial_CAF_ratio_Q4_vs_Q1.png")
img <- readPNG(png_file, native = TRUE)
record("publication_PNG_exact_dimensions", identical(dim(img), c(3120L, 4320L)),
       paste(dim(img), collapse = "x"))
file_info <- system2("file", shQuote(png_file), stdout = TRUE)
record("publication_PNG_RGB_noninterlaced",
       grepl("8-bit/color RGB, non-interlaced", file_info, fixed = TRUE), file_info)
sips_info <- system2("sips", c("-g", "pixelWidth", "-g", "pixelHeight", "-g", "dpiWidth",
                                "-g", "dpiHeight", "-g", "profile", shQuote(png_file)), stdout = TRUE)
record("publication_PNG_600dpi_sRGB",
       any(grepl("dpiWidth: 600.000", sips_info, fixed = TRUE)) &&
         any(grepl("dpiHeight: 600.000", sips_info, fixed = TRUE)) &&
         any(grepl("profile: sRGB IEC61966-2.1", sips_info, fixed = TRUE)),
       paste(sips_info, collapse = " | "))
xattr_info <- system2("xattr", c("-l", shQuote(png_file)), stdout = TRUE, stderr = TRUE)
flags_info <- system2("ls", c("-lO", shQuote(png_file)), stdout = TRUE)
record("publication_PNG_clean_filesystem_metadata",
       length(xattr_info) == 0L && !grepl("hidden", flags_info, fixed = TRUE),
       if (length(xattr_info)) paste(xattr_info, collapse = "; ") else "no xattrs; no hidden flag")

fwrite(checks, file.path(table_dir, "09_independent_verification_checks.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "02_verification_sessionInfo.txt"))
if (any(!checks$passed)) {
  print(checks[passed == FALSE])
  stop(sum(!checks$passed), " independent verification checks failed")
}
message("ALL INDEPENDENT VERIFICATION CHECKS PASSED: ", nrow(checks), "/", nrow(checks))
print(checks)
