#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
base_dir <- normalizePath(file.path(analysis_dir, ".."), mustWork = TRUE)
analysis_root <- normalizePath(file.path(base_dir, ".."), mustWork = TRUE)
prior_dir <- file.path(analysis_root, "02_lineage_reannotation")

.libPaths(c(file.path(base_dir, "R_libs"), file.path(prior_dir, "R_libs"), .libPaths()))
suppressPackageStartupMessages(library(data.table))

table_dir <- file.path(analysis_dir, "results", "tables")
log_dir <- file.path(analysis_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) {
  sub(" .*", "", system2("shasum", c("-a", "256", shQuote(path)), stdout = TRUE))
}
check <- function(name, condition, detail = "") {
  data.table(check = name, status = if (isTRUE(condition)) "PASS" else "FAIL", detail = detail)
}
near_equal <- function(x, y, tolerance = 1e-12) {
  length(x) == length(y) && all(is.na(x) == is.na(y)) &&
    all(abs(x[!is.na(x)] - y[!is.na(y)]) <= tolerance)
}

score_file <- file.path(base_dir, "intermediate", "01_cell_paper_style_scores.rds")
counts_file <- file.path(prior_dir, "intermediate", "01_full_counts_mt20_qc.rds")
lineage_file <- file.path(prior_dir, "results", "tables", "09_final_cell_lineages_FROZEN.csv")
expected_hashes <- c(
  counts = "7552ba882c306fa68380390a8172e0fc96af2c7be7cb3f5778345204e9fd61b2",
  lineages = "97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1",
  scores = "44a1c56f947910c2930ec0f1ce1e27a3b36cf059facd1835a7ae29ea05068262"
)
observed_hashes <- c(
  counts = sha256(counts_file),
  lineages = sha256(lineage_file),
  scores = sha256(score_file)
)

score <- as.data.table(readRDS(score_file))
membership <- fread(file.path(table_dir, "01_selected_cell_quartile_membership.csv.gz"))
thresholds <- fread(file.path(table_dir, "01_quartile_thresholds_and_counts.csv"))
deg <- fread(file.path(table_dir, "02_lineage_quartile_FindMarkers_all.csv"))
audit <- fread(file.path(table_dir, "02_DEG_enrichment_input_audit.csv"))
go_all <- fread(file.path(table_dir, "03_lineage_quartile_GO_BP_all_terms.csv"))
go_long <- fread(file.path(table_dir, "03_lineage_quartile_GO_BP_directional_all.csv"))
plot_data <- fread(file.path(table_dir, "plotdata_top_GO_BP_terms.csv"))

target <- copy(score[final_lineage %chin% c("Epithelial", "CAF") &
                       is.finite(production_degradation_ratio)])
target[, `:=`(
  q1_recomputed = as.numeric(quantile(production_degradation_ratio, 0.25, type = 7)),
  q3_recomputed = as.numeric(quantile(production_degradation_ratio, 0.75, type = 7))
), by = final_lineage]
target[, group_recomputed := fifelse(
  production_degradation_ratio <= q1_recomputed, "Bottom 25%",
  fifelse(production_degradation_ratio >= q3_recomputed, "Top 25%", "Middle 50%")
)]
selected_recomputed <- target[group_recomputed != "Middle 50%"]

thresholds_recomputed <- target[, .(
  finite_ratio_cells = .N,
  q1 = unique(q1_recomputed),
  q3 = unique(q3_recomputed),
  bottom_cells = sum(group_recomputed == "Bottom 25%"),
  middle_cells = sum(group_recomputed == "Middle 50%"),
  top_cells = sum(group_recomputed == "Top 25%"),
  q1_ties = sum(production_degradation_ratio == unique(q1_recomputed)),
  q3_ties = sum(production_degradation_ratio == unique(q3_recomputed)),
  total_patients = uniqueN(Patient),
  bottom_patients = uniqueN(Patient[group_recomputed == "Bottom 25%"]),
  top_patients = uniqueN(Patient[group_recomputed == "Top 25%"])
), by = final_lineage]
setorder(thresholds, final_lineage)
setorder(thresholds_recomputed, final_lineage)

membership_recomputed <- selected_recomputed[, .(
  cell_id, final_lineage, production_degradation_ratio,
  quartile_group = group_recomputed
)]
membership_observed <- membership[, .(
  cell_id, final_lineage, production_degradation_ratio,
  quartile_group = as.character(quartile_group)
)]
setkey(membership_recomputed, cell_id)
setkey(membership_observed, cell_id)

expected_direction <- fifelse(
  deg$p_val_adj < 0.05 & deg$avg_log2FC >= 0.25, "Top 25%",
  fifelse(deg$p_val_adj < 0.05 & deg$avg_log2FC <= -0.25, "Bottom 25%", "Not selected")
)

go_recomputed <- rbindlist(lapply(unique(go_all$lineage), function(lin) {
  z <- go_all[lineage == lin]
  rbindlist(list(
    z[, .(lineage, GO_ID, direction = "Top 25%", q_recomputed = p.adjust(P.Top, "BH"))],
    z[, .(lineage, GO_ID, direction = "Bottom 25%", q_recomputed = p.adjust(P.Bottom, "BH"))]
  ))
}))
setkey(go_recomputed, lineage, GO_ID, direction)
go_long_verify <- merge(
  go_long,
  go_recomputed,
  by = c("lineage", "GO_ID", "direction"),
  all.x = TRUE,
  sort = FALSE
)
go_long_verify[, `:=`(
  gene_fraction_recomputed = DE / N,
  background_fraction_recomputed = selected_entrez / universe_entrez,
  fold_enrichment_recomputed = (DE / N) / (selected_entrez / universe_entrez)
)]

eligible <- go_long[
  N >= 10 & N <= 500 & DE > 0 & is.finite(q) & q < 0.05 & is.finite(fold_enrichment)
]
setorder(eligible, lineage, direction, q, p, -fold_enrichment)
expected_plot <- eligible[, head(.SD, 10L), by = .(lineage, direction)]
setkey(expected_plot, lineage, direction, GO_ID)
setkey(plot_data, lineage, direction, GO_ID)

pngs <- list.files(file.path(analysis_dir, "results", "figures_publication"),
                   pattern = "\\.png$", full.names = TRUE)
sips_ok <- vapply(pngs, function(path) {
  out <- system2("sips", c("-g", "pixelWidth", "-g", "pixelHeight", "-g", "dpiWidth",
                           "-g", "dpiHeight", "-g", "profile", shQuote(path)),
                 stdout = TRUE, stderr = TRUE)
  any(grepl("dpiWidth: 600", out, fixed = TRUE)) &&
    any(grepl("dpiHeight: 600", out, fixed = TRUE)) &&
    any(grepl("profile: sRGB IEC61966-2.1", out, fixed = TRUE))
}, logical(1))

checks <- rbindlist(list(
  check("Frozen input SHA-256", identical(observed_hashes, expected_hashes),
        paste(names(observed_hashes), substr(observed_hashes, 1, 12), collapse = "; ")),
  check("Unique membership cell IDs", uniqueN(membership$cell_id) == nrow(membership),
        sprintf("n=%d", nrow(membership))),
  check("Recomputed membership cell-ID set", setequal(membership_recomputed$cell_id, membership_observed$cell_id)),
  check("Recomputed membership lineage", identical(membership_recomputed$final_lineage, membership_observed$final_lineage)),
  check("Recomputed membership ratio", near_equal(membership_recomputed$production_degradation_ratio,
                                                    membership_observed$production_degradation_ratio)),
  check("Recomputed membership quartile", identical(membership_recomputed$quartile_group,
                                                      membership_observed$quartile_group)),
  check("Recomputed quartile thresholds", near_equal(thresholds$q1, thresholds_recomputed$q1) &&
          near_equal(thresholds$q3, thresholds_recomputed$q3)),
  check("Recomputed quartile counts", identical(
    thresholds[, .(finite_ratio_cells, bottom_cells, middle_cells, top_cells, q1_ties, q3_ties,
                   total_patients, bottom_patients, top_patients)],
    thresholds_recomputed[, .(finite_ratio_cells, bottom_cells, middle_cells, top_cells, q1_ties, q3_ties,
                              total_patients, bottom_patients, top_patients)]
  )),
  check("No top/bottom cell overlap", !anyDuplicated(membership[, .(cell_id)])),
  check("DEG direction labels", identical(expected_direction, deg$direction)),
  check("BH q recomputation", near_equal(go_long_verify$q, go_long_verify$q_recomputed, tolerance = 1e-14)),
  check("GO gene fraction", near_equal(go_long_verify$gene_fraction,
                                        go_long_verify$gene_fraction_recomputed)),
  check("GO background fraction", near_equal(go_long_verify$background_fraction,
                                              go_long_verify$background_fraction_recomputed)),
  check("GO fold enrichment", near_equal(go_long_verify$fold_enrichment,
                                          go_long_verify$fold_enrichment_recomputed)),
  check("Displayed GO rows obey filters", all(plot_data$N >= 10 & plot_data$N <= 500 &
                                                plot_data$DE > 0 & plot_data$q < 0.05)),
  check("Displayed GO top-10 selection", identical(
    expected_plot[, .(lineage, direction, GO_ID)],
    plot_data[, .(lineage, direction, GO_ID)]
  )),
  check("Displayed GO max 10 per group", all(plot_data[, .N, by = .(lineage, direction)]$N <= 10L)),
  check("Publication PNG existence", length(pngs) == 3L && all(file.exists(pngs)),
        paste(basename(pngs), collapse = "; ")),
  check("Publication PNG 600 dpi sRGB", length(sips_ok) == 3L && all(sips_ok)),
  check("No missing key table values", !anyNA(thresholds) && !anyNA(audit$value) &&
          !anyNA(plot_data[, .(lineage, direction, GO_ID, Term, N, DE, p, q, fold_enrichment)]))
))

if (any(checks$status != "PASS")) {
  fwrite(checks, file.path(table_dir, "04_verification_summary.csv"))
  stop("One or more independent verification checks failed")
}

fwrite(checks, file.path(table_dir, "04_verification_summary.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "02_verification_sessionInfo.txt"))
message("Independent verification complete: ", nrow(checks), "/", nrow(checks), " checks PASS")
