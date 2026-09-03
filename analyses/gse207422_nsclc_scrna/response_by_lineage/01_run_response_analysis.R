#!/usr/bin/env Rscript

args_all <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args_all[grepl("^--file=", args_all)])
if (length(script_path) != 1L) stop("Could not resolve script path")
analysis_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
parent_dir <- normalizePath(file.path(analysis_dir, ".."), mustWork = TRUE)
analysis_root <- normalizePath(file.path(parent_dir, ".."), mustWork = TRUE)
prior_dir <- file.path(analysis_root, "02_lineage_reannotation")
.libPaths(c(file.path(parent_dir, "R_libs"), file.path(prior_dir, "R_libs"), .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

options(stringsAsFactors = FALSE, warn = 1)
MASTER_SEED <- 260825L
set.seed(MASTER_SEED)

table_dir <- file.path(analysis_dir, "results", "tables")
pub_dir <- file.path(analysis_dir, "results", "figures_publication")
log_dir <- file.path(analysis_dir, "logs")
for (d in c(table_dir, pub_dir, log_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

sha256 <- function(path) sub(" .*", "", system2("shasum", c("-a", "256", path), stdout = TRUE))
plan_file <- file.path(analysis_dir, "ANALYSIS_PLAN_FROZEN.md")
expected_plan_sha <- "a20203c72f28c344530f1a729885e37c51c32cdb00dd694157ac36c08b976bfe"
if (!file.exists(plan_file) || sha256(plan_file) != expected_plan_sha) stop("Frozen plan hash mismatch")

score_file <- file.path(parent_dir, "results", "tables", "01_cell_paper_style_scores.csv.gz")
lineage_file <- file.path(prior_dir, "results", "tables", "09_final_cell_lineages_FROZEN.csv")
expected_input_sha <- c(
  scores = "17724b2b1fcbad436e6ce4cee91a34bda851361aac7b14c26a3557e5d006fb61",
  lineages = "97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1"
)
observed_input_sha <- c(scores = sha256(score_file), lineages = sha256(lineage_file))
if (!identical(observed_input_sha, expected_input_sha)) stop("Verified input SHA-256 mismatch")

lineage_order <- c(
  "Epithelial", "CAF", "B cell", "Plasma cell", "CD4 T cell",
  "CD8 T cell", "Cycling T cell", "NK cell", "Mast cell",
  "Neutrophil", "Monocyte", "Macrophage", "Conventional DC", "pDC"
)
score_order <- c("Production", "Degradation", "Production/degradation ratio")
focus_lineages <- c("Epithelial", "CAF", "CD8 T cell", "NK cell")
COL_MPR <- "#2E5F8A"
COL_NMPR <- "#C47B3B"
COL_TEXT <- "#202020"
COL_GREY <- "#BDBDBD"
group_cols <- c("MPR/pCR" = COL_MPR, "NMPR" = COL_NMPR)

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
      strip.background = element_blank(),
      strip.text = element_text(size = 6.7, face = "bold", colour = COL_TEXT),
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6),
      plot.margin = margin(3, 3, 3, 3, unit = "pt")
    )
}

save_plot <- function(stem, plot, width, height) {
  png_file <- file.path(pub_dir, paste0(stem, ".png"))
  pdf_file <- file.path(pub_dir, paste0(stem, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, units = "in", dpi = 600,
         device = "png", bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, units = "in",
         device = cairo_pdf, bg = "white")
  if (!all(file.exists(c(png_file, pdf_file)))) stop("Figure write failed: ", stem)
}

q1 <- function(x) as.numeric(quantile(x, 0.25, na.rm = TRUE, type = 7))
q3 <- function(x) as.numeric(quantile(x, 0.75, na.rm = TRUE, type = 7))

exact_rank_sum <- function(x_mpr, x_nmpr) {
  x <- c(x_mpr, x_nmpr)
  n1 <- length(x_mpr)
  n <- length(x)
  if (n1 < 1L || n1 >= n || any(!is.finite(x))) return(NA_real_)
  ranks <- rank(x, ties.method = "average")
  observed <- sum(ranks[seq_len(n1)])
  expected <- n1 * (n + 1) / 2
  allocations <- combn(n, n1)
  permuted <- colSums(matrix(ranks[allocations], nrow = n1))
  tolerance <- 1e-12
  mean(abs(permuted - expected) >= abs(observed - expected) - tolerance)
}

hl_shift <- function(x_mpr, x_nmpr) {
  median(as.vector(outer(x_mpr, x_nmpr, "-")))
}

bootstrap_hl <- function(x_mpr, x_nmpr, seed, B = 10000L) {
  set.seed(seed)
  boot <- replicate(B, {
    a <- sample(x_mpr, length(x_mpr), replace = TRUE)
    b <- sample(x_nmpr, length(x_nmpr), replace = TRUE)
    hl_shift(a, b)
  })
  as.numeric(quantile(boot, c(0.025, 0.975), type = 7, names = FALSE))
}

analyze_long <- function(metric_long, do_bootstrap = FALSE, seed_offset = 0L) {
  metric_long[, {
    x_mpr <- value[response_group == "MPR/pCR" & is.finite(value)]
    x_nmpr <- value[response_group == "NMPR" & is.finite(value)]
    n_mpr <- length(x_mpr)
    n_nmpr <- length(x_nmpr)
    tested <- n_mpr >= 3L && n_nmpr >= 3L
    effect <- if (tested) hl_shift(x_mpr, x_nmpr) else NA_real_
    ci <- c(NA_real_, NA_real_)
    if (tested && do_bootstrap) {
      lineage_index <- match(as.character(final_lineage[1]), lineage_order)
      score_index <- match(as.character(score[1]), score_order)
      ci <- bootstrap_hl(
        x_mpr, x_nmpr,
        seed = MASTER_SEED + seed_offset + 100L * score_index + lineage_index
      )
    }
    list(
      n_MPR_pCR = n_mpr,
      n_NMPR = n_nmpr,
      MPR_pCR_median = if (n_mpr) median(x_mpr) else NA_real_,
      MPR_pCR_q1 = if (n_mpr) q1(x_mpr) else NA_real_,
      MPR_pCR_q3 = if (n_mpr) q3(x_mpr) else NA_real_,
      NMPR_median = if (n_nmpr) median(x_nmpr) else NA_real_,
      NMPR_q1 = if (n_nmpr) q1(x_nmpr) else NA_real_,
      NMPR_q3 = if (n_nmpr) q3(x_nmpr) else NA_real_,
      HL_shift_MPR_minus_NMPR = effect,
      HL_bootstrap_CI_low = ci[1],
      HL_bootstrap_CI_high = ci[2],
      exact_P = if (tested) exact_rank_sum(x_mpr, x_nmpr) else NA_real_,
      tested = tested
    )
  }, by = .(final_lineage, score)]
}

message("Reading frozen score values only after plan freeze")
scores <- fread(score_file)
lin_check <- fread(
  lineage_file,
  select = c("cell_id", "Patient", "Sample", "Resource", "Pathologic.Response", "final_lineage", "analysis_eligible")
)
if (nrow(scores) != 92053L || nrow(lin_check) != 92053L) stop("Unexpected input rows")
if (anyDuplicated(scores$cell_id) || anyDuplicated(lin_check$cell_id)) stop("Duplicate cell IDs")
setkey(scores, cell_id)
setkey(lin_check, cell_id)
if (!identical(scores$cell_id, lin_check$cell_id)) stop("Cell-ID order mismatch after keying")
for (nm in c("Patient", "Sample", "Resource", "Pathologic.Response", "final_lineage", "analysis_eligible")) {
  if (!identical(scores[[nm]], lin_check[[nm]])) stop("Metadata mismatch: ", nm)
}
required_values <- c(
  "production_module_shifted", "degradation_module_shifted",
  "ratio_defined", "production_degradation_ratio"
)
if (!all(required_values %in% names(scores))) stop("Required score column missing")
if (any(!is.finite(scores$production_module_shifted)) ||
    any(!is.finite(scores$degradation_module_shifted)) ||
    any(scores$production_module_shifted < 0) ||
    any(scores$degradation_module_shifted < 0)) {
  stop("Invalid shifted module score")
}
if (sum(!scores$ratio_defined) != 1L || sum(!is.finite(scores$production_degradation_ratio)) != 1L) {
  stop("Unexpected ratio-defined structure")
}

cells <- scores[
  analysis_eligible == TRUE &
    Resource == "Post-treatment surgery" &
    Pathologic.Response %chin% c("MPR", "pCR", "NMPR")
]
cells[, response_group := fifelse(Pathologic.Response %chin% c("MPR", "pCR"), "MPR/pCR", "NMPR")]
cells[, response_group := factor(response_group, levels = c("MPR/pCR", "NMPR"))]
cells[, final_lineage := factor(final_lineage, levels = lineage_order)]
expected_patients <- c("P03", "P06", "P11", "P14", "P02", "P04", "P07", "P09", "P10", "P12", "P13", "P15")
if (!setequal(unique(cells$Patient), expected_patients)) stop("Primary patient set mismatch")

patient_summary <- cells[, {
  ratio_values <- production_degradation_ratio[ratio_defined & is.finite(production_degradation_ratio)]
  list(
    n_cells = .N,
    n_finite_ratio_cells = length(ratio_values),
    production_mean = mean(production_module_shifted),
    production_median = median(production_module_shifted),
    production_q1 = q1(production_module_shifted),
    production_q3 = q3(production_module_shifted),
    degradation_mean = mean(degradation_module_shifted),
    degradation_median = median(degradation_module_shifted),
    degradation_q1 = q1(degradation_module_shifted),
    degradation_q3 = q3(degradation_module_shifted),
    ratio_median = if (length(ratio_values)) median(ratio_values) else NA_real_,
    ratio_mean = if (length(ratio_values)) mean(ratio_values) else NA_real_,
    ratio_q1 = if (length(ratio_values)) q1(ratio_values) else NA_real_,
    ratio_q3 = if (length(ratio_values)) q3(ratio_values) else NA_real_,
    ratio_of_means = if (mean(degradation_module_shifted) > 0) {
      mean(production_module_shifted) / mean(degradation_module_shifted)
    } else NA_real_
  )
}, by = .(Patient, Sample, Pathologic.Response, response_group, final_lineage)]
setorder(patient_summary, final_lineage, response_group, Patient)

make_primary_long <- function(summary_dt, min_cells = 10L) {
  s <- summary_dt[n_cells >= min_cells]
  out <- rbindlist(list(
    s[, .(Patient, response_group, final_lineage, n_cells, score = "Production", value = production_mean)],
    s[, .(Patient, response_group, final_lineage, n_cells, score = "Degradation", value = degradation_mean)],
    s[, .(Patient, response_group, final_lineage, n_cells, score = "Production/degradation ratio", value = ratio_median)]
  ))
  out[, score := factor(score, levels = score_order)]
  out[, final_lineage := factor(final_lineage, levels = lineage_order)]
  setorder(out, score, final_lineage, response_group, Patient)
  out
}

primary_long <- make_primary_long(patient_summary, 10L)
primary_results <- analyze_long(primary_long, do_bootstrap = TRUE)
primary_results[tested == TRUE, BH_q := p.adjust(exact_P, method = "BH")]
primary_results[tested == FALSE, BH_q := NA_real_]
primary_results[, family_tests := sum(tested)]
if (unique(primary_results$family_tests) != 39L) stop("Expected 39 primary tests")
setorder(primary_results, score, final_lineage)

threshold_results <- rbindlist(lapply(c(1L, 10L, 50L, 100L), function(min_cells) {
  d <- make_primary_long(patient_summary, min_cells)
  z <- analyze_long(d, do_bootstrap = FALSE)
  z[tested == TRUE, BH_q := p.adjust(exact_P, method = "BH")]
  z[tested == FALSE, BH_q := NA_real_]
  z[, `:=`(min_cells = min_cells, family_tests = sum(tested))]
  z
}), use.names = TRUE)
setorder(threshold_results, min_cells, score, final_lineage)

aggregation_long <- rbindlist(list(
  patient_summary[n_cells >= 10L, .(
    Patient, response_group, final_lineage, n_cells,
    score = "Production", aggregation = "Cell median", value = production_median
  )],
  patient_summary[n_cells >= 10L, .(
    Patient, response_group, final_lineage, n_cells,
    score = "Degradation", aggregation = "Cell median", value = degradation_median
  )],
  patient_summary[n_cells >= 10L, .(
    Patient, response_group, final_lineage, n_cells,
    score = "Production/degradation ratio", aggregation = "Ratio of cell-score means", value = ratio_of_means
  )]
))
aggregation_long[, score := factor(score, levels = score_order)]
aggregation_results <- aggregation_long[, {
  z <- analyze_long(.SD, do_bootstrap = FALSE)
  z
}, by = aggregation]
aggregation_results[tested == TRUE, BH_q := p.adjust(exact_P, method = "BH"), by = aggregation]
aggregation_results[tested == FALSE, BH_q := NA_real_]
aggregation_results[, family_tests := sum(tested), by = aggregation]
setorder(aggregation_results, aggregation, score, final_lineage)

loo <- primary_long[, {
  patients <- unique(Patient)
  rbindlist(lapply(patients, function(drop_patient) {
    d <- .SD[Patient != drop_patient]
    x_mpr <- d[response_group == "MPR/pCR", value]
    x_nmpr <- d[response_group == "NMPR", value]
    data.table(
      dropped_patient = drop_patient,
      dropped_group = as.character(.SD[Patient == drop_patient, response_group][1]),
      n_MPR_pCR = length(x_mpr),
      n_NMPR = length(x_nmpr),
      HL_shift_MPR_minus_NMPR = if (length(x_mpr) >= 2L && length(x_nmpr) >= 2L) {
        hl_shift(x_mpr, x_nmpr)
      } else NA_real_
    )
  }))
}, by = .(final_lineage, score)]
setorder(loo, score, final_lineage, dropped_patient)

fwrite(patient_summary, file.path(table_dir, "01_patient_lineage_score_summaries.csv"))
fwrite(primary_long, file.path(table_dir, "01_primary_plot_values.csv"))
fwrite(primary_results, file.path(table_dir, "02_primary_exact_results.csv"))
fwrite(threshold_results, file.path(table_dir, "03_cell_threshold_sensitivity.csv"))
fwrite(aggregation_results, file.path(table_dir, "04_aggregation_sensitivity.csv"))
fwrite(loo, file.path(table_dir, "05_leave_one_patient_out_HL.csv"))

format_q <- function(q, tested) {
  out <- rep("NT", length(q))
  ok <- tested & is.finite(q)
  out[ok & q < 0.001] <- "q<0.001"
  out[ok & q >= 0.001] <- sprintf("q=%.3f", q[ok & q >= 0.001])
  out
}
primary_results[, q_label := format_q(BH_q, tested)]

score_titles <- c(
  "Production" = "Production module score",
  "Degradation" = "Degradation module score",
  "Production/degradation ratio" = "Production/degradation module-score ratio"
)

make_score_plot <- function(score_name, show_legend = TRUE) {
  d <- primary_long[score == score_name]
  r <- primary_results[score == score_name]
  y_rng <- range(d$value, na.rm = TRUE)
  y_span <- diff(y_rng)
  if (!is.finite(y_span) || y_span == 0) y_span <- max(abs(y_rng), 1)
  r[, annotation_y := y_rng[2] + 0.10 * y_span]
  ggplot(d, aes(final_lineage, value, colour = response_group, fill = response_group)) +
    geom_boxplot(
      aes(group = interaction(final_lineage, response_group)),
      position = position_dodge(width = 0.72), width = 0.58,
      outlier.shape = NA, alpha = 0.16, linewidth = 0.30
    ) +
    geom_point(
      position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.72, seed = MASTER_SEED),
      size = 1.20, alpha = 0.90, stroke = 0
    ) +
    geom_text(
      data = r,
      aes(x = final_lineage, y = annotation_y, label = q_label),
      inherit.aes = FALSE, size = 1.65, colour = COL_TEXT, vjust = 0
    ) +
    scale_colour_manual(values = group_cols, drop = FALSE) +
    scale_fill_manual(values = group_cols, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.04, 0.18))) +
    labs(
      title = unname(score_titles[score_name]),
      subtitle = "Post-treatment surgery; each point is one patient; exact rank-sum q across 39 primary tests",
      x = NULL, y = unname(score_titles[score_name]),
      colour = "Pathologic response", fill = "Pathologic response"
    ) +
    theme_nc() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.35),
      legend.position = if (show_legend) "top" else "none",
      legend.justification = "left"
    )
}

p_prod <- make_score_plot("Production", TRUE)
p_deg <- make_score_plot("Degradation", FALSE)
p_ratio <- make_score_plot("Production/degradation ratio", FALSE)
save_plot("Fig_response_A_production_by_lineage", p_prod, 7.20, 3.25)
save_plot("Fig_response_B_degradation_by_lineage", p_deg, 7.20, 3.25)
save_plot("Fig_response_C_ratio_by_lineage", p_ratio, 7.20, 3.25)

p_full <- p_prod / p_deg / p_ratio +
  plot_layout(heights = c(1, 0.95, 0.95)) +
  plot_annotation(
    title = "Sarcosine module scores by pathologic response across tumour-microenvironment lineages",
    subtitle = "MPR/pCR versus NMPR after neoadjuvant anti-PD-1 plus chemotherapy; patient is the independent unit",
    theme = theme(
      plot.title = element_text(family = "Arial", size = 8, face = "bold", colour = COL_TEXT),
      plot.subtitle = element_text(family = "Arial", size = 6.5, colour = "#555555")
    )
  )
save_plot("Fig_response_D_three_scores_all_lineages", p_full, 7.20, 8.75)

focus_values <- primary_long[final_lineage %chin% focus_lineages]
focus_values[, final_lineage := factor(final_lineage, levels = focus_lineages)]
focus_results <- primary_results[final_lineage %chin% focus_lineages]
focus_results[, final_lineage := factor(final_lineage, levels = focus_lineages)]
focus_y <- focus_values[, {
  rr <- range(value, na.rm = TRUE)
  span <- diff(rr)
  if (!is.finite(span) || span == 0) span <- max(abs(rr), 1)
  .(annotation_y = rr[2] + 0.11 * span)
}, by = score]
focus_results <- merge(focus_results, focus_y, by = "score", all.x = TRUE)
focus_results[, stat_label := fifelse(
  tested,
  paste0("P=", formatC(exact_P, format = "g", digits = 2), "\nq=", formatC(BH_q, format = "g", digits = 2)),
  "NT"
)]

p_focus <- ggplot(focus_values, aes(response_group, value, colour = response_group, fill = response_group)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, alpha = 0.16, linewidth = 0.30) +
  geom_point(position = position_jitter(width = 0.08, seed = MASTER_SEED), size = 1.35, alpha = 0.92) +
  geom_text(
    data = focus_results,
    aes(x = 1.5, y = annotation_y, label = stat_label),
    inherit.aes = FALSE, size = 1.70, colour = COL_TEXT, lineheight = 0.90
  ) +
  facet_grid(score ~ final_lineage, scales = "free_y", switch = "y") +
  scale_colour_manual(values = group_cols, drop = FALSE) +
  scale_fill_manual(values = group_cols, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.20))) +
  labs(
    title = "Lineage-resolved sarcosine module scores by pathologic response",
    subtitle = "Post-treatment surgery; points are patients; exact P and 39-test BH q",
    x = NULL, y = NULL, colour = "Pathologic response", fill = "Pathologic response"
  ) +
  theme_nc() +
  theme(
    legend.position = "top", legend.justification = "left",
    axis.text.x = element_text(angle = 25, hjust = 1, size = 5.6),
    strip.placement = "outside", panel.spacing = unit(2.3, "mm")
  )
save_plot("Fig_response_E_focus_epithelial_CAF_CD8_NK", p_focus, 7.20, 6.30)

matrix_data <- copy(primary_results)
matrix_data[, score := factor(score, levels = score_order)]
matrix_data[, final_lineage := factor(final_lineage, levels = rev(lineage_order))]
matrix_data[, scale_denominator := {
  z <- IQR(HL_shift_MPR_minus_NMPR[tested], na.rm = TRUE)
  if (!is.finite(z) || z == 0) z <- sd(HL_shift_MPR_minus_NMPR[tested], na.rm = TRUE)
  if (!is.finite(z) || z == 0) z <- 1
  z
}, by = score]
matrix_data[, standardized_HL := HL_shift_MPR_minus_NMPR / scale_denominator]
matrix_data[, display_HL := pmax(-2.5, pmin(2.5, standardized_HL))]
matrix_data[, size_value := fifelse(tested, pmin(4, -log10(pmax(BH_q, 1e-4))), NA_real_)]
size_max <- max(matrix_data$size_value, na.rm = TRUE)
size_breaks <- pretty(c(0, size_max), n = 4)
size_breaks <- size_breaks[size_breaks >= 0 & size_breaks <= size_max]

p_matrix <- ggplot(matrix_data[tested == TRUE], aes(score, final_lineage)) +
  geom_point(
    aes(size = size_value, fill = display_HL), shape = 21,
    colour = "#666666", stroke = 0.25
  ) +
  geom_point(
    data = matrix_data[tested == TRUE & BH_q < 0.05],
    aes(size = size_value, fill = display_HL), shape = 21,
    colour = COL_TEXT, stroke = 0.95
  ) +
  geom_point(
    data = matrix_data[tested == FALSE],
    shape = 4, size = 2.2, stroke = 0.65, colour = COL_GREY
  ) +
  scale_fill_gradient2(
    low = COL_NMPR, mid = "white", high = COL_MPR,
    midpoint = 0, limits = c(-2.5, 2.5),
    name = "MPR/pCR − NMPR\nrobust-scaled HL shift"
  ) +
  scale_size_continuous(
    range = c(1.5, 5.2), limits = c(0, size_max), breaks = size_breaks,
    name = expression(-log[10](q))
  ) +
  labs(
    title = "Pathologic-response association matrix",
    subtitle = "Black outline: q < 0.05; cross: insufficient patients; display scaling does not alter tests",
    x = NULL, y = NULL
  ) +
  theme_nc() +
  theme(
    axis.text.x = element_text(angle = 25, hjust = 1, size = 5.8),
    axis.text.y = element_text(size = 5.8),
    legend.position = "right", panel.grid = element_blank()
  )
save_plot("Fig_response_F_effect_matrix", p_matrix, 5.60, 4.60)

figure_spec <- data.table(
  figure = c(
    "Fig_response_A_production_by_lineage", "Fig_response_B_degradation_by_lineage",
    "Fig_response_C_ratio_by_lineage", "Fig_response_D_three_scores_all_lineages",
    "Fig_response_E_focus_epithelial_CAF_CD8_NK", "Fig_response_F_effect_matrix"
  ),
  width_in = c(7.20, 7.20, 7.20, 7.20, 7.20, 5.60),
  height_in = c(3.25, 3.25, 3.25, 8.75, 6.30, 4.60),
  dpi = 600L
)
figure_spec[, `:=`(
  expected_width_px = as.integer(round(width_in * dpi)),
  expected_height_px = as.integer(round(height_in * dpi))
)]
fwrite(primary_results, file.path(table_dir, "02_primary_exact_results.csv"))
fwrite(matrix_data, file.path(table_dir, "06_effect_matrix_plotdata.csv"))
fwrite(figure_spec, file.path(table_dir, "07_figure_specifications.csv"))
writeLines(capture.output(sessionInfo()), file.path(log_dir, "01_analysis_sessionInfo.txt"))
message("Frozen response-by-lineage analysis completed")
