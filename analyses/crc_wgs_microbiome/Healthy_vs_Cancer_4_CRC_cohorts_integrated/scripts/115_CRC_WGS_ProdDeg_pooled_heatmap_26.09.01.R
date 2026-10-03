#!/usr/bin/env Rscript

# Pooled sample-level CRC WGS analysis of species associations with community
# sarcosine Production and Degradation functional-potential scores.
#
# Statistical design:
#   1. Pool the four analysis-selected CRC WGS cohorts at the Run-ID level.
#   2. Use the frozen sample-level Production_sum and Degradation_sum scores.
#   3. Build the union species matrix, filling a species absent from a cohort's
#      exported matrix with zero, matching integrated_analysis_pooled.R.
#   4. Retain species detected in >=10% of all pooled score-matched samples.
#   5. Compute pooled Spearman correlations and BH correction separately for
#      Production and Degradation.
#   6. Display the 15 strongest positive BH-significant associations per axis.
#      No significance symbols are drawn in the figure.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggtext)
  library(ragg)
  library(magick)
  library(png)
})

options(stringsAsFactors = FALSE, scipen = 999)

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- sub("^--file=", "", script_arg[grepl("^--file=", script_arg)])
if (length(script_file) != 1L) {
  stop("Run this analysis with Rscript so the script path can be resolved.")
}

integrated_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
crc_root <- normalizePath(file.path(integrated_root, ".."), mustWork = TRUE)
project_root <- normalizePath(file.path(integrated_root, "..", "..", ".."), mustWork = TRUE)
use_root <- file.path(crc_root, "사용데이터_모음")

score_file <- file.path(
  use_root,
  "Derived_Sarcosine_Functional_Scores",
  "CRC_WGS_sarcosine_functional_scores_all_cohorts.csv"
)

cohorts <- c("PRJEB10878", "PRJEB27928", "PRJEB6070", "PRJNA429097")
species_files <- setNames(
  file.path(
    use_root,
    cohorts,
    paste0(cohorts, "_species_relative_abundance_matrix.csv")
  ),
  cohorts
)

cancer_file <- file.path(
  integrated_root,
  "results_integrated",
  "cross_cohort",
  "cross_cohort_CRC_enriched_membership_matrix.csv"
)
healthy_file <- file.path(
  integrated_root,
  "results_integrated",
  "SupFig3de_HEALTHY_ENRICHED_26.08.27",
  "healthy_enriched_membership_matrix.csv"
)

out_dir <- file.path(
  integrated_root,
  "results_integrated",
  "CRC_WGS_PROD_DEG_POOLED_HEATMAP_26.09.01"
)
use_out_dir <- file.path(use_root, "CRC_WGS_ProdDeg_PooledHeatmap_26.09.01")
manuscript_candidates <- list.dirs(project_root, recursive = FALSE, full.names = TRUE)
manuscript_candidates <- manuscript_candidates[
  grepl("^Manuscript", basename(manuscript_candidates))
]
if (length(manuscript_candidates) != 1L) {
  stop("Could not resolve the unique Manuscript work directory.")
}
transport_dir <- file.path(manuscript_candidates, "FigDesign&Manuscript")
transport_file <- file.path(
  transport_dir,
  "PPT_INSERT_CRC_WGS_ProdDeg_POOLED_Heatmap.png"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(use_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(transport_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(score_file, unname(species_files), cancer_file, healthy_file)
if (!all(file.exists(required_files))) {
  stop(
    "Missing required input(s): ",
    paste(required_files[!file.exists(required_files)], collapse = "; ")
  )
}

scores <- read.csv(score_file, check.names = FALSE)
required_score_cols <- c(
  "Run_ID", "Cohort", "Analysis_Group", "Production_sum", "Degradation_sum"
)
if (!all(required_score_cols %in% colnames(scores))) {
  stop("Score table is missing required columns.")
}
scores <- scores |>
  select(all_of(required_score_cols)) |>
  arrange(match(Cohort, cohorts), Run_ID)

stopifnot(
  nrow(scores) == 1647L,
  !anyDuplicated(scores$Run_ID),
  identical(sort(unique(scores$Cohort)), sort(cohorts)),
  all(is.finite(scores$Production_sum)),
  all(is.finite(scores$Degradation_sum)),
  all(scores$Production_sum >= 0),
  all(scores$Degradation_sum >= 0)
)

species_tables <- list()
species_union <- character()

for (cohort in cohorts) {
  dat <- read.csv(species_files[[cohort]], check.names = FALSE)
  if (!("Run_ID" %in% colnames(dat))) stop("Run_ID missing from ", species_files[[cohort]])
  if (anyDuplicated(dat$Run_ID)) stop("Duplicate Run_ID in ", species_files[[cohort]])
  if (anyDuplicated(colnames(dat)[-1])) stop("Duplicate species columns in ", species_files[[cohort]])

  cohort_score_ids <- scores$Run_ID[scores$Cohort == cohort]
  missing_species_profiles <- setdiff(cohort_score_ids, dat$Run_ID)
  if (length(missing_species_profiles) > 0L) {
    stop(
      "Score-matched Run IDs missing from species matrix for ", cohort, ": ",
      paste(missing_species_profiles, collapse = ", ")
    )
  }

  dat <- dat[match(cohort_score_ids, dat$Run_ID), , drop = FALSE]
  if (!identical(dat$Run_ID, cohort_score_ids)) stop("Run-ID ordering failed for ", cohort)
  species_tables[[cohort]] <- dat
  species_union <- union(species_union, colnames(dat)[-1])
}

pooled_species <- matrix(
  0,
  nrow = nrow(scores),
  ncol = length(species_union),
  dimnames = list(scores$Run_ID, species_union)
)

for (cohort in cohorts) {
  dat <- species_tables[[cohort]]
  ids <- dat$Run_ID
  taxa <- colnames(dat)[-1]
  values <- as.matrix(dat[, taxa, drop = FALSE])
  storage.mode(values) <- "double"
  if (any(!is.finite(values)) || any(values < 0)) {
    stop("Non-finite or negative species abundance in ", cohort)
  }
  pooled_species[ids, taxa] <- values
}

if (!identical(rownames(pooled_species), scores$Run_ID)) {
  stop("Pooled species matrix and score table are not identically ordered.")
}

prevalence <- colMeans(pooled_species > 0)
prevalence_table <- tibble(
  Species_full = names(prevalence),
  Species = sub("^.*\\|s__", "", names(prevalence)),
  n_detected = colSums(pooled_species > 0),
  n_samples = nrow(pooled_species),
  prevalence = as.numeric(prevalence),
  pass_prevalence_10pct = prevalence >= 0.10
) |>
  arrange(desc(prevalence), Species)

kept_species <- prevalence_table |>
  filter(pass_prevalence_10pct) |>
  pull(Species_full)

if (length(kept_species) == 0L) stop("No species passed the pooled 10% prevalence filter.")

correlate_axis <- function(score_values, axis_name) {
  if (sd(score_values) == 0) stop(axis_name, " score has zero variance.")

  rows <- lapply(kept_species, function(species_full) {
    abundance <- pooled_species[, species_full]
    if (sd(abundance) == 0) return(NULL)
    test <- suppressWarnings(
      cor.test(abundance, score_values, method = "spearman", exact = FALSE)
    )
    tibble(
      Axis = axis_name,
      Species_full = species_full,
      Species = sub("^.*\\|s__", "", species_full),
      rho = unname(test$estimate),
      p_value = test$p.value,
      n_samples = length(score_values),
      n_detected = sum(abundance > 0),
      prevalence = mean(abundance > 0)
    )
  })

  bind_rows(rows) |>
    mutate(q_value = p.adjust(p_value, method = "BH")) |>
    arrange(desc(rho), q_value, Species)
}

cor_all <- bind_rows(
  correlate_axis(scores$Production_sum, "Production"),
  correlate_axis(scores$Degradation_sum, "Degradation")
)

stopifnot(
  nrow(cor_all) == 2L * length(kept_species),
  !anyDuplicated(cor_all[c("Axis", "Species_full")]),
  all(is.finite(cor_all$rho)),
  all(is.finite(cor_all$p_value)),
  all(is.finite(cor_all$q_value))
)

top_per_axis <- cor_all |>
  filter(rho > 0, q_value < 0.05) |>
  group_by(Axis) |>
  arrange(desc(rho), q_value, Species, .by_group = TRUE) |>
  slice_head(n = 15L) |>
  mutate(axis_rank = row_number()) |>
  ungroup()

selection_counts <- top_per_axis |>
  count(Axis, name = "n_selected")
if (!identical(sort(selection_counts$n_selected), c(15L, 15L))) {
  stop("Fewer than 15 positive BH-significant species were available for one or both axes.")
}

selected_species <- top_per_axis |>
  group_by(Species_full, Species) |>
  arrange(desc(rho), axis_rank, .by_group = TRUE) |>
  summarise(
    selected_in_axes = paste(sort(unique(Axis)), collapse = ";"),
    primary_axis = first(Axis),
    primary_axis_rank = first(axis_rank),
    primary_axis_rho = first(rho),
    primary_axis_q = first(q_value),
    .groups = "drop"
  ) |>
  arrange(factor(primary_axis, levels = c("Production", "Degradation")), primary_axis_rank, Species)

plot_long <- selected_species |>
  select(
    Species_full, Species, selected_in_axes, primary_axis,
    primary_axis_rank
  ) |>
  crossing(Axis = c("Production", "Degradation")) |>
  left_join(
    cor_all |>
      select(
        Axis, Species_full, rho, p_value, q_value,
        n_samples, n_detected, prevalence
      ),
    by = c("Axis", "Species_full")
  )

cancer_membership <- read.csv(cancer_file, check.names = FALSE) |>
  filter(n_cohorts_significant >= 3L) |>
  transmute(Species_full, Cancer_recurrent_n = n_cohorts_significant)

healthy_membership <- read.csv(healthy_file, check.names = FALSE) |>
  filter(n_cohorts_strict >= 3L) |>
  transmute(Species_full, Healthy_recurrent_n = n_cohorts_strict)

disease_overlap <- selected_species |>
  left_join(cancer_membership, by = "Species_full") |>
  left_join(healthy_membership, by = "Species_full") |>
  mutate(
    Cancer_recurrent_n = replace_na(Cancer_recurrent_n, 0L),
    Healthy_recurrent_n = replace_na(Healthy_recurrent_n, 0L),
    disease_enrichment = case_when(
      Cancer_recurrent_n >= 3L & Healthy_recurrent_n >= 3L ~ "Both (unexpected)",
      Cancer_recurrent_n >= 3L ~ "Cancer-enriched in >=3/4 cohorts",
      Healthy_recurrent_n >= 3L ~ "Healthy-enriched in >=3/4 cohorts",
      TRUE ~ "Not recurrently disease-enriched"
    )
  ) |>
  arrange(factor(primary_axis, levels = c("Production", "Degradation")), primary_axis_rank, Species)

if (any(disease_overlap$disease_enrichment == "Both (unexpected)")) {
  stop("A selected species was recurrently enriched in both directions.")
}

selected_abundance <- as.data.frame(
  pooled_species[, selected_species$Species_full, drop = FALSE],
  check.names = FALSE
)
short_abundance_names <- paste0("Species__", selected_species$Species)
if (anyDuplicated(short_abundance_names)) {
  stop("Selected species short names are not unique; cannot create concise sample table.")
}
colnames(selected_abundance) <- short_abundance_names

sample_table <- bind_cols(scores, selected_abundance)
taxonomy_map <- selected_species |>
  transmute(
    sample_table_column = paste0("Species__", Species),
    Species,
    Species_full,
    selected_in_axes,
    primary_axis,
    primary_axis_rank
  )

sample_counts <- scores |>
  count(Cohort, Analysis_Group, name = "n_samples") |>
  arrange(match(Cohort, cohorts), Analysis_Group)

# Self-contained pooled analysis inputs for downstream reuse and exact auditing.
# The all-union table reproduces prevalence filtering; the prevalence-filtered
# table is the exact wide matrix used for the 411-species correlation tests.
pooled_input_all_union <- bind_cols(
  scores,
  as.data.frame(pooled_species, check.names = FALSE)
)
pooled_input_prevalence10 <- bind_cols(
  scores,
  as.data.frame(pooled_species[, kept_species, drop = FALSE], check.names = FALSE)
)

stopifnot(
  nrow(pooled_input_all_union) == 1647L,
  ncol(pooled_input_all_union) == 5L + ncol(pooled_species),
  nrow(pooled_input_prevalence10) == 1647L,
  ncol(pooled_input_prevalence10) == 5L + length(kept_species),
  identical(pooled_input_all_union$Run_ID, scores$Run_ID),
  identical(pooled_input_prevalence10$Run_ID, scores$Run_ID)
)

prod_rows <- selected_species |>
  filter(primary_axis == "Production") |>
  arrange(desc(primary_axis_rank))
deg_rows <- selected_species |>
  filter(primary_axis == "Degradation") |>
  arrange(desc(primary_axis_rank))
row_order <- c(deg_rows$Species, prod_rows$Species)

label_map <- selected_species |>
  mutate(label = paste0("<i>", gsub("_", " ", Species), "</i>")) |>
  select(Species, label)

plot_long <- plot_long |>
  mutate(
    Axis = factor(Axis, levels = c("Production", "Degradation")),
    Species = factor(Species, levels = row_order)
  )

group_boundary <- nrow(deg_rows) + 0.5
fill_limit <- max(0.50, ceiling(max(abs(plot_long$rho)) * 10) / 10)

p <- ggplot(plot_long, aes(x = Axis, y = Species, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.75, width = 0.96, height = 0.96) +
  geom_hline(yintercept = group_boundary, linewidth = 0.65, color = "#414141") +
  scale_fill_gradient2(
    low = "#2B6CB0", mid = "#F7F7F7", high = "#C43C39",
    midpoint = 0, limits = c(-fill_limit, fill_limit),
    oob = scales::squish,
    name = expression("Spearman " * rho)
  ) +
  scale_y_discrete(
    labels = setNames(label_map$label, label_map$Species),
    expand = c(0, 0)
  ) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  labs(
    title = "Pooled CRC WGS sarcosine functional associations",
    subtitle = paste0(
      "Sample-level Spearman correlations across four pooled cohorts (n=",
      nrow(scores), ")\nTop 15 positive BH-significant associations per functional axis"
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "Species prevalence ≥10%  •  Production: K00315 + K08688\n",
      "Degradation: K00301 + K00302 + K00303 + K00305 + K00306"
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 13, base_family = "Arial") +
  theme(
    plot.title = element_text(
      size = 18, face = "bold", color = "#111111", margin = margin(b = 5)
    ),
    plot.subtitle = element_text(
      size = 11.5, color = "#4D4D4D", lineheight = 1.15, margin = margin(b = 15)
    ),
    plot.caption = element_text(
      size = 9.5, color = "#555555", hjust = 0, margin = margin(t = 12)
    ),
    axis.text.x = element_text(
      size = 13, face = "bold", color = "#111111", margin = margin(b = 8)
    ),
    axis.text.y = ggtext::element_markdown(
      size = 10.5, color = "#1A1A1A", margin = margin(r = 8)
    ),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    legend.key.height = grid::unit(3.2, "cm"),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 20, 15, 15)
  )

png_file <- file.path(
  out_dir,
  "Fig_CRC_WGS_species_Production_Degradation_POOLED_heatmap.png"
)
pdf_file <- file.path(
  out_dir,
  "Fig_CRC_WGS_species_Production_Degradation_POOLED_heatmap.pdf"
)

DPI <- 400
ragg::agg_png(
  filename = png_file,
  width = 9.5,
  height = 12.0,
  units = "in",
  res = DPI,
  background = "white",
  scaling = 1
)
print(p)
dev.off()

srgb_profile <- "/System/Library/ColorSync/Profiles/sRGB Profile.icc"
if (!file.exists(srgb_profile)) stop("System sRGB profile not found.")

standardize_png <- function(path, dpi) {
  img <- magick::image_read(path)
  img <- magick::image_background(img, "white", flatten = TRUE)
  img <- magick::image_convert(img, format = "png", colorspace = "sRGB", depth = 8)
  img <- magick::image_strip(img)
  magick::image_write(img, path = path, format = "png", depth = 8)
  before <- png::readPNG(path, native = FALSE, info = TRUE)

  tmp <- tempfile(fileext = ".png")
  sips_output <- system2(
    "/usr/bin/sips",
    args = c(
      "-s", "dpiWidth", as.character(dpi),
      "-s", "dpiHeight", as.character(dpi),
      "-e", shQuote(srgb_profile), shQuote(path), "--out", shQuote(tmp)
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  if (!file.exists(tmp) || file.info(tmp)$size <= 5000) {
    stop("sips failed: ", paste(sips_output, collapse = "\n"))
  }
  if (!file.copy(tmp, path, overwrite = TRUE)) stop("Could not install standardized PNG.")
  unlink(tmp)

  after <- png::readPNG(path, native = FALSE, info = TRUE)
  if (!identical(dim(before), dim(after)) || !identical(as.numeric(before), as.numeric(after))) {
    stop("Embedding sRGB/dpi metadata changed decoded PNG pixels.")
  }
}

standardize_png(png_file, DPI)

if (!file.copy(png_file, transport_file, overwrite = TRUE)) {
  stop("Could not create PowerPoint transport PNG.")
}

for (path in c(png_file, transport_file)) {
  invisible(system2("/usr/bin/xattr", args = c("-c", shQuote(path)), stdout = TRUE, stderr = TRUE))
  invisible(system2("/bin/chmod", args = c("644", shQuote(path)), stdout = TRUE, stderr = TRUE))
  invisible(system2("/usr/bin/chflags", args = c("nohidden", shQuote(path)), stdout = TRUE, stderr = TRUE))
}

if (!identical(unname(tools::md5sum(png_file)), unname(tools::md5sum(transport_file)))) {
  stop("Transport PNG is not byte-identical to analytical PNG.")
}

grDevices::cairo_pdf(pdf_file, width = 9.5, height = 12.0, family = "Arial")
print(p)
dev.off()

output_tables <- list(
  pooled_analysis_input_all_union_species = pooled_input_all_union,
  pooled_analysis_input_prevalence10_species = pooled_input_prevalence10,
  pooled_species_prevalence = prevalence_table,
  pooled_species_Production_Degradation_correlations_all = cor_all,
  selected_top15_per_axis = top_per_axis,
  heatmap_plot_values_long = plot_long,
  selected_species_disease_enrichment_overlap = disease_overlap,
  selected_species_taxonomy_map = taxonomy_map,
  pooled_sample_scores_selected_species = sample_table,
  sample_counts_by_cohort_group = sample_counts
)

for (nm in names(output_tables)) {
  analysis_path <- file.path(out_dir, paste0(nm, ".csv"))
  use_path <- file.path(use_out_dir, paste0(nm, ".csv"))
  write.csv(output_tables[[nm]], analysis_path, row.names = FALSE, na = "")
  if (!file.copy(analysis_path, use_path, overwrite = TRUE)) {
    stop("Could not copy CSV to use-data folder: ", use_path)
  }
  if (!identical(unname(tools::md5sum(analysis_path)), unname(tools::md5sum(use_path)))) {
    stop("Use-data CSV copy does not match analysis output: ", basename(use_path))
  }
}

input_manifest <- data.frame(
  file = required_files,
  md5 = unname(tools::md5sum(required_files)),
  stringsAsFactors = FALSE
)
write.csv(input_manifest, file.path(out_dir, "input_file_md5.csv"), row.names = FALSE)

validation_lines <- c(
  "CRC WGS pooled Production/Degradation heatmap validation",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Design: raw sample-level pooled Spearman analysis; no cohort-level meta-analysis used",
  paste0("Score-matched pooled samples: ", nrow(scores)),
  paste0("Cohorts: ", paste(cohorts, collapse = ", ")),
  paste0("Union species features: ", ncol(pooled_species)),
  paste0("Species passing >=10% pooled prevalence: ", length(kept_species)),
  paste0("Correlation rows: ", nrow(cor_all)),
  paste0(
    "Selected rows per axis: ",
    paste(selection_counts$Axis, selection_counts$n_selected, sep = "=", collapse = "; ")
  ),
  paste0("Unique selected species: ", nrow(selected_species)),
  paste0(
    "Recurrent Cancer-enriched overlaps (>=3/4): ",
    sum(disease_overlap$Cancer_recurrent_n >= 3L)
  ),
  paste0(
    "Recurrent Healthy-enriched overlaps (>=3/4): ",
    sum(disease_overlap$Healthy_recurrent_n >= 3L)
  ),
  paste0("PNG exists and non-empty: ", file.exists(png_file) && file.info(png_file)$size > 0),
  paste0(
    "PowerPoint transport copy byte-identical: ",
    identical(unname(tools::md5sum(png_file)), unname(tools::md5sum(transport_file)))
  ),
  "No significance symbols are drawn; P and BH q values are retained in CSV outputs."
)
writeLines(validation_lines, file.path(out_dir, "validation_report.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

readme <- c(
  "CRC WGS pooled species–sarcosine functional-potential heatmap",
  "",
  "Design",
  "- This is a raw sample-level pooled analysis, not a meta-analysis.",
  "- Four analysis-selected CRC WGS cohorts were concatenated by unique Run ID.",
  "- Species absent from a cohort-specific exported matrix were zero-filled, matching the frozen integrated workflow.",
  "- Species detected in >=10% of 1,647 score-matched samples were tested.",
  "- Spearman correlation was computed against Production_sum and Degradation_sum; BH correction was performed separately for each axis.",
  "- The figure shows the 15 strongest positive BH-significant associations per axis without significance glyphs.",
  "",
  "Interpretation guardrail",
  "- Scores are WGS community functional-potential abundances, not sarcosine concentration, enzyme activity, metabolic flux or causality.",
  "- Because samples are pooled directly, the displayed correlations are not adjusted for cohort.",
  "",
  "PowerPoint transport",
  paste0("- ", transport_file),
  "",
  "CSV delivery folder",
  paste0("- ", use_out_dir),
  "- pooled_analysis_input_all_union_species.csv: all pooled union species used to derive prevalence",
  "- pooled_analysis_input_prevalence10_species.csv: exact 411-species matrix used for pooled correlations",
  "",
  "Reproduction",
  paste0("Rscript ", script_file)
)
writeLines(readme, file.path(out_dir, "README.txt"))

cat("Pooled analysis complete.\n")
cat("Samples:", nrow(scores), "\n")
cat("Prevalence-filtered species:", length(kept_species), "\n")
cat("Selected unique species:", nrow(selected_species), "\n")
cat("Figure:", png_file, "\n")
cat("CSV delivery:", use_out_dir, "\n")
