#!/usr/bin/env Rscript

# Four-cohort CRC WGS integration of species associations with community-level
# sarcosine production and degradation functional potentials.
#
# Frozen inputs:
#   results_integrated/cross_cohort/cross_cohort_prod_assoc_correlations.csv
#   results_integrated/cross_cohort/cross_cohort_deg_assoc_correlations.csv
#
# Statistical design:
#   1. The frozen files contain within-cohort Spearman correlations.
#   2. For each species and functional axis, correlations observed in >=3 cohorts
#      are integrated with a random-effects model on Fisher-z transformed rho.
#   3. BH correction is applied across all meta-analysed species within each axis.
#   4. The display contains the 15 strongest positive associations per axis among
#      species with BH q < 0.05. No significance symbols are drawn in the figure.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(metafor)
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
root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

cross_dir <- file.path(root, "results_integrated", "cross_cohort")
prod_file <- file.path(cross_dir, "cross_cohort_prod_assoc_correlations.csv")
deg_file <- file.path(cross_dir, "cross_cohort_deg_assoc_correlations.csv")
cancer_file <- file.path(cross_dir, "cross_cohort_CRC_enriched_membership_matrix.csv")
healthy_file <- file.path(
  root,
  "results_integrated",
  "SupFig3de_HEALTHY_ENRICHED_26.08.27",
  "healthy_enriched_membership_matrix.csv"
)

out_dir <- file.path(
  root,
  "results_integrated",
  "CRC_WGS_PROD_DEG_META_HEATMAP_26.08.31"
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

project_root <- normalizePath(file.path(root, "..", "..", ".."), mustWork = TRUE)
transport_dir <- file.path(project_root, "Manuscript작업", "FigDesign&Manuscript")
dir.create(transport_dir, recursive = TRUE, showWarnings = FALSE)
transport_file <- file.path(transport_dir, "PPT_INSERT_CRC_WGS_ProdDeg_Heatmap.png")
srgb_profile <- "/System/Library/ColorSync/Profiles/sRGB Profile.icc"
if (!file.exists(srgb_profile)) stop("System sRGB profile not found: ", srgb_profile)

required_files <- c(prod_file, deg_file, cancer_file, healthy_file)
if (!all(file.exists(required_files))) {
  stop("Missing required frozen input(s): ", paste(required_files[!file.exists(required_files)], collapse = "; "))
}

prod <- read.csv(prod_file, check.names = FALSE) |>
  mutate(Axis = "Production", KO_set = prod_kos_summed) |>
  select(Cohort, Species_full, Species, rho, p_value, p_adj, KO_set, n_samples, Axis)

deg <- read.csv(deg_file, check.names = FALSE) |>
  mutate(Axis = "Degradation", KO_set = deg_kos_summed) |>
  select(Cohort, Species_full, Species, rho, p_value, p_adj, KO_set, n_samples, Axis)

cor_long <- bind_rows(prod, deg)
expected_cohorts <- c("PRJEB10878", "PRJEB27928", "PRJEB6070", "PRJNA429097")

stopifnot(
  identical(sort(unique(cor_long$Cohort)), sort(expected_cohorts)),
  all(is.finite(cor_long$rho)),
  all(cor_long$rho >= -1 & cor_long$rho <= 1),
  all(is.finite(cor_long$n_samples)),
  all(cor_long$n_samples > 3),
  nrow(cor_long) == nrow(distinct(cor_long, Axis, Cohort, Species_full))
)

sample_n_check <- cor_long |>
  distinct(Axis, Cohort, n_samples) |>
  count(Cohort, n_samples, name = "n_axes")
if (nrow(sample_n_check) != 4L || any(sample_n_check$n_axes != 2L)) {
  stop("Production and degradation inputs do not agree on cohort sample counts.")
}

meta_one <- function(dat) {
  dat <- dat |>
    filter(is.finite(rho), is.finite(n_samples), n_samples > 3) |>
    arrange(Cohort)

  if (nrow(dat) < 3L) {
    return(tibble(
      meta_z = NA_real_, meta_se_z = NA_real_, meta_rho = NA_real_,
      ci_low_rho = NA_real_, ci_high_rho = NA_real_, meta_p = NA_real_,
      tau2 = NA_real_, I2 = NA_real_, Q_p = NA_real_
    ))
  }

  rho_clamped <- pmin(pmax(dat$rho, -0.999999), 0.999999)
  yi <- atanh(rho_clamped)
  vi <- 1 / (dat$n_samples - 3)
  fit <- metafor::rma.uni(yi = yi, vi = vi, method = "REML")

  tibble(
    meta_z = as.numeric(fit$b),
    meta_se_z = as.numeric(fit$se),
    meta_rho = tanh(as.numeric(fit$b)),
    ci_low_rho = tanh(as.numeric(fit$ci.lb)),
    ci_high_rho = tanh(as.numeric(fit$ci.ub)),
    meta_p = as.numeric(fit$pval),
    tau2 = as.numeric(fit$tau2),
    I2 = as.numeric(fit$I2),
    Q_p = as.numeric(fit$QEp)
  )
}

meta_all <- cor_long |>
  group_by(Axis, Species_full, Species) |>
  summarise(
    k_cohorts = n_distinct(Cohort),
    total_n = sum(n_samples),
    cohorts = paste(sort(unique(Cohort)), collapse = ";"),
    KO_sets = paste(sort(unique(KO_set)), collapse = ";"),
    .groups = "drop"
  ) |>
  filter(k_cohorts >= 3L) |>
  left_join(
    cor_long |>
      group_by(Axis, Species_full, Species) |>
      filter(n_distinct(Cohort) >= 3L) |>
      group_modify(~meta_one(.x)) |>
      ungroup(),
    by = c("Axis", "Species_full", "Species")
  ) |>
  group_by(Axis) |>
  mutate(meta_q = p.adjust(meta_p, method = "BH")) |>
  ungroup() |>
  arrange(Axis, desc(meta_rho), meta_q, Species)

if (anyDuplicated(meta_all[c("Axis", "Species_full")])) {
  stop("Duplicate axis/species rows after meta-analysis.")
}

top_per_axis <- meta_all |>
  filter(meta_rho > 0, meta_q < 0.05) |>
  group_by(Axis) |>
  arrange(desc(meta_rho), meta_q, Species, .by_group = TRUE) |>
  slice_head(n = 15L) |>
  mutate(axis_rank = row_number()) |>
  ungroup()

selection_counts <- top_per_axis |>
  count(Axis, name = "n_selected")
if (!identical(sort(selection_counts$n_selected), c(15L, 15L))) {
  stop("Fewer than 15 positive BH-significant species were available for one or both axes.")
}

selected_species <- top_per_axis |>
  arrange(factor(Axis, levels = c("Production", "Degradation")), axis_rank) |>
  distinct(Species_full, .keep_all = TRUE) |>
  transmute(
    Species_full,
    Species,
    selected_by = Axis,
    selected_axis_rank = axis_rank,
    selected_axis_meta_rho = meta_rho,
    selected_axis_meta_q = meta_q
  )

plot_long <- selected_species |>
  select(Species_full, Species, selected_by, selected_axis_rank) |>
  crossing(Axis = c("Production", "Degradation")) |>
  left_join(
    meta_all |>
      select(Axis, Species_full, meta_rho, meta_q, k_cohorts, total_n, I2),
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
  arrange(factor(selected_by, levels = c("Production", "Degradation")), selected_axis_rank)

if (any(disease_overlap$disease_enrichment == "Both (unexpected)")) {
  stop("A selected species was recurrently enriched in both directions; inspect identifiers.")
}

row_order <- selected_species |>
  arrange(factor(selected_by, levels = c("Degradation", "Production")), desc(selected_axis_rank)) |>
  pull(Species)

label_map <- selected_species |>
  mutate(label = paste0("<i>", gsub("_", " ", Species), "</i>")) |>
  select(Species, label)

plot_long <- plot_long |>
  mutate(
    Axis = factor(Axis, levels = c("Production", "Degradation")),
    Species = factor(Species, levels = row_order)
  )

group_boundary <- 15.5

p <- ggplot(plot_long, aes(x = Axis, y = Species, fill = meta_rho)) +
  geom_tile(color = "white", linewidth = 0.75, width = 0.96, height = 0.96) +
  geom_hline(yintercept = group_boundary, linewidth = 0.65, color = "#414141") +
  scale_fill_gradient2(
    low = "#2B6CB0", mid = "#F7F7F7", high = "#C43C39",
    midpoint = 0, limits = c(-0.60, 0.60), oob = scales::squish,
    name = expression("Meta-analytic Spearman " * rho)
  ) +
  scale_y_discrete(labels = setNames(label_map$label, label_map$Species), expand = c(0, 0)) +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  labs(
    title = "Sarcosine functional associations in CRC WGS",
    subtitle = paste0(
      "Four-cohort random-effects integration of within-cohort Spearman correlations\n",
      "Top 15 positive BH-significant associations per functional axis"
    ),
    x = NULL,
    y = NULL,
    caption = paste0(
      "Production KO panel: K00315/K00552/K08688  •  ",
      "Degradation KO panel: K00301–K00306"
    )
  ) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 13, base_family = "Arial") +
  theme(
    plot.title = element_text(size = 18, face = "bold", color = "#111111", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11.5, color = "#4D4D4D", lineheight = 1.15, margin = margin(b = 15)),
    plot.caption = element_text(size = 9.5, color = "#555555", hjust = 0, margin = margin(t = 12)),
    axis.text.x = element_text(size = 13, face = "bold", color = "#111111", margin = margin(b = 8)),
    axis.text.y = ggtext::element_markdown(size = 10.5, color = "#1A1A1A", margin = margin(r = 8)),
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

png_file <- file.path(out_dir, "Fig_CRC_WGS_species_Production_Degradation_meta_heatmap.png")
pdf_file <- file.path(out_dir, "Fig_CRC_WGS_species_Production_Degradation_meta_heatmap.pdf")
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
    stop("sips failed to create standardized PNG: ", paste(sips_output, collapse = "\n"))
  }
  if (!file.copy(tmp, path, overwrite = TRUE)) stop("Could not install standardized PNG")
  unlink(tmp)

  after <- png::readPNG(path, native = FALSE, info = TRUE)
  if (!identical(dim(before), dim(after)) || !identical(as.numeric(before), as.numeric(after))) {
    stop("Embedding sRGB/dpi metadata changed decoded PNG pixels")
  }

  metadata <- system2(
    "/usr/bin/sips",
    args = c("-g", "profile", "-g", "dpiWidth", "-g", "dpiHeight", shQuote(path)),
    stdout = TRUE,
    stderr = TRUE
  )
  if (!any(grepl("profile: sRGB IEC61966-2.1", metadata, fixed = TRUE)) ||
      !any(grepl(paste0("dpiWidth: ", dpi, ".000"), metadata, fixed = TRUE)) ||
      !any(grepl(paste0("dpiHeight: ", dpi, ".000"), metadata, fixed = TRUE))) {
    stop("PNG profile/resolution validation failed: ", paste(metadata, collapse = "\n"))
  }
}

standardize_png(png_file, DPI)

if (!file.copy(png_file, transport_file, overwrite = TRUE)) {
  stop("Could not create PowerPoint transport copy: ", transport_file)
}
for (path in c(png_file, transport_file)) {
  invisible(system2("/usr/bin/xattr", args = c("-c", shQuote(path)), stdout = TRUE, stderr = TRUE))
  invisible(system2("/bin/chmod", args = c("644", shQuote(path)), stdout = TRUE, stderr = TRUE))
  invisible(system2("/usr/bin/chflags", args = c("nohidden", shQuote(path)), stdout = TRUE, stderr = TRUE))
  remaining_xattrs <- system2("/usr/bin/xattr", args = c("-l", shQuote(path)), stdout = TRUE, stderr = TRUE)
  if (length(remaining_xattrs) > 0L) {
    stop("Extended attributes remain on PNG: ", path, "\n", paste(remaining_xattrs, collapse = "\n"))
  }
}
if (!identical(unname(tools::md5sum(png_file)), unname(tools::md5sum(transport_file)))) {
  stop("PowerPoint transport copy is not byte-identical to the analytical PNG")
}

grDevices::cairo_pdf(pdf_file, width = 9.5, height = 12.0, family = "Arial")
print(p)
dev.off()

write.csv(meta_all, file.path(out_dir, "meta_correlations_all_species.csv"), row.names = FALSE, na = "")
write.csv(top_per_axis, file.path(out_dir, "selected_top15_per_axis.csv"), row.names = FALSE, na = "")
write.csv(plot_long, file.path(out_dir, "heatmap_plot_values_long.csv"), row.names = FALSE, na = "")
write.csv(disease_overlap, file.path(out_dir, "selected_species_disease_enrichment_overlap.csv"), row.names = FALSE, na = "")
write.csv(sample_n_check, file.path(out_dir, "cohort_sample_counts.csv"), row.names = FALSE)

input_md5 <- data.frame(
  file = required_files,
  md5 = unname(tools::md5sum(required_files)),
  stringsAsFactors = FALSE
)
write.csv(input_md5, file.path(out_dir, "input_file_md5.csv"), row.names = FALSE)

validation_lines <- c(
  "CRC WGS Production/Degradation meta-heatmap validation",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Input rows: Production=", nrow(prod), "; Degradation=", nrow(deg)),
  paste0("Cohorts: ", paste(expected_cohorts, collapse = ", ")),
  paste0("Cohort sample sizes: ", paste(sample_n_check$Cohort, sample_n_check$n_samples, sep = "=", collapse = "; ")),
  paste0("Meta-analysis eligibility: species observed in >=3 cohorts"),
  paste0("Meta-analysed rows: ", nrow(meta_all)),
  paste0("Selected rows per axis: ", paste(selection_counts$Axis, selection_counts$n_selected, sep = "=", collapse = "; ")),
  paste0("Unique selected species: ", nrow(selected_species)),
  paste0("Recurrent Cancer-enriched overlaps (>=3/4): ", sum(disease_overlap$Cancer_recurrent_n >= 3L)),
  paste0("Recurrent Healthy-enriched overlaps (>=3/4): ", sum(disease_overlap$Healthy_recurrent_n >= 3L)),
  paste0("PNG exists and non-empty: ", file.exists(png_file) && file.info(png_file)$size > 0),
  paste0("PowerPoint transport copy byte-identical: ", identical(unname(tools::md5sum(png_file)), unname(tools::md5sum(transport_file)))),
  "No significance symbols are drawn; q values are retained in CSV outputs."
)
writeLines(validation_lines, file.path(out_dir, "validation_report.txt"))
capture.output(sessionInfo(), file = file.path(out_dir, "sessionInfo.txt"))

readme <- c(
  "CRC WGS species–sarcosine functional-potential heatmap",
  "",
  "Purpose",
  "- Integrate species associations with community sarcosine Production and Degradation potentials across four CRC WGS cohorts.",
  "",
  "Frozen inputs",
  paste0("- ", prod_file),
  paste0("- ", deg_file),
  paste0("- ", cancer_file),
  paste0("- ", healthy_file),
  "",
  "Method",
  "- Within-cohort Spearman rho values were taken from the frozen outputs; correlations were not recomputed from raw abundance matrices.",
  "- For each species and functional axis, rho was Fisher-z transformed and integrated across >=3 cohorts with a REML random-effects model using variance 1/(n-3).",
  "- BH correction was applied within each functional axis across all meta-analysed species.",
  "- The plot displays the 15 strongest positive associations per axis among species with BH q < 0.05.",
  "- The heatmap contains no significance symbols, as requested; exact p/q values remain in the CSV files.",
  "",
  "Interpretation guardrail",
  "- Production and Degradation are community metagenomic functional potentials derived from cohort-detected members of the predefined KO panels (Production: K00315/K00552/K08688; Degradation: K00301-K00306), not measured sarcosine flux or metabolite concentration.",
  "- Disease-enrichment overlap uses the strict recurrent criterion in the frozen memberships: q < 0.05 and |log2FC| > 1 in >=3 of 4 cohorts.",
  "",
  "PowerPoint transport copy",
  paste0("- ", transport_file),
  "- The transport PNG is byte-identical to the analytical PNG and uses an embedded sRGB IEC61966-2.1 profile with exact 400 dpi metadata.",
  "",
  "Independent verification",
  paste0("- Rscript ", file.path(dirname(script_file), "114_verify_CRC_WGS_ProdDeg_meta_heatmap_26.08.31.R")),
  paste0("- Report: ", file.path(out_dir, "independent_verification_report.txt")),
  "",
  "Reproduction",
  paste0("Rscript ", script_file)
)
writeLines(readme, file.path(out_dir, "README.txt"))

message("Completed: ", out_dir)
message("Figure: ", png_file)
