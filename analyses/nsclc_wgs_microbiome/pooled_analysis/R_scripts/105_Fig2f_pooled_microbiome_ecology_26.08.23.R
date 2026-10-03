#!/usr/bin/env Rscript
# Main Figure 2 candidate: pooled NSCLC gut microbiome ecology.
#
# Biological question
#   Across four independent NSCLC ICI WGS cohorts, do responders (R) and
#   non-responders (NR) differ in within-sample Shannon diversity and in overall
#   Bray-Curtis community composition?
#
# Statistical design
#   Unit: one pretreatment faecal WGS sample/patient.
#   Design: two independent response groups across four cohorts.
#   Alpha: lm(Shannon ~ cohort + response); the displayed P value is the
#          prespecified response coefficient after cohort adjustment. Because
#          only Shannon is promoted, no multiple-testing correction is applied.
#   Beta: adonis2(Bray-Curtis ~ cohort + response, by = "terms", 999
#         permutations); response is tested after cohort. A pooled ordination is
#         descriptive; the model, not visual separation, supplies inference.
#   QC: recomputed values must match the previously saved alpha and PERMANOVA
#       tables; sample counts, ID alignment, row sums and PNG properties are
#       asserted. PERMDISP is reported as an interpretation diagnostic.
#
# Workflow order
#   1. Read metadata and species relative-abundance tables.
#   2. Filter to response-labelled NSCLC WGS samples and align IDs.
#   3. Recompute Shannon on each cohort's species table.
#   4. Build the union species matrix, zero-fill, re-TSS, then calculate
#      Bray-Curtis and PCoA.
#   5. Fit the prespecified cohort-adjusted models and verify prior results.
#   6. Render the display-only composite and export its source data/provenance.
#
# Reproducibility
#   Run from the project root:
#   Rscript \
#     '공공_Metabolomics&Metagenomics_분석모음/HGMT_NSCLC_ICI_RvsNR_WGS/pooled_analysis/R_scripts/105_Fig2f_pooled_microbiome_ecology_26.08.23.R'

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(here)
})

## ---- paths and constants ---------------------------------------------------
project_root <- here::here()
analysis_rel <- file.path(
  "공공_Metabolomics&Metagenomics_분석모음",
  "HGMT_NSCLC_ICI_RvsNR_WGS"
)
analysis_dir <- file.path(project_root, analysis_rel)
pooled_dir <- file.path(analysis_dir, "pooled_analysis")
prior_dir <- file.path(pooled_dir, "results")
out_dir <- file.path(prior_dir, "Fig2f_pooled_microbiome_ecology_26.08.23")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cohort_dirs <- c(
  PRJNA751792 = "NSCLC_PRJNA751792",
  PRJNA1023797 = "NSCLC_PRJNA1023797",
  PRJEB22863 = "NSCLC_RCC_PRJEB22863",
  PRJEB26531 = "NSCLC_PRJEB26531"
)
nsclc_label <- "Carcinoma, Non-Small-Cell Lung"
response_cols <- c(R = "#2E5F8A", NR = "#C47B3B")
expected_counts <- c(R = 445L, NR = 404L)
expected_union_species <- 2598L
expected_alpha_p <- 0.0372575021751629
expected_beta_r2 <- 0.00326905870004063
expected_beta_p <- 0.001

input_files <- character()
record_input <- function(path) {
  stopifnot(file.exists(path))
  input_files <<- unique(c(input_files, path))
  path
}

## ---- loaders: retained from the validated pooled-diversity workflow --------
parse_one_response <- function(description) {
  if (is.na(description)) return(NA_character_)
  matched <- regmatches(description, regexpr("response_group:[^;]*", description))
  if (length(matched) == 0L || !nzchar(matched)) {
    matched <- regmatches(description, regexpr("response:[^;]*", description))
  }
  if (length(matched) == 0L) return(NA_character_)
  value <- trimws(sub("^response(_group)?:[[:space:]]*", "", matched))
  if (value %in% c("R", "NR")) value else NA_character_
}

parse_response <- function(x) {
  vapply(as.character(x), parse_one_response, character(1), USE.NAMES = FALSE)
}

read_metadata <- function(cohort_dir) {
  files <- list.files(cohort_dir, "^selected_project_.*\\.txt$", full.names = TRUE)
  stopifnot(length(files) == 1L)
  path <- record_input(files)
  lines <- sub("\r$", "", readLines(path, warn = FALSE))
  header <- strsplit(lines[2], "\t", fixed = TRUE)[[1]]
  rows <- lines[-(1:2)]
  rows <- rows[nzchar(rows)]
  matrix_rows <- do.call(
    rbind,
    lapply(strsplit(rows, "\t", fixed = TRUE), function(x) x[seq_len(length(header))])
  )
  result <- as.data.table(matrix_rows)
  setnames(result, header)
  result
}

read_species <- function(cohort_dir) {
  files <- list.files(cohort_dir, "^Bacteria_.*\\.txt$", full.names = TRUE)
  stopifnot(length(files) == 1L)
  path <- record_input(files)
  long <- fread(path, sep = "\t", header = TRUE, quote = "", showProgress = FALSE)
  stopifnot(ncol(long) == 3L)
  setnames(long, c("Taxa", "RunID", "Abundance"))
  stopifnot(is.numeric(long$Abundance), all(long$Abundance >= 0, na.rm = TRUE))
  species <- long[
    grepl("|s__", Taxa, fixed = TRUE) & !grepl("|t__", Taxa, fixed = TRUE)
  ]
  wide <- dcast(species, RunID ~ Taxa, value.var = "Abundance", fill = 0,
                fun.aggregate = sum)
  matrix_values <- as.matrix(wide[, -1L, with = FALSE])
  rownames(matrix_values) <- wide$RunID
  matrix_values
}

## ---- load, inspect and align ------------------------------------------------
alpha_parts <- list()
species_parts <- list()
metadata_parts <- list()

for (cohort in names(cohort_dirs)) {
  cohort_dir <- file.path(analysis_dir, cohort_dirs[[cohort]])
  metadata <- read_metadata(cohort_dir)
  metadata <- metadata[
    `Assay type` == "WGS" & `Phenotype name` == nsclc_label
  ]
  metadata[, response := parse_response(`Sample description`)]
  metadata <- metadata[!is.na(response)]
  metadata[, response := factor(response, levels = c("NR", "R"))]

  species_matrix <- read_species(cohort_dir)
  run_ids <- metadata[["Run ID"]]
  stopifnot(!anyDuplicated(run_ids), all(run_ids %in% rownames(species_matrix)))
  species_matrix <- species_matrix[run_ids, , drop = FALSE]
  stopifnot(identical(rownames(species_matrix), run_ids))

  sample_ids <- paste0(cohort, "::", run_ids)
  alpha_parts[[cohort]] <- data.table(
    sample = sample_ids,
    cohort = cohort,
    response = as.character(metadata$response),
    Shannon = diversity(species_matrix, index = "shannon", MARGIN = 1)
  )
  rownames(species_matrix) <- sample_ids
  species_parts[[cohort]] <- species_matrix
  metadata_parts[[cohort]] <- data.table(
    sample = sample_ids,
    cohort = cohort,
    response = as.character(metadata$response)
  )
}

alpha_data <- rbindlist(alpha_parts)
ordination_metadata <- rbindlist(metadata_parts)
alpha_data[, cohort := factor(cohort, levels = names(cohort_dirs))]
alpha_data[, response := factor(response, levels = c("NR", "R"))]
ordination_metadata[, cohort := factor(cohort, levels = names(cohort_dirs))]
ordination_metadata[, response := factor(response, levels = c("NR", "R"))]

stopifnot(
  nrow(alpha_data) == 849L,
  all(as.integer(table(alpha_data$response)[names(expected_counts)]) == expected_counts),
  !anyNA(alpha_data$Shannon),
  all(is.finite(alpha_data$Shannon)),
  identical(alpha_data$sample, ordination_metadata$sample)
)

## ---- Shannon: cohort-adjusted response effect ------------------------------
alpha_model <- lm(Shannon ~ cohort + response, data = alpha_data)
alpha_coef <- summary(alpha_model)$coefficients
alpha_beta <- unname(alpha_coef["responseR", "Estimate"])
alpha_p <- unname(alpha_coef["responseR", "Pr(>|t|)"])

prior_alpha_path <- record_input(file.path(prior_dir, "pooled_alpha_combined_stats.csv"))
prior_alpha <- fread(prior_alpha_path)[index == "Shannon"]
stopifnot(
  nrow(prior_alpha) == 1L,
  isTRUE(all.equal(alpha_beta, prior_alpha$beta_adj_RvsNR, tolerance = 1e-12)),
  isTRUE(all.equal(alpha_p, prior_alpha$p_adj_cohort, tolerance = 1e-12)),
  isTRUE(all.equal(alpha_p, expected_alpha_p, tolerance = 1e-12))
)

## ---- Bray-Curtis, PCoA and cohort-adjusted PERMANOVA -----------------------
all_species <- sort(unique(unlist(lapply(species_parts, colnames))))
community <- matrix(
  0,
  nrow = nrow(ordination_metadata),
  ncol = length(all_species),
  dimnames = list(ordination_metadata$sample, all_species)
)
for (cohort in names(species_parts)) {
  block <- species_parts[[cohort]]
  community[rownames(block), colnames(block)] <- block
}
stopifnot(
  identical(rownames(community), ordination_metadata$sample),
  ncol(community) == expected_union_species,
  all(rowSums(community) > 0)
)
community <- community / rowSums(community)
stopifnot(max(abs(rowSums(community) - 1)) < 1e-9)

bray <- vegdist(community, method = "bray")
pcoa_fit <- cmdscale(bray, k = 2, eig = TRUE)
variance_explained <- pcoa_fit$eig / sum(pcoa_fit$eig[pcoa_fit$eig > 0]) * 100
pcoa_data <- data.table(
  sample = rownames(pcoa_fit$points),
  PCoA1 = pcoa_fit$points[, 1],
  PCoA2 = pcoa_fit$points[, 2]
)
pcoa_data <- merge(pcoa_data, ordination_metadata, by = "sample", sort = FALSE)
pcoa_data <- pcoa_data[match(ordination_metadata$sample, sample)]
stopifnot(identical(pcoa_data$sample, ordination_metadata$sample))

set.seed(42)
permanova <- adonis2(
  bray ~ cohort + response,
  data = ordination_metadata,
  by = "terms",
  permutations = 999
)
beta_r2 <- unname(permanova["response", "R2"])
beta_p <- unname(permanova["response", "Pr(>F)"])

prior_beta_path <- record_input(file.path(prior_dir, "pooled_beta_combined_permanova.csv"))
prior_beta <- fread(prior_beta_path)[term == "response|cohort"]
stopifnot(
  nrow(prior_beta) == 1L,
  isTRUE(all.equal(beta_r2, prior_beta$R2, tolerance = 1e-12)),
  isTRUE(all.equal(beta_p, prior_beta$p, tolerance = 1e-12)),
  isTRUE(all.equal(beta_r2, expected_beta_r2, tolerance = 1e-12)),
  isTRUE(all.equal(beta_p, expected_beta_p, tolerance = 1e-12))
)

# PERMDISP is a diagnostic only; it does not replace the prespecified PERMANOVA.
dispersion <- betadisper(bray, ordination_metadata$response)
set.seed(42)
dispersion_test <- permutest(dispersion, permutations = 999)
dispersion_f <- unname(dispersion_test$tab[1, "F"])
dispersion_p <- unname(dispersion_test$tab[1, "Pr(>F)"])

## ---- export exact panel source data and statistics -------------------------
alpha_export <- copy(alpha_data)
alpha_export[, cohort := as.character(cohort)]
alpha_export[, response := as.character(response)]
fwrite(alpha_export, file.path(out_dir, "Fig2f_source_data_Shannon.csv"))

pcoa_export <- copy(pcoa_data)
pcoa_export[, cohort := as.character(cohort)]
pcoa_export[, response := as.character(response)]
fwrite(pcoa_export, file.path(out_dir, "Fig2f_source_data_PCoA.csv"))

stats_export <- data.table(
  component = c("Shannon", "Bray-Curtis", "Bray-Curtis dispersion"),
  model_or_test = c(
    "lm(Shannon ~ cohort + response)",
    "adonis2(Bray ~ cohort + response), response after cohort, 999 permutations",
    "PERMDISP by response, 999 permutations"
  ),
  effect_name = c("beta R vs NR", "partial/sequential R2 response|cohort", "F"),
  effect = c(alpha_beta, beta_r2, dispersion_f),
  p = c(alpha_p, beta_p, dispersion_p),
  n = nrow(alpha_data)
)
fwrite(stats_export, file.path(out_dir, "Fig2f_statistics.csv"))

counts_export <- alpha_data[, .N, by = .(cohort, response)]
counts_export[, cohort := as.character(cohort)]
counts_export[, response := as.character(response)]
fwrite(counts_export, file.path(out_dir, "Fig2f_sample_counts.csv"))

input_checksums <- data.table(
  file = sub(paste0("^", project_root, "/?"), "", input_files),
  md5 = unname(tools::md5sum(input_files))
)
fwrite(input_checksums, file.path(out_dir, "Fig2f_input_md5.tsv"), sep = "\t")

## ---- publication styling ---------------------------------------------------
theme_main <- function() {
  theme_classic(base_family = "Arial", base_size = 10.5) +
    theme(
      plot.title = element_text(size = 11.5, face = "bold", hjust = 0),
      axis.title = element_text(size = 10.5),
      axis.text = element_text(size = 9.2, colour = "black"),
      legend.title = element_text(size = 9.2),
      legend.text = element_text(size = 9.2),
      legend.key.width = grid::unit(0.9, "lines"),
      legend.key.height = grid::unit(0.8, "lines"),
      plot.margin = margin(6, 8, 6, 8)
    )
}

# Display order is R then NR; model reference remains NR.
alpha_plot_data <- copy(alpha_data)
alpha_plot_data[, response_display := factor(as.character(response), levels = c("R", "NR"))]
alpha_range <- diff(range(alpha_plot_data$Shannon))
bracket_y <- max(alpha_plot_data$Shannon) + 0.055 * alpha_range
tick_y <- bracket_y - 0.025 * alpha_range
label_y <- bracket_y + 0.075 * alpha_range

alpha_plot <- ggplot(alpha_plot_data, aes(response_display, Shannon)) +
  geom_boxplot(
    aes(fill = response_display, colour = response_display),
    width = 0.56,
    alpha = 0.24,
    linewidth = 0.7,
    outlier.shape = NA
  ) +
  geom_jitter(
    aes(colour = response_display),
    width = 0.13,
    height = 0,
    size = 0.8,
    alpha = 0.42,
    stroke = 0
  ) +
  stat_summary(
    fun = median,
    geom = "point",
    shape = 21,
    size = 2.2,
    stroke = 0.65,
    fill = "white",
    colour = "black"
  ) +
  annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y,
           linewidth = 0.5, colour = "black") +
  annotate("segment", x = 1, xend = 1, y = tick_y, yend = bracket_y,
           linewidth = 0.5, colour = "black") +
  annotate("segment", x = 2, xend = 2, y = tick_y, yend = bracket_y,
           linewidth = 0.5, colour = "black") +
  annotate(
    "text", x = 1.5, y = label_y,
    label = sprintf("Cohort-adjusted P = %.3f", alpha_p),
    family = "Arial", size = 3.2
  ) +
  scale_colour_manual(values = response_cols, guide = "none") +
  scale_fill_manual(values = response_cols, guide = "none") +
  scale_x_discrete(labels = c(
    R = sprintf("R\n(n = %d)", expected_counts[["R"]]),
    NR = sprintf("NR\n(n = %d)", expected_counts[["NR"]])
  )) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.22))) +
  labs(title = "Shannon diversity", x = "ICI response", y = "Shannon index") +
  theme_main()

pcoa_plot_data <- copy(pcoa_data)
pcoa_plot_data[, response_display := factor(as.character(response), levels = c("R", "NR"))]

pcoa_plot <- ggplot(pcoa_plot_data, aes(PCoA1, PCoA2, colour = response_display)) +
  stat_ellipse(
    aes(group = response_display), type = "norm", level = 0.95,
    linewidth = 0.8, alpha = 0.95, show.legend = FALSE
  ) +
  geom_point(size = 0.8, alpha = 0.48, stroke = 0) +
  scale_colour_manual(
    values = response_cols, breaks = c("R", "NR"),
    labels = c("Responder (R)", "Non-responder (NR)"), name = NULL
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1))) +
  labs(
    title = "Bray–Curtis composition",
    subtitle = sprintf("Response | cohort: R² = %.4f, P = %.3f", beta_r2, beta_p),
    x = sprintf("PCoA 1 (%.1f%%)", variance_explained[1]),
    y = sprintf("PCoA 2 (%.1f%%)", variance_explained[2])
  ) +
  theme_main() +
  theme(
    plot.subtitle = element_text(size = 9.2, colour = "black", margin = margin(b = 3)),
    legend.position = "top", legend.justification = "left",
    legend.margin = margin(0, 0, 1, 0)
  )

composite <- alpha_plot + pcoa_plot +
  plot_layout(widths = c(0.82, 1.18)) +
  plot_annotation(
    title = "Pooled gut microbiome ecology",
    theme = theme(
      plot.title = element_text(
        family = "Arial", face = "bold", size = 13, hjust = 0,
        margin = margin(0, 0, 3, 4)
      ),
      plot.margin = margin(4, 5, 4, 5)
    )
  )

## ---- render: PowerPoint RGB PNG, high-resolution PNG and vector SVG --------
save_ragg_png <- function(filename, dpi) {
  ggsave(
    file.path(out_dir, filename), composite, device = ragg::agg_png,
    width = 7.8, height = 3.35, units = "in", dpi = dpi,
    background = "white", scaling = 1
  )
}

save_ragg_png("Fig2f_pooled_microbiome_ecology_PPT.png", 300)
save_ragg_png("Fig2f_pooled_microbiome_ecology_600dpi.png", 600)
# Cairo PDF is deliberately not used: this Mac R build lacks its XQuartz-linked
# Cairo libraries. Remove a stale partial PDF from an interrupted earlier run;
# SVG is the verified vector deliverable.
stale_pdf <- file.path(out_dir, "Fig2f_pooled_microbiome_ecology.pdf")
if (file.exists(stale_pdf)) unlink(stale_pdf)
ggsave(
  file.path(out_dir, "Fig2f_pooled_microbiome_ecology.svg"), composite,
  device = svglite::svglite, width = 7.8, height = 3.35, units = "in", bg = "white"
)

## ---- output checks ---------------------------------------------------------
png_paths <- file.path(out_dir, c(
  "Fig2f_pooled_microbiome_ecology_PPT.png",
  "Fig2f_pooled_microbiome_ecology_600dpi.png"
))
png_dims <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  stopifnot(identical(
    readBin(con, "raw", n = 8L),
    as.raw(c(137, 80, 78, 71, 13, 10, 26, 10))
  ))
  invisible(readBin(con, "integer", n = 1L, size = 4L, endian = "big"))
  stopifnot(identical(rawToChar(readBin(con, "raw", n = 4L)), "IHDR"))
  c(
    width = readBin(con, "integer", n = 1L, size = 4L, endian = "big"),
    height = readBin(con, "integer", n = 1L, size = 4L, endian = "big")
  )
}
dimensions <- t(vapply(png_paths, png_dims, integer(2)))
stopifnot(
  identical(unname(dimensions[1, ]), c(2340L, 1005L)),
  identical(unname(dimensions[2, ]), c(4680L, 2010L))
)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
log_lines <- c(
  sprintf("Generated: %s", format(Sys.time(), tz = "Asia/Seoul", usetz = TRUE)),
  "Purpose: display-only Main Figure 2 composite from the validated pooled diversity workflow.",
  sprintf("Samples: n=%d; R=%d; NR=%d; cohorts=%d", nrow(alpha_data),
          expected_counts[["R"]], expected_counts[["NR"]], length(cohort_dirs)),
  sprintf("Species matrix: %d samples x %d union species; re-TSS row sums verified.",
          nrow(community), ncol(community)),
  sprintf("Shannon: beta(R vs NR)=%.12f; cohort-adjusted P=%.12f; prior table match=TRUE.",
          alpha_beta, alpha_p),
  sprintf("Bray-Curtis response|cohort: R2=%.12f; P=%.3f; prior table match=TRUE.",
          beta_r2, beta_p),
  sprintf("PERMDISP diagnostic: F=%.6f; P=%.3f.", dispersion_f, dispersion_p),
  sprintf("PCoA variance displayed: axis1=%.3f%%; axis2=%.3f%%.",
          variance_explained[1], variance_explained[2]),
  "Palette: R=#2E5F8A; NR=#C47B3B.",
  sprintf("PPT PNG dimensions: %d x %d px; high-resolution PNG: %d x %d px.",
          dimensions[1, 1], dimensions[1, 2], dimensions[2, 1], dimensions[2, 2])
)
writeLines(log_lines, file.path(out_dir, "Fig2f_analysis_log.txt"))

# Write this last so every other final deliverable, including the final log, is
# covered by the checksum manifest from the same run.
output_files <- list.files(out_dir, full.names = TRUE)
output_files <- output_files[basename(output_files) != "Fig2f_output_md5.tsv"]
output_checksums <- data.table(
  file = basename(output_files), md5 = unname(tools::md5sum(output_files))
)
fwrite(output_checksums, file.path(out_dir, "Fig2f_output_md5.tsv"), sep = "\t")
cat(paste(log_lines, collapse = "\n"), "\n")
