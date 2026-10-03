#!/usr/bin/env Rscript

# Corrected cohort-adjusted PERMANOVA and publication render for current
# Supplementary Figure 3b (pooled CRC WGS Bray-Curtis PCoA).
#
# Statistical question:
#   Does disease group explain species-level Bray-Curtis composition after
#   removing cohort effects, with labels exchangeable only within cohort?
#
# Primary model:
#   adonis2(distance ~ Group + Cohort, by = "margin",
#           strata = Cohort, permutations = 999)
#
# This script does not overwrite an earlier result. Outputs are written to a
# new dated folder under results_integrated.

suppressPackageStartupMessages({
  library(vegan)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})
if (!requireNamespace("ragg", quietly = TRUE)) {
  stop("Package 'ragg' is required for 300-dpi PNG output.")
}

SEED <- 42L
N_PERM <- 999L

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine this script's path.")
SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_arg)))
BASE_DIR <- normalizePath(file.path(SCRIPT_DIR, ".."))
SOURCE_CSV <- normalizePath(file.path(
  BASE_DIR, "..", "..", "..", "Figure_Panel_Source_Data_26.07.23",
  "SupFig2b.csv"
))
OUT_DIR <- file.path(
  BASE_DIR, "results_integrated",
  "For_Publication_PERMANOVA_CORRECTED_26.08.20"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

COHORT_DIRS <- c(
  PRJEB6070 = file.path(BASE_DIR, "PRJEB6070_CRC_AdenomatousPolyps"),
  PRJEB10878 = file.path(BASE_DIR, "PRJEB10878_CRC"),
  PRJEB27928 = file.path(BASE_DIR, "PRJEB27928_CRC"),
  PRJNA429097 = file.path(BASE_DIR, "PRJNA429097_CRC")
)
stopifnot(all(dir.exists(COHORT_DIRS)), file.exists(SOURCE_CSV))

read_hgmt_metadata <- function(cohort_id, cohort_dir) {
  metadata_files <- sort(list.files(
    cohort_dir, pattern = "^selected_project_.*\\.txt$", full.names = TRUE
  ), decreasing = TRUE)
  if (!length(metadata_files)) stop("No metadata file for ", cohort_id)
  metadata_file <- metadata_files[[1L]]

  metadata_lines <- readLines(metadata_file)
  if (length(metadata_lines) < 3L) stop("Malformed metadata: ", metadata_file)
  header <- strsplit(metadata_lines[[2L]], "\t", fixed = TRUE)[[1L]]
  data_lines <- metadata_lines[seq.int(3L, length(metadata_lines))]
  data_lines <- data_lines[nzchar(data_lines)]
  split_rows <- strsplit(data_lines, "\t", fixed = TRUE)
  split_rows <- lapply(split_rows, function(x) {
    if (length(x) >= length(header)) x[seq_along(header)] else
      c(x, rep(NA_character_, length(header) - length(x)))
  })
  metadata <- as.data.frame(do.call(rbind, split_rows), stringsAsFactors = FALSE)
  names(metadata) <- gsub(" ", ".", header, fixed = TRUE)

  required <- c("Run.ID", "Phenotype.name")
  missing <- setdiff(required, names(metadata))
  if (length(missing)) {
    stop("Metadata missing columns: ", paste(missing, collapse = ", "))
  }
  phenotypes <- unique(metadata$Phenotype.name)
  cancer_phenotypes <- setdiff(phenotypes, c("Health", "Adenomatous Polyps"))
  if (length(cancer_phenotypes) != 1L) {
    stop(cohort_id, ": expected one cancer phenotype, found ",
         paste(cancer_phenotypes, collapse = ", "))
  }

  metadata <- metadata %>%
    filter(Phenotype.name %in% c("Health", cancer_phenotypes[[1L]])) %>%
    mutate(
      Group = if_else(Phenotype.name == "Health", "Healthy", "Cancer"),
      Cohort = cohort_id
    )
  if ("Assay.type" %in% names(metadata)) {
    metadata <- metadata %>% filter(Assay.type == "WGS")
  }
  attr(metadata, "input_file") <- metadata_file
  metadata
}

message("Loading cohort metadata...")
metadata_list <- Map(read_hgmt_metadata, names(COHORT_DIRS), COHORT_DIRS)
meta_pooled <- bind_rows(metadata_list)
stopifnot(!anyDuplicated(meta_pooled$Run.ID))

message("Loading species-level HGMT profiles...")
bacteria_list <- Map(function(cohort_id, cohort_dir) {
  bacteria_files <- sort(list.files(
    cohort_dir, pattern = "^Bacteria_.*\\.txt$", full.names = TRUE
  ), decreasing = TRUE)
  if (length(bacteria_files) != 1L) {
    stop(cohort_id, ": expected one Bacteria file, found ",
         length(bacteria_files))
  }
  bacteria <- read.delim(bacteria_files[[1L]], stringsAsFactors = FALSE)
  names(bacteria) <- trimws(names(bacteria))
  required <- c("Taxa", "Run.ID", "Abundance")
  missing <- setdiff(required, names(bacteria))
  if (length(missing)) {
    stop("Bacteria file missing columns: ", paste(missing, collapse = ", "))
  }
  cohort_runs <- meta_pooled$Run.ID[meta_pooled$Cohort == cohort_id]
  out <- bacteria %>%
    filter(grepl("\\|s__", Taxa), !grepl("\\|t__", Taxa), Run.ID %in% cohort_runs) %>%
    select(Taxa, Run.ID, Abundance)
  attr(out, "input_file") <- bacteria_files[[1L]]
  out
}, names(COHORT_DIRS), COHORT_DIRS)

bacteria_long <- bind_rows(bacteria_list)
bacteria_wide <- bacteria_long %>%
  pivot_wider(names_from = Taxa, values_from = Abundance, values_fill = 0)
bacteria_matrix <- as.data.frame(bacteria_wide)
rownames(bacteria_matrix) <- bacteria_matrix$Run.ID
bacteria_matrix$Run.ID <- NULL
bacteria_matrix <- as.matrix(bacteria_matrix)
storage.mode(bacteria_matrix) <- "double"

sample_order <- intersect(rownames(bacteria_matrix), meta_pooled$Run.ID)
bacteria_matrix <- bacteria_matrix[sample_order, , drop = FALSE]
meta_matched <- meta_pooled[match(sample_order, meta_pooled$Run.ID), , drop = FALSE]
meta_matched$Group <- factor(meta_matched$Group, levels = c("Healthy", "Cancer"))
meta_matched$Cohort <- factor(meta_matched$Cohort, levels = names(COHORT_DIRS))

pcoa_source <- read.csv(SOURCE_CSV, stringsAsFactors = FALSE, check.names = FALSE)
required_pcoa <- c("Run.ID", "PC1", "PC2", "Group", "Cohort")
missing_pcoa <- setdiff(required_pcoa, names(pcoa_source))
if (length(missing_pcoa)) {
  stop("PCoA source missing columns: ", paste(missing_pcoa, collapse = ", "))
}

expected_counts <- matrix(
  c(476L, 592L, 54L, 74L, 120L, 140L, 95L, 98L),
  nrow = 4L, byrow = TRUE,
  dimnames = list(names(COHORT_DIRS), c("Healthy", "Cancer"))
)
observed_counts <- with(meta_matched, table(Cohort, Group))
print(observed_counts)
stopifnot(
  nrow(bacteria_matrix) == 1649L,
  nrow(meta_matched) == 1649L,
  ncol(bacteria_matrix) > 0L,
  identical(dim(observed_counts), dim(expected_counts)),
  identical(rownames(observed_counts), rownames(expected_counts)),
  identical(colnames(observed_counts), colnames(expected_counts)),
  all(as.matrix(observed_counts) == expected_counts),
  setequal(meta_matched$Run.ID, pcoa_source$Run.ID),
  !anyNA(bacteria_matrix),
  all(is.finite(bacteria_matrix)),
  all(bacteria_matrix >= 0)
)
message("Verified 1,649 aligned samples (745 Healthy; 904 Cancer).")

message("Computing Bray-Curtis distances...")
bray <- vegan::vegdist(bacteria_matrix, method = "bray")

set.seed(SEED)
perm_naive <- vegan::adonis2(
  bray ~ Group, data = meta_matched, permutations = N_PERM
)
set.seed(SEED)
perm_full_omnibus <- vegan::adonis2(
  bray ~ Group + Cohort, data = meta_matched, permutations = N_PERM
)
set.seed(SEED)
perm_adjusted <- vegan::adonis2(
  bray ~ Group + Cohort,
  data = meta_matched,
  permutations = N_PERM,
  by = "margin",
  strata = meta_matched$Cohort
)
set.seed(SEED)
perm_sequential <- vegan::adonis2(
  bray ~ Cohort + Group,
  data = meta_matched,
  permutations = N_PERM,
  by = "terms",
  strata = meta_matched$Cohort
)

extract_row <- function(result, row_name, test_name, design) {
  if (!row_name %in% rownames(result)) {
    stop("Missing row '", row_name, "' in ", test_name)
  }
  data.frame(
    Test = test_name,
    Row = row_name,
    Df = unname(result[row_name, "Df"]),
    SumOfSqs = unname(result[row_name, "SumOfSqs"]),
    R2 = unname(result[row_name, "R2"]),
    F = unname(result[row_name, "F"]),
    P = unname(result[row_name, "Pr(>F)"]),
    Permutations = N_PERM,
    Seed = SEED,
    Design = design,
    stringsAsFactors = FALSE
  )
}

results <- bind_rows(
  extract_row(
    perm_naive, "Model", "Group_unadjusted",
    "bray ~ Group; free permutations"
  ),
  extract_row(
    perm_full_omnibus, "Model", "Group_plus_Cohort_full_model",
    "bray ~ Group + Cohort; omnibus model (not Group-adjusted R2)"
  ),
  extract_row(
    perm_adjusted, "Group", "Group_cohort_adjusted_primary",
    "bray ~ Group + Cohort; marginal Group effect; permutations within Cohort"
  ),
  extract_row(
    perm_sequential, "Group", "Group_after_Cohort_sensitivity",
    "bray ~ Cohort + Group; sequential Group effect; permutations within Cohort"
  )
)

primary <- results %>% filter(Test == "Group_cohort_adjusted_primary")
sensitivity <- results %>% filter(Test == "Group_after_Cohort_sensitivity")
stopifnot(
  nrow(primary) == 1L,
  nrow(sensitivity) == 1L,
  isTRUE(all.equal(primary$R2, sensitivity$R2, tolerance = 1e-12)),
  isTRUE(all.equal(primary$F, sensitivity$F, tolerance = 1e-12))
)

write.csv(
  results,
  file.path(OUT_DIR, "SupFig3b_PERMANOVA_results.csv"),
  row.names = FALSE
)
write.csv(
  as.data.frame.matrix(observed_counts),
  file.path(OUT_DIR, "SupFig3b_sample_counts.csv")
)

input_files <- c(
  vapply(metadata_list, attr, character(1L), which = "input_file"),
  vapply(bacteria_list, attr, character(1L), which = "input_file"),
  PCoA_source = SOURCE_CSV
)
checksum_table <- data.frame(
  Input = names(input_files),
  Path = unname(input_files),
  MD5 = unname(tools::md5sum(input_files)),
  stringsAsFactors = FALSE
)
write.csv(
  checksum_table,
  file.path(OUT_DIR, "SupFig3b_input_checksums.csv"),
  row.names = FALSE
)

format_p <- function(x) {
  if (x < 0.001) "p < 0.001" else paste0("p = ", sprintf("%.3f", x))
}
annotation <- paste0(
  "Cohort-adjusted PERMANOVA: R² = ",
  sprintf("%.4f", primary$R2), ", ", format_p(primary$P)
)

pcoa_source <- pcoa_source %>%
  mutate(
    Group = factor(Group, levels = c("Healthy", "Cancer")),
    Cohort = factor(Cohort, levels = names(COHORT_DIRS))
  )

cohort_shapes <- c(
  PRJEB6070 = 16, PRJEB10878 = 15, PRJEB27928 = 18, PRJNA429097 = 17
)
group_colors <- c(Healthy = "#1B9E8F", Cancer = "#C43C3C")

plot_obj <- ggplot(pcoa_source, aes(PC1, PC2, colour = Group, shape = Cohort)) +
  geom_point(size = 2, alpha = 0.65) +
  stat_ellipse(
    data = pcoa_source,
    aes(x = PC1, y = PC2, colour = Group, group = Group),
    inherit.aes = FALSE, level = 0.95, linetype = 2,
    linewidth = 0.85, show.legend = FALSE
  ) +
  scale_colour_manual(values = group_colors) +
  scale_shape_manual(values = cohort_shapes) +
  labs(
    title = "Bray–Curtis PCoA", subtitle = annotation,
    x = "PCoA1", y = "PCoA2"
  ) +
  guides(
    shape = guide_legend(nrow = 1, order = 1),
    colour = guide_legend(nrow = 1, order = 2)
  ) +
  theme_classic(base_family = "Arial", base_size = 15) +
  theme(
    plot.title = element_text(
      family = "Arial", size = 20, face = "bold", hjust = 0,
      margin = margin(b = 3)
    ),
    plot.subtitle = element_text(
      family = "Arial", size = 15, hjust = 0, margin = margin(b = 8)
    ),
    axis.title = element_text(family = "Arial", size = 17),
    axis.text = element_text(family = "Arial", size = 15, colour = "black"),
    legend.title = element_blank(),
    legend.text = element_text(family = "Arial", size = 15),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.box.just = "left",
    legend.spacing.y = grid::unit(0.04, "in"),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    plot.margin = margin(10, 12, 10, 10)
  )

png_path <- file.path(OUT_DIR, "SupFig3b_pooled_pcoa_corrected.png")
ragg::agg_png(
  filename = png_path, width = 9.5, height = 6.4, units = "in",
  res = 300, background = "white"
)
print(plot_obj)
grDevices::dev.off()
if (!file.exists(png_path) || file.info(png_path)$size <= 0) {
  stop("Corrected PNG was not written: ", png_path)
}

capture.output(
  sessionInfo(),
  file = file.path(OUT_DIR, "sessionInfo_SupFig3b_PERMANOVA.txt")
)

message("\nCorrected PERMANOVA results:")
print(results)
message("\nSaved: ", png_path)
