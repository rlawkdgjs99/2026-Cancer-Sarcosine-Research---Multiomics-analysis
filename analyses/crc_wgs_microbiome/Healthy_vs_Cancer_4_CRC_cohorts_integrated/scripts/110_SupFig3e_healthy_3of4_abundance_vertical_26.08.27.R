#!/usr/bin/env Rscript

# Supplementary Figure 3e: pooled Healthy-versus-Cancer abundance distributions
# for the eight species that met the strict Healthy-enriched criterion in
# exactly three of four CRC WGS cohorts.
#
# The species set is frozen by script 108. This script reads the same taxonomic
# sample roster used by the former Cancer-enriched Sup. Fig. 3e (n = 1,649),
# reconstructs zero-complete per-sample abundances for the eight selected taxa,
# and renders vertical Healthy/Cancer boxplots with raw points. Pooled Wilcoxon
# tests are descriptive; the selection/formal evidence remains the four frozen
# per-cohort differential-abundance analyses.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggtext)
  library(readr)
})

required_packages <- c("digest", "magick", "png")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "))
}

cmd_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", cmd_args, value = TRUE)
if (length(file_arg) != 1) stop("Run with Rscript so the script path can be resolved.")
SCRIPT_PATH <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
SCRIPT_DIR <- dirname(SCRIPT_PATH)
INTEGRATED_DIR <- dirname(SCRIPT_DIR)
CRC_ROOT <- dirname(INTEGRATED_DIR)
ANALYSIS_ROOT <- dirname(CRC_ROOT)

THEME_FILE <- file.path(ANALYSIS_ROOT, "_shared", "theme_nc_26.08.18.R")
if (!file.exists(THEME_FILE)) stop("Shared theme not found: ", THEME_FILE)
source(THEME_FILE)

OUTDIR <- file.path(
  INTEGRATED_DIR, "results_integrated",
  "SupFig3e_HEALTHY_3OF4_ABUNDANCE_26.08.27"
)
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

TARGET_FILE <- file.path(
  INTEGRATED_DIR, "results_integrated",
  "SupFig3de_HEALTHY_ENRICHED_26.08.27",
  "healthy_enriched_exactly_3of4_species.csv"
)

BACTERIA_FILES <- c(
  PRJEB10878 = file.path(CRC_ROOT, "PRJEB10878_CRC", "Bacteria_1771697876393.txt"),
  PRJEB27928 = file.path(CRC_ROOT, "PRJEB27928_CRC", "Bacteria_1771698146062.txt"),
  PRJEB6070 = file.path(
    CRC_ROOT, "PRJEB6070_CRC_AdenomatousPolyps", "Bacteria_1771697637952.txt"
  ),
  PRJNA429097 = file.path(CRC_ROOT, "PRJNA429097_CRC", "Bacteria_1771698356470.txt")
)

METADATA_FILES <- c(
  PRJEB10878 = file.path(
    CRC_ROOT, "PRJEB10878_CRC", "selected_project_1771697867394.txt"
  ),
  PRJEB27928 = file.path(
    CRC_ROOT, "PRJEB27928_CRC", "selected_project_1771698138019.txt"
  ),
  PRJEB6070 = file.path(
    CRC_ROOT, "PRJEB6070_CRC_AdenomatousPolyps",
    "selected_project_1771697628182.txt"
  ),
  PRJNA429097 = file.path(
    CRC_ROOT, "PRJNA429097_CRC", "selected_project_1771698343076.txt"
  )
)

EXPECTED_SHA256 <- c(
  target = "b07ca1bea99a5883b38c9ed1182cfffd7544087f031e549b3be0695c8bc8d8e3",
  bacteria_PRJEB10878 = "a48f57980e32af46fe8fe99b54c6114ec04958f23894f0831e85071755a16761",
  metadata_PRJEB10878 = "96899d0827e46464ab5d22ef1f9fe5190be8210b8e83bf9abd68c6690da1b8da",
  bacteria_PRJEB27928 = "92a21e2d2236561674853dc52b99a22539655dde62f99ce359ba915ef66a7b6f",
  metadata_PRJEB27928 = "19c43c5d25069b7dd33013eb649b63519cd5385af511548e5431d86e3ec39596",
  bacteria_PRJEB6070 = "af6cd42a0ecb3195d508a75e71f2089109d66143931cd738785cd7df0410eddc",
  metadata_PRJEB6070 = "4d536ed64feb75a945e29deebbea518b450c048d029f01beea78caaf2eb2f8a9",
  bacteria_PRJNA429097 = "2bde43412105bc8ace1c459d928f28d737528a14906489ad376e30a9b24a8f1e",
  metadata_PRJNA429097 = "970b7761ead98f55ec8b52d2e17b73a6412c7938027ad6214b9bb2f4c97f77b4"
)

DPI <- 600
WIDTH_IN <- 3.372
HEIGHT_IN <- 2.090
PSEUDO_DIV <- 5
SRGB_PROFILE <- "/System/Library/ColorSync/Profiles/sRGB Profile.icc"
if (!file.exists(SRGB_PROFILE)) stop("System sRGB profile not found")

INPUT_FILES <- c(
  target = TARGET_FILE,
  setNames(BACTERIA_FILES, paste0("bacteria_", names(BACTERIA_FILES))),
  setNames(METADATA_FILES, paste0("metadata_", names(METADATA_FILES)))
)
if (!all(file.exists(INPUT_FILES))) {
  stop("Missing input(s): ", paste(INPUT_FILES[!file.exists(INPUT_FILES)], collapse = "; "))
}
input_sha <- vapply(
  INPUT_FILES, digest::digest, character(1),
  algo = "sha256", file = TRUE, serialize = FALSE
)
if (!identical(unname(input_sha), unname(EXPECTED_SHA256[names(input_sha)]))) {
  bad <- names(input_sha)[input_sha != EXPECTED_SHA256[names(input_sha)]]
  stop("Input SHA-256 mismatch: ", paste(bad, collapse = ", "))
}

targets <- readr::read_csv(TARGET_FILE, show_col_types = FALSE, progress = FALSE)
required_target_columns <- c("Species_full", "Species", "n_cohorts_strict", "mean_log2FC")
if (!all(required_target_columns %in% names(targets))) stop("Target table columns changed")
if (nrow(targets) != 8L || any(targets$n_cohorts_strict != 3L) ||
    anyDuplicated(targets$Species_full)) {
  stop("Expected eight unique species, each strict in exactly three cohorts")
}
target_taxa <- targets$Species_full

read_groups <- function(path, cohort) {
  lines <- readLines(path, warn = FALSE)
  if (length(lines) < 3L) stop(cohort, ": malformed metadata file")
  header <- strsplit(lines[2], "\t", fixed = FALSE)[[1]]
  body <- lines[seq.int(3L, length(lines))]
  body <- body[nzchar(body)]
  rows <- lapply(strsplit(body, "\t", fixed = FALSE), function(x) {
    if (length(x) >= length(header)) x[seq_along(header)] else
      c(x, rep(NA_character_, length(header) - length(x)))
  })
  meta <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
  names(meta) <- gsub(" ", ".", header, fixed = TRUE)
  required <- c("Run.ID", "Phenotype.name")
  if (!all(required %in% names(meta))) stop(cohort, ": metadata columns changed")
  cancer_labels <- setdiff(
    unique(meta$Phenotype.name[!is.na(meta$Phenotype.name)]),
    c("Health", "Adenomatous Polyps")
  )
  if (length(cancer_labels) != 1L) {
    stop(cohort, ": expected one cancer phenotype; found ",
         paste(cancer_labels, collapse = ", "))
  }
  meta <- meta[meta$Phenotype.name %in% c("Health", cancer_labels), , drop = FALSE]
  if ("Assay.type" %in% names(meta)) {
    meta <- meta[meta$Assay.type == "WGS", , drop = FALSE]
  }
  out <- tibble(
    Cohort = cohort,
    RunID = meta$Run.ID,
    Group = ifelse(meta$Phenotype.name == "Health", "Healthy", "Cancer")
  ) %>%
    filter(!is.na(RunID), RunID != "")
  if (anyDuplicated(out$RunID)) stop(cohort, ": duplicated Run.ID in selected metadata")
  out
}

metadata <- bind_rows(lapply(names(METADATA_FILES), function(cohort) {
  read_groups(METADATA_FILES[[cohort]], cohort)
}))
if (anyDuplicated(metadata$RunID)) stop("Run IDs overlap across cohorts")

observed_rows <- lapply(names(BACTERIA_FILES), function(cohort) {
  b <- read.delim(
    BACTERIA_FILES[[cohort]], stringsAsFactors = FALSE,
    check.names = FALSE, quote = "", comment.char = ""
  )
  names(b) <- trimws(names(b))
  required <- c("Taxa", "Run ID", "Abundance")
  if (!all(required %in% names(b))) stop(cohort, ": bacteria columns changed")
  b <- b[
    grepl("\\|s__", b$Taxa) & !grepl("\\|t__", b$Taxa) &
      b$Taxa %in% target_taxa,
    required, drop = FALSE
  ]
  out <- tibble(
    Cohort = cohort,
    RunID = b[["Run ID"]],
    Species_full = b$Taxa,
    Abundance = as.numeric(b$Abundance)
  )
  if (any(!is.finite(out$Abundance)) || any(out$Abundance < 0)) {
    stop(cohort, ": invalid abundance values")
  }
  if (anyDuplicated(out[c("RunID", "Species_full")])) {
    stop(cohort, ": duplicated sample-species rows")
  }
  out
})
names(observed_rows) <- names(BACTERIA_FILES)
observed <- bind_rows(observed_rows)

present_samples <- unique(observed$RunID)
# Use every selected taxonomic sample, not merely samples with a non-zero target.
all_bacteria_samples <- unique(unlist(lapply(names(BACTERIA_FILES), function(cohort) {
  b <- read.delim(
    BACTERIA_FILES[[cohort]], stringsAsFactors = FALSE,
    check.names = FALSE, quote = "", comment.char = "",
    colClasses = c("character", "character", "numeric")
  )
  b[["Run ID"]]
})))
sample_meta <- metadata %>% filter(RunID %in% all_bacteria_samples)
if (nrow(sample_meta) != 1649L ||
    sum(sample_meta$Group == "Healthy") != 745L ||
    sum(sample_meta$Group == "Cancer") != 904L) {
  stop("Taxonomic roster changed; expected n=1,649 (Healthy 745, Cancer 904)")
}

abundance_long <- tidyr::crossing(
  RunID = sample_meta$RunID,
  Species_full = target_taxa
) %>%
  left_join(observed %>% select(RunID, Species_full, Abundance),
            by = c("RunID", "Species_full")) %>%
  mutate(Abundance = replace_na(Abundance, 0)) %>%
  left_join(sample_meta, by = "RunID")

if (nrow(abundance_long) != 1649L * 8L || any(is.na(abundance_long$Group)) ||
    anyDuplicated(abundance_long[c("RunID", "Species_full")])) {
  stop("Zero-complete sample-by-species table failed structural checks")
}
if (!all(vapply(split(abundance_long$Abundance, abundance_long$Species_full),
                  max, numeric(1)) > 0)) {
  stop("At least one target species is never detected")
}

stats <- bind_rows(lapply(seq_len(nrow(targets)), function(i) {
  sp <- targets$Species_full[i]
  x <- abundance_long %>% filter(Species_full == sp)
  h <- x$Abundance[x$Group == "Healthy"]
  c <- x$Abundance[x$Group == "Cancer"]
  tibble(
    Species_full = sp,
    Species = targets$Species[i],
    n_Healthy = length(h), n_Cancer = length(c),
    prevalence_Healthy_pct = mean(h > 0) * 100,
    prevalence_Cancer_pct = mean(c > 0) * 100,
    median_Healthy = median(h), median_Cancer = median(c),
    mean_Healthy = mean(h), mean_Cancer = mean(c),
    pooled_log2FC = log2((mean(c) + 1e-6) / (mean(h) + 1e-6)),
    pooled_Wilcoxon_p = suppressWarnings(
      wilcox.test(h, c, exact = FALSE, correct = TRUE)$p.value
    )
  )
})) %>%
  mutate(pooled_Wilcoxon_BH_q = p.adjust(pooled_Wilcoxon_p, method = "BH"))

q_label <- function(q) {
  if (is.na(q)) return(NA_character_)
  if (q < 0.001) return("q<0.001")
  if (q < 0.01) return(sprintf("q=%.3f", q))
  sprintf("q=%.2f", q)
}

species_labels <- c(
  Blautia_stercoris = "<i>B. stercoris</i>",
  Clostridium_sp_AF34_13 = "<i>Clostridium</i><br>sp. AF34-13",
  Anaerostipes_hadrus = "<i>A. hadrus</i>",
  Eubacterium_ventriosum = "<i>E. ventriosum</i>",
  Faecalibacillus_intestinalis = "<i>F. intestinalis</i>",
  Adlercreutzia_equolifaciens = "<i>A. equolifaciens</i>",
  Anaerobutyricum_hallii = "<i>A. hallii</i>",
  Blautia_wexlerae = "<i>B. wexlerae</i>"
)
if (!all(targets$Species %in% names(species_labels))) stop("Missing species display label")

species_order <- targets %>% arrange(mean_log2FC, Species) %>% pull(Species)
abundance_long <- abundance_long %>%
  left_join(targets %>% select(Species_full, Species), by = "Species_full") %>%
  mutate(
    Species = factor(Species, levels = species_order),
    Group = factor(Group, levels = c("Healthy", "Cancer"))
  ) %>%
  group_by(Species) %>%
  mutate(
    floor_value = min(Abundance[Abundance > 0], na.rm = TRUE) / PSEUDO_DIV,
    Plot_abundance = ifelse(Abundance <= 0, floor_value, Abundance)
  ) %>%
  ungroup()

annotation <- stats %>%
  mutate(
    Species = factor(Species, levels = species_order),
    label = vapply(pooled_Wilcoxon_BH_q, q_label, character(1))
  ) %>%
  left_join(
    abundance_long %>% group_by(Species) %>%
      summarise(y = max(Plot_abundance) * 1.45, .groups = "drop"),
    by = "Species"
  )

set.seed(42)
p <- ggplot(abundance_long, aes(x = Group, y = Plot_abundance, fill = Group)) +
  geom_boxplot(
    outlier.shape = NA, width = 0.62, linewidth = 0.28,
    colour = "#202020", alpha = 1
  ) +
  geom_jitter(
    position = position_jitter(width = 0.13, height = 0, seed = 42),
    shape = 21, size = 0.22, stroke = 0.10,
    colour = "#202020", fill = "#202020", alpha = 0.16
  ) +
  geom_boxplot(
    outlier.shape = NA, width = 0.62, linewidth = 0.36,
    colour = "#202020", fill = NA
  ) +
  geom_text(
    data = annotation,
    aes(x = 1.5, y = y, label = label),
    inherit.aes = FALSE, family = "Arial", size = 1.55,
    colour = "#202020"
  ) +
  facet_wrap(~Species, scales = "free_y", ncol = 4,
             labeller = as_labeller(species_labels)) +
  scale_y_log10(
    breaks = scales::breaks_log(n = 3),
    labels = scales::label_log(),
    expand = expansion(mult = c(0.02, 0.17))
  ) +
  scale_fill_manual(values = c(Healthy = COL_HEALTHY, Cancer = COL_CANCER)) +
  labs(
    title = "Healthy-enriched species",
    x = NULL, y = "Relative abundance (%)"
  ) +
  theme_minimal(base_family = "Arial", base_size = 6) +
  theme(
    plot.title = element_text(
      family = "Arial", face = "bold", size = 7.3,
      hjust = 0, margin = margin(b = 1.0)
    ),
    strip.text = ggtext::element_markdown(
      family = "Arial", size = 4.5, colour = "#111111",
      lineheight = 0.88, margin = margin(b = 0.7)
    ),
    strip.background = element_blank(),
    axis.text.x = element_text(
      family = "Arial", size = 4.3, colour = "#222222",
      angle = 24, hjust = 1, vjust = 1
    ),
    axis.text.y = element_text(family = "Arial", size = 4.1, colour = "#222222"),
    axis.title.y = element_text(family = "Arial", size = 5.1, margin = margin(r = 1.3)),
    axis.ticks = element_line(colour = "#202020", linewidth = 0.28),
    axis.line = element_line(colour = "#202020", linewidth = 0.35),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(colour = "#E3E3E3", linewidth = 0.18),
    panel.spacing.x = unit(3.0, "pt"),
    panel.spacing.y = unit(3.0, "pt"),
    legend.position = "none",
    plot.margin = margin(1.4, 2.0, 1.2, 1.5)
  )

standardize_png <- function(path) {
  img <- magick::image_read(path)
  img <- magick::image_background(img, "white", flatten = TRUE)
  img <- magick::image_convert(img, format = "png", colorspace = "sRGB", depth = 8)
  img <- magick::image_strip(img)
  magick::image_write(img, path = path, format = "png", depth = 8)
  before <- png::readPNG(path, native = FALSE, info = TRUE)

  tmp <- tempfile(fileext = ".png")
  output <- system2(
    "/usr/bin/sips",
    args = c(
      "-s", "dpiWidth", as.character(DPI),
      "-s", "dpiHeight", as.character(DPI),
      "-e", shQuote(SRGB_PROFILE), shQuote(path), "--out", shQuote(tmp)
    ),
    stdout = TRUE, stderr = TRUE
  )
  if (!file.exists(tmp) || file.info(tmp)$size <= 5000) {
    stop("sips failed: ", paste(output, collapse = "\n"))
  }
  if (!file.copy(tmp, path, overwrite = TRUE)) stop("Could not install standardized PNG")
  unlink(tmp)
  after <- png::readPNG(path, native = FALSE, info = TRUE)
  if (!identical(dim(before), dim(after)) ||
      !identical(as.numeric(before), as.numeric(after))) {
    stop("PNG metadata insertion changed decoded pixels")
  }
  metadata_check <- system2(
    "/usr/bin/sips",
    args = c("-g", "profile", "-g", "dpiWidth", "-g", "dpiHeight", shQuote(path)),
    stdout = TRUE, stderr = TRUE
  )
  if (!any(grepl("profile: sRGB IEC61966-2.1", metadata_check, fixed = TRUE)) ||
      !any(grepl("dpiWidth: 600.000", metadata_check, fixed = TRUE)) ||
      !any(grepl("dpiHeight: 600.000", metadata_check, fixed = TRUE))) {
    stop("PNG profile/resolution validation failed")
  }
}

png_file <- file.path(OUTDIR, "SupFig3e_healthy_enriched_3of4_abundance_vertical.png")
grDevices::png(
  png_file,
  width = as.integer(round(WIDTH_IN * DPI)),
  height = as.integer(round(HEIGHT_IN * DPI)),
  units = "px", res = DPI, type = "cairo-png", bg = "white"
)
print(p)
grDevices::dev.off()
standardize_png(png_file)

ggsave(
  file.path(OUTDIR, "SupFig3e_healthy_enriched_3of4_abundance_vertical.pdf"),
  p, width = WIDTH_IN, height = HEIGHT_IN, device = cairo_pdf, bg = "white"
)
readr::write_csv(
  abundance_long %>%
    select(Cohort, RunID, Group, Species_full, Species, Abundance),
  file.path(OUTDIR, "per_sample_abundance_zero_complete.csv")
)
readr::write_csv(stats, file.path(OUTDIR, "pooled_descriptive_wilcoxon_stats.csv"))
readr::write_tsv(
  tibble(Input = names(INPUT_FILES), File = unname(INPUT_FILES), SHA256 = input_sha),
  file.path(OUTDIR, "input_sha256.tsv")
)

writeLines(
  c(
    "Supplementary Figure 3e vertical Healthy-versus-Cancer abundance plot (2026-08-27)",
    "Target set: eight species meeting BH q < 0.05 and log2FC(Cancer/Healthy) < -1 in exactly 3 of 4 cohorts.",
    "Sample roster: n=1,649 taxonomic samples (Healthy 745; Cancer 904).",
    "Per-sample missing target taxa were reconstructed as zero, matching the former Cancer-enriched Sup. Fig. 3e workflow.",
    "Plot: vertical boxplots plus raw points; y-axis log10; zeros displayed at one-fifth the smallest non-zero value within species.",
    "Statistics: pooled Wilcoxon rank-sum tests on raw abundance; BH correction across 8 species; descriptive only.",
    "Formal selection evidence remains the four frozen per-cohort differential-abundance analyses.",
    paste0("PNG: ", round(WIDTH_IN * DPI), " x ", round(HEIGHT_IN * DPI),
           " px; 600 dpi; 8-bit RGB; embedded sRGB profile.")
  ),
  file.path(OUTDIR, "VALIDATION_REPORT.txt"), useBytes = TRUE
)
capture.output(sessionInfo(), file = file.path(OUTDIR, "sessionInfo.txt"))

manifest_files <- list.files(OUTDIR, full.names = TRUE)
manifest_files <- manifest_files[basename(manifest_files) != "output_sha256.tsv"]
output_sha <- vapply(
  manifest_files, digest::digest, character(1),
  algo = "sha256", file = TRUE, serialize = FALSE
)
readr::write_tsv(
  tibble(File = basename(manifest_files), SHA256 = output_sha),
  file.path(OUTDIR, "output_sha256.tsv")
)

cat("Vertical Healthy-versus-Cancer abundance panel complete.\n")
cat(png_file, "\n")
print(stats %>%
        select(Species, n_Healthy, n_Cancer, pooled_log2FC,
               pooled_Wilcoxon_p, pooled_Wilcoxon_BH_q))
