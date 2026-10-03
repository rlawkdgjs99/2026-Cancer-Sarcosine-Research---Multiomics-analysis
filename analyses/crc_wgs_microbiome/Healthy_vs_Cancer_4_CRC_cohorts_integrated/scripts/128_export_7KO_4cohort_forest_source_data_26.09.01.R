#!/usr/bin/env Rscript

# Export the exact plotting/statistics tables and sample-level KO values for
# the seven cohort-stratified CRC-WGS sarcosine-KO forest plots.

suppressPackageStartupMessages(library(digest))

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Could not determine script path")
script_path <- normalizePath(sub("^--file=", "", script_arg))
base_dir <- normalizePath(file.path(dirname(script_path), ".."))
analysis_root <- normalizePath(file.path(base_dir, ".."))

input_path <- file.path(
  base_dir, "results_integrated", "sarcosine",
  "sarcosine_KO_per_sample_pooled.csv"
)
dest_dir <- file.path(
  analysis_root, "사용데이터_모음",
  "CRC_WGS_7KO_4cohort_forest_source_data_26.09.01"
)
dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

five_dir <- file.path(
  base_dir, "results_integrated",
  "publication_style_plots_REMAINING_5KO_4COHORT_26.09.01"
)
creatinase_dir <- file.path(
  base_dir, "results_integrated",
  "publication_style_plots_CREATINASE_K08688_FOREST_26.09.01"
)
soxb_dir <- file.path(analysis_root, "Manuscript_Final_Panels_26.08.22")

info <- data.frame(
  KO = c("K00301", "K00302", "K00303", "K00305", "K00306", "K00315", "K08688"),
  Label = c("Sarcosine_oxidase", "soxA", "soxB", "soxG", "PIPOX", "DMGDH", "Creatinase"),
  Plot_source = c(
    file.path(five_dir, "K00301_Sarcosine_oxidase_4cohort_forest_plot_data.csv"),
    file.path(five_dir, "K00302_soxA_4cohort_forest_plot_data.csv"),
    file.path(soxb_dir, "Fig1j_soxB_cross_cohort_forest_statistics.csv"),
    file.path(five_dir, "K00305_soxG_4cohort_forest_plot_data.csv"),
    file.path(five_dir, "K00306_PIPOX_4cohort_forest_plot_data.csv"),
    file.path(five_dir, "K00315_DMGDH_4cohort_forest_plot_data.csv"),
    file.path(creatinase_dir, "Creatinase_K08688_cross_cohort_forest_statistics.csv")
  ),
  stringsAsFactors = FALSE
)

if (!file.exists(input_path)) stop("Missing pooled KO input")
if (!all(file.exists(info$Plot_source))) stop("Missing one or more plot/statistics sources")

dat <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("Run.ID", "Cohort", "Group", info$KO)
if (!all(required %in% names(dat))) stop("Pooled input lacks required columns")
if (nrow(dat) != 1647L || anyDuplicated(dat$Run.ID)) stop("Expected 1,647 unique Runs")
if (!setequal(unique(dat$Cohort), c("PRJEB6070", "PRJEB10878", "PRJEB27928", "PRJNA429097"))) {
  stop("Unexpected cohort labels")
}
if (!setequal(unique(dat$Group), c("Healthy", "Cancer"))) stop("Unexpected group labels")

manifest <- vector("list", nrow(info))
for (i in seq_len(nrow(info))) {
  ko <- info$KO[i]
  label <- info$Label[i]
  plot_dest <- file.path(
    dest_dir,
    paste0("CRC_WGS_4COHORT_", ko, "_", label, "_forest_plot_data.csv")
  )
  raw_dest <- file.path(
    dest_dir,
    paste0("CRC_WGS_4COHORT_", ko, "_", label, "_raw_values.csv")
  )

  if (!file.copy(info$Plot_source[i], plot_dest, overwrite = TRUE)) {
    stop("Could not copy plot source: ", ko)
  }
  if (!identical(
    digest(info$Plot_source[i], algo = "sha256", file = TRUE),
    digest(plot_dest, algo = "sha256", file = TRUE)
  )) stop("Plot-source hash mismatch: ", ko)

  values <- as.numeric(dat[[ko]])
  if (any(!is.finite(values)) || any(values < 0)) stop("Invalid KO abundance: ", ko)
  raw <- data.frame(
    Run.ID = dat$Run.ID,
    Cohort = dat$Cohort,
    Group = dat$Group,
    KO = ko,
    KO_relative_abundance = values,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  write.csv(raw, raw_dest, row.names = FALSE, quote = TRUE)

  manifest[[i]] <- data.frame(
    KO = ko,
    Label = label,
    Plot_data_file = basename(plot_dest),
    Plot_data_rows = nrow(read.csv(plot_dest, check.names = FALSE)),
    Plot_data_source = info$Plot_source[i],
    Plot_data_SHA256 = digest(plot_dest, algo = "sha256", file = TRUE),
    Raw_values_file = basename(raw_dest),
    Raw_values_rows = nrow(raw),
    Raw_input_source = input_path,
    Raw_input_SHA256 = digest(input_path, algo = "sha256", file = TRUE),
    stringsAsFactors = FALSE
  )
}

manifest <- do.call(rbind, manifest)
write.csv(
  manifest,
  file.path(dest_dir, "CRC_WGS_7KO_4cohort_forest_SOURCE_PROVENANCE.csv"),
  row.names = FALSE,
  quote = TRUE
)

readme <- c(
  "CRC WGS seven-KO four-cohort forest source data",
  "",
  "For each KO:",
  "- *_forest_plot_data.csv: exact byte-for-byte copy of the table used to draw the four-cohort forest.",
  "- *_raw_values.csv: all 1,647 Run-level KO relative-abundance values used to calculate the cohort estimates.",
  "",
  "Sign convention: r_rb = P(Healthy > CRC) - P(Healthy < CRC).",
  "The five newly drawn KO tables retain explicit all-zero/non-estimable annotations where applicable.",
  "K00303/soxB and K08688/Creatinase retain their original frozen four-cohort statistics tables."
)
writeLines(readme, file.path(dest_dir, "README.txt"), useBytes = TRUE)

for (f in list.files(dest_dir, full.names = TRUE)) {
  Sys.chmod(f, "0644")
  system2("/usr/bin/xattr", c("-c", shQuote(f)), stdout = FALSE, stderr = FALSE)
  system2("/usr/bin/chflags", c("nohidden", shQuote(f)), stdout = FALSE, stderr = FALSE)
}

cat("Created seven plot-data and seven raw-value CSV pairs in:\n", dest_dir, "\n", sep = "")
