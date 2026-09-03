## 00_setup.R — shared configuration for the PRJEB26531 Healthy vs NSCLC 16S analysis
## Project : HGMT 16S, Healthy vs Cancer, NSCLC cohort PRJEB26531 (Korea)
## Data    : HGMT portal export — genus-level 16S relative abundance (TSS, compositional)
## Author  : analysis pipeline (reproducible); R 4.6.0
## NOTE    : Each numbered script sources this file, then runs independently.

## ===== The ONLY path to edit if the project is relocated =====
BASE_DIR <- "."

## ----- Derived paths -----
ANALYSIS_DIR <- file.path(BASE_DIR, "analysis")
DATA_DIR     <- BASE_DIR                       # raw HGMT files live directly here
DERIVED_DIR  <- file.path(ANALYSIS_DIR, "data_derived")
RESULTS_DIR  <- file.path(ANALYSIS_DIR, "results")
FIG_DIR      <- file.path(ANALYSIS_DIR, "figures")
LOG_DIR      <- file.path(ANALYSIS_DIR, "logs")
for (d in c(DERIVED_DIR, RESULTS_DIR, FIG_DIR, LOG_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ----- Raw input files (exact names as downloaded from HGMT) -----
BACTERIA_FILE <- file.path(DATA_DIR, "Bacteria_1781090611263.txt")
META_FILE     <- file.path(DATA_DIR, "selected_project_1781090596057.txt")

## ----- Reproducibility -----
set.seed(42)                  # fixed, arbitrary — documented per faithful-coder
options(stringsAsFactors = FALSE)

## ----- Documented analysis constants / thresholds -----
PREVALENCE_MIN <- 0.10        # genus must be present (>0) in >=10% of samples for DA (HGMT paper)
ALPHA_SIG      <- 0.05        # significance level
N_PERM         <- 999         # PERMANOVA / permutation tests (HGMT paper)
GROUP_LEVELS   <- c("Health", "NSCLC")  # Health = reference; positive effect => enriched in NSCLC

suppressPackageStartupMessages({
  library(data.table)
  library(vegan)
  library(ggplot2)
})

## ----- Shared plotting aesthetics -----
theme_hvc    <- theme_bw(base_size = 16) + theme(panel.grid.minor = element_blank())
GROUP_COLORS <- c(Health = "#1B7837", NSCLC = "#B2182B")

## ggplot2's default PNG backend (ragg) fails on this non-ASCII (Korean) project path;
## cairo works, so route all PNGs through it. PDFs use grDevices pdf() which also works.
save_png <- function(plot, file, width, height, dpi = 200) {
  grDevices::png(file, width = width, height = height, units = "in", res = dpi, type = "cairo")
  print(plot); grDevices::dev.off(); invisible(file)
}

cat("[00_setup] BASE_DIR =", BASE_DIR, "\n")
cat("[00_setup] versions: data.table", as.character(packageVersion("data.table")),
    "| vegan", as.character(packageVersion("vegan")),
    "| Maaslin2", as.character(packageVersion("Maaslin2")),
    "| R", as.character(getRversion()), "\n")
