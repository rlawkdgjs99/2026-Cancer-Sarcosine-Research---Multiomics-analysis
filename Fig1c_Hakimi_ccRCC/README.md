# Sarcosine in clear cell RCC: tumour vs normal (Hakimi et al., Cancer Cell 2016)

R code accompanying the analysis of tissue sarcosine abundance in clear cell
renal cell carcinoma (ccRCC), based on the published metabolomics supplementary
data of Hakimi et al.

## Figure -> script mapping

- **Figure 1c** (tumour-tissue sarcosine elevated in ccRCC, 138 matched pairs).
  Single self-contained script reads the Hakimi et al. 2016 supplementary xlsx
  and writes the boxplot.
  - `R_scripts/01_sarcosine_tumor_vs_normal.R`

## Data source / accessions

- Hakimi et al. (Cancer Cell 2016) supplementary table (`mmc2.xlsx`); GEO **GSE74734**.

The script reads the supplementary table file
`1-s2.0-S1535610815004687-mmc2.xlsx` (sheets `Median Normalized` and `Raw Data`)
expected at the analysis-folder root.

**Raw data are NOT included in this deposit.** Obtain the supplementary table
and associated data from the accessions above (the Hakimi et al. 2016 Cancer Cell
supplementary materials and GEO GSE74734) and place
`1-s2.0-S1535610815004687-mmc2.xlsx` at the analysis-folder root before running.

## Run order

1. `R_scripts/01_sarcosine_tumor_vs_normal.R`

Outputs (a stats CSV and two PNG figures) are written to `Outputs/`, created by
the script if it does not already exist.

## Conventions

- **Set the R working directory to the analysis-folder root before running;
  paths in the scripts are relative to it.** In each script `BASE <- "."` denotes
  that root, and all input/output paths are derived from it via `file.path(BASE, ...)`.

## Requirements

R with the packages used by the script: `readxl`, `dplyr`, `tidyr`, `ggplot2`.
