# ccRCC paired tumor–normal metabolomics

This module tests sarcosine abundance in 138 matched ccRCC tumor/normal pairs using the Hakimi et al. (Cancer Cell 2016, doi:10.1016/j.ccell.2015.12.004) supplementary metabolomics workbook (`1-s2.0-S1535610815004687-mmc2.xlsx`, median-normalized sheet). No separate repository accession is used for these metabolomics data.

Run `R_scripts/01_sarcosine_tumor_vs_normal.R` from a directory containing the workbook; it writes statistics and plots to `Outputs/`. `R_scripts/101_ROC_stackedbar_26.07.27.R` produces the ROC and paired-classification panels. Raw and supplementary data are not redistributed.
