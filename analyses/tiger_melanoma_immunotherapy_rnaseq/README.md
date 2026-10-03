# TIGER melanoma immunotherapy RNA-seq

This module analyzes melanoma cohort PRJEB23709 from the TIGER portal. The core pipeline relates SARDH/PIPOX degradation and GNMT/DMGDH production programs to pretreatment ICI response and survival. Extended scripts test cytotoxicity and candidate upstream immune pathways.

## Core run order

1. `analysis/01_data_prep.R`
2. `analysis/02_response_analysis.R`
3. `analysis/03_survival_analysis.R`
4. `analysis/05_tiger_all_samples.R`
5. `analysis/06_generate_summary_outputs.R`

`analysis/00_setup.R` supplies shared plotting and path helpers. Scripts `07` and `08` provide gene-resolved cytotoxicity analyses. `analysis/shared_axis_upstream/` contains the degradation-group GSEA, deconvolution, TME integration, IFNG-axis, TF-activity, T-cell-context, focused SPI1/NF-κB, route-narrowing, and cell-context sensitivity analyses that feed the cross-cohort shared-axis module.

Raw expression and clinical tables are not distributed. The primary discovery set contains 73 pretreatment biopsies. All results are observational bulk-tumor associations.

## Figure 5 scripts (`results/<branch>/code/`)

The current Figure 5 and Supplementary Figures 20–24 use the 73 pretreatment tumors and the Balance score, (zSARDH + zPIPOX − zGNMT − zDMGDH)/2, with z scores from log₂(FPKM + 1); panels g–j use 16 matched pretreatment/on-treatment pairs. The scripts below are kept in their original relative layout, `results/<branch>/code/` (see [LAYOUT.md](../../LAYOUT.md)).

| Branch | Content |
|---|---|
| `Fig5_PRE73_Balance_HighLow_OS_26.09.22` | Overall survival of Balance-High versus Balance-Low tumors (Figure 5d) |
| `Fig5_PRE73_Balance_HighLow_CD8_CYT_26.09.22` | quanTIseq CD8⁺ T-cell percentage and cytolytic expression (CYT) by Balance group (Figure 5e,f) |
| `Fig5_PRE73_Balance_CD8_CYT_Correlations_26.09.22` | Balance versus CD8 and CYT (Figure 5e,f) |
| `Fig5_PRE73_FourGenes_CD8_CYT_26.09.22` | SARDH, PIPOX, GNMT and DMGDH versus CD8 and CYT (Supplementary Figure 23c) |
| `Fig5_PRE73_Ratio_Coordinates_CD8_CYT_26.09.22` | Log-expression ratio and two-gene coordinates versus CD8 and CYT (Supplementary Figure 23c) |
| `Fig5_Matched16_PRE_EDT_Immune_Genes_26.09.22` | quanTIseq (IOBR) and CYT for the 16 matched pairs: inputs of Figure 5g–j and Supplementary Figure 24 |

What these scripts need and what is not here:

- They read author-prepared workbooks (`73Pre.xlsx`, `16Matched.xlsx`) that tabulate the TIGER expression values, the scores and the clinical covariates for the 73 patients and the 16 pairs, and they compare their inputs with tables of earlier analyses of the project. Neither is distributed, so the scripts document the calculations but cannot be run from this repository alone.
- For the matched pairs, `02_quantiseq.R` documents the quanTIseq and CYT calculation (IOBR 2.2.3; CYT from GZMA and PRF1 TPM), which uses the same settings as the earlier 73-tumor calculation.
- The hazard ratio and confidence interval of Figure 5d (Prism log-rank estimate, as stated in the legend) and the repeated-measures ANOVA of Figure 5g and Supplementary Figure 24a,b were run in GraphPad Prism.
- Scripts for the data ellipses and centroid regions (Figure 5c,h; Supplementary Figure 21), the paired Hotelling T² tests (Figure 5i,j; Supplementary Figure 24c) and the ten-measure response and survival summaries (Supplementary Figure 22) are not part of this repository.
