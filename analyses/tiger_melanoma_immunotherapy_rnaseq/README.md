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
