# CRC fecal metabolomics

This module reproduces the CRC-versus-healthy fecal-metabolomics analysis from MetaboLights MTBLS10232, including sarcosine, ROC and median-split summaries, an anatomical-site check, differential tables, volcano summaries, and quality-control outputs.

Pipeline (from the analysis-data root): `R_scripts/01_install_and_inventory.R` and `02_download_and_parse.R` (retrieve and parse the deposited MAF files), `03_final_merge.R`, `04_parsing_fixes.R` (parsing corrections applied before the merge is used downstream), `05_characterize_matrices.R`, `13_consolidated_regeneration.R` (normalization and regeneration of the analysis tables), `101_ROC_stackedbar_26.07.27.R` (ROC and stacked bars), `15_sarcosine_anatomical_site_consistency_26.08.23.R`, and `14_dissertation_plots.R`. Required inputs are the positive- and negative-mode MAF tables, merged metadata, and the source-publication supplementary table used for replication. Raw data are not included.
