# CRC shotgun-metagenomic sarcosine analysis

This module analyzes four CRC stool WGS cohorts: PRJEB6070, PRJEB10878, PRJEB27928, and PRJNA429097. It compares healthy and CRC microbiomes, quantifies microbial sarcosine-production/degradation functions, and evaluates cross-cohort consistency.

Run each cohort's `analysis_healthy_vs_cancer.R`, then the scripts in `Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/`. Regeneration scripts consume the saved cohort and pooled tables. Raw sequences and HGMT taxonomic/KO profiles are not included.

The numbered scripts `100`–`129` in the integrated folder produce the pooled and cohort-specific panels of the current manuscript (for example volcano, forest, heatmap, lollipop and species panels); scripts named `verify_*` or `*_verify_*` check the inputs or outputs of the script they follow. Input matrices are written by `사용데이터_모음/scripts/export_analysis_used_inputs.py` (`사용데이터_모음` = "collection of the data actually used"). Many of these scripts find their inputs relative to their own location; see `../../LAYOUT.md`. A panel-level overview is in `../../FIGURE_MAP.md`.
