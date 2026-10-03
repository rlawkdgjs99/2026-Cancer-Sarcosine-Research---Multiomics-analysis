# NSCLC ICI shotgun-metagenomic sarcosine analysis

This module analyzes the NSCLC ICI shotgun-metagenomic projects PRJNA751792, PRJNA1023797 and PRJEB22863 for microbiome diversity, species abundance, microbial sarcosine functions, ICI response, and survival. `NSCLC_PRJEB26531/` holds an analysis of the PRJEB26531 WGS cohort that is not part of the manuscript analyses (the manuscript uses PRJEB26531 only for the 16S analysis in `../nsclc_16s_microbiome/`).

Run `sarcosine_KO_derivation.R`, each cohort's analysis scripts, and then the pooled and cross-cohort scripts. Raw sequences and HGMT taxonomic/KO profiles are not included.

The panels of the current manuscript are produced by the scripts in `pooled_analysis/R_scripts/` and by three September pipelines under `pooled_analysis/results/` (`CRC_Method_7KO_NSCLC_26.09.07`, `CRC_Method_Species_DA_26.09.09`, `NSCLC_7KO_CRC_style_panels_26.09.09`; each folder holds its scripts, a frozen analysis plan was kept with the original run, and `verify` scripts check the exported numbers). Input matrices are written by `사용데이터_모음/scripts/export_analysis_used_inputs.py` and checked by `validate_exported_inputs.py`. These pipelines find their inputs relative to their own location; see `../../LAYOUT.md`. A panel-level overview is in `../../FIGURE_MAP.md`.

Important provenance boundary: PRJNA751792 and PRJNA1023797 contain 283 overlapping biological sample names, including 45 discordant response labels. Until a BioSample-level reconciliation is completed, pooled cross-project estimates involving both projects are sensitivity results and project-run counts must not be described as independent-patient counts.
