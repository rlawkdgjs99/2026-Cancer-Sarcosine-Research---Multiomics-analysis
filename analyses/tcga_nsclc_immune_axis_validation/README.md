# TCGA NSCLC validation of the shared immune axis

This module tests the four pathways selected before the TCGA analysis: APC cross-presentation, TCR signaling, the Reactome TNFR2-related noncanonical NF-κB expression program, and IFNG response. The manuscript validation cohort is NSCLC, combining LUAD/LUSC primary tumors and splitting them at the median RNA-defined degradation score `mean(z(SARDH), z(PIPOX))`.

## Result

- NSCLC: four of four exact pathways have positive NES and collection-wide BH q<0.05. This is the validation result selected for the manuscript.
- The compact public-release table audit passes all prespecified scope and result checks; see `results/logs/public_release_table_check.txt`.

NSCLC is the full exact-set validation used in the manuscript.

## Run order

1. `scripts/00_download_tcga_xena.R`
2. `scripts/01_prepare_sarcosine_scores.R`
3. `scripts/02_validate_shared_immune_axis.R`
4. `scripts/03_verify_results.R`
Run scripts from this module directory. Raw Xena downloads and large RDS objects are not included. The selected audit tables and reports are under `results/`.

The full verification script is intended to run after steps 1–2 regenerate the local RDS inputs and step 3 regenerates the sample-level score table. That sample-level table is not redistributed in this compact release.

The ordered arrows used in pathway schematics are a candidate topology. The adjacent-node models support coordinated expression after adjustment; they do not establish temporal causality or same-cell signaling.
