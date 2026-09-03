# Melanoma–ccRCC shared immune-axis discovery

This module compares exact MSigDB enrichment results from two independently analyzed bulk-tumor cohorts:

- TIGER melanoma: Degradation-High versus Low in 73 pretreatment tumors.
- Matched ccRCC: measured Sarcosine-Low versus High in 100 tumors.

Positive harmonized values always denote these target states.

## Result

The shared four-program intersection is APC cross-presentation, TCR signaling, the Reactome TNFR2-related noncanonical NF-κB expression program, and IFNG response. Each exact set has positive NES and BH q<0.05 in both cohorts, with cross-cohort leading-edge overlap.

Gene-disjoint sample-level models retain APC → TCR → noncanonical NF-κB → IFNG as the direction-consistent downstream branch candidate. The complete target-state-linked serial route does not replicate in ccRCC, and the TF analysis does not demonstrate common NFKB2/RELB activation. IL-12/STAT4 and specific CD28 costimulation failed the two-cohort exact-set gate.

## Scripts

- `scripts/01_compare_common_pathways.R`: joins the same exact MSigDB sets, harmonizes direction, applies the dual-cohort BH gate, and calculates leading-edge intersections.
- `scripts/02_narrow_shared_immune_axis.R`: constructs gene-disjoint scores, tests group/continuous effects, adjacent nodes, serial products, branch competition, and TF diagnostics.

Selected result tables and reports are under `results/pathway_screen/` and `results/mechanism_narrowing/`. The full 4,649-row joined pathway table and expression matrices are not distributed; the scripts regenerate them from the source cohort outputs.

By default the scripts read documented filenames under `inputs/`. Every input can instead be supplied explicitly with the `SARCO_CCRCC_*`, `SARCO_TIGER_*`, and `SARCO_COMMON_PATHWAY_JOIN` environment variables defined at the top of the scripts. `SARCO_GSEA_R_LIB` is an optional package-library override, not a required machine-specific path.

This is a bulk-tumor pathway-enrichment result, not a cell-type-specific or causal signaling claim.
