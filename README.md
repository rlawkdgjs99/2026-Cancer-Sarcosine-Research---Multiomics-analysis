# Sarcosine–cancer multi-omics analyses

This repository contains analysis code and selected machine-readable result tables for the public-cohort components of a study of sarcosine metabolism, tumor immunity, and cancer-associated microbiomes.

The repository is organized by **cohort and analysis**, not by manuscript panel number. This keeps the code stable when the manuscript layout changes.

## Central result adopted for the manuscript

Two independent bulk-tumor cohorts were used for discovery:

- TIGER melanoma (PRJEB23709): pretreatment tumors with high versus low RNA-defined sarcosine-degradation score, `mean(z(SARDH), z(PIPOX))`.
- Matched ccRCC metabolome–RNA-seq cohort: tumors with low versus high measured sarcosine intensity using the prespecified within-cohort median split.

The same four exact MSigDB programs were enriched in the target direction in both cohorts:

`APC cross-presentation → TCR signaling → TNFR2-related noncanonical NF-κB expression program → IFNG response`

The arrows denote a biologically ordered **candidate pathway-enrichment axis**. The shared result is a bulk-tumor program-level association, not proof of temporal causality, receptor engagement, phosphorylation, or signaling within one cell type. The complete exposure-linked serial route did not replicate in ccRCC, and common NFKB2/RELB transcription-factor activation was not demonstrated.

Prespecified RNA-only validation used TCGA NSCLC primary tumors (LUAD/LUSC), where all four exact pathways replicated. This NSCLC result, together with the melanoma–ccRCC discovery analysis, is the immune-pathway analysis selected for the manuscript.

## Repository map

| Directory | Cohort / data type | Purpose |
|---|---|---|
| `analyses/melanoma_ccrcc_shared_immune_axis/` | TIGER melanoma + matched ccRCC bulk tumor | Exact-set intersection, leading-edge overlap, gene-disjoint route models, TF diagnostics, and final narrowing |
| `analyses/tcga_nsclc_immune_axis_validation/` | TCGA LUAD/LUSC bulk tumor | Prespecified RNA-only validation of the shared immune axis |
| `analyses/tiger_melanoma_immunotherapy_rnaseq/` | PRJEB23709 bulk RNA-seq | ICI response/survival analyses and upstream immune-pathway analyses |
| `analyses/ccrcc_matched_metabolome_rnaseq/` | matched tumor metabolome + RNA-seq | Sarcosine grouping, GSEA, deconvolution, TF inference, continuous/quartile sensitivities, and pathway narrowing |
| `analyses/gse207422_nsclc_scrna/` | GSE207422 single-cell RNA-seq | Sarcosine functional-score analyses with frozen lineage annotation and patient-aware lineage summaries |
| `analyses/crc_fecal_metabolomics/` | MTBLS10232 | CRC-versus-control fecal metabolomics |
| `analyses/ccrcc_tumor_metabolomics/` | Hakimi et al. ccRCC tissue metabolomics | Paired tumor-versus-normal sarcosine analysis |
| `analyses/crc_wgs_microbiome/` | four CRC shotgun-metagenomic cohorts | Microbial sarcosine-degradation capacity and cross-cohort consistency |
| `analyses/nsclc_wgs_microbiome/` | four NSCLC ICI shotgun-metagenomic projects | Microbial sarcosine functions, taxa, ICI response, and survival |
| `analyses/nsclc_16s_microbiome/` | PRJEB26531 16S | Independent NSCLC microbiome corroboration |
| `analyses/crc_matched_16s_metabolomics/` | PRJNA763023 + MTBLS10232 | Matched genus-level microbiome–sarcosine integration |

## Data and reproducibility

Raw public data and large processed matrices are not redistributed. See [DATASETS.md](DATASETS.md) for accessions and specimen definitions, [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for the execution contract, and [ANALYSIS_STATUS.md](ANALYSIS_STATUS.md) for the evidence level of each module.

Selected result tables are included where they are needed to audit the reported conclusions. Machine-specific source paths have been removed from the public manifests while retaining file hashes and dataset identities.

## Important interpretation notes

- “TNFR2-related noncanonical NF-κB” is the name of the exact Reactome expression program used in the analysis. Enrichment of that set does not by itself demonstrate TNFR2 ligation or NFKB2/RELB protein activation.
- Melanoma and ccRCC use different sarcosine-related exposures: an RNA-derived degradation score versus measured metabolite intensity.
- All discovery and TCGA validation results are from bulk tumor tissue unless a module is explicitly labeled single-cell.
- No license is asserted here; users should obtain permission before redistributing or reusing code beyond ordinary scholarly inspection.
