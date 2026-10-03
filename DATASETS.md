# Public datasets and specimen definitions

| Analysis | Public source | Unit used in the analysis | Key definition |
|---|---|---|---|
| TIGER melanoma bulk RNA-seq | PRJEB23709; TIGER portal | 73 pretreatment tumor biopsies | Degradation-High/Low from the median of `mean(z(SARDH), z(PIPOX))` |
| Matched ccRCC metabolome–RNA-seq | Zenodo record 8063124; source article DOI 10.1038/s41588-024-01662-5 | 100 tumors with matched transcriptome and metabolomics | Sarcosine-Low/High from normalized tumor sarcosine intensity; median split |
| TCGA NSCLC | LUAD + LUSC through UCSC Xena | primary tumors, one sample per patient | RNA-defined Degradation-High/Low within NSCLC |
| NSCLC single-cell RNA-seq | GSE207422 | cells nested within patients; frozen lineage map | Production, degradation, and production/degradation functional-score analyses |
| CRC fecal metabolomics | MTBLS10232 | fecal metabolomic profiles | CRC versus healthy controls |
| ccRCC tissue metabolomics | Hakimi et al., Cancer Cell 2016 (doi:10.1016/j.ccell.2015.12.004), supplementary table S2; no repository accession | 138 matched tumor/normal pairs | Paired tumor-versus-normal sarcosine abundance |
| CRC shotgun metagenomics | PRJEB6070, PRJEB10878, PRJEB27928, PRJNA429097 | cohort-specific stool metagenomes | Healthy versus CRC; microbial sarcosine functions |
| NSCLC ICI shotgun metagenomics | PRJNA751792, PRJNA1023797, PRJEB22863 | project-specific stool metagenomes | ICI responder versus non-responder; microbial sarcosine functions |
| NSCLC 16S | PRJEB26531 (the PRJEB26531 WGS cohort is not used in the manuscript analyses) | 16S profiles | Healthy versus NSCLC |
| CRC matched 16S/metabolomics | PRJNA763023 + MTBLS10232 | matched genus abundance and fecal metabolite records | genus–sarcosine integration |

## Not distributed and not in this repository

- The NSCLC plasma sarcosine measurements (Figure 1b–d) come from a previously reported laboratory cohort (Lee et al., Drug Resist. Updat. 2024) and are not distributed in this repository; the Data availability statement of the manuscript describes how the deidentified values are provided.
- The mouse and T-cell experiments (Figure 7, Supplementary Figures 27–28) were analyzed in GraphPad Prism; there is no analysis code for them.
- The author-prepared workbooks `73Pre.xlsx` and `16Matched.xlsx` (TIGER-derived expression values, scores and clinical covariates of the 73 pretreatment patients and 16 matched pairs) are read by the Figure 5 scripts and are not distributed.
- Panels marked "Prism only" in [FIGURE_MAP.md](FIGURE_MAP.md) were drawn in GraphPad Prism on a separate computer; no script exists for them.

## Data boundaries

- The matched ccRCC metabolomics are normalized untargeted intensities, not absolute concentrations.
- The TIGER score is a transcriptional proxy for degradation capacity; it is not a metabolite measurement or flux assay.
- TCGA validation is RNA-only and therefore validates association with the prespecified transcriptional state, not measured sarcosine.
- The NSCLC WGS project metadata require an unresolved independent-unit audit: PRJNA751792 and PRJNA1023797 contain 283 overlapping biological sample names, including 45 discordant response labels. Until resolved, pooled cross-project estimates should be treated as sensitivity results rather than counts of independent patients.
- Raw and controlled-access data are not included. Users must follow the terms of the original repositories and publications.
