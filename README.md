# Analysis code — sarcosine, the gut microbiome, and cancer immunotherapy

Custom analysis code for the public-data analyses in the manuscript
**"[MANUSCRIPT TITLE — to be filled in]"** (authors, year).

This repository contains the R scripts that reproduce the public-data figures
(Figure 1–3 and Supplementary Figures 1–15). Each top-level folder corresponds
to one analysis and to the manuscript figure(s) it produces.

> **Raw data are not included here.** They are public and must be obtained from
> the accessions listed below. Each folder's `README.md` documents its inputs and
> run order. Set the R working directory to the relevant analysis folder before
> running; all paths in the scripts are relative to it.

## Repository structure

| Folder                                | Manuscript figure(s)                                                                                        | Data source (accession)                                                                       |
| ------------------------------------- | ----------------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------- |
| `Fig1b_MTBLS10232_faecal_metabolome/` | Figure 1b — faecal sarcosine elevated in CRC                                                                | MetaboLights **MTBLS10232**                                                                   |
| `Fig1c_Hakimi_ccRCC/`                 | Figure 1c — tumour-tissue sarcosine elevated in ccRCC                                                       | Hakimi et al. _Cancer Cell_ 2016 supplementary table; GEO **GSE74734**                        |
| `Fig1dg_TIGER_melanoma/`              | Figure 1d–g + Supplementary Figure 1b–d — host sarcosine-metabolic gene expression vs ICI response/survival | RNA-seq cohort **PRJEB23709** (TIGER portal, http://tiger.canceromics.org)                    |
| `Fig2_CRC_WGS/`                       | Figure 2 + CRC Supplementary Figures — gut-microbial sarcosine-degradation capacity in CRC                  | WGS BioProjects **PRJEB6070, PRJEB10878, PRJEB27928, PRJNA429097** (via the HGMT database)    |
| `Fig3_NSCLC_WGS/`                     | Figure 3 + NSCLC Supplementary Figures — sarcosine-degrading gut microbiome marks ICI response in NSCLC     | WGS BioProjects **PRJNA751792, PRJNA1023797, PRJEB22863, PRJEB26531** (via the HGMT database) |
| `SupFig13_NSCLC_16S/`                 | Supplementary Figure 13 — independent 16S corroboration (NSCLC)                                             | 16S BioProject **PRJEB26531** (via the HGMT database)                                         |
| `SupFig14-15_Liang_multiomics/`       | Supplementary Figures 14–15 — matched 16S + faecal metabolome                                               | 16S BioProject **PRJNA763023** + faecal metabolome MetaboLights **MTBLS10232**                |

## Not in this repository

- **Figure 4–5 (mouse experiments):** analysed in GraphPad Prism (not R); source data reported with the manuscript.
- **Figure 1a / Supplementary Figure 1a (pathway schematics):** drawn from KEGG pathway maps (map00260, map00670, map00330) and assembled in BioRender — no code.

## Software

Analyses were run in **R (≥ 4.5)** with Bioconductor; key packages are cited in the
manuscript Methods and listed at the top of each script (e.g. MetaPhlAn4-derived
profiles, MaAsLin2, ANCOM-BC2, metafor, vegan, limma, DADA2, phyloseq, survival).
Each script lists the specific package versions used.

## How to use

1. Pick the folder for the figure you want to reproduce.
2. Read that folder's `README.md` (run order + inputs).
3. Download the listed data from its accession(s) into the layout described there.
4. Set the R working directory to the analysis folder and run the scripts in order.

## License

No license file is included yet. To allow reuse, consider adding a `LICENSE`
(e.g. MIT for permissive reuse). Until then, all rights are reserved by default.

## Citation

If you use this code, please cite the manuscript:
_[Authors]. [Title]. [Journal] [Year]. [DOI]._ — to be completed on acceptance.
