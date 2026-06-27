# Gut-microbial sarcosine-degradation capacity in colorectal cancer (Healthy vs Cancer, 4 integrated CRC WGS cohorts)

R code deposit for the colorectal-cancer (CRC) shotgun-metagenome (WGS) analysis of this study.
The analysis compares Healthy vs Cancer gut microbiomes across four CRC BioProjects and
characterises microbial sarcosine metabolism — in particular the sarcosine-**degradation**
capacity (sarcosine oxidase / SOX panel and related KEGG orthologues) — both pooled across
cohorts and per cohort.

## Figure → script mapping

**Figure 2 (gut-microbial sarcosine-degradation capacity in CRC, 4 integrated cohorts) + its
Supplementary Figures.** Pipeline:
`integrated_analysis_pooled.R` (upstream stats / CSVs) →
`regenerate_*.R` (figure PNGs) +
`cross_cohort_*.R` +
`supp_*boxplot.R` +
per-cohort scripts.

Pooled (4-cohort integrated) scripts live in
`Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/`:

| Script                                                       | Role                                                                                                                                          |
| ------------------------------------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------- |
| `integrated_analysis_pooled.R`                               | **Upstream**: pools the 4 cohorts and computes the integrated statistics / result CSVs that the figure scripts below consume. Run this first. |
| `regenerate_plots.R`                                         | Pooled differential-abundance / sarcosine-bacteria / pathway-balance figure PNGs (from the upstream CSVs).                                    |
| `regenerate_diversity_plots.R`                               | Pooled alpha-diversity and beta-diversity (PCoA) figures.                                                                                     |
| `regenerate_volcano_KO_paper.R`                              | Pooled KEGG-KO differential-abundance volcano (sarcosine KOs highlighted).                                                                    |
| `regenerate_correlation_heatmap_paper.R`                     | Pooled species × sarcosine-KO correlation heatmap.                                                                                            |
| `regenerate_sarcosine_boxplot_log10.R`                       | Pooled sarcosine-metabolism boxplot (log10 scale).                                                                                            |
| `scatter_sarcosine_species.R` / `regenerate_scatter_plots.R` | Sarcosine vs species scatter plots.                                                                                                           |
| `cross_cohort_consistency.R`                                 | Cross-cohort consistency of the bacterial / KO signals.                                                                                       |
| `cross_cohort_degradation_consistency.R`                     | Cross-cohort consistency of the sarcosine-**degradation** signal specifically.                                                                |
| `supp_deg_assoc_5species_boxplot.R`                          | Supplementary: degradation-associated 5-species boxplots.                                                                                     |
| `supp_CRC_enriched_3species_boxplot.R`                       | Supplementary: CRC-enriched 3-species boxplots.                                                                                               |
| `supp_hhathewayi_boxplot.R`                                  | Supplementary: _Hungatella hathewayi_ boxplot.                                                                                                |
| `regenerate_diversity_plots_percohort.R`                     | Per-cohort analogue of the pooled diversity figures (loops over the 4 cohorts).                                                               |
| `regenerate_volcano_KO_paper_percohort.R`                    | Per-cohort analogue of the pooled KO volcano (loops over the 4 cohorts).                                                                      |
| `regenerate_correlation_heatmap_paper_percohort.R`           | Per-cohort analogue of the pooled correlation heatmap (loops over the 4 cohorts).                                                             |

Per-cohort scripts live in each cohort folder
(`PRJEB6070_CRC_AdenomatousPolyps/scripts/`, `PRJEB10878_CRC/scripts/`,
`PRJEB27928_CRC/scripts/`, `PRJNA429097_CRC/scripts/`):

| Script                         | Role                                                                                                                                                                                                                                                        |
| ------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `analysis_healthy_vs_cancer.R` | **Per-cohort upstream**: the full Healthy-vs-Cancer pipeline for that one cohort (diversity, bacterial DA, KEGG-KO DA, sarcosine-enzyme analysis, species–sarcosine correlation, production-vs-degradation balance). Produces the cohort `results_*/` CSVs. |
| `regenerate_plots.R`           | Per-cohort figure PNGs re-rendered from that cohort's result CSVs.                                                                                                                                                                                          |

## Data source / accessions

WGS BioProjects **PRJEB6070, PRJEB10878, PRJEB27928, PRJNA429097**
(uniformly processed via the HGMT database).

> **Raw data are NOT included in this deposit.** The scripts read processed per-cohort inputs
> (taxonomic / KEGG-KO relative-abundance tables and the derived `results_*` CSVs). To reproduce
> from scratch, obtain the raw sequencing data from the four accessions above and process them via
> the HGMT pipeline to regenerate those inputs. The figure-regeneration scripts (`regenerate_*.R`)
> can be run from the saved result CSVs without re-running the upstream analyses.

## Run order

Set the R working directory to the **analysis-folder root** before running anything (see below).

1. **Per-cohort upstream** — for each cohort, run
   `<cohort>/scripts/analysis_healthy_vs_cancer.R`.
   This generates each cohort's `results_bacteria/`, `results_kegg/`, `results_sarcosine/` CSVs.
2. **Pooled upstream** —
   `Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/integrated_analysis_pooled.R`.
   This pools the four cohorts and writes the integrated result CSVs.
3. **Figures (run in any order after the upstream steps that produce their inputs):**
   - Pooled figures: `regenerate_plots.R`, `regenerate_diversity_plots.R`,
     `regenerate_volcano_KO_paper.R`, `regenerate_correlation_heatmap_paper.R`,
     `regenerate_sarcosine_boxplot_log10.R`, `scatter_sarcosine_species.R`,
     `regenerate_scatter_plots.R`
   - Cross-cohort: `cross_cohort_consistency.R`, `cross_cohort_degradation_consistency.R`
   - Supplementary boxplots: `supp_deg_assoc_5species_boxplot.R`,
     `supp_CRC_enriched_3species_boxplot.R`, `supp_hhathewayi_boxplot.R`
   - Per-cohort figure re-renders: `<cohort>/scripts/regenerate_plots.R`, and the per-cohort
     loop scripts `regenerate_diversity_plots_percohort.R`,
     `regenerate_volcano_KO_paper_percohort.R`,
     `regenerate_correlation_heatmap_paper_percohort.R`

## Path convention (important)

All paths in these scripts are **relative to the analysis-folder root** — the directory that
contains `Healthy_vs_Cancer_4_CRC_cohorts_integrated/` and the four `PRJ*` cohort folders.
**Set your R working directory to that root before running any script** (e.g.
`setwd("/path/to/analysis_root")`, or launch R from there). The original absolute
machine-specific paths have been replaced with these portable relative paths for this deposit;
no analysis logic, parameters, seeds, or statistics were changed.
