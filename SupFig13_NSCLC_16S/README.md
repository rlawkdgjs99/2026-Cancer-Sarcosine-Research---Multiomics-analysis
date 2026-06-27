# 16S rRNA corroboration in the Korean NSCLC cohort (HGMT, PRJEB26531)

R code for the independent 16S rRNA gene survey of healthy vs. non-small-cell lung
cancer (NSCLC) gut microbiomes in a Korean cohort, used as an orthogonal corroboration
of the metagenomic findings. Analyses cover genus-level alpha/beta diversity and
differential abundance (Health vs. NSCLC).

## Figure -> script mapping

These scripts produce **Supplementary Figure 13** (independent 16S rRNA corroboration
in the Korean NSCLC cohort).

- `analysis/00_setup.R` and `analysis/01_load_prepare.R` are shared upstream
  (configuration + data loading/preparation).
- `analysis/02_alpha_diversity.R`, `analysis/03_beta_diversity.R`,
  `analysis/04_differential_abundance.R`, `analysis/05_DA_plots.R`,
  `analysis/regenerate_alpha_pvals.R`, and `analysis/07_suppfig_lachnospira_boxplot.R`
  produce the panels.

## Data source / accessions

- 16S BioProject **PRJEB26531** (via the HGMT database).

**Raw data are NOT included in this deposit.** The 16S relative-abundance export and the
sample metadata must be obtained from the accession above (PRJEB26531, via the HGMT
portal). The pipeline expects the two HGMT export files referenced in `00_setup.R`
(`Bacteria_1781090611263.txt` and `selected_project_1781090596057.txt`) to be present at
the analysis-folder root before `01_load_prepare.R` is run.

## Run order

Run the scripts in this order:

1. `analysis/00_setup.R` — shared configuration (paths, constants, plotting theme).
   Sourced automatically by the numbered scripts below; no need to run it on its own.
2. `analysis/01_load_prepare.R` — load HGMT export, subset to 16S, build the genus-level
   TSS relative-abundance matrix (writes `analysis/data_derived/prepared_16S.rds`).
3. `analysis/02_alpha_diversity.R` — within-sample (alpha) diversity, Health vs. NSCLC.
4. `analysis/regenerate_alpha_pvals.R` — re-render the 4-metric alpha boxplot with
   per-facet Wilcoxon P labels (read-only; reuses `02_alpha_diversity.R` outputs, recomputes
   no statistics). Run after `02_alpha_diversity.R`.
5. `analysis/03_beta_diversity.R` — between-sample (beta) diversity / PERMANOVA.
6. `analysis/04_differential_abundance.R` — genus-level differential abundance (MaAsLin2).
7. `analysis/05_DA_plots.R` — differential-abundance plots (reuses `04` outputs).
8. `analysis/07_suppfig_lachnospira_boxplot.R` — supplementary _Lachnospira_ boxplot
   (reuses prepared data and `04` outputs; recomputes no statistics).

Scripts 2, 3, and 6 each source `00_setup.R` and read the prepared `.rds`; scripts 4, 5,
7, and `regenerate_alpha_pvals.R` consume outputs written by earlier steps, so keep the
order above.

## Path convention

**Set the R working directory to the analysis-folder root before running; all paths are
relative to it.** The original absolute paths have been replaced with portable relative
paths under this convention (`BASE_DIR`/`BASE` is `.`; the `source()` calls and
derived paths resolve from the working directory). For example, from R:

```r
setwd("/path/to/HGMT_16S_Healthy_vs_Cancer_NSCLC_PRJEB26531")
source("analysis/01_load_prepare.R")
```

## Environment

R 4.6.0 with `data.table`, `vegan`, `ggplot2`, and `Maaslin2` (see `00_setup.R` for the
versions printed at setup). A random seed of 42 is fixed in `00_setup.R` for
reproducibility.
