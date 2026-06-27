# Sarcosine-metabolism gene expression vs ICI response/survival in melanoma (TIGER, PRJEB23709)

R analysis code for the melanoma RNA-seq arm of the sarcosine / ICI-responsiveness
study. The pipeline tests host sarcosine-metabolic enzyme genes
(degradation: `SARDH`, `PIPOX`; production: `GNMT`, `DMGDH`) and two module
z-scores against immune-checkpoint-inhibitor (ICI) response and overall survival.

## Figure -> script mapping

- **Figure 1d-g** (host sarcosine-metabolic gene expression vs ICI response /
  survival in melanoma) and **Supplementary Figure 1b-d** (production enzymes)
  are produced by the numbered pipeline `00` -> `06`.
- `06_manuscript_figures.R` writes the final figure PNGs into
  `results/manuscript_figures/` (`Main/` = degradation story = Fig 1d-g;
  `Supplementary/` = production arm = Sup Fig 1b-d).

## Data source / accessions

- RNA-seq cohort **PRJEB23709**, obtained via the **TIGER portal**
  (http://tiger.canceromics.org).
- **Raw data are NOT included** in this deposit. The expression and clinical
  tables must be obtained from the accession/portal above:
  - `Melanoma-PRJEB23709_ExpressionData.tsv`
  - `Melanoma-PRJEB23709_ClinicalData.tsv`

  `01_data_prep.R` expects both files in the analysis-folder root (i.e. the
  working directory; see "Running" below).

## Run order

Run the scripts in numbered order:

1. `analysis/01_data_prep.R` — load + verify expression/clinical TSVs, build the
   PRE-only (n=73) and all-91 analysis data frames, write
   `results/analysis_data_PRE.rds`, `results/analysis_data_all91.rds`, and
   `results/tables/analysis_data_PRE.csv`.
2. `analysis/02_response_analysis.R` — ICI responder (R) vs non-responder (NR)
   Wilcoxon tests + boxplots (PRE-only cohort).
3. `analysis/03_survival_analysis.R` — overall-survival KM / Cox by
   median-split expression (PRE-only cohort).
4. `analysis/05_tiger_all_samples.R` — TIGER-portal-style all-samples analysis
   (pre + on-treatment, repeats included) with repeated-measures-corrected p
   (mixed model / cluster-robust Cox) alongside the naive p.
5. `analysis/06_manuscript_figures.R` — curate the paper-ready figure set into
   `results/manuscript_figures/{Main,Supplementary}/`.

`analysis/00_setup.R` is not run directly; it is `source()`d by scripts 02, 03,
05, and 06 for shared colors, the publication ggplot theme, and save helpers.

Note: the sensitivity / outlier-removed step `04_sensitivity_outlier_removed.R`
is part of the full project but is **not** included in this deposit; scripts
`02`/`03`/`05`/`06` do not depend on its outputs.

## Running

Set the **R working directory to the analysis-folder root** before running, then
launch each script, e.g.:

```sh
# working directory = the analysis-folder root (the folder containing this
# code_for_publication/ directory, the results/ directory, and the two .tsv files)
Rscript code_for_publication/analysis/01_data_prep.R
Rscript code_for_publication/analysis/02_response_analysis.R
Rscript code_for_publication/analysis/03_survival_analysis.R
Rscript code_for_publication/analysis/05_tiger_all_samples.R
Rscript code_for_publication/analysis/06_manuscript_figures.R
```

All paths in the scripts are **relative to the analysis-folder root**:

- input TSVs: `Melanoma-PRJEB23709_ExpressionData.tsv`,
  `Melanoma-PRJEB23709_ClinicalData.tsv` (place in the root);
- `source("analysis/00_setup.R")` and outputs under `results/` resolve from the
  root.

The scripts as deposited live under `code_for_publication/analysis/`, but they
read `analysis/00_setup.R` and write `results/...` relative to the working
directory; either copy `analysis/00_setup.R` (and run outputs) into the root's
own `analysis/` and `results/` layout, or run the scripts from a checkout whose
root mirrors that layout. The intent is unchanged from the original project:
**cwd = analysis-folder root, paths relative to it.**
