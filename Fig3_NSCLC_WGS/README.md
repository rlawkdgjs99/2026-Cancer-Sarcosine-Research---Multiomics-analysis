# Sarcosine-degrading gut microbiome and ICI response in NSCLC — analysis code

R code that reproduces **Figure 3** (a sarcosine-degrading gut microbiome marks
immune-checkpoint-inhibitor [ICI] response in non-small-cell lung cancer, NSCLC)
and its Supplementary Figures, from shotgun metagenomic (WGS) profiles of four
ICI cohorts.

This is a public code deposit. **No raw or processed data are included** — only
the analysis scripts. See _Data source_ below for how to obtain the inputs.

---

## Data source / accessions

WGS BioProjects (obtained via the HGMT database):

| Accession    | Folder                 | Role       |
| ------------ | ---------------------- | ---------- |
| PRJNA751792  | `NSCLC_PRJNA751792`    | discovery  |
| PRJNA1023797 | `NSCLC_PRJNA1023797`   | discovery  |
| PRJEB22863   | `NSCLC_RCC_PRJEB22863` | discovery  |
| PRJEB26531   | `NSCLC_PRJEB26531`     | validation |

Raw sequence data and the MetaPhlAn4 species / KEGG-Ortholog (KO) functional
profiles are **NOT distributed here** and must be obtained from the accessions
above (via the HGMT database). Each per-cohort folder is expected to contain the
HGMT export files the scripts read (a `selected_project_*.txt` metadata table, a
`Bacteria_*.txt` MetaPhlAn4 species table, and `KO_relative_abundance.tsv`); the
per-cohort and pooled scripts write their own `results/` subfolders, which the
downstream scripts then read.

---

## Running convention (important)

**Set the R working directory to the analysis-folder root before running each
script; every path in the scripts is relative to that root.** For example:

```r
setwd("/path/to/HGMT_NSCLC_ICI_RvsNR_WGS")
source("NSCLC_PRJNA751792/R_scripts/run_analysis.R")
```

In the original working copy each script located itself from its own file path
(`commandArgs("--file=")`); for this deposit that self-location was replaced with
an explicit relative base (`.` = root, `pooled_analysis/`, `<cohort>/`) so the
scripts run correctly from the root regardless of where this `code_for_publication/`
folder sits. Nothing else (statistics, seeds, thresholds, parameters) was changed.

---

## Run order (4-tier pipeline)

Outputs of each tier are inputs to the next, so run them in order.

**Tier 0 — sarcosine KO set (run once)**

- `sarcosine_KO_derivation.R` — derives the sarcosine degradation/production KO
  set directly from KEGG via the `KEGGREST` package and writes
  `sarcosine_KO_set.csv` (consumed by the per-cohort and pooled sarcosine scripts).

**Tier 1 — per cohort** (run for each of the four cohort folders)

- `<cohort>/R_scripts/run_analysis.R` — alpha/beta diversity + species
  differential abundance (MaAsLin2), R vs NR.
- `<cohort>/R_scripts/run_sarcosine.R` — sarcosine pathway scores
  (degradation / production) and individual sarcosine-KO differential abundance.

**Tier 2 — pooled / cross-cohort** (`pooled_analysis/R_scripts/`)

- `pooled_taxa_diversity.R` — species DA + diversity meta-analysis.
- `pooled_sarcosine.R` — sarcosine KO + pathway-score meta-analysis.
- `pooled_sarcosine_combined.R`, `pooled_diversity_combined.R`,
  `pooled_diversity_plots.R` — combined/per-cohort pooled summaries and plots.
- `cross_cohort_degradation_venn.R` — cross-cohort degradation-associated
  species (Venn).
- `species_survival_km.R`, `score_survival_km.R` — Kaplan–Meier survival
  (per-species and per-score).

**Tier 3 — figure assembly / extra analyses** (analysis-folder root)

- `sarcosine_KO_derivation.R` (Tier 0, listed above).
- `sarcosine_species_assoc_dualfilter_bar.R` — sarcosine-associated bacteria,
  dual-filter (DA × correlation) bar.
- `sarcosine_species_scatter_grids.R`, `sarcosine_taxa_function_scatter.R` —
  taxa↔function scatter grids / integration.
- `replot_DA_volcano_byICIgroup.R`, `plot_DA_lollipop_byICIgroup.R` — DA volcano
  / lollipop by ICI-response group.
- `figure_hungatella_main.R` — _Hungatella hathewayi_ main figure panel.
- `suppl_hungatella_production.R`, `suppl_hungatella_abundance_boxplot.R` —
  _Hungatella_ supplementary panels.

Some Tier-2/Tier-3 scripts read intermediate CSVs written by earlier Tier-2
scripts into `pooled_analysis/results/` (e.g. `pooled_species_meta.csv`,
`sarcosine_species_assoc_meta.csv`), so run the pooled scripts before the figure
scripts.

---

## Figure → script mapping

**Figure 3** (sarcosine-degrading gut microbiome marks ICI response in NSCLC)
**and its Supplementary Figures.** 4-tier pipeline:

per-cohort `run_analysis.R` + `run_sarcosine.R` → `sarcosine_KO_derivation.R` →
`pooled_analysis/*` → figure scripts.

---

## Notes / caveats

- `sarcosine_KO_derivation.R` reads optional, pre-extracted per-cohort KO-count
  files from `/tmp/koc_*.txt` (lines 53–56). These are only used to annotate
  cohort presence and are guarded by `if (file.exists())`; the KO _set_ and its
  roles come entirely from live KEGG/KEGGREST queries, so the script runs and
  writes `sarcosine_KO_set.csv` even when those `/tmp` files are absent. The
  `/tmp` paths were left unchanged because they are scratch locations, not part
  of the analysis-folder data tree.
- The scripts render PNGs via an ASCII temp file and then `file.copy()` to the
  final destination (a workaround for the non-ASCII project path on macOS); this
  behaviour is preserved unchanged.
- Required R packages (from the scripts): `data.table`, `vegan`, `ggplot2`,
  `ggrepel`, `Maaslin2`, `metafor`, `survival`, `broom`, `scales`, `cowplot`,
  `KEGGREST`. Each script also writes a `sessionInfo_*.txt` for provenance.
- Not part of this deposit but present in the original analysis folder:
  `sarcosine_KO_verify.R` (an independent KO-set verification) and
  `final_audit.R` (a reproducibility audit). They are referenced in some script
  headers but are not required to reproduce the figures.
