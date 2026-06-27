# Faecal metabolome analysis — CRC vs healthy (MTBLS10232)

R scripts that reproduce the faecal-metabolomics differential analysis and the
Figure 1b sarcosine boxplot (faecal sarcosine elevated in CRC vs healthy
controls), plus the supporting differential tables, volcano plots and QC
figures.

## Data source / accessions

- **Faecal metabolome:** MetaboLights **MTBLS10232** (LC-MS positive- and
  negative-mode reverse-phase metabolite profiling; `m_*_maf.tsv` MAF tables).
- Replication reference: Li 2024 supplementary Table S4
  (`Other_files/Supp info from the Article/mmc3.xlsx`).

**Raw data are NOT included in this deposit.** Download the MTBLS10232 study
files from MetaboLights and obtain the article's supplementary file separately,
then place them at the relative paths the scripts expect (see "Inputs" below)
before running.

## Run order

Run the scripts in numeric order, each from the same R session/working
directory:

1. `R_scripts/05_characterize_matrices.R`
   Loads the POS/NEG MAF intensity matrices, characterizes them, and writes
   `analysis/intensity_matrices.rds`.
2. `R_scripts/13_consolidated_regeneration.R`
   Re-runs the differential analysis (limma + Wilcoxon + MaAsLin2, age-adjusted),
   writes the differential tables `analysis/main/diff_results_POS.tsv` and
   `analysis/main/diff_results_NEG.tsv`, `significant_features_annotated.tsv`,
   the Li 2024 replication check, and volcano plots.
3. `R_scripts/14_dissertation_plots.R`
   Renders the dissertation-quality figures, including the headline Fig 1b
   boxplot `analysis/main/sarcosine_CRC_vs_CTRL.png`.

## Figure → script mapping

**Figure 1b** (faecal sarcosine elevated in CRC vs healthy):
`05_characterize_matrices.R` (builds `intensity_matrices.rds`)
→ `13_consolidated_regeneration.R` (differential tables `diff_results_*.tsv`)
→ `14_dissertation_plots.R` (writes the Fig 1b boxplot
`analysis/main/sarcosine_CRC_vs_CTRL.png`).

## Inputs the scripts expect (relative to the analysis-folder root)

- `MTBLS10232_inventory/m_MTBLS10232_LC-MS_positive_reverse-phase_metabolite_profiling_v2_maf.tsv`
- `MTBLS10232_inventory/m_MTBLS10232_LC-MS_negative_reverse-phase_metabolite_profiling_v2_maf.tsv`
- `metadata_merged_619samples.tsv`
- `Other_files/Supp info from the Article/mmc3.xlsx` (Li 2024 Table S4; replication step in script 13)

## Convention

Set the R working directory to the **analysis-folder root** before running each
script; all paths in the scripts are relative to it. Outputs are written to
`analysis/` (with `analysis/main/` and `analysis/supplementary/`
subdirectories) and `sessionInfo_*.txt` snapshots to `R_scripts/`.

> Note: these scripts live under `code_for_publication/R_scripts/` in this
> deposit, but they are written to be run with the working directory set to the
> analysis-folder root (one level above `code_for_publication/`), where the
> input data and `analysis/` output tree reside. They were copied byte-for-byte
> from the canonical analysis scripts; no path edits were required because the
> originals already use root-relative paths.
