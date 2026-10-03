# Matched CRC 16S–fecal-metabolome integration

This module integrates PRJNA763023 16S profiles with MTBLS10232 fecal metabolomics (matched subset, n = 308). It reprocesses the 16S reads (cutadapt, DADA2, SILVA v138.1), analyzes CRC-versus-healthy microbiome structure, and tests genus-level associations with fecal sarcosine intensity.

The scripts keep their original relative layout and are run from this module's folder (`16S_데이터_분석` = "16S data analysis"):

| Step | Script in `16S_데이터_분석/raw_download/pe_pipeline/` | Purpose |
|---|---|---|
| 1 | `run_pe_pipeline.sh` | cutadapt primer removal and DADA2 driver; set `CUTADAPT` if the executable is not on `PATH`, `SAMPLES`/`WORK` are optional |
| 2 | `02_dada2_pe.R` | DADA2 denoising and SILVA taxonomy (reads `WORK` and `BASE` from the environment) |
| 3 | `03_downstream_analysis.R` | genus-level diversity, PERMANOVA and differential abundance, healthy versus CRC |
| 4 | `04_sarcosine_analysis.R` | core analysis: merge with fecal sarcosine, degrader/producer genus abundance and correlations |
| 4b | `04b_sarcosine_analysis_eligens.R` | variant of step 4 (see its header); the original is not modified |
| 5–8 | `05_degrader_highlow.R`, `06_degrader_individual_highlow.R`, `07_individual_highlow_no_outliers.R`, `08_producer_individual_highlow.R` | High/Low comparisons of sarcosine; step 7 is an exploratory outlier-removal sensitivity analysis |
| 9 | `09_genus_sarcosine_correlation_heatmap.R` | genus–sarcosine Spearman correlations overall and within healthy/CRC |

Not included: FASTQ files (`16S_데이터_분석/raw_download/pe/`), `pe_samples.tsv`, the SILVA training set (`ref/`), the MetaboLights negative-mode MAF table (`Metabolomics_Data_MetaboLights/Result_files/`), and sample metadata. The scripts that build the genus-level effect table of Figure 1l / Supplementary Figure 5a and the Supplementary Figure 4 mixture model are not part of this release.
