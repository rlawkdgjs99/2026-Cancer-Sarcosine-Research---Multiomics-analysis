# Gut-microbial sarcosine metabolism in colorectal cancer — 16S + faecal metabolome analysis code

R scripts for the matched 16S rRNA amplicon + faecal metabolome analysis of the
Liang/Ma colorectal cancer (CRC) cohort. The pipeline reconstructs the faecal
microbiome from paired-end 16S reads (DADA2), characterises microbiome structure
in CRC vs. healthy controls, and integrates genus-level abundances of
sarcosine-degrading / -producing commensals with faecal sarcosine intensity.

## Figure → script mapping

Supplementary Figures 14-15 (matched 16S + faecal metabolome).
pe_pipeline:

- `02_dada2_pe.R` — DADA2 (paired-end denoising, ASV inference, Silva v138.1
  taxonomy → genus table)
- `03_downstream_analysis.R` — microbiome structure (alpha/beta diversity,
  composition, differential abundance: CRC vs. healthy)
- `04_sarcosine_analysis.R` — sarcosine integration (faecal sarcosine by group;
  degrader/producer abundance vs. sarcosine correlation; degrader–sarcosine–cancer triangle)
- `09_genus_sarcosine_correlation_heatmap.R` — genus × sarcosine Spearman
  correlation heatmap (Overall / within-Healthy / within-CRC)

## Data source / accessions

- 16S rRNA amplicon: BioProject **PRJNA763023**
- Faecal metabolome: MetaboLights **MTBLS10232**

Raw data are **NOT** included in this repository. They must be obtained from the
accessions above (16S FASTQs from PRJNA763023; the negative-mode metabolite
assignment file `m_MTBLS10232_LC-MS_negative_reverse-phase_metabolite_profiling_v2_maf.tsv`
and associated metadata from MTBLS10232).

## Run order

Run the scripts in numeric order:

1. `02_dada2_pe.R` — produces `genus_table.tsv` (and ASV / tracking outputs).
   Requires cutadapt-trimmed paired-end FASTQs and the Silva v138.1 training set;
   driven by the `WORK`, `BASE`, and `THREADS` environment variables (see the
   header of the script).
2. `03_downstream_analysis.R` — microbiome structure; consumes `genus_table.tsv`
   and `pe_samples.tsv`.
3. `04_sarcosine_analysis.R` — sarcosine integration; consumes `genus_table.tsv`,
   `pe_samples.tsv`, and the MTBLS10232 negative-mode MAF. Writes
   `sarcosine/merged_per_sample.tsv`.
4. `09_genus_sarcosine_correlation_heatmap.R` — heatmap; consumes
   `sarcosine/merged_per_sample.tsv` produced by step 3.

## Path convention

Set the R working directory to the **analysis-folder root** before running; all
paths in the scripts are relative to it. In these scripts:

- `03_downstream_analysis.R`: `BASE <- "16S_데이터_분석/raw_download"`
- `04_sarcosine_analysis.R`: `BASE <- "."` (the analysis-folder root)
- `09_genus_sarcosine_correlation_heatmap.R`: `OUT <- "16S_데이터_분석/raw_download/pe_pipeline/sarcosine"`

`02_dada2_pe.R` resolves its inputs/outputs from the `WORK` and `BASE`
environment variables rather than relative paths, so it is unaffected by the
working directory and contains no absolute paths.

The directory substructure under `code_for_publication/`
(`16S_데이터_분석/raw_download/pe_pipeline/`) mirrors the original analysis-folder
layout so the relative paths resolve correctly when the working directory is the
analysis-folder root.
