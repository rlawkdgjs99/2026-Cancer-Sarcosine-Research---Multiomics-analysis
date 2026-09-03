# Matched CRC 16S–fecal-metabolome integration

This module integrates PRJNA763023 16S profiles with MTBLS10232 fecal metabolomics. It performs DADA2 denoising, CRC-versus-healthy microbiome analyses, and genus-level associations with fecal sarcosine intensity.

Run `02_dada2_pe.R`, `03_downstream_analysis.R`, `04_sarcosine_analysis.R`, and `09_genus_sarcosine_correlation_heatmap.R` in numeric order. FASTQs, the SILVA training set, metabolomics MAF files, and sample metadata are not included.
