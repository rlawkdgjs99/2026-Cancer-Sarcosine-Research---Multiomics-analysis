# NSCLC 16S microbiome corroboration

This module analyzes 16S data from PRJEB26531 for healthy-versus-NSCLC genus-level alpha diversity, beta diversity, differential abundance, and selected taxa.

Run the scripts under `analysis/` in numeric order after placing the HGMT abundance and metadata exports at the expected analysis root. `00_setup.R` defines paths and reproducibility settings; downstream scripts consume the prepared RDS and earlier result tables. `06_query_genera.R` looks up the selected genera in the differential-abundance tables. Raw data are not included.
