# NSCLC ICI shotgun-metagenomic sarcosine analysis

This module analyzes PRJNA751792, PRJNA1023797, PRJEB22863, and PRJEB26531 for microbiome diversity, species abundance, microbial sarcosine functions, ICI response, and survival.

Run `sarcosine_KO_derivation.R`, each cohort's analysis scripts, and then the pooled and cross-cohort scripts. Raw sequences and HGMT taxonomic/KO profiles are not included.

Important provenance boundary: PRJNA751792 and PRJNA1023797 contain 283 overlapping biological sample names, including 45 discordant response labels. Until a BioSample-level reconciliation is completed, pooled cross-project estimates involving both projects are sensitivity results and project-run counts must not be described as independent-patient counts.
