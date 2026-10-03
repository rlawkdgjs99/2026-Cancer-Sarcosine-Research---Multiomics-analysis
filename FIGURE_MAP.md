# Figure map

Figure and panel numbers follow the current manuscript (Figures 1–7, Supplementary Figures 1–28). Several script names keep the earlier numbering used while the figures were being laid out (for example `Fig1h`, `Fig1j`, `SupFig3` in the CRC WGS scripts). The mapping below is therefore given at the level of figure groups, based on the analysis each script performs; confirm a specific panel against the script output before relying on a name.

Everything below uses public data unless stated (the Figure 5 scripts read author-prepared workbooks derived from the public TIGER data). Data are not distributed (see `DATASETS.md`); several pipelines assume the original workspace layout (see `LAYOUT.md`). Panels marked "Prism only" were drawn in GraphPad Prism outside this repository; no script is included for them.

## Figure 1 and Supplementary Figures 1–5: sarcosine in human metabolomes and 16S microbiomes

| Panels | Analysis | Module and scripts |
|---|---|---|
| 1a | Overview schematic | none |
| 1b–d | NSCLC plasma sarcosine, early on-treatment change and ROC | **Not in this repository.** Laboratory-shared Biocrates measurements from the previously reported cohort (Lee et al., Drug Resist. Updat. 2024); the deidentified paired values are described in the Data availability statement of the manuscript. |
| 1e–g, Sup 1a | CRC fecal sarcosine (MTBLS10232): levels, median split, ROC, anatomical-site check | `analyses/crc_fecal_metabolomics/R_scripts/`: `01`–`05` (download, parsing, merging), `13_consolidated_regeneration.R` (normalization and regeneration), `101_ROC_stackedbar_26.07.27.R` (ROC, stacked bars), `15_sarcosine_anatomical_site_consistency_26.08.23.R` |
| 1h–j, Sup 1b | ccRCC tumor versus adjacent-normal sarcosine (Hakimi et al.) | `analyses/ccrcc_tumor_metabolomics/R_scripts/`: `01_sarcosine_tumor_vs_normal.R`, `101_ROC_stackedbar_26.07.27.R` |
| 1k, Sup 2 | Matched CRC 16S (PRJNA763023) genus abundances and diversity | `analyses/crc_matched_16s_metabolomics/16S_데이터_분석/raw_download/pe_pipeline/`: `run_pe_pipeline.sh`, `02_dada2_pe.R`, `03_downstream_analysis.R` |
| 1k, Sup 3 | NSCLC stool 16S (PRJEB26531) | `analyses/nsclc_16s_microbiome/analysis/` |
| 1l, Sup 4–5 | Genus–sarcosine associations in the matched CRC subset | `…/pe_pipeline/`: `04_sarcosine_analysis.R`, `04b_sarcosine_analysis_eligens.R`, `05`–`08` (high/low comparisons), `09_genus_sarcosine_correlation_heatmap.R`. The Figure 1l / Supplementary Figure 5a effect-size panels and the Supplementary Figure 4 mixture-model panels are Prism only. |

## Figures 2–3 and Supplementary Figures 6–14: CRC shotgun metagenomes (four cohorts)

Folder: `analyses/crc_wgs_microbiome/`. Per-cohort analyses are in `PRJEB6070_CRC_AdenomatousPolyps/`, `PRJEB10878_CRC/`, `PRJEB27928_CRC/`, `PRJNA429097_CRC/`; pooled and cross-cohort analyses are in `Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/`.

| Panels | Analysis | Script family (in `Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/`) |
|---|---|---|
| 2a, 2b | Pooled cohort profiles and species differential abundance | `117_export_CRC_WGS_4cohort_pooled_inputs_26.09.01.R`, `lollipop_diff_abundance_pooled.R`, `regenerate_diff_abundance_plots.R` |
| 2d, 2e | Degradation, production and balance scores; cohort-specific effects | `integrated_analysis_pooled.R`, `113_CRC_WGS_ProdDeg_meta_heatmap_26.08.31.R`, `115_CRC_WGS_ProdDeg_pooled_heatmap_26.09.01.R` |
| 2f | Genome-wide KO differential abundance | `102_Fig1h_pooled_KO_volcano_main_26.08.22.R`, `regenerate_volcano_KO_paper.R` |
| 2g–i | soxB and creatinase: pooled and cohort-specific effects | `103_Fig1j_soxB_cross_cohort_forest_26.08.22.R`, `104_creatinase_K08688_cross_cohort_forest_26.09.01.R`, `126_remaining_5KO_pooled_forests_26.09.01.R`, `127_remaining_5KO_4cohort_forests_26.09.01.R`, `128`–`129` (source-data export and legends) |
| 2j | Candidate species selected by species–KO correlation | `123_CRC_WGS_pooled_sarcosine_species_pathway_tags_26.09.01.R`, `125`–`126_CRC_WGS_species_*_26.09.06.R` |
| 3a, 3b | Cohort membership of selected species | `plot_Fig3ab_cohort_membership_26.09.11.R` |
| 3c–g | Species abundance, KO correlations | `119_export_CRC_WGS_species_KO_correlation_inputs_26.09.01.R`, `scatter_sarcosine_species.R`, `regenerate_correlation_heatmap_paper.R` |
| Sup 6–14 | Diversity, cohort-specific differential abundance, KO sets, cross-cohort overlap | `regenerate_diversity_plots*.R`, `recompute_Fig2d_SupFig8_taxonomic_complete_case_26.08.12.R`, `cross_cohort_*consistency.R`, `100_correct_SupFig3b_PERMANOVA_26.08.20.R`, `108_SupFig3de_healthy_enriched_26.08.27.R` |

Scripts whose name contains `verify` or `validate` check the exported inputs or outputs of the analysis they follow. The input matrices are produced by `사용데이터_모음/scripts/export_analysis_used_inputs.py`.

## Figure 4 and Supplementary Figures 15–19: NSCLC ICI shotgun metagenomes

Folder: `analyses/nsclc_wgs_microbiome/` (cohort folders `NSCLC_PRJNA751792/`, `NSCLC_PRJNA1023797/`, `NSCLC_RCC_PRJEB22863/`; pooled and September panel pipelines in `pooled_analysis/`).

| Panels | Analysis | Where |
|---|---|---|
| 4a–d | Pooled diversity, species differential abundance, scores | `pooled_analysis/R_scripts/` (`pooled_*.R`, `105_Fig2f_pooled_microbiome_ecology_26.08.23.R`), `pooled_analysis/results/CRC_Method_Species_DA_26.09.09/` |
| 4e–g | Score effects and seven-KO panels | `pooled_analysis/R_scripts/CRC_method_seven_KO_NSCLC_26.09.07.R`, `redraw_CRC_seven_KO_NSCLC_26.09.07.R`, `pooled_analysis/results/NSCLC_7KO_CRC_style_panels_26.09.09/` |
| 4h | Species versus KO correlations | `sarcosine_taxa_function_scatter.R`, `sarcosine_species_scatter_grids.R` |
| 4i, 4j | Cross-disease KO and species–degradation concordance | `pooled_analysis/R_scripts/106_candidate_Fig2_genomewide_KO_cross_disease_concordance_26.08.23.R`, `102_candidate_Fig2j_cross_disease_degradation_concordance_26.08.23.R` |
| 4k, Sup 19 | Overall survival, soxB × *L. eligens* groups | `pooled_analysis/R_scripts/species_survival_km*.R`, `score_survival_km.R` |
| Sup 15–18 | Cohort-specific diversity, differential abundance, KOs, species | cohort folders, `plot_DA_lollipop_byICIgroup.R`, `replot_DA_volcano_byICIgroup.R`, `sarcosine_species_assoc_dualfilter_bar.R` |

## Figure 5 and Supplementary Figures 20–24: melanoma (TIGER PRJEB23709)

Module: `analyses/tiger_melanoma_immunotherapy_rnaseq/`. The current Figure 5 uses the 73 pretreatment tumors and the Balance score; panels g–j use 16 matched pretreatment/on-treatment pairs. The Figure 5 scripts are kept in their original layout, `results/<branch>/code/` (see `LAYOUT.md`). Tests marked "Prism" were run in GraphPad Prism.

| Panels | Analysis | Where |
|---|---|---|
| 5a | Host sarcosine production and degradation enzymes | schematic, no code |
| 5b | Balance by response, exact Mann–Whitney | Balance values are in the patient table written by `results/Fig5_PRE73_Balance_CD8_CYT_Correlations_26.09.22/code/01_prepare.py`; the P value in the legend reproduces from this table with an exact test |
| 5c | Overall four-gene expression versus Balance (data ellipses, centroids) | Prism only |
| 5d | Overall survival by median Balance | `results/Fig5_PRE73_Balance_HighLow_OS_26.09.22/code/`; hazard ratio and confidence interval: Prism log-rank estimate |
| 5e, 5f | CD8 and CYT by Balance group and versus Balance | `results/Fig5_PRE73_Balance_HighLow_CD8_CYT_26.09.22/code/`, `results/Fig5_PRE73_Balance_CD8_CYT_Correlations_26.09.22/code/`; the P values and Spearman correlations in the text reproduce from the patient table |
| 5g–5j | 16 matched pairs: Balance by time and response, centroid shifts, Balance versus CD8 and CYT | inputs (quanTIseq with IOBR, CYT): `results/Fig5_Matched16_PRE_EDT_Immune_Genes_26.09.22/code/`; the panels, including the repeated-measures ANOVA of 5g, the centroid shifts, confidence regions and paired Hotelling T² tests, are Prism only |
| Sup 20 | Responders versus non-responders: four genes, scores, ratio | `analysis/01_data_prep.R`, `analysis/02_response_analysis.R` (the gene-level P values in the figure reproduce from their output table); the score-level panels (standardization across all biopsies) are Prism only |
| Sup 21 | Joint expression distributions (data ellipses, Hotelling regions) | Prism only |
| Sup 22 | Ten measures: response odds ratios, AUCs, Cox hazard ratios | Prism only |
| Sup 23 | CD8 and CYT by response; Spearman and High–Low heatmap for ten measures | `results/Fig5_PRE73_FourGenes_CD8_CYT_26.09.22/code/`, `results/Fig5_PRE73_Ratio_Coordinates_CD8_CYT_26.09.22/code/`, `results/Fig5_PRE73_Balance_CD8_CYT_Correlations_26.09.22/code/`; the Spearman values of 23c and the exact P values of 23a,b reproduce from the patient table |
| Sup 24 | Paired CD8 and CYT changes; paired Hotelling tests | inputs: `results/Fig5_Matched16_PRE_EDT_Immune_Genes_26.09.22/code/`; the panels, including the repeated-measures ANOVA and the paired Hotelling tests, are Prism only |

## Figure 6, Supplementary Figures 25–26: shared immune programs

| Panels | Analysis | Where |
|---|---|---|
| 6a, 6b | Melanoma and ccRCC NES concordance; running-enrichment curves | `analyses/melanoma_ccrcc_shared_immune_axis/` (`01_compare_common_pathways.R`, `02_narrow_shared_immune_axis.R`). The renderer of panel 6a is not included; the panel 6b curves are Prism only. |
| 6c | NES and BH q across melanoma, ccRCC and TCGA NSCLC | `analyses/tcga_nsclc_immune_axis_validation/` |
| 6d | Post-treatment GSE207422 lineage UMAP and degradation scores | `analyses/gse207422_nsclc_scrna/02_lineage_reannotation/` (frozen annotation), `03_lee_fig3_style_FINAL/` (UMAP renders), `14_AllCell_Degradation_CD8_cDC1_26.09.08/scripts/01_prepare.R` (patient exposure and groups) |
| 6e | Whole-tumor pseudobulk, four exact programs | `15_AllCell_Degradation_Pathway_Overview_26.09.14/` |
| 6f | Conventional DCs, cross-presentation | `17_cDC_AllCell_Degradation_26.09.16/` |
| 6g | CD8⁺ T cells, TCR / TNFR2-related / IFN-γ programs | `14_AllCell_Degradation_CD8_cDC1_26.09.08/` |
| Sup 25 | Module maps and lineage composition | `03_lee_fig3_style_FINAL/`, `14_…/results/Fig6_ImmuneComposition_26.09.21/code/` |
| Sup 26 | Hallmark and Reactome overviews by compartment | `15_…/scripts/13`–`15`, `18_Epithelial_CAF_AllCell_Degradation_26.09.20/` |

## Figure 7 and Supplementary Figures 27–28: mouse and T-cell experiments

Prism only: GraphPad Prism analyses of the experimental data; there is no analysis code in this repository. The experimental data are reported in the manuscript figures and source files.
