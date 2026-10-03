# Workspace layout expected by some scripts

Most scripts are run from the folder of their analysis module (see the module README). A subset of the later scripts instead finds its inputs by walking up from its own file location and naming sibling folders of the original analysis workspace. These scripts are included without changes to their logic (only machine-specific absolute paths were removed), so that they stay identical to the versions that produced the manuscript results.

To run them, recreate the original layout by **copying** (not symlinking, because the scripts resolve real paths) the repository folders:

| Repository folder | Folder name expected by the scripts (relative to a workspace root) |
|---|---|
| `analyses/crc_wgs_microbiome/` | `공공_Metabolomics&Metagenomics_분석모음/HGMT_CRC_WGS-Healthy_vs_Cancer/` |
| `analyses/nsclc_wgs_microbiome/` | `공공_Metabolomics&Metagenomics_분석모음/HGMT_NSCLC_ICI_RvsNR_WGS/` |
| `analyses/_shared/` | `공공_Metabolomics&Metagenomics_분석모음/_shared/` |
| `analyses/gse207422_nsclc_scrna/<NN_branch>/` | `2024_Drug_Res_Updates_NSCLC/RNA-seq공공데이터_GSE207422/analysis_sarcosine_FINAL_26.08.25/<NN_branch>/` (branch folders keep their original names) |
| `analyses/tiger_melanoma_immunotherapy_rnaseq/` | `공공_BulkRNAseq_ICI반응성관련_분석모음/Sarcosine-TIGER상_ICI반응성_RNAseq분석/Melanoma-PRJEB23709/` (the Figure 5 scripts sit in `results/<branch>/code/`) |

Several Figure 5 scripts of the melanoma module locate the workspace root by looking upward for a file named `PROJECT_HANDOFF.md`; create an empty file with that name in the workspace root (an ancestor of the copied module folder). The matched-pair scripts also read the TIGER FPKM matrix from `사용데이터_모음/Source_Input/` next to the module folder and the IOBR 2.2.3 reference data from the module's `reference_data/`. The author-prepared workbooks `73Pre.xlsx` and `16Matched.xlsx` are not distributed (see `DATASETS.md`).

Translation of the folder names:`공공_Metabolomics&Metagenomics_분석모음` = "public metabolomics and metagenomics analyses"; `사용데이터_모음` = "collection of the data actually used" (matrices written by `export_analysis_used_inputs.py` from the HGMT downloads); `공공_BulkRNAseq_ICI반응성관련_분석모음` = "public bulk RNA-seq analyses related to ICI response"; `Sarcosine-TIGER상_ICI반응성_RNAseq분석` = "sarcosine, TIGER ICI-response RNA-seq analysis".

Input data are not distributed. The exported input matrices, HGMT downloads, GSE207422 matrices and TCGA/TIGER downloads must be obtained as described in `DATASETS.md`.

## Scripts that depend on this layout (86)

**crc_wgs_microbiome** (25)

- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/results_integrated/CRC_WGS_6KO_input_export_26.09.07/verify_6KO_inputs.py`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/102_Fig1h_pooled_KO_volcano_main_26.08.22.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/105_Fig1j_boxplot_BioSample_26.08.25.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/106_Fig1j_display_variants_26.08.25.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/107_Fig1j_boxplot_runlevel_26.08.25.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/115_CRC_WGS_ProdDeg_pooled_heatmap_26.09.01.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/116_verify_CRC_WGS_ProdDeg_pooled_heatmap_26.09.01.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/117_export_CRC_WGS_4cohort_pooled_inputs_26.09.01.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/118_verify_CRC_WGS_4cohort_pooled_inputs_26.09.01.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/119_export_CRC_WGS_species_KO_correlation_inputs_26.09.01.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/120_verify_CRC_WGS_species_KO_correlation_inputs_26.09.01.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/121_CRC_WGS_sarcosine_7KO_lollipop_26.09.01.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/122_verify_CRC_WGS_sarcosine_7KO_lollipop_26.09.01.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/123_CRC_WGS_pooled_sarcosine_species_pathway_tags_26.09.01.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/124_verify_CRC_WGS_pooled_sarcosine_species_pathway_tags_26.09.01.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/125_CRC_WGS_species_readable_26.09.06.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/126_CRC_WGS_species_bar_style_26.09.06.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/126_remaining_5KO_pooled_forests_26.09.01.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/127_remaining_5KO_4cohort_forests_26.09.01.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/128_export_7KO_4cohort_forest_source_data_26.09.01.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/129_7KO_4cohort_forests_explicit_legend_26.09.01.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/129_verify_7KO_4cohort_forest_source_data_26.09.01.R`
- `analyses/crc_wgs_microbiome/Healthy_vs_Cancer_4_CRC_cohorts_integrated/scripts/plot_Fig3ab_cohort_membership_26.09.11.R`
- `analyses/crc_wgs_microbiome/사용데이터_모음/scripts/export_analysis_used_inputs.py`
- `analyses/crc_wgs_microbiome/사용데이터_모음/scripts/export_sarcosine_functional_scores.py`

**gse207422_nsclc_scrna** (24)

- `analyses/gse207422_nsclc_scrna/03_lee_fig3_style_FINAL/results/figures_publication/Fig6_Post12_Group_Module_UMAP_26.09.29/code/export_four_csv.py`
- `analyses/gse207422_nsclc_scrna/03_lee_fig3_style_FINAL/results/figures_publication/UMAP_Contrast_26.09.20/qa/verify_outputs.py`
- `analyses/gse207422_nsclc_scrna/14_AllCell_Degradation_CD8_cDC1_26.09.08/results/Fig6_Degradation_CellFraction_Correlations_26.09.29/code/analyze.py`
- `analyses/gse207422_nsclc_scrna/14_AllCell_Degradation_CD8_cDC1_26.09.08/results/Fig6_ImmuneComposition_26.09.21/code/02_verify.py`
- `analyses/gse207422_nsclc_scrna/14_AllCell_Degradation_CD8_cDC1_26.09.08/results/Fig6_ImmuneComposition_26.09.21/code/04_finalize.py`
- `analyses/gse207422_nsclc_scrna/14_AllCell_Degradation_CD8_cDC1_26.09.08/results/Fig6_ImmuneComposition_26.09.21/code/07_export_statistics_csv.py`
- `analyses/gse207422_nsclc_scrna/14_AllCell_Degradation_CD8_cDC1_26.09.08/results/Fig8_CellType_Pathways_26.09.14/scripts/04_artifact_qa.py`
- `analyses/gse207422_nsclc_scrna/14_AllCell_Degradation_CD8_cDC1_26.09.08/results/Fig8_CellType_Pathways_26.09.14/scripts/06_finalize.py`
- `analyses/gse207422_nsclc_scrna/14_AllCell_Degradation_CD8_cDC1_26.09.08/scripts/05_artifact_qa.py`
- `analyses/gse207422_nsclc_scrna/14_AllCell_Degradation_CD8_cDC1_26.09.08/scripts/06_finalize.py`
- `analyses/gse207422_nsclc_scrna/14_AllCell_Degradation_CD8_cDC1_26.09.08/scripts/08_finalize_Fig8a_Cell_Context_26.09.17.py`
- `analyses/gse207422_nsclc_scrna/14_AllCell_Degradation_CD8_cDC1_26.09.08/scripts/10_verify_CD8_cDC_Proportions_26.09.17.py`
- `analyses/gse207422_nsclc_scrna/15_AllCell_Degradation_Pathway_Overview_26.09.14/Fig6e_Assigned_26.09.20/draw_fig6e_split.R`
- `analyses/gse207422_nsclc_scrna/15_AllCell_Degradation_Pathway_Overview_26.09.14/scripts/06_scRNA_GSEA_heatmap.R`
- `analyses/gse207422_nsclc_scrna/15_AllCell_Degradation_Pathway_Overview_26.09.14/scripts/07_verify_scRNA_heatmap.py`
- `analyses/gse207422_nsclc_scrna/15_AllCell_Degradation_Pathway_Overview_26.09.14/scripts/08_scRNA_two_direction_heatmap.R`
- `analyses/gse207422_nsclc_scrna/15_AllCell_Degradation_Pathway_Overview_26.09.14/scripts/09_verify_scRNA_two_directions.py`
- `analyses/gse207422_nsclc_scrna/15_AllCell_Degradation_Pathway_Overview_26.09.14/scripts/11_plot_exact_GSEA_curves.R`
- `analyses/gse207422_nsclc_scrna/15_AllCell_Degradation_Pathway_Overview_26.09.14/scripts/15_verify_Hallmark_Reactome.py`
- `analyses/gse207422_nsclc_scrna/18_Epithelial_CAF_AllCell_Degradation_26.09.20/scripts/05_verify_export.py`
- `analyses/gse207422_nsclc_scrna/20_Fig6_Aligned_Immune_Programs_26.09.21/scripts/01_plot.R`
- `analyses/gse207422_nsclc_scrna/20_Fig6_Aligned_Immune_Programs_26.09.21/scripts/02_verify_export.py`
- `analyses/gse207422_nsclc_scrna/20_Fig6_Aligned_Immune_Programs_26.09.21/scripts/03_plot_large.R`
- `analyses/gse207422_nsclc_scrna/20_Fig6_Aligned_Immune_Programs_26.09.21/scripts/04_verify_large.py`

**nsclc_wgs_microbiome** (23)

- `analyses/nsclc_wgs_microbiome/pooled_analysis/R_scripts/101_candidate_Fig2g_pathway_score_meta_26.08.22.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/R_scripts/102_candidate_Fig2j_cross_disease_degradation_concordance_26.08.23.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/R_scripts/103_validate_CRC_KO_cache_against_raw_26.08.23.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/R_scripts/104_Fig2i_individual_sarcosine_gene_meta_26.08.23.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/R_scripts/105_Fig2f_pooled_microbiome_ecology_26.08.23.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/R_scripts/106_candidate_Fig2_genomewide_KO_cross_disease_concordance_26.08.23.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/R_scripts/CRC_method_seven_KO_NSCLC_26.09.07.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/R_scripts/cross_disease_deg_core_overlap_26.08.19.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/R_scripts/redraw_CRC7_direction_labels_26.09.07.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/R_scripts/redraw_CRC_seven_KO_NSCLC_26.09.07.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/R_scripts/validate_CRC_method_seven_KO_NSCLC_26.09.07.py`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/results/CRC_Method_7KO_NSCLC_26.09.07/Fig3_eijk_Update_26.09.07/scripts/01_analysis.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/results/CRC_Method_7KO_NSCLC_26.09.07/Fig3_eijk_Update_26.09.07/scripts/02_plot.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/results/CRC_Method_7KO_NSCLC_26.09.07/Fig3_eijk_Update_26.09.07/scripts/03_verify_numerics.py`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/results/CRC_Method_Species_DA_26.09.09/per_cohort/scripts/analysis.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/results/CRC_Method_Species_DA_26.09.09/per_cohort/scripts/plot.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/results/CRC_Method_Species_DA_26.09.09/per_cohort/scripts/verify.py`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/results/CRC_Method_Species_DA_26.09.09/scripts/analysis.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/results/CRC_Method_Species_DA_26.09.09/scripts/plot.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/results/CRC_Method_Species_DA_26.09.09/scripts/verify.py`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/results/NSCLC_7KO_CRC_style_panels_26.09.09/analysis.R`
- `analyses/nsclc_wgs_microbiome/pooled_analysis/results/NSCLC_7KO_CRC_style_panels_26.09.09/verify_cohort_membership.py`
- `analyses/nsclc_wgs_microbiome/사용데이터_모음/scripts/export_analysis_used_inputs.py`

**tiger_melanoma_immunotherapy_rnaseq** (14)

- `analyses/tiger_melanoma_immunotherapy_rnaseq/results/Fig5_Matched16_PRE_EDT_Immune_Genes_26.09.22/code/01_prepare.py`
- `analyses/tiger_melanoma_immunotherapy_rnaseq/results/Fig5_Matched16_PRE_EDT_Immune_Genes_26.09.22/code/02_quantiseq.R`
- `analyses/tiger_melanoma_immunotherapy_rnaseq/results/Fig5_PRE73_Balance_CD8_CYT_Correlations_26.09.22/code/01_prepare.py`
- `analyses/tiger_melanoma_immunotherapy_rnaseq/results/Fig5_PRE73_Balance_CD8_CYT_Correlations_26.09.22/code/02_analyse_plot.R`
- `analyses/tiger_melanoma_immunotherapy_rnaseq/results/Fig5_PRE73_Balance_CD8_CYT_Correlations_26.09.22/code/03_verify.py`
- `analyses/tiger_melanoma_immunotherapy_rnaseq/results/Fig5_PRE73_Balance_HighLow_OS_26.09.22/code/01_prepare.py`
- `analyses/tiger_melanoma_immunotherapy_rnaseq/results/Fig5_PRE73_Balance_HighLow_OS_26.09.22/code/02_analyse_plot.R`
- `analyses/tiger_melanoma_immunotherapy_rnaseq/results/Fig5_PRE73_Balance_HighLow_OS_26.09.22/code/03_verify.py`
- `analyses/tiger_melanoma_immunotherapy_rnaseq/results/Fig5_PRE73_FourGenes_CD8_CYT_26.09.22/code/01_prepare.py`
- `analyses/tiger_melanoma_immunotherapy_rnaseq/results/Fig5_PRE73_FourGenes_CD8_CYT_26.09.22/code/02_analyse_plot.R`
- `analyses/tiger_melanoma_immunotherapy_rnaseq/results/Fig5_PRE73_FourGenes_CD8_CYT_26.09.22/code/03_verify_export.py`
- `analyses/tiger_melanoma_immunotherapy_rnaseq/results/Fig5_PRE73_Ratio_Coordinates_CD8_CYT_26.09.22/code/01_prepare.py`
- `analyses/tiger_melanoma_immunotherapy_rnaseq/results/Fig5_PRE73_Ratio_Coordinates_CD8_CYT_26.09.22/code/02_analyse_plot.R`
- `analyses/tiger_melanoma_immunotherapy_rnaseq/results/Fig5_PRE73_Ratio_Coordinates_CD8_CYT_26.09.22/code/03_verify_export.py`

