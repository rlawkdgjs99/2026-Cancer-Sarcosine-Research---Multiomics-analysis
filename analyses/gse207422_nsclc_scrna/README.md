# GSE207422 NSCLC single-cell RNA-seq

This module analyzes the tumors of GSE207422 (Hu et al., Genome Medicine 2023) with a frozen, marker-audited lineage annotation and patients as the unit of replication. A patient-level degradation score (mean of the z-scored SARDH and PIPOX expression in whole-tumor pseudobulk of all 15 patients) defines degradation-High and -Low groups. Whole-tumor pseudobulk and lineage-specific pseudobulk profiles (CD8⁺ T cells, conventional dendritic cells, epithelial cells) are compared by histology-adjusted limma-voom, and the ranked genes are tested by preranked GSEA. Differential expression and enrichment treat patients, not cells, as replicates.

The deposited inputs yield 92,053 cells at the reported mitochondrial-content threshold, whereas the source paper reports 92,031. Because the exact unpublished downstream object and 11-label mapping are unavailable, this workflow retains all objectively eligible cells and uses the independently audited frozen 14-lineage map. It does not claim exact object-level reproduction of the source figure.

## Branches used for the current manuscript

Folder names are kept as in the original analysis; branches 14–20 are numbered in the order they were run. The branches read the frozen annotation and the outputs of earlier branches, so run them in numeric order (see `../../LAYOUT.md`).

| Folder | Content | Manuscript use |
|---|---|---|
| `02_lineage_reannotation/`, `lineage_reannotation/` | frozen 14-lineage annotation; `references/` holds the source article's QC, doublet and clustering scripts as reference | all cell-type analyses |
| `03_lee_fig3_style_FINAL/` | UMAP and module-score renders and panel-data exports | Figure 6d; Supplementary Figure 25a,b |
| `14_AllCell_Degradation_CD8_cDC1_26.09.08/` | patient exposure and groups (`scripts/01_prepare.R`), CD8⁺ and cDC1-like GSEA with CAMERA and patient-level score checks, immune composition (`results/Fig6_ImmuneComposition_26.09.21/code/`) | Figure 6g; Supplementary Figure 25c |
| `15_AllCell_Degradation_Pathway_Overview_26.09.14/` | whole-tumor pseudobulk GSEA, exact programs, Hallmark/Reactome overviews | Figure 6e; Supplementary Figure 26 |
| `17_cDC_AllCell_Degradation_26.09.16/` | conventional-DC GSEA with fixed sensitivity scopes | Figure 6f |
| `18_Epithelial_CAF_AllCell_Degradation_26.09.20/` | epithelial GSEA and CAF eligibility audit | Supplementary Figure 26 |
| `20_Fig6_Aligned_Immune_Programs_26.09.21/` | display of the aligned programs | Figure 6e–g |

The folders `scripts/`, `hallmark_gsea/`, `lineage_gobp/`, `response_by_lineage/` and `provenance/` come from the first (August) version of this analysis and are not used for the current Figure 6.

Most scripts take the analysis directory as a command-line argument. Large expression objects and raw matrices are not included.
