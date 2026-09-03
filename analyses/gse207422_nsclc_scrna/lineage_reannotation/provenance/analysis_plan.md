# GSE207422 single-cell reannotation plan (frozen before rerun)

Date: 2026-08-24 (Asia/Seoul)

## Objective

Reanalyse the GSE207422 processed scRNA-seq count matrix and define cell lineages using the published workflows rather than the previous custom relative marker-module winner rule. Only after the lineage annotation has passed the audits below may SARDH/PIPOX or sarcosine-degradation results be recomputed by cell lineage.

## Immutable inputs

- `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz`
  - expected SHA-256: `aba15960fc7ee6a2443511bce5177e4d71b964131b6e98597e7a85d0a213ba36`
- `GSE207422_NSCLC_scRNAseq_metadata.xlsx`
  - expected SHA-256: `d098a750c7ebc595994c929b666177f25c4da6fe3d3e38bb795269a1fb21053e`

No file outside this new analysis directory will be edited.

## Method authorities

1. Lee et al., Drug Resistance Updates 77 (2024), 101159, Supplementary Methods:
   - exclude cells with mitochondrial gene content >20%;
   - split the 15 patient datasets, normalize and select features;
   - `FindIntegrationAnchors` and `IntegrateData`;
   - PCA, `FindNeighbors`, Louvain clustering at resolution 0.1, UMAP;
   - identify cluster markers with `FindAllMarkers` and annotate from marker-expression profiles.
2. Hu et al., Genome Medicine 15 (2023), 14, the source GSE207422 study:
   - the processed matrix follows upstream QC and Scrublet doublet removal;
   - CCA integration with 3,000 variable features and dimensions 1:20;
   - all-cell validation clustering at resolution 0.6;
   - `FindAllMarkers(only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)` (modern Seurat name for the original `thresh.use=0.25`);
   - manual marker panel: `PTPRC, CD3E, LYZ, CD79A, MS4A1, IGHG1, CSF3R, FGFBP2, KIT, LILRA4, VWF, COL1A1, EPCAM, MKI67`;
   - cycling immune clusters require lineage re-clustering rather than forced assignment.
3. Source-code commit fixed at `e4c837d72b9726f5dde6a7c1b42f6cdc22b981dd` from `Junjie-Hu/NSCLC-immunotherapy`.

## Primary and validation analyses

- Primary clustering: Drug Resistance Updates workflow, Louvain resolution 0.1.
- Validation clustering: source-study workflow, resolution 0.6.
- The integrated assay is used only for dimensionality reduction and clustering. Marker expression, differential expression and target-gene summaries use the unintegrated RNA counts/data.
- Random seed is fixed at 260824. This records the requested seed but does not imply cross-version bitwise identity with the authors' older Seurat environment.

## Cell-lineage decisions

- Published marker genes and cluster marker profiles are the evidence. A label is not assigned solely because it has the largest relative score.
- The source authors' broad marker panel and their published cluster identities are the primary biological reference.
- The 2024 figure categories are: B cell, CAF, CD4 T cell, CD8 T cell, epithelial cell, mast cell, macrophage, neutrophil, NK cell, pDC and plasma cell.
- T/NK separation must show concordant evidence: T requires CD3-complex/TCR expression; NK requires NK markers and lack of broad CD3/TCR expression. A mixed T/NK cluster remains `T/NK unresolved`.
- A cluster with incompatible lineage programs, marked single-sample dominance, or suspected doublets remains `Mixed/Unresolved`.
- Cycling immune cells are re-clustered and retain `Cycling immune` if lineage cannot be resolved.
- The final mapping is frozen in a separate audited CSV before any sarcosine-gene result is run.

## Hard gates

1. Matrix dimensions, cell IDs, gene IDs and checksums must match the inputs above.
2. Mitochondrial-content QC is recomputed from the full count matrix. The expected 2024-paper count is 92,031 cells; any mismatch is reported and explained before clustering.
3. No duplicated cell or gene identifiers; no negative/non-integer counts; no zero-library cells.
4. Every final lineage must have positive-marker support and absence of a contradictory dominant program in cluster-level marker tables/plots.
5. T/NK, pDC/myeloid, B/myeloid and cycling clusters receive dedicated audits.
6. Cell-lineage response comparisons use patient x lineage pseudobulk/summary units; cells are not treated as biological replicates.
7. The completed workflow is rerun from clean derived outputs and machine-readable tables are compared.

## Interpretation boundary

The exact per-cell labels from the authors are not provided with the downloaded GEO files. Therefore the goal is a source-code-grounded, auditable reannotation, not a claim of perfect identity to the unpublished author Seurat object. Any unresolved discrepancy remains explicit and is not replaced by an invented label.
