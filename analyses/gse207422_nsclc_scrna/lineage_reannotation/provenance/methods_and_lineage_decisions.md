# Source-method and lineage decisions before sarcosine analysis

Date: 2026-08-25 (Asia/Seoul)

This addendum does not change the frozen analysis plan. It records information
that was verified after the original research team's full public code archive
was downloaded and audited.

## Verified source-method parameters

- Hu et al., *Genome Medicine* 15:14 (2023), PDF pp. 5–6:
  - sample-wise `NormalizeData`;
  - CCA integration using 3,000 variable features for T/myeloid analyses;
  - integration anchors from dimensions 1:20;
  - PCA followed by `FindNeighbors` using the first 15 PCs;
  - `FindClusters` resolution 0.4 for T and myeloid cells;
  - UMAP using the same 15 PCs;
  - `FindAllMarkers(min.pct = 0.25, thresh.use = 0.25)` for cluster markers;
  - manual removal of clusters co-expressing two or more incompatible major
    lineage markers as residual doublets.
- Hu et al., PDF p. 11:
  - T/NK cells were re-clustered into 14 annotated subtypes: two NK, five CD8,
    four conventional CD4, two Treg, and one proliferating subtype.
- Lee et al., *Drug Resistance Updates* 77:101159 (2024), Supplementary
  Methods, paragraphs 25 and 28:
  - the sarcosine signature was `BHMT, BHMT2, DMGDH, ETFB, GNMT, PIPOX,
    SARDH, SHMT1, SHMT2`;
  - their scRNA analysis used the deposited GSE207422 data, mitochondrial
    content <=20%, CCA integration, and cell-level `AddModuleScore`.

## What the public code does and does not contain

- Commit `e4c837d72b9726f5dde6a7c1b42f6cdc22b981dd` contains the exact all-cell
  workflow and the downstream T/NK annotation script.
- `05_T_analysis_new.R` starts by reading `Tcell.rds`; the code that created
  that object and `Tmem.rds` is not present. The public code therefore does not
  expose the exact cell-selection vector, manual doublet removals, or the
  separate memory-T split.
- The source paper nevertheless reports the integration, PC, and resolution
  parameters above. The present workflow reconstructs those reported
  parameters but does **not** claim identity to the unavailable author objects.

## T/NK input decision for the reconstruction

The 23-cluster all-cell validation solution (resolution 0.6) is used only to
select the broad T/NK compartment before source-method re-clustering. The
selected clusters are 0, 1, 2, 3, 10, and 14 (38,730 cells total):

- 0: CD3E/CD8A/GZMK/NKG7-positive cytotoxic lymphocytes;
- 1: mixed cytotoxic CD8/NK compartment retained for re-clustering rather than
  forced to either lineage;
- 2: CD3E/CD4/IL7R-positive T cells;
- 3: activated CD4/Treg-like T cells;
- 10: CD3E/CD8A/GZMB/CXCL13-positive cytotoxic/exhausted T cells;
- 14: MKI67/TOP2A-positive cycling T cells.

The pDC cluster (22), despite GZMB expression, is excluded because it expresses
the pDC program (LILRA4/GZMB/IRF8/GPR183/TCF4) rather than a T/NK program.

No sarcosine result will be calculated until the reconstructed T/NK clusters
have been audited against the authors' complete T/NK marker panel and the final
annotation mapping has been frozen.

## Clinical-group reconstruction

- Pre-treatment biopsy samples: `TN` (n = 3 patients).
- Post-treatment samples with MPR or pCR: `MPR` (n = 4 patients, pCR included).
- Post-treatment samples with NMPR: `NMPR` (n = 8 patients).

This reproduces the group counts reported by Hu et al. The later inferential
unit will be patient x lineage pseudobulk; individual cells will not be treated
as biological replicates.
