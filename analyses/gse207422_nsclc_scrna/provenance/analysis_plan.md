# Lee et al. Figure 3-style sarcosine functional-score analysis

Date frozen: 2026-08-25 (Asia/Seoul), before any result from this workflow was examined.

## Objective

Reanalyse the GSE207422 single-cell RNA-seq data using the scRNA-seq logic reported for Figure 3C-F of Lee et al., *Drug Resistance Updates* 77 (2024) 101159, while changing only the biological modules to the author's prespecified sarcosine functions:

1. Production: `GNMT`, `DMGDH`
2. Degradation: `SARDH`, `PIPOX`
3. Production/degradation ratio: the paper-style ratio of non-negative Seurat module scores

The prior direct-expression analysis is retained unchanged in `analysis_sarcosine_functional_scores_26.08.25`. This is a separate method-adaptation workflow and does not overwrite it.

## Method authorities and reproducibility boundary

- Lee et al. main paper, Figure 3C-F and Methods section 4.13.
- Lee et al. Supplementary Methods paragraphs 25 and 28.
- Hu et al., *Genome Medicine* 15:14 (2023), deposited GSE207422 data and downloaded public source code, for integration and lineage provenance.

Lee et al. report sample-wise normalization, Seurat CCA integration, PCA, Louvain clustering at resolution 0.1, UMAP, marker-based annotation, `AddModuleScore`, minimum-to-zero transformation, a module-score ratio, median high/low cells, `FindMarkers`, and pathway enrichment. The exact downstream code and final Seurat object used for the published Figure 3 score panels are not public. Exact pixel- or object-level reproduction is therefore not claimed.

The frozen source-method reconstruction contains 92,053 cells meeting the reported mitochondrial-content threshold, whereas Figure 3 reports 92,031. The deposited inputs do not identify the 22 additional exclusions. All 92,053 objectively eligible cells are retained and the discrepancy remains disclosed.

The published 11-label cluster-to-lineage mapping is not recoverable exactly from the public downstream code. To avoid inventing mappings, this workflow uses the independently marker-audited frozen 14-lineage assignment from `analysis_sarcosine_degradation_REANNOTATION_26.08.24`; 90,512 cells are analysis-eligible for lineage-stratified displays. All 92,053 cells are retained for module scoring and score UMAPs. This is an intentional integrity-preserving deviation from the published 11-label display.

## Immutable inputs

- `GSE207422_NSCLC_scRNAseq_UMI_matrix.txt.gz`: `aba15960fc7ee6a2443511bce5177e4d71b964131b6e98597e7a85d0a213ba36`
- `GSE207422_NSCLC_scRNAseq_metadata.xlsx`: `d098a750c7ebc595994c929b666177f25c4da6fe3d3e38bb795269a1fb21053e`
- Full mt<=20% count matrix RDS: `7552ba882c306fa68380390a8172e0fc96af2c7be7cb3f5778345204e9fd61b2`
- Integrated clusters/UMAP table: `6bee4b491ba9f621a1df0231b91ec3dbf467e404b2ac4b0143d2d59787d5f569`
- Frozen lineage table: `97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1`
- Lee et al. main PDF: `7a7e11c364d850b54ae6f083887757fd7def0ad3e17f37556d66e8e511fe475f`
- Lee et al. Supplementary Methods DOCX: `67b301ce1712a4ce407b05d1fee87f627961125f37e273382a4486389d1fe2a6`

Raw inputs, the prior reannotation folder and the prior functional-score folder are read-only.

## Paper-style score calculation

1. Build a Seurat object from the verified 24,292-gene x 92,053-cell sparse UMI matrix.
2. Apply `NormalizeData(normalization.method="LogNormalize", scale.factor=10000)` to the RNA assay.
3. Run one `AddModuleScore` call with the two modules above, RNA `data` slot, `nbin=24`, `ctrl=100`, `k=FALSE`, `search=FALSE`, and seed 260825.
4. Let raw module scores be `P` and `D`. Following the paper's reported minimum-to-zero transformation, calculate:
   - `P0 = P - min(P)`
   - `D0 = D - min(D)`
5. Calculate the paper-style functional ratio as `P0 / D0`.

The paper does not specify what was done when the shifted denominator equals zero. No arbitrary epsilon will be invented. Cells with `D0 == 0` have an undefined ratio and are excluded only from ratio-specific displays, median classification, differential expression and enrichment. Their exact count and identifiers will be exported. Production and Degradation module-score displays retain them.

## Figure 3-style outputs

1. Frozen lineage UMAP, corresponding to Figure 3C's annotated cellular landscape.
2. Violin plots by frozen lineage for Production, Degradation and the finite Production/Degradation ratio, extending the Figure 3D ratio violin while preserving the author's three requested axes.
3. Continuous UMAPs for Production, Degradation and finite Production/Degradation ratio, corresponding to Figure 3E's continuous score UMAP.
4. Median high/low UMAP for the finite ratio, corresponding to Figure 3E's high/low UMAP.
5. High-versus-low cell differential expression using Seurat `FindMarkers` with explicit current defaults: Wilcoxon, RNA `data` slot, `logfc.threshold=0.1`, `min.pct=0.01`, `only.pos=FALSE`, no downsampling, seed 260825.
6. GO Biological Process over-representation for significant high- and low-ratio genes, corresponding in purpose to Figure 3F. Entrez mapping uses the installed `org.Hs.eg.db`; `limma::goana` provides the GO mapping and hypergeometric tests. BH correction is applied separately to the high- and low-direction GO families.

Genes entering over-representation analysis must have Seurat Bonferroni-adjusted `p_val_adj < 0.05` and absolute `avg_log2FC >= 0.25`, with the sign defining the High or Low direction. The tested-gene universe from `FindMarkers` is used as the enrichment background after Entrez mapping. Publication enrichment plots are restricted to Biological Process terms containing 10-500 represented universe genes and show at most the ten smallest BH-adjusted terms per direction; all tested terms remain in the output table.

For the median split, finite-ratio cells strictly above the median are High and those strictly below are Low. Cells exactly equal to the median are excluded from high/low DEG and their count is reported. This follows the paper's wording of above versus below rather than silently assigning ties.

## Interpretation boundary

- Cell-level violin, UMAP, high/low DEG and enrichment reproduce the paper's descriptive analysis unit and are exploratory. Individual cells are not independent patients; cell-level P values will not be presented as patient-level evidence.
- `AddModuleScore` is a relative expression score after control-gene subtraction. It is not enzyme activity, metabolite concentration, pathway flux or microbial metabolism.
- The ratio is a ratio of shifted module scores, not a biochemical production/degradation ratio.
- High/low enrichment can be driven by cell-lineage composition. Results will be cross-tabulated by lineage and patient before interpretation.
- The two-gene modules are much smaller than the gene signatures used in the source paper. Component expression and module-score/direct-expression concordance will be exported.

## Verification requirements

1. Recheck every declared input checksum.
2. Assert exact cell-ID alignment among counts, metadata, UMAP and frozen lineage tables.
3. Verify all four target genes are present exactly once.
4. Independently reproduce LogNormalize values for a deterministic cell subset.
5. Re-run `AddModuleScore` from a clean object and require maximum score difference <=1e-12.
6. Verify the non-negative shift, ratio formula, undefined denominator count, median and high/low membership.
7. Verify plot-data row counts, DEG direction, GO BH values, exact PNG dimensions and PDF existence.
8. Re-run the derived workflow and compare machine-readable results and PNG hashes.

Random seed: 260825.
