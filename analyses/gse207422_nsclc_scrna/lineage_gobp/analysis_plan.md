# Frozen analysis plan: lineage-restricted ratio-quartile GO BP

Date frozen: 2026-08-26 (Asia/Seoul)

## Question

Within Epithelial cells and within CAFs separately, which transcriptional programs distinguish cells in the top versus bottom quartile of the previously validated cell-level sarcosine Production/Degradation ratio?

## Locked inputs

- Frozen sparse counts: `../../02_lineage_reannotation/intermediate/01_full_counts_mt20_qc.rds`
- Frozen lineage table: `../../02_lineage_reannotation/results/tables/09_final_cell_lineages_FROZEN.csv`
- Previously validated Lee-Figure-3-style scores: `../intermediate/01_cell_paper_style_scores.rds`

The workflow must verify SHA-256 identities before analysis. No normalization, lineage, module-score, ratio, UMAP or patient metadata definition may be changed.

## Cell selection

1. Analyse `Epithelial` and `CAF` separately.
2. Retain only cells with a finite previously calculated Production/Degradation ratio.
3. Within each lineage, calculate type-7 25th and 75th percentiles.
4. Assign ratios `<= Q1` to `Bottom 25%` and ratios `>= Q3` to `Top 25%`; exclude the middle 50% from differential expression. Retain every boundary tie and report its count.
5. Export cell and patient contributions for every lineage/quartile group.

## Differential expression and enrichment

For each lineage independently:

1. Reconstruct the RNA assay from the frozen raw counts for selected Top/Bottom cells and apply Seurat LogNormalize with scale factor 10,000. Cell-wise normalization is unchanged by subsetting cells.
2. Compare Top 25% against Bottom 25% using Seurat `FindMarkers`: Wilcoxon, RNA `data` layer, `logfc.threshold=0.1`, `min.pct=0.01`, `only.pos=FALSE`, no downsampling, seed 260826.
3. Define Top-direction genes by Bonferroni-adjusted `p_val_adj < 0.05` and `avg_log2FC >= 0.25`; define Bottom-direction genes by `p_val_adj < 0.05` and `avg_log2FC <= -0.25`.
4. Map tested symbols to Entrez IDs using the installed `org.Hs.eg.db`. The mapped genes returned by that lineage's `FindMarkers` call are the enrichment universe.
5. Run `limma::goana` hypergeometric over-representation for the Top and Bottom gene sets. Restrict publication displays to GO Biological Process terms with 10-500 represented-universe genes, at least one directional DE gene and BH q < 0.05. Apply BH correction separately to Top and Bottom within each lineage.
6. Display at most the ten smallest-q terms in each lineage/direction. Export every tested GO term.

## Display

The x axis is fold enrichment, `(DE / N) / (directional DE genes / mapped universe)`, rather than raw `DE/N`, so directions with different total DEG-set sizes are more honestly comparable. Point size is the directional DE count, and colour is `-log10(BH q)`. Export Epithelial, CAF and combined 600-dpi PNG plus vector PDF. Embed sRGB in PNG files and clear macOS extended/hidden attributes for PowerPoint compatibility.

## Interpretation boundary

This is an exploratory cell-level analysis designed to remove broad Epithelial-versus-CAF lineage composition as an explanation for the earlier all-cell GO result. Cells remain nested within 15 patients and patient contributions are unequal; cell-level P/q values are not patient-level evidence. The ratio is a ratio of minimum-shifted two-gene AddModuleScores, not biochemical flux or sarcosine concentration. Associations must not be described as causal.
