# Frozen analysis plan: cell-level Hallmark GSEA for ratio Q4 versus Q1

Date frozen: 2026-08-26 (Asia/Seoul)

## Question

Within `Epithelial` cells and within `CAF` cells separately, which MSigDB Human Hallmark transcriptional programs distinguish cells in the upper versus lower quartile of the previously validated sarcosine Production/Degradation transcriptional-score ratio?

This is **not** an MPR-versus-NMPR analysis. Patient response labels will not enter cell selection, differential expression, gene ranking, gene-set testing or plotting.

## Locked biological definitions

- Production score genes: `GNMT`, `DMGDH`.
- Degradation score genes: `SARDH`, `PIPOX`.
- Ratio: the previously validated finite minimum-shifted Production/Degradation AddModuleScore ratio. It is a transcriptional-score balance, not sarcosine concentration, enzyme activity or metabolic flux.
- Lineages: the frozen reannotation labels `Epithelial` and `CAF`; analyse them separately.
- Within each lineage, calculate R type-7 Q1 and Q3 over all finite-ratio cells.
- Q1 group: ratio `<= Q1`; Q4 group: ratio `>= Q3`; exclude the middle 50%. Retain boundary ties and report them.
- Contrast orientation: `Q4 - Q1`. Positive rank/NES means Q4-enriched; negative rank/NES means Q1-enriched.

Expected locked selection from the prior verified analysis:

| Lineage | Finite-ratio cells | Q1 | Q3 | Q1 cells | Q4 cells |
|---|---:|---:|---:|---:|---:|
| Epithelial | 11,912 | 0.376174483868376 | 0.536637221540880 | 2,978 | 2,978 |
| CAF | 1,102 | 0.304701212119105 | 0.501039946395598 | 276 | 276 |

## Locked inputs and SHA-256

- Counts: `../../02_lineage_reannotation/intermediate/01_full_counts_mt20_qc.rds` — `7552ba882c306fa68380390a8172e0fc96af2c7be7cb3f5778345204e9fd61b2`.
- Frozen lineages: `../../02_lineage_reannotation/results/tables/09_final_cell_lineages_FROZEN.csv` — `97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1`.
- Validated scores: `../intermediate/01_cell_paper_style_scores.rds` — `44a1c56f947910c2930ec0f1ce1e27a3b36cf059facd1835a7ae29ea05068262`.

The workflow must abort if any input identity, dimension, cell-ID alignment, lineage count, quartile threshold or selected-cell count differs.

## Hallmark resource

- Use `msigdbr` package version `26.1.1`, database release `2026.1.Hs`, human collection `H` (`Hallmark`).
- Expected resource: exactly 50 gene sets and 7,331 gene-set membership rows before deduplication.
- Freeze the exact resource as a sorted CSV and GMT under `resources/`, record SHA-256, and use only the frozen GMT/CSV for the main and reproducibility runs.
- Use human HGNC gene symbols. Export the exact intersection and membership counts for every pathway.

## Cell-level gene ranking

For each lineage independently:

1. Select its frozen Q4 and Q1 cells only.
2. Reconstruct an RNA assay from frozen raw counts and apply Seurat `LogNormalize` with scale factor 10,000.
3. Run Seurat `FindMarkers` with `ident.1 = Q4`, `ident.2 = Q1`, normalized RNA `data`, two-sided Wilcoxon, `logfc.threshold = 0`, `min.pct = 0`, `min.diff.pct = -Inf`, `only.pos = FALSE`, no downsampling and seed `260826`.
4. Retain all returned genes with finite `avg_log2FC`; do not impose a DEG P/q or fold-change threshold before GSEA.
5. Primary rank metric: `avg_log2FC`, oriented Q4 minus Q1. Resolve only exact numeric ties with a deterministic gene-symbol-ordered perturbation no larger than `1e-12` of the ranking scale; export the unperturbed and final ranks plus a tie audit.
6. Exclude `GNMT`, `DMGDH`, `SARDH`, and `PIPOX` from the primary rank to reduce direct circularity from the score used to define Q1/Q4. Run a prespecified sensitivity with those genes retained.

This is a cell-level exploratory association. Cells are nested within 15 patients and patient contributions are unequal; cell-level P/q values are not patient-level replication.

## GSEA

- Primary engine: `fgsea::fgseaMultilevel` version `1.38.0`.
- Input: all 50 frozen Hallmark sets and the full named primary rank for each lineage.
- Parameters: `minSize = 15`, `maxSize = 500`, `eps = 0`, `scoreType = "std"`, seed `260826`.
- Export ES, NES, nominal P, pathway size, `log2err`, leading-edge genes and BH-adjusted q.
- Primary multiplicity: BH separately over all tested Hallmark sets within each lineage.
- Prespecified global sensitivity: BH over the combined Epithelial-plus-CAF test family.
- Prespecified rank sensitivity: rerun with the four score-defining genes retained.
- Independent method sensitivity: `limma::cameraPR` over the same primary rank and frozen Hallmark membership, with its directional result and BH q exported. This is a robustness check, not a replacement for fgsea.

## Display

- Create a two-facet publication plot with Epithelial and CAF separated.
- Horizontal axis: NES; positive/right is Q4-enriched and negative/left is Q1-enriched.
- Colour: direction; point size: leading-edge count; colour intensity or outline: BH q.
- Display all Hallmark pathways with primary within-lineage BH q < 0.05, capped at the ten smallest-q pathways per direction and lineage. If none pass in a direction, state that explicitly rather than filling the panel with non-significant pathways.
- Export the full 50-set results for each lineage regardless of display filtering.
- Export 600-dpi, embedded-sRGB PNG and vector PDF. Clear macOS extended attributes and hidden flags.
- Prepare a pixel-identical short ASCII PowerPoint transport PNG under `Manuscript작업/FigDesign&Manuscript/` and perform scratch-PPTX insert/reopen/render validation before reporting the raster deliverable complete.

## Verification and go/no-go

The independent verifier must recompute or assert:

1. input SHA-256 values, dimensions and cell-ID identity;
2. exact quartile thresholds, group membership, zero overlap and expected cell counts;
3. absence of response-group use in the analysis metadata/design;
4. frozen Hallmark version, 50-set count, resource hashes, unique pathway/gene membership and mapping counts;
5. Q4-minus-Q1 rank direction, finite unique named rank values and deterministic tie handling;
6. exclusion of the four score-defining genes from the primary rank and their permitted presence only in the sensitivity rank;
7. independent BH recomputation within lineage and globally;
8. NES/display direction consistency and exact displayed top-term selection;
9. complete leading-edge membership and figure-data correspondence;
10. reversed-rank sensitivity gives the expected NES sign reversal within numerical tolerance;
11. publication PNG dimensions, dpi, RGB/non-interlaced status, sRGB profile and PowerPoint transport identity;
12. a second clean full run reproduces the core resource, rank, result, plot-data and PNG files byte-identically where the file format is deterministic.

No result is guaranteed in advance. If no Hallmark set passes BH q < 0.05, report the null result. Do not select only favourable pathways, imply cancer-versus-normal biology, infer biochemical flux or make causal claims.
