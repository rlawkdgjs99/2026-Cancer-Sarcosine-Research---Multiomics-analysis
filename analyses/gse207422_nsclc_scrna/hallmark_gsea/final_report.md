# Hallmark GSEA of sarcosine production/degradation-ratio quartiles

Date: 2026-08-26

## Question and frozen contrast

The analysis tests, separately within Epithelial cells and CAFs, which MSigDB Hallmark gene sets are associated with the upper versus lower quartile of the single-cell sarcosine Production/Degradation score ratio.

- Production score genes: `GNMT` and `DMGDH`
- Degradation score genes: `SARDH` and `PIPOX`
- Contrast: Q4 minus Q1, calculated independently within each lineage
- Epithelial: 11,912 finite-ratio cells; Q1 threshold 0.376174483868376; Q3 threshold 0.536637221540880; 2,978 Q1 and 2,978 Q4 cells; 15 patients represented in both groups
- CAF: 1,102 finite-ratio cells; Q1 threshold 0.304701212119105; Q3 threshold 0.501039946395598; 276 Q1 and 276 Q4 cells; 15 Q1 and 14 Q4 patients represented
- No cells were tied exactly at either quartile boundary
- MPR/NMPR response labels were not used

The frozen analysis-plan SHA-256 is `5ca4942bb26e7b45f7a0ab4c7f46f8718a34e705ca3a1a581012bae0d40f0cb4`.

## Method

Counts for the 6,508 selected cells were LogNormalized in Seurat with scale factor 10,000. Within each lineage, Q4 versus Q1 differential expression was ranked by Seurat `FindMarkers` Wilcoxon `avg_log2FC`, with `logfc.threshold = 0` and `min.pct = 0`. No DEG significance threshold was applied before GSEA.

The primary rank excluded the four genes used to define the score (`GNMT`, `DMGDH`, `SARDH`, `PIPOX`) to avoid circular self-enrichment. `fgseaMultilevel` 1.38.0 tested all 50 gene-symbol Hallmark sets from MSigDB Human release 2026.1.Hs using `minSize = 15`, `maxSize = 500`, `eps = 0`, `scoreType = "std"`, serial execution and seed 260826. The primary multiplicity correction was BH within each lineage; pooled BH over both lineages was retained as a sensitivity output.

## Primary results

### Epithelial cells

Eleven Hallmark sets passed within-lineage BH q < 0.05.

Q4/high-ratio enrichment:

- E2F targets: NES 1.801239, q 3.360647e-05
- G2M checkpoint: NES 1.721035, q 1.569780e-04

Q1/low-ratio enrichment:

- KRAS signaling up: NES -2.180957, q 2.541486e-09
- IL6/JAK/STAT3 signaling: NES -2.134354, q 1.334073e-06
- Allograft rejection: NES -2.101580, q 2.765551e-08
- Inflammatory response: NES -2.071485, q 2.565747e-08
- Complement: NES -1.801284, q 5.442205e-05
- Interferon-gamma response: NES -1.779594, q 6.393664e-05
- Oxidative phosphorylation: NES -1.750674, q 1.602384e-04
- IL2/STAT5 signaling: NES -1.540438, q 7.350658e-03
- Coagulation: NES -1.498255, q 2.654303e-02

Thus, in this exploratory cell-level analysis, epithelial Q4/high-ratio cells were associated mainly with cell-cycle programs, whereas epithelial Q1/low-ratio cells were associated with immune/inflammatory and related Hallmark programs.

### CAFs

No Hallmark set passed the primary within-lineage BH q < 0.05 threshold. The smallest primary q value was 0.3476477. Therefore the primary CAF Hallmark result is null.

## Sensitivity analyses

- Including `GNMT`, `DMGDH`, `SARDH` and `PIPOX` in the rank preserved the enrichment direction for all 100 lineage-pathway tests and preserved the Epithelial count of 2 Q4- and 9 Q1-enriched sets; CAF remained null.
- Reversing the rank reversed all NES signs. NES correlations were -0.999986 for CAF and -0.999985 for Epithelial; maximum absolute reverse-direction discrepancies were 0.0145 and 0.0281, respectively.
- `cameraPR` agreed with most significant Epithelial patterns. For CAF, `cameraPR` alone identified epithelial-mesenchymal transition as Q4-up (q = 0.006014), whereas the primary `fgsea` result for that set did not pass BH correction (NES 1.340765, q = 0.347648). This method-dependent CAF signal must not be presented as a robust primary finding.

## Interpretation boundary

These results describe transcriptional associations with a shifted AddModuleScore ratio, not sarcosine concentration, enzymatic flux, causal pathway activity, cancer-versus-normal differences, or patient-level treatment response. The analysis is cell-level and exploratory; cells are nested within patients, so it is not a patient-level replicated causal test.

## Verification

- Frozen input and Hallmark-resource SHA-256 values were rechecked.
- The full statistical workflow was run twice. Five core output tables were byte-identical across runs.
- An independent verifier passed 26/26 checks, including quartile membership, rank integrity, BH recomputation, pathway intersections, leading-edge membership, reversed-rank behavior, display filtering, and PNG properties.
- Publication PNG: 4,320 x 3,120 px, 7.2 x 5.2 in at 600 dpi, 8-bit RGB, non-interlaced, sRGB IEC61966-2.1.
- A short-path PowerPoint transport PNG was pixel-identical to the analytical PNG and passed scratch PPTX insertion, reopen and render validation; the pre- and post-reopen slide renders had identical SHA-256 values.

## Principal outputs

- Publication figure: `results/figures_publication/Fig_Hallmark_GSEA_Epithelial_CAF_ratio_Q4_vs_Q1.png`
- Vector figure: `results/figures_publication/Fig_Hallmark_GSEA_Epithelial_CAF_ratio_Q4_vs_Q1.pdf`
- Full primary GSEA table: `results/tables/05_fgsea_primary_all_50_sets.csv`
- Primary summary: `results/tables/07_primary_result_summary.csv`
- Independent verification: `results/tables/09_independent_verification_checks.csv`
- Frozen plan: `ANALYSIS_PLAN_FROZEN.md`
- Main analysis code: `scripts/01_run_hallmark_GSEA.R`
- Independent verifier: `scripts/02_verify_hallmark_GSEA.R`
- Presentation transport copies are not included; regenerate the PNG from the module script when needed.

No manuscript or author-owned PowerPoint file was modified by this analysis.
