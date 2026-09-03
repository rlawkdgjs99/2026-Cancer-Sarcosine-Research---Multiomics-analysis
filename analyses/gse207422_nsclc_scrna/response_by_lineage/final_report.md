# Patient-level sarcosine functional scores by pathologic response in GSE207422

## Verdict

This analysis does **not** demonstrate a lineage-specific association between the three prespecified sarcosine scores and pathologic response after neoadjuvant anti-PD-1 plus chemotherapy. Five of 39 exact patient-level comparisons had nominal `P < 0.05`, but **none survived BH correction**; the smallest adjusted value was `q = 0.1576`.

The result is therefore a negative/exploratory response analysis. It should not be described as evidence that Production, Degradation or their ratio predicts MPR/pCR. If shown at all, it is better suited to a complete Supplementary panel than a selectively focused main-text panel.

## Frozen design

- Frozen plan: `ANALYSIS_PLAN_FROZEN.md`
- Plan SHA-256: `a20203c72f28c344530f1a729885e37c51c32cdb00dd694157ac36c08b976bfe`
- Scores were not inspected until after the plan was written and hashed.
- Production module: `GNMT + DMGDH`
- Degradation module: `SARDH + PIPOX`
- Ratio: the already verified finite cell-level shifted Production/Degradation score ratio.
- Clinical endpoint: pathologic `MPR/pCR` versus `NMPR` after neoadjuvant anti-PD-1 plus chemotherapy. This is not an ICI-monotherapy endpoint.
- Primary cohort: post-treatment surgical samples only.
  - MPR/pCR: 4 patients (`P03`, `P06`, `P11`, `P14`)
  - NMPR: 8 patients (`P02`, `P04`, `P07`, `P09`, `P10`, `P12`, `P13`, `P15`)
  - Excluded from the primary cohort: two pretreatment NMPR biopsies (`P05`, `P08`) and one not-evaluable patient (`P01`).
- Independent unit: patient, never cell.
- Patient × lineage eligibility: at least 10 cells.
- Tested lineage requirement: at least 3 eligible patients in each response group.
- Test: exact two-sided permutation Wilcoxon rank-sum test, enumerating every group allocation with midranks.
- Effect: Hodges–Lehmann shift, MPR/pCR minus NMPR.
- Multiplicity: BH across all 39 primary tests (13 testable lineages × 3 scores).
- Primary patient summary:
  - Production: mean cell shifted module score
  - Degradation: mean cell shifted module score
  - Ratio: median finite cell-level score ratio
- Bootstrap intervals: 10,000 resamples with frozen master seed `260825`.

All 14 frozen lineages were retained. pDC was displayed but not tested because only 2 MPR/pCR and 5 NMPR patients met the 10-cell rule. The other 13 lineages were tested for all three scores.

## Primary results

No comparison had `q < 0.05`. The five nominal comparisons were:

| Lineage | Score | Eligible patients, MPR/pCR vs NMPR | HL shift, MPR/pCR − NMPR | Bootstrap 95% CI | Exact P | BH q |
|---|---|---:|---:|---:|---:|---:|
| Macrophage | Production | 4 vs 8 | -0.001696 | -0.002871 to -0.000763 | 0.00404 | 0.1576 |
| Epithelial | Degradation | 4 vs 8 | -0.019614 | -0.031895 to -0.005818 | 0.02828 | 0.3545 |
| NK cell | Production/degradation ratio | 4 vs 8 | +0.011056 | +0.003838 to +0.021235 | 0.02828 | 0.3545 |
| Mast cell | Production/degradation ratio | 4 vs 7 | -0.011330 | -0.019246 to -0.004313 | 0.04242 | 0.3545 |
| Monocyte | Production | 4 vs 8 | -0.002842 | -0.006249 to -0.000349 | 0.04848 | 0.3545 |

These are nominal patterns only. They do not support a corrected response association and must not be presented as discoveries.

For the lineages specifically raised before analysis:

- **Epithelial:** lower Degradation in MPR/pCR was nominal (`P = 0.0283`, `q = 0.3545`); the ratio was higher but not nominally significant (`P = 0.0727`, `q = 0.3545`).
- **CAF:** none of the three scores differed (`q = 0.5143–0.9004`).
- **CD8 T cell:** none differed (`q = 0.5965–0.8877`).
- **NK cell:** the ratio was nominally higher in MPR/pCR (`P = 0.0283`, `q = 0.3545`), but Production and Degradation were not significant.

## Sensitivity analyses

No sensitivity family produced a BH-FDR discovery:

| Analysis family | Tests | Minimum q | q < 0.05 |
|---|---:|---:|---:|
| Cell-count threshold ≥1 | 42 | 0.1697 | 0 |
| Primary threshold ≥10 | 39 | 0.1576 | 0 |
| Threshold ≥50 | 27 | 0.1636 | 0 |
| Threshold ≥100 | 24 | 0.1455 | 0 |
| Cell-median module-score summaries | 26 | 0.1051 | 0 |
| Ratio of patient mean module scores | 13 | 0.6303 | 0 |

Leave-one-patient-out analysis retained the direction of several nominal patterns, including macrophage Production, epithelial Degradation, epithelial and NK-cell ratios, and mast-cell ratio. This improves confidence that those directions are not caused by one patient, but it does **not** rescue statistical significance or establish predictive value.

## Figure recommendation

- Preferred full audit figure: `results/figures_publication/Fig_response_D_three_scores_all_lineages.png`
- Preferred compact overview: `results/figures_publication/Fig_response_F_effect_matrix.png`
- Focused Epithelial/CAF/CD8/NK figure: `results/figures_publication/Fig_response_E_focus_epithelial_CAF_CD8_NK.png`

The focused figure should not be used alone unless the text explicitly states that these lineages were biologically prespecified and the full 39-test result is supplied alongside it. Otherwise, it can look selectively filtered. For a Supplementary Figure, the complete all-lineage plot or complete effect matrix is the more defensible choice.

## Reproducibility and verification

- Locked score input SHA-256: `17724b2b1fcbad436e6ce4cee91a34bda851361aac7b14c26a3557e5d006fb61`
- Locked frozen-lineage input SHA-256: `97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1`
- Independent recomputation passed for:
  - every patient × lineage summary;
  - all 39 primary exact tests, HL shifts, bootstrap intervals and BH values;
  - cell-count threshold sensitivities;
  - alternative aggregation sensitivities;
  - all leave-one-patient-out HL shifts;
  - all PNG dimensions.
- The complete workflow plus independent verifier was executed twice. All 16 targeted CSV/PNG outputs were byte-identical between runs.
- Verification tables:
  - `results/tables/08_independent_verification.csv`
  - `results/tables/08_png_dimension_verification.csv`
  - `results/tables/09_double_run_reproducibility.csv`

## Interpretation boundary

These transcript-derived module scores are not measurements of sarcosine concentration, enzyme activity or metabolic flux. The cohort is small at the patient level, many lineages have fewer than 12 eligible patients, and treatment combined anti-PD-1 with chemotherapy. The analysis is associative and post-treatment; it is neither a pretreatment predictive biomarker analysis nor evidence of causality.
