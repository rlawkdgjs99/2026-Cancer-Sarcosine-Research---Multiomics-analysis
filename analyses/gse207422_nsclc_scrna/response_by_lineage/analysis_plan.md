# Frozen analysis plan — GSE207422 sarcosine module scores by pathologic response and lineage

Frozen before loading or examining Production, Degradation or ratio values for this response analysis: 2026-08-25 22:16 KST.

## Scientific question

Among post-treatment surgical NSCLC samples collected after neoadjuvant anti-PD-1 plus chemotherapy, do patient-level sarcosine Production, Degradation and Production/Degradation module-score summaries differ between favourable pathologic responders (MPR or pCR) and non-major pathologic responders (NMPR) within each frozen cell lineage?

This is an exploratory association analysis. It is not an analysis of ICI monotherapy, not a pretreatment predictive-biomarker test, not a measurement of sarcosine abundance or flux, and not causal.

## Locked inputs and provenance

1. Verified cell-level score table from the final Lee Figure 3-style workflow:
   `../results/tables/01_cell_paper_style_scores.csv.gz`
   SHA-256 `17724b2b1fcbad436e6ce4cee91a34bda851361aac7b14c26a3557e5d006fb61`.
2. Frozen cell-lineage and clinical metadata:
   `../../02_lineage_reannotation/results/tables/09_final_cell_lineages_FROZEN.csv`
   SHA-256 `97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1`.
3. The structure-only preflight loaded no score-value column. It verified exact cell-ID and metadata identity across the two inputs, 92,053 unique cells, 12 evaluable post-treatment surgical patients, MPR/pCR `n=4` and NMPR `n=8`, and the frozen 14-lineage set.

No score will be recomputed. Production remains shifted Seurat AddModuleScore for `GNMT + DMGDH`; Degradation remains shifted AddModuleScore for `SARDH + PIPOX`; the ratio remains the already computed cell-level shifted Production/Degradation ratio. One globally undefined denominator remains excluded only from ratio calculations.

## Cohort and groups

- Include only `Resource == "Post-treatment surgery"`.
- Include only `Pathologic.Response` in `MPR`, `pCR`, or `NMPR`.
- Define favourable pathologic response as `MPR/pCR`; pCR is retained as the strongest favourable response rather than discarded.
- Define comparison group as `NMPR`.
- Exclude the not-evaluable patient and all pre-treatment biopsy samples from the primary and sensitivity analyses.
- Expected patients: MPR/pCR `P03, P06, P11, P14`; NMPR `P02, P04, P07, P09, P10, P12, P13, P15`.

The treatment contains anti-PD-1 plus chemotherapy and uses three PD-1 antibodies and two main chemotherapy backbones. The small `n=12` dataset does not support reliable multivariable adjustment. Results will therefore be labelled unadjusted exploratory post-treatment associations.

## Frozen cell lineages

`Epithelial`, `CAF`, `B cell`, `Plasma cell`, `CD4 T cell`, `CD8 T cell`, `Cycling T cell`, `NK cell`, `Mast cell`, `Neutrophil`, `Monocyte`, `Macrophage`, `Conventional DC`, and `pDC`.

## Independent unit and cell-count eligibility

- The independent biological unit is the patient, never the cell.
- Primary eligibility: at least 10 cells in a patient×lineage compartment.
- A lineage/score is inferentially tested only when at least three eligible patients remain in each response group.
- Preflight consequence fixed before score inspection: 13 lineages are testable at the primary threshold; pDC has only two eligible MPR/pCR patients and will be shown descriptively without P or q.
- Sensitivity thresholds: at least 1, 50, and 100 cells per patient×lineage. Tests still require at least three patients per group.
- Patients are equally weighted regardless of cell count after eligibility is met.

## Primary patient-level summaries

For every eligible patient×lineage:

1. Production: arithmetic mean of the cell-level shifted Production module score.
2. Degradation: arithmetic mean of the cell-level shifted Degradation module score.
3. Production/Degradation ratio: median of finite cell-level ratios, chosen a priori because the ratio has a long upper tail and a zero-denominator boundary.

Record cell count, finite-ratio cell count, mean, median, Q1 and Q3 for every score. Do not replace absent/undefined values with zero or epsilon.

## Primary estimand, test and multiplicity

- Direction is always MPR/pCR minus NMPR.
- Report each group's patient-level median and IQR.
- Report the Hodges–Lehmann location shift: median of every pairwise `MPR/pCR - NMPR` difference.
- Obtain a 95% percentile bootstrap interval for the Hodges–Lehmann shift by resampling patients independently within groups 10,000 times, fixed master seed `260825`, with deterministic score/lineage-specific sub-seeds.
- Primary two-sided P value: exact permutation Wilcoxon rank-sum test implemented by enumerating every allocation of the observed group sizes, using midranks for ties and the absolute deviation of rank sum from its permutation expectation. Include the observed allocation in the exact tail count.
- Apply Benjamini–Hochberg correction once across all testable primary lineage×score comparisons. Expected family size is 13 lineages × 3 scores = 39.
- Interpret `q < 0.05` as FDR-significant. Report exact P and q values regardless of significance. Do not reinterpret nominal P values as confirmed evidence.

## Prespecified sensitivity analyses

1. Repeat the same patient-level summaries/tests at minimum cell thresholds 1, 50 and 100, applying BH separately within each threshold's testable family.
2. At the primary 10-cell threshold, repeat Production and Degradation with the patient×lineage cell median rather than mean.
3. At the primary threshold, repeat the balance metric as the ratio of patient×lineage mean shifted Production to mean shifted Degradation; do not add epsilon and mark a zero denominator undefined.
4. Compute leave-one-patient-out Hodges–Lehmann shifts for all primary comparisons as an influence diagnostic. These are not additional significance tests.

Sensitivity analyses may qualify robustness but cannot replace the frozen primary result.

## Figures fixed before results

1. Full three-row grouped patient-dot/box figure across all 14 lineages. Rows: Production, Degradation, ratio. Every point is one eligible patient. pDC is retained descriptively if points exist; no inferential label if the ≥3-per-group rule fails.
2. Focused publication candidate containing Epithelial, CAF, CD8 T cell and NK cell across the three scores, again showing every patient point.
3. Compact result matrix across the 14 lineages and three scores: colour shows the signed Hodges–Lehmann shift after within-score robust scaling for display only; point size shows `-log10(q)`; an outline denotes `q < 0.05`; untested cells are grey crossed symbols. The exact unscaled effects remain in tables.

Use responder blue `#2E5F8A` for MPR/pCR and non-responder amber `#C47B3B` for NMPR, Arial typography, white background, all patient points, and no cell-level violin. The manuscript/deck will not be edited unless the author subsequently selects a panel.

## Verification

1. Assert both frozen input hashes and exact cell-ID/metadata alignment.
2. Independently recompute patient summaries, exact permutation P values, Hodges–Lehmann effects, BH q values and bootstrap intervals in a separate script that does not read the primary result tables.
3. Run the complete analysis twice and compare all machine-readable CSVs and PNGs byte-for-byte. PDF metadata may contain timestamps and will be checked structurally/visually rather than by byte hash.
4. Inspect every publication PNG at full resolution and verify physical dimensions, labels, group mapping, sample counts and clipping.
5. Confirm that no raw input, manuscript, PowerPoint, Supplementary Information or Source Data workbook changed.

## Interpretation boundaries

- AddModuleScore is a relative control-gene-adjusted transcriptional score, not enzyme activity.
- The shifted score ratio is not a biochemical production/degradation flux ratio.
- Association with pathologic response does not demonstrate prediction or causation.
- Cell-state composition, post-treatment sampling, antibody/chemotherapy regimen, pathology and other patient covariates may confound results.
- Small group sizes, especially MPR/pCR `n=4`, limit precision. A visually large effect is not strong evidence without multiplicity-adjusted support and sensitivity consistency.
