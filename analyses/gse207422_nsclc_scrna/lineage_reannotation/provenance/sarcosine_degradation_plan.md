# Sarcosine-degradation statistical plan

Frozen: 2026-08-25 (Asia/Seoul), before inspecting any lineage-specific
SARDH/PIPOX result from the reconstructed annotation.

## Primary biological question

Among post-treatment NSCLC specimens, is direct sarcosine-degradation capacity
within a defined cell lineage associated with major pathologic response?

## Analysis units and groups

- Biological replicate: patient (one deposited sample per patient).
- Expression unit: raw-count pseudobulk for each patient x final cell lineage.
- Primary contrast: post-treatment MPR/pCR (n = 4) versus post-treatment NMPR
  (n = 8), coefficient oriented as MPR minus NMPR.
- The three treatment-naive biopsy patients are displayed descriptively. They
  are not used to infer a treatment effect because treatment status and tissue
  acquisition (biopsy versus resection) are confounded and the samples are not
  paired.
- Individual cells are never treated as biological replicates.

## Targets fixed before analysis

Primary target family:

1. `SARDH`;
2. `PIPOX`;
3. direct degradation composite count: `SARDH + PIPOX` raw pseudobulk counts.

Secondary contextual family:

- the nine-gene sarcosine-metabolism signature reported by Lee et al.:
  `BHMT, BHMT2, DMGDH, ETFB, GNMT, PIPOX, SARDH, SHMT1, SHMT2`.

The nine-gene panel is not called a direct degradation module because it spans
more than the two direct sarcosine-degradation enzymes central to this
manuscript.

## Pseudobulk eligibility

- Primary inclusion threshold: at least 10 cells from a patient in a lineage,
  matching the source study's minimum for patient-level cluster signature
  summaries.
- A lineage is inferentially eligible only when at least 3 of 4 MPR patients and
  6 of 8 NMPR patients meet that threshold.
- A prespecified sensitivity analysis repeats the model using at least 50 cells
  per patient-lineage.
- A target is tested only if it passes `edgeR::filterByExpr` within that
  lineage/contrast; otherwise it is reported as `UNTESTABLE_LOW_EXPRESSION`.

## Normalization and model

For each eligible lineage:

1. sum the full raw gene-count matrix by patient;
2. construct an edgeR `DGEList` from the full pseudobulk transcriptome;
3. calculate TMM normalization factors;
4. estimate dispersion and fit a robust quasi-likelihood negative-binomial GLM
   with design `~ response_group`;
5. test MPR minus NMPR for each primary target and the prespecified secondary
   panel/score.

Raw counts, library sizes, TMM factors, cell numbers, log2 fold changes,
confidence intervals where estimable, nominal P values, and adjusted P values
will all be exported.

## Multiplicity and stability

- Benjamini–Hochberg adjustment is applied jointly across all eligible
  lineage x target tests in the primary three-target family.
- Secondary nine-gene results form a separately labelled family.
- Leave-one-patient-out estimates and the >=50-cell threshold analysis are used
  to evaluate direction stability; they do not replace the primary result.

## Main-figure decision boundary

A cell-level P value is never sufficient. A candidate main-figure panel requires
all of the following:

1. a patient-level pseudobulk effect relevant to the manuscript's T-cell or
   immune mechanism;
2. directionally compatible SARDH/PIPOX evidence rather than a result driven by
   an unrelated member of the nine-gene contextual signature;
3. no reversal under the >=50-cell sensitivity analysis;
4. no single patient determining the direction;
5. transparent distinction between post-treatment response association and
   the purely descriptive treatment-naive samples.

If these conditions fail, the analysis will not be placed in a main figure and
will not be rescued by cell-level testing.
