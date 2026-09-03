# GSE207422 sarcosine-degradation reanalysis: final report

Date: 2026-08-25 (Asia/Seoul)

## Bottom line

The patient-level analysis does **not** support a response-associated main-
figure panel. Across reconstructed cell lineages, neither `SARDH`, `PIPOX`, nor
the direct-degradation composite count (`SARDH + PIPOX`) was associated with
MPR/pCR after the prespecified multiple-testing correction. The minimum primary
family FDR was 0.798541 at the >=10-cell threshold.

A separate descriptive result is reproducible: `SARDH` is preferentially
detected in T-cell lineages, especially cycling T cells, whereas `PIPOX` is
preferentially detected in macrophages, conventional dendritic cells, and
monocytes. This is cell-lineage localization, not evidence of treatment
response or sarcosine-mediated T-cell function.

## Data and source-method audit

- Deposited count matrix: 24,292 genes x 92,330 cells before the mitochondrial
  filter.
- Cells with mitochondrial fraction <=20%: 92,053.
- The published 2024 report states 92,031 cells. The 22-cell excess in the
  deposited-data reconstruction remains unexplained. No undocumented deletion
  was introduced to force agreement.
- Input UMI SHA-256:
  `aba15960fc7ee6a2443511bce5177e4d71b964131b6e98597e7a85d0a213ba36`.
- Input metadata SHA-256:
  `d098a750c7ebc595994c929b666177f25c4da6fe3d3e38bb795269a1fb21053e`.
- The downloaded original-team code contains 27 files at commit
  `e4c837d72b9726f5dde6a7c1b42f6cdc22b981dd`.
- The exact all-cell workflow and downstream T/NK annotation script are
  present. The code creating the authors' intermediate `Tcell.rds` and
  `Tmem.rds` objects is absent. Therefore, this workflow reconstructs the
  published method but does not claim identity to the unavailable objects.

The reconstructed workflow followed the reported source parameters: sample-
wise normalization, 3,000 variable features, CCA integration using dimensions
1:20, PCA/neighbors/UMAP using the first 15 PCs, and resolution 0.4 for the
T/NK and memory-T re-clustering. Final labels were frozen before the sarcosine
association analysis.

## Final cell-lineage reconstruction

The frozen annotation contains 92,053 cells and 16 labels:

- 90,512 cells in 14 analysis-eligible lineages;
- 1,332 unresolved NK/gamma-delta T cells excluded from lineage inference;
- 209 residual T-myeloid doublets excluded according to the source study's
  incompatible-lineage-marker rule.

Major eligible lineage counts were CD8 T 17,614; CD4 T 14,858; macrophage
12,473; epithelial 11,913; neutrophil 10,009; B 7,292; monocyte 4,984; plasma
3,421; NK 2,643; cycling T 2,074; CAF 1,102; conventional DC 1,007; mast 845;
and pDC 277.

Frozen lineage-table SHA-256:
`97f1438f1f80082b5f1c4977288adc005bc3616a3d8a49c04f34cf1a8301e1a1`.

## Statistical design

- Biological replicate: patient, never an individual cell.
- Primary contrast: four post-treatment MPR/pCR patients minus eight
  post-treatment NMPR patients.
- Three treatment-naive biopsy patients were retained for descriptive
  summaries only because treatment status and sampling procedure are
  confounded and samples are not paired.
- Expression unit: full-transcriptome patient x lineage raw-count pseudobulk.
- Primary inclusion threshold: >=10 cells per patient-lineage; prespecified
  sensitivity threshold: >=50 cells.
- A lineage required at least 3 MPR and 6 NMPR patients.
- Model: edgeR robust quasi-likelihood negative-binomial GLM with TMM
  normalization and coefficient MPR minus NMPR.
- Primary targets fixed before result inspection: `SARDH`, `PIPOX`, and raw
  pseudobulk `SARDH + PIPOX` composite count.
- BH correction was applied jointly over tested lineage x target combinations
  in the primary family.
- Low-expression targets failing `filterByExpr` were reported as untestable,
  not treated as zero effects.

The statistical plan was frozen before lineage-specific target results were
inspected. Plan SHA-256:
`7443fdae621896ecc5d05b89d7be603e47386e4fc0af5969ca69497bb146fab6`.

## Primary response-association results

Eleven lineage-target combinations passed the expression gate at the primary
>=10-cell threshold. None survived BH correction.

| Lineage | Target | log2FC MPR-NMPR | 95% CI | P | BH q |
|---|---:|---:|---:|---:|---:|
| Epithelial | SARDH+PIPOX | -1.515 | -3.186 to 0.156 | 0.07259 | 0.79854 |
| Macrophage | SARDH | -0.704 | -2.738 to 1.330 | 0.47039 | 0.98234 |
| CD4 T cell | SARDH | 0.234 | -0.598 to 1.066 | 0.55622 | 0.98234 |
| CD8 T cell | SARDH | -0.292 | -1.334 to 0.750 | 0.55750 | 0.98234 |
| CD4 T cell | SARDH+PIPOX | 0.218 | -0.588 to 1.023 | 0.57200 | 0.98234 |
| CD8 T cell | SARDH+PIPOX | -0.256 | -1.260 to 0.748 | 0.59272 | 0.98234 |
| Monocyte | SARDH+PIPOX | -0.212 | -1.229 to 0.805 | 0.66148 | 0.98234 |
| Macrophage | SARDH+PIPOX | -0.123 | -0.831 to 0.586 | 0.71443 | 0.98234 |
| Cycling T cell | SARDH | -0.023 | -0.944 to 0.898 | 0.95879 | 0.98678 |
| Cycling T cell | SARDH+PIPOX | -0.023 | -0.944 to 0.898 | 0.95879 | 0.98678 |
| Macrophage | PIPOX | 0.005 | -0.692 to 0.703 | 0.98678 | 0.98678 |

At the >=50-cell sensitivity threshold, the epithelial composite was nominally
negative (log2FC -1.628, 95% CI -3.163 to -0.093, P=0.03915) but remained
non-significant after primary-family correction (q=0.35234) and was based on
4 MPR versus 7 NMPR patients. It is not a validated response signal.

Leave-one-patient-out directions were stable for some null estimates (CD8
`SARDH` was negative in 12/12 folds; CD4 `SARDH` positive in 11/12), but all
full-model confidence intervals crossed zero and all primary q values were
null. Directional stability alone is not evidence of association.

## Secondary published signature

The Lee et al. signature declares nine genes, but `BHMT` and `BHMT2` are absent
from the deposited 24,292-gene matrix. They were explicitly reported as
`NOT_AVAILABLE_IN_DEPOSITED_MATRIX` and were not replaced with zeros. Seven of
nine genes were analyzed as a separately labelled contextual family. Nominal
ETFB/SHMT2 findings did not survive BH correction; the smallest secondary q was
0.19288.

## Descriptive cell-lineage localization

For each patient x lineage with >=10 cells, the fraction of cells with at least
one target UMI and the mean log1p(CP10k) across all cells were calculated. The
reported values below are medians across patients. Within-patient ranks were
also computed to check that localization was not driven by a single deeply
sequenced specimen.

### SARDH

- Cycling T: 10.000% median detection (IQR 5.516-13.566%); detection top-three
  lineage in 13/15 patients (86.7%).
- CD4 T: 3.504% (IQR 3.011-7.250%); top-three in 13/15 (86.7%).
- CD8 T: 2.784% (IQR 1.692-5.983%); top-three in 9/15 (60.0%).

### PIPOX

- Macrophage: 5.497% median detection (IQR 4.699-9.059%); top-three in 15/15
  patients (100%).
- Conventional DC: 5.778% (IQR 4.720-6.667%); top-three in 12/14 (85.7%).
- Monocyte: 1.737% (IQR 1.009-2.972%); top-three in 9/15 (60.0%).

This separation is biologically descriptive but does not show that sarcosine
caused the lineage states, nor that the localization predicts MPR/pCR. The
diagnostic dot plot should therefore not be relabelled as an inferential
response panel.

## Figure decision

**Response-association panel: NO-GO for a main figure.** The prespecified
patient-level targets are null after correction, and T-cell directions do not
form a coherent MPR-associated degradation signal.

**Localization dot plot: at most a descriptive/contextual or supplementary
panel.** It can honestly show that SARDH and PIPOX occupy different NSCLC cell
compartments. It does not rescue the response result and, by itself, is weaker
than a mechanistic Figure 4 panel.

## Reproducibility and verification

- Original input SHA-256 values remained unchanged.
- Count and frozen-lineage cell IDs aligned exactly; unmatched cells: 0.
- All 14 full-transcriptome lineage pseudobulks were independently checked.
- Independent target-count aggregation had maximum difference 0.
- BH values, confidence-interval containment, leave-one-out summaries, and
  explicit missing-gene handling passed 35 checks.
- The localization table passed six independent sparse-matrix aggregation
  checks. Maximum differences were 4.44e-14 percentage points for detection
  and 1.55e-15 for mean log1p(CP10k).
- The localization script was run twice; the summary-table and PNG SHA-256
  values were byte-identical between runs.

Key output SHA-256 values:

- Primary model results:
  `364778cbef6ab97dbdf4774761d57fb0ec38c318fb31e6f965df0b81c04e200c`
- Primary verification checks:
  `4739d0dd6e9adcb07de7c67b147b92c997320680cc95997197021357981b3b1e`
- Localization summary:
  `0ddde0d45aa79a006a7058abe78e1dc12e91c05b41454eab1fc822813c38a7d6`
- Localization PNG:
  `8477e50269d104b8a3b4d9e18adbd51ea8d6ca09b01f1b53a3e7100ff9251ea8`

No manuscript, PowerPoint deck, handoff file, or other analysis directory was
modified by this workflow.
