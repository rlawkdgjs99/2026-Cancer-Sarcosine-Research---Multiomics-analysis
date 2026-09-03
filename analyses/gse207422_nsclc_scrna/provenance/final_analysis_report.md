# Lee et al. Figure 3-style sarcosine functional-score analysis

Date completed: 2026-08-25 (Asia/Seoul)

## Bottom line

The requested method adaptation is complete. The only biological substitutions were:

- Production: `GNMT + DMGDH`
- Degradation: `SARDH + PIPOX`

The downstream workflow follows the logic reported for Lee et al. Figure 3C-F: Seurat `LogNormalize`, `AddModuleScore`, minimum-to-zero score transformation, Production/Degradation ratio, finite-score median High/Low split, cell-level `FindMarkers`, and pathway enrichment.

This is not an exact reproduction of the published figure. The authors did not release the exact downstream Figure 3 score script or final Seurat object, their figure reports 92,031 cells whereas the deposited inputs objectively yield 92,053 cells at the reported mitochondrial threshold, and the exact published 11-label cluster mapping cannot be recovered from the public code. This workflow therefore retains all 92,053 objectively eligible cells and uses the independent marker-audited frozen 14-lineage map rather than inventing labels.

## Exact result

- Cells scored: 92,053.
- Raw Production score minimum: -0.0458731238487056.
- Raw Degradation score minimum: -0.148574315360367.
- Both minima were shifted to exactly zero as reported in the source method.
- Cells with shifted Degradation equal to zero: 1. No epsilon was added; this cell was excluded only from ratio-specific analyses.
- Finite ratios: 92,052.
- Finite-score median: 0.343967676527376.
- High-ratio cells: 46,026.
- Low-ratio cells: 46,026.
- Exact median ties: 0.

The ratio is strongly structured by cell lineage. Median ratios were highest in epithelial cells (0.4505) and CAFs (0.3937), whereas the lowest medians were in mast cells (0.3199), pDCs (0.3283), neutrophils (0.3315), CD4 T cells (0.3324) and CD8 T cells (0.3334). Consistently, 84.1% of epithelial cells and 64.8% of CAFs were above the global median, compared with 42.3% of CD4 T cells, 43.5% of CD8 T cells, 32.8% of mast cells and 25.9% of neutrophils.

Cell-level High-versus-Low `FindMarkers` tested 11,515 genes. Under the frozen selection rule (Bonferroni-adjusted `p < 0.05` and absolute average log2 fold change >=0.25), 5,625 genes were High-direction and 1,185 were Low-direction. GO Biological Process over-representation was dominated by epithelial/ciliary and morphogenesis terms in High-ratio cells and immune/T-cell activation and cytokine-response terms in Low-ratio cells.

This enrichment contrast is primarily a cell-lineage-composition result, not evidence that the ratio cell-intrinsically suppresses immune activation. The High group is enriched for epithelial/CAF cells, while multiple immune lineages are more frequent in the Low group. Cell-level P values also do not provide patient-level replication.

## Interpretation boundary

The analysis supports a descriptive statement:

> Sarcosine Production/Degradation transcript-module balance differs across tumour cell lineages, with relatively higher ratios in epithelial and stromal compartments and lower ratios in several immune compartments.

It does not establish:

- measured sarcosine production or degradation;
- enzyme activity or metabolic flux;
- a causal effect of the ratio on immune function;
- immunotherapy response;
- patient-level association.

`AddModuleScore` is control-gene-subtracted relative expression. The ratio is a ratio of two minimum-shifted module scores, not a biochemical Production/Degradation flux ratio. Production is particularly sparse: 99% of cells have the same shifted Production value (0.0458731), so continuous Production gradients should not be over-interpreted.

## Publication outputs

All PNGs are 600 dpi and have matching vector PDFs in `results/figures_publication/`.

- `Fig3_style_C_lineage_UMAP`: annotated cell landscape.
- `Fig3_style_D_ratio_violin`: preferred Figure 3D-style ratio-by-lineage panel; upper 1% is winsorized for display only and raw values remain unchanged in the tables.
- `Fig3_style_D_three_axis_violins`: three-axis diagnostic; less visually efficient because Production is extremely sparse.
- `Fig3_style_E_continuous_UMAPs`: Production, Degradation and ratio continuous UMAPs.
- `Fig3_style_E_ratio_high_low_UMAP`: finite-ratio median High/Low UMAP.
- `Fig3_style_F_high_low_GO_BP`: exploratory High/Low GO Biological Process comparison.
- `Fig_scRNA_Lee2024_Fig3_style_compact`: compact Figure 3C-E analogue.

If one concise descriptive panel is needed, prefer `Fig3_style_D_ratio_violin` together with either the lineage UMAP or the High/Low UMAP. Do not use the GO panel alone as mechanistic evidence.

## Verification

- Input SHA-256 values and exact cell-ID sets were rechecked.
- Count and lineage files contained 92,053 unique cell IDs each; there were zero unmatched IDs.
- Independent raw-count LogNormalize verification had a maximum absolute difference of 5.55e-17.
- A clean Seurat object reproduced both AddModuleScore vectors with maximum absolute difference exactly 0.
- Ratio transformation, denominator-zero handling, median and all High/Low memberships were independently reproduced.
- The complete cell-level Wilcoxon run was performed twice. The 11,515-row DEG table was byte-identical both times (SHA-256 `af5b154e2d16434be435204ab511449a1713206773fc4d62087a8325e669dc09`).
- GO BH corrections, every declared PNG dimension and every companion PDF were verified.
- Fifty-three independent assertions passed in `results/tables/05_verification_summary.csv`.
- A final derived-workflow rerun produced 30 machine-readable CSV/CSV.GZ and PNG outputs byte-identically.

## Reproducibility entry points

- Frozen plan: `ANALYSIS_PLAN_FROZEN.md`
- Preflight: `scripts/00_preflight.R`
- Main workflow: `scripts/01_run_fig3_style_analysis.R`
- Independent verification: `scripts/02_verify_fig3_style.R`
- Byte-level reproducibility check: `scripts/03_repro_hash_check.R`
- Exact result tables: `results/tables/`
- Session information and reproducibility manifests: `logs/`

No manuscript, PowerPoint, Supplementary Information, Source Data workbook, upstream downloaded data or prior analysis folder was edited by this workflow.
