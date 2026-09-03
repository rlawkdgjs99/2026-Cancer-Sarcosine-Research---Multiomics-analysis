# Epithelial- and CAF-restricted sarcosine-ratio quartile analysis

Date: 26 August 2026
Dataset: GSE207422 single-cell RNA-seq
Status: Analysis complete and independently verified; not installed in the manuscript or figure deck.

## Question

Within Epithelial cells and CAFs separately, which GO Biological Process terms are enriched among genes expressed more highly in cells in the upper versus lower quartile of the sarcosine production/degradation transcriptional score ratio?

Production was defined as `GNMT + DMGDH`; degradation was defined as `SARDH + PIPOX`. The ratio is the previously validated shifted AddModuleScore ratio. It is a transcriptional score balance, not a measurement of sarcosine concentration, enzyme activity, or metabolic flux.

## Frozen design

- Only cells with a finite ratio were eligible.
- Q1 and Q3 were computed separately within each lineage using R type-7 quantiles.
- Bottom group: ratio <= Q1; top group: ratio >= Q3; the middle 50% was excluded.
- Differential expression: Seurat `FindMarkers`, top versus bottom, Wilcoxon test, normalized RNA data, `logfc.threshold = 0.1`, `min.pct = 0.01`, no downsampling.
- Directional enrichment inputs: BH-adjusted DEG q < 0.05 and absolute average log2 fold change >= 0.25.
- GO Biological Process over-representation: `limma::goana`; BH correction was performed independently for each lineage and DEG direction.
- Displayed terms: GO BP size 10–500 genes, DE overlap > 0, BH q < 0.05, top 10 by q for each lineage/direction.
- The plotted x-axis is fold enrichment, `(DE/N)/(directional mapped DEG/universe)`. Point size is the number of directional DE genes in the GO term; colour is `-log10(BH q)`.

The full locked specification and frozen hashes are in `ANALYSIS_PLAN_FROZEN.md`.

## Cell selection

| Lineage | Finite-ratio cells | Q1 | Q3 | Bottom Q1 cells | Top Q4 cells | Patients in Q1 | Patients in Q4 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Epithelial | 11,912 | 0.376174 | 0.536637 | 2,978 | 2,978 | 15 | 15 |
| CAF | 1,102 | 0.304701 | 0.501040 | 276 | 276 | 15 | 14 |

No cells were tied exactly at either quartile threshold.

## DEG inputs to GO enrichment

| Lineage | Tested genes | Mapped universe | Q4-up symbols / Entrez | Q1-up symbols / Entrez |
|---|---:|---:|---:|---:|
| Epithelial | 14,226 | 12,315 | 5,486 / 4,997 | 1,554 / 1,366 |
| CAF | 12,900 | 11,650 | 134 / 132 | 58 / 57 |

## Main findings

### Epithelial cells

- The high-ratio Q4 program was enriched for chromatin remodelling, cilium assembly/organization, microtubule cytoskeleton and cell-projection organization, and DNA replication. The strongest displayed term was chromatin remodelling (265 overlapping DE genes; fold enrichment 1.464; BH q = 1.04e-13).
- The low-ratio Q1 program was enriched for translation and protein biosynthesis, together with several cell-fusion/syncytium-related terms. Translation and protein biosynthetic process each included 118 directional DE genes (fold enrichment 2.235; BH q = 4.39e-14).

### CAFs

- The high-ratio Q4 program was strongly enriched for extracellular-matrix and extracellular-structure organization and collagen-related processes. Extracellular matrix organization, extracellular structure organization, and external encapsulating structure organization each included 28 directional DE genes (fold enrichment 12.673; BH q = 8.79e-20).
- The low-ratio Q1 program was enriched for angiogenesis, blood-vessel/vasculature development and morphogenesis, endothelial differentiation/migration, and sprouting angiogenesis. Angiogenesis included 19 directional DE genes (fold enrichment 10.847; BH q = 2.63e-11).

## Correct interpretation

The defensible descriptive conclusion is:

> Within the reannotated GSE207422 data, epithelial and CAF cells at opposite extremes of the sarcosine production/degradation transcriptional-score ratio showed different gene-program enrichments. High-ratio epithelial cells were associated mainly with chromatin/ciliary/cell-projection programs, whereas low-ratio epithelial cells showed translational programs. In CAFs, high-ratio cells showed extracellular-matrix/collagen programs, while low-ratio cells showed angiogenic and vascular-development programs.

This does **not** demonstrate that sarcosine metabolism caused these states, that the cells had higher or lower sarcosine flux, or that the findings generalize at the patient level.

## Critical limitation

Cells are nested within patients and are not independent biological replicates. Patient contributions were markedly unequal:

- Epithelial Q4: one patient contributed 1,495/2,978 cells (50.2%).
- Epithelial Q1: the largest patient contributed 711/2,978 cells (23.9%).
- CAF Q4: one patient contributed 138/276 cells (50.0%).
- CAF Q1: the largest patient contributed 107/276 cells (38.8%).

Therefore the very small cell-level DEG and GO q-values must not be presented as patient-level evidence. Before use as a main or supplementary manuscript panel, a patient-aware sensitivity analysis should be added, ideally using within-patient contrasts or a pseudobulk/mixed-effects design. The current figure is suitable for exploratory review, not yet for a causal or patient-generalized claim.

## Verification

- Frozen input SHA-256 values matched before analysis.
- An independent verification script passed 20/20 checks, including quartile membership, DEG-direction rules, BH q-values, enrichment fractions, fold enrichment, and displayed top-term selection.
- A second full clean run reproduced 14/14 core CSV and PNG files byte-for-byte (identical SHA-256).
- Publication PNGs are 600 dpi with embedded sRGB profiles.

## Key outputs

- Combined figure: `results/figures_publication/Fig_Epithelial_CAF_ratio_Q4_vs_Q1_GO_BP_combined.png`
- Epithelial figure: `results/figures_publication/Fig_Epithelial_ratio_Q4_vs_Q1_GO_BP.png`
- CAF figure: `results/figures_publication/Fig_CAF_ratio_Q4_vs_Q1_GO_BP.png`
- All DEG results: `results/tables/02_lineage_quartile_FindMarkers_all.csv`
- All directional GO results: `results/tables/03_lineage_quartile_GO_BP_directional_all.csv`
- Displayed GO terms: `results/tables/plotdata_top_GO_BP_terms.csv`
- Independent verification: `results/tables/04_verification_summary.csv`
- Reproducibility comparison: `logs/03_reproducibility_comparison.csv`
