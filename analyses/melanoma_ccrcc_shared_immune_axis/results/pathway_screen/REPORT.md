# TIGER melanoma and ccRCC: exact shared-pathway comparison

## Comparison fixed before interpretation

- TIGER melanoma: PRE-treatment tumors, sarcosine-degradation **High versus Low**; positive original NES denotes Degradation-High.
- ccRCC: tumors split at the median of **measured tumor Sarcosine**; the original contrast is High minus Low, so its NES was multiplied by -1 for this comparison.
- Therefore, positive harmonized NES denotes **TIGER Degradation-High** or **ccRCC measured Sarcosine-Low**.
- Only primary adjusted grouped GSEA results were used. Exact pathway IDs were matched within MSigDB 2026.1.Hs; BH q values remain the original within-collection values.
- This is an intersection analysis, not a pooled meta-analysis. The two exposures are biologically related hypotheses but are not the same measurement.

## Result counts

- Exact sets compared: 4649 (Hallmark 50, Reactome 1,027, GO:BP 3,572).
- Same target direction and BH q<0.05 in both cohorts: 62 exact sets (Hallmark 9, Reactome 17, GO:BP 36).

## IL-12-specific answer

Six exact IL-12-related sets were common to the two result universes. All six had the expected directional alignment, but **none was BH q<0.05 in both cohorts**. This means that the current data do not provide a replicated IL-12-specific pathway hit.

- GO:BP IL-12 production: TIGER NES 2.457, q=5.2e-10; ccRCC harmonized NES 1.249, q=0.476.
- Reactome IL-12-family signaling: TIGER NES 1.438, q=0.095; ccRCC harmonized NES 1.321, q=0.261.
- Reactome IL-12 signaling: TIGER NES 1.367, q=0.156; ccRCC harmonized NES 1.517, q=0.099.
- Reactome JAK/STAT expression after IL-12 stimulation is asymmetric: TIGER NES 1.008, q=0.606; ccRCC harmonized NES 1.791, q=0.018, with zero shared leading-edge genes.

The shared downstream signal is **type-II-IFN/IFNG**, not an IL-12-specific chain. Existing targeted analyses also show low IL12A/B abundance, absent IL-12 identity after family competition, or failed serial/context checks; therefore IL-12 should remain a secondary hypothesis rather than the lead common mechanism.

## Most defensible common pathway to pursue

### 1. Antigen presentation / TCR activation -> type-II-IFN/IFNG response (lead common axis)

This is the most coherent nonredundant cross-cohort axis because several distinct exact sets cover successive biological steps and pass BH q<0.05 in both cohorts:

- Antigen processing and presentation: TIGER NES 2.910, q=7.2e-23; ccRCC NES 1.661, q=0.019.
- Antigen cross-presentation: TIGER NES 2.101, q=2.1e-06; ccRCC NES 2.136, q=1e-05.
- TCR signaling: TIGER NES 2.399, q=8.9e-11; ccRCC NES 1.993, q=6.3e-05.
- IFN-gamma response: TIGER NES 3.222, q=1.6e-48; ccRCC NES 2.141, q=5.3e-09; 48 shared leading-edge genes.

This pattern supports prioritizing an **APC/antigen-presentation–TCR activation–IFNG immune-inflamed axis**. It does not prove a temporal signaling chain, and it does not establish that the same cell population generates every component.

### 2. TNF/NF-kappaB inflammatory program (strong secondary axis)

Hallmark TNF-alpha signaling via NF-kappaB is strongly aligned in both cohorts: TIGER NES 2.312, q=4.2e-13; ccRCC NES 2.186, q=8.3e-10; 42 shared leading-edge genes. This is a robust common program, but is less cell- and receptor-specific than the antigen/TCR axis.

### 3. Broad immune-inflamed state (supporting context)

IFN-alpha response, inflammatory response, and complement are also strict common pathways. They support a shared immune-inflamed state, but are broader and less useful as a single mechanistic route.

## Important exclusions and interpretation limits

- Do not describe the result as a replicated IL-12 mechanism; the exact IL-12 sets fail the two-cohort significance criterion.
- Do not describe CD28 as the shared cross-cohort route. CD28 was a leading TIGER candidate, but the ccRCC median analysis did not support CD28-specific enrichment.
- Bulk RNA cannot resolve APC, T-cell, and tumor-cell sources, and shared immune pathways may partly reflect immune-cell abundance.
- Cross-sectional associations do not establish causality, ligand secretion, receptor engagement, or phospho-STAT/NF-kappaB activation.

## Reader-facing outputs

- `figures/Fig_TIGER_Melanoma_ccRCC_Common_Pathway_Priorities.png`
- `tables/02_strict_aligned_both_BH_significant_pathways.csv`
- `tables/03_IL12_exact_pathway_comparison.csv`
- `tables/04_reader_priority_nonredundant_pathways.csv`
- `tables/05_priority_shared_leading_edge_genes.csv`

## Reproducibility

Input paths and MD5 hashes are recorded in `tables/07_input_manifest_md5.csv`. Session information is in `sessionInfo.txt`.
