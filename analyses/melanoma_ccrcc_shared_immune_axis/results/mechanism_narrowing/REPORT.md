# Cross-cohort common immune-program final narrowing

## Bottom line

The reproducible cross-cohort result is a **shared four-program enrichment intersection**: APC cross-presentation, TCR signaling, the TNFR2-related noncanonical NF-kappaB expression program, and IFNG response are each enriched in the target direction in both cohorts.
Within that shared program, **APC → TCR → noncanonical NF-kappaB → IFNG** is a direction-consistent exploratory branch association, not a replicated serial mechanism.
No complete exposure-linked serial route was reproduced in both cohorts: the first target-state-to-APC link was absent in ccRCC, and common NFKB2/RELB TF activation was not demonstrated.

## Target-state definitions

- TIGER melanoma: sarcosine-degradation High versus Low among 73 PRE-treatment tumors.
- ccRCC: measured tumor Sarcosine-Low versus High by the prespecified median split among 100 matched tumors.
- Positive harmonized NES and positive model effects denote those target states.

## Broad screen

- Exact common-pathway gates were evaluated for 11 prespecified candidates.
- 6 candidates were direction-aligned and BH q<0.05 in both cohorts.
- The displayed shared downstream axis passed the exact-set gate at all four nodes: APC cross-presentation (shared leading edge n=12), TCR signaling (n=14), noncanonical NF-kappaB (n=3; NFKB2, RELB, TRAF3), and IFNG response (n=48).
- Canonical NF-kappaB: TIGER NES 2.31, q 0.000000000000423; ccRCC NES 2.19, q 0.000000000834; shared leading-edge genes 42.
- Noncanonical NF-kappaB: TIGER NES 1.97, q 0.000142; ccRCC NES 2.12, q 0.0000463; shared leading-edge genes 3.
- Specific CD28 co-stimulation failed the two-cohort gate (TIGER q 0.0111; ccRCC q 0.69). The broad CD28-family set was retained only as a mixed regulatory sensitivity node.
- IL-12/STAT4 failed the two-cohort gate (TIGER q 0.0954; ccRCC q 0.261).

## Gene-disjoint route design

APC cross-presentation, TCR, broad CD28-family regulation, canonical NF-kappaB, noncanonical NF-kappaB, and the IFNG-response outcome were scored from MSigDB 2026.1.Hs. IFNG-response genes were removed from every predictor, and genes shared by two or more predictor modules were removed from every predictor. All retained predictor/outcome scores therefore have zero pairwise gene overlap.

## Direct sample-level group separation

- Melanoma gene-disjoint route scores all showed positive adjusted Degradation-High effects: APC beta 0.38, TCR beta 0.427, noncanonical NF-kappaB beta 0.516, and IFNG-response beta 0.497; all BH q<0.002.
- ccRCC gene-disjoint route scores did not separate Sarcosine-Low from High: APC beta -0.00765, TCR beta 0.0169, noncanonical NF-kappaB beta -0.0474, and IFNG-response beta 0.0702; all BH q=0.947.
This explains the weak visual group separation in the ccRCC heatmap. The cross-cohort overlap is restricted to full exact-set ranked-list GSEA and downstream branch-level covariance; it is not a replicated ccRCC sample-level route contrast.

## Primary grouped serial products for the leading route

- Melanoma: 0.025 [0.00768, 0.0499].
- ccRCC: -0.000975 [-0.0385, 0.026].
- Complete forward serial replication in both cohorts: **FAIL**.
- Reverse-order products are retained in the source table and are treated as an ambiguity diagnostic, not as proof of reverse biology.

## Branch competition

- Melanoma leading-branch beta: 0.176, BH q 0.0536.
- ccRCC leading-branch beta: 0.363, BH q 0.000226.
- Two-cohort fixed-effect branch estimate: beta 0.268 [0.151, 0.385], p 0.00000741, I2 59.2%.
- Canonical NF-kappaB was direction-discordant after joint competition; its fixed-effect beta was -0.0154, p 0.709, I2 95.7%. This is why the strong canonical enrichment was not promoted to the final common branch.

## TF evidence

SPI1, NFKB1, RELA, and NFKB2 activities were inferred from cross-cohort-consensus signed DoRothEA A/B target edges after excluding every route-score and IFNG-response gene. NFKB2 retained 10 targets and met the prespecified minimum of 10; RELB retained only 5 targets and its activity was therefore not estimated.
- NFKB2 inferred-activity target-state effect: melanoma beta -0.316, q 0.0202; ccRCC beta -0.148, q 0.897.
- NFKB2 activity-to-noncanonical-program association: melanoma beta -0.00313, q 0.977; ccRCC beta 0.0376, q 0.74.
- NFKB2 transcript target-state effect: melanoma beta 0.484, q 0.000000918; ccRCC beta 0.132, q 0.233.
- RELB transcript target-state effect: melanoma beta 0.508, q 0.000000918; ccRCC beta 0.207, q 0.137.
SPI1 remains an APC-context check, while NFKB1/RELA remain the prespecified canonical-branch check. TF activities and focal transcripts are corroborative diagnostics and are not inserted as causal nodes.

## Interpretation boundary

The two exposures are not identical measurements: TIGER uses a SARDH/PIPOX transcriptional degradation score, whereas ccRCC uses measured tumor Sarcosine. Bulk RNA cannot identify ligand secretion, receptor engagement, phosphorylation, cell of origin, or temporal ordering. Named T-cell/APC context-adjusted models are sensitivity analyses because adjustment can remove either confounding or true immune-composition mediation.

## Reader-facing outputs

Use only the following two figures as the connected standalone reader sequence:

1. `figures_final_reader/Fig_01_What_Is_Shared_Across_Cohorts.png`: directly presents the selected common candidate axis with paired cohort NES/q values and shared leading-edge evidence at every node.
2. `figures_final_reader/Fig_02_Cohort_Specific_Group_Contrasts.png`: shows why the shared enrichment cannot be promoted to a replicated sample-level serial mechanism.

The legacy `figures/` directory contains extended diagnostics and provenance. It is not a standalone narrative sequence. Full clinical, continuous, reverse-order, and named-cell-context results remain in the tables.
