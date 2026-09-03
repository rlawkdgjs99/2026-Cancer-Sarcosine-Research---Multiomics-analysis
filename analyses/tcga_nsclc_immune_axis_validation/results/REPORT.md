# TCGA NSCLC validation result

The prespecified four-program axis was evaluated in 1,017 LUAD/LUSC primary tumors. The covariate-adjusted model retained 988 tumors. Degradation-High/Low was defined within NSCLC using the median of `mean(z(SARDH), z(PIPOX))`.

## Exact-set GSEA

Positive NES denotes Degradation-High.

| Program | NES | Collection-wide BH q |
|---|---:|---:|
| APC cross-presentation | 1.435 | 0.0281 |
| TCR signaling | 1.977 | 1.04×10⁻⁸ |
| TNFR2-related noncanonical NF-κB expression program | 1.440 | 0.0311 |
| IFNG response | 2.051 | 1.42×10⁻¹⁷ |

All four exact pathways therefore replicate in the prespecified direction.

## Sample-level support

The adjusted High-minus-Low effects for the same exact-pathway scores were positive and BH-significant for all four programs. The three gene-disjoint adjacent-node associations were also positive and BH-significant: APC–TCR β=0.763, TCR–noncanonical NF-κB β=0.804, and noncanonical NF-κB–IFNG β=0.760.

These associations support a coordinated bulk-tumor immune program. They do not establish temporal ordering, TNFR2 receptor engagement, NFKB2/RELB protein activation, or same-cell signaling.
