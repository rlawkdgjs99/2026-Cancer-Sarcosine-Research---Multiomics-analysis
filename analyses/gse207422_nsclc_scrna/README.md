# GSE207422 NSCLC single-cell RNA-seq

This module reconstructs a sarcosine functional-score analysis on GSE207422 with a frozen, marker-audited lineage annotation. It contains the primary cell-state workflow, Hallmark GSEA, lineage GO:BP analysis, response-by-lineage models, and full lineage-reannotation provenance.

The deposited inputs yield 92,053 cells at the reported mitochondrial-content threshold, whereas the source paper reports 92,031. Because the exact unpublished downstream object and 11-label mapping are unavailable, this workflow retains all objectively eligible cells and uses the independently audited frozen 14-lineage map. It does not claim exact object-level reproduction of the source figure.

Run the numbered scripts inside each submodule according to its analysis plan. Large expression objects and raw matrices are not included.
