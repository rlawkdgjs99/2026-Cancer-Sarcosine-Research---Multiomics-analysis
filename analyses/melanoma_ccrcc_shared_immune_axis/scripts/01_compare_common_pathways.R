#!/usr/bin/env Rscript

# Cross-cohort comparison of the primary adjusted pathway results from:
#   1) TIGER melanoma PRE tumors, sarcosine-degradation High vs Low
#   2) ccRCC tumors, measured Sarcosine High vs Low
#
# Direction harmonization used throughout:
#   positive harmonized NES = TIGER Degradation-High or ccRCC Sarcosine-Low.
#
# This is an exact MSigDB gene-set comparison, not a pooled meta-analysis.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_arg)) stop("Run this file with Rscript.")
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
module_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
input_dir <- Sys.getenv("SARCO_SHARED_AXIS_INPUT_DIR", unset = file.path(module_dir, "inputs"))

input_path <- function(env_name, default_name) {
  value <- Sys.getenv(env_name, unset = file.path(input_dir, default_name))
  if (!file.exists(value)) stop("Missing input ", env_name, ": ", value)
  normalizePath(value, mustWork = TRUE)
}

ccrcc_file <- input_path(
  "SARCO_CCRCC_GSEA",
  "ccrcc_GSEA_Hallmark_Reactome_GOBP_adjusted_High_vs_Low.csv"
)
tiger_hallmark_file <- input_path(
  "SARCO_TIGER_HALLMARK_GSEA", "tiger_fgsea_Hallmark_adjusted_High_vs_Low.csv"
)
tiger_reactome_file <- input_path(
  "SARCO_TIGER_REACTOME_GSEA", "tiger_fgsea_Reactome_adjusted_High_vs_Low.csv"
)
tiger_gobp_file <- input_path(
  "SARCO_TIGER_GOBP_GSEA", "tiger_fgsea_GO_BP_adjusted_High_vs_Low.csv"
)

out_dir <- file.path(module_dir, "results", "pathway_screen")
table_dir <- file.path(out_dir, "tables")
figure_dir <- file.path(out_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

required_columns <- c(
  "pathway", "collection", "pval", "padj", "ES", "NES", "size",
  "leadingEdge", "model", "direction", "description", "MSigDB_version"
)

read_checked <- function(path, source_label) {
  x <- read_csv(path, show_col_types = FALSE, progress = FALSE)
  missing_columns <- setdiff(required_columns, names(x))
  if (length(missing_columns)) {
    stop(source_label, " is missing columns: ", paste(missing_columns, collapse = ", "))
  }
  if (anyDuplicated(paste(x$collection, x$pathway, sep = "::"))) {
    stop(source_label, " contains duplicated collection/pathway IDs")
  }
  if (any(!is.finite(x$NES)) || any(!is.finite(x$padj))) {
    stop(source_label, " contains non-finite NES or adjusted P values")
  }
  x
}

ccrcc <- read_checked(ccrcc_file, "ccRCC")
tiger_hallmark <- read_checked(tiger_hallmark_file, "TIGER Hallmark")
tiger_reactome <- read_checked(tiger_reactome_file, "TIGER Reactome")
tiger_gobp <- read_checked(tiger_gobp_file, "TIGER GO:BP")
tiger <- bind_rows(tiger_hallmark, tiger_reactome, tiger_gobp)

stopifnot(
  nrow(tiger_hallmark) == 50L,
  nrow(tiger_reactome) == 1035L,
  nrow(tiger_gobp) == 3630L,
  all(unique(tiger$MSigDB_version) == "2026.1.Hs"),
  all(unique(ccrcc$MSigDB_version) == "2026.1.Hs"),
  any(grepl("High", tiger$model, fixed = TRUE)),
  any(grepl("Sarcosine", ccrcc$model, fixed = TRUE))
)

split_genes <- function(x) {
  if (is.na(x) || !nzchar(x)) return(character())
  unique(trimws(strsplit(x, ";", fixed = TRUE)[[1]]))
}

common <- tiger %>%
  select(
    pathway, collection,
    tiger_NES = NES, tiger_pval = pval, tiger_q = padj,
    tiger_size = size, tiger_leading_edge = leadingEdge,
    tiger_model = model, tiger_direction = direction,
    description, MSigDB_version
  ) %>%
  inner_join(
    ccrcc %>%
      select(
        pathway, collection,
        ccrcc_NES_high_minus_low = NES,
        ccrcc_pval = pval, ccrcc_q = padj,
        ccrcc_size = size, ccrcc_leading_edge = leadingEdge,
        ccrcc_model = model, ccrcc_direction = direction
      ),
    by = c("pathway", "collection")
  ) %>%
  rowwise() %>%
  mutate(
    tiger_harmonized_NES = tiger_NES,
    ccrcc_harmonized_NES = -ccrcc_NES_high_minus_low,
    target_direction_aligned = tiger_harmonized_NES > 0 & ccrcc_harmonized_NES > 0,
    opposite_direction_aligned = tiger_harmonized_NES < 0 & ccrcc_harmonized_NES < 0,
    both_BH_q_lt_0_05 = tiger_q < 0.05 & ccrcc_q < 0.05,
    strict_target_common = target_direction_aligned & both_BH_q_lt_0_05,
    max_q = max(tiger_q, ccrcc_q),
    min_abs_harmonized_NES = min(abs(tiger_harmonized_NES), abs(ccrcc_harmonized_NES)),
    shared_leading_edge_genes = paste(
      sort(intersect(split_genes(tiger_leading_edge), split_genes(ccrcc_leading_edge))),
      collapse = ";"
    ),
    shared_leading_edge_n = length(
      intersect(split_genes(tiger_leading_edge), split_genes(ccrcc_leading_edge))
    ),
    leading_edge_union_n = length(
      union(split_genes(tiger_leading_edge), split_genes(ccrcc_leading_edge))
    ),
    leading_edge_jaccard = if_else(
      leading_edge_union_n > 0,
      shared_leading_edge_n / leading_edge_union_n,
      0
    )
  ) %>%
  ungroup() %>%
  arrange(collection, pathway)

stopifnot(
  nrow(common) == 4649L,
  sum(common$strict_target_common) == 62L,
  all(common$tiger_harmonized_NES == common$tiger_NES),
  all(common$ccrcc_harmonized_NES == -common$ccrcc_NES_high_minus_low)
)

strict_common <- common %>%
  filter(strict_target_common) %>%
  arrange(max_q, desc(leading_edge_jaccard), desc(min_abs_harmonized_NES))

il12_pattern <- "INTERLEUKIN_12|IL_12|IL12"
il12 <- common %>%
  filter(grepl(il12_pattern, pathway)) %>%
  mutate(
    tiger_significance = if_else(tiger_q < 0.05, "BH q<0.05", "BH q>=0.05"),
    ccrcc_significance = if_else(ccrcc_q < 0.05, "BH q<0.05", "BH q>=0.05")
  ) %>%
  arrange(collection, pathway)

stopifnot(nrow(il12) == 6L, sum(il12$strict_target_common) == 0L)

priority_key <- tribble(
  ~pathway, ~theme, ~display_label, ~plot_order,
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION", "Antigen presentation/T-cell interface", "Antigen processing and presentation", 1,
  "REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION", "Antigen presentation/T-cell interface", "Antigen cross-presentation", 2,
  "REACTOME_TCR_SIGNALING", "Antigen presentation/T-cell interface", "TCR signaling", 3,
  "GOBP_REGULATION_OF_T_CELL_ACTIVATION", "Antigen presentation/T-cell interface", "Regulation of T-cell activation", 4,
  "GOBP_TYPE_II_INTERFERON_PRODUCTION", "Type-II-IFN/IFNG", "Type-II IFN production", 5,
  "REACTOME_INTERFERON_GAMMA_SIGNALING", "Type-II-IFN/IFNG", "IFN-gamma signaling", 6,
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "Type-II-IFN/IFNG", "IFN-gamma response", 7,
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "TNF/NF-kappaB inflammatory", "TNF-alpha signaling via NF-kappaB", 8,
  "HALLMARK_INFLAMMATORY_RESPONSE", "TNF/NF-kappaB inflammatory", "Inflammatory response", 9,
  "HALLMARK_INTERFERON_ALPHA_RESPONSE", "Broad immune-inflamed state", "IFN-alpha response", 10,
  "HALLMARK_COMPLEMENT", "Broad immune-inflamed state", "Complement", 11
)

priority <- priority_key %>%
  left_join(common, by = "pathway") %>%
  arrange(plot_order)

if (any(is.na(priority$collection)) || !all(priority$strict_target_common)) {
  stop("At least one prespecified reader-priority pathway is missing or not strict-common")
}

counts <- common %>%
  group_by(collection) %>%
  summarise(
    exact_sets_compared = n(),
    target_direction_aligned = sum(target_direction_aligned),
    both_BH_q_lt_0_05 = sum(both_BH_q_lt_0_05),
    strict_target_common = sum(strict_target_common),
    .groups = "drop"
  ) %>%
  bind_rows(
    common %>%
      summarise(
        collection = "All collections",
        exact_sets_compared = n(),
        target_direction_aligned = sum(target_direction_aligned),
        both_BH_q_lt_0_05 = sum(both_BH_q_lt_0_05),
        strict_target_common = sum(strict_target_common)
      )
  )

priority_gene_rows <- priority %>%
  select(pathway, theme, display_label, shared_leading_edge_genes) %>%
  separate_longer_delim(shared_leading_edge_genes, delim = ";") %>%
  filter(shared_leading_edge_genes != "") %>%
  rename(shared_leading_edge_gene = shared_leading_edge_genes)

input_manifest <- tibble(
  source = c("ccRCC combined GSEA", "TIGER Hallmark", "TIGER Reactome", "TIGER GO:BP"),
  path = c(ccrcc_file, tiger_hallmark_file, tiger_reactome_file, tiger_gobp_file),
  md5 = unname(tools::md5sum(path)),
  rows = c(nrow(ccrcc), nrow(tiger_hallmark), nrow(tiger_reactome), nrow(tiger_gobp)),
  MSigDB_version = "2026.1.Hs"
)

write_csv(common, file.path(table_dir, "01_all_exact_common_pathway_sets.csv"))
write_csv(strict_common, file.path(table_dir, "02_strict_aligned_both_BH_significant_pathways.csv"))
write_csv(il12, file.path(table_dir, "03_IL12_exact_pathway_comparison.csv"))
write_csv(priority, file.path(table_dir, "04_reader_priority_nonredundant_pathways.csv"))
write_csv(priority_gene_rows, file.path(table_dir, "05_priority_shared_leading_edge_genes.csv"))
write_csv(counts, file.path(table_dir, "06_comparison_counts.csv"))
write_csv(input_manifest, file.path(table_dir, "07_input_manifest_md5.csv"))

format_q <- function(x) {
  ifelse(x < 0.001, format(x, scientific = TRUE, digits = 2), sprintf("%.3f", x))
}

priority_long <- priority %>%
  transmute(
    pathway, theme, display_label, plot_order,
    shared_leading_edge_n, leading_edge_jaccard,
    `TIGER melanoma\nDegradation-High` = tiger_harmonized_NES,
    `ccRCC\nSarcosine-Low` = ccrcc_harmonized_NES,
    tiger_q, ccrcc_q
  ) %>%
  pivot_longer(
    cols = c(`TIGER melanoma\nDegradation-High`, `ccRCC\nSarcosine-Low`),
    names_to = "cohort_target_group", values_to = "harmonized_NES"
  ) %>%
  mutate(
    q = if_else(grepl("TIGER", cohort_target_group), tiger_q, ccrcc_q),
    cell_label = paste0("NES ", sprintf("%.2f", harmonized_NES), "\nq ", format_q(q)),
    display_label = factor(display_label, levels = rev(priority_key$display_label)),
    cohort_target_group = factor(
      cohort_target_group,
      levels = c("TIGER melanoma\nDegradation-High", "ccRCC\nSarcosine-Low")
    )
  )

il12_labels <- c(
  GOBP_INTERLEUKIN_12_PRODUCTION = "IL-12 production",
  GOBP_NEGATIVE_REGULATION_OF_INTERLEUKIN_12_PRODUCTION = "Negative regulation of IL-12 production",
  GOBP_POSITIVE_REGULATION_OF_INTERLEUKIN_12_PRODUCTION = "Positive regulation of IL-12 production",
  REACTOME_GENE_AND_PROTEIN_EXPRESSION_BY_JAK_STAT_SIGNALING_AFTER_INTERLEUKIN_12_STIMULATION = "JAK/STAT expression after IL-12 stimulation",
  REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING = "IL-12-family signaling",
  REACTOME_INTERLEUKIN_12_SIGNALING = "IL-12 signaling"
)

il12_long <- il12 %>%
  mutate(display_label = unname(il12_labels[pathway])) %>%
  transmute(
    pathway, display_label, shared_leading_edge_n, leading_edge_jaccard,
    `TIGER melanoma\nDegradation-High` = tiger_harmonized_NES,
    `ccRCC\nSarcosine-Low` = ccrcc_harmonized_NES,
    tiger_q, ccrcc_q
  ) %>%
  pivot_longer(
    cols = c(`TIGER melanoma\nDegradation-High`, `ccRCC\nSarcosine-Low`),
    names_to = "cohort_target_group", values_to = "harmonized_NES"
  ) %>%
  mutate(
    q = if_else(grepl("TIGER", cohort_target_group), tiger_q, ccrcc_q),
    cell_label = paste0("NES ", sprintf("%.2f", harmonized_NES), "\nq ", format_q(q)),
    display_label = factor(display_label, levels = rev(unname(il12_labels))),
    cohort_target_group = factor(
      cohort_target_group,
      levels = c("TIGER melanoma\nDegradation-High", "ccRCC\nSarcosine-Low")
    )
  )

fill_scale <- scale_fill_gradient2(
  low = "#2A9D8F", mid = "#F7F7F7", high = "#C83E3E",
  midpoint = 0, limits = c(-3.3, 3.3), oob = squish,
  name = "Harmonized\nNES"
)

base_theme <- theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_text(colour = "black", size = 10),
    axis.text.x = element_text(colour = "black", face = "bold", size = 10),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, colour = "#4D4D4D"),
    legend.position = "right"
  )

p_priority <- ggplot(priority_long, aes(cohort_target_group, display_label, fill = harmonized_NES)) +
  geom_tile(colour = "white", linewidth = 1.1) +
  geom_text(aes(label = cell_label), size = 3.1, lineheight = 0.95) +
  fill_scale +
  labs(
    title = "a  Nonredundant pathways significant in both cohorts",
    subtitle = "Each cell uses the original within-collection BH q; all displayed sets pass q<0.05 in both cohorts"
  ) +
  base_theme

p_il12 <- ggplot(il12_long, aes(cohort_target_group, display_label, fill = harmonized_NES)) +
  geom_tile(aes(colour = q < 0.05), linewidth = 1.2) +
  geom_text(aes(label = cell_label), size = 3.0, lineheight = 0.95) +
  fill_scale +
  scale_colour_manual(
    values = c(`TRUE` = "#111111", `FALSE` = "#D0D0D0"),
    labels = c(`TRUE` = "BH q<0.05", `FALSE` = "BH q>=0.05"),
    name = "Cell border"
  ) +
  labs(
    title = "b  IL-12-related exact gene sets",
    subtitle = "All six lean in the target direction, but none reaches BH q<0.05 in both cohorts"
  ) +
  base_theme

p_overlap <- ggplot(
  priority,
  aes(leading_edge_jaccard, factor(display_label, levels = rev(priority_key$display_label)))
) +
  geom_segment(aes(x = 0, xend = leading_edge_jaccard, yend = factor(display_label, levels = rev(priority_key$display_label))), colour = "#BBBBBB") +
  geom_point(aes(size = shared_leading_edge_n, colour = theme), alpha = 0.9) +
  geom_text(aes(label = paste0("n=", shared_leading_edge_n)), hjust = -0.25, size = 3) +
  scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 0.52), expand = expansion(mult = c(0, 0.12))) +
  scale_size_continuous(range = c(3, 9), name = "Shared\nleading-edge n") +
  scale_colour_manual(
    values = c(
      "Antigen presentation/T-cell interface" = "#7A5195",
      "Type-II-IFN/IFNG" = "#D45087",
      "TNF/NF-kappaB inflammatory" = "#F95D6A",
      "Broad immune-inflamed state" = "#2A9D8F"
    ),
    name = "Theme"
  ) +
  labs(
    title = "c  Shared leading-edge genes",
    subtitle = "Jaccard overlap of the two cohort-specific leading edges",
    x = "Leading-edge Jaccard overlap",
    y = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(colour = "black", size = 9),
    axis.title.x = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, colour = "#4D4D4D"),
    legend.position = "right"
  )

combined_figure <- p_priority / (p_il12 | p_overlap) +
  plot_layout(heights = c(1.35, 1), guides = "collect") +
  plot_annotation(
    title = "Shared immune programs: melanoma Degradation-High and ccRCC measured Sarcosine-Low",
    subtitle = paste0(
      "Exact MSigDB 2026.1.Hs set comparison of primary adjusted GSEA results. ",
      "Positive harmonized NES denotes the named target group; this is not a pooled meta-analysis."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 19),
      plot.subtitle = element_text(size = 11, colour = "#4D4D4D"),
      plot.margin = margin(10, 15, 10, 10)
    )
  ) & theme(legend.position = "right")

ggsave(
  file.path(figure_dir, "Fig_TIGER_Melanoma_ccRCC_Common_Pathway_Priorities.png"),
  combined_figure, width = 16, height = 15, units = "in", dpi = 300, bg = "white"
)
ggsave(
  file.path(figure_dir, "Fig_TIGER_Melanoma_ccRCC_Common_Pathway_Priorities.pdf"),
  combined_figure, width = 16, height = 15, units = "in", device = cairo_pdf, bg = "white"
)

lookup <- function(id, field) priority[[field]][priority$pathway == id][1]
lookup_il12 <- function(id, field) il12[[field]][il12$pathway == id][1]

report_lines <- c(
  "# TIGER melanoma and ccRCC: exact shared-pathway comparison",
  "",
  "## Comparison fixed before interpretation",
  "",
  "- TIGER melanoma: PRE-treatment tumors, sarcosine-degradation **High versus Low**; positive original NES denotes Degradation-High.",
  "- ccRCC: tumors split at the median of **measured tumor Sarcosine**; the original contrast is High minus Low, so its NES was multiplied by -1 for this comparison.",
  "- Therefore, positive harmonized NES denotes **TIGER Degradation-High** or **ccRCC measured Sarcosine-Low**.",
  "- Only primary adjusted grouped GSEA results were used. Exact pathway IDs were matched within MSigDB 2026.1.Hs; BH q values remain the original within-collection values.",
  "- This is an intersection analysis, not a pooled meta-analysis. The two exposures are biologically related hypotheses but are not the same measurement.",
  "",
  "## Result counts",
  "",
  paste0("- Exact sets compared: ", nrow(common), " (Hallmark 50, Reactome 1,027, GO:BP 3,572)."),
  paste0("- Same target direction and BH q<0.05 in both cohorts: ", nrow(strict_common), " exact sets ",
         "(Hallmark ", sum(strict_common$collection == "Hallmark"), ", Reactome ",
         sum(strict_common$collection == "Reactome"), ", GO:BP ", sum(strict_common$collection == "GO:BP"), ")."),
  "",
  "## IL-12-specific answer",
  "",
  paste0("Six exact IL-12-related sets were common to the two result universes. All six had the expected directional alignment, but **none was BH q<0.05 in both cohorts**. ",
         "This means that the current data do not provide a replicated IL-12-specific pathway hit."),
  "",
  paste0("- GO:BP IL-12 production: TIGER NES ", sprintf("%.3f", lookup_il12("GOBP_INTERLEUKIN_12_PRODUCTION", "tiger_harmonized_NES")),
         ", q=", format_q(lookup_il12("GOBP_INTERLEUKIN_12_PRODUCTION", "tiger_q")),
         "; ccRCC harmonized NES ", sprintf("%.3f", lookup_il12("GOBP_INTERLEUKIN_12_PRODUCTION", "ccrcc_harmonized_NES")),
         ", q=", format_q(lookup_il12("GOBP_INTERLEUKIN_12_PRODUCTION", "ccrcc_q")), "."),
  paste0("- Reactome IL-12-family signaling: TIGER NES ", sprintf("%.3f", lookup_il12("REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING", "tiger_harmonized_NES")),
         ", q=", format_q(lookup_il12("REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING", "tiger_q")),
         "; ccRCC harmonized NES ", sprintf("%.3f", lookup_il12("REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING", "ccrcc_harmonized_NES")),
         ", q=", format_q(lookup_il12("REACTOME_INTERLEUKIN_12_FAMILY_SIGNALING", "ccrcc_q")), "."),
  paste0("- Reactome IL-12 signaling: TIGER NES ", sprintf("%.3f", lookup_il12("REACTOME_INTERLEUKIN_12_SIGNALING", "tiger_harmonized_NES")),
         ", q=", format_q(lookup_il12("REACTOME_INTERLEUKIN_12_SIGNALING", "tiger_q")),
         "; ccRCC harmonized NES ", sprintf("%.3f", lookup_il12("REACTOME_INTERLEUKIN_12_SIGNALING", "ccrcc_harmonized_NES")),
         ", q=", format_q(lookup_il12("REACTOME_INTERLEUKIN_12_SIGNALING", "ccrcc_q")), "."),
  paste0("- Reactome JAK/STAT expression after IL-12 stimulation is asymmetric: TIGER NES ",
         sprintf("%.3f", lookup_il12("REACTOME_GENE_AND_PROTEIN_EXPRESSION_BY_JAK_STAT_SIGNALING_AFTER_INTERLEUKIN_12_STIMULATION", "tiger_harmonized_NES")),
         ", q=", format_q(lookup_il12("REACTOME_GENE_AND_PROTEIN_EXPRESSION_BY_JAK_STAT_SIGNALING_AFTER_INTERLEUKIN_12_STIMULATION", "tiger_q")),
         "; ccRCC harmonized NES ",
         sprintf("%.3f", lookup_il12("REACTOME_GENE_AND_PROTEIN_EXPRESSION_BY_JAK_STAT_SIGNALING_AFTER_INTERLEUKIN_12_STIMULATION", "ccrcc_harmonized_NES")),
         ", q=", format_q(lookup_il12("REACTOME_GENE_AND_PROTEIN_EXPRESSION_BY_JAK_STAT_SIGNALING_AFTER_INTERLEUKIN_12_STIMULATION", "ccrcc_q")),
         ", with zero shared leading-edge genes."),
  "",
  "The shared downstream signal is **type-II-IFN/IFNG**, not an IL-12-specific chain. Existing targeted analyses also show low IL12A/B abundance, absent IL-12 identity after family competition, or failed serial/context checks; therefore IL-12 should remain a secondary hypothesis rather than the lead common mechanism.",
  "",
  "## Most defensible common pathway to pursue",
  "",
  "### 1. Antigen presentation / TCR activation -> type-II-IFN/IFNG response (lead common axis)",
  "",
  "This is the most coherent nonredundant cross-cohort axis because several distinct exact sets cover successive biological steps and pass BH q<0.05 in both cohorts:",
  "",
  paste0("- Antigen processing and presentation: TIGER NES ", sprintf("%.3f", lookup("GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION", "tiger_harmonized_NES")),
         ", q=", format_q(lookup("GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION", "tiger_q")),
         "; ccRCC NES ", sprintf("%.3f", lookup("GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION", "ccrcc_harmonized_NES")),
         ", q=", format_q(lookup("GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION", "ccrcc_q")), "."),
  paste0("- Antigen cross-presentation: TIGER NES ", sprintf("%.3f", lookup("REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION", "tiger_harmonized_NES")),
         ", q=", format_q(lookup("REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION", "tiger_q")),
         "; ccRCC NES ", sprintf("%.3f", lookup("REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION", "ccrcc_harmonized_NES")),
         ", q=", format_q(lookup("REACTOME_ANTIGEN_PROCESSING_CROSS_PRESENTATION", "ccrcc_q")), "."),
  paste0("- TCR signaling: TIGER NES ", sprintf("%.3f", lookup("REACTOME_TCR_SIGNALING", "tiger_harmonized_NES")),
         ", q=", format_q(lookup("REACTOME_TCR_SIGNALING", "tiger_q")),
         "; ccRCC NES ", sprintf("%.3f", lookup("REACTOME_TCR_SIGNALING", "ccrcc_harmonized_NES")),
         ", q=", format_q(lookup("REACTOME_TCR_SIGNALING", "ccrcc_q")), "."),
  paste0("- IFN-gamma response: TIGER NES ", sprintf("%.3f", lookup("HALLMARK_INTERFERON_GAMMA_RESPONSE", "tiger_harmonized_NES")),
         ", q=", format_q(lookup("HALLMARK_INTERFERON_GAMMA_RESPONSE", "tiger_q")),
         "; ccRCC NES ", sprintf("%.3f", lookup("HALLMARK_INTERFERON_GAMMA_RESPONSE", "ccrcc_harmonized_NES")),
         ", q=", format_q(lookup("HALLMARK_INTERFERON_GAMMA_RESPONSE", "ccrcc_q")),
         "; 48 shared leading-edge genes."),
  "",
  "This pattern supports prioritizing an **APC/antigen-presentation–TCR activation–IFNG immune-inflamed axis**. It does not prove a temporal signaling chain, and it does not establish that the same cell population generates every component.",
  "",
  "### 2. TNF/NF-kappaB inflammatory program (strong secondary axis)",
  "",
  paste0("Hallmark TNF-alpha signaling via NF-kappaB is strongly aligned in both cohorts: TIGER NES ",
         sprintf("%.3f", lookup("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "tiger_harmonized_NES")),
         ", q=", format_q(lookup("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "tiger_q")),
         "; ccRCC NES ", sprintf("%.3f", lookup("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "ccrcc_harmonized_NES")),
         ", q=", format_q(lookup("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "ccrcc_q")),
         "; 42 shared leading-edge genes. This is a robust common program, but is less cell- and receptor-specific than the antigen/TCR axis."),
  "",
  "### 3. Broad immune-inflamed state (supporting context)",
  "",
  "IFN-alpha response, inflammatory response, and complement are also strict common pathways. They support a shared immune-inflamed state, but are broader and less useful as a single mechanistic route.",
  "",
  "## Important exclusions and interpretation limits",
  "",
  "- Do not describe the result as a replicated IL-12 mechanism; the exact IL-12 sets fail the two-cohort significance criterion.",
  "- Do not describe CD28 as the shared cross-cohort route. CD28 was a leading TIGER candidate, but the ccRCC median analysis did not support CD28-specific enrichment.",
  "- Bulk RNA cannot resolve APC, T-cell, and tumor-cell sources, and shared immune pathways may partly reflect immune-cell abundance.",
  "- Cross-sectional associations do not establish causality, ligand secretion, receptor engagement, or phospho-STAT/NF-kappaB activation.",
  "",
  "## Reader-facing outputs",
  "",
  "- `figures/Fig_TIGER_Melanoma_ccRCC_Common_Pathway_Priorities.png`",
  "- `tables/02_strict_aligned_both_BH_significant_pathways.csv`",
  "- `tables/03_IL12_exact_pathway_comparison.csv`",
  "- `tables/04_reader_priority_nonredundant_pathways.csv`",
  "- `tables/05_priority_shared_leading_edge_genes.csv`",
  "",
  "## Reproducibility",
  "",
  "Input paths and MD5 hashes are recorded in `tables/07_input_manifest_md5.csv`. Session information is in `sessionInfo.txt`."
)

writeLines(report_lines, file.path(out_dir, "REPORT.md"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))

cat("Output directory:", out_dir, "\n")
cat("Exact common sets:", nrow(common), "\n")
cat("Strict target-direction common sets:", nrow(strict_common), "\n")
cat("IL-12 exact sets:", nrow(il12), "; strict common:", sum(il12$strict_target_common), "\n")
