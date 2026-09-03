#!/usr/bin/env Rscript

# Descriptive preflight for the exact GC-MS Sarcosine feature.
# This script consumes outputs from 002_matched_multiomics_preflight.py and
# does not modify the public source matrices.

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
if (length(script_arg) != 1L) stop("Could not resolve script path")
script_path <- normalizePath(sub("^--file=", "", script_arg))
project_dir <- dirname(dirname(script_path))
out_dir <- file.path(project_dir, "results", "00_input_audit_26.09.02")

paired_path <- file.path(out_dir, "paired_tumor_normal_sarcosine.csv")
group_path <- file.path(out_dir, "tumor_sarcosine_group_map.csv")
stopifnot(file.exists(paired_path), file.exists(group_path))

paired <- read.csv(paired_path, check.names = FALSE)
groups <- read.csv(group_path, check.names = FALSE)
stopifnot(nrow(paired) == 50L, nrow(groups) == 100L)
stopifnot(all(is.finite(paired$log2_tumor_normal_ratio)))

wilcox_result <- wilcox.test(
  paired$log2_tumor_normal_ratio,
  mu = 0,
  alternative = "two.sided",
  exact = FALSE,
  conf.int = TRUE
)

nonzero <- paired$log2_tumor_normal_ratio[paired$log2_tumor_normal_ratio != 0]
ranks <- rank(abs(nonzero), ties.method = "average")
w_pos <- sum(ranks[nonzero > 0])
w_neg <- sum(ranks[nonzero < 0])
paired_rank_biserial <- (w_pos - w_neg) / (w_pos + w_neg)

stats_table <- data.frame(
  analysis = "Paired tumor vs NAT log2 normalized Sarcosine intensity",
  n_pairs = nrow(paired),
  median_log2_tumor_normal_ratio = median(paired$log2_tumor_normal_ratio),
  hodges_lehmann_pseudomedian = unname(wilcox_result$estimate),
  confidence_low = wilcox_result$conf.int[1],
  confidence_high = wilcox_result$conf.int[2],
  wilcoxon_W = unname(wilcox_result$statistic),
  two_sided_p = wilcox_result$p.value,
  paired_rank_biserial = paired_rank_biserial,
  tumor_higher_pairs = sum(nonzero > 0),
  tumor_lower_pairs = sum(nonzero < 0),
  stringsAsFactors = FALSE
)
write.csv(
  stats_table,
  file.path(out_dir, "paired_tumor_normal_sarcosine_statistics.csv"),
  row.names = FALSE,
  quote = TRUE
)

paired_long <- rbind(
  data.frame(patient_id = paired$patient_id, tissue = "NAT", intensity = paired$normal_sarcosine),
  data.frame(patient_id = paired$patient_id, tissue = "Tumor", intensity = paired$tumor_sarcosine)
)
paired_long$tissue <- factor(paired_long$tissue, levels = c("NAT", "Tumor"))

colors <- c("NAT" = "#1F9E89", "Tumor" = "#C83E3E", "Low" = "#1F9E89", "High" = "#C83E3E")

p1 <- ggplot(paired_long, aes(tissue, log2(intensity), group = patient_id)) +
  geom_line(color = "grey65", linewidth = 0.35, alpha = 0.6) +
  geom_point(aes(color = tissue), size = 1.8, alpha = 0.85) +
  scale_color_manual(values = colors[c("NAT", "Tumor")]) +
  labs(
    title = "Paired tissue Sarcosine",
    subtitle = sprintf("n = 50 pairs; Wilcoxon P = %.3g", wilcox_result$p.value),
    x = NULL,
    y = expression(log[2]~"normalized GC-MS intensity")
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

groups$sarcosine_group <- factor(groups$sarcosine_group, levels = c("Low", "High"))
p2 <- ggplot(groups, aes(sarcosine_group, log2_sarcosine_intensity, color = sarcosine_group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.55) +
  geom_jitter(width = 0.14, height = 0, size = 1.5, alpha = 0.65) +
  scale_color_manual(values = colors[c("Low", "High")]) +
  labs(
    title = "Tumor median split",
    subtitle = sprintf("cutpoint = %.4f; n = 50 per group", unique(groups$tumor_median_cutpoint)),
    x = "Tissue Sarcosine group",
    y = expression(log[2]~"normalized GC-MS intensity")
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none", plot.title = element_text(face = "bold"))

png_path <- file.path(out_dir, "Fig_Sarcosine_Matched_Data_Preflight.png")
png(png_path, width = 3300, height = 1500, res = 300)
gridExtra::grid.arrange(p1, p2, ncol = 2, top = grid::textGrob(
  "TJ-RCC matched tissue Sarcosine preflight",
  gp = grid::gpar(fontsize = 17, fontface = "bold")
))
dev.off()

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_R.txt"))
print(stats_table)
cat("Figure:", png_path, "\n")
