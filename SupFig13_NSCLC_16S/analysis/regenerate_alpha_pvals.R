#!/usr/bin/env Rscript
# ============================================================================
# Re-render the 16S 4-metric alpha-diversity boxplot WITH per-facet Wilcoxon P
# labels. READ-ONLY: reads the already-computed alpha_diversity_per_sample.csv
# and alpha_diversity_tests.csv (from 02_alpha_diversity.R); NO statistic is
# recomputed. Mirrors the original plot (00_setup.R aesthetics) and only adds
# the P labels. Path-robust (resolves results/figures relative to this script,
# so it works on the backup where 00_setup.R's BASE_DIR is stale).
# Output: analysis/figures/alpha_all_metrics_pval.png  (original preserved)
# ============================================================================
options(bitmapType = "quartz")
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

## Run with the R working directory set to the analysis-folder root; the script's
## inputs/outputs live under analysis/ (relative to that root).
adir <- "analysis"
RES <- file.path(adir, "results"); FIG <- file.path(adir, "figures")

alpha <- fread(file.path(RES, "alpha_diversity_per_sample.csv"))
tests <- fread(file.path(RES, "alpha_diversity_tests.csv"))

metrics      <- c("Shannon", "Simpson", "InvSimpson", "Observed")
GROUP_COLORS <- c(Health = "#1B7837", NSCLC = "#B2182B")   # verbatim from 00_setup.R
theme_hvc    <- theme_bw(base_size = 16) + theme(panel.grid.minor = element_blank())

alpha$group <- factor(alpha$group, levels = c("Health", "NSCLC"))
long <- melt(alpha, id.vars = c("run_id", "group"), measure.vars = metrics,
             variable.name = "metric", value.name = "value")
long$metric <- factor(long$metric, levels = metrics)

fmtp <- function(p) ifelse(p < 0.0001, "P < 0.0001", sprintf("P = %.3g", p))
labs_df <- data.frame(metric = factor(tests$metric, levels = metrics),
                      label  = fmtp(tests$p_value), stringsAsFactors = FALSE)

set.seed(42)   # reproducible jitter
p <- ggplot(long, aes(group, value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = GROUP_COLORS) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.13))) +
  geom_text(data = labs_df, aes(x = 1.5, label = label), y = Inf, vjust = 1.5,
            size = 4.0, fontface = "bold", inherit.aes = FALSE) +
  labs(x = NULL, y = "Alpha diversity",
       title = "PRJEB26531 16S — alpha diversity (Health vs NSCLC)",
       subtitle = "Wilcoxon rank-sum, Health vs NSCLC") +
  theme_hvc + theme(legend.position = "none")

out <- file.path(FIG, "alpha_all_metrics_pval.png")
grDevices::png(out, width = 10, height = 3.4, units = "in", res = 200, type = "quartz")
print(p); grDevices::dev.off()
cat("wrote", out, "\n"); print(labs_df)
