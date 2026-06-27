# =============================================================================
# 02_response_analysis.R  --  ICI Responder (R) vs Non-responder (NR)
#
# STATISTICAL DESIGN
#   Unit of observation : patient (PRE/baseline biopsy, one per patient; n=73)
#   Groups              : R (CR+PR, n=40) vs NR (SD+PD, n=33) -- independent, unpaired
#   Variables tested    : 4 genes [log2(FPKM+1)] + 2 module z-scores (6 features)
#   Test                : Wilcoxon rank-sum (Mann-Whitney), two-sided.
#                         Chosen because (a) FPKM is continuous & non-normal,
#                         (b) it matches TIGER's own DE method, (c) robust at this n.
#   Cohorts             : pooled (n=73); anti-PD-1 mono (n=41); combo (n=32)
#                         -> user asked to compare the two therapies separately.
#   Multiple testing    : Benjamini-Hochberg FDR across the 6 features, within cohort.
#   The p-value shown on each plot is the SAME value stored in the stats table.
# =============================================================================
set.seed(42)
source("analysis/00_setup.R")
suppressPackageStartupMessages({ library(cowplot) })

dat <- readRDS("results/analysis_data_PRE.rds")
stopifnot(nrow(dat) == 73L)

# ---- Feature definitions ----------------------------------------------------
vars <- list(
  list(key = "SARDH",             label = "SARDH",             ylab = "log2(FPKM + 1)", grp = "Degradation"),
  list(key = "PIPOX",             label = "PIPOX",             ylab = "log2(FPKM + 1)", grp = "Degradation"),
  list(key = "GNMT",              label = "GNMT",              ylab = "log2(FPKM + 1)", grp = "Production"),
  list(key = "DMGDH",             label = "DMGDH",             ylab = "log2(FPKM + 1)", grp = "Production"),
  list(key = "Degradation_score", label = "Degradation score", ylab = "Module score (z)", grp = "Module"),
  list(key = "Production_score",  label = "Production score",  ylab = "Module score (z)", grp = "Module")
)
cohorts <- list(
  list(key = "pooled",  label = "All (PRE)",  sub = dat),
  list(key = "antiPD1", label = "anti-PD-1",  sub = dat[dat$therapy_short == "antiPD1", ]),
  list(key = "combo",   label = "Combo",      sub = dat[dat$therapy_short == "combo", ])
)

# ---- Boxplot builder (annotates the EXACT p passed in) ----------------------
make_box <- function(d, var, ylab, title, pval) {
  ns <- table(d$response_group)
  xlabs <- c(R  = sprintf("R\n(n=%d)",  ns[["R"]]),
             NR = sprintf("NR\n(n=%d)", ns[["NR"]]))
  ggplot(d, aes(x = response_group, y = .data[[var]], fill = response_group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, linewidth = 0.5) +
    geom_jitter(width = 0.15, height = 0, size = 1.8, alpha = 0.55, color = "grey20") +
    scale_fill_manual(values = RESP_COLORS) +
    scale_x_discrete(labels = xlabs) +
    labs(title = title, subtitle = paste0("Wilcoxon  ", fmt_p(pval)),
         x = NULL, y = ylab) +
    theme_pub()
}

# ---- Run: stats + plots -----------------------------------------------------
res_rows <- list()
plot_store <- list()   # key: paste(var,cohort)

for (co in cohorts) {
  d <- co$sub
  for (v in vars) {
    x <- d[[v$key]]
    g <- d$response_group
    wt <- wilcox.test(x ~ g)            # two-sided, default (exact if no ties & small n)
    med <- tapply(x, g, median)
    res_rows[[length(res_rows) + 1]] <- data.frame(
      feature   = v$label,
      feature_grp = v$grp,
      cohort    = co$label,
      n_R       = sum(g == "R"),
      n_NR      = sum(g == "NR"),
      median_R  = round(med[["R"]],  4),
      median_NR = round(med[["NR"]], 4),
      W         = unname(wt$statistic),
      p_raw     = wt$p.value,
      stringsAsFactors = FALSE
    )
    ttl <- v$label
    plot_store[[paste(v$key, co$key, sep = "__")]] <-
      make_box(d, v$key, v$ylab, ttl, wt$p.value)
  }
}
res <- do.call(rbind, res_rows)

# BH adjustment across the 6 features, within each cohort
res$p_BH <- NA_real_
for (ck in unique(res$cohort)) {
  idx <- res$cohort == ck
  res$p_BH[idx] <- p.adjust(res$p_raw[idx], method = "BH")
}
res$p_raw <- signif(res$p_raw, 4)
res$p_BH  <- signif(res$p_BH, 4)

write.csv(res, "results/tables/response_RvsNR_stats.csv", row.names = FALSE)
cat("=== R vs NR Wilcoxon results ===\n"); print(res, row.names = FALSE)

# ---- Save individual plots --------------------------------------------------
for (nm in names(plot_store)) {
  save_plot(plot_store[[nm]], paste0("box_", nm), width = 4.6, height = 5.2)
}

# ---- Combined panels --------------------------------------------------------
gene_keys <- c("SARDH", "PIPOX", "GNMT", "DMGDH")
mod_keys  <- c("Degradation_score", "Production_score")

# Pooled: 4-gene 2x2 panel + 2-module panel
panel_genes_pooled <- plot_grid(plotlist = lapply(gene_keys,
  function(k) plot_store[[paste(k, "pooled", sep = "__")]]),
  ncol = 2, labels = "AUTO", label_size = 18)
save_plot(panel_genes_pooled, "panel_box_genes_pooled", width = 9.2, height = 10.4)

panel_mods_pooled <- plot_grid(plotlist = lapply(mod_keys,
  function(k) plot_store[[paste(k, "pooled", sep = "__")]]),
  ncol = 2, labels = "AUTO", label_size = 18)
save_plot(panel_mods_pooled, "panel_box_modules_pooled", width = 9.2, height = 5.4)

# Therapy comparison: per feature, anti-PD-1 | combo side by side
for (v in vars) {
  pp <- plot_grid(
    plot_store[[paste(v$key, "antiPD1", sep = "__")]] + labs(title = paste0(v$label, " — anti-PD-1")),
    plot_store[[paste(v$key, "combo",   sep = "__")]] + labs(title = paste0(v$label, " — Combo")),
    ncol = 2, labels = NULL)
  save_plot(pp, paste0("panel_box_byTherapy_", v$key), width = 9.2, height = 5.4)
}

cat("\nSaved individual boxplots + panels to results/figures/\n")
cat("Saved stats table: results/tables/response_RvsNR_stats.csv\n")
