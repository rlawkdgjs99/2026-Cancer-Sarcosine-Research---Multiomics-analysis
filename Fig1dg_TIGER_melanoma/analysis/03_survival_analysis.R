# =============================================================================
# 03_survival_analysis.R  --  Overall survival by expression High vs Low
#
# STATISTICAL DESIGN
#   Endpoint  : Overall survival (days); event = Dead(1)/Alive(0). (No PFS available.)
#   Grouping  : per feature, median split WITHIN the analyzed cohort.
#               High = value > cohort median ; Low = value <= cohort median.
#   Cohorts   : pooled PRE (n=73, 29 events) and anti-PD-1 mono (n=41, 24 events).
#               Combo (n=32) is EXCLUDED from survival: only 5 events -> KM/Cox
#               would be statistically unreliable (user-approved decision).
#   Tests     : log-rank (survdiff) for KM difference; univariate Cox PH for
#               HR(High vs Low) with 95% CI (Low = reference).
#   Colors    : Low = GREEN, High = RED (user preference).
#   The log-rank P and HR shown on each KM plot are the SAME values in the table.
# =============================================================================
set.seed(42)
source("analysis/00_setup.R")
suppressPackageStartupMessages({ library(cowplot) })

dat <- readRDS("results/analysis_data_PRE.rds")
stopifnot(nrow(dat) == 73L)

vars <- list(
  list(key = "SARDH",             label = "SARDH"),
  list(key = "PIPOX",             label = "PIPOX"),
  list(key = "GNMT",              label = "GNMT"),
  list(key = "DMGDH",             label = "DMGDH"),
  list(key = "Degradation_score", label = "Degradation score"),
  list(key = "Production_score",  label = "Production score")
)
cohorts <- list(
  list(key = "pooled",  label = "All (PRE)", sub = dat),
  list(key = "antiPD1", label = "anti-PD-1", sub = dat[dat$therapy_short == "antiPD1", ])
)

# ---- One feature x cohort: stats + KM plot ----------------------------------
km_one <- function(d, var, title) {
  cut <- median(d[[var]])
  d$expr_group <- factor(ifelse(d[[var]] > cut, "High", "Low"), levels = c("Low", "High"))
  # Guard: both groups must exist
  stopifnot(all(c("Low", "High") %in% levels(droplevels(d$expr_group))))

  # Use column-name formulas with data = d (required by survminer::ggsurvplot,
  # which re-evaluates the survfit call in the data environment).
  sd  <- survival::survdiff(Surv(os_days, event) ~ expr_group, data = d)
  lr_p <- pchisq(sd$chisq, df = length(sd$n) - 1, lower.tail = FALSE)
  cox <- survival::coxph(Surv(os_days, event) ~ expr_group, data = d)  # Low = reference
  cf  <- summary(cox)$coefficients
  ci  <- summary(cox)$conf.int
  hr  <- cf["expr_groupHigh", "exp(coef)"]
  hr_lo <- ci["expr_groupHigh", "lower .95"]
  hr_hi <- ci["expr_groupHigh", "upper .95"]
  cox_p <- cf["expr_groupHigh", "Pr(>|z|)"]

  ev <- tapply(d$event, d$expr_group, sum)
  nn <- table(d$expr_group)
  stats <- data.frame(
    feature = title, n_low = nn[["Low"]], n_high = nn[["High"]],
    events_low = ev[["Low"]], events_high = ev[["High"]],
    median_cut = round(cut, 4),
    logrank_p = signif(lr_p, 4),
    HR_high_vs_low = round(hr, 3), HR_lower = round(hr_lo, 3), HR_upper = round(hr_hi, 3),
    cox_p = signif(cox_p, 4), stringsAsFactors = FALSE
  )

  fit <- survival::survfit(Surv(os_days, event) ~ expr_group, data = d)
  pval_txt <- sprintf("Log-rank %s\nHR = %.2f (%.2f-%.2f)", fmt_p(lr_p), hr, hr_lo, hr_hi)
  gg <- survminer::ggsurvplot(
    fit, data = d,
    palette = c(COL_GREEN, COL_RED),          # Low, High (matches factor level order)
    legend.title = "Expression", legend.labs = c("Low", "High"), legend = "top",
    pval = pval_txt, pval.size = 5, pval.coord = c(0.02 * max(d$os_days), 0.10),
    risk.table = TRUE, risk.table.height = 0.27, risk.table.fontsize = 4.6,
    risk.table.title = "Number at risk", tables.y.text = FALSE,
    conf.int = FALSE, censor.shape = "|", censor.size = 4,
    xlab = "Time (days)", ylab = "Overall survival", title = title,
    font.main = c(22, "bold"), font.x = c(17), font.y = c(17),
    font.tickslab = c(15), font.legend = c(15),
    ggtheme = theme_classic(base_size = 16)
  )
  gg$plot <- gg$plot + theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"))
  list(stats = stats, gg = gg)
}

# ---- Run --------------------------------------------------------------------
res_rows <- list(); plot_only <- list()
for (co in cohorts) {
  for (v in vars) {
    out <- km_one(co$sub, v$key, v$label)
    out$stats$cohort <- co$label
    res_rows[[length(res_rows) + 1]] <- out$stats
    # full KM (curve + risk table) saved individually
    full <- cowplot::plot_grid(out$gg$plot, out$gg$table, ncol = 1,
                               rel_heights = c(3.2, 1.15), align = "v", axis = "lr")
    save_plot(full, paste0("km_", v$key, "__", co$key), width = 6, height = 6.8)
    plot_only[[paste(v$key, co$key, sep = "__")]] <- out$gg$plot
  }
}
res <- do.call(rbind, res_rows)
res <- res[, c("feature", "cohort", "n_low", "n_high", "events_low", "events_high",
               "median_cut", "logrank_p", "HR_high_vs_low", "HR_lower", "HR_upper", "cox_p")]
write.csv(res, "results/tables/survival_KM_stats.csv", row.names = FALSE)
cat("=== Survival (median-split) results ===\n"); print(res, row.names = FALSE)

# ---- Combined KM panels (curves only) per cohort ----------------------------
gene_keys <- c("SARDH", "PIPOX", "GNMT", "DMGDH")
mod_keys  <- c("Degradation_score", "Production_score")
for (co in cohorts) {
  pg <- cowplot::plot_grid(plotlist = lapply(gene_keys,
        function(k) plot_only[[paste(k, co$key, sep = "__")]]),
        ncol = 2, labels = "AUTO", label_size = 18)
  save_plot(pg, paste0("panel_km_genes_", co$key), width = 11, height = 10)
  pm <- cowplot::plot_grid(plotlist = lapply(mod_keys,
        function(k) plot_only[[paste(k, co$key, sep = "__")]]),
        ncol = 2, labels = "AUTO", label_size = 18)
  save_plot(pm, paste0("panel_km_modules_", co$key), width = 11, height = 5.4)
}

cat("\nSaved KM plots + panels to results/figures/\n")
cat("Saved stats table: results/tables/survival_KM_stats.csv\n")
