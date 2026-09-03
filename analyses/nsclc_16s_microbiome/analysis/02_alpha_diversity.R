## 02_alpha_diversity.R
## Within-sample (alpha) diversity: Health vs NSCLC.
## Test: Wilcoxon rank-sum (two independent groups, non-parametric; HGMT paper).
## NOTE on data type: HGMT provides genus relative abundance (TSS), NOT counts.
##  - Shannon/Simpson/InvSimpson are valid on proportions (vegan normalizes by row sums).
##  - Observed richness = #genera with abundance>0; it is detection-based (depends on HGMT
##    profiling, not raw read depth) so it is reported with that caveat.

source("analysis/00_setup.R")
log_con <- file(file.path(LOG_DIR, "02_alpha.log"), open = "wt"); sink(log_con, split = TRUE)
cat("==== 02_alpha_diversity ====", as.character(Sys.time()), "\n")

prep <- readRDS(file.path(DERIVED_DIR, "prepared_16S.rds"))
meta <- prep$meta; tss <- prep$genus_tss
stopifnot(all(meta$run_id == rownames(tss)))           # ID alignment guard

## ---- compute alpha metrics on genus relative abundance ----
alpha <- data.table(
  run_id     = rownames(tss),
  group      = meta$group,
  Shannon    = diversity(tss, index = "shannon"),
  Simpson    = diversity(tss, index = "simpson"),
  InvSimpson = diversity(tss, index = "invsimpson"),
  Observed   = rowSums(tss > 0)
)
fwrite(alpha, file.path(RESULTS_DIR, "alpha_diversity_per_sample.csv"))

## ---- Wilcoxon rank-sum per metric (Health vs NSCLC) ----
metrics <- c("Shannon", "Simpson", "InvSimpson", "Observed")
res <- rbindlist(lapply(metrics, function(m) {
  x <- alpha[[m]]; g <- alpha$group
  w <- wilcox.test(x ~ g, exact = FALSE)
  data.table(metric = m,
             median_Health = median(x[g == "Health"]),
             median_NSCLC  = median(x[g == "NSCLC"]),
             W = unname(w$statistic), p_value = w$p.value)
}))
res[, p_adj_BH := p.adjust(p_value, "BH")]
fwrite(res, file.path(RESULTS_DIR, "alpha_diversity_tests.csv"))
cat("\n[Wilcoxon rank-sum: Health vs NSCLC]\n"); print(res)

## ---- sensitivity: covariate-adjusted Shannon (age + sex + BMI) ----
df <- merge(alpha, meta[, .(run_id, age, sex, bmi)], by = "run_id")
lm_sh <- lm(Shannon ~ group + age + sex + bmi, data = df)
co <- summary(lm_sh)$coefficients
adj <- data.table(term = rownames(co), estimate = co[, 1],
                  std_error = co[, 2], t = co[, 3], p_value = co[, 4])
fwrite(adj, file.path(RESULTS_DIR, "alpha_shannon_adjusted_lm.csv"))
cat("\n[Shannon ~ group + age + sex + BMI] group term:\n")
print(adj[grepl("group", term)])

## ---- plots ----
alpha_long <- melt(alpha, id.vars = c("run_id", "group"),
                   measure.vars = metrics, variable.name = "metric", value.name = "value")
p <- ggplot(alpha_long, aes(group, value, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.6, alpha = 0.4) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(x = NULL, y = "Alpha diversity", title = "PRJEB26531 16S — alpha diversity (Health vs NSCLC)") +
  theme_hvc + theme(legend.position = "none")
save_png(p, file.path(FIG_DIR, "alpha_all_metrics.png"), 10, 3.2)

pS <- ggplot(alpha, aes(group, Shannon, fill = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.45) +
  scale_fill_manual(values = GROUP_COLORS) +
  labs(x = NULL, y = "Shannon index",
       subtitle = sprintf("Wilcoxon P = %.3g", res[metric == "Shannon", p_value])) +
  theme_hvc + theme(legend.position = "none")
save_png(pS, file.path(FIG_DIR, "alpha_shannon_boxplot.png"), 3.6, 3.6)
ggsave(file.path(FIG_DIR, "alpha_shannon_boxplot.pdf"), pS, width = 3.6, height = 3.6)

cat("\n==== done ====\n"); sink(); close(log_con)
