## 04_differential_abundance.R
## Genus-level differential abundance: NSCLC vs Health (reference = Health).
## Method: MaAsLin2 (general linear model on relative-abundance microbiome data).
##   - Primary    : normalization=TSS, transform=LOG
##   - Sensitivity: normalization=CLR, transform=NONE (compositional)
## WHY NOT ANCOM-BC2 (the HGMT paper's tool): ANCOM-BC2 requires COUNT data to estimate
##   sampling fractions; HGMT provides only relative abundance, so feeding proportions would
##   violate its model. MaAsLin2 is designed for relative-abundance input. (ANCOM-BC2 can be
##   added as a caveated sensitivity if count data are ever obtained.)
## Confounders: per HGMT, test age/sex/BMI vs group; include those with P<0.05 as covariates.
## Prevalence filter: genus must be present (>0) in >=10% of samples (MaAsLin2 min_prevalence).

source("analysis/00_setup.R")
suppressPackageStartupMessages(library(Maaslin2))
log_con <- file(file.path(LOG_DIR, "04_DA.log"), open = "wt"); sink(log_con, split = TRUE)
cat("==== 04_differential_abundance ====", as.character(Sys.time()), "\n")

prep <- readRDS(file.path(DERIVED_DIR, "prepared_16S.rds"))
meta <- prep$meta; raw <- prep$genus_raw; lookup <- prep$feature_lookup
stopifnot(all(meta$run_id == rownames(raw)))

## ---- confounder testing (HGMT approach) ----
sex_p <- fisher.test(table(meta$sex, meta$group))$p.value
age_p <- wilcox.test(age ~ group, data = meta, exact = FALSE)$p.value
bmi_p <- wilcox.test(bmi ~ group, data = meta, exact = FALSE)$p.value
conf <- data.table(variable = c("sex", "age", "bmi"),
                   test = c("Fisher", "Wilcoxon", "Wilcoxon"),
                   p_value = c(sex_p, age_p, bmi_p))
conf[, is_confounder := p_value < ALPHA_SIG]
fwrite(conf, file.path(RESULTS_DIR, "confounder_tests.csv"))
cat("\n[confounder tests vs group]\n"); print(conf)
covars <- conf[is_confounder == TRUE, variable]
fixed_effects <- c(covars, "group")
cat("fixed effects for DA:", paste(fixed_effects, collapse = " + "), "\n")

## ---- metadata + feature inputs for MaAsLin2 (samples as rows) ----
md <- as.data.frame(meta[, c("age", "sex", "bmi", "group"), with = FALSE])
rownames(md) <- meta$run_id
feat <- as.data.frame(raw); rownames(feat) <- rownames(raw)   # samples x genera (raw % )

run_maaslin <- function(norm, trans, outdir) {
  Maaslin2(input_data = feat, input_metadata = md, output = outdir,
           normalization = norm, transform = trans, analysis_method = "LM",
           correction = "BH", min_prevalence = PREVALENCE_MIN, min_abundance = 0,
           fixed_effects = fixed_effects, reference = "group,Health",
           plot_heatmap = FALSE, plot_scatter = FALSE, max_significance = 0.25, cores = 1)
  invisible(NULL)
}
run_maaslin("TSS", "LOG", file.path(RESULTS_DIR, "maaslin2_TSS_LOG"))
run_maaslin("CLR", "NONE", file.path(RESULTS_DIR, "maaslin2_CLR"))

## ---- collect the 'group' (NSCLC vs Health) effect from each model ----
read_group_da <- function(dir, label) {
  r <- fread(file.path(dir, "all_results.tsv"))
  r <- r[metadata == "group"]                # the group term; value == "NSCLC"
  r[, method := label]
  merge(r, lookup, by = "feature", all.x = TRUE)
}
da_tss <- read_group_da(file.path(RESULTS_DIR, "maaslin2_TSS_LOG"), "TSS_LOG")
da_clr <- read_group_da(file.path(RESULTS_DIR, "maaslin2_CLR"), "CLR")
fwrite(da_tss, file.path(RESULTS_DIR, "DA_maaslin2_TSS_LOG.csv"))
fwrite(da_clr, file.path(RESULTS_DIR, "DA_maaslin2_CLR.csv"))

cat("\n[primary TSS_LOG] features tested:", nrow(da_tss),
    "| q<0.05:", sum(da_tss$qval < 0.05), "| q<0.25:", sum(da_tss$qval < 0.25), "\n")
cat("[sens.  CLR]     features tested:", nrow(da_clr),
    "| q<0.05:", sum(da_clr$qval < 0.05), "| q<0.25:", sum(da_clr$qval < 0.25), "\n")

cat("\n[Top 15 by qval — primary TSS_LOG] (coef>0 => higher in NSCLC)\n")
print(head(da_tss[order(qval), .(genus, coef, pval, qval)], 15))

## ---- concordance between the two models (q<0.05, same direction) ----
sig_tss <- da_tss[qval < 0.05, .(feature, genus, coef_tss = coef, q_tss = qval)]
sig_clr <- da_clr[qval < 0.05, .(feature, coef_clr = coef, q_clr = qval)]
conc <- merge(sig_tss, sig_clr, by = "feature")
conc[, same_direction := sign(coef_tss) == sign(coef_clr)]
fwrite(conc, file.path(RESULTS_DIR, "DA_concordance_TSSLOG_vs_CLR.csv"))
cat("\n[concordance] q<0.05 in BOTH:", nrow(conc),
    "| same direction:", sum(conc$same_direction), "\n")

## ---- volcano (primary) ----
da_tss[, sig := ifelse(qval < 0.05, "q<0.05", ifelse(qval < 0.25, "q<0.25", "ns"))]
v <- ggplot(da_tss, aes(coef, -log10(qval), color = sig)) +
  geom_point(size = 1.6, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey50") +
  geom_vline(xintercept = 0, linetype = 3, color = "grey60") +
  scale_color_manual(values = c("q<0.05" = "#D7191C", "q<0.25" = "#FDAE61", "ns" = "grey75")) +
  labs(x = "MaAsLin2 coefficient (NSCLC vs Health; >0 = higher in NSCLC)",
       y = "-log10(q)", color = NULL,
       title = "Genus DA — NSCLC vs Health (TSS+LOG, adjusted)") +
  theme_hvc
save_png(v, file.path(FIG_DIR, "DA_volcano_TSS_LOG.png"), 6, 4.5)
ggsave(file.path(FIG_DIR, "DA_volcano_TSS_LOG.pdf"), v, width = 6, height = 4.5)

cat("\n==== done ====\n"); sink(); close(log_con)
