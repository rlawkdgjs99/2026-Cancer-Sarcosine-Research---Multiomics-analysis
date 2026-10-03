## 06_query_genera.R
## Look up specific genera (all GTDB variants) in the Healthy-vs-NSCLC DA results.
## Distinguishes: tested & significant / tested & ns / NOT tested (filtered by <10% prevalence).

source("analysis/00_setup.R")   # relative to the analysis-folder root (working dir)

prep <- readRDS(file.path(DERIVED_DIR, "prepared_16S.rds"))
raw  <- prep$genus_raw; look <- as.data.table(prep$feature_lookup); meta <- prep$meta
stopifnot(all(meta$run_id == rownames(raw)))
da  <- fread(file.path(RESULTS_DIR, "DA_maaslin2_TSS_LOG.csv"))   # primary
dac <- fread(file.path(RESULTS_DIR, "DA_maaslin2_CLR.csv"))       # sensitivity

targets <- c("Alistipes", "Lachnospira", "Roseburia", "Hungatella")
hits <- look[grepl(paste(targets, collapse = "|"), genus, ignore.case = TRUE)]
cat("matched GTDB variants:", nrow(hits), "\n")

g <- meta$group
hits[, `:=`(
  prevalence_pct = round(colMeans(raw > 0)[feature] * 100, 1),
  present_Health_pct = round(colMeans(raw[g == "Health", , drop = FALSE] > 0)[feature] * 100, 1),
  present_NSCLC_pct  = round(colMeans(raw[g == "NSCLC",  , drop = FALSE] > 0)[feature] * 100, 1),
  mean_Health_pct = round(colMeans(raw[g == "Health", , drop = FALSE])[feature], 3),
  mean_NSCLC_pct  = round(colMeans(raw[g == "NSCLC",  , drop = FALSE])[feature], 3)
)]
hits <- merge(hits, da[, .(feature, coef_tss = coef, p_tss = pval, q_tss = qval)], by = "feature", all.x = TRUE)
hits <- merge(hits, dac[, .(feature, q_clr = qval)], by = "feature", all.x = TRUE)
hits[, tested_in_DA := !is.na(q_tss)]                        # appears in MaAsLin2 results?
hits[, direction := fifelse(is.na(coef_tss), "NA",
                     fifelse(coef_tss > 0, "NSCLC-higher", "Health-higher"))]
hits[, sig_q05 := fifelse(is.na(q_tss), NA, q_tss < 0.05)]

setorder(hits, genus)
print(hits[, .(genus, prevalence_pct, present_Health_pct, present_NSCLC_pct,
               mean_Health_pct, mean_NSCLC_pct, tested_in_DA, direction,
               coef_tss = round(coef_tss, 3), q_tss = signif(q_tss, 3),
               q_clr = signif(q_clr, 3), sig_q05)])
fwrite(hits, file.path(RESULTS_DIR, "query_genera_lookup.csv"))
cat("\nNote: prevalence filter for DA = >=", PREVALENCE_MIN * 100, "% of samples.\n")
