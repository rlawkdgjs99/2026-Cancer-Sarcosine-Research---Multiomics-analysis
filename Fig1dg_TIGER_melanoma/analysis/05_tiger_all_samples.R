# =============================================================================
# 05_tiger_all_samples.R
# Reproduce the TIGER-website-style analysis: ALL samples (NO PRE-only filter,
# NO patient deduplication), naive tests — this is what TIGER's portal shows.
#
# INTEGRITY GUARDRAILS:
#   - These pool pre- AND on-treatment (EDT) biopsies and count repeated patients
#     more than once (pseudoreplication). Every figure is labelled
#     "all samples (incl. on-treatment & repeated biopsies)" so it is NOT confused
#     with the rigorous PRE-only analysis (figures/.../primary_PRE_only).
#   - A companion table reports the repeated-measures-CORRECTED p (mixed model for
#     R vs NR; cluster-robust Cox for survival) next to the naive p, so the reader
#     sees whether significance survives proper handling of the repeats.
#
# Run with cwd = Melanoma-PRJEB23709/ :  Rscript analysis/05_tiger_all_samples.R
# =============================================================================
set.seed(42)
source("analysis/00_setup.R")
suppressPackageStartupMessages({ library(cowplot); library(lmerTest) })

d <- readRDS("results/analysis_data_all91.rds")     # all 91 samples; genes = log2(FPKM+1)
stopifnot(nrow(d) == 91L)

# Module z-scores recomputed WITHIN the all-91 cohort (parallel to PRE-only build)
z <- function(x) (x - mean(x)) / sd(x)
d$Degradation_score <- rowMeans(sapply(c("SARDH","PIPOX"), function(g) z(d[[g]])))
d$Production_score  <- rowMeans(sapply(c("GNMT","DMGDH"),  function(g) z(d[[g]])))

vars <- list(
  list(key="SARDH",             label="SARDH",             ylab="log2(FPKM + 1)"),
  list(key="PIPOX",             label="PIPOX",             ylab="log2(FPKM + 1)"),
  list(key="GNMT",              label="GNMT",              ylab="log2(FPKM + 1)"),
  list(key="DMGDH",             label="DMGDH",             ylab="log2(FPKM + 1)"),
  list(key="Degradation_score", label="Degradation score", ylab="Module score (z)"),
  list(key="Production_score",  label="Production score",  ylab="Module score (z)")
)
# 2-line ASCII caption (em/en dashes break under cairo; keep it plain & wrapped)
CAP <- "All samples: pre + on-treatment biopsies, repeats\nincluded (NOT baseline-only)"

# ---- builders (label all-samples basis; annotate NAIVE p, matching TIGER) ----
make_box <- function(dd, var, ylab, title, pval, nrep) {
  ns <- table(dd$response_group)
  xlabs <- c(R=sprintf("R\n(n=%d)",ns[["R"]]), NR=sprintf("NR\n(n=%d)",ns[["NR"]]))
  cap <- sprintf("All samples: pre + on-treatment, repeats included\nn=%d (%d repeat rows); NOT baseline-only", nrow(dd), nrep)
  ggplot(dd, aes(x=response_group, y=.data[[var]], fill=response_group)) +
    geom_boxplot(width=0.6, outlier.shape=NA, alpha=0.85, linewidth=0.5) +
    geom_jitter(width=0.15, height=0, size=1.8, alpha=0.55, color="grey20") +
    scale_fill_manual(values=RESP_COLORS) + scale_x_discrete(labels=xlabs) +
    labs(title=title, subtitle=sprintf("Wilcoxon %s", fmt_p(pval)),
         x=NULL, y=ylab, caption=cap) +
    theme_pub() + theme(plot.caption=element_text(size=10.5, hjust=0, color="grey35"))
}
km_one <- function(dd, var, title, nrep) {
  cut <- median(dd[[var]])
  dd$expr_group <- factor(ifelse(dd[[var]]>cut,"High","Low"), levels=c("Low","High"))
  sd <- survival::survdiff(Surv(os_days,event)~expr_group, data=dd)
  lr <- pchisq(sd$chisq, length(sd$n)-1, lower.tail=FALSE)
  cox <- survival::coxph(Surv(os_days,event)~expr_group, data=dd)
  cf <- summary(cox)$coefficients; ci <- summary(cox)$conf.int
  hr <- cf["expr_groupHigh","exp(coef)"]; lo<-ci["expr_groupHigh","lower .95"]; hi<-ci["expr_groupHigh","upper .95"]
  fit <- survival::survfit(Surv(os_days,event)~expr_group, data=dd)
  ptxt <- sprintf("Log-rank %s\nHR = %.2f (%.2f-%.2f)\nall samples n=%d (%d repeats)",
                  fmt_p(lr), hr, lo, hi, nrow(dd), nrep)
  gg <- survminer::ggsurvplot(fit, data=dd, palette=c(COL_GREEN,COL_RED),
    legend.title="Expression", legend.labs=c("Low","High"), legend="top",
    pval=ptxt, pval.size=4.4, pval.coord=c(0.02*max(dd$os_days),0.14),
    risk.table=TRUE, risk.table.height=0.27, risk.table.fontsize=4.6,
    risk.table.title="Number at risk", tables.y.text=FALSE, conf.int=FALSE,
    censor.shape="|", censor.size=4, xlab="Time (days)", ylab="Overall survival", title=title,
    font.main=c(22,"bold"), font.x=17, font.y=17, font.tickslab=15, font.legend=15,
    ggtheme=theme_classic(base_size=16))
  gg$plot <- gg$plot + theme(plot.title=element_text(hjust=0.5, size=22, face="bold")) +
    labs(caption=CAP) + theme(plot.caption=element_text(size=10.5, hjust=0, color="grey35"))
  list(lr=lr, hr=hr, lo=lo, hi=hi, gg=gg,
       n_low=sum(dd$expr_group=="Low"), n_high=sum(dd$expr_group=="High"),
       ev_low=sum(dd$event[dd$expr_group=="Low"]), ev_high=sum(dd$event[dd$expr_group=="High"]))
}
nrep_of <- function(dd) sum(duplicated(dd$patient_name))   # how many rows are repeat-patient

# ============================ R vs NR (all samples) =========================
resp_cohorts <- list(
  list(key="pooled",  label="All (91)", sub=d),
  list(key="antiPD1", label="anti-PD-1", sub=d[d$therapy_short=="antiPD1",]),
  list(key="combo",   label="Combo",     sub=d[d$therapy_short=="combo",])
)
resp_rows <- list(); resp_plots <- list()
for (co in resp_cohorts) for (v in vars) {
  dd <- co$sub; nrep <- nrep_of(dd)
  pw <- suppressWarnings(wilcox.test(dd[[v$key]] ~ dd$response_group)$p.value)   # naive (TIGER)
  # repeated-measures corrected: linear mixed model with random patient effect
  pm <- tryCatch({
    m <- suppressWarnings(lmerTest::lmer(reformulate(c("response_group","(1|patient_name)"), v$key), data=dd))
    summary(m)$coefficients["response_groupNR","Pr(>|t|)"]
  }, error=function(e) NA_real_)
  med <- tapply(dd[[v$key]], dd$response_group, median)
  resp_rows[[length(resp_rows)+1]] <- data.frame(
    feature=v$label, cohort=co$label, n=nrow(dd), n_repeat_rows=nrep,
    n_R=sum(dd$response_group=="R"), n_NR=sum(dd$response_group=="NR"),
    median_R=round(med[["R"]],4), median_NR=round(med[["NR"]],4),
    p_naive_Wilcox=signif(pw,4), p_mixedmodel_corrected=signif(pm,4), stringsAsFactors=FALSE)
  resp_plots[[paste(v$key,co$key,sep="__")]] <- make_box(dd, v$key, v$ylab, v$label, pw, nrep)
}
resp <- do.call(rbind, resp_rows)
write.csv(resp, "results/tables/TIGER_response_RvsNR_stats.csv", row.names=FALSE)
cat("=== TIGER-style R vs NR (all samples): naive vs mixed-model-corrected p ===\n")
print(resp[,c("feature","cohort","n","p_naive_Wilcox","p_mixedmodel_corrected")], row.names=FALSE)

for (nm in names(resp_plots)) save_plot(resp_plots[[nm]], paste0("tiger_box_", nm), 5.3, 6.0)
gk <- c("SARDH","PIPOX","GNMT","DMGDH"); mk <- c("Degradation_score","Production_score")
save_plot(plot_grid(plotlist=lapply(gk,function(k) resp_plots[[paste(k,"pooled",sep="__")]]),
          ncol=2, labels="AUTO", label_size=18), "tiger_panel_box_genes_pooled", 10.6, 12)
save_plot(plot_grid(plotlist=lapply(mk,function(k) resp_plots[[paste(k,"pooled",sep="__")]]),
          ncol=2, labels="AUTO", label_size=18), "tiger_panel_box_modules_pooled", 10.6, 6.4)

# ============================ Survival (all samples) ========================
surv_cohorts <- list(
  list(key="pooled",  label="All (91)", sub=d),
  list(key="antiPD1", label="anti-PD-1", sub=d[d$therapy_short=="antiPD1",])
)
surv_rows <- list(); surv_plots <- list()
for (co in surv_cohorts) for (v in vars) {
  dd <- co$sub; nrep <- nrep_of(dd)
  out <- km_one(dd, v$key, v$label, nrep)
  # cluster-robust Cox (same data, corrects within-patient correlation)
  cut <- median(dd[[v$key]]); dd$g <- factor(ifelse(dd[[v$key]]>cut,"High","Low"),levels=c("Low","High"))
  pcr <- tryCatch(summary(survival::coxph(Surv(os_days,event)~g+cluster(patient_name),data=dd))$coefficients["gHigh","Pr(>|z|)"],
                  error=function(e) NA_real_)
  surv_rows[[length(surv_rows)+1]] <- data.frame(
    feature=v$label, cohort=co$label, n=nrow(dd), n_repeat_rows=nrep,
    n_low=out$n_low, n_high=out$n_high, events_low=out$ev_low, events_high=out$ev_high,
    p_naive_logrank=signif(out$lr,4), HR_high_vs_low=round(out$hr,3),
    HR_lower=round(out$lo,3), HR_upper=round(out$hi,3),
    p_clusterrobust_corrected=signif(pcr,4), stringsAsFactors=FALSE)
  full <- plot_grid(out$gg$plot, out$gg$table, ncol=1, rel_heights=c(3.2,1.15), align="v", axis="lr")
  save_plot(full, paste0("tiger_km_", v$key, "__", co$key), 6, 7.1)
  surv_plots[[paste(v$key,co$key,sep="__")]] <- out$gg$plot
}
surv <- do.call(rbind, surv_rows)
write.csv(surv, "results/tables/TIGER_survival_KM_stats.csv", row.names=FALSE)
cat("\n=== TIGER-style survival (all samples): naive vs cluster-robust-corrected p ===\n")
print(surv[,c("feature","cohort","n","p_naive_logrank","HR_high_vs_low","p_clusterrobust_corrected")], row.names=FALSE)

for (co in surv_cohorts) {
  save_plot(plot_grid(plotlist=lapply(gk,function(k) surv_plots[[paste(k,co$key,sep="__")]]),
            ncol=2, labels="AUTO", label_size=18), paste0("tiger_panel_km_genes_", co$key), 11, 10.4)
  save_plot(plot_grid(plotlist=lapply(mk,function(k) surv_plots[[paste(k,co$key,sep="__")]]),
            ncol=2, labels="AUTO", label_size=18), paste0("tiger_panel_km_modules_", co$key), 11, 5.8)
}

cat("\nSaved TIGER-style figures (tiger_*) -> figures/*/TIGER_all_samples/\n")
cat("Saved tables: TIGER_response_RvsNR_stats.csv, TIGER_survival_KM_stats.csv (naive + corrected p)\n")
