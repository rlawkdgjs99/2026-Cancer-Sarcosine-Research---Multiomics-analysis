# =============================================================================
# 06_manuscript_figures.R
# Curate a small, paper-ready figure set into results/manuscript_figures/.
#   Main/          = TIGER-style (all-samples) significant DEGRADATION story
#   Supplementary/ = production arm (NS/trend, honest) + PRE-only baseline robustness
#
# Paper-ready = cohort shown in subtitle + the REPEATED-MEASURES-CORRECTED p
#   (mixed-model for R-vs-NR, cluster-robust for survival) read straight from the
#   verified stats tables, so each figure matches its table exactly.
# Originals in results/figures/ are left untouched.
#
# Run with cwd = Melanoma-PRJEB23709/ :  Rscript analysis/06_manuscript_figures.R
# =============================================================================
set.seed(42)
source("analysis/00_setup.R")          # COL_GREEN/RED, RESP_COLORS, theme_pub
suppressPackageStartupMessages({ library(cowplot) })

# ---- data -------------------------------------------------------------------
pre   <- readRDS("results/analysis_data_PRE.rds")     # 73; module scores present
all91 <- readRDS("results/analysis_data_all91.rds")   # 91; genes only
z <- function(x) (x - mean(x)) / sd(x)
all91$Degradation_score <- rowMeans(sapply(c("SARDH","PIPOX"), function(g) z(all91[[g]])))
all91$Production_score  <- rowMeans(sapply(c("GNMT","DMGDH"),  function(g) z(all91[[g]])))

# ---- stats tables (source of truth for the p-values shown) ------------------
tigResp <- read.csv("results/tables/TIGER_response_RvsNR_stats.csv")
tigSurv <- read.csv("results/tables/TIGER_survival_KM_stats.csv")
preResp <- read.csv("results/tables/response_RvsNR_stats.csv")
preSurv <- read.csv("results/tables/survival_KM_stats.csv")
pget <- function(tab, feat, coh, col) {
  v <- tab[tab$feature == feat & tab$cohort == coh, col]
  if (length(v) != 1) stop(sprintf("lookup failed: %s | %s | %s", feat, coh, col))
  v
}

# ---- data subset ------------------------------------------------------------
subset_data <- function(source, cohort_key) {
  d <- if (source == "tiger") all91 else pre
  if (cohort_key == "antiPD1") d <- d[d$therapy_short == "antiPD1", ]
  else if (cohort_key == "combo") d <- d[d$therapy_short == "combo", ]
  d
}

# ---- save helper (manuscript folder; PNG only) ------------------------------
save_ms <- function(plot, sub, stem, w = 5.3, h = 5.9) {
  d <- file.path("results/manuscript_figures", sub)
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  ggsave(file.path(d, paste0(stem, ".png")), plot, width = w, height = h,
         dpi = 300, bg = "white", device = png, type = "cairo")
}
pstar <- function(p) paste0(sprintf("P = %.3f", p), ifelse(p < 0.05, " *", ""))

# ---- renderers --------------------------------------------------------------
box_paper <- function(d, var, feat_label, coh_disp, pval, note) {
  ns <- table(d$response_group)
  xlabs <- c(R = sprintf("R\n(n=%d)", ns[["R"]]), NR = sprintf("NR\n(n=%d)", ns[["NR"]]))
  ylab <- if (grepl("_score$", var)) "Module score (z)" else "log2(FPKM + 1)"
  ggplot(d, aes(response_group, .data[[var]], fill = response_group)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.85, linewidth = 0.5) +
    geom_jitter(width = 0.15, height = 0, size = 1.8, alpha = 0.55, color = "grey20") +
    scale_fill_manual(values = RESP_COLORS) + scale_x_discrete(labels = xlabs) +
    labs(title = feat_label, subtitle = sprintf("%s | %s", coh_disp, pstar(pval)),
         x = NULL, y = ylab, caption = note) +
    theme_pub() + theme(plot.caption = element_text(size = 10.5, hjust = 0, color = "grey35"))
}
km_paper <- function(d, var, feat_label, coh_disp, pval, hr, lo, hi, note) {
  cut <- median(d[[var]])
  d$expr_group <- factor(ifelse(d[[var]] > cut, "High", "Low"), levels = c("Low", "High"))
  fit <- survival::survfit(Surv(os_days, event) ~ expr_group, data = d)
  ptxt <- sprintf("%s\nHR = %.2f (%.2f-%.2f)", pstar(pval), hr, lo, hi)
  gg <- survminer::ggsurvplot(fit, data = d, palette = c(COL_GREEN, COL_RED),
    legend.title = "Expression", legend.labs = c("Low", "High"), legend = "top",
    pval = ptxt, pval.size = 4.6, pval.coord = c(0.02 * max(d$os_days), 0.12),
    risk.table = TRUE, risk.table.height = 0.27, risk.table.fontsize = 4.6,
    risk.table.title = "Number at risk", tables.y.text = FALSE, conf.int = FALSE,
    censor.shape = "|", censor.size = 4, xlab = "Time (days)", ylab = "Overall survival",
    title = feat_label, font.main = c(22, "bold"), font.x = 17, font.y = 17,
    font.tickslab = 15, font.legend = 15, ggtheme = theme_classic(base_size = 16))
  gg$plot <- gg$plot +
    labs(subtitle = coh_disp, caption = note) +
    theme(plot.title = element_text(hjust = 0.5, size = 22, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 15),
          plot.caption = element_text(size = 10.5, hjust = 0, color = "grey35"))
  list(curve = gg$plot,
       full = cowplot::plot_grid(gg$plot, gg$table, ncol = 1,
                                 rel_heights = c(3.3, 1.1), align = "v", axis = "lr"))
}
NOTE_TIGER  <- "All samples (pre + on-treatment); P corrected for repeated\nbiopsies (mixed-model / cluster-robust)"
NOTE_BASE   <- "Pre-treatment baseline, 1 sample/patient (rigorous comparator)"

# =============================================================================
# MAIN  (all TIGER-style, all significant after correction)
# =============================================================================
mainA <- box_paper(subset_data("tiger","pooled"), "Degradation_score", "Degradation score",
                   "All ICI (n=91)", pget(tigResp,"Degradation score","All (91)","p_mixedmodel_corrected"), NOTE_TIGER)
mainB <- box_paper(subset_data("tiger","antiPD1"), "SARDH", "SARDH",
                   "anti-PD-1 (n=50)", pget(tigResp,"SARDH","anti-PD-1","p_mixedmodel_corrected"), NOTE_TIGER)
mainC <- box_paper(subset_data("tiger","pooled"), "PIPOX", "PIPOX",
                   "All ICI (n=91)", pget(tigResp,"PIPOX","All (91)","p_mixedmodel_corrected"), NOTE_TIGER)
mainD <- km_paper(subset_data("tiger","pooled"), "SARDH", "SARDH", "All ICI (n=91)",
                  pget(tigSurv,"SARDH","All (91)","p_clusterrobust_corrected"),
                  pget(tigSurv,"SARDH","All (91)","HR_high_vs_low"),
                  pget(tigSurv,"SARDH","All (91)","HR_lower"),
                  pget(tigSurv,"SARDH","All (91)","HR_upper"), NOTE_TIGER)
# Panel letters follow the user's Figure 1 layout (A,B already used -> start at C)
save_ms(mainA, "Main", "Main_C_Degradation_score_RvsNR_allICI", 5.3, 5.9)
save_ms(mainB, "Main", "Main_D_SARDH_RvsNR_antiPD1",           5.3, 5.9)
save_ms(mainC, "Main", "Main_E_PIPOX_RvsNR_allICI",            5.3, 5.9)
save_ms(mainD$full, "Main", "Main_F_SARDH_OS_allICI",          6.0, 7.2)
# combined 2x2 main panel (letters C-F; D as curve-only to fit grid)
save_ms(cowplot::plot_grid(mainA, mainB, mainC, mainD$curve, ncol = 2,
        labels = c("C","D","E","F"), label_size = 20), "Main", "Main_panel_CDEF", 11, 12.4)

# =============================================================================
# SUPPLEMENTARY
# =============================================================================
# (1) Production arm, R vs NR (TIGER, pooled) — honestly NS
sP1 <- box_paper(subset_data("tiger","pooled"), "Production_score", "Production score",
                 "All ICI (n=91)", pget(tigResp,"Production score","All (91)","p_mixedmodel_corrected"), NOTE_TIGER)
sP2 <- box_paper(subset_data("tiger","pooled"), "GNMT", "GNMT",
                 "All ICI (n=91)", pget(tigResp,"GNMT","All (91)","p_mixedmodel_corrected"), NOTE_TIGER)
sP3 <- box_paper(subset_data("tiger","pooled"), "DMGDH", "DMGDH",
                 "All ICI (n=91)", pget(tigResp,"DMGDH","All (91)","p_mixedmodel_corrected"), NOTE_TIGER)
# Supplementary letters start at B (panel A already used in the user's Sup Fig 1)
save_ms(sP1, "Supplementary", "Sup_B_Production_score_RvsNR_allICI")
save_ms(sP2, "Supplementary", "Sup_C_GNMT_RvsNR_allICI")
save_ms(sP3, "Supplementary", "Sup_D_DMGDH_RvsNR_allICI")
save_ms(cowplot::plot_grid(sP1, sP2, sP3, ncol = 3, labels = c("B","C","D"), label_size = 18),
        "Supplementary", "Sup_panel_BCD_production_RvsNR", 15, 5.9)

# NOTE (per user curation): the supplementary SURVIVAL-TREND panels (Production/
# Degradation/PIPOX OS) and the BASELINE-robustness panels (PRE-only) were removed
# from the manuscript set to keep Supplementary = the production arm only (B-D).
# Those analyses still exist in results/figures/ (TIGER_all_samples & primary_PRE_only)
# and the baseline sensitivity remains described in the Methods/Results text.

cat("Manuscript figures written to results/manuscript_figures/{Main,Supplementary}/\n")
cat("Main:", length(list.files("results/manuscript_figures/Main")), "files | ",
    "Supplementary:", length(list.files("results/manuscript_figures/Supplementary")), "files\n")
