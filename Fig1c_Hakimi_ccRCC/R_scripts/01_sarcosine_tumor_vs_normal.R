#!/usr/bin/env Rscript
# =====================================================================
# Sarcosine: Tumor vs Normal in clear cell RCC
# Source : Hakimi et al., Cancer Cell 29:104-116 (2016), Suppl. Table S2
#          File: 1-s2.0-S1535610815004687-mmc2.xlsx, sheet "Median Normalized"
# Design : 138 matched tumor/normal pairs (Metabolon GC/MS; sarcosine = GC/MS)
# Tests  : (primary) paired Wilcoxon signed-rank + paired t-test on log2
#          (paper replication) unpaired Mann-Whitney U  [Fig 1B / Table S3]
# Integrity: data extracted by LABEL matching (not hard-coded positions);
#            computed values cross-checked against published S3
#            (log2FC = +1.487, adj.P = 4.53e-18).
# =====================================================================

suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2)
})
set.seed(1)  # reproducible jitter

BASE    <- "."  # run with the R working directory set to the analysis-folder root
xlsx    <- file.path(BASE, "1-s2.0-S1535610815004687-mmc2.xlsx")
out_dir <- file.path(BASE, "Outputs")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- 1. Read 'Median Normalized' sheet as a raw character grid ------
raw <- read_excel(xlsx, sheet = "Median Normalized",
                  col_names = FALSE, col_types = "text", .name_repair = "minimal")
M <- as.matrix(raw)

find_one_row <- function(M, label) {
  hits <- which(apply(M, 1, function(r) any(trimws(r) == label, na.rm = TRUE)))
  if (length(hits) != 1) stop(sprintf("Expected exactly 1 row for '%s', found %d", label, length(hits)))
  hits
}
tt_row    <- find_one_row(M, "TISSUE TYPE")
match_row <- find_one_row(M, "MATCHING INDEX")

sar_hits <- which(apply(M, 1, function(r) any(grepl("sarcosine", r, ignore.case = TRUE))))
if (length(sar_hits) != 1) stop(sprintf("Expected exactly 1 sarcosine row, found %d", length(sar_hits)))
sar_row  <- sar_hits
sar_col  <- which(grepl("sarcosine", M[sar_row, ], ignore.case = TRUE))[1]
sar_name <- trimws(M[sar_row, sar_col])

## ---- 2. Identify sample columns + extract sarcosine ----------------
tissue_all  <- toupper(trimws(M[tt_row, ]))
sample_cols <- which(tissue_all %in% c("TUMOR", "NORMAL"))

df <- data.frame(
  matching  = trimws(M[match_row, sample_cols]),
  tissue    = factor(tissue_all[sample_cols], levels = c("NORMAL", "TUMOR")),
  abundance = as.numeric(M[sar_row, sample_cols]),
  stringsAsFactors = FALSE
)

## ---- 3. Integrity checks -------------------------------------------
cat("Metabolite identified :", sar_name, " (sheet row", sar_row, ")\n")
cat("Tissue counts         :", paste(names(table(df$tissue)), table(df$tissue), collapse = " | "), "\n")
stopifnot(sum(df$tissue == "TUMOR") == 138, sum(df$tissue == "NORMAL") == 138)
cat("Missing sarcosine vals:", sum(is.na(df$abundance)), "(expected 0)\n")
stopifnot(sum(is.na(df$abundance)) == 0, all(df$abundance > 0))
dup <- df %>% count(matching, tissue) %>% filter(n > 1)
stopifnot(nrow(dup) == 0)                       # each (pair,tissue) unique
df$log2abund <- log2(df$abundance)

## ---- 4. Matched pairs by MATCHING INDEX ----------------------------
wide <- df %>% select(matching, tissue, abundance) %>%
  pivot_wider(names_from = tissue, values_from = abundance) %>%
  filter(!is.na(TUMOR), !is.na(NORMAL))
n_pairs <- nrow(wide)
cat("Complete tumor/normal pairs:", n_pairs, "\n")
wide$log2ratio <- log2(wide$TUMOR / wide$NORMAL)

## ---- 5. Statistics --------------------------------------------------
pw <- wilcox.test(wide$TUMOR, wide$NORMAL, paired = TRUE)               # primary
pt <- t.test(log2(wide$TUMOR), log2(wide$NORMAL), paired = TRUE)        # primary (parametric)
mw <- wilcox.test(abundance ~ tissue, data = df)                       # paper: unpaired Mann-Whitney
ut <- t.test(log2abund ~ tissue, data = df)                           # unpaired t on log2

mean_T <- mean(df$abundance[df$tissue == "TUMOR"]);  mean_N <- mean(df$abundance[df$tissue == "NORMAL"])
med_T  <- median(df$abundance[df$tissue == "TUMOR"]); med_N <- median(df$abundance[df$tissue == "NORMAL"])
log2FC_mean   <- log2(mean_T / mean_N)
log2FC_median <- log2(med_T / med_N)
log2FC_paired <- mean(wide$log2ratio)

# format a p-value, flagging R's numerical floor (~2.2e-16) honestly
fmt_p <- function(p) if (p <= 2.23e-16) "p < 2.2e-16" else paste0("p = ", format(p, scientific = TRUE, digits = 2))

## ---- 5b. Sensitivity: complete-case (exclude min-imputed values) ----
# In "Median Normalized" the missing sarcosine values were imputed with the
# per-metabolite minimum; missingness is differential (more in NORMAL), which
# can mildly inflate the fold change. Re-run on pairs where BOTH tumor and
# normal sarcosine were actually observed in the "Raw Data" sheet.
rawg <- read_excel(xlsx, sheet = "Raw Data", col_names = FALSE, col_types = "text", .name_repair = "minimal")
Rm   <- as.matrix(rawg)
tt_r  <- find_one_row(Rm, "TISSUE TYPE")
mt_r  <- find_one_row(Rm, "MATCHING INDEX")
sar_r <- which(apply(Rm, 1, function(r) any(grepl("sarcosine", r, ignore.case = TRUE))))
stopifnot(length(sar_r) == 1)
tis_r  <- toupper(trimws(Rm[tt_r, ]))
scol_r <- which(tis_r %in% c("TUMOR", "NORMAL", "T", "N"))
raw_df <- data.frame(
  matching = trimws(Rm[mt_r, scol_r]),
  tissue   = ifelse(tis_r[scol_r] %in% c("TUMOR", "T"), "TUMOR", "NORMAL"),
  present  = !is.na(suppressWarnings(as.numeric(Rm[sar_r, scol_r]))),
  stringsAsFactors = FALSE)
raw_miss_T <- sum(raw_df$tissue == "TUMOR"  & !raw_df$present)
raw_miss_N <- sum(raw_df$tissue == "NORMAL" & !raw_df$present)
comp_idx <- raw_df %>% group_by(matching) %>%
  summarise(ok = (n() == 2) && all(present), .groups = "drop") %>% filter(ok) %>% pull(matching)
wide_cc  <- wide %>% filter(matching %in% comp_idx)
n_cc     <- nrow(wide_cc)
cc_log2FC <- log2(mean(wide_cc$TUMOR) / mean(wide_cc$NORMAL))
cc_pt     <- t.test(log2(wide_cc$TUMOR), log2(wide_cc$NORMAL), paired = TRUE)
cc_up     <- sum(wide_cc$TUMOR > wide_cc$NORMAL)
cat(sprintf("Raw missing: tumor=%d, normal=%d | complete-case pairs=%d, log2FC(mean)=%.4f, paired-t %s, tumor>normal=%d/%d\n",
            raw_miss_T, raw_miss_N, n_cc, cc_log2FC, fmt_p(cc_pt$p.value), cc_up, n_cc))

## ---- 6. Results table ----------------------------------------------
res <- data.frame(
  metric = c("metabolite", "n_pairs", "n_tumor", "n_normal",
             "mean_tumor", "mean_normal", "median_tumor", "median_normal",
             "log2FC_mean_ratio", "log2FC_median_ratio", "log2FC_paired_mean_log2ratio",
             "paired_wilcoxon_p", "paired_t_test_p",
             "unpaired_mannwhitney_p", "unpaired_t_test_p",
             "raw_missing_tumor", "raw_missing_normal",
             "complete_case_n_pairs", "complete_case_log2FC_mean", "complete_case_paired_t_p",
             "complete_case_pairs_tumor_gt_normal",
             "PUBLISHED_S3_log2FC", "PUBLISHED_S3_adjP"),
  value = c(sar_name, n_pairs, sum(df$tissue == "TUMOR"), sum(df$tissue == "NORMAL"),
            round(mean_T, 4), round(mean_N, 4), round(med_T, 4), round(med_N, 4),
            round(log2FC_mean, 4), round(log2FC_median, 4), round(log2FC_paired, 4),
            signif(pw$p.value, 4), signif(pt$p.value, 4),
            signif(mw$p.value, 4), signif(ut$p.value, 4),
            raw_miss_T, raw_miss_N,
            n_cc, round(cc_log2FC, 4), signif(cc_pt$p.value, 4),
            paste0(cc_up, "/", n_cc),
            1.487, "4.53e-18"),
  stringsAsFactors = FALSE)
write.csv(res, file.path(out_dir, "sarcosine_tumor_vs_normal_stats.csv"), row.names = FALSE)
cat("\n==== RESULTS ====\n"); print(res, row.names = FALSE)

## ---- 7. Figures (PNG, large fonts, concise titles) -----------------
big_theme <- theme_classic(base_size = 18) +
  theme(plot.title    = element_text(size = 24, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5, color = "grey25"),
        axis.text     = element_text(size = 18, color = "black"),
        axis.title    = element_text(size = 18),
        plot.margin   = margin(16, 16, 16, 16))

# Fig 1: boxplot tumor vs normal
p1 <- ggplot(df, aes(tissue, log2abund, fill = tissue)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.65) +
  geom_jitter(width = 0.12, height = 0, size = 1.7, alpha = 0.45) +
  scale_fill_manual(values = c(NORMAL = "#55A868", TUMOR = "#C44E52"), guide = "none") +  # Control=green, Cancer=red (thesis convention)
  labs(title = "Sarcosine: Tumor vs Normal",
       subtitle = sprintf("Paired Wilcoxon %s    |    log2 FC (T/N) = %+.2f",
                          fmt_p(pw$p.value), log2FC_mean),
       x = NULL, y = expression(log[2]*" median-normalized abundance")) +
  big_theme
ggsave(file.path(out_dir, "sarcosine_tumor_vs_normal_boxplot.png"),
       p1, width = 7, height = 6, dpi = 300, bg = "white")

# Fig 2: paired tumor/normal ratio
p2 <- ggplot(wide, aes(x = "", y = log2ratio)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.7) +
  geom_violin(fill = "#C44E52", alpha = 0.30, width = 0.7) +
  geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.08, height = 0, size = 1.6, alpha = 0.45) +
  labs(title = "Sarcosine: Paired Tumor/Normal Ratio",
       subtitle = sprintf("n = %d pairs   |   median log2(T/N) = %+.2f   |   paired %s",
                          n_pairs, median(wide$log2ratio), fmt_p(pw$p.value)),
       x = NULL, y = expression(log[2]*"(Tumor / Normal)")) +
  big_theme
ggsave(file.path(out_dir, "sarcosine_paired_ratio.png"),
       p2, width = 7, height = 6, dpi = 300, bg = "white")

cat("\nSaved -> Outputs/: sarcosine_tumor_vs_normal_stats.csv, ",
    "sarcosine_tumor_vs_normal_boxplot.png, sarcosine_paired_ratio.png\n", sep = "")
cat("R sessionInfo platform:", R.version.string, "\n")
