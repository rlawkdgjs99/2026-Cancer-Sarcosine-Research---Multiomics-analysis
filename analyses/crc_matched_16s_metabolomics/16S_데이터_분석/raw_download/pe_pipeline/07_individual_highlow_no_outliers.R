#!/usr/bin/env Rscript
# ============================================================
#  07_individual_highlow_no_outliers.R
#   06과 동일(균종별 High/Low -> sarcosine)하되 sarcosine 이상치 제거 버전
#   이상치 정의: log10(sarcosine) 전체 분포의 Tukey IQR(1.5x), 한 번만 적용(비교와 독립)
#   주의: Wilcoxon은 순위기반이라 이상치에 강건 -> 큰 변화 기대 어려움. 탐색용.
#  출력: 07a/07b PNG, 07_highlow_outlier_comparison.tsv (원본 vs 제거 나란히)
# ============================================================
suppressMessages(library(ggplot2)); set.seed(42)
OUT <- "16S_데이터_분석/raw_download/pe_pipeline/sarcosine"   # relative to the analysis-folder root (working dir)
df0 <- read.delim(file.path(OUT, "merged_per_sample.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
df0$group <- factor(df0$group, levels = c("Healthy", "CRC"))
degraders <- c("Lachnospira", "Faecalibacterium", "Coprococcus", "Fusicatenibacter", "Roseburia", "Blautia")
COLHL <- c("Low" = "#D95F02", "High" = "#2CA25F")
theme_pub <- theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 20, margin = margin(b = 3)),
        plot.subtitle = element_text(size = 14, colour = "grey30", margin = margin(b = 8)),
        axis.title = element_text(face = "bold"), strip.text = element_text(face = "bold", size = 14),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "none", plot.margin = margin(12, 16, 10, 12))
savepng <- function(p, f, w, h) ggsave(file.path(OUT, f), p, width = w, height = h, dpi = 300, bg = "white", device = grDevices::png)

## ---------- 이상치 제거 (log10 sarcosine Tukey IQR) ----------
ls <- log10(df0$sarcosine); q <- quantile(ls, c(.25, .75)); iqr <- q[2] - q[1]
lo <- q[1] - 1.5 * iqr; hi <- q[2] + 1.5 * iqr
is_out <- ls < lo | ls > hi
df <- df0[!is_out, ]
cat(sprintf("[이상치] log10(sarcosine) Tukey IQR로 %d개 제거 -> n: %d -> %d (Healthy %d / CRC %d)\n",
            sum(is_out), nrow(df0), nrow(df), sum(df$group == "Healthy"), sum(df$group == "CRC")))

hl_test <- function(x, y) { hl <- ifelse(x >= median(x), "High", "Low"); if (length(unique(hl)) < 2) return(NA_real_); suppressWarnings(wilcox.test(y ~ factor(hl))$p.value) }
mkres <- function(d) sapply(degraders, function(g) c(
  pooled  = hl_test(d[[g]], d$sarcosine),
  Healthy = hl_test(d[[g]][d$group == "Healthy"], d$sarcosine[d$group == "Healthy"]),
  CRC     = hl_test(d[[g]][d$group == "CRC"], d$sarcosine[d$group == "CRC"])))

## ---------- 원본 vs 이상치제거 비교 ----------
r0 <- mkres(df0); r1 <- mkres(df)
comp <- data.frame(genus = degraders,
  pooled_orig = round(r0["pooled", ], 4),  pooled_noOut = round(r1["pooled", ], 4),
  withinCRC_orig = round(r0["CRC", ], 4),  withinCRC_noOut = round(r1["CRC", ], 4),
  withinHealthy_orig = round(r0["Healthy", ], 4), withinHealthy_noOut = round(r1["Healthy", ], 4))
write.table(comp, file.path(OUT, "07_highlow_outlier_comparison.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\n=== 균종별 Wilcoxon p: 원본 vs 이상치제거 ===\n"); print(comp, row.names = FALSE)

## ---------- 플롯 (이상치 제거 데이터) ----------
mklong <- function(d) do.call(rbind, lapply(degraders, function(g) {
  x <- d[[g]]; data.frame(genus = g, highlow = factor(ifelse(x >= median(x), "High", "Low"), levels = c("Low", "High")), sarcosine = d$sarcosine) }))
mkplot <- function(d, pcol, ptitle, psub, fn) {
  l <- mklong(d)
  lab <- setNames(sprintf("%s\n(p = %.2g)", degraders, pcol), degraders)
  l$facet <- factor(lab[l$genus], levels = lab[degraders])
  g <- ggplot(l, aes(highlow, sarcosine, fill = highlow)) +
    geom_boxplot(outlier.shape = NA, alpha = .75, width = .6) +
    geom_jitter(width = .14, size = .5, alpha = .3) +
    facet_wrap(~facet, ncol = 3) + scale_y_log10() + scale_fill_manual(values = COLHL) +
    labs(title = ptitle, subtitle = psub, x = "Abundance of that genus (median split)",
         y = "Sarcosine intensity (a.u., log scale)") + theme_pub
  savepng(g, fn, 12, 7)
}
mkplot(df, r1["pooled", ],
       "Fecal sarcosine by each degrader genus -- pooled, outliers removed",
       sprintf("log10-sarcosine Tukey-IQR outliers removed (n=%d). CAUTION pooled split tracks cancer status.", nrow(df)),
       "07a_each_degrader_pooled_noOutliers.png")
crc <- df[df$group == "CRC", ]
mkplot(crc, r1["CRC", ],
       "Fecal sarcosine by each degrader genus -- within CRC, outliers removed",
       sprintf("Within-CRC median split, outliers removed (n=%d); controls for cancer status.", nrow(crc)),
       "07b_each_degrader_within_CRC_noOutliers.png")
cat("\nNO_OUTLIER_DONE\n")
