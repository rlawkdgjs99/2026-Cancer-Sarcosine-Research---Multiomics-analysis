#!/usr/bin/env Rscript
# ============================================================
#  05_degrader_highlow.R  —  분해균 High vs Low 에 따른 sarcosine 비교
#   (a) pooled: 전체에서 분해균 High/Low(중앙값 분할) -> sarcosine  [교란 주의]
#   (b) within-group: Healthy/CRC 각각에서 High/Low -> sarcosine    [교란 통제]
#   + 교란 가시화: High/Low 그룹의 암/대조 구성비
#  입력: sarcosine/merged_per_sample.tsv  (04에서 생성)
#  출력: sarcosine/05a_*.png, 05b_*.png, 05_highlow_summary.txt, 05_highlow_stats.tsv
# ============================================================
suppressMessages(library(ggplot2))
set.seed(42)

BASE <- "."   # analysis-folder root (working dir)
OUT  <- file.path(BASE, "16S_데이터_분석/raw_download/pe_pipeline/sarcosine")
df <- read.delim(file.path(OUT, "merged_per_sample.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
df$group <- factor(df$group, levels = c("Healthy", "CRC"))

COLHL <- c("Low degraders" = "#D95F02", "High degraders" = "#2CA25F")  # 적=낮음, 녹=높음
theme_pub <- theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 20, margin = margin(b = 3)),
        plot.subtitle = element_text(size = 14, colour = "grey30", margin = margin(b = 8)),
        axis.title = element_text(face = "bold"), strip.text = element_text(face = "bold", size = 16),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "none", plot.margin = margin(12, 16, 10, 12))
savepng <- function(p, f, w, h) ggsave(file.path(OUT, f), p, width = w, height = h, dpi = 300, bg = "white", device = grDevices::png)

## ---------- (a) pooled 중앙값 분할 ----------
med <- median(df$degrader_sum)
df$deg_pool <- factor(ifelse(df$degrader_sum >= med, "High degraders", "Low degraders"),
                      levels = c("Low degraders", "High degraders"))
p_pool <- wilcox.test(sarcosine ~ deg_pool, data = df)$p.value
mLo <- median(df$sarcosine[df$deg_pool == "Low degraders"])
mHi <- median(df$sarcosine[df$deg_pool == "High degraders"])
# 교란 가시화: 각 분해균 그룹의 암 비율
ct <- table(df$deg_pool, df$group)
pctCRC <- round(100 * ct[, "CRC"] / rowSums(ct), 0)

ga <- ggplot(df, aes(deg_pool, sarcosine, fill = deg_pool)) +
  geom_boxplot(outlier.shape = NA, alpha = .75, width = .6) +
  geom_jitter(width = .15, size = .8, alpha = .4) +
  scale_y_log10() + scale_fill_manual(values = COLHL) +
  labs(title = "Fecal sarcosine by degrader abundance -- pooled",
       subtitle = sprintf("Degrader_sum split; n=%d. Wilcoxon p = %.3g; Low=%.0f, High=%.0f\nCAUTION confounded: Low group %d%% CRC, High group %d%% CRC",
                          nrow(df), p_pool, mLo, mHi,
                          pctCRC[["Low degraders"]], pctCRC[["High degraders"]]),
       x = NULL, y = "Sarcosine intensity (a.u., log scale)") + theme_pub
savepng(ga, "05a_sarcosine_by_degrader_pooled.png", 8, 6.2)

## ---------- (b) within-group 중앙값 분할 (교란 통제) ----------
df$deg_within <- NA_character_
for (g in levels(df$group)) {
  idx <- df$group == g
  df$deg_within[idx] <- ifelse(df$degrader_sum[idx] >= median(df$degrader_sum[idx]), "High", "Low")
}
df$deg_within <- factor(df$deg_within, levels = c("Low", "High"))
p_within <- sapply(levels(df$group), function(g) {
  d <- df[df$group == g, ]; wilcox.test(d$sarcosine ~ d$deg_within)$p.value
})
lab.g <- setNames(sprintf("%s  (Wilcoxon p = %.2g)", levels(df$group), p_within), levels(df$group))
df$facet <- factor(lab.g[as.character(df$group)], levels = lab.g[levels(df$group)])

gb <- ggplot(df, aes(deg_within, sarcosine, fill = deg_within)) +
  geom_boxplot(outlier.shape = NA, alpha = .75, width = .6) +
  geom_jitter(width = .15, size = .8, alpha = .4) +
  facet_wrap(~facet) + scale_y_log10() +
  scale_fill_manual(values = c("Low" = "#D95F02", "High" = "#2CA25F")) +
  labs(title = "Fecal sarcosine by degrader abundance, by disease status",
       subtitle = sprintf("Within-group median split, controlling for cancer status (n=%d)", nrow(df)),
       x = "Degrader abundance (within-group median split)", y = "Sarcosine intensity (a.u., log scale)") +
  theme_pub
savepng(gb, "05b_sarcosine_by_degrader_within_group.png", 9, 6)

## ---------- 통계 테이블 + 요약 ----------
stat <- data.frame(
  comparison = c("pooled High vs Low", "within Healthy High vs Low", "within CRC High vs Low"),
  median_sarcosine_Low = c(mLo, median(df$sarcosine[df$group=="Healthy" & df$deg_within=="Low"]),
                           median(df$sarcosine[df$group=="CRC" & df$deg_within=="Low"])),
  median_sarcosine_High = c(mHi, median(df$sarcosine[df$group=="Healthy" & df$deg_within=="High"]),
                            median(df$sarcosine[df$group=="CRC" & df$deg_within=="High"])),
  wilcoxon_p = c(p_pool, p_within[["Healthy"]], p_within[["CRC"]]))
write.table(stat, file.path(OUT, "05_highlow_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sink(file.path(OUT, "05_highlow_summary.txt"))
cat("=== 분해균 High vs Low 에 따른 sarcosine ===\n")
cat(sprintf("페어드 n=%d (Healthy %d / CRC %d)\n\n", nrow(df), sum(df$group=="Healthy"), sum(df$group=="CRC")))
cat("[교란 가시화] pooled 중앙값 분할 시 각 분해균 그룹의 암 비율:\n")
print(addmargins(ct)); cat(sprintf("  -> Low-degrader=%d%% CRC, High-degrader=%d%% CRC (분해균 분할이 사실상 암 상태를 가름)\n\n",
    pctCRC[["Low degraders"]], pctCRC[["High degraders"]]))
cat("[통계]\n"); print(stat, row.names = FALSE)
cat("\n[해석] pooled에서 차이가 나도 그건 대부분 암 상태 교란.",
    "그룹 내(within Healthy/CRC)에서도 High<Low(유의)면 직접 관계 시사; 비유의면 교란으로 해석.\n")
sink()
cat("HIGHLOW_DONE\n")
print(stat)
