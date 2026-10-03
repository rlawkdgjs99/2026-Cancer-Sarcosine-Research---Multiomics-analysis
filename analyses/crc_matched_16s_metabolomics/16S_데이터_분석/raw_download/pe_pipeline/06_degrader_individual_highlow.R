#!/usr/bin/env Rscript
# ============================================================
#  06_degrader_individual_highlow.R
#   각 degrader 균종을 *개별적으로* High/Low(중앙값 분할) -> sarcosine 비교
#   (합산 degrader_sum 아님! 균종마다 따로)
#   - 06a: pooled (전체) 균종별 High vs Low            [교란 주의]
#   - 06b: within-CRC (암 환자 내) 균종별 High vs Low  [교란 통제]
#   - stats: 균종별 pooled / within-Healthy / within-CRC p
#  입력: sarcosine/merged_per_sample.tsv
# ============================================================
suppressMessages(library(ggplot2)); set.seed(42)
OUT <- "16S_데이터_분석/raw_download/pe_pipeline/sarcosine"   # relative to the analysis-folder root (working dir)
df <- read.delim(file.path(OUT, "merged_per_sample.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
df$group <- factor(df$group, levels = c("Healthy", "CRC"))
degraders <- c("Lachnospira", "Faecalibacterium", "Coprococcus", "Fusicatenibacter", "Roseburia", "Blautia")
COLHL <- c("Low" = "#D95F02", "High" = "#2CA25F")

theme_pub <- theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 20, margin = margin(b = 3)),
        plot.subtitle = element_text(size = 14, colour = "grey30", margin = margin(b = 8)),
        axis.title = element_text(face = "bold"), strip.text = element_text(face = "bold", size = 14),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "none", plot.margin = margin(12, 16, 10, 12))
savepng <- function(p, f, w, h) ggsave(file.path(OUT, f), p, width = w, height = h, dpi = 300, bg = "white", device = grDevices::png)

# 중앙값 분할 후 sarcosine Wilcoxon (한 벡터의 High/Low로 다른 벡터 비교)
hl_test <- function(x, y) {
  hl <- ifelse(x >= median(x), "High", "Low")
  if (length(unique(hl)) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(y ~ factor(hl))$p.value)
}

## ---------- 균종별 통계 (pooled / within Healthy / within CRC) ----------
res <- do.call(rbind, lapply(degraders, function(g) {
  x <- df[[g]]
  data.frame(genus = g,
    p_pooled = hl_test(x, df$sarcosine),
    p_within_Healthy = hl_test(df[[g]][df$group == "Healthy"], df$sarcosine[df$group == "Healthy"]),
    p_within_CRC = hl_test(df[[g]][df$group == "CRC"], df$sarcosine[df$group == "CRC"]),
    median_sarc_Low = median(df$sarcosine[x < median(x)]),
    median_sarc_High = median(df$sarcosine[x >= median(x)]))
}))
write.table(res, file.path(OUT, "06_individual_highlow_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---------- 롱 데이터 (균종별 High/Low) ----------
mklong <- function(d) do.call(rbind, lapply(degraders, function(g) {
  x <- d[[g]]
  data.frame(genus = g, highlow = factor(ifelse(x >= median(x), "High", "Low"), levels = c("Low", "High")),
             sarcosine = d$sarcosine)
}))

## ---------- 06a: pooled ----------
la <- mklong(df)
labA <- setNames(sprintf("%s\n(pooled p = %.2g)", res$genus, res$p_pooled), res$genus)
la$facet <- factor(labA[la$genus], levels = labA[degraders])
ga <- ggplot(la, aes(highlow, sarcosine, fill = highlow)) +
  geom_boxplot(outlier.shape = NA, alpha = .75, width = .6) +
  geom_jitter(width = .14, size = .5, alpha = .3) +
  facet_wrap(~facet, ncol = 3) + scale_y_log10() + scale_fill_manual(values = COLHL) +
  labs(title = "Fecal sarcosine by High vs Low of EACH degrader genus (pooled, all 308)",
       subtitle = "Per-genus median split; orange=Low, green=High. CAUTION: pooled split tracks cancer status (confounded).",
       x = "Abundance of that genus (median split)", y = "Sarcosine intensity (a.u., log scale)") +
  theme_pub
savepng(ga, "06a_sarcosine_by_each_degrader_pooled.png", 12, 7)

## ---------- 06b: within CRC (암 환자 내, 교란 통제) ----------
crc <- df[df$group == "CRC", ]
lb <- mklong(crc)
labB <- setNames(sprintf("%s\n(within-CRC p = %.2g)", res$genus, res$p_within_CRC), res$genus)
lb$facet <- factor(labB[lb$genus], levels = labB[degraders])
gb <- ggplot(lb, aes(highlow, sarcosine, fill = highlow)) +
  geom_boxplot(outlier.shape = NA, alpha = .75, width = .6) +
  geom_jitter(width = .14, size = .5, alpha = .3) +
  facet_wrap(~facet, ncol = 3) + scale_y_log10() + scale_fill_manual(values = COLHL) +
  labs(title = "Fecal sarcosine by High vs Low of each degrader genus -- within CRC only",
       subtitle = sprintf("Median split within CRC (n=%d) -> controls for cancer status; tests a direct genus->sarcosine effect.", nrow(crc)),
       x = "Abundance of that genus (within-CRC median split)", y = "Sarcosine intensity (a.u., log scale)") +
  theme_pub
savepng(gb, "06b_sarcosine_by_each_degrader_within_CRC.png", 12, 7)

cat("=== 균종별 High vs Low -> sarcosine (Wilcoxon p) ===\n")
print(res, row.names = FALSE)
cat("\nINDIVIDUAL_HIGHLOW_DONE\n")
