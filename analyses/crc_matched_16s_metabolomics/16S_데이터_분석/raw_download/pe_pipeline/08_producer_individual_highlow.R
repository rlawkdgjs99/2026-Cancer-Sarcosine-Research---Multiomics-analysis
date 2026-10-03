#!/usr/bin/env Rscript
# ============================================================
#  08_producer_individual_highlow.R
#   생성(producer) 균종 개별 High/Low -> sarcosine (06과 동일 방식)
#   producers: Hungatella, Clostridium sensu stricto 1
#   적응적 분할: median>0 이면 High/Low(중앙값), 아니면(희박) Present/Absent(검출)
#      - Hungatella prev 43% -> Present/Absent
#      - Clostridium s.s.1 prev 76% -> 중앙값 분할
#   가설: 생성균 많을수록 sarcosine 높음(Higher>Lower) — 분해균과 반대
#   주의: Clostridium s.s.1 != 반드시 C.symbiosum (Silva 매핑 불확실)
# ============================================================
suppressMessages(library(ggplot2)); set.seed(42)
OUT <- "16S_데이터_분석/raw_download/pe_pipeline/sarcosine"   # relative to the analysis-folder root (working dir)
df <- read.delim(file.path(OUT, "merged_per_sample.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
df$group <- factor(df$group, levels = c("Healthy", "CRC"))
producers <- c("Hungatella")   # Clostridium s.s.1 제거: C.symbiosum 아님
COLHL <- c("Lower" = "#2C7FB8", "Higher" = "#7B3294")
theme_pub <- theme_bw(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 19, margin = margin(b = 3)),
        plot.subtitle = element_text(size = 14, colour = "grey30", margin = margin(b = 8)),
        axis.title = element_text(face = "bold"), strip.text = element_text(face = "bold", size = 14),
        strip.background = element_rect(fill = "grey92", colour = NA),
        legend.position = "none", plot.margin = margin(12, 16, 10, 12))
savepng <- function(p, f, w, h) ggsave(file.path(OUT, f), p, width = w, height = h, dpi = 300, bg = "white", device = grDevices::png)

# 적응적 분할: list(factor Lower/Higher, type)
split_info <- function(x) {
  if (median(x) > 0) list(f = factor(ifelse(x >= median(x), "Higher", "Lower"), levels = c("Lower", "Higher")), type = "median split")
  else               list(f = factor(ifelse(x > 0, "Higher", "Lower"), levels = c("Lower", "Higher")), type = "present/absent")
}
hl_test <- function(x, y) {
  s <- split_info(x)$f
  if (nlevels(droplevels(s)) < 2 || min(table(s)) < 5) return(NA_real_)
  suppressWarnings(wilcox.test(y ~ s)$p.value)
}

## ---------- 통계 (pooled / within Healthy / within CRC) + 분할유형 ----------
res <- do.call(rbind, lapply(producers, function(g) {
  x <- df[[g]]
  data.frame(genus = g, split_type = split_info(x)$type,
    prevalence_pct = round(100 * mean(x > 0), 1),
    p_pooled = hl_test(x, df$sarcosine),
    p_within_Healthy = hl_test(df[[g]][df$group == "Healthy"], df$sarcosine[df$group == "Healthy"]),
    p_within_CRC = hl_test(df[[g]][df$group == "CRC"], df$sarcosine[df$group == "CRC"]))
}))
write.table(res, file.path(OUT, "08_producer_highlow_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("=== producer individual High/Low -> sarcosine ===\n"); print(res, row.names = FALSE)

## ---------- 플롯 ----------
mklong <- function(d) do.call(rbind, lapply(producers, function(g) {
  s <- split_info(d[[g]])
  data.frame(genus = g, split = s$f, sarcosine = d$sarcosine) }))
mkplot <- function(d, pvec, ptitle, psub, fn) {
  l <- mklong(d)
  lab <- setNames(sprintf("%s  [%s]\n(p = %.2g)", res$genus, res$split_type, pvec), producers)
  l$facet <- factor(lab[l$genus], levels = lab[producers])
  g <- ggplot(l, aes(split, sarcosine, fill = split)) +
    geom_boxplot(outlier.shape = NA, alpha = .75, width = .55) +
    geom_jitter(width = .14, size = .6, alpha = .35) +
    facet_wrap(~facet, ncol = 2) + scale_y_log10() + scale_fill_manual(values = COLHL) +
    labs(title = ptitle, subtitle = psub,
         x = "Producer abundance (Lower vs Higher)", y = "Sarcosine intensity (a.u., log scale)") + theme_pub
  savepng(g, fn, 9, 5.5)
}
mkplot(df, res$p_pooled,
       "Fecal sarcosine by producer genus -- pooled (n=308)",
       "Hypothesis: producers raise sarcosine (Higher>Lower). CAUTION: pooled tracks cancer.",
       "08a_each_producer_pooled.png")
crc <- df[df$group == "CRC", ]
mkplot(crc, res$p_within_CRC,
       "Fecal sarcosine by producer genus -- within CRC only",
       sprintf("Within-CRC split (n=%d) controls for cancer status.", nrow(crc)),
       "08b_each_producer_within_CRC.png")
cat("\nPRODUCER_HIGHLOW_DONE\n")
