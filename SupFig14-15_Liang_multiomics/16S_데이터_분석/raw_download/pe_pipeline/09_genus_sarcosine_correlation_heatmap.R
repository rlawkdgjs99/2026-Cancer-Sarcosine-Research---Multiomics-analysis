#!/usr/bin/env Rscript
# ============================================================
#  09_genus_sarcosine_correlation_heatmap.R
#   분해/생성 균(genus) <-> 분변 sarcosine Spearman 상관 heatmap
#   대상 대사체 = sarcosine 하나만 (사용자 주관심사)
#   각 균 x {Overall(pooled) / within-Healthy / within-CRC} 3개 맥락
#  입력: sarcosine/merged_per_sample.tsv
#  NOTE: 폴더명 16S_데이터_분석 (구 HGMT_Data-16S 에서 변경됨) 기준 경로
# ============================================================
suppressMessages(library(ggplot2)); set.seed(42)
OUT <- "16S_데이터_분석/raw_download/pe_pipeline/sarcosine"
df <- read.delim(file.path(OUT, "merged_per_sample.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
df$group <- factor(df$group, levels = c("Healthy", "CRC"))

genera <- data.frame(
  genus = c("Lachnospira", "Faecalibacterium", "Coprococcus", "Fusicatenibacter", "Roseburia",
            "[Eubacterium] eligens group", "Hungatella", "Blautia"),
  disp  = c("Lachnospira", "Faecalibacterium", "Coprococcus", "Fusicatenibacter", "Roseburia",
            "[Eub.] eligens grp", "Hungatella", "Blautia"),
  # Blautia = 'Mixed' (분해속 단정 불가: CRC에서 증가, within-CRC만 음의 상관);
  # [Eubacterium] eligens group = L.eligens 특이(SILVA), 분해속
  role  = c(rep("Degraders (expect negative)", 6), "Producers (expect positive)", "Mixed (genus-level)"),
  stringsAsFactors = FALSE)

contexts <- list(Overall = rep(TRUE, nrow(df)), Healthy = df$group == "Healthy", CRC = df$group == "CRC")

res <- do.call(rbind, lapply(seq_len(nrow(genera)), function(i) {
  g <- genera$genus[i]
  do.call(rbind, lapply(names(contexts), function(cx) {
    idx <- contexts[[cx]]
    ct <- suppressWarnings(cor.test(df[[g]][idx], df$sarcosine[idx], method = "spearman"))
    data.frame(genus = g, disp = genera$disp[i], role = genera$role[i], context = cx,
               n = sum(idx), rho = unname(ct$estimate), p = ct$p.value)
  }))
}))
res$stars <- ifelse(res$p < 0.001, "***", ifelse(res$p < 0.01, "**", ifelse(res$p < 0.05, "*", "")))
write.table(res, file.path(OUT, "09_genus_sarcosine_correlation.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("=== genus <-> sarcosine Spearman ===\n"); print(res[, c("disp","context","n","rho","p","stars")], row.names = FALSE)

## ---------- heatmap ----------
res$disp <- factor(res$disp, levels = rev(genera$disp))
res$context <- factor(res$context, levels = c("Overall", "Healthy", "CRC"))
res$role <- factor(res$role, levels = c("Degraders (expect negative)", "Producers (expect positive)", "Mixed (genus-level)"))
lim <- max(abs(res$rho))

g <- ggplot(res, aes(context, disp, fill = rho)) +
  geom_tile(colour = "white", linewidth = .6) +
  geom_text(aes(label = sprintf("%.2f%s", rho, stars)), size = 5.6) +
  facet_grid(role ~ ., scales = "free_y", space = "free", switch = "y") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                       limits = c(-lim, lim), name = "Spearman\nrho") +
  scale_x_discrete(position = "top", labels = c("Overall\n(pooled, n=308)", "within\nHealthy (121)", "within\nCRC (187)")) +
  labs(title = "Sarcosine-related genera vs. fecal sarcosine (Spearman rho)",
       subtitle = "* p<0.05  ** p<0.01  *** p<0.001.   'Overall' = pooled (confounded by cancer); Healthy/CRC = within-group.",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 16) +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0),
        plot.subtitle = element_text(size = 13, colour = "grey30", hjust = 0, margin = margin(b = 8)),
        axis.text.y = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 14, lineheight = .9),
        strip.text.y.left = element_text(face = "bold", size = 14, angle = 0),
        strip.placement = "outside", panel.grid = element_blank(),
        legend.title = element_text(face = "bold"), plot.margin = margin(12, 14, 12, 12))
ggsave(file.path(OUT, "09_genus_sarcosine_correlation_heatmap.png"), g,
       width = 13.5, height = 6, dpi = 300, bg = "white", device = grDevices::png)
cat("\nCORR_HEATMAP_DONE\n")
